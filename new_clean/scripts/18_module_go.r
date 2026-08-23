# =============================================================================
# 18_module_go.r — one GO enrichment per nitrogen-responsive module
#
# 09_go_enrichment.r asks what the CONSERVED GENE SET is for, one test per
# species. It cannot say what any individual module is for. This does: every
# responsive module is its own gene set, tested against the same network-node
# background, so a module can be named -- "Module_026 is photosynthesis" --
# rather than only counted.
#
# THE MODULE SET IS EVERY RESPONSIVE MODULE. With the response called by Spearman
# alone there are no response classes to pool -- what used to be `both` /
# `mi_only` / `pearson_only` is one set. `direction` (rises or falls with
# nitrogen) rides through as a column so it can be split afterwards, but it does
# not partition the run: the question here is what responsive modules do, not
# what separates the two halves.
#
# ONE topGOdata OBJECT, REUSED. 09 builds a fresh one per gene set, which re-runs
# the DAG mapping every time; at ~650 modules that is hours. `updateGenes()`
# swaps the gene list and keeps the graph. The reason to prefer it is not only
# speed: every module is then scored against an IDENTICAL term universe, so the
# per-module results are comparable to each other. Rebuilding per module would
# let the tested term set drift with the gene set.
#
# SCORES, NOT GenTable. `score(runTest(...))` gives the p-value for every tested
# term directly. GenTable formats them into strings ("< 1e-30", 4 significant
# digits) which 09 then has to parse back, and it is the slowest call in the
# loop. termStat() supplies Annotated/Significant/Expected for the terms that
# survive. Same numbers, exactly, without the round trip through text.
#
# THRESHOLD is the RAW weight01 p, per 09's reasoning: weight01 conditions each
# term on its DAG neighbours, so the terms are not an exchangeable family and BH
# does not apply to them. Two BH columns are written for reference and neither
# selects -- see the note above p.adj_global below.
#
# TWO GRAINS OF OUTPUT. The joined table across all modules is what you read to
# see whether the responsive set has a coherent function; the per-module
# directories are what you read to decide whether ONE module is worth following
# up. Those are different questions and they do not fit in one file, so each gets
# its own artefact and its own plot:
#
#   module_go/module_GO_<ONT>_<study>.tsv        joined, every module
#   module_go/..._global.{png,pdf}               GLOBAL: which terms recur
#   module_go/modules/<Module>/GO_..._<mod>.tsv  that module alone
#   module_go/modules/<Module>/GO_..._<mod>.png  GRAIN: that module's terms
#
# RUN: through run.sh  ->  ./run.sh modulego sugarcane
# =============================================================================

suppressMessages({
  library(topGO)
  library(data.table)
  library(ggplot2)
})

source(file.path(dirname(sub("--file=", "",
       grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "lib", "common.R"))

STUDY      <- env_req("CLEAN_STUDY")
PROFILE    <- env_req("CLEAN_MODULE_PROFILE")
MEMBERSHIP <- env_req("CLEAN_MEMBERSHIP")
NODES      <- env_req("CLEAN_NODE_METRICS")
EMAPPER    <- env_req("CLEAN_EMAPPER")
ONTOLOGY   <- env_opt("CLEAN_ONTOLOGY", "BP")
GO_P       <- env_num("CLEAN_GO_P", 0.05)
MIN_ANN    <- as.integer(env_num("CLEAN_MODULE_GO_MIN_ANNOTATED", 3))
NTOP       <- as.integer(env_num("CLEAN_GO_NTOP", 20))   # terms drawn per figure
CORES      <- as.integer(env_num("CLEAN_MODULE_GO_CORES", 1))
LIMIT      <- as.integer(env_num("CLEAN_MODULE_GO_LIMIT", 0))   # 0 = all; smoke tests
OUT_DIR    <- ensure_dir(env_req("CLEAN_OUT_DIR"))
setDTthreads(as.integer(env_num("CLEAN_CORES", 100)))

banner(paste0("per-module GO ", ONTOLOGY, ": ", STUDY))

for (f in c(PROFILE, MEMBERSHIP, NODES, EMAPPER))
  if (!file.exists(f)) stop("missing ", basename(f), call. = FALSE)

# --- annotation and universe -------------------------------------------------
# Both helpers are carried over from 09_go_enrichment.r unchanged. The eggNOG
# layout (col 1 = query id, col 10 = comma-separated GOs) and the isoform merge
# are the same file and the same problem, so this must not be a second
# implementation that can drift from it.
parse_eggnog <- function(annotation_file) {
  raw <- read.table(
    annotation_file, sep = "\t", header = FALSE,
    comment.char = "", quote = "", fill = TRUE, stringsAsFactors = FALSE
  )
  raw <- raw[!grepl("^##",     raw[[1]]), ]
  raw <- raw[!grepl("^#query", raw[[1]]), ]

  gene_col <- raw[[1]]
  go_col   <- raw[[10]]

  go_list <- strsplit(go_col, ",", fixed = TRUE)
  go_list <- lapply(go_list, function(x) {
    x <- trimws(x)
    x[x != "-" & x != "" & grepl("^GO:", x)]
  })

  has_go         <- sapply(go_list, length) > 0
  gene2GO        <- go_list[has_go]
  names(gene2GO) <- gene_col[has_go]

  gene2GO_merged <- tapply(
    seq_along(gene2GO), names(gene2GO),
    function(idx) unique(unlist(gene2GO[idx]))
  )
  as.list(gene2GO_merged)
}

say("parsing eggNOG annotation ...")
gene2GO_all <- parse_eggnog(EMAPPER)
say("  genome GO-annotated genes: ", fmt_n(length(gene2GO_all)))

# Universe = GO-annotated genes of THIS network, identical to 09's background and
# to the TF hypergeometric universe in 15_module_profile.r. Keeping all three on
# the same denominator is what makes "this module is TF-rich" and "this module is
# enriched for photosynthesis" statements about the same population.
node_ids <- read.table(NODES, sep = "\t", header = TRUE,
                       stringsAsFactors = FALSE, quote = "")[[1]]
node_ids     <- unique(trimws(node_ids))
geneUniverse <- intersect(node_ids, names(gene2GO_all))
gene2GO      <- gene2GO_all[geneUniverse]
say("  network nodes: ", fmt_n(length(node_ids)),
    "  |  GO-annotated (= background): ", fmt_n(length(geneUniverse)))
if (!length(geneUniverse)) stop("empty GO background", call. = FALSE)

# --- which modules -----------------------------------------------------------
prof <- fread(PROFILE)
sel  <- prof[responsive == TRUE]
setorder(sel, padj)
say("responsive modules: ", fmt_n(nrow(sel)), "  (",
    paste(sprintf("%s %d", sel[, .N, by = direction]$direction,
                  sel[, .N, by = direction]$N), collapse = " | "), ")")
if (!nrow(sel)) { say("nothing to test"); quit(save = "no", status = 0) }

mem <- fread(MEMBERSHIP, select = c("gene", "module_name"))
mem[, gene := strip_version(gene)]
setnames(mem, "module_name", "module")
mem <- mem[module %chin% sel$module]
by_module <- split(mem$gene, mem$module)

# Annotated members decide whether a module is testable at all. A 3-gene module
# CAN reach p < 0.05 against a 25,000-gene universe, so this gate is not a
# judgement on small modules -- it only skips ones with too few annotated members
# for the test to be defined. n_annotated is reported per module so a reader can
# filter harder without a re-run.
sel[, n_annotated := vapply(module, function(m)
      sum(by_module[[m]] %in% geneUniverse), integer(1))]
testable <- sel[n_annotated >= MIN_ANN]
say(sprintf("annotated members >= %d: %s of %s modules testable  (%s gated out)",
            MIN_ANN, fmt_n(nrow(testable)), fmt_n(nrow(sel)),
            fmt_n(nrow(sel) - nrow(testable))))
say(sprintf("  annotated members: median %.0f  max %.0f",
            median(sel$n_annotated), max(sel$n_annotated)))
if (LIMIT > 0) {
  testable <- head(testable, LIMIT)
  say("CLEAN_MODULE_GO_LIMIT set — testing only the first ", nrow(testable))
}
if (!nrow(testable)) { say("no module clears the annotation gate"); quit(save = "no", status = 0) }

# --- the shared topGOdata object ---------------------------------------------
# Seeded with the union of every testable module's genes purely so the factor has
# both levels; updateGenes replaces it before anything is scored. nodeSize stays
# at topGO's default, as in 09 -- raising it here would make the module-level and
# gene-level term universes different and the two analyses incomparable.
seed_genes <- unique(unlist(by_module[testable$module], use.names = FALSE))
mk_list <- function(genes) {
  gl <- factor(as.integer(geneUniverse %in% genes), levels = c(0L, 1L))
  names(gl) <- geneUniverse
  gl
}
say("building the topGO graph once (", ONTOLOGY, ") ...")
GOdata <- suppressMessages(
  new("topGOdata", ontology = ONTOLOGY, allGenes = mk_list(seed_genes),
      annot = annFUN.gene2GO, gene2GO = gene2GO))
all_go <- usedGO(GOdata)
say("  terms in the shared universe: ", fmt_n(length(all_go)))

# --- per-module test ---------------------------------------------------------
go_terms <- suppressMessages(AnnotationDbi::Term(GO.db::GOTERM[all_go]))

# --- direction-level test ----------------------------------------------------
# Three grains, not two. The per-module test says what ONE module is for; the
# global figure says which terms recur across modules. Neither answers "what do
# the modules that go UP with nitrogen do, as a set, and is it different from the
# ones that go DOWN" -- and that is the question a reader asks first, because the
# two directions are different biology rather than two ways of detecting the
# same thing.
#
# Pooling also breaks the constraint that dominates the per-module analysis. Most
# responsive modules are too small to carry MIN_ANN annotated genes and are gated
# out individually; pooled by direction they all contribute, so this test sees
# the whole responsive set rather than the large-module tail of it.
#
# Same GOdata object, same background, same weight01 statistic and same raw-p
# threshold as every other GO test in this project, so the three grains are
# directly comparable.
say("")
banner(paste0("direction-level ", ONTOLOGY, " — ", STUDY))
dir_sets <- split(sel$module, sel$direction)
dir_rows <- rbindlist(lapply(names(dir_sets), function(d) {
  genes <- unique(unlist(by_module[dir_sets[[d]]], use.names = FALSE))
  n_ann <- sum(genes %in% geneUniverse)
  say(sprintf("  %-8s %s modules, %s genes, %s GO-annotated",
              d, fmt_n(length(dir_sets[[d]])), fmt_n(length(genes)), fmt_n(n_ann)))
  if (n_ann < MIN_ANN) { say("    below the annotation gate — skipped"); return(NULL) }
  gd <- suppressMessages(updateGenes(GOdata, mk_list(genes)))
  rt <- suppressMessages(runTest(gd, algorithm = "weight01", statistic = "fisher"))
  p  <- score(rt)
  padj <- p.adjust(p, method = "BH")
  keep <- names(p)[p <= GO_P]
  say(sprintf("    %s terms at raw p <= %.2g  (%s clearing BH within direction)",
              fmt_n(length(keep)), GO_P, fmt_n(sum(padj <= 0.05))))
  if (!length(keep)) return(NULL)
  st <- suppressMessages(termStat(gd, keep))
  data.table(direction = d, n_modules = length(dir_sets[[d]]),
             n_genes = length(genes), n_annotated = n_ann,
             GO.ID = keep, Term = unname(go_terms[keep]),
             Annotated = st$Annotated, Significant = st$Significant,
             Expected = round(st$Expected, 3),
             pvalue = unname(p[keep]),
             p.adj = signif(unname(padj[keep]), 4))
}))

dir_tag <- file.path(OUT_DIR, sprintf("module_GO_%s_%s_bydirection", ONTOLOGY, STUDY))
if (nrow(dir_rows)) {
  setorder(dir_rows, direction, pvalue)
  write_tsv(dir_rows, paste0(dir_tag, ".tsv"))
  # Which terms are direction-SPECIFIC and which are shared is the whole point of
  # splitting, so it is written as its own small table rather than left for a
  # reader to derive.
  wide <- dcast(dir_rows, GO.ID + Term ~ direction, value.var = "pvalue")
  dirs_present <- setdiff(names(wide), c("GO.ID", "Term"))
  wide[, class := if (length(dirs_present) < 2) "single direction tested" else
        fifelse(!is.na(get(dirs_present[1])) & !is.na(get(dirs_present[2])), "shared",
        fifelse(!is.na(get(dirs_present[1])), paste0(dirs_present[1], " only"),
                paste0(dirs_present[2], " only")))]
  setorder(wide, class, Term)
  write_tsv(wide, paste0(dir_tag, "_comparison.tsv"))
  say("")
  print(wide[, .N, by = class], row.names = FALSE)
} else {
  say("no direction-level term cleared p <= ", GO_P)
  write_tsv(data.table(direction = character(), n_modules = integer(),
                       n_genes = integer(), n_annotated = integer(),
                       GO.ID = character(), Term = character(),
                       Annotated = integer(), Significant = integer(),
                       Expected = numeric(), pvalue = numeric(), p.adj = numeric()),
            paste0(dir_tag, ".tsv"))
}
say("")

test_one <- function(i) {
  mod   <- testable$module[i]
  genes <- by_module[[mod]]
  res <- tryCatch({
    gd  <- suppressMessages(updateGenes(GOdata, mk_list(genes)))
    rt  <- suppressMessages(runTest(gd, algorithm = "weight01", statistic = "fisher"))
    p   <- score(rt)
    keep <- names(p)[p <= GO_P]
    list(p = p, gd = gd, keep = keep)
  }, error = function(e) { say("  ", mod, ": topGO error — ", conditionMessage(e)); NULL })
  if (is.null(res)) return(NULL)

  p <- res$p
  # BH within the module, over EVERY tested term. 09's note applies verbatim:
  # correcting a pre-filtered subset shrinks m to terms already known to be small
  # and makes the adjusted values anti-conservative.
  padj_local <- p.adjust(p, method = "BH")

  rows <- NULL
  if (length(res$keep)) {
    st <- suppressMessages(termStat(res$gd, res$keep))
    rows <- data.table(
      module      = mod,
      GO.ID       = res$keep,
      Term        = unname(go_terms[res$keep]),
      Annotated   = st$Annotated,
      Significant = st$Significant,
      Expected    = round(st$Expected, 3),
      pvalue      = unname(p[res$keep]),
      p.adj       = signif(unname(padj_local[res$keep]), 4),
      # position of each retained term inside this module's full p-vector, so the
      # cross-module BH below can be computed on the complete set of tests
      local_idx   = match(res$keep, names(p)))
  }
  list(mod = mod, p = unname(p), n_terms = length(p), rows = rows)
}

t0 <- Sys.time()
idx <- seq_len(nrow(testable))
if (CORES > 1) {
  say("testing ", fmt_n(length(idx)), " modules on ", CORES, " cores ...")
  out <- parallel::mclapply(idx, test_one, mc.cores = CORES, mc.preschedule = FALSE)
} else {
  say("testing ", fmt_n(length(idx)), " modules ...")
  out <- vector("list", length(idx))
  for (i in idx) {
    out[[i]] <- test_one(i)
    if (i %% 50 == 0)
      say(sprintf("  %s/%s  (%.1f min)", fmt_n(i), fmt_n(length(idx)),
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}
failed <- vapply(out, is.null, logical(1))
if (any(failed)) say("modules that errored: ", sum(failed))
out <- out[!failed]
say(sprintf("scored %s modules in %.1f min", fmt_n(length(out)),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# --- cross-module BH ---------------------------------------------------------
# Per-module testing introduces a second burden 09 never faced: ~650 modules x
# ~thousands of terms. p.adj_global is BH over ALL of those tests at once. It is
# reported, NOT used to select -- weight01 p-values are not an exchangeable
# family, and pooling them across modules only compounds that -- but a term that
# survives it is on much firmer ground than one that only clears the raw p.
#
# The full p-vector of every module is kept for this, not just the retained
# terms: BH takes a cumulative minimum from the largest p downwards, so the
# non-significant tail genuinely changes the adjusted values of the head.
lens    <- vapply(out, function(o) o$n_terms, integer(1))
offsets <- c(0L, head(cumsum(lens), -1L))
all_p   <- unlist(lapply(out, function(o) o$p), use.names = FALSE)
all_padj <- p.adjust(all_p, method = "BH")
say("cross-module BH over ", fmt_n(length(all_p)), " module x term tests")

rows <- rbindlist(lapply(seq_along(out), function(j) {
  r <- out[[j]]$rows
  if (is.null(r)) return(NULL)
  r[, p.adj_global := signif(all_padj[offsets[j] + local_idx], 4)]
  r[, local_idx := NULL]
  r
}))

# --- plotting ----------------------------------------------------------------
# Both figures are dot plots in 09_go_enrichment.r's idiom, so a module figure and
# the conserved-set figure can be read side by side without relearning the axes.
# The p = 0 flooring is carried over from there too: -log10(0) is infinite and
# silently drops the most significant term off the panel.
floor_zero <- function(p) {
  z <- p == 0
  if (any(z)) {
    nz <- p[!z]
    p[z] <- if (length(nz)) min(nz) / 2 else GO_P / 1000
  }
  p
}

# Counts of genes and of modules are integers; pretty() happily returns 1.5 and a
# legend reading "1.5 genes" is nonsense.
int_breaks <- function(x) { b <- unique(round(pretty(x))); b[b >= 1] }

# Long subtitles are silently clipped at the device edge, and they carry the
# denominators, so they get wrapped rather than trimmed.
wrap_sub <- function(x, width = 78) paste(strwrap(x, width), collapse = "\n")

save_plot <- function(gg, base, w, h) {
  suppressMessages(ggsave(paste0(base, ".png"), gg, device = "png",
                          width = w, height = h, units = "cm", dpi = 300, limitsize = FALSE))
  suppressMessages(ggsave(paste0(base, ".pdf"), gg, device = "pdf",
                          width = w, height = h, units = "cm", limitsize = FALSE))
}

# GRAIN: one module. -log10(p) per term, sized by how many of the module's genes
# carry it, and filled by whether the term also survives the cross-module BH --
# that fill is the honest part of the figure, because most terms do not.
plot_module <- function(d, mod, base) {
  d <- head(d[order(pvalue)], NTOP)
  d[, pv := floor_zero(pvalue)]
  d[, Term := factor(Term, levels = rev(unique(Term)))]
  d[, robust := fifelse(p.adj_global <= 0.05, "survives cross-module BH", "raw p only")]
  r <- testable[module == mod]
  gg <- ggplot(d, aes(x = Term, y = -log10(pv), size = Significant, colour = robust)) +
    geom_point() +
    scale_size(range = c(2.5, 9), breaks = int_breaks) +
    scale_colour_manual(values = c(`survives cross-module BH` = "#B2182B",
                                   `raw p only` = "grey55"), drop = FALSE) +
    coord_flip() +
    labs(title = sprintf("GO %s - %s (%s)", ONTOLOGY, mod, STUDY),
         subtitle = wrap_sub(sprintf(
           "%s, %s genes (%s GO-annotated), PC1 %.0f%% var  |  %d terms at raw p <= %.2g",
           r$direction, fmt_n(r$n_genes), fmt_n(r$n_annotated),
           r$pc1_var_pct, nrow(d), GO_P)),
         x = NULL, y = expression(-log[10](p)),
         size = "genes in\nthe module", colour = NULL) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          plot.title.position = "plot",
          legend.position = "bottom", legend.box = "horizontal",
          axis.text.y = element_text(size = 9))
  # Floor is per-row, not a fixed minimum: MF and CC often return one or two terms
  # and a 9 cm device stretches a single point across an empty panel.
  save_plot(gg, base, 24, max(6.5, 4 + 0.65 * nrow(d)))
}

# GLOBAL: all modules at once. Terms are RANKED by how many distinct modules they
# are enriched in -- a term found once is a lead, a term found in six independent
# modules is a pathway -- and recurrence gets the size channel, reinforced by
# colour. It is deliberately NOT the x axis: when no term recurs, as in purple,
# that axis collapses onto a single value and the figure says nothing. -log10(p)
# on x is informative in both cases, so one layout serves both studies.
plot_global <- function(rows, base) {
  top <- rows[, .(modules = uniqueN(module), best_p = min(pvalue),
                  genes = sum(Significant)), by = .(GO.ID, Term)]
  setorder(top, -modules, best_p)
  n_terms <- nrow(top)
  top <- head(top, NTOP)
  top[, bp := floor_zero(best_p)]
  top[, Term := factor(Term, levels = rev(unique(Term)))]
  recurs <- top[modules > 1, .N]
  # When nothing recurs -- purple, where 3 modules share no term -- the module
  # channel is constant. A legend showing a single value is noise, so it goes and
  # the points fall back to one size and colour; the subtitle still states it.
  gg <- if (recurs == 0)
    ggplot(top, aes(x = Term, y = -log10(bp))) +
      geom_point(size = 4, colour = "#7B3294")
  else
    ggplot(top, aes(x = Term, y = -log10(bp), size = modules, colour = modules)) +
      geom_point() +
      scale_size(range = c(2.5, 10), breaks = int_breaks) +
      scale_colour_viridis_c(option = "magma", end = 0.85, direction = -1,
                             breaks = int_breaks) +
      guides(colour = guide_colourbar(order = 1), size = guide_legend(order = 2))
  gg <- gg +
    coord_flip() +
    labs(title = sprintf("GO %s across nitrogen-responsive modules - %s", ONTOLOGY, STUDY),
         subtitle = wrap_sub(sprintf(
           "%s of %s responsive modules testable (>= %d GO-annotated genes); %s of those have >= 1 term at raw p <= %.2g. Top %d of %s terms, ranked by how many modules share them (%d shared by more than one).",
           fmt_n(nrow(testable)), fmt_n(nrow(sel)), MIN_ANN,
           fmt_n(uniqueN(rows$module)), GO_P, nrow(top), fmt_n(n_terms), recurs), 95),
         x = NULL, y = expression(-log[10](best~p)),
         size = "modules", colour = "modules") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          plot.title.position = "plot",
          axis.text.y = element_text(size = 9))
  save_plot(gg, base, 28, max(7, 4 + 0.65 * nrow(top)))
}

# --- write -------------------------------------------------------------------
info <- testable[, .(module, direction, n_genes, pc1_var_pct, n_annotated)]
tag  <- file.path(OUT_DIR, sprintf("module_GO_%s_%s", ONTOLOGY, STUDY))

if (nrow(rows)) {
  rows <- merge(rows, info, by = "module", all.x = TRUE)
  setcolorder(rows, c("module", "direction", "n_genes", "n_annotated", "pc1_var_pct",
                      "GO.ID", "Term", "Annotated", "Significant", "Expected",
                      "pvalue", "p.adj", "p.adj_global"))
  setorder(rows, pvalue)
  write_tsv(rows, paste0(tag, ".tsv"))
} else {
  say("NOTE: no term cleared p <= ", GO_P, " in any module")
  write_tsv(data.table(module = character(), direction = character(),
                       n_genes = integer(), n_annotated = integer(),
                       pc1_var_pct = numeric(),
                       GO.ID = character(), Term = character(),
                       Annotated = integer(), Significant = integer(),
                       Expected = numeric(), pvalue = numeric(),
                       p.adj = numeric(), p.adj_global = numeric()),
            paste0(tag, ".tsv"))
}

# One row per module, INCLUDING the ones that returned nothing and the ones the
# annotation gate excluded. 09 drops skipped networks from its summary; here the
# denominator is the result -- "8 of 647 modules have any enriched term" is only
# readable if all 647 are in the file.
best <- if (nrow(rows)) rows[order(pvalue), .SD[1L], by = module,
                             .SDcols = c("GO.ID", "Term", "pvalue")] else
        data.table(module = character(), GO.ID = character(),
                   Term = character(), pvalue = numeric())
setnames(best, c("GO.ID", "Term", "pvalue"), c("top_GO", "top_Term", "top_pvalue"))
nsig <- if (nrow(rows)) rows[, .(n_sig_terms = .N), by = module] else
        data.table(module = character(), n_sig_terms = integer())

summ <- sel[, .(module, direction, n_genes, pc1_var_pct, n_annotated)]
summ[, tested := module %chin% testable$module]
summ <- merge(summ, nsig, by = "module", all.x = TRUE)
summ <- merge(summ, best, by = "module", all.x = TRUE)
summ[is.na(n_sig_terms), n_sig_terms := 0L]
summ[tested == FALSE, n_sig_terms := NA_integer_]
setorder(summ, -n_sig_terms, -n_annotated, na.last = TRUE)
write_tsv(summ, paste0(tag, "_summary.tsv"))

# --- per-module directories --------------------------------------------------
# One directory per TESTED module, table always, figure only when there is
# something to draw. Every tested module gets a directory even when it returned
# nothing, so the presence of a directory means "this was tested" and the absence
# of a figure inside it means "nothing cleared the threshold" -- neither has to be
# looked up in the summary table.
mod_root <- ensure_dir(file.path(OUT_DIR, "modules"))
n_dirs <- 0L; n_figs <- 0L
for (mod in testable$module) {
  d <- if (nrow(rows)) rows[module == mod] else rows
  md <- ensure_dir(file.path(mod_root, mod))
  base <- file.path(md, sprintf("GO_%s_%s_%s", ONTOLOGY, mod, STUDY))
  fwrite(d, paste0(base, ".tsv"), sep = "\t", quote = FALSE)
  n_dirs <- n_dirs + 1L
  if (nrow(d)) { plot_module(copy(d), mod, base); n_figs <- n_figs + 1L }
}
say("wrote ", fmt_n(n_dirs), " per-module directories under ", basename(mod_root),
    "/  (", fmt_n(n_figs), " with a figure)")

# --- the global figure -------------------------------------------------------
if (nrow(rows)) {
  plot_global(rows, paste0(tag, "_global"))
  say("wrote ", basename(tag), "_global.{png,pdf}")
} else {
  say("no global figure — nothing enriched anywhere")
}

# --- readout -----------------------------------------------------------------
banner(paste0("per-module GO ", ONTOLOGY, " — ", STUDY))
n_hit <- summ[n_sig_terms > 0, .N]
say(sprintf("%s of %s tested modules have >= 1 enriched %s term (raw p <= %.2g)",
            fmt_n(n_hit), fmt_n(nrow(testable)), ONTOLOGY, GO_P))
say(sprintf("  %s of %s responsive modules cleared the annotation gate",
            fmt_n(nrow(testable)), fmt_n(nrow(sel))))
if (nrow(rows)) {
  say(sprintf("  terms written: %s  |  clearing cross-module BH 0.05: %s",
              fmt_n(nrow(rows)), fmt_n(rows[p.adj_global <= 0.05, .N])))
  say("")
  say("most recurrent terms across modules:")
  top <- rows[, .(modules = uniqueN(module), best_p = min(pvalue)), by = .(GO.ID, Term)]
  setorder(top, -modules, best_p)
  print(head(top, 15), row.names = FALSE)
  say("")
  say("modules by enriched-term count:")
  print(head(summ[n_sig_terms > 0,
                  .(module, direction, n_genes, n_annotated, n_sig_terms, top_Term)], 15),
        row.names = FALSE)
}
say("")
say("done: ", STUDY)
