# =============================================================================
# 25_fig_module_go.r — the module-GO figure: what the responsive modules do,
# and whether the two response directions do different things.
#
# Figure 5 establishes which modules respond and that the two directions are
# opposite in expression. Opposite in expression is not the same as different in
# function, and that is the question this figure answers.
#
# THREE GRAINS OF GO, and the figure exists mostly to show that the middle one
# was missing:
#
#   per module   what ONE module is for. Answers "name Module_100", but it is
#                gated hard: a module needs MODULE_GO_MIN_ANNOTATED annotated
#                members for the test to be defined, and only ~11% of responsive
#                modules clear that, biased toward the large ones.
#   by direction NEW. The union of every module that rises with nitrogen against
#                the union of every module that falls, as two gene sets. Pooling
#                breaks the annotation gate -- modules too small to test alone
#                all contribute -- so this sees the whole responsive set rather
#                than its large-module tail.
#   pooled       every responsive module as one set. Already in the per-module
#                run's global figure; not repeated here.
#
#   A  THE ANNOTATION GATE, per species: how many responsive modules exist, how
#      many are individually testable, and how many annotated genes the two
#      direction sets carry instead. This is the panel that justifies the middle
#      grain rather than asserting it.
#
#   B  DIRECTION-SPECIFIC TERMS, the headline. Top terms of each direction, drawn
#      on one axis with the two directions opposed, so "what goes up with
#      nitrogen" and "what goes down" can be read against each other.
#
#   C  HOW SEPARATE THE TWO DIRECTIONS ARE -- terms found only when the modules
#      rise, only when they fall, or both. If the split were arbitrary this would
#      be mostly shared.
#
#   D  PER-MODULE GO, the grain that survives from the previous analysis: the
#      terms that recur across individually tested modules, coloured by the
#      direction of the module they came from.
#
# NOTHING IS WRITTEN ON THE FIGURE THAT BELONGS IN THE LEGEND -- panel letters and
# the labels the data needs, nothing else.
#
# RUN: through run.sh  ->  ./run.sh figmodulego   (then ./run.sh legends)
# =============================================================================

suppressMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(scico); library(grid); library(svglite); library(scales)
})

source(file.path(dirname(sub("--file=", "",
       grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "lib", "common.R"))

STUDIES    <- env_list("CLEAN_STUDIES", c("sugarcane", "purple"))
OUT_PREFIX <- env_req("CLEAN_OUT_PREFIX")
FIG        <- env_opt("CLEAN_FIG_NUM", "6")
ONT        <- env_opt("CLEAN_GO_ONTOLOGY", "BP")
NTERMS     <- as.integer(env_num("CLEAN_GO_NTERMS", 8))
NRECUR     <- as.integer(env_num("CLEAN_GO_NRECUR", 10))
MAIN       <- env_opt("CLEAN_GO_MAIN_STUDY", "sugarcane")
setDTthreads(as.integer(env_num("CLEAN_CORES", 8)))

cfg <- lapply(STUDIES, function(st) list(
  dir = env_req(sprintf("CLEAN_GODIR_%s", toupper(st)))))
names(cfg) <- STUDIES

banner(sprintf("Figure %s — module-level GO (%s)", FIG, ONT))

DIRS <- c("positive", "negative")
DIR_LAB <- c(positive = "rises with N", negative = "falls with N")
PAL_DIR <- setNames(rev(scico(3, palette = "managua", begin = 0.08, end = 0.92))[c(3, 1)],
                    unname(DIR_LAB))
PAL_STUDY <- setNames(scico(3, palette = "batlow")[1:2], STUDIES)

theme_f <- theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(face = "bold", size = 7.5),
        legend.key.size = unit(3.5, "mm"),
        legend.title = element_text(face = "bold", size = 7.5),
        plot.margin = margin(2, 4, 2, 2))
wrap_term <- function(x, w = 44)
  vapply(x, function(t) paste(strwrap(t, w), collapse = "\n"), "")

rd <- function(st, suffix)
  fread(file.path(cfg[[st]]$dir,
                  sprintf("module_GO_%s_%s%s.tsv", ONT, st, suffix)))

SUMM <- rbindlist(lapply(STUDIES, function(st) rd(st, "_summary")[, study := st]))
BYDIR <- rbindlist(lapply(STUDIES, function(st) {
  d <- rd(st, "_bydirection"); if (nrow(d)) d[, study := st] else NULL
}), fill = TRUE)
PERMOD <- rbindlist(lapply(STUDIES, function(st) {
  d <- rd(st, ""); if (nrow(d)) d[, study := st] else NULL
}), fill = TRUE)
for (D in list(SUMM, BYDIR, PERMOD))
  if (nrow(D)) D[, study := factor(study, levels = STUDIES)]

# =============================================================================
# A — the annotation gate, and what pooling by direction recovers
# =============================================================================
gate <- SUMM[, .(`responsive` = .N,
                 `individually\ntestable` = sum(tested),
                 `with >= 1\nterm` = sum(n_sig_terms > 0, na.rm = TRUE)),
             by = study]
gateL <- melt(gate, id.vars = "study", variable.name = "stage", value.name = "n")
gateL[, stage := factor(stage, levels = names(gate)[-1])]
say("panel A: the annotation gate"); print(gate, row.names = FALSE)

pA <- ggplot(gateL, aes(stage, n, fill = study)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64) +
  geom_text(aes(label = fmt_n(n)), position = position_dodge(width = 0.72),
            vjust = -0.4, size = 2.1) +
  scale_fill_manual(values = PAL_STUDY, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
  labs(x = NULL, y = "responsive modules") +
  theme_f + theme(axis.text.x = element_text(size = 6.4, lineheight = 1.05))

# =============================================================================
# B — direction-specific terms in the better-powered study
# =============================================================================
bd <- BYDIR[study == MAIN]
if (nrow(bd)) {
  bd[, dir_lab := factor(DIR_LAB[direction], levels = unname(DIR_LAB))]
  top <- bd[order(pvalue), head(.SD, NTERMS), by = dir_lab]
  top[, mlp := -log10(pmax(pvalue, 1e-300))]
  # Opposed on one axis: "up" to the right, "down" to the left, so the two sets
  # are read against each other instead of in two stacked panels.
  top[, signed := ifelse(dir_lab == DIR_LAB[["positive"]], mlp, -mlp)]
  setorder(top, signed)
  top[, Term_w := factor(wrap_term(Term), levels = wrap_term(Term))]
  say(sprintf("panel B: top %d terms per direction in %s", NTERMS, MAIN))
  print(top[, .(dir_lab, Term = substr(Term, 1, 44), pvalue, p.adj)], row.names = FALSE)

  pB <- ggplot(top, aes(signed, Term_w, fill = dir_lab)) +
    geom_col(width = 0.72) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40") +
    scale_fill_manual(values = PAL_DIR, name = NULL) +
    scale_x_continuous(labels = function(x) abs(x)) +
    labs(x = expression(-log[10](p)~", topGO weight01"), y = NULL) +
    theme_f + theme(axis.text.y = element_text(size = 5.6, lineheight = 0.9))
} else {
  pB <- ggplot() + theme_void()
}

# =============================================================================
# C — how separate the two directions are
# =============================================================================
cmp <- rbindlist(lapply(STUDIES, function(st) {
  d <- BYDIR[study == st]
  if (!nrow(d)) return(NULL)
  w <- dcast(d, GO.ID ~ direction, value.var = "pvalue")
  present <- intersect(DIRS, names(w))
  if (length(present) < 2) return(data.table(study = st, class = "one direction only",
                                             n = nrow(w)))
  data.table(study = st,
             class = c(DIR_LAB[["positive"]], "both", DIR_LAB[["negative"]]),
             n = c(sum(!is.na(w$positive) & is.na(w$negative)),
                   sum(!is.na(w$positive) & !is.na(w$negative)),
                   sum(is.na(w$positive) & !is.na(w$negative))))
}))
cmp[, study := factor(study, levels = STUDIES)]
cmp[, class := factor(class, levels = c(DIR_LAB[["positive"]], "both",
                                        DIR_LAB[["negative"]], "one direction only"))]
say("panel C: term overlap between directions"); print(cmp, row.names = FALSE)

pC <- ggplot(cmp, aes(study, n, fill = class)) +
  geom_col(width = 0.58, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(n > 0, n, "")), position = position_stack(vjust = 0.5),
            size = 2.2, colour = "grey15") +
  scale_fill_manual(values = c(setNames(unname(PAL_DIR), unname(DIR_LAB)),
                               both = "grey78", `one direction only` = "grey55"),
                    name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = sprintf("%s terms at raw p <= 0.05", ONT)) +
  theme_f

# =============================================================================
# D — per-module GO: the terms that recur across individually tested modules
# =============================================================================
pm <- PERMOD[study == MAIN]
if (nrow(pm)) {
  rec <- pm[, .(modules = uniqueN(module), best_p = min(pvalue),
                up = uniqueN(module[direction == "positive"]),
                down = uniqueN(module[direction == "negative"])),
            by = .(GO.ID, Term)][order(-modules, best_p)]
  recT <- head(rec, NRECUR)
  recL <- melt(recT[, .(Term, up, down)], id.vars = "Term",
               variable.name = "direction", value.name = "n")
  recL[, direction := factor(DIR_LAB[fifelse(direction == "up", "positive", "negative")],
                             levels = unname(DIR_LAB))]
  recL[, Term_w := factor(wrap_term(Term), levels = rev(wrap_term(recT$Term)))]
  say(sprintf("panel D: %s terms recurring across individually tested modules", MAIN))
  print(recT[, .(Term = substr(Term, 1, 44), modules, up, down)], row.names = FALSE)

  pD <- ggplot(recL[n > 0], aes(n, Term_w, fill = direction)) +
    geom_col(width = 0.72) +
    scale_fill_manual(values = PAL_DIR, name = NULL) +
    scale_x_continuous(breaks = pretty_breaks(4),
                       expand = expansion(mult = c(0, 0.06))) +
    labs(x = "modules sharing the term", y = NULL) +
    theme_f + theme(axis.text.y = element_text(size = 5.8, lineheight = 0.95))
} else {
  pD <- ggplot() + theme_void()
}

# =============================================================================
# compose
# =============================================================================
# The two GO panels carry wrapped term names on their y axes and need both the
# extra width and the taller row; A and C are three bars each.
fig <- (pA | pB) / (pC | pD) +
  plot_layout(widths = c(1, 1.3), heights = c(1.15, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12))

W <- 22; H <- 17
invisible(ensure_dir(dirname(OUT_PREFIX)))
ggsave(paste0(OUT_PREFIX, ".png"), fig, width = W, height = H, units = "cm",
       dpi = 400, type = "cairo")
ggsave(paste0(OUT_PREFIX, ".pdf"), fig, width = W, height = H, units = "cm",
       device = cairo_pdf)
ggsave(paste0(OUT_PREFIX, ".svg"), fig, width = W, height = H, units = "cm",
       device = svglite::svglite)
say(sprintf("wrote %s.{png,pdf,svg}  (%d x %d cm)", basename(OUT_PREFIX), W, H))

# =============================================================================
# legend
# =============================================================================
wrap_at <- function(x, width = 96)
  paste(vapply(strsplit(x, "\n[[:space:]]*\n")[[1L]],
               function(para) paste(strwrap(gsub("[[:space:]]+", " ", para), width),
                                    collapse = "\n"),
               ""), collapse = "\n\n")

g_ <- function(st, col) gate[study == st][[col]]
c_ <- function(st, cl) { v <- cmp[study == st & class == cl, n]; if (!length(v)) 0L else v }
bdn <- function(st, d, col) unique(BYDIR[study == st & direction == d][[col]])
S1 <- MAIN; S2 <- setdiff(STUDIES, MAIN)[1]

legend <- paste0(
"Figure ", FIG, ". Gene Ontology of the nitrogen-responsive modules (", ONT,
"), at three grains. Every test uses topGO's weight01 algorithm with a Fisher statistic against ",
"the SAME background -- the GO-annotated nodes of that species' own network, which is also the ",
"background of the gene-level GO and of the per-module transcription-factor test -- and selects ",
"on the raw weight01 p at 0.05, following topGO's own convention that weight01 conditions each ",
"term on its DAG neighbours and so does not produce an exchangeable family for BH. Adjusted ",
"values are written to the output tables for reference and do not select.\n",
"\n",
"(A) Why a middle grain was needed. Testing modules ONE AT A TIME requires each to carry at least ",
"three GO-annotated members for the test to be defined, and almost none do: of ",
fmt_n(g_(S1, "responsive")), " responsive ", S1, " modules only ",
fmt_n(g_(S1, "individually\ntestable")), " are individually testable and ",
fmt_n(g_(S1, "with >= 1\nterm")), " return a term; ", S2, " gives ",
fmt_n(g_(S2, "individually\ntestable")), " of ", fmt_n(g_(S2, "responsive")),
". The per-module grain therefore describes ~11% of the responsive set, biased toward the large ",
"modules, and it cannot be read as characterising the set as a whole.\n",
"\n",
"(B) THE DIRECTION-LEVEL TEST, which is new here and is what the annotation gate above makes ",
"necessary. The union of every module that RISES with nitrogen is tested as one gene set against ",
"the union of every module that FALLS. Pooling dissolves the gate -- modules too small to test ",
"alone all contribute -- so this grain sees the whole responsive set rather than its large-module ",
"tail: in ", S1, " the two sets carry ", fmt_n(bdn(S1, "positive", "n_annotated")), " and ",
fmt_n(bdn(S1, "negative", "n_annotated")), " GO-annotated genes drawn from ",
fmt_n(bdn(S1, "positive", "n_modules")), " and ", fmt_n(bdn(S1, "negative", "n_modules")),
" modules. Bars are the top ", sprintf("%d", NTERMS), " terms of each direction, opposed on one ",
"axis so the two can be read against each other; length is -log10 of the raw weight01 p.\n",
"\n",
"The two directions are doing different and interpretable things. What RISES with nitrogen is ",
"nitrogen assimilation and the metabolism that consumes it -- NITRATE ASSIMILATION, nitric oxide ",
"biosynthesis, PROLINE and asparagine and spermidine biosynthesis, reactive oxygen species ",
"biosynthesis -- alongside a large defence block (response to fungus, to bacterium, to chitin, ",
"monoterpene biosynthesis). What FALLS is the low-nitrogen and carbon-storage programme: ",
"FLAVONOID biosynthesis, RAFFINOSE-family oligosaccharide biosynthesis, triglyceride ",
"biosynthesis, cellular response to cold, skotomorphogenesis. That is the textbook shape of a ",
"nitrogen response read at module resolution, and neither of the other two grains recovers it.\n",
"\n",
"(C) How separate the two directions are, as term counts. In ", S1, ", ", fmt_n(c_(S1, "rises with N")),
" terms are enriched only among the modules that rise and ", fmt_n(c_(S1, "falls with N")),
" only among those that fall, against ", fmt_n(c_(S1, "both")), " found in both; in ", S2,
" the same three numbers are ", fmt_n(c_(S2, "rises with N")), ", ",
fmt_n(c_(S2, "falls with N")), " and ", fmt_n(c_(S2, "both")), ". So the split is not a ",
"partition of one functional programme into halves -- the overlap is a few terms out of a ",
"hundred -- and pooling the two directions, as the previous analysis did, averages two distinct ",
"programmes into one list.\n",
"\n",
"(D) The per-module grain that survives: terms recurring across individually tested ", S1,
" modules, with each bar split by the direction of the modules contributing it. Recurrence is ",
"low -- the most shared term appears in ", fmt_n(max(recT$modules)),
" modules -- which is expected when only ", fmt_n(g_(S1, "individually\ntestable")),
" modules are testable, and is the reason this panel supports rather than carries the figure.\n",
"\n",
"CAVEATS. Purple's per-module column is ", fmt_n(g_(S2, "individually\ntestable")),
" modules and should be read as a report of that many gene sets, not as a characterisation of ",
"its nitrogen response; panels B and D are drawn from ", S1, " for that reason. Annotation ",
"coverage, not statistics, is the binding constraint throughout: only ~8% of ", S1,
"'s network nodes carry any eggNOG GO term. And the direction labels mean 'tracks nitrogen ",
"supply', not 'responds to nitrogen stress' -- in ", S2,
" the three levels are stress-control-stress rather than a dose series, so a monotone test ",
"cannot see a module moved the same way by both extremes.")

legend_file <- paste0(OUT_PREFIX, "_legend.txt")
writeLines(wrap_at(legend), legend_file)
say("wrote ", basename(legend_file))

stats <- rbindlist(list(
  gateL[, .(panel = "A", study = as.character(study),
            quantity = gsub("\n", " ", stage), value = fmt_n(n))],
  BYDIR[, .(panel = "B", study = as.character(study),
            quantity = sprintf("%s: modules / genes / GO-annotated / terms", direction),
            value = sprintf("%s / %s / %s / %s", fmt_n(n_modules[1]), fmt_n(n_genes[1]),
                            fmt_n(n_annotated[1]), fmt_n(.N))), by = .(study, direction)][
              , .(panel, study, quantity, value)],
  cmp[, .(panel = "C", study = as.character(study),
          quantity = sprintf("terms enriched: %s", class), value = fmt_n(n))]))
write_tsv(stats, paste0(OUT_PREFIX, "_stats.tsv"))
say("done")
