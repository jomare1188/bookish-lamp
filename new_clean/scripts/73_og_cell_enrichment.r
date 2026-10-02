#!/usr/bin/env Rscript
# =============================================================================
# 73_og_cell_enrichment.r -- what the two orthogroup cells are FOR: GO and
# InterPro, on orthogroups rather than genes.
#
# THE CONTRAST THAT ISOLATES DIVERGENCE is divergent/separable against
# uniform/separable, not either against the genome. Both cells are multi-copy
# separable orthogroups, so testing one against all genes would mostly recover
# "is in a multi-copy family" and "is expressed". 09_go_enrichment.r makes the
# same argument for conserved/nonconserved: two sets that partition one universe
# are a real contrast; a set against a universe that contains it is not.
#
# THE BACKGROUND IS COPY-NUMBER MATCHED, AND MATCHED WITHIN THE ANNOTATED SET.
# 72 did the matching and wrote the assignments; this stage reads them rather than
# re-deriving, so the two stages cannot disagree about which orthogroups are in
# which set. Two reasons the matching is not optional:
#   copy number  uniform groups have MORE copies (3.55 vs 2.82, p = 6e-236),
#                because concordance is easier to prove with more members
#   annotation   GO coverage is 57.5% vs 67.5% and the gap SURVIVES copy-number
#                matching -- ~10 points inside every size band, and InterPro shows
#                the same gap at a higher baseline (79.1% vs 88.2%). It cannot be
#                matched away, so matching happens among annotated orthogroups and
#                the gap is reported as a finding by 72 instead.
#
# AN ORTHOGROUP'S TERMS ARE THE UNION OF ITS MEMBERS'. 72 built og2go_<study>.tsv
# for this, with its first column named `gene` though it holds an orthogroup id --
# which is what lets parse_gene2go() and annFUN.gene2GO take it unchanged. The
# union is the right rule here: a term held by one copy and not its sibling IS the
# copy divergence under study, so requiring a fraction of copies to agree would
# discard the signal.
#
# GO AND INTERPRO NEED OPPOSITE MULTIPLE-TESTING TREATMENTS, and this is the one
# place in this project where BH selects.
#   GO       weight01 conditions each term on its neighbours in the DAG, so the
#            p-values are deliberately NOT an exchangeable family and BH's
#            assumptions do not hold. Selection is on the RAW p, as in 09, 18, 63
#            and 07. p.adj is written and carried but does not select.
#   InterPro has no DAG. There is nothing to condition on, the ~11,800 accessions
#            ARE exchangeable, the test is a plain Fisher, and BH applies properly.
#            So for InterPro the selecting column is p.adj. That inversion is
#            deliberate; without this paragraph it reads as an oversight.
#
# BASE R FOR THE topGO PATH. data.table's auto-indexing segfaults in topGO_env --
# `d[gene %chin% names(gene2GO)]` dies in forderv -> setkeyv -> setindexv
# (63_degree_go.r:51). fread for I/O is safe and 18_module_go.r relies on it, but
# nothing here filters on a large character key with data.table.
#
# RUN: through run.sh -> ./run.sh ogcellgo <study> [BP|MF|CC]
# =============================================================================

suppressMessages({
  library(topGO)
  library(GO.db)
})

env_req <- function(k) {
  v <- Sys.getenv(k)
  if (!nzchar(v)) stop("required environment variable ", k, " is unset -- ",
                       "this stage is meant to be launched through run.sh", call. = FALSE)
  v
}
env_opt <- function(k, d = "") { v <- Sys.getenv(k); if (nzchar(v)) v else d }
fmt_n   <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
say     <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

STUDY    <- env_req("CLEAN_STUDY")
SETS     <- env_req("CLEAN_CELL_SETS")
OG2GO    <- env_req("CLEAN_OG2GO")
IPS_FILE <- env_req("CLEAN_IPS_TSV")
ORTHO    <- env_req("CLEAN_ORTHOGROUPS")
OG_COL   <- env_req("CLEAN_OG_SPECIES")
OUT_DIR  <- env_req("CLEAN_OUT_DIR")
ONT      <- env_opt("CLEAN_ONTOLOGY", "BP")
P_THR    <- as.numeric(env_opt("CLEAN_GO_P", "0.05"))
NODESIZE <- as.integer(env_opt("CLEAN_GO_NODESIZE", "10"))
MIN_ANN  <- as.integer(env_opt("CLEAN_OG_CELL_MIN_ANNOTATED", "10"))
# THE SAME MINIMUM topGO APPLIES TO GO TERMS, APPLIED TO IPR ACCESSIONS. Without
# it this stage tested 6,435 sugarcane accessions when only 421 occur in >= 10
# orthogroups: a domain present in one or two groups cannot reach significance,
# and 6,000 such tests inflate the BH denominator 15-fold. Measured 2026-10-02,
# that buried the transporter signal on the uniform side -- MFS transporter at
# p = 1.4e-04 came back p.adj = 0.174. nodeSize already does this for GO
# (CLEAN_GO_NODESIZE); leaving InterPro without it was an unjustified asymmetry.
IPR_MIN  <- as.integer(env_opt("CLEAN_IPR_MIN_COUNT", "10"))
MAX_COP  <- as.integer(env_opt("CLEAN_MAX_OG_COPIES", "20"))
# "" tests every duplication class; "homeolog" restricts to copies the recent
# polyploidy retained, which is the control that separates a polyploidy result
# from a transposon-family artefact. 72 pre-matched both, so this only selects a
# column -- the matching itself never happens here.
CLASS    <- env_opt("CLEAN_CELL_CLASS", "")
if (!CLASS %in% c("", "homeolog"))
  stop("CLEAN_CELL_CLASS must be empty or 'homeolog' (got '", CLASS, "')", call. = FALSE)
SUFFIX   <- if (nzchar(CLASS)) paste0("_", CLASS) else ""

if (!ONT %in% c("BP", "MF", "CC"))
  stop("CLEAN_ONTOLOGY must be BP, MF or CC (got '", ONT, "')", call. = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat(strrep("=", 70), "\n", STUDY, " cell enrichment: ", ONT, "\n",
    strrep("=", 70), "\n", sep = "")

# --- the sets, as 72 assigned them -------------------------------------------
s <- read.delim(SETS, colClasses = "character")
for (cl in c("og_2sp", "divergent", "matched_go", "matched_ipr", "n_members"))
  if (!cl %in% names(s))
    stop(basename(SETS), " needs a `", cl, "` column -- re-run ogcellcontrast",
         call. = FALSE)
s$divergent   <- s$divergent   == "TRUE"
gcol_m <- paste0("matched_go", SUFFIX); icol_m <- paste0("matched_ipr", SUFFIX)
for (cl in c(gcol_m, icol_m))
  if (!cl %in% names(s))
    stop(basename(SETS), " has no `", cl, "` column -- re-run ogcellcontrast",
         call. = FALSE)
s$matched_go  <- s[[gcol_m]] == "TRUE"
s$matched_ipr <- s[[icol_m]] == "TRUE"
s$n_members   <- as.integer(s$n_members)
say(sprintf("%s orthogroups in the two cells | %s matched for GO | %s for IPR%s",
            fmt_n(nrow(s)), fmt_n(sum(s$matched_go)), fmt_n(sum(s$matched_ipr)),
            if (nzchar(CLASS)) paste0("  [", CLASS, " class only]") else ""))

# --- the orthogroup -> GO map -------------------------------------------------
g2 <- read.delim(OG2GO, colClasses = "character")
if (!all(c("gene", "go_id") %in% names(g2)))
  stop(basename(OG2GO), " needs `gene` and `go_id` columns", call. = FALSE)
g2 <- g2[grepl("^GO:", g2$go_id), ]
og2GO <- split(g2$go_id, g2$gene)
say(sprintf("og2go: %s orthogroups annotated", fmt_n(length(og2GO))))

# The all-orthogroup background for the secondary contrast, CAPPED at MAX_COP
# copies. Uncapped it would include OG0000000's 1,233 purple copies, whose union
# of 1,233 genes' annotations would dominate every term count. 67 caps the
# classified cells at the same 20, so the two universes are then comparable.
oall <- read.delim(ORTHO, colClasses = "character", check.names = FALSE)
gcol <- grep(OG_COL, names(oall), fixed = TRUE)[1]
if (is.na(gcol)) stop("species column '", OG_COL, "' not in ", basename(ORTHO), call. = FALSE)
nmem <- vapply(strsplit(oall[[gcol]], ","), function(x) sum(nzchar(trimws(x))), integer(1))
og_all <- oall[[1]][nmem >= 1L & nmem <= MAX_COP]
say(sprintf("all-orthogroup background: %s groups with 1-%d expressed-capable copies",
            fmt_n(length(og_all)), MAX_COP))

# --- one topGOdata per universe, both directions off the same object ---------
# A control computed on a different background controls nothing: the two
# directions of a contrast are one contingency seen twice, so they must share the
# graph. updateGenes() is what makes that cheap (18_module_go.r:216).
mk <- function(universe, interesting) {
  gl <- factor(as.integer(universe %in% interesting), levels = c(0L, 1L))
  names(gl) <- universe
  gl
}

run_contrast <- function(label, fg, bg_extra, universe_name) {
  universe <- intersect(unique(c(fg, bg_extra)), names(og2GO))
  fg_a <- intersect(fg, universe)
  bg_a <- setdiff(universe, fg_a)
  if (length(fg_a) < MIN_ANN || length(bg_a) < MIN_ANN) {
    say(sprintf("  %-28s SKIP -- %s foreground / %s background annotated (min %d)",
                label, fmt_n(length(fg_a)), fmt_n(length(bg_a)), MIN_ANN))
    return(NULL)
  }
  GOd <- suppressMessages(new("topGOdata", ontology = ONT,
                              allGenes = mk(universe, fg_a),
                              annot = annFUN.gene2GO, gene2GO = og2GO[universe],
                              nodeSize = NODESIZE))
  rt <- suppressMessages(runTest(GOd, algorithm = "weight01", statistic = "fisher"))
  p  <- score(rt)                       # p directly: no GenTable string round-trip
  st <- suppressMessages(termStat(GOd, names(p)))
  tt <- suppressMessages(AnnotationDbi::Term(GO.db::GOTERM[names(p)]))
  out <- data.frame(
    study = STUDY, ontology = ONT, dup_class = if (nzchar(CLASS)) CLASS else "all",
    contrast = label, universe = universe_name,
    n_foreground = length(fg_a), n_background = length(bg_a),
    n_terms_tested = length(p),
    GO.ID = names(p), Term = unname(tt[names(p)]),
    Annotated = st$Annotated, Significant = st$Significant,
    Expected = round(st$Expected, 3),
    enrichment = round(st$Significant / pmax(st$Expected, 1e-9), 3),
    pvalue = unname(p), stringsAsFactors = FALSE)
  out$p.adj <- signif(p.adjust(out$pvalue, method = "BH"), 4)
  n_tested <- length(p)
  out <- out[out$pvalue <= P_THR, ]
  out <- out[order(out$pvalue), ]
  # THE CHANCE EXPECTATION, BESIDE THE COUNT. Selecting on the raw p is this
  # project's convention because weight01's p-values are not exchangeable -- but
  # that convention assumes the list is enriched for signal in the first place. If
  # the number clearing the threshold is at or below P_THR * n_tested, the raw-p
  # list is noise and only the BH survivors mean anything. Printing both is what
  # lets a reader tell the two situations apart.
  exp_chance <- P_THR * n_tested
  n_bh <- sum(out$p.adj <= 0.05)
  say(sprintf("  %-28s %s fg / %s bg | %s of %s terms at raw p <= %.2g (chance ~%.0f) | %s survive BH",
              label, fmt_n(length(fg_a)), fmt_n(length(bg_a)), fmt_n(nrow(out)),
              fmt_n(n_tested), P_THR, exp_chance, fmt_n(n_bh)))
  if (nrow(out) <= exp_chance)
    say("      ^ at or below chance: read the BH survivors, not the raw-p list")
  out
}

div_m <- s$og_2sp[s$divergent & s$matched_go]
uni_m <- s$og_2sp[!s$divergent & s$matched_go]
div_a <- s$og_2sp[s$divergent]
uni_a <- s$og_2sp[!s$divergent]

say("GO, weight01 + fisher, selecting on the raw p")
GO_ROWS <- list()
# 1 primary: the two cells against each other, copy-number matched
GO_ROWS$dm <- run_contrast("divergent vs uniform", div_m, uni_m, "matched cells")
GO_ROWS$um <- run_contrast("uniform vs divergent", uni_m, div_m, "matched cells")
# 2 context: each against every orthogroup. Reported as background, not a result:
#    it largely recovers multi-copy-ness, which is what both cells are.
GO_ROWS$da <- run_contrast("divergent vs all orthogroups", div_a, og_all, "all orthogroups")
GO_ROWS$ua <- run_contrast("uniform vs all orthogroups", uni_a, og_all, "all orthogroups")
go <- do.call(rbind, GO_ROWS[!vapply(GO_ROWS, is.null, TRUE)])

# --- the two directions are ONE contingency: check the sign, never infer it --
# A term enriched in the divergent run must be depleted in the uniform run on the
# same universe. Taking the direction from "which run produced it" would invert
# the biology wherever a term happens to clear threshold in both.
if (!is.null(GO_ROWS$dm) && !is.null(GO_ROWS$um)) {
  bothd <- intersect(GO_ROWS$dm$GO.ID, GO_ROWS$um$GO.ID)
  if (length(bothd)) {
    a <- GO_ROWS$dm[match(bothd, GO_ROWS$dm$GO.ID), ]
    b <- GO_ROWS$um[match(bothd, GO_ROWS$um$GO.ID), ]
    bad <- sum(sign(a$enrichment - 1) == sign(b$enrichment - 1))
    say(sprintf("  %s terms clear threshold in both directions; %s have the same sign",
                fmt_n(length(bothd)), fmt_n(bad)))
    if (bad > 0)
      say("  NOTE: same-sign terms in both directions mean the term is enriched ",
          "relative to different sub-universes; read `enrichment`, not the run name")
  }
}

gf <- file.path(OUT_DIR, sprintf("og_cell_GO_%s_%s%s.tsv", ONT, STUDY, SUFFIX))
write.table(go, gf, sep = "\t", quote = FALSE, row.names = FALSE)
say(sprintf("wrote %s (%s rows)", basename(gf), fmt_n(nrow(go))))

# --- InterPro: flat terms, Fisher, and BH genuinely applies -----------------
if (ONT == "BP") {           # IPR has no ontology; run it once, not three times
  say("InterPro, Fisher per accession, selecting on BH p.adj")
  cmd <- sprintf("awk -F'\\t' 'NF>=13 && $12!=\"-\" && $12!=\"\" {print $1\"\\t\"$12\"\\t\"$13}' %s | sort -u",
                 shQuote(IPS_FILE))
  ip <- read.delim(pipe(cmd), header = FALSE, colClasses = "character",
                   col.names = c("gene", "ipr", "desc"))
  # Map genes -> orthogroups without data.table, since this runs in topGO_env.
  gl <- strsplit(oall[[gcol]], ",")
  og_of <- rep(oall[[1]], lengths(gl))
  gn    <- trimws(unlist(gl, use.names = FALSE))
  keep  <- nzchar(gn)
  og_of <- og_of[keep]; gn <- gn[keep]
  g2og  <- og_of[match(ip$gene, gn)]
  ok    <- !is.na(g2og)
  pair  <- unique(data.frame(og = g2og[ok], ipr = ip$ipr[ok],
                             desc = ip$desc[ok], stringsAsFactors = FALSE))
  say(sprintf("  %s (orthogroup, IPR) pairs over %s accessions",
              fmt_n(nrow(pair)), fmt_n(length(unique(pair$ipr)))))

  fg <- s$og_2sp[s$divergent & s$matched_ipr]
  bg <- s$og_2sp[!s$divergent & s$matched_ipr]
  pf <- pair[pair$og %in% c(fg, bg), ]
  tab <- table(pf$ipr, pf$og %in% fg)
  if (ncol(tab) == 2L) {
    occ <- rowSums(tab)
    say(sprintf("  %s accessions present | %s occur in >= %d orthogroups and are tested",
                fmt_n(nrow(tab)), fmt_n(sum(occ >= IPR_MIN)), IPR_MIN))
    tab <- tab[occ >= IPR_MIN, , drop = FALSE]
    if (!nrow(tab)) { say("  no accession clears the minimum"); tab <- NULL }
  }
  if (!is.null(tab) && ncol(tab) == 2L) {
    nf <- length(fg); nb <- length(bg)
    res <- data.frame(
      study = STUDY, dup_class = if (nzchar(CLASS)) CLASS else "all",
      ipr = rownames(tab),
      desc = pair$desc[match(rownames(tab), pair$ipr)],
      n_divergent = as.integer(tab[, "TRUE"]), n_uniform = as.integer(tab[, "FALSE"]),
      stringsAsFactors = FALSE)
    res$pct_divergent <- round(100 * res$n_divergent / nf, 3)
    res$pct_uniform   <- round(100 * res$n_uniform / nb, 3)
    pv <- vapply(seq_len(nrow(res)), function(i) {
      m <- matrix(c(res$n_divergent[i], nf - res$n_divergent[i],
                    res$n_uniform[i],   nb - res$n_uniform[i]), 2)
      ft <- fisher.test(m)
      c(ft$p.value, unname(ft$estimate))
    }, numeric(2))
    res$odds_ratio <- round(pv[2, ], 4)
    res$pvalue <- pv[1, ]
    res$p.adj  <- signif(p.adjust(res$pvalue, method = "BH"), 4)
    res$direction <- ifelse(res$odds_ratio > 1, "divergent", "uniform")
    res <- res[order(res$p.adj, res$pvalue), ]
    ipf <- file.path(OUT_DIR, sprintf("og_cell_InterPro_%s%s.tsv", STUDY, SUFFIX))
    write.table(res, ipf, sep = "\t", quote = FALSE, row.names = FALSE)
    sig <- sum(res$p.adj <= 0.05)
    say(sprintf("  %s accessions tested | %s at BH p.adj <= 0.05 (%s divergent, %s uniform)",
                fmt_n(nrow(res)), fmt_n(sig),
                fmt_n(sum(res$p.adj <= 0.05 & res$direction == "divergent")),
                fmt_n(sum(res$p.adj <= 0.05 & res$direction == "uniform"))))
    say(sprintf("wrote %s", basename(ipf)))
  }
}
say("done")
