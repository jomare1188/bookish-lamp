#!/usr/bin/env Rscript
# =============================================================================
# 63_degree_go.r -- what HUBS are for, and what the PERIPHERY is for.
#
# Ranks every GO-annotated network node by degree and asks which GO terms sit
# systematically at one end of that ranking. Two readings of one ranking:
#
#     score = rank(-degree)   ->  terms concentrated among HUBS
#     score = rank(+degree)   ->  terms concentrated at the PERIPHERY
#
# WHY A KS TEST ON THE RANKING AND NOT A CUT. Degree spans four orders of magnitude
# -- sugarcane's 10th percentile is 2, its median 53, its 90th 6,677 -- so any
# decile boundary is arbitrary and throws away the middle 80% of the genes. The KS
# statistic uses the whole ranking, which is the same thing GSEA does.
#
# WHY topGO AND NOT fgsea. weight01 decorrelates the GO DAG: without it a parent and
# its children all score on the same genes and the top of the table fills with
# "translation", "cytoplasmic translation", "peptide biosynthetic process". fgsea
# 1.32.2 is installed but treats GO terms as independent sets, and it lives in a
# different conda env. weight01 + ks gives the threshold-free test AND the DAG
# handling every other GO panel in this paper uses. Measured: 26 s per direction.
#
# ---------------------------------------------------------------------------
# THE RESULT THIS EXISTS TO QUALIFY
#
# Degree and edge-conservation are badly confounded, and partly by arithmetic.
# Measured on the current graphs:
#
#     median degree, sugarcane   on a conserved edge  228   no conserved edge  18
#     median degree, purple      on a conserved edge 6916   no conserved edge 376
#     share on a conserved edge  top-degree decile 63.3%    bottom decile 7.1%
#
# "On AT LEAST ONE conserved edge" is nearly automatic for a gene with 6,677 edges
# at a 10.4% per-edge conservation rate, and unlikely for a gene with 2. So figure
# 4 panel B's conserved/non-conserved GO contrast is entangled with this one, and
# the two must be read together. That is why this runs on the SAME universe panel B
# uses -- network nodes with GO -- which is asserted below rather than assumed.
# ---------------------------------------------------------------------------
#
# MULTIPLE TESTING, following 09_go_enrichment.r:177-181 exactly: selection is on
# the RAW weight01 p. weight01's p-values are not an exchangeable family -- the
# algorithm deliberately makes a term's score depend on its neighbours' -- so BH's
# assumptions do not hold on them. A BH column is written anyway, over every tested
# term, because it is the number a reader will ask for, but it does not select.
#
# AN EFFECT SIZE IS WRITTEN BESIDE EVERY p. At n = 64,178 a KS test calls terms
# significant on small shifts, and a p-value does not say how much higher a term's
# genes sit. `degree_ratio` is the median degree of the term's annotated genes over
# the background median, which does.
#
# NOTE ON THE ENVIRONMENT: data.table's auto-indexing SEGFAULTS in topGO_env --
# `d[gene %chin% names(gene2GO)]` dies inside forderv -> setkeyv -> setindexv,
# reproducibly. Everything here is base R for that reason. Do not "modernise" it.
#
# RUN: through run.sh  ->  ./run.sh degreego sugarcane [BP|MF|CC]
# =============================================================================

suppressMessages(library(topGO))

env_req <- function(k) {
  v <- Sys.getenv(k)
  if (!nzchar(v)) stop("required environment variable ", k, " is unset -- ",
                       "this stage is meant to be launched through run.sh", call. = FALSE)
  v
}
env_opt <- function(k, d = "") { v <- Sys.getenv(k); if (nzchar(v)) v else d }
say <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
fmt_n <- function(x) format(x, big.mark = ",", scientific = FALSE)

STUDY    <- env_req("CLEAN_STUDY")
NODES    <- env_req("CLEAN_NODE_METRICS")
GENE2GO  <- env_req("CLEAN_GENE2GO")
OUT_FILE <- env_req("CLEAN_OUT_FILE")
ONT      <- env_opt("CLEAN_ONTOLOGY", "BP")
NODESIZE <- as.integer(env_opt("CLEAN_GO_NODESIZE", "10"))
P_THR    <- as.numeric(env_opt("CLEAN_GO_P", "0.05"))
EXPECT_N <- as.integer(env_opt("CLEAN_EXPECT_ANNOTATED", "0"))

if (!ONT %in% c("BP", "MF", "CC"))
  stop("CLEAN_ONTOLOGY must be BP, MF or CC (got '", ONT, "')", call. = FALSE)

cat("\n", strrep("=", 70), "\n",
    sprintf("degree-ranked GO (%s): %s   [weight01 + KS on the full ranking]", ONT, STUDY),
    "\n", strrep("=", 70), "\n", sep = "")

# --- inputs, base R throughout (see the environment note above) ---------------
nm <- read.delim(NODES, colClasses = "character")
if (!all(c("gene", "degree") %in% names(nm)))
  stop(basename(NODES), " needs `gene` and `degree` columns", call. = FALSE)
nm$degree <- as.numeric(nm$degree)

g2 <- read.delim(GENE2GO, colClasses = "character")
if (!all(c("gene", "go_id") %in% names(g2)))
  stop(basename(GENE2GO), " needs `gene` and `go_id` columns", call. = FALSE)
g2 <- g2[grepl("^GO:", g2$go_id), ]
gene2GO <- split(g2$go_id, g2$gene)

say(fmt_n(nrow(nm)), " network nodes; ", fmt_n(length(gene2GO)),
    " genes carry GO genome-wide")

keep <- nm$gene %in% names(gene2GO)
nm <- nm[keep, ]
say("universe = network nodes WITH GO: ", fmt_n(nrow(nm)))

# THE UNIVERSE MUST BE THE ONE FIGURE 4 PANEL B USES. This panel is a control on
# that one -- conserved-edge genes are hubs -- and a control computed on a different
# background controls nothing.
if (EXPECT_N > 0 && nrow(nm) != EXPECT_N)
  stop("universe is ", fmt_n(nrow(nm)), " genes but the conserved-set GO background ",
       "is ", fmt_n(EXPECT_N), ".\n  The two panels would not be comparable.",
       call. = FALSE)
if (EXPECT_N > 0) say("  matches the conserved-set GO background (", fmt_n(EXPECT_N), ")")

BG_MED <- median(nm$degree)
say("degree: min ", min(nm$degree), "  median ", BG_MED, "  max ", fmt_n(max(nm$degree)))

# The annotation-coverage gradient, measured and reported rather than left to be
# discovered: if the periphery were far less annotated its terms would be less well
# powered and the contrast partly an artefact of that.
allg <- read.delim(NODES, colClasses = "character")
allg$degree <- as.numeric(allg$degree)
allg <- allg[order(allg$degree), ]
k <- nrow(allg) %/% 10
cov_lo <- mean(allg$gene[seq_len(k)] %in% names(gene2GO))
cov_hi <- mean(allg$gene[seq.int(nrow(allg) - k + 1, nrow(allg))] %in% names(gene2GO))
say(sprintf("GO coverage: bottom degree decile %.1f%%, top %.1f%%  (gradient %.1f pts)",
            100 * cov_lo, 100 * cov_hi, 100 * (cov_hi - cov_lo)))

# --- one ranking, two directions ---------------------------------------------
# topGO's KS path treats SMALL scores as the interesting end, so the hub test ranks
# on -degree and the periphery test on +degree. Same genes, same annotation, only
# the sign differs -- which is what makes these two readings of one ranking rather
# than two analyses.
run_dir <- function(direction) {
  sgn <- if (direction == "hub") -1 else 1
  sc <- setNames(as.numeric(rank(sgn * nm$degree, ties.method = "average")), nm$gene)

  GOd <- new("topGOdata", ontology = ONT, allGenes = sc,
             # geneSelectionFun is unused by the KS statistic but topGOdata requires
             # one; the decile it names never selects anything here.
             geneSelectionFun = function(x) x <= quantile(x, 0.1),
             annot = annFUN.gene2GO, gene2GO = gene2GO, nodeSize = NODESIZE)

  t0 <- Sys.time()
  res <- runTest(GOd, algorithm = "weight01", statistic = "ks")
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  n_tested <- length(score(res))
  tab <- GenTable(GOd, pvalue = res, topNodes = n_tested, numChar = 200)
  # GenTable returns p as character, with "< 1e-30" for the floor.
  tab$pvalue <- as.numeric(sub("^\\s*<\\s*", "", tab$pvalue))
  tab$p.adj <- signif(p.adjust(tab$pvalue, method = "BH"), 4)

  # Effect size: how much higher (or lower) does this term's ranking actually sit?
  gl <- genesInTerm(GOd, tab$GO.ID)
  dmed <- vapply(tab$GO.ID, function(id) {
    g <- intersect(gl[[id]], nm$gene)
    if (!length(g)) NA_real_ else median(nm$degree[match(g, nm$gene)])
  }, numeric(1))

  out <- data.frame(
    study = STUDY, ontology = ONT, direction = direction,
    GO.ID = tab$GO.ID, Term = tab$Term, Annotated = tab$Annotated,
    pvalue = tab$pvalue, p.adj = tab$p.adj,
    median_degree_term = as.numeric(dmed),
    median_degree_background = BG_MED,
    degree_ratio = signif(as.numeric(dmed) / BG_MED, 4),
    stringsAsFactors = FALSE)
  out <- out[order(out$pvalue), ]
  say(sprintf("  %-9s %s terms tested, %s at raw p <= %.2g  (%.0f s)",
              direction, fmt_n(n_tested), fmt_n(sum(out$pvalue <= P_THR)), P_THR, secs))
  out
}

say("")
res <- do.call(rbind, lapply(c("hub", "periphery"), run_dir))

dir.create(dirname(OUT_FILE), showWarnings = FALSE, recursive = TRUE)
write.table(res, OUT_FILE, sep = "\t", quote = FALSE, row.names = FALSE)
say("")
say("wrote ", basename(OUT_FILE), "  (", fmt_n(nrow(res)), " rows)")

for (d in c("hub", "periphery")) {
  say("")
  say("top ", ONT, " terms at the ", toupper(d), " end:")
  h <- head(res[res$direction == d, ], 8)
  for (i in seq_len(nrow(h)))
    say(sprintf("  %-52s p %-9.3g  n %-6s  median degree %s (%.2fx bg)",
                substr(h$Term[i], 1, 52), h$pvalue[i], fmt_n(h$Annotated[i]),
                fmt_n(round(h$median_degree_term[i])), h$degree_ratio[i]))
}
say("")
say("Selection is on the RAW weight01 p: weight01 scores are not an exchangeable")
say("family, so BH's assumptions do not hold on them. p.adj is carried, not used.")
say("done: ", STUDY, " ", ONT)
