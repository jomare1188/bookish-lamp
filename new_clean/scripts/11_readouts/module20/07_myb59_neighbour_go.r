#!/usr/bin/env Rscript
# =============================================================================
# 07_myb59_neighbour_go.r -- what the AtMYB59 copies' NEIGHBOURHOODS are for.
#
# WHY THIS IS THE INTERESTING QUESTION NOW. Step 06 found that no sugarcane copy
# reproduces a meaningful part of the purple gene's neighbourhood: the best recovers
# 1.89% of its 635 neighbours, nothing survives BH across the nine copies, and the
# 1:1 ortholog (06Ag012300, OG0088735) is SPECIFICALLY indistinguishable from chance
# -- z = 0.53, p = 0.29, its 11 shared genes being exactly what a 2,540-gene
# neighbourhood yields against a degree-matched purple gene.
#
# So the gene-level answer is no. This asks the weaker and more interesting version:
# do the neighbourhoods do the same KIND of thing without sharing genes?
#
# One topGO weight01 Fisher run per neighbourhood per ontology. BP, MF and CC: MF
# cleared 272 terms against BP's 164 at module level, so on a domain-derived
# annotation it often says more about what a gene set DOES.
#
# THE BACKGROUND IS THAT NETWORK'S GO-ANNOTATED NODES -- the same universe 09, 18
# and 63 use (64,178 and 109,591) -- so these results are comparable with every other
# GO panel in the paper rather than being their own island.
#
# SELECTION IS ON THE RAW weight01 p, per 09_go_enrichment.r:177-181: weight01 makes
# a term's score depend on its neighbours', so its p-values are not an exchangeable
# family and BH's assumptions do not hold. A BH column is carried and does not select.
#
# POWER IS NOT EQUAL ACROSS SETS and the output must not hide it: neighbourhoods run
# from 76 to 2,198 genes. n_annotated is written per set, so a thin result is
# visibly thin rather than looking like an absence of function.
#
# > Base R throughout: data.table's auto-indexing SEGFAULTS in topGO_env
# > (forderv -> setkeyv -> setindexv). See 63_degree_go.r.
#
# RUN: ./run_all.sh 07
# =============================================================================

suppressMessages(library(topGO))

BASE <- "/dados04/jorge/comparative_saccharum"
RES  <- file.path(BASE, "new_clean/results")
OUT  <- file.path(RES, "readouts/module20")
NBR  <- file.path(OUT, "myb59_neighbourhoods")
ONTS <- c("BP", "MF", "CC")
NODESIZE <- 10
P_THR <- 0.05
EXPECT_BG <- c(sugarcane = 64178, purple = 109591)

say <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
fmt <- function(x) format(x, big.mark = ",", scientific = FALSE)

files <- list.files(NBR, pattern = "\\.txt$", full.names = TRUE)
if (!length(files)) stop("no neighbour lists -- run step 06 first", call. = FALSE)
say(length(files), " neighbourhoods to test x ", length(ONTS), " ontologies")

# --- per-species annotation and universe, built once -------------------------
G2 <- list(); UNIV <- list()
for (sp in c("sugarcane", "purple")) {
  g2 <- read.delim(file.path(BASE, "annotation", sp, sprintf("gene2go_%s.tsv", sp)),
                   colClasses = "character")
  g2 <- g2[grepl("^GO:", g2$go_id), ]
  G2[[sp]] <- split(g2$go_id, g2$gene)
  nm <- read.delim(file.path(RES, sp, sprintf("network_%s_node_metrics.tsv", sp)),
                   colClasses = "character")
  u <- intersect(nm$gene, names(G2[[sp]]))
  # The background must be the one every other GO stage tests against, or these
  # results cannot be set beside them.
  if (length(u) != EXPECT_BG[[sp]])
    stop(sp, " background is ", fmt(length(u)), " genes but 09/18/63 use ",
         fmt(EXPECT_BG[[sp]]), call. = FALSE)
  UNIV[[sp]] <- u
  say("  ", sp, " background: ", fmt(length(u)), " GO-annotated network nodes")
}

rows <- list()
for (f in files) {
  base <- sub("\\.txt$", "", basename(f))
  sp   <- sub("__.*$", "", base)
  gene <- sub("^[^_]*__", "", base)
  nb   <- readLines(f, warn = FALSE)
  nb   <- nb[nzchar(nb)]
  ann  <- intersect(nb, UNIV[[sp]])
  say("")
  say(gene, " (", sp, "): ", fmt(length(nb)), " neighbours, ",
      fmt(length(ann)), " GO-annotated")
  if (length(ann) < NODESIZE) {
    say("  too few annotated neighbours to test")
    next
  }
  sel <- factor(as.integer(UNIV[[sp]] %in% ann), levels = c(0, 1))
  names(sel) <- UNIV[[sp]]
  for (ont in ONTS) {
    GOd <- new("topGOdata", ontology = ont, allGenes = sel,
               annot = annFUN.gene2GO, gene2GO = G2[[sp]], nodeSize = NODESIZE)
    r <- runTest(GOd, algorithm = "weight01", statistic = "fisher")
    n_tested <- length(score(r))
    tb <- GenTable(GOd, pvalue = r, topNodes = n_tested, numChar = 200)
    tb$pvalue <- as.numeric(sub("^\\s*<\\s*", "", tb$pvalue))
    tb$p.adj <- signif(p.adjust(tb$pvalue, method = "BH"), 4)
    keep <- tb[tb$pvalue <= P_THR, , drop = FALSE]
    say(sprintf("    %s: %s terms tested, %s at raw p <= %.2g",
                ont, fmt(n_tested), fmt(nrow(keep)), P_THR))
    if (nrow(keep))
      rows[[length(rows) + 1]] <- data.frame(
        species = sp, gene = gene, ontology = ont,
        n_neighbours = length(nb), n_annotated = length(ann),
        GO.ID = keep$GO.ID, Term = keep$Term,
        Annotated = keep$Annotated, Significant = keep$Significant,
        Expected = keep$Expected, pvalue = keep$pvalue, p.adj = keep$p.adj,
        stringsAsFactors = FALSE)
  }
}

res <- do.call(rbind, rows)
f_out <- file.path(OUT, "myb59_neighbour_go.tsv")
write.table(res, f_out, sep = "\t", quote = FALSE, row.names = FALSE)
say("")
say("wrote ", basename(f_out), "  (", fmt(nrow(res)), " rows)")

# --- do the neighbourhoods converge on FUNCTION where they do not on genes? ---
# Step 06 says the gene sets barely overlap. If the term sets do, that is functional
# convergence without gene-level conservation; if they do not, the negative is
# complete. Either way it is measured here rather than asserted.
pu <- unique(res$gene[res$species == "purple"])
conv <- list()
for (ont in ONTS) {
  for (pg in pu) {
    pt <- res$GO.ID[res$species == "purple" & res$gene == pg & res$ontology == ont]
    for (sg in unique(res$gene[res$species == "sugarcane"])) {
      st <- res$GO.ID[res$species == "sugarcane" & res$gene == sg & res$ontology == ont]
      sh <- length(intersect(pt, st))
      un <- length(union(pt, st))
      conv[[length(conv) + 1]] <- data.frame(
        ontology = ont, purple_gene = pg, sugarcane_gene = sg,
        n_purple_terms = length(pt), n_sugarcane_terms = length(st),
        shared_terms = sh, jaccard = if (un) round(sh / un, 4) else 0,
        stringsAsFactors = FALSE)
    }
  }
}
cv <- do.call(rbind, conv)
write.table(cv, file.path(OUT, "myb59_neighbour_go_convergence.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n", strrep("=", 78), "\n",
    "TERM-LEVEL CONVERGENCE of neighbourhoods (they barely share GENES)\n",
    strrep("=", 78), "\n", sep = "")
for (ont in ONTS) {
  d <- cv[cv$ontology == ont & cv$purple_gene == "Soffic.09G0001580-9H", ]
  if (!nrow(d)) next
  d <- d[order(-d$shared_terms), ]
  cat(sprintf("\n%s  (purple focus gene has %d enriched terms)\n",
              ont, d$n_purple_terms[1]))
  for (i in seq_len(min(5, nrow(d))))
    cat(sprintf("   %-30s %4d terms, %3d shared, Jaccard %.3f\n",
                d$sugarcane_gene[i], d$n_sugarcane_terms[i],
                d$shared_terms[i], d$jaccard[i]))
}
say("")
say("done")
