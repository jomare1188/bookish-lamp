#!/usr/bin/env Rscript
# =============================================================================
# 60_normalise_gene2go.r — reduce every gene's GO set to its MOST SPECIFIC terms.
#
# WHY THIS IS NOT OPTIONAL. Annotation sources disagree about whether they emit
# direct terms or terms already propagated up the GO DAG, and mixing the two
# breaks both things we want to do with them. Measured on sugarcane:
#
#                       terms/gene   carries GO:0008150 (biological_process)
#   InterPro + Pfam            2.3   0 genes
#   eggNOG (--go_evidence all) 69.9  75,696 genes
#
# eggNOG ships the full ancestor closure, including the three ontology ROOTS.
#
# It breaks the coherence judge: pairwise Dice similarity between gene GO sets is
# inflated for a propagated source because every pair shares the roots and most
# high-level ancestors, so that source would "win" on shared generalities rather
# than on shared function.
#
# It breaks topGO: topGO builds the ontology graph and propagates annotations
# ITSELF. Feeding it pre-propagated terms counts each gene again at every ancestor
# node, inflating the counts that the enrichment p-values are computed from.
#
# So every source is reduced here to the same convention -- keep a term only if no
# OTHER term assigned to that gene is a descendant of it -- and all downstream
# comparisons are then like-for-like. Ancestors come from GO.db, which ships with
# topGO, so no ontology file has to be fetched or version-matched separately.
#
# RUN: through run.sh  ->  ./run.sh gonorm <in.tsv> <out.tsv>
# =============================================================================
suppressMessages({ library(data.table); library(GO.db) })

args <- commandArgs(trailingOnly = TRUE)
IN  <- if (length(args) >= 1) args[1] else Sys.getenv("CLEAN_IN_FILE")
OUT <- if (length(args) >= 2) args[2] else Sys.getenv("CLEAN_OUT_FILE")
if (!nzchar(IN) || !nzchar(OUT)) stop("usage: 60_normalise_gene2go.r <in.tsv> <out.tsv>")
if (!file.exists(IN)) stop("no input at ", IN)

cat(sprintf("[%s] normalising %s\n", format(Sys.time(), "%H:%M:%S"), basename(IN)))
d <- fread(IN, sep = "\t", header = TRUE, colClasses = "character")
if (!all(c("gene", "go_id") %in% names(d))) stop("need gene and go_id columns")
d <- unique(d[grepl("^GO:", go_id)])
cat(sprintf("  in : %s pairs, %s genes, %.1f terms/gene\n",
            format(nrow(d), big.mark = ","), format(uniqueN(d$gene), big.mark = ","),
            nrow(d) / uniqueN(d$gene)))

# term -> its ancestors, across all three ontologies
anc <- c(as.list(GOBPANCESTOR), as.list(GOMFANCESTOR), as.list(GOCCANCESTOR))
roots <- c("GO:0008150", "GO:0003674", "GO:0005575")

used <- unique(d$go_id)
cat(sprintf("  %s distinct terms, %s of them known to GO.db\n",
            format(length(used), big.mark = ","),
            format(sum(used %in% names(anc)), big.mark = ",")))

# For each gene: drop any term that is an ancestor of another term it carries.
# Built once as a flat ancestor table, then anti-joined -- a per-gene loop over
# hundreds of thousands of genes would not finish.
A <- rbindlist(lapply(intersect(used, names(anc)), function(t) {
  a <- anc[[t]]; a <- a[a != "all"]
  if (!length(a)) NULL else data.table(go_id = t, ancestor = a)
}))
setkey(A, go_id)
da <- merge(d[, .(gene, go_id)], A, by = "go_id", allow.cartesian = TRUE)
drop <- unique(da[, .(gene, go_id = ancestor)])       # these are implied, not specific
out <- fsetdiff(unique(d[, .(gene, go_id)]), drop)
out <- out[!go_id %chin% roots]

# carry the extra columns back
extra <- setdiff(names(d), c("gene", "go_id"))
if (length(extra)) out <- merge(out, unique(d), by = c("gene", "go_id"), all.x = TRUE)

fwrite(out, OUT, sep = "\t", quote = FALSE)
cat(sprintf("  out: %s pairs, %s genes, %.1f terms/gene  (removed %s implied)\n",
            format(nrow(out), big.mark = ","), format(uniqueN(out$gene), big.mark = ","),
            nrow(out) / uniqueN(out$gene), format(nrow(d) - nrow(out), big.mark = ",")))
for (r in roots) {
  n <- sum(out$go_id == r)
  if (n > 0) stop("root ", r, " survived normalisation on ", n, " genes")
}
cat(sprintf("[%s] done\n", format(Sys.time(), "%H:%M:%S")))
