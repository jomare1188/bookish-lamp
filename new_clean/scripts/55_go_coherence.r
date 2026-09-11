#!/usr/bin/env Rscript
# =============================================================================
# 55_go_coherence.r — does an annotation make co-expressed modules more
# functionally coherent? The judge that decides which GO sources get adopted.
#
# WHY THIS EXISTS. Coverage is the wrong target on its own, and we measured that
# directly: moving GO from 8% to 56% of genes made 4.6x more modules testable and
# yet produced FEWER significant terms (sugarcane 147 -> 143, purple 14 -> 8),
# because the test family grew from 229,140 to 996,268 pairs and the old 8%
# background had been biased toward well-studied genes. So "more GO" cannot be the
# criterion for adopting an annotation source.
#
# THE CRITERION USED INSTEAD is the one 37_cluster_homogeneity.r already applies to
# PFAM domains, and the maths here is copied from it unchanged: the mean pairwise
# Sorensen-Dice similarity between the GO sets of genes in a module, measured
# against a SIZE-MATCHED RANDOM null. It asks whether genes that are co-expressed
# share function more than chance -- which is the property an annotation has to
# have to be worth anything downstream.
#
# THE NULL IS NOT OPTIONAL. Dice similarity rises trivially as modules shrink and
# as the annotation gets denser, so raw H would reward any annotation that simply
# adds terms. Only the excess over a random partition of the same module sizes,
# drawn from the same annotated gene pool, is evidence.
#
# SCORES SEVERAL ANNOTATIONS ON ONE PARTITION, which is the point: the modules are
# held fixed and only the annotation varies, so the comparison is like-for-like.
#
#   CLEAN_GENE2GO_SET="current=/path/a.tsv eggnog_auto=/path/b.tsv ..."
#
# RUN: through run.sh  ->  ./run.sh gocoherence sugarcane
# =============================================================================
suppressMessages({ library(data.table); library(Matrix) })
source(file.path(dirname(sub("--file=", "",
       grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "lib", "common.R"))

STUDY      <- env_req("CLEAN_STUDY")
MEMBERSHIP <- env_req("CLEAN_MEMBERSHIP")
SETSPEC    <- env_req("CLEAN_GENE2GO_SET")
OUT_FILE   <- env_req("CLEAN_OUT_FILE")
MAX_GENES  <- as.integer(env_num("CLEAN_MAX_GENES", 2000))
N_PERM     <- as.integer(env_num("CLEAN_PERM", 100))
SEED       <- as.integer(env_num("CLEAN_SEED", 1188))
MIN_MOD    <- as.integer(env_num("CLEAN_MIN_MODULE", 3))
setDTthreads(as.integer(env_num("CLEAN_CORES", 32)))
set.seed(SEED)

banner(paste("GO coherence:", STUDY))

mem <- fread(MEMBERSHIP, select = c("gene", "module_name"))
mem <- mem[module_name != "Unassigned"]
setnames(mem, c("gene", "module"))
keep <- mem[, .N, by = module][N >= MIN_MOD, module]
mem <- mem[module %chin% keep]
say(fmt_n(uniqueN(mem$module)), " modules of >= ", MIN_MOD, " genes, ",
    fmt_n(nrow(mem)), " genes")

sets <- strsplit(trimws(SETSPEC), "[[:space:]]+")[[1]]
sets <- sets[nzchar(sets)]
res <- list(); per_mod <- list()

for (spec in sets) {
  nm   <- sub("=.*$", "", spec)
  path <- sub("^[^=]*=", "", spec)
  if (!file.exists(path)) { say("SKIP ", nm, ": no file at ", path); next }

  d <- fread(path, sep = "\t", header = TRUE, select = c("gene", "go_id"),
             colClasses = "character")
  d <- unique(d[grepl("^GO:", go_id)])
  if (!nrow(d)) { say("SKIP ", nm, ": no GO rows"); next }

  # gene x term incidence, sparse. Same construction as 37.
  gi <- data.table(gene = sort(unique(d$gene)))[, i := .I]
  ti <- data.table(go_id = sort(unique(d$go_id)))[, j := .I]
  d  <- merge(merge(d, gi, by = "gene"), ti, by = "go_id")
  B  <- sparseMatrix(i = d$i, j = d$j, x = 1, dims = c(nrow(gi), nrow(ti)))
  rownames(B) <- gi$gene
  SETSZ <- rowSums(B)
  ALL <- rownames(B)

  dice_mean <- function(gv) {
    gv <- gv[gv %chin% ALL]; m <- length(gv)
    if (m < 2L) return(c(NA_real_, m, 0))
    sub <- FALSE
    if (m > MAX_GENES) { gv <- sample(gv, MAX_GENES); m <- MAX_GENES; sub <- TRUE }
    Bm <- B[gv, , drop = FALSE]
    I  <- as.matrix(tcrossprod(Bm)); s <- SETSZ[gv]
    D  <- 2 * I / outer(s, s, "+")
    c(mean(D[upper.tri(D)]), m, as.numeric(sub))
  }
  # Size-matched null, memoised by module size: the expensive part is the m x m
  # block, and many modules share a size.
  cache <- new.env(hash = TRUE)
  null_for <- function(m) {
    vapply(m, function(k) {
      key <- as.character(k)
      if (!is.null(cache[[key]])) return(cache[[key]])
      v <- mean(vapply(seq_len(N_PERM),
                       function(i) dice_mean(sample(ALL, min(k, length(ALL))))[1],
                       numeric(1)), na.rm = TRUE)
      assign(key, v, envir = cache); v
    }, numeric(1))
  }

  s <- mem[gene %chin% ALL][, { r <- dice_mean(gene)
         .(n_annotated = as.integer(r[2]), subsampled = r[3] > 0, H = r[1]) }, by = module]
  ok <- s[!is.na(H)]
  if (!nrow(ok)) { say("SKIP ", nm, ": no scorable module"); next }
  ok[, H_null := null_for(n_annotated)]
  # Gene-weighted, so an annotation is not rewarded for scoring many tiny modules.
  obs  <- ok[, sum(H * n_annotated) / sum(n_annotated)]
  nul  <- ok[, sum(H_null * n_annotated) / sum(n_annotated)]
  covg <- uniqueN(intersect(mem$gene, ALL))

  say(sprintf("  %-18s genes %7s  modules scored %6s  H %.4f  null %.4f  excess %+.4f",
              nm, fmt_n(covg), fmt_n(nrow(ok)), obs, nul, obs - nul))
  res[[nm]] <- data.table(study = STUDY, annotation = nm,
                          genes_annotated_in_modules = covg,
                          terms = nrow(ti), modules_scored = nrow(ok),
                          n_subsampled = ok[subsampled == TRUE, .N],
                          H = obs, H_null = nul, H_excess = obs - nul)
  per_mod[[nm]] <- cbind(annotation = nm, ok)
}

if (!length(res)) stop("no annotation set could be scored", call. = FALSE)
out <- rbindlist(res)
setorder(out, -H_excess)
write_tsv(out, OUT_FILE)
write_tsv(rbindlist(per_mod), sub("\\.tsv$", "_per_module.tsv", OUT_FILE))

say("")
say("ranked by excess over a size-matched null (the adoption criterion):")
print(out[, .(annotation, genes = genes_annotated_in_modules, terms,
              modules = modules_scored, H = round(H, 4),
              null = round(H_null, 4), excess = round(H_excess, 4))], row.names = FALSE)
say("")
say("coverage and coherence can move in opposite directions; that is the point.")
say("done: ", STUDY)
