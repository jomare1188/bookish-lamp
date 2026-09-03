#!/usr/bin/env Rscript
# =============================================================================
# 37_cluster_homogeneity.r -- are the modules functionally coherent, or just small?
#
# clm info says how well a partition captures edge mass. It cannot say whether
# the modules mean anything biologically, and raising inflation until the giant
# module shatters would score well on area fraction while producing rubble.
# This is the independent criterion.
#
# THE METRIC is cogeqc's: per group, the mean pairwise Sorensen-Dice similarity
# between the genes' annotation sets -- here PFAM domains from eggNOG-mapper. A
# group whose genes all share a domain scores 1; a group with nothing in common
# scores 0.
#
# THE PACKAGE IS NOT USED, deliberately, and this is not a style preference:
#   * cogeqc::calculate_H has max_size = 200 and returns NOTHING for larger
#     groups. Purple's Module_001 has 47,887 genes. It would silently refuse to
#     score exactly the object being diagnosed, and subtract a constant from
#     everything else for the privilege.
#   * it enumerates pairs with combn() inside an R lapply. 47,887 genes is 1.1e9
#     pairs for one module of one clustering, and there are dozens of cells.
#   * it is installed in none of the three envs this pipeline uses.
# Here the same score is computed from a sparse binary gene x domain matrix --
# intersections are B %*% t(B), set sizes are rowSums -- and large modules are
# SUBSAMPLED rather than skipped, with the subsampling recorded per module.
# 38 checks this implementation against cogeqc on modules small enough for it.
#
# THE NULL IS NOT OPTIONAL. Homogeneity rises trivially as modules shrink, and
# shrinking modules is exactly what higher inflation does -- so raw H would
# reward the most fragmented clustering no matter what it contained. Every
# clustering is scored against a SIZE-MATCHED random partition: the same module
# size distribution, genes shuffled between them. The reported quantity is the
# gene-weighted mean of (H_observed - H_null), which is comparable across
# clusterings of different granularity in a way raw H is not.
#
# RUN: through run.sh  ->  ./run.sh clusterhomog sugarcane
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(Matrix)})
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE),
                                                 value = TRUE)[1])), "lib", "common.R"))

STUDY     <- env_req("CLEAN_STUDY")
WORK      <- env_req("CLEAN_WORK_DIR")
OUT_DIR   <- ensure_dir(env_req("CLEAN_OUT_DIR"))
EMAPPER   <- env_req("CLEAN_EMAPPER")
MAX_GENES <- as.integer(env_num("CLEAN_MAX_GENES", 2000))
N_PERM    <- as.integer(env_num("CLEAN_PERM", 100))
SEED      <- as.integer(env_num("CLEAN_SEED", 1188))
setDTthreads(as.integer(env_num("CLEAN_CORES", 32)))
set.seed(SEED)

banner(paste("annotation homogeneity of each clustering:", STUDY))

# --- annotations --------------------------------------------------------------
# Column 21 of emapper.annotations is PFAMs, comma-separated, "-" for none.
say("reading ", basename(EMAPPER))
ann <- fread(cmd = sprintf("grep -v '^##' %s", shQuote(EMAPPER)), sep = "\t",
             header = TRUE, quote = "")
setnames(ann, 1, "query")
if (!"PFAMs" %in% names(ann))
  stop("no PFAMs column in ", basename(EMAPPER), call. = FALSE)
a <- ann[PFAMs != "-" & !is.na(PFAMs) & nzchar(PFAMs), .(gene = strip_version(query), PFAMs)]
a <- a[, .(dom = trimws(strsplit(PFAMs, ",", fixed = TRUE)[[1]])), by = gene]
a <- unique(a[nzchar(dom) & dom != "-"])
say("  ", fmt_n(uniqueN(a$gene)), " genes carry ", fmt_n(uniqueN(a$dom)), " distinct PFAM domains")

genes_idx <- data.table(gene = sort(unique(a$gene)))[, gi := .I]
doms_idx  <- data.table(dom  = sort(unique(a$dom)))[,  di := .I]
a <- merge(merge(a, genes_idx, by = "gene"), doms_idx, by = "dom")
B <- sparseMatrix(i = a$gi, j = a$di, x = 1,
                  dims = c(nrow(genes_idx), nrow(doms_idx)))
rownames(B) <- genes_idx$gene
SETSZ <- rowSums(B)

# --- the score ----------------------------------------------------------------
# Mean over the strict upper triangle of 2*|Di ^ Dj| / (|Di| + |Dj|).
# sum(I) over all ordered pairs is sum(crossprod), and the diagonal is sum(|Di|),
# so the upper triangle is (total - diagonal)/2 -- no pair is ever enumerated.
# The denominator varies per pair, so the ratio cannot be reduced the same way and
# the m x m intersection block is formed; m is capped by MAX_GENES for that reason.
dice_mean <- function(gv) {
  gv <- gv[gv %chin% rownames(B)]
  m <- length(gv)
  if (m < 2L) return(c(NA_real_, m, 0))
  sub <- FALSE
  if (m > MAX_GENES) { gv <- sample(gv, MAX_GENES); m <- MAX_GENES; sub <- TRUE }
  Bm <- B[gv, , drop = FALSE]
  I  <- as.matrix(tcrossprod(Bm))            # m x m intersection counts
  s  <- SETSZ[gv]
  D  <- 2 * I / outer(s, s, "+")
  c(mean(D[upper.tri(D)]), m, as.numeric(sub))
}

score_partition <- function(memb) {
  # memb: data.table(gene, module). Returns per-module scores.
  memb <- memb[gene %chin% rownames(B)]
  memb[, {
    r <- dice_mean(gene)
    .(n_annotated = as.integer(r[2]), subsampled = r[3] > 0, H = r[1])
  }, by = module]
}

# --- the size-matched null ----------------------------------------------------
# Under a random partition the expected homogeneity of a module depends on ONE
# thing: how many annotated genes it holds. So it is tabulated once per study
# against a log-spaced grid of sizes and interpolated, instead of re-scoring
# every module of every clustering 100 times over. That is not an approximation
# of a different null -- it is the same null, computed where it actually varies,
# and it turns hours into seconds while giving a smoother estimate.
ALL_GENES <- rownames(B)
null_curve <- function(reps) {
  grid <- unique(pmin(c(2:9, round(10 ^ seq(1, log10(MAX_GENES), length.out = 18))),
                      length(ALL_GENES)))
  rbindlist(lapply(grid, function(m) {
    v <- vapply(seq_len(reps), function(i) dice_mean(sample(ALL_GENES, m))[1], numeric(1))
    data.table(m = m, H_null = mean(v, na.rm = TRUE), H_null_sd = sd(v, na.rm = TRUE))
  }))
}
say("tabulating the random-partition expectation over ", N_PERM, " draws per size")
NC <- null_curve(N_PERM)
print(NC[m <= 200], row.names = FALSE)

# Linear interpolation on log size; sizes past the grid use the last value, where
# the curve has already flattened.
null_for <- function(m) {
  m <- pmin(pmax(m, min(NC$m)), max(NC$m))
  approx(log(NC$m), NC$H_null, xout = log(m), rule = 2)$y
}

# --- the clusterings to score --------------------------------------------------
# Every cell the sweep produced, plus the partition the pipeline currently ships,
# so the two are on one scale.
tab <- fread(file.path(WORK, sprintf("%s.tab", STUDY)), header = FALSE,
             col.names = c("idx", "gene"))
tab[, gene := strip_version(gene)]

read_native_cls <- function(path) {
  # Native cluster format: records run from a cluster index to a '$' terminator
  # and may wrap across lines, so tokens are accumulated rather than split by row.
  ln <- readLines(path, warn = FALSE)
  i0 <- grep("^begin", ln); i1 <- grep("^\\)", ln)
  ln <- ln[(i0 + 1):(i1[i1 > i0][1] - 1)]
  tok <- unlist(strsplit(paste(ln, collapse = " "), "[ \t]+"))
  tok <- tok[nzchar(tok)]
  ends <- which(tok == "$")
  starts <- c(1L, head(ends, -1L) + 1L)
  rbindlist(lapply(seq_along(ends), function(c) {
    m <- tok[(starts[c] + 1L):(ends[c] - 1L)]
    if (!length(m)) return(NULL)
    data.table(idx = as.integer(m), module = sprintf("C%06d", c))
  }))
}

cells <- list()
sw <- file.path(WORK, sprintf("sweep_%s", STUDY))
if (dir.exists(sw)) {
  for (f in sort(list.files(sw, pattern = "^cls\\.k.*", full.names = TRUE))) {
    if (grepl("\\.secs$", f)) next
    nm <- sub("^cls\\.", "", basename(f))
    cells[[nm]] <- f
  }
}
say("clusterings to score: ", length(cells), " sweep cells + the shipped partition")

prod_file <- file.path(OUT_DIR, sprintf("mcl_%s_membership.tsv", STUDY))
res <- list(); per_module <- list()
score_one <- function(label, memb) {
  s <- score_partition(memb)
  ok <- s[!is.na(H)]
  if (!nrow(ok)) { say("  ", label, ": no scorable module"); return(NULL) }
  # Gene-weighted, so a clustering is not rewarded for having many tiny modules.
  obs   <- ok[, sum(H * n_annotated) / sum(n_annotated)]
  ok[, H_null := null_for(n_annotated)]
  nullw <- ok[, sum(H_null * n_annotated) / sum(n_annotated)]
  say(sprintf("  %-22s modules %6s  scored %6s  H %.4f  null %.4f  excess %+.4f",
              label, fmt_n(uniqueN(memb$module)), fmt_n(nrow(ok)), obs, nullw, obs - nullw))
  per_module[[label]] <<- cbind(clustering = label, s)
  # Carry the grid coordinates as their own columns so 38 can join on them
  # rather than re-parsing a label format that could drift.
  kk <- sub("^k([^.]+)\\..*$", "\\1", label)
  ii <- sub("^.*\\.I([0-9]+)$", "\\1", label)
  data.table(study = STUDY, clustering = label,
             knn = if (grepl("^k", label)) kk else NA_character_,
             inflation = if (grepl("\\.I[0-9]+$", label))
                           as.numeric(sub("^(.)(.*)$", "\\1.\\2", ii)) else NA_real_,
             n_modules = uniqueN(memb$module), n_modules_scored = nrow(ok),
             n_subsampled = ok[subsampled == TRUE, .N],
             H = obs, H_null_mean = nullw,
             H_excess = obs - nullw)
}

if (file.exists(prod_file)) {
  m <- fread(prod_file, select = c("gene", "module_name"))
  setnames(m, c("gene", "module"))
  m[, gene := strip_version(gene)]
  res[["shipped"]] <- score_one("shipped (-I 2, no knn)", m[module != "Unassigned"])
}
for (nm in names(cells)) {
  cl <- read_native_cls(cells[[nm]])
  m  <- merge(cl, tab, by = "idx")[, .(gene, module)]
  res[[nm]] <- score_one(nm, m)
}

out <- rbindlist(res, fill = TRUE)
write_tsv(out, file.path(OUT_DIR, sprintf("cluster_homogeneity_%s.tsv", STUDY)))
write_tsv(rbindlist(per_module, fill = TRUE),
          file.path(OUT_DIR, sprintf("cluster_homogeneity_%s_per_module.tsv", STUDY)))
say("")
print(out[order(-H_excess)], row.names = FALSE)
say("")
say("H_excess is the readout, not H: raw homogeneity rises as modules shrink,")
say("so only the excess over a size-matched random partition is evidence.")
say("done: ", STUDY)
