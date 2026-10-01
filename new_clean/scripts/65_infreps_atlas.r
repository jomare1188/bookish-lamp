#!/usr/bin/env Rscript
# =============================================================================
# 65_infreps_atlas.r -- what the quantifier was UNSURE about, per gene and per
# orthogroup, from the 30 Gibbs replicates that have been sitting unused on disk.
#
# WHY THIS EXISTS. Every "the networks diverge" result in this project is
# entangled with read-assignment ambiguity, because in an 8-12x polyploid reads
# cannot be confidently assigned among near-identical copies. The stage's own
# control says so: in 29_dnds/12_decoupling_test.r, frac_unique_min carries
# coefficient -1.53 (p 6.8e-08) in sugarcane and -1.85 (p 1.6e-31) in purple --
# LESS UNIQUELY MAPPABLE MEANS MORE APPARENT NETWORK DIVERGENCE, which is the
# signature of an artefact rather than of biology.
#
# WHY frac_unique IS NOT ENOUGH. It is a sequence proxy, computed from canonical
# 31-mers shared with SIBLING COPIES ONLY (29_dnds/10_kmer_uniqueness.py), so it
# is blind to ambiguity against a non-sibling paralog, and it exists for just 29%
# of sugarcane and 49% of purple network nodes. What is here is the quantifier's
# OWN uncertainty, propagated from the actual reads, genome-wide, for every gene.
#
# THE THREE MEASUREMENTS
#
#   infrv            fishpond's inferential relative variance, per gene. How much
#                    of this gene's count variance is read assignment rather than
#                    biology. The published definition, via computeInfRV() -- not
#                    reimplemented here.
#
#   boot_r_within    per gene PAIR inside an orthogroup: the correlation of the
#                    two genes' counts ACROSS THE 30 REPLICATES WITHIN A SAMPLE,
#                    averaged over samples. This is the direct read-stealing
#                    diagnostic. When the EM cannot tell two copies apart it moves
#                    reads between them: one rises exactly as the other falls
#                    while their sum stays put, so STRONGLY NEGATIVE MEANS NOT
#                    SEPARATELY ESTIMABLE. Near zero means the copies are
#                    distinguishable and any divergence between them is real.
#
#   infrv_gain       infrv averaged over an orthogroup's members, divided by the
#                    infrv of the orthogroup's SUMMED profile. How much
#                    uncertainty collapsing to the orthogroup removes. This is the
#                    number that justifies -- or refuses to justify -- treating
#                    the orthogroup as the unit of analysis.
#
# THE NULL IS NOT OPTIONAL. boot_r_within is a correlation of 30 numbers; it has
# a distribution even between unrelated genes. The reference is
# expression-decile-matched RANDOM gene pairs, drawn exactly as
# 29_dnds/12_decoupling_test.r:149-172 draws its own, and it must centre on ~0.
# If it does not, the statistic is broken and nothing below it means anything.
#
# WHAT THIS IS NOT. Not a differential-expression test (that is 66), not a filter
# (nothing here removes a gene from any analysis), and not a replacement for
# frac_unique -- the two measure different things and the crosswalk carries both.
#
# THE COUNTS ARE THE PIPELINE'S OWN. tximport is run with
# countsFromAbundance="lengthScaledTPM" over ALL samples at once, which
# reproduces nf-core's salmon.merged.gene_counts_length_scaled.tsv to 4.7e-10 and
# counts(dds) exactly after rounding -- asserted below. The scaling depends on
# mean transcript length across the whole dataset, so a subset of samples gives
# DIFFERENT counts; that is why this cannot be chunked by sample.
#
# RUN: through run.sh -> ./run.sh infreps <study>
# =============================================================================

suppressMessages({
  library(tximport)
  library(fishpond)
  library(SummarizedExperiment)
  library(data.table)
})

source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE),
                                                 value = TRUE)[1])), "lib", "common.R"))

STUDY      <- env_req("CLEAN_STUDY")
QUANT_DIR  <- env_req("CLEAN_QUANT_DIR")
META_FILE  <- env_req("CLEAN_META")
ORTHO_FILE <- env_req("CLEAN_ORTHOGROUPS")
OG_COL     <- env_req("CLEAN_OG_SPECIES")
OUT_DIR    <- env_req("CLEAN_OUT_DIR")
STRIP_VER  <- env_flag("CLEAN_STRIP_VERSION")
MAX_COPIES <- as.integer(env_num("CLEAN_MAX_OG_COPIES", 20))
NULL_PAIRS <- as.integer(env_num("CLEAN_NULL_PAIRS", 200000))
SEED       <- as.integer(env_num("CLEAN_SEED", 1))

banner(paste("inferential replicates:", STUDY))
ensure_dir(OUT_DIR)
set.seed(SEED)

# --- load -------------------------------------------------------------------
meta <- fread(META_FILE, header = TRUE)
setnames(meta, tolower(names(meta)))
samples <- as.character(meta$sample)
qf <- file.path(QUANT_DIR, samples, "quant.sf")
names(qf) <- samples
if (!all(file.exists(qf)))
  stop("missing quant.sf for: ",
       paste(samples[!file.exists(qf)], collapse = ", "), call. = FALSE)

t2g_f <- file.path(QUANT_DIR, "salmon.merged.tx2gene.tsv")
if (!file.exists(t2g_f))
  stop("tx2gene not found: ", t2g_f, call. = FALSE)
t2g <- fread(t2g_f, header = TRUE)
say(sprintf("%d samples | tx2gene %s rows", length(samples), fmt_n(nrow(t2g))))

say("tximport with inferential replicates (the slow step)")
txi <- tximport(qf, type = "salmon", tx2gene = as.data.frame(t2g[, 1:2]),
                countsFromAbundance = "lengthScaledTPM", dropInfReps = FALSE)
if (is.null(txi$infReps))
  stop("tximport returned no inferential replicates -- was salmon run with ",
       "--numGibbsSamples/--numBootstraps?", call. = FALSE)
N_REP <- ncol(txi$infReps[[1]])
say(sprintf("counts %s x %d | %d replicates per sample",
            fmt_n(nrow(txi$counts)), ncol(txi$counts), N_REP))

# THE RECONCILIATION. Everything downstream describes the counts the networks
# were built from only if this holds. Asserted against the merged table rather
# than the dds so no DESeq2 is needed; the two are equal after rounding.
mg_f <- file.path(QUANT_DIR, "salmon.merged.gene_counts_length_scaled.tsv")
if (!file.exists(mg_f))
  stop("cannot reconcile: ", mg_f, " is missing", call. = FALSE)
mg <- as.data.frame(fread(mg_f), stringsAsFactors = FALSE)
rownames(mg) <- mg[[1]]
mg <- mg[, intersect(colnames(mg), samples), drop = FALSE]
gg <- intersect(rownames(txi$counts), rownames(mg))
dev <- max(abs(round(txi$counts[gg, colnames(mg), drop = FALSE]) -
               round(as.matrix(mg[gg, , drop = FALSE]))))
if (length(gg) < 0.99 * nrow(txi$counts) || dev != 0)
  stop("tximport does not reproduce ", basename(mg_f), ": ", fmt_n(length(gg)),
       " genes matched, max rounded deviation ", dev,
       "\n  countsFromAbundance is wrong, or the quantification moved", call. = FALSE)
say(sprintf("reconciled against %s: %s genes, rounded deviation 0",
            basename(mg_f), fmt_n(length(gg))))
rm(mg); invisible(gc())

genes <- rownames(txi$counts)
if (STRIP_VER) {
  bare <- strip_version(genes)
  if (anyDuplicated(bare))
    stop("stripping versions collides ", fmt_n(sum(duplicated(bare))),
         " gene ids", call. = FALSE)
  genes <- bare
}

# Mean expression level, for the decile-matched null. Same construction as
# 29_dnds/12_decoupling_test.r:34-38 so the two nulls are comparable.
mean_log <- rowMeans(log2(txi$counts + 1))
names(mean_log) <- genes

# --- the orthogroup partition -----------------------------------------------
og <- as.data.table(read_orthogroups_tsv(ORTHO_FILE))
og[, Gene := strip_version(sub("\\.p[0-9]+$", "", Gene))]
og <- og[Species == OG_COL & Gene %chin% genes, .(Orthogroup, Gene)]
og <- unique(og)
sz <- og[, .N, by = Orthogroup]
say(sprintf("orthogroups with >=1 expressed gene: %s | multi-copy: %s",
            fmt_n(nrow(sz)), fmt_n(sz[N >= 2, .N])))

# Oversized groups are excluded from the PAIR statistic only -- they are
# quadratically expensive and, past ~20 copies, lumped superfamilies rather than
# one locus retained several times (29_dnds/08_homeolog_families.py:77-80). They
# keep their per-gene infrv and their collapsed profile.
keep_og <- sz[N >= 2 & N <= MAX_COPIES, Orthogroup]
over <- sz[N > MAX_COPIES, .N]
say(sprintf("pairs drawn from %s groups (<= %d copies); %s oversized groups skipped",
            fmt_n(length(keep_og)), MAX_COPIES, fmt_n(over)))

gidx <- setNames(seq_along(genes), genes)
pair_dt <- og[Orthogroup %chin% keep_og][order(Orthogroup, Gene)][
  , {
      k <- .N
      i <- rep(seq_len(k - 1L), times = rev(seq_len(k - 1L)))
      j <- unlist(lapply(seq_len(k - 1L), function(a) (a + 1L):k), use.names = FALSE)
      .(gene_a = Gene[i], gene_b = Gene[j])
    }, by = Orthogroup]
say(sprintf("within-orthogroup pairs: %s", fmt_n(nrow(pair_dt))))

# --- boot_r_within, vectorised one sample at a time -------------------------
# For a sample, infReps[[s]] is genes x replicates. The correlation of two genes
# across replicates is then a dot product of two row-standardised rows, so every
# pair in the table is one elementwise multiply and one rowSums -- no loop over
# pairs, which at 169k-247k pairs x 48 samples would not finish.
row_standardise <- function(M) {
  M <- M - rowMeans(M)
  nrm <- sqrt(rowSums(M * M))
  bad <- !is.finite(nrm) | nrm < 1e-12     # a gene the sampler never moved
  M <- M / nrm
  M[bad, ] <- NA_real_
  M
}

#
# A pair with n_usable == 0 is NOT missing data. It means that in every sample at
# least one of the two genes had IDENTICAL counts in all 30 replicates -- the
# sampler never moved it, because its reads were never in contention. That is the
# strongest evidence of unambiguous quantification the data can offer, so it is
# carried as its own category rather than dropped as an NA.
pair_cor_over_samples <- function(ia, ib) {
  acc <- matrix(0, nrow = length(ia), ncol = 2)   # [, 1] sum, [, 2] n usable
  for (s in seq_along(samples)) {
    Z <- row_standardise(txi$infReps[[s]])
    v <- rowSums(Z[ia, , drop = FALSE] * Z[ib, , drop = FALSE])
    ok <- is.finite(v)
    acc[ok, 1] <- acc[ok, 1] + v[ok]
    acc[ok, 2] <- acc[ok, 2] + 1
  }
  list(r = ifelse(acc[, 2] > 0, acc[, 1] / pmax(acc[, 2], 1), NA_real_),
       n_usable = as.integer(acc[, 2]))
}

say("boot_r_within over the within-orthogroup pairs")
ia <- gidx[pair_dt$gene_a]; ib <- gidx[pair_dt$gene_b]
pc <- pair_cor_over_samples(ia, ib)
pair_dt[, boot_r_within := pc$r]
pair_dt[, n_samples_usable := pc$n_usable]
pair_dt[, verdict := fifelse(n_samples_usable == 0L, "determined",
                     fifelse(boot_r_within < -0.5, "read_stealing", "separable"))]
pair_dt[, mean_log_a := mean_log[gene_a]]
pair_dt[, mean_log_b := mean_log[gene_b]]
n_det <- pair_dt[n_samples_usable == 0L, .N]
n_measured <- nrow(pair_dt) - n_det
say(sprintf("  median %.4f over %s measurable pairs | %s read-stealing (%.1f%%) | %s determined",
            median(pair_dt$boot_r_within, na.rm = TRUE), fmt_n(n_measured),
            fmt_n(pair_dt[verdict == "read_stealing", .N]),
            100 * pair_dt[verdict == "read_stealing", .N] / max(n_measured, 1),
            fmt_n(n_det)))

# --- the null: expression-decile-matched random pairs ------------------------
say(sprintf("null: %s expression-decile-matched random pairs", fmt_n(NULL_PAIRS)))
dec <- cut(rank(mean_log, ties.method = "first"),
           breaks = quantile(rank(mean_log, ties.method = "first"),
                             probs = seq(0, 1, 0.1)),
           include.lowest = TRUE, labels = FALSE)
by_dec <- split(seq_along(genes), dec)
take <- sample(nrow(pair_dt), min(NULL_PAIRS, nrow(pair_dt)), replace = NULL_PAIRS > nrow(pair_dt))
na_idx <- vapply(dec[gidx[pair_dt$gene_a[take]]], function(d) sample(by_dec[[d]], 1L), 1L)
nb_idx <- vapply(dec[gidx[pair_dt$gene_b[take]]], function(d) sample(by_dec[[d]], 1L), 1L)
same_og <- og[, .(Orthogroup, Gene)][match(genes[na_idx], Gene), Orthogroup] ==
           og[, .(Orthogroup, Gene)][match(genes[nb_idx], Gene), Orthogroup]
null_r <- pair_cor_over_samples(na_idx, nb_idx)$r
null_r[which(same_og)] <- NA_real_      # a random draw that landed on a sibling
say(sprintf("  null median %.4f | mean %.4f | %s drawn, %s usable",
            median(null_r, na.rm = TRUE), mean(null_r, na.rm = TRUE),
            fmt_n(length(null_r)), fmt_n(sum(is.finite(null_r)))))
if (abs(median(null_r, na.rm = TRUE)) > 0.15)
  stop("the null median is ", sprintf("%.3f", median(null_r, na.rm = TRUE)),
       ", not ~0 -- boot_r_within is measuring something other than read ",
       "assignment and nothing downstream of it can be trusted", call. = FALSE)

# --- per-gene InfRV ----------------------------------------------------------
# computeInfRV wants infRep1..K as genes x samples assays; tximport hands back
# genes x replicates per SAMPLE, so the orientation is transposed here once.
say("per-gene InfRV")
mk_assays <- function(mats, ids) {
  out <- vector("list", N_REP)
  for (k in seq_len(N_REP)) {
    A <- vapply(mats, function(M) M[, k], numeric(nrow(mats[[1]])))
    dimnames(A) <- list(ids, samples)
    out[[k]] <- A
  }
  setNames(out, paste0("infRep", seq_len(N_REP)))
}
cm <- txi$counts; dimnames(cm) <- list(genes, samples)
se <- SummarizedExperiment(assays = c(list(counts = cm), mk_assays(txi$infReps, genes)))
se <- computeInfRV(se)
gene_dt <- data.table(gene = genes,
                      infrv = as.numeric(mcols(se)$meanInfRV),
                      mean_log_count = as.numeric(mean_log))
say(sprintf("  infrv: median %.4f | 90th pct %.4f | max %.3f",
            median(gene_dt$infrv, na.rm = TRUE),
            quantile(gene_dt$infrv, 0.9, na.rm = TRUE),
            max(gene_dt$infrv, na.rm = TRUE)))
rm(se); invisible(gc())

# --- collapsed orthogroup profiles, and the gain ----------------------------
say("InfRV of the summed orthogroup profiles")
og_ids <- sort(unique(og$Orthogroup))
oidx <- setNames(seq_along(og_ids), og_ids)
row_of <- oidx[og$Orthogroup]
col_of <- gidx[og$Gene]
agg <- function(M) {
  A <- matrix(0, nrow = length(og_ids), ncol = ncol(M))
  for (k in seq_len(ncol(M))) A[, k] <- rowsum(M[col_of, k], row_of, reorder = FALSE)[, 1]
  A
}
og_counts <- agg(cm)
dimnames(og_counts) <- list(og_ids, samples)
se_og <- SummarizedExperiment(
  assays = c(list(counts = og_counts),
             mk_assays(lapply(txi$infReps, function(M) agg(M)), og_ids)))
se_og <- computeInfRV(se_og)
og_infrv <- setNames(as.numeric(mcols(se_og)$meanInfRV), og_ids)
rm(se_og); invisible(gc())

member <- merge(og, gene_dt[, .(Gene = gene, infrv)], by = "Gene")
og_dt <- member[, .(n_members = .N,
                    infrv_member_mean = mean(infrv, na.rm = TRUE),
                    infrv_member_max = max(infrv, na.rm = TRUE)),
                by = Orthogroup]
og_dt[, infrv_collapsed := og_infrv[Orthogroup]]
og_dt[, infrv_gain := infrv_member_mean / infrv_collapsed]
pr <- pair_dt[, .(boot_r_within = mean(boot_r_within, na.rm = TRUE),
                  n_pairs = .N), by = Orthogroup]
og_dt <- merge(og_dt, pr, by = "Orthogroup", all.x = TRUE)
setnames(og_dt, "Orthogroup", "og_2sp")

m <- og_dt[n_members >= 2 & is.finite(infrv_gain)]
say(sprintf("  multi-copy groups: median infrv_gain %.3f | %.1f%% above 1 (collapsing helps)",
            median(m$infrv_gain), 100 * mean(m$infrv_gain > 1)))

# --- write -------------------------------------------------------------------
write_tsv(gene_dt[order(-infrv)], file.path(OUT_DIR, sprintf("%s_infrv.tsv", STUDY)))
setnames(pair_dt, "Orthogroup", "og_2sp")
write_tsv(pair_dt[order(boot_r_within)],
          file.path(OUT_DIR, sprintf("%s_readsteal_pairs.tsv", STUDY)))
write_tsv(og_dt[order(-infrv_gain)],
          file.path(OUT_DIR, sprintf("%s_og_uncertainty.tsv", STUDY)))
write_tsv(data.table(boot_r_within = null_r[is.finite(null_r)]),
          file.path(OUT_DIR, sprintf("%s_readsteal_null.tsv", STUDY)))

summ <- data.table(
  metric = c("study", "n_samples", "n_replicates", "n_genes", "n_orthogroups",
             "n_multicopy_orthogroups", "n_oversized_skipped", "max_og_copies",
             "n_pairs", "n_pairs_determined", "pair_boot_r_median",
             "pair_pct_read_stealing",
             "null_boot_r_median", "null_boot_r_mean", "n_null_pairs",
             "infrv_median", "infrv_p90", "infrv_gain_median",
             "pct_gain_above_1"),
  value = c(STUDY, length(samples), N_REP, length(genes), nrow(sz),
            sz[N >= 2, .N], over, MAX_COPIES, nrow(pair_dt), n_det,
            sprintf("%.4f", median(pair_dt$boot_r_within, na.rm = TRUE)),
            sprintf("%.2f", 100 * pair_dt[verdict == "read_stealing", .N] /
                      max(n_measured, 1)),
            sprintf("%.4f", median(null_r, na.rm = TRUE)),
            sprintf("%.4f", mean(null_r, na.rm = TRUE)),
            sum(is.finite(null_r)),
            sprintf("%.4f", median(gene_dt$infrv, na.rm = TRUE)),
            sprintf("%.4f", quantile(gene_dt$infrv, 0.9, na.rm = TRUE)),
            sprintf("%.4f", median(m$infrv_gain)),
            sprintf("%.2f", 100 * mean(m$infrv_gain > 1))))
write_tsv(summ, file.path(OUT_DIR, sprintf("%s_infreps_summary.tsv", STUDY)))

say("done")
