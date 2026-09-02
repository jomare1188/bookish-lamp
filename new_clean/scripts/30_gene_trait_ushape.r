#!/usr/bin/env Rscript
# =============================================================================
# 30_gene_trait_ushape.r — the U-shape contrast, genome-wide, at gene level
#
# WHY THIS EXISTS. Purple's design is stress-control-stress: 0 mM and 6 mM are
# deficiency and excess, 2 mM is the control. Every OTHER genome-wide test in
# this pipeline is monotone -- Pearson (07), mutual information (12), Spearman at
# module level (19) -- and a gene moved the SAME WAY by both stresses is
# structurally invisible to all of them.
#
# The contrast c(+1,-2,+1) already exists in the repo, but only ever on a handful
# of named genes (11_readouts/myb61/08, 11_readouts/module20/04) and once at
# module level. So the headline negative -- "of purple's 32 responsive genes only
# 2 are orthologs of a sugarcane-responsive gene" -- is at present a statement
# about MONOTONE responses only, and Soffic.09G0001580-9H is an uncalibrated
# anecdote: its p-value has never been placed against a genome-wide distribution.
#
# This script produces that distribution. Downstream, CLEAN_SELECTION=ushape in
# 08_conserved_cor_genes.r runs the existing conservation funnel on it.
#
# TWO MODELS, both reported:
#   blocked       expr ~ genotype + nlev, contrast on nlev. n=18, resid df=14.
#                 Genotype is a BLOCK, not a pool: 51NG3 is S. robustum and TAGZ
#                 is S. officinarum -- different species. This is the version
#                 comparable with the existing genome-wide negative.
#   per-genotype  n=9, resid df=6, identical to what myb61/08 and module20/04 do,
#                 so the anecdotes sit on the same scale as the distribution they
#                 are being calibrated against.
#
# DIRECTION IS NOT A DETAIL. A significant contrast means either:
#   est > 0  trough at the control -> UP at both 0N and 6N -> stress-induced
#   est < 0  peak at the control   -> DOWN at both extremes
# and the two known cases point opposite ways: Soffic.09G0001580-9H is a trough,
# every significant MYB61 copy is a peak. Both are carried; neither is discarded.
#
# NO EFFECT FLOOR is applied, unlike the Pearson rule's |r| >= 0.6, because the
# contrast has no natural standardised equivalent on the log2 scale. That makes
# this rule MORE permissive: a null under it is a stronger null, but any positive
# needs u_est inspected before it is believed. u_est is carried for that purpose.
# =============================================================================
suppressPackageStartupMessages({library(data.table)})
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE),
                                                 value = TRUE)[1])), "lib", "common.R"))

STUDY      <- env_opt("CLEAN_STUDY", "purple")
VST_PREFIX <- env_req("CLEAN_VST_PREFIX")
META_FILE  <- env_req("CLEAN_META")
OUT_FILE   <- env_req("CLEAN_OUT_FILE")
TRAIT      <- env_opt("CLEAN_SELECT_TRAIT", "treatment")
LEVELS     <- strsplit(env_opt("CLEAN_USHAPE_LEVELS", "0N,2N,6N"), ",")[[1]]
BLOCK      <- env_opt("CLEAN_USHAPE_BLOCK", "genotype")
FOCUS      <- env_opt("CLEAN_USHAPE_FOCUS", "Soffic.09G0001580-9H")
CVEC       <- as.numeric(strsplit(env_opt("CLEAN_USHAPE_CONTRAST", "1,-2,1"), ",")[[1]])

banner(sprintf("U-shape contrast c(%s) over %s — %s",
               paste(sprintf("%+g", CVEC), collapse = ","),
               paste(LEVELS, collapse = "/"), STUDY))

# --- the contrast, vectorised over genes -------------------------------------
# 11_readouts/myb61/08_myb61_expression_test.r::contrast_test does this one gene
# at a time with tapply. Over 170,740 genes that will not finish, so the same
# arithmetic is done once for all genes via the QR of the model matrix. It is
# EXACT, not an approximation -- 30_verify below requires it to reproduce that
# function's published values to full precision before anything here is trusted.
contrast_matrix <- function(Y, lev, block = NULL, cvec) {
  # Y: genes x samples. lev: factor of nitrogen level. block: optional factor.
  stopifnot(nlevels(lev) == length(cvec))
  X <- if (is.null(block)) model.matrix(~ 0 + lev) else model.matrix(~ 0 + lev + block)
  n <- ncol(Y); k <- qr(X)$rank
  df <- n - k
  if (df <= 0) stop("no residual degrees of freedom", call. = FALSE)

  qrX  <- qr(X)
  B    <- qr.coef(qrX, t(Y))                 # coefficients x genes
  FIT  <- X %*% B
  RES  <- t(Y) - FIT
  s2   <- colSums(RES^2) / df                # per gene

  # the contrast acts on the level means only; block terms get zero weight
  cfull <- numeric(ncol(X))
  cfull[seq_len(nlevels(lev))] <- cvec
  # var(c'b) = s2 * c' (X'X)^-1 c
  XtXi  <- chol2inv(qr.R(qrX))
  vscale <- as.numeric(t(cfull) %*% XtXi %*% cfull)

  est <- as.numeric(crossprod(cfull, B))
  se  <- sqrt(s2 * vscale)
  tt  <- est / se
  list(est = est, se = se, t = tt,
       p = 2 * pt(-abs(tt), df = df), df = df)
}

# --- inputs -------------------------------------------------------------------
say("reading VST from ", basename(VST_PREFIX), ".f32")
vst <- read_vst(VST_PREFIX)
say("  ", fmt_n(nrow(vst)), " genes x ", ncol(vst), " samples")

meta <- fread(META_FILE, header = TRUE)
setnames(meta, tolower(names(meta)))
m <- meta[match(colnames(vst), sample)]
if (any(is.na(m$sample)))
  stop("VST samples missing from ", basename(META_FILE), call. = FALSE)
for (cn in c(TRAIT, BLOCK))
  if (!cn %in% names(m)) stop("column '", cn, "' not in ", basename(META_FILE), call. = FALSE)

lev <- factor(m[[TRAIT]], levels = LEVELS)
if (any(is.na(lev)))
  stop("samples with a ", TRAIT, " outside {", paste(LEVELS, collapse = ","), "}: ",
       paste(unique(m[[TRAIT]][is.na(lev)]), collapse = ", "), call. = FALSE)
blk <- factor(m[[BLOCK]])
say("design: ", ncol(vst), " samples, ", nlevels(blk), " ", BLOCK, " x ",
    nlevels(lev), " levels")
print(table(setNames(data.frame(lev, blk), c(TRAIT, BLOCK))))

# --- blocked model, the primary -----------------------------------------------
say("")
say("blocked model: expr ~ ", BLOCK, " + ", TRAIT)
# A gene with zero variance has no estimable contrast; drop it here rather than
# letting it become an NaN downstream.
row_var <- rowSums((vst - rowMeans(vst))^2)
ok <- which(is.finite(row_var) & row_var > 0)
say("  ", fmt_n(length(ok)), " of ", fmt_n(nrow(vst)), " genes have non-zero variance")
res <- contrast_matrix(vst[ok, , drop = FALSE], lev, blk, CVEC)
say("  residual df = ", res$df)

out <- data.table(gene = rownames(vst)[ok], trait = TRAIT,
                  u_est = res$est, u_se = res$se, u_t = res$t, u_pval = res$p)
out[, u_padj := p.adjust(u_pval, method = "BH")]
out[, direction := fifelse(u_est > 0, "trough_at_control", "peak_at_control")]

# --- per-genotype, to match the anecdotes' scale -------------------------------
for (g in levels(blk)) {
  idx <- which(blk == g)
  say("per-genotype model: ", g, " (n = ", length(idx), ")")
  sub <- vst[ok, idx, drop = FALSE]
  r <- contrast_matrix(sub, droplevels(lev[idx]), NULL, CVEC)
  out[, (paste0("u_est_", g)) := r$est]
  out[, (paste0("u_p_", g))   := r$p]
  say("  residual df = ", r$df)
}

setorder(out, u_pval)
fwrite(out, OUT_FILE, sep = "\t")
say("")
say("wrote ", fmt_n(nrow(out)), " genes -> ", basename(OUT_FILE))

# --- what the distribution looks like -----------------------------------------
say("")
say("BH over ", fmt_n(nrow(out)), " tested genes:")
for (a in c(0.05, 0.1, 0.2))
  say(sprintf("  padj <= %.2f : %s  (trough %s, peak %s)", a,
              fmt_n(out[u_padj <= a, .N]),
              fmt_n(out[u_padj <= a & direction == "trough_at_control", .N]),
              fmt_n(out[u_padj <= a & direction == "peak_at_control", .N])))
say(sprintf("  raw p <= 0.05 (UNCORRECTED, the most generous reading): %s",
            fmt_n(out[u_pval <= 0.05, .N])))
say(sprintf("  expected by chance at raw p <= 0.05: %s", fmt_n(round(0.05 * nrow(out)))))

# Null calibration: a p-value distribution far from uniform, with nothing
# surviving BH, points at unmodelled structure rather than signal.
br <- table(cut(out$u_pval, breaks = seq(0, 1, 0.1)))
say("  p-value histogram (should be ~flat under the null):")
say("    ", paste(sprintf("%.0f%%", 100 * as.numeric(br) / nrow(out)), collapse = " "))

# --- the focus gene, calibrated ------------------------------------------------
if (nzchar(FOCUS) && FOCUS %chin% out$gene) {
  f <- out[gene == FOCUS]
  say("")
  say("FOCUS ", FOCUS, " — the anecdote, placed against the distribution:")
  say(sprintf("  blocked : u_est %+.4f  p %.4g  padj %.4g  [%s]",
              f$u_est, f$u_pval, f$u_padj, f$direction))
  say(sprintf("  rank %s of %s by raw p  (percentile %.2f)",
              fmt_n(which(out$gene == FOCUS)), fmt_n(nrow(out)),
              100 * which(out$gene == FOCUS) / nrow(out)))
  for (g in levels(blk))
    say(sprintf("  %-8s: u_est %+.4f  p %.4g", g,
                f[[paste0("u_est_", g)]], f[[paste0("u_p_", g)]]))
} else if (nzchar(FOCUS)) {
  say("FOCUS ", FOCUS, " is not in the tested set")
}
# --- the summary the argument actually needs ----------------------------------
# Not "did anything reach significance" but "how much non-monotone signal is
# there, and how much of it can this design resolve" -- the two are different and
# only the pair of them bounds the negative result.
n_raw <- out[u_pval <= 0.05, .N]
n_exp <- 0.05 * nrow(out)
foc   <- if (nzchar(FOCUS) && FOCUS %chin% out$gene) out[gene == FOCUS] else NULL
summ <- data.table(
  metric = c("genes_tested", "residual_df_blocked",
             "padj_le_0.05", "padj_le_0.05_trough", "padj_le_0.05_peak",
             "raw_p_le_0.05_UNCORRECTED", "raw_p_le_0.05_expected_by_chance",
             "raw_p_le_0.05_fold_over_chance", "excess_genes_over_chance",
             "focus_gene", "focus_rank", "focus_percentile",
             "focus_pval", "focus_padj", "focus_direction"),
  value = c(nrow(out), res$df,
            out[u_padj <= 0.05, .N],
            out[u_padj <= 0.05 & direction == "trough_at_control", .N],
            out[u_padj <= 0.05 & direction == "peak_at_control", .N],
            n_raw, round(n_exp), round(n_raw / n_exp, 3), round(n_raw - n_exp),
            if (is.null(foc)) NA else FOCUS,
            if (is.null(foc)) NA else which(out$gene == FOCUS),
            if (is.null(foc)) NA else round(100 * which(out$gene == FOCUS) / nrow(out), 3),
            if (is.null(foc)) NA else signif(foc$u_pval, 4),
            if (is.null(foc)) NA else signif(foc$u_padj, 4),
            if (is.null(foc)) NA else foc$direction))
sfile <- sub("\\.tsv$", "_summary.tsv", OUT_FILE)
fwrite(summ, sfile, sep = "\t")
say("wrote ", basename(sfile))
say("done")
