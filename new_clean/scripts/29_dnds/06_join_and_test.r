#!/usr/bin/env Rscript
# ============================================================================
# 06_join_and_test.r -- does sequence evolution follow network conservation?
#
# Joins the per-gene dN/dS values onto the network predictors and runs the
# tests in the order that makes the answer interpretable. The order matters:
# the marginal correlation (test 1) is NOT the result. Degree, expression
# level and edge conservation are all mutually correlated, and expression
# level is the strongest known predictor of omega in every organism it has
# been measured in -- so the PARTIAL effect of conservation, after those are
# in the model (test 2), is what the hypothesis actually claims.
#
# The sugarcane<->purple pair is never used per gene. It is 95-99% identical,
# ~40% of genes have zero synonymous differences, and a per-gene ratio there
# is undefined or noise. It enters only through test 5, where NG86 COUNTS are
# summed within a bin and one ratio is formed per bin.
# ============================================================================
suppressPackageStartupMessages({library(data.table)})

env <- function(k, d = NA) { v <- Sys.getenv(k); if (nzchar(v)) v else d }
OUTDIR  <- env("OUTDIR")
DS_MIN  <- as.numeric(env("DNDS_DS_MIN"))
DS_MAX  <- as.numeric(env("DNDS_DS_MAX"))
OM_FLAG <- as.numeric(env("DNDS_OMEGA_FLAG"))

say <- function(...) cat(sprintf(...), "\n", sep = "")

# ---------------------------------------------------------------- inputs ----
trip <- fread(env("TRIPLETS"))
raw  <- fread(file.path(OUTDIR, "dnds_raw.tsv"))
d    <- merge(trip, raw, by = "orthogroup")
say("== %d triplets with dN/dS", nrow(d))

cons_sc <- fread(file.path(OUTDIR, "conserved_degree_sugarcane.tsv"))
setnames(cons_sc, c("gene", "n_edges", "n_conserved", "frac_conserved"),
         c("sugarcane_gene", "sc_n_edges", "sc_n_conserved", "sc_frac_conserved"))
d <- merge(d, cons_sc, by = "sugarcane_gene", all.x = TRUE)

# The purple pass over a 37 GB edge table takes ~40 min, so it may legitimately
# still be running. An empty or partial file is treated as absent rather than
# joined -- a half-written conservation column would silently corrupt the
# replicate test rather than fail it.
pu_path <- file.path(OUTDIR, "conserved_degree_purple.tsv")
have_purple <- file.exists(pu_path) && file.size(pu_path) > 0 &&
  identical(names(fread(pu_path, nrows = 1)),
            c("gene", "n_edges", "n_conserved", "frac_conserved"))
if (have_purple) {
  cons_pu <- fread(pu_path)
  setnames(cons_pu, c("gene", "n_edges", "n_conserved", "frac_conserved"),
           c("purple_gene", "pu_n_edges", "pu_n_conserved", "pu_frac_conserved"))
  d <- merge(d, cons_pu, by = "purple_gene", all.x = TRUE)
} else {
  say("!! conserved_degree_purple.tsv absent -- the purple replicate is skipped")
}

nm <- fread(env("NODE_METRICS_sugarcane"), select = c("gene", "degree", "strength"))
setnames(nm, c("sugarcane_gene", "sc_degree", "sc_strength"))
d <- merge(d, nm, by = "sugarcane_gene", all.x = TRUE)

hubs <- fread(env("HUBS_sugarcane"), select = "gene")
d[, sc_is_hub := sugarcane_gene %in% hubs$gene]

# expression: the confounder that must be in the model, not a nicety
mean_log_tpm <- function(path, strip_version) {
  tpm <- fread(path)
  id <- tpm[[1]]
  if (strip_version) id <- sub("\\.v[0-9.]+$", "", id)
  m <- as.matrix(tpm[, -c(1, 2), with = FALSE])
  data.table(gene = id, mean_log_tpm = rowMeans(log2(m + 1)))
}
e_sc <- mean_log_tpm(env("TPM_sugarcane"), TRUE)
setnames(e_sc, c("sugarcane_gene", "sc_mean_log_tpm"))
d <- merge(d, unique(e_sc, by = "sugarcane_gene"), by = "sugarcane_gene", all.x = TRUE)

# ------------------------------------------- a better constraint readout ----
# omega = dN/dS assumes amino-acid divergence scales 1:1 with synonymous
# divergence. It does not here: regressing log dN on log dS ALONE gives a slope
# of 0.511 (CI 0.492-0.530), and 0.725 once GC3, expression and length are also
# in the model -- the gap between the two is how much of the dS signal is really
# base composition. Either way the slope is below 1, so omega OVER-CORRECTS, and
# because GC3 inflates dS about threefold (spearman(GC3, log dS) = +0.744 against
# +0.176 for log dN) dividing by dS injects composition into the constraint
# measure instead of removing the neutral rate.
#
# constraint_score is the residual of log(dN) once the neutral rate and
# composition are regressed out with FITTED rather than assumed coefficients:
# "more or less amino-acid change than expected for a gene with this synonymous
# rate, base composition, expression level and length". Negative = constrained.
#
# It is built here, not in a later stage, so every downstream analysis can use it.
d[, constraint_score := NA_real_]
fit_rows <- d[, which(sc_sb_dN > 0 & sc_sb_dS >= DS_MIN & sc_sb_dS <= DS_MAX &
                        !is.na(gc3) & !is.na(sc_mean_log_tpm) & !is.na(aln_codons))]
if (length(fit_rows) > 100) {
  fit <- lm(log(sc_sb_dN) ~ log(sc_sb_dS) + gc3 + sc_mean_log_tpm + aln_codons,
            data = d[fit_rows])
  d[fit_rows, constraint_score := residuals(fit)]
  b <- coef(fit)[["log(sc_sb_dS)"]]
  say("== constraint_score: log(dN) ~ log(dS) + GC3 + expression + length")
  say("   fitted dS slope = %.3f (omega assumes 1.000) over %d genes; R2 = %.3f",
      b, length(fit_rows), summary(fit)$r.squared)
  # Both are oriented so that LOW = more constrained (low omega; negative
  # residual = less amino-acid change than expected), so they should agree
  # POSITIVELY. A strong but imperfect agreement is the point: they are two
  # different corrections of the same quantity, and the disagreement is exactly
  # the dS/GC3 distortion omega carries and this does not.
  say("   spearman(constraint_score, log omega) = %+.3f  (both: low = constrained)",
      cor(d[fit_rows, constraint_score], log(d[fit_rows, sc_sb_omega]),
          method = "spearman", use = "complete.obs"))
}

# --------------------------------------------------------------- filters ----
# Every filter is counted and written out. None is applied silently.
n0 <- nrow(d)
d[, ds_ok_sc  := sc_sb_dS >= DS_MIN & sc_sb_dS <= DS_MAX]
d[, ds_ok_pu  := pu_sb_dS >= DS_MIN & pu_sb_dS <= DS_MAX]
d[, omega_outlier := sc_sb_omega > OM_FLAG]
d[, in_network := !is.na(sc_degree) & sc_degree > 0]

filt <- data.table(
  step = c("triplets with dN/dS", "sugarcane gene in the network",
           sprintf("dS(sugarcane,sorghum) in [%g,%g]", DS_MIN, DS_MAX),
           sprintf("omega <= %g (outliers reported separately)", OM_FLAG),
           "has conservation + degree + expression"),
  n = c(n0,
        d[in_network == TRUE, .N],
        d[in_network == TRUE & ds_ok_sc == TRUE, .N],
        d[in_network == TRUE & ds_ok_sc == TRUE & omega_outlier == FALSE, .N],
        d[in_network == TRUE & ds_ok_sc == TRUE & omega_outlier == FALSE &
            !is.na(sc_frac_conserved) & !is.na(sc_mean_log_tpm), .N]))
say("== filters"); print(filt)
fwrite(filt, file.path(OUTDIR, "dnds_filter_report.tsv"), sep = "\t")

fwrite(d, file.path(OUTDIR, "dnds_gene_table.tsv"), sep = "\t")
say("== wrote dnds_gene_table.tsv (%d rows, %d cols)", nrow(d), ncol(d))

a <- d[in_network == TRUE & ds_ok_sc == TRUE & omega_outlier == FALSE &
         !is.na(sc_frac_conserved) & !is.na(sc_mean_log_tpm) & sc_sb_omega > 0]
say("== analysis set: %d genes", nrow(a))

res <- list()
add <- function(test, statistic, value, n, note = "")
  res[[length(res) + 1]] <<- data.table(test, statistic, value, n, note)

# --- 1. the marginal claim --------------------------------------------------
ct <- suppressWarnings(cor.test(a$sc_sb_omega, a$sc_frac_conserved, method = "spearman"))
add("1 marginal: omega(sc,sorghum) ~ frac_conserved", "spearman_rho",
    unname(ct$estimate), nrow(a), sprintf("p = %.3g", ct$p.value))
w <- wilcox.test(sc_sb_omega ~ (sc_n_conserved > 0), data = a)
add("1 marginal: omega, >=1 conserved edge vs none", "wilcox_p", w$p.value, nrow(a),
    sprintf("median omega: %.4f conserved vs %.4f not",
            a[sc_n_conserved > 0, median(sc_sb_omega)],
            a[sc_n_conserved == 0, median(sc_sb_omega)]))

if (have_purple) {
  b <- d[ds_ok_pu == TRUE & pu_sb_omega > 0 & pu_sb_omega <= OM_FLAG &
           !is.na(pu_frac_conserved)]
  ct2 <- suppressWarnings(cor.test(b$pu_sb_omega, b$pu_frac_conserved, method = "spearman"))
  add("1 replicate: omega(purple,sorghum) ~ purple frac_conserved", "spearman_rho",
      unname(ct2$estimate), nrow(b), sprintf("p = %.3g", ct2$p.value))
}

# --- 2. the partial effect: THIS is the result ------------------------------
m <- lm(log(sc_sb_omega) ~ sc_frac_conserved + log(sc_degree) + sc_mean_log_tpm +
          aln_codons + gc3, data = a)
cf <- summary(m)$coefficients
ci <- confint(m)
for (v in rownames(cf)[-1])
  add("2 partial: log(omega) ~ conservation + degree + expression + length + GC3",
      v, cf[v, "Estimate"], nrow(a),
      sprintf("CI [%.4f, %.4f], p = %.3g", ci[v, 1], ci[v, 2], cf[v, "Pr(>|t|)"]))
add("2 model", "adj_r_squared", summary(m)$adj.r.squared, nrow(a), "")

m0 <- lm(log(sc_sb_omega) ~ log(sc_degree) + sc_mean_log_tpm + aln_codons + gc3, data = a)
add("2 model", "delta_adj_r2_from_conservation",
    summary(m)$adj.r.squared - summary(m0)$adj.r.squared, nrow(a),
    "how much conservation adds once the confounders are already in")

# --- 3. within degree deciles ----------------------------------------------
a[, deg_decile := cut(rank(sc_degree, ties.method = "first"),
                      breaks = quantile(rank(sc_degree, ties.method = "first"),
                                        probs = seq(0, 1, 0.1)),
                      include.lowest = TRUE, labels = 1:10)]
strat <- a[, {
  ct <- if (.N >= 10)
          suppressWarnings(cor.test(sc_sb_omega, sc_frac_conserved, method = "spearman"))
        else list(estimate = NA_real_, p.value = NA_real_)
  .(n = .N, rho = as.numeric(unname(ct$estimate)), p = as.numeric(ct$p.value),
    median_omega = as.numeric(median(sc_sb_omega)),
    median_degree = as.numeric(median(sc_degree)))
}, by = deg_decile][order(deg_decile)]
say("== 3. within degree deciles"); print(strat)
fwrite(strat, file.path(OUTDIR, "dnds_by_degree_decile.tsv"), sep = "\t")
add("3 stratified", "deciles_with_rho<0", sum(strat$rho < 0, na.rm = TRUE), nrow(a),
    "consistency of sign across the 10 degree deciles")

# --- 4. hubs: the positive control -----------------------------------------
wh <- wilcox.test(sc_sb_omega ~ sc_is_hub, data = a)
add("4 control: omega, hub vs non-hub", "wilcox_p", wh$p.value, nrow(a),
    sprintf("median omega: %.4f hub (n=%d) vs %.4f non-hub",
            a[sc_is_hub == TRUE, median(sc_sb_omega)], a[sc_is_hub == TRUE, .N],
            a[sc_is_hub == FALSE, median(sc_sb_omega)]))

# --- 5. the Saccharum pair, aggregated (never per gene) --------------------
# omega for a BIN = (sum Nd / sum N) / (sum Sd / sum S). Ratios cannot be
# averaged across genes when the denominator is often zero; counts can be summed.
binned_omega <- function(dt) {
  pn <- sum(dt$ng86_Nd) / sum(dt$ng86_N)
  ps <- sum(dt$ng86_Sd) / sum(dt$ng86_S)
  if (ps <= 0) return(NA_real_)
  pn / ps
}
boot_ci <- function(dt, R = 2000) {
  if (nrow(dt) < 5) return(c(NA_real_, NA_real_))
  set.seed(1)
  v <- replicate(R, binned_omega(dt[sample.int(.N, .N, replace = TRUE)]))
  unname(quantile(v, c(0.025, 0.975), na.rm = TRUE))
}
p <- d[in_network == TRUE & !is.na(sc_frac_conserved)]
p[, cons_bin := cut(sc_frac_conserved, breaks = c(-Inf, 0, 0.05, 0.1, 0.2, Inf),
                    labels = c("0", "(0,0.05]", "(0.05,0.1]", "(0.1,0.2]", ">0.2"))]
pairbins <- p[, {
  ci <- boot_ci(.SD)
  .(n_genes = .N, Nd = sum(ng86_Nd), Sd = sum(ng86_Sd),
    N = sum(ng86_N), S = sum(ng86_S),
    omega = binned_omega(.SD), lo = ci[1], hi = ci[2])
}, by = cons_bin][order(cons_bin)]
say("== 5. sugarcane<->LA purple, NG86 counts summed within conservation bins")
print(pairbins)
fwrite(pairbins, file.path(OUTDIR, "dnds_saccharum_pair_binned.tsv"), sep = "\t")

# --- 5b. is that trend conservation, or is it expression? ------------------
# The conservation bins are NOT otherwise comparable: median log2 TPM runs from
# 0.45 in the unconserved bin to 2.16 in the most conserved, and GC3 from 0.77
# to 0.54. Expression is the strongest known predictor of omega in any organism,
# so a pooled trend across these bins is confounded by construction. Re-form the
# SAME summed-count ratio within expression tertiles and within degree tertiles:
# if the trend is about conservation it survives; if it is about expression it
# collapses. This is the same discipline test 2 applies to the outgroup omega,
# and it is not optional.
tertile <- function(x) cut(rank(x, ties.method = "first"),
                           breaks = quantile(rank(x, ties.method = "first"),
                                             c(0, 1/3, 2/3, 1)),
                           include.lowest = TRUE, labels = c("low", "mid", "high"))
p[, tpm_tertile := tertile(sc_mean_log_tpm)]
p[, deg_tertile := tertile(sc_degree)]

strat_pair <- rbindlist(list(
  p[!is.na(tpm_tertile), .(stratum = "expression", level = as.character(tpm_tertile),
                           n_genes = .N, omega = binned_omega(.SD)),
    by = .(cons_bin, tpm_tertile)][, tpm_tertile := NULL],
  p[!is.na(deg_tertile), .(stratum = "degree", level = as.character(deg_tertile),
                           n_genes = .N, omega = binned_omega(.SD)),
    by = .(cons_bin, deg_tertile)][, deg_tertile := NULL]))
setcolorder(strat_pair, c("stratum", "level", "cons_bin", "n_genes", "omega"))
say("== 5b. the same ratio, within expression and degree tertiles")
print(dcast(strat_pair[stratum == "expression"], cons_bin ~ level, value.var = "omega"))
fwrite(strat_pair, file.path(OUTDIR, "dnds_saccharum_pair_stratified.tsv"), sep = "\t")

# Does the pooled trend survive stratification? Compare the drop from the least
# to the most conserved bin, pooled vs within each stratum level.
drop_of <- function(dt) {
  a <- dt[cons_bin == levels(cons_bin)[1], omega]
  b <- dt[cons_bin == levels(cons_bin)[nlevels(cons_bin)], omega]
  if (length(a) && length(b) && is.finite(a) && is.finite(b)) b / a else NA_real_
}
add("5b pooled", "omega_ratio_mostcons_over_leastcons", drop_of(pairbins), nrow(p),
    "the pooled trend, which the confounders make uninterpretable")
for (lv in c("low", "mid", "high"))
  add("5b within expression tertile", sprintf("omega_ratio_mostcons_over_leastcons_%s", lv),
      drop_of(strat_pair[stratum == "expression" & level == lv]), p[tpm_tertile == lv, .N],
      "if these do not reproduce the pooled ratio, the trend was expression")

conf <- p[, .(n = .N, med_degree = as.numeric(median(sc_degree)),
              med_log_tpm = as.numeric(median(sc_mean_log_tpm, na.rm = TRUE)),
              med_gc3 = as.numeric(median(gc3)),
              med_omega_outgroup = as.numeric(median(sc_sb_omega, na.rm = TRUE))),
          by = cons_bin][order(cons_bin)]
say("== 5b. how the conservation bins differ in the confounders")
print(conf)
fwrite(conf, file.path(OUTDIR, "dnds_cons_bin_confounders.tsv"), sep = "\t")

# --- 6. effect size, which is what should be reported ----------------------
q <- quantile(a$sc_frac_conserved, c(0.1, 0.9))
top <- a[sc_frac_conserved >= q[2], median(sc_sb_omega)]
bot <- a[sc_frac_conserved <= q[1], median(sc_sb_omega)]
add("6 effect size", "median_omega_top_decile_conservation", top, a[sc_frac_conserved >= q[2], .N], "")
add("6 effect size", "median_omega_bottom_decile_conservation", bot, a[sc_frac_conserved <= q[1], .N], "")
add("6 effect size", "ratio_top_over_bottom", top / bot, nrow(a),
    "below 1 = more conserved neighbourhoods are under stronger constraint")

out <- rbindlist(res)
fwrite(out, file.path(OUTDIR, "dnds_tests.tsv"), sep = "\t")
say("")
print(out)
say("")
say("== wrote dnds_tests.tsv, dnds_by_degree_decile.tsv, dnds_saccharum_pair_binned.tsv")
