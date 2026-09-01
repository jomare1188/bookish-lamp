#!/usr/bin/env Rscript
# ============================================================================
# 15_centrality_test.r -- does local clustering say anything degree did not?
#
# THE SCREEN COMES FIRST, AND IT CAN END THE ANALYSIS. The worry going in was
# that clustering would be degree wearing a different label -- in many networks
# C(k) ~ 1/k, and then a clustering-vs-constraint correlation just restates the
# degree result. Measured here it is not: rho(clustering, degree) = +0.48, and
# 73% of clustering's spread survives within degree deciles. The screen passes,
# so clustering is a genuinely second description of position. The reported
# result is still always the PARTIAL effect after degree.
#
# Two readouts, because omega is distorted here (it assumes dN scales 1:1 with
# dS; the fitted slope is 0.511 alone, 0.725 with covariates, and GC3 inflates dS
# threefold):
#   log(omega)        the published readout, kept for comparability
#   constraint_score  residual of log(dN) after log(dS), GC3, expression, length
#
# Their nuisance structure differs, so their models must differ. constraint_score
# ALREADY has composition, expression, length and the synonymous rate regressed
# out; putting them back in would be double-counting. omega has not.
# ============================================================================
suppressPackageStartupMessages({library(data.table)})

env <- function(k, d = NA) { v <- Sys.getenv(k); if (nzchar(v)) v else d }
OUTDIR <- env("OUTDIR")
say <- function(...) cat(sprintf(...), "\n", sep = "")

res <- list()
add <- function(test, readout, statistic, value, n, note = "", p = NA_real_,
                dR2 = NA_real_, std_coef = NA_real_)
  res[[length(res) + 1]] <<- data.table(test, readout, statistic, value, n, p,
                                        dR2, std_coef, note)

d <- fread(file.path(OUTDIR, "dnds_gene_table.tsv"))
cen <- fread(file.path(OUTDIR, "centrality_sugarcane.tsv"))
setnames(cen, "gene", "sugarcane_gene")
d <- merge(d, cen[, .(sugarcane_gene, clustering, coreness)],
           by = "sugarcane_gene", all.x = TRUE)

a <- d[in_network == TRUE & ds_ok_sc == TRUE & omega_outlier == FALSE &
         !is.na(sc_frac_conserved) & !is.na(sc_mean_log_tpm) &
         sc_sb_omega > 0 & !is.na(clustering)]
a[, `:=`(lomega = log(sc_sb_omega), ldeg = log(sc_degree))]
say("== analysis set: %d genes with clustering and a constraint value", nrow(a))

# ---------------------------------------------------------------- screen ----
rho_cd <- cor(a$clustering, a$sc_degree, method = "spearman")
add("0 redundancy", "-", "spearman_clustering_vs_degree", rho_cd, nrow(a),
    "if |rho| > 0.95 clustering is a relabelling of degree")
say("== redundancy screen")
say("   spearman(clustering, degree) = %+.3f", rho_cd)

a[, deg_decile := cut(rank(sc_degree, ties.method = "first"),
                      breaks = quantile(rank(sc_degree, ties.method = "first"),
                                        probs = seq(0, 1, 0.1)),
                      include.lowest = TRUE, labels = 1:10)]
ck <- a[, .(n = .N, median_degree = as.numeric(median(sc_degree)),
            median_clustering = as.numeric(median(clustering)),
            sd_clustering = as.numeric(sd(clustering)),
            iqr_clustering = as.numeric(IQR(clustering))),
        by = deg_decile][order(deg_decile)]
# NOTE, measured rather than assumed: C(k) RISES here, from 0.00 at degree 1 to
# 0.81 at degree 2,637 (turning over slightly to 0.71 at the top). That is the
# opposite of the C(k) ~ 1/k of preferential-attachment networks, and it is why
# clustering is NOT redundant with degree in this graph (rho = +0.48). A
# thresholded correlation network is not a growth network: a hub sits inside a
# dense co-expressed module, so its neighbours are correlated with each other
# too and its neighbourhood is near-complete. Corroborated by the two global
# statistics -- mean local clustering 0.464 against global transitivity 0.690;
# transitivity is triangle-weighted and so dominated by exactly those hubs.
say("== C(k): clustering by degree decile (measured; see the note above)")
print(ck)
fwrite(ck, file.path(OUTDIR, "centrality_ck_curve.tsv"), sep = "\t")

# Is there anything left to test once degree is fixed? If clustering has no
# spread within a degree decile, no partial test can work.
add("0 redundancy", "-", "median_within_decile_sd_clustering",
    median(ck$sd_clustering), nrow(a),
    sprintf("overall sd = %.4f; the ratio is the usable signal", sd(a$clustering)))
say("   clustering sd overall %.4f, median within degree decile %.4f (%.0f%% retained)",
    sd(a$clustering), median(ck$sd_clustering),
    100 * median(ck$sd_clustering) / sd(a$clustering))

GATE_SHUT <- abs(rho_cd) > 0.95 || median(ck$sd_clustering) / sd(a$clustering) < 0.1
if (GATE_SHUT)
  say("   ** SCREEN FAILED: clustering is not an independent description of\n      position in this network. Results below are reported but must NOT be\n      presented as evidence separate from degree. **")

# ----------------------------------------------------------------- tests ----
READOUTS <- list(
  omega = list(col = "lomega",
               form = "%s ~ sc_frac_conserved + ldeg + clustering + sc_mean_log_tpm + aln_codons + gc3"),
  # constraint_score is already a residual on dS, GC3, expression and length --
  # re-adding them would be double-counting, so its model carries only the
  # network terms.
  constraint = list(col = "constraint_score",
                    form = "%s ~ sc_frac_conserved + ldeg + clustering"))

for (rn in names(READOUTS)) {
  spec <- READOUTS[[rn]]
  b <- a[!is.na(get(spec$col))]
  if (nrow(b) < 100) { say("!! %s: too few rows, skipped", rn); next }
  say("")
  say("================ readout: %s (n = %d) ================", rn, nrow(b))

  # 3. marginal -- confounded by construction, reported for completeness
  ct <- suppressWarnings(cor.test(b$clustering, b[[spec$col]], method = "spearman"))
  add("3 marginal", rn, "spearman_clustering", unname(ct$estimate), nrow(b),
      "confounded with degree; not the result", ct$p.value)

  # 4. partial -- THE result
  m <- lm(as.formula(sprintf(spec$form, spec$col)), data = b)
  cf <- summary(m)$coefficients
  ci <- confint(m)
  R2 <- summary(m)$r.squared
  for (v in c("clustering", "ldeg", "sc_frac_conserved")) {
    if (!v %in% rownames(cf)) next
    drop <- R2 - summary(update(m, as.formula(paste(". ~ . -", v))))$r.squared
    sc <- cf[v, "Estimate"] * sd(b[[v]]) / sd(b[[spec$col]])
    add("4 partial", rn, v, cf[v, "Estimate"], nrow(b),
        sprintf("CI [%.4f, %.4f]", ci[v, 1], ci[v, 2]), cf[v, "Pr(>|t|)"],
        dR2 = drop, std_coef = sc)
  }
  add("4 partial", rn, "model_r2", R2, nrow(b), "")
  say("   partial coefficients -- dR2 is the comparable quantity across terms:")
  print(rbindlist(res)[test == "4 partial" & readout == rn,
                       .(statistic, coef = round(value, 4), dR2 = signif(dR2, 3),
                         std_coef = round(std_coef, 3), p = signif(p, 3))])

  # 5. within degree deciles -- is the sign even stable?
  strat <- b[, {
    ct <- if (.N >= 10) suppressWarnings(cor.test(clustering, get(spec$col), method = "spearman"))
          else list(estimate = NA_real_, p.value = NA_real_)
    .(n = .N, rho = as.numeric(unname(ct$estimate)), p = as.numeric(ct$p.value))
  }, by = deg_decile][order(deg_decile)]
  fwrite(strat, file.path(OUTDIR, sprintf("centrality_by_degree_decile_%s.tsv", rn)), sep = "\t")
  add("5 stratified", rn, "deciles_with_negative_rho", sum(strat$rho < 0, na.rm = TRUE), nrow(b),
      "out of 10; consistency of sign is the evidence, not any single bin")

  # 6. effect size, overall and inside one degree decile
  q <- quantile(b$clustering, c(0.1, 0.9))
  hi <- b[clustering >= q[2], median(get(spec$col))]
  lo <- b[clustering <= q[1], median(get(spec$col))]
  add("6 effect size", rn, "median_top_minus_bottom_clustering_decile", hi - lo, nrow(b),
      sprintf("top %.4f vs bottom %.4f", hi, lo))
  mid <- b[deg_decile == 5]
  if (nrow(mid) > 50) {
    qm <- quantile(mid$clustering, c(0.1, 0.9))
    add("6 effect size", rn, "same_within_degree_decile_5",
        mid[clustering >= qm[2], median(get(spec$col))] -
          mid[clustering <= qm[1], median(get(spec$col))], nrow(mid),
        "degree held roughly fixed")
  }
}

# 7. did adding clustering disturb what was already concluded?
m_before <- lm(lomega ~ sc_frac_conserved + ldeg + sc_mean_log_tpm + aln_codons + gc3, data = a)
m_after  <- lm(lomega ~ sc_frac_conserved + ldeg + clustering + sc_mean_log_tpm + aln_codons + gc3, data = a)
for (v in c("sc_frac_conserved", "ldeg")) {
  d1 <- summary(m_before)$r.squared -
    summary(update(m_before, as.formula(paste(". ~ . -", v))))$r.squared
  d2 <- summary(m_after)$r.squared -
    summary(update(m_after, as.formula(paste(". ~ . -", v))))$r.squared
  add("7 stability", "omega", sprintf("dR2_%s_before_vs_after", v), d2 - d1, nrow(a),
      sprintf("before %.5f, after %.5f -- should barely move", d1, d2))
}

out <- rbindlist(res)
# Holm over the family of constraint tests actually claimed (marginal + partial
# clustering, both readouts). n = 12k makes significance cheap, so effect sizes
# lead and this is a guard, not the headline.
fam <- out[(test == "3 marginal") | (test == "4 partial" & statistic == "clustering")]
out[, p_holm := NA_real_]
if (nrow(fam) > 0) {
  idx <- which((out$test == "3 marginal") |
                 (out$test == "4 partial" & out$statistic == "clustering"))
  out[idx, p_holm := p.adjust(p, method = "holm")]
}
out[, screen_passed := !GATE_SHUT]
fwrite(out, file.path(OUTDIR, "centrality_tests.tsv"), sep = "\t")
say("")
print(out[, .(test, readout, statistic, value = round(value, 5),
              dR2 = signif(dR2, 3), p = signif(p, 3), p_holm = signif(p_holm, 3))])
say("")
say("== wrote centrality_tests.tsv, centrality_ck_curve.tsv, centrality_by_degree_decile_*.tsv")
if (GATE_SHUT) say("== SCREEN FAILED -- see the redundancy rows before quoting anything above")
