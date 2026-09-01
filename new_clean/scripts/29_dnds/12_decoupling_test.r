#!/usr/bin/env Rscript
# ============================================================================
# 12_decoupling_test.r -- do near-identical gene copies hold the same network job?
#
# The claim this tests: sugarcane's polyploid copies are ~98% identical in
# sequence yet sit in completely different places in the co-expression network.
# If that survives its controls it is the mirror image of the 1:1 null -- there,
# sequence constraint was blind to network conservation; here, network position
# moves while sequence stands still.
#
# TWO CONTROLS DECIDE WHETHER IT IS REAL, and both are capable of killing it:
#
#   quantification  copies sharing nearly all their 31-mers cannot be told apart
#                   by salmon, so their "different expression" -- and hence their
#                   different degree -- may be EM arithmetic rather than biology.
#                   19.6% of sugarcane copies and 45.5% of purple copies have
#                   ZERO unique 31-mers. Test 4 splits the headline on this.
#
#   expression      a copy that is barely expressed cannot acquire high degree in
#                   a correlation network. So the question is never "do copies
#                   differ in degree" but "do they differ MORE than their
#                   expression difference accounts for". Test 3, same discipline
#                   step 06 uses.
# ============================================================================
suppressPackageStartupMessages({library(data.table)})

env <- function(k, d = NA) { v <- Sys.getenv(k); if (nzchar(v)) v else d }
OUTDIR <- env("OUTDIR")
FLOOR  <- as.numeric(env("MIN_FRAC_UNIQUE"))
NEAR   <- as.numeric(env("PID_NEAR_IDENTICAL"))
NNULL  <- as.integer(env("NULL_PAIRS"))
say <- function(...) cat(sprintf(...), "\n", sep = "")

mean_log_tpm <- function(path, strip) {
  tpm <- fread(path); id <- tpm[[1]]
  if (strip) id <- sub("\\.v[0-9.]+$", "", id)
  data.table(gene = id, tpm = rowMeans(log2(as.matrix(tpm[, -c(1, 2), with = FALSE]) + 1)))
}

res <- list()
add <- function(study, test, statistic, value, n, note = "")
  res[[length(res) + 1]] <<- data.table(study, test, statistic, value, n, note)

for (study in c("sugarcane", "purple")) {
  pf <- file.path(OUTDIR, sprintf("family_pairs_%s.tsv", study))
  if (!file.exists(pf) || file.size(pf) == 0) { say("!! %s: no pairs, skipped", study); next }
  say("================ %s ================", study)
  p  <- fread(pf)
  ku <- fread(file.path(OUTDIR, sprintf("kmer_uniqueness_%s.tsv", study)))
  nm <- fread(env(sprintf("NODE_METRICS_%s", study)), select = c("gene", "degree"))
  cd <- fread(file.path(OUTDIR, sprintf("conserved_degree_%s.tsv", study)),
              select = c("gene", "frac_conserved"))
  md <- fread(env(sprintf("MODULES_%s", study)), select = c("gene", "module_name"))
  ex <- mean_log_tpm(env(sprintf("TPM_%s", study)), study == "sugarcane")
  ex <- unique(ex, by = "gene")
  fam <- fread(file.path(OUTDIR, sprintf("families_%s.tsv", study)),
               select = c("orthogroup", "n_copies"))

  attach_side <- function(dt, col, suffix) {
    for (tb in list(list(nm, "degree"), list(cd, "frac_conserved"),
                    list(md, "module_name"), list(ex, "tpm"),
                    list(ku[, .(gene, frac_unique)], "frac_unique"))) {
      s <- copy(tb[[1]]); setnames(s, "gene", col)
      s <- unique(s, by = col)
      dt <- merge(dt, s, by = col, all.x = TRUE)
      setnames(dt, tb[[2]], paste0(tb[[2]], suffix))
    }
    dt
  }
  p <- attach_side(p, "gene_a", "_a")
  p <- attach_side(p, "gene_b", "_b")
  p <- merge(p, fam, by = "orthogroup", all.x = TRUE)

  # both copies must be IN the network for a network comparison to mean anything
  p <- p[!is.na(degree_a) & !is.na(degree_b) & degree_a > 0 & degree_b > 0]
  p[, net_div     := abs(log2(degree_a / degree_b))]
  p[, cons_div    := abs(frac_conserved_a - frac_conserved_b)]
  p[, expr_div    := abs(tpm_a - tpm_b)]
  p[, frac_unique_min := pmin(frac_unique_a, frac_unique_b, na.rm = FALSE)]
  p[, diff_module := module_name_a != module_name_b &
       module_name_a != "Unassigned" & module_name_b != "Unassigned"]
  p[, distinguishable := !is.na(frac_unique_min) & frac_unique_min >= FLOOR]

  fwrite(p, file.path(OUTDIR, sprintf("decoupling_pairs_%s.tsv", study)), sep = "\t")
  say("pairs with both copies in the network: %d (%d homeolog, %d tandem, %d dispersed)",
      nrow(p), p[pair_class == "homeolog", .N], p[pair_class == "tandem", .N],
      p[pair_class == "dispersed", .N])

  h <- p[pair_class == "homeolog"]

  # --- 1. the headline, as a magnitude ------------------------------------
  ni <- h[pid_cds >= NEAR]
  add(study, "1 headline", "n_near_identical_pairs", nrow(ni), nrow(h),
      sprintf("CDS identity >= %.2f", NEAR))
  if (nrow(ni) > 0) {
    add(study, "1 headline", "median_fold_degree_difference", 2^median(ni$net_div), nrow(ni), "")
    add(study, "1 headline", "pct_pairs_over_10x", 100 * mean(ni$net_div > log2(10)), nrow(ni), "")
    add(study, "1 headline", "pct_in_different_modules", 100 * mean(ni$diff_module, na.rm = TRUE),
        sum(!is.na(ni$diff_module)), "")
  }

  # --- 2. does identity predict network similarity at all? ----------------
  for (v in c("net_div", "cons_div", "expr_div")) {
    ct <- suppressWarnings(cor.test(h$pid_cds, h[[v]], method = "spearman"))
    add(study, "2 identity vs divergence", sprintf("spearman_rho_%s", v),
        unname(ct$estimate), nrow(h), sprintf("p = %.3g", ct$p.value))
  }

  # --- 3. the expression control, which decides it ------------------------
  m <- lm(net_div ~ expr_div + pid_cds + n_copies + frac_unique_min, data = h)
  cf <- summary(m)$coefficients; ci <- confint(m)
  for (v in rownames(cf)[-1])
    add(study, "3 expression control", v, cf[v, "Estimate"], nrow(m$model),
        sprintf("CI [%.3f, %.3f], p = %.3g", ci[v, 1], ci[v, 2], cf[v, "Pr(>|t|)"]))
  m0 <- lm(net_div ~ expr_div, data = h)
  add(study, "3 expression control", "adj_r2_expression_alone", summary(m0)$adj.r.squared,
      nrow(m0$model), "how much of network divergence is just expression divergence")
  add(study, "3 expression control", "adj_r2_full", summary(m)$adj.r.squared, nrow(m$model), "")

  h[, expr_stratum := cut(expr_div, breaks = quantile(expr_div, c(0, 1/3, 2/3, 1), na.rm = TRUE),
                          include.lowest = TRUE, labels = c("low", "mid", "high"))]
  h[, id_bin := cut(pid_cds, breaks = c(-Inf, 0.95, 0.98, 0.99, 0.999, Inf),
                    labels = c("<0.95", "0.95-0.98", "0.98-0.99", "0.99-0.999", ">=0.999"))]
  strat <- h[!is.na(expr_stratum), .(n = .N, median_fold = 2^median(net_div),
                                     pct_diff_module = 100 * mean(diff_module, na.rm = TRUE)),
             by = .(expr_stratum, id_bin)][order(expr_stratum, id_bin)]
  say("== 3. median fold degree difference, by identity bin within expression strata")
  print(dcast(strat, id_bin ~ expr_stratum, value.var = "median_fold"))
  fwrite(strat, file.path(OUTDIR, sprintf("decoupling_by_identity_bin_%s.tsv", study)), sep = "\t")

  # --- 4. the mapping-artefact control ------------------------------------
  mc <- h[, .(n = .N, median_fold = 2^median(net_div),
              pct_diff_module = 100 * mean(diff_module, na.rm = TRUE),
              median_expr_div = median(expr_div, na.rm = TRUE)),
          by = .(distinguishable, near_identical = pid_cds >= NEAR)][order(-distinguishable, near_identical)]
  say("== 4. split on whether salmon can tell the copies apart (frac_unique >= %.2f)", FLOOR)
  print(mc)
  fwrite(mc, file.path(OUTDIR, sprintf("decoupling_mapping_control_%s.tsv", study)), sep = "\t")
  for (d in c(TRUE, FALSE)) {
    s <- h[distinguishable == d & pid_cds >= NEAR]
    if (nrow(s) > 10)
      add(study, "4 mapping control",
          sprintf("median_fold_degree_diff_%s", ifelse(d, "distinguishable", "ambiguous")),
          2^median(s$net_div), nrow(s), "near-identical pairs only")
  }

  # --- 5. the null: what does 'far apart' even look like? -----------------
  # Random gene pairs matched on expression decile, so the comparison is against
  # unrelated genes that are equally detectable -- not against nothing.
  pool <- merge(nm[degree > 0], ex, by = "gene")
  pool <- pool[!is.na(tpm)]
  pool[, edec := cut(rank(tpm, ties.method = "first"),
                     quantile(rank(tpm, ties.method = "first"), seq(0, 1, 0.1)),
                     include.lowest = TRUE, labels = FALSE)]
  set.seed(1)
  bydec <- split(pool$degree, pool$edec)
  hh <- h[!is.na(tpm_a) & !is.na(tpm_b)]
  idx <- sample.int(nrow(hh), min(NNULL, nrow(hh)), replace = TRUE)
  ra <- hh[idx]
  ea <- pool[match(ra$gene_a, gene), edec]; eb <- pool[match(ra$gene_b, gene), edec]
  drawn <- function(e) vapply(e, function(i)
    if (is.na(i)) NA_real_ else bydec[[as.character(i)]][sample.int(length(bydec[[as.character(i)]]), 1)],
    numeric(1))
  drawn_gene <- function(e) vapply(e, function(i)
    if (is.na(i)) NA_character_ else {
      g <- byg[[as.character(i)]]; g[sample.int(length(g), 1)] },
    character(1))
  byg <- split(pool$gene, pool$edec)
  nd <- abs(log2(drawn(ea) / drawn(eb)))
  nd <- nd[is.finite(nd)]

  # The module baseline matters as much as the degree one: with ~10,000 MCL
  # modules, two unrelated genes land in different modules almost always, so
  # "77% of copy pairs are in different modules" is meaningless until it is read
  # against that. Draw the null pairs' module-difference rate the same way.
  mmap <- setNames(md$module_name, md$gene)
  ga <- drawn_gene(ea); gb <- drawn_gene(eb)
  ma <- mmap[ga]; mb <- mmap[gb]
  keep <- !is.na(ma) & !is.na(mb) & ma != "Unassigned" & mb != "Unassigned"
  add(study, "5 null", "pct_diff_module_random_pairs", 100 * mean((ma != mb)[keep]),
      sum(keep), "chance rate; compare with the headline's pct_in_different_modules")

  add(study, "5 null", "median_fold_real_homeolog_pairs", 2^median(hh$net_div), nrow(hh), "")
  add(study, "5 null", "median_fold_expression_matched_random_pairs", 2^median(nd), length(nd),
      "if random pairs are NOT more divergent than copies, the matching is broken")
  add(study, "5 null", "ratio_random_over_real", 2^median(nd) / 2^median(hh$net_div), length(nd), "")

  # --- 6. contrasts --------------------------------------------------------
  cls <- p[, .(n = .N, median_fold = 2^median(net_div),
               pct_diff_module = 100 * mean(diff_module, na.rm = TRUE),
               median_pid_cds = median(pid_cds)), by = pair_class][order(-n)]
  say("== 6. by pair class"); print(cls)
  fwrite(cls, file.path(OUTDIR, sprintf("decoupling_by_class_%s.tsv", study)), sep = "\t")
}

out <- rbindlist(res)
fwrite(out, file.path(OUTDIR, "decoupling_tests.tsv"), sep = "\t")
say(""); print(out)
