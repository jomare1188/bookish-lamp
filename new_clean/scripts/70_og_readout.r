#!/usr/bin/env Rscript
# =============================================================================
# 70_og_readout.r -- the keystone result re-tested against quantification
# uncertainty, and the three ICCs side by side.
#
# WHY THIS EXISTS. results/dnds/decoupling_tests.tsv reports this project's most
# consequential negative: gene copies at CDS identity >= 0.99, in ONE genome and
# ONE set of samples, differ 8.2-fold in degree (sugarcane) and sit in different
# MCL modules 77.6% of the time. Everything that follows from it -- that wiring is
# not a conserved property, that the orthogroup is the better unit -- rests on
# that number not being a multi-mapping artefact.
#
# THE OLD CONTROL WAS TOO THIN TO SETTLE IT. 12_decoupling_test.r splits on
# `distinguishable`, a k-mer proxy scoped WITHIN the family, and only 300 of 8,562
# sugarcane and 714 of 53,264 purple near-identical pairs qualify. In those small
# subsets the effect shrinks from 8.4x to 6.3x and from 7.1x to 3.7x -- in purple
# it nearly halves, which is exactly the shape of an artefact and exactly too few
# pairs to be sure.
#
# WHAT REPLACES IT. 65's boot_r_within is the quantifier's own answer, genome-wide:
# strongly negative when the EM is trading reads between two copies. Joining it to
# the decoupling pairs ON THE GENE PAIR -- which is namespace-free, so the 2sp/3sp
# orthogroup mismatch cannot intrude -- lets the keystone be recomputed separately
# for pairs the quantifier CAN and CANNOT separate.
#
# THE PREDICTION THE ARTEFACT HYPOTHESIS MAKES IS FALSIFIABLE. If multi-mapping
# manufactures the decoupling, then the read-stealing pairs carry the inflated
# divergence and the separable pairs carry less. Measured 2026-10-01 it is the
# other way round in both species, so the artefact hypothesis is refuted rather
# than merely unsupported. This stage writes that comparison rather than asserting
# it, so a later rebuild can overturn it.
#
# THE THREE ICCs. One statistic (icc_oneway, common.R), one grouping (the
# orthogroup), three quantities -- and they do not agree, which is the point:
# omega is a property of the orthogroup, degree is not.
#
# RUN: through run.sh -> ./run.sh ogreadout
# =============================================================================

suppressMessages(library(data.table))
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE),
                                                 value = TRUE)[1])), "lib", "common.R"))

STUDIES    <- env_list("CLEAN_STUDIES", c("sugarcane", "purple"))
RESULTS    <- env_req("CLEAN_RESULTS")
DNDS_DIR   <- env_req("CLEAN_DNDS_DIR")
OUT_DIR    <- env_req("CLEAN_OUT_DIR")
STEAL_THR  <- env_num("CLEAN_STEAL_THR", -0.5)

banner("the keystone, re-tested against quantification uncertainty")
ensure_dir(OUT_DIR)

pair_key <- function(a, b) fifelse(a < b, paste(a, b, sep = "|"), paste(b, a, sep = "|"))
ID_BINS <- c(-Inf, 0.95, 0.98, 0.99, 0.999, Inf)

keystone <- list(); idcurve <- list(); iccs <- list()

for (s in STUDIES) {
  dp_f <- file.path(DNDS_DIR, sprintf("decoupling_pairs_%s.tsv", s))
  rs_f <- file.path(RESULTS, s, "infreps", sprintf("%s_readsteal_pairs.tsv", s))
  for (f in c(dp_f, rs_f))
    if (!file.exists(f))
      stop(basename(f), " missing\n  ", if (identical(f, rs_f))
           paste0("run  ./run.sh infreps ", s) else "run the dN/dS stage",
           "  first", call. = FALSE)

  dp <- fread(dp_f)
  rs <- fread(rs_f, select = c("gene_a", "gene_b", "boot_r_within",
                               "n_samples_usable", "verdict"))
  dp[, k := pair_key(gene_a, gene_b)]
  rs[, k := pair_key(gene_a, gene_b)]
  m <- merge(dp, rs[, .(k, boot_r_within, n_samples_usable, verdict)], by = "k")
  say(sprintf("%-10s %s decoupling pairs | %s carry a read-stealing verdict (%.1f%%)",
              s, fmt_n(nrow(dp)), fmt_n(nrow(m)), 100 * nrow(m) / nrow(dp)))
  # The unjoined majority are pairs whose two genes sit in DIFFERENT 2sp
  # orthogroups: the families are 3sp groups, which merge 2sp ones. Near-identical
  # pairs -- the ones that matter here -- are ~90% within one 2sp group, so the
  # join covers the question even though it does not cover the table.
  m[, near := pid_cds >= 0.99]
  m[, id_bin := cut(pid_cds, ID_BINS)]

  idcurve[[s]] <- m[, .(study = s, n = .N,
                        pct_read_stealing = round(100 * mean(verdict == "read_stealing"), 3),
                        pct_determined = round(100 * mean(verdict == "determined"), 3),
                        median_boot_r = round(median(boot_r_within, na.rm = TRUE), 4)),
                    by = id_bin][order(id_bin)]

  for (scope in c("near_identical", "all_pairs")) {
    mm <- if (scope == "near_identical") m[near == TRUE] else m
    keystone[[paste(s, scope)]] <- mm[is.finite(net_div), .(
      study = s, scope = scope, n = .N,
      median_fold_degree = round(2^median(net_div), 3),
      pct_different_module = round(100 * mean(diff_module, na.rm = TRUE), 2),
      median_pid_cds = round(median(pid_cds, na.rm = TRUE), 5)),
      by = verdict][order(verdict)]
  }

  # The three ICCs, on this species' orthogroups.
  og_f <- file.path(RESULTS, s, sprintf("og_uniformity_%s_summary.tsv", s))
  om_f <- file.path(DNDS_DIR, sprintf("percopy_omega_icc_%s.tsv", s))
  got <- data.table(study = s, quantity = character(), icc = numeric(),
                    source = character())[0]
  if (file.exists(om_f)) {
    o <- fread(om_f)
    got <- rbind(got, data.table(study = s, quantity = "omega_vs_sorghum",
                                 icc = as.numeric(o[metric == "icc", value]),
                                 source = basename(om_f)))
  }
  if (file.exists(og_f)) {
    u <- fread(og_f)
    for (q in c("nitrogen_r_partial", "log_degree")) {
      v <- u[metric == paste0("icc_", q), value]
      if (length(v))
        got <- rbind(got, data.table(study = s, quantity = q,
                                     icc = as.numeric(v), source = basename(og_f)))
    }
  }
  iccs[[s]] <- got
}

ks <- rbindlist(keystone)
ic <- rbindlist(idcurve)
ii <- rbindlist(iccs)

banner("read-stealing follows CDS identity, as it must if the statistic works")
for (s in STUDIES) {
  say(sprintf("%s:", s))
  for (i in seq_len(nrow(ic[study == s])))
    with(ic[study == s][i], say(sprintf(
      "  %-14s %9s pairs | %6.2f%% read-stealing | %5.2f%% determined | median boot_r %7.4f",
      as.character(id_bin), fmt_n(n), pct_read_stealing, pct_determined, median_boot_r)))
}

banner("THE KEYSTONE, split by whether the copies are separable at all")
for (s in STUDIES) {
  for (sc in c("near_identical", "all_pairs")) {
    r <- ks[study == s & scope == sc]
    if (!nrow(r)) next
    say(sprintf("%s, %s:", s, sc))
    for (i in seq_len(nrow(r)))
      with(r[i], say(sprintf("  %-14s %8s pairs | median degree fold %6.2fx | %5.1f%% different module",
                             verdict, fmt_n(n), median_fold_degree, pct_different_module)))
  }
}

# The falsifiable comparison, stated as a verdict rather than left to the reader.
banner("verdict")
for (s in STUDIES) {
  r <- ks[study == s & scope == "near_identical"]
  sep <- r[verdict == "separable", median_fold_degree]
  stl <- r[verdict == "read_stealing", median_fold_degree]
  if (!length(sep) || !length(stl)) next
  say(sprintf("%-10s separable %.2fx vs read-stealing %.2fx -> %s",
              s, sep, stl,
              if (sep > stl)
                "the decoupling is LARGER where copies are separable; the artefact hypothesis is refuted"
              else
                "the decoupling is larger where reads are traded; IT MAY BE AN ARTEFACT"))
}

banner("one statistic, one grouping, three quantities")
for (i in seq_len(nrow(ii)))
  with(ii[i], say(sprintf("  %-10s %-20s ICC %.4f", study, quantity, icc)))

write_tsv(ks, file.path(OUT_DIR, "keystone_by_separability.tsv"))
write_tsv(ic, file.path(OUT_DIR, "readsteal_by_cds_identity.tsv"))
write_tsv(ii, file.path(OUT_DIR, "icc_three_quantities.tsv"))
say("done")
