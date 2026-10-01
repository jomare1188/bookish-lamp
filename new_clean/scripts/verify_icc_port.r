#!/usr/bin/env Rscript
# =============================================================================
# verify_icc_port.r -- prove common.R's icc_oneway() against a published answer.
#
# WHY THIS EXISTS. icc_oneway() is a port of 29_dnds/11_percopy_omega.py:117-130,
# and a port is worth nothing unless it reproduces the original. That script
# already wrote its answer to results/dnds/percopy_omega_icc_<study>.tsv on the
# per-copy omega data, so the port has a fixture with a known value rather than a
# self-consistency check:
#
#     sugarcane  0.940211        purple  0.938674
#
# The filtering has to match too, or the port would be scored on different
# numbers: dS within [DNDS_DS_MIN, DNDS_DS_MAX] and 0 < omega <= DNDS_OMEGA_FLAG,
# grouped by orthogroup (11_percopy_omega.py:169-172).
#
# A MISMATCH MEANS THE PORT IS WRONG, not that the data moved. Tolerance is PER
# METRIC, because the python wrote them at different precisions (11_percopy_omega.py
# :186-190): icc and the mean squares with %.6f, mean_family_size with %.4f, the
# counts as integers. One blanket 1e-6 fails on mean_family_size for that reason
# alone, which is a defect in the check and not in the port.
#
# RUN: $RSCRIPT_NET scripts/verify_icc_port.r
# =============================================================================

suppressMessages(library(data.table))
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE),
                                                 value = TRUE)[1])), "lib", "common.R"))

DNDS_DIR   <- env_opt("CLEAN_DNDS_DIR",
                      "/dados04/jorge/comparative_saccharum/new_clean/results/dnds")
DS_MIN     <- env_num("DNDS_DS_MIN", 0.01)
DS_MAX     <- env_num("DNDS_DS_MAX", 2.0)
OMEGA_FLAG <- env_num("DNDS_OMEGA_FLAG", 10)
TOL <- c(icc = 1e-6, ms_between = 1e-6, ms_within = 1e-6,
         mean_family_size = 1e-4, n_families = 0, n_copies = 0)

banner("icc_oneway(): R port vs the published python answer")
fail <- 0L

for (study in c("sugarcane", "purple")) {
  om_f  <- file.path(DNDS_DIR, sprintf("percopy_omega_%s.tsv", study))
  fix_f <- file.path(DNDS_DIR, sprintf("percopy_omega_icc_%s.tsv", study))
  if (!file.exists(om_f) || !file.exists(fix_f)) {
    say(sprintf("%-10s SKIP -- %s not on disk", study,
                basename(if (file.exists(om_f)) fix_f else om_f)))
    next
  }

  om <- fread(om_f)
  om <- om[dS >= DS_MIN & dS <= DS_MAX & omega > 0 & omega <= OMEGA_FLAG]
  got <- icc_oneway(om$omega, om$orthogroup)
  if (is.null(got)) {
    say(sprintf("%-10s FAIL -- port returned NULL on %s rows", study, fmt_n(nrow(om))))
    fail <- fail + 1L
    next
  }

  fix <- fread(fix_f)
  want <- setNames(fix$value, fix$metric)
  num <- function(k) suppressWarnings(as.numeric(want[[k]]))

  checks <- list(
    icc              = c(got$icc,             num("icc")),
    ms_between       = c(got$ms_between,      num("ms_between")),
    ms_within        = c(got$ms_within,       num("ms_within")),
    mean_family_size = c(got$mean_group_size, num("mean_family_size")),
    n_families       = c(got$n_groups,        num("n_families")),
    n_copies         = c(got$n_values,        num("n_copies")))

  bad <- names(which(vapply(names(checks), function(nm) {
    p <- checks[[nm]]
    !is.finite(p[1]) || !is.finite(p[2]) || abs(p[1] - p[2]) > TOL[[nm]]
  }, TRUE)))

  say(sprintf("%-10s icc %.6f  (published %.6f)  delta %.2g",
              study, got$icc, num("icc"), abs(got$icc - num("icc"))))
  say(sprintf("           ms_b %.6f / %.6f | ms_w %.6f / %.6f | k0 %.4f / %.4f | %s fam / %s copies",
              got$ms_between, num("ms_between"), got$ms_within, num("ms_within"),
              got$mean_group_size, num("mean_family_size"),
              fmt_n(got$n_groups), fmt_n(got$n_values)))
  if (length(bad)) {
    say(sprintf("           FAIL on: %s", paste(bad, collapse = ", ")))
    fail <- fail + 1L
  } else {
    say("           PASS -- every component matches to its written precision")
  }
}

if (fail > 0L)
  stop(fail, " study/studies disagree with the published ICC -- the port is wrong",
       call. = FALSE)
say("icc_oneway() reproduces 29_dnds/11_percopy_omega.py to its written precision")
