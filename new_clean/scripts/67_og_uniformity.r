#!/usr/bin/env Rscript
# =============================================================================
# 67_og_uniformity.r -- is an orthogroup one signal or several? And when its
# copies disagree, is that biology or is it the quantifier?
#
# WHY THIS EXISTS. This project's keystone negative is that gene copies at CDS
# identity >= 0.99, in ONE genome and ONE set of samples, differ 8.2-fold in
# degree and sit in different MCL modules 77.6% of the time
# (results/dnds/decoupling_tests.tsv). If network position is not a conserved
# property of a gene at that distance, then cross-species conservation of wiring
# was never findable at the gene level -- which is the argument for moving to the
# orthogroup. But before collapsing anything, two questions have to be answered
# per orthogroup:
#
#   does it behave as ONE unit?    -> expr_concordance, sign concordance
#   can its copies be told apart?  -> boot_r_within, from 65
#
# Crossing those two gives the only cell in this project that is a POSITIVE
# finding rather than a negative:
#
#                   separable copies          inseparable copies
#   uniform         redundant: collapse       collapse, mandatory
#                   loses nothing
#   divergent       REAL post-polyploidy      artefactual divergence --
#                   regulatory divergence     the set that was polluting
#                                             the earlier results
#
# WHAT IS PER GROUP AND WHAT IS GLOBAL. expr_concordance, sign concordance and
# the degree spread are per orthogroup and go in the main table. ICC is a
# BETWEEN-vs-WITHIN variance decomposition over all groups at once -- one number
# per species, not per group -- so it goes in the summary. Reporting an ICC per
# group would be a category error.
#
# THE CONCORDANCE HAS A NULL, and needs one: a 2-member group of unrelated genes
# does not score 0. The reference is SIZE-MATCHED random gene sets, drawn from the
# same expression matrix, so "uniform" means uniform above what group size alone
# buys.
#
# THE TRAIT IS LEFT IN THE RESIDUALS. Only the design blocks come out
# (blocked_residuals(), common.R). Removing nitrogen as well would measure
# agreement in the leftover noise, when the question is agreement in the biology.
#
# WHAT THIS IS NOT. Not a filter -- no gene or group is dropped from anything
# downstream on the strength of this table. The classification is a label to
# stratify by, and the plan's whole point is that BOTH sides of it get reported.
#
# RUN: through run.sh -> ./run.sh oguniformity <study>
# =============================================================================

suppressMessages(library(data.table))
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE),
                                                 value = TRUE)[1])), "lib", "common.R"))

STUDY      <- env_req("CLEAN_STUDY")
VST_PREFIX <- env_req("CLEAN_VST_PREFIX")
META_FILE  <- env_req("CLEAN_META")
ORTHO_FILE <- env_req("CLEAN_ORTHOGROUPS")
OG_COL     <- env_req("CLEAN_OG_SPECIES")
TRAIT_FILE <- env_req("CLEAN_GENETRAIT")
INFREP_DIR <- env_req("CLEAN_INFREP_DIR")
CROSSWALK  <- env_req("CLEAN_CROSSWALK")
OUT_DIR    <- env_req("CLEAN_OUT_DIR")
NODES      <- env_opt("CLEAN_NODE_METRICS")
MEMBERSHIP <- env_opt("CLEAN_MEMBERSHIP")
BLOCK      <- env_list("CLEAN_BLOCK", "genotype")
TRAIT      <- env_opt("CLEAN_SELECT_TRAIT", "treatment")
MAX_COPIES <- as.integer(env_num("CLEAN_MAX_OG_COPIES", 20))
NULL_SETS  <- as.integer(env_num("CLEAN_NULL_SETS", 200))
STEAL_THR  <- env_num("CLEAN_STEAL_THR", -0.5)
UNIFORM_Q  <- env_num("CLEAN_UNIFORM_QUANTILE", 0.95)
SEED       <- as.integer(env_num("CLEAN_SEED", 1))

banner(paste("orthogroup uniformity:", STUDY))
ensure_dir(OUT_DIR)
set.seed(SEED)

# --- expression, residualised on the design ---------------------------------
vst <- read_vst(VST_PREFIX)
meta <- fread(META_FILE, header = TRUE)
setnames(meta, tolower(names(meta)))
m <- meta[match(colnames(vst), sample)]
if (anyNA(m$sample))
  stop("VST samples missing from ", basename(META_FILE), call. = FALSE)
for (b in BLOCK)
  if (!b %in% names(m))
    stop("block column '", b, "' not in ", basename(META_FILE),
         "; have: ", paste(names(m), collapse = ", "), call. = FALSE)
blocks <- setNames(lapply(BLOCK, function(b) factor(m[[b]])), BLOCK)
say(sprintf("%s genes x %d samples | blocks: %s",
            fmt_n(nrow(vst)), ncol(vst), paste(BLOCK, collapse = " + ")))

RES <- blocked_residuals(vst, blocks)          # trait deliberately left in
RK  <- row_midranks(RES)                       # Spearman, as the purple trait test uses

# Row-standardise once. Then for a group of n unit-norm centred rows, the MEAN
# PAIRWISE correlation is (||sum of rows||^2 - n) / (n(n-1)) exactly -- one
# colSums per group instead of an n x n cor() per group, which is what makes
# 34,000-40,000 groups and a 200-replicate null affordable.
Z <- RK - rowMeans(RK)
nrm <- sqrt(rowSums(Z * Z))
bad <- !is.finite(nrm) | nrm < 1e-12
Z <- Z / nrm
Z[bad, ] <- 0
say(sprintf("%s genes have no residual variance and score as uncorrelated",
            fmt_n(sum(bad))))

mean_pairwise <- function(idx_list) {
  vapply(idx_list, function(ii) {
    n <- length(ii)
    if (n < 2L) return(NA_real_)
    s <- colSums(Z[ii, , drop = FALSE])
    (sum(s * s) - n) / (n * (n - 1))
  }, numeric(1))
}

# --- the orthogroup partition -----------------------------------------------
og <- as.data.table(read_orthogroups_tsv(ORTHO_FILE))
og[, Gene := strip_version(sub("\\.p[0-9]+$", "", Gene))]
og <- unique(og[Species == OG_COL & Gene %chin% rownames(vst), .(Orthogroup, Gene)])
gidx <- setNames(seq_len(nrow(vst)), rownames(vst))
members <- og[, .(idx = list(gidx[Gene]), genes = list(Gene), n = .N), by = Orthogroup]
say(sprintf("orthogroups with >=1 expressed gene: %s | multi-copy: %s",
            fmt_n(nrow(members)), fmt_n(members[n >= 2, .N])))

members[, expr_concordance := mean_pairwise(idx)]

# --- the size-matched null ---------------------------------------------------
sizes <- sort(unique(members[n >= 2 & n <= MAX_COPIES, n]))
say(sprintf("null: %d random sets per group size, sizes %d-%d",
            NULL_SETS, min(sizes), max(sizes)))
null_dt <- rbindlist(lapply(sizes, function(k) {
  draws <- replicate(NULL_SETS, sample.int(nrow(vst), k), simplify = FALSE)
  data.table(n = k, v = mean_pairwise(draws))
}))
null_q <- null_dt[, .(null_mean = mean(v, na.rm = TRUE),
                      null_cut = as.numeric(quantile(v, UNIFORM_Q, na.rm = TRUE))),
                  by = n]
members <- merge(members, null_q, by = "n", all.x = TRUE)
say(sprintf("  null mean at n=2: %.4f | %.0fth pct cut: %.4f",
            null_q[n == 2, null_mean], 100 * UNIFORM_Q, null_q[n == 2, null_cut]))

# --- nitrogen response per member -------------------------------------------
tr <- fread(TRAIT_FILE, select = c("gene", "r_partial", "padj", "responsive"))
tr[, gene := strip_version(gene)]
mm <- merge(og, tr, by.x = "Gene", by.y = "gene")
nres <- mm[, {
  s <- sign(r_partial[is.finite(r_partial)])
  .(n_tested = sum(is.finite(r_partial)),
    n_r_mean = mean(r_partial, na.rm = TRUE),
    n_r_sd   = if (sum(is.finite(r_partial)) >= 2) sd(r_partial, na.rm = TRUE) else NA_real_,
    sign_concordance = if (length(s)) max(sum(s > 0), sum(s < 0)) / length(s) else NA_real_,
    n_responsive = sum(responsive %in% c(TRUE, "TRUE")))
}, by = Orthogroup]
members <- merge(members, nres, by = "Orthogroup", all.x = TRUE)

# --- network position per member (where the gene is a node) -----------------
if (nzchar(NODES) && file.exists(NODES)) {
  nm <- fread(NODES, select = c("gene", "degree"))
  nm[, gene := strip_version(gene)]
  dd <- merge(og, nm, by.x = "Gene", by.y = "gene")
  dg <- dd[, .(n_nodes = .N,
               degree_median = as.numeric(median(degree)),
               degree_fold_range = if (.N >= 2 && min(degree) > 0)
                 max(degree) / min(degree) else NA_real_), by = Orthogroup]
  members <- merge(members, dg, by = "Orthogroup", all.x = TRUE)
}
if (nzchar(MEMBERSHIP) && file.exists(MEMBERSHIP)) {
  mb <- fread(MEMBERSHIP, select = c("gene", "module_name"))
  mb[, gene := strip_version(gene)]
  mj <- merge(og, mb, by.x = "Gene", by.y = "gene")
  md <- mj[module_name != "Unassigned",
           .(n_distinct_modules = uniqueN(module_name), n_in_modules = .N),
           by = Orthogroup]
  members <- merge(members, md, by = "Orthogroup", all.x = TRUE)
}

# --- the uncertainty side, from 65 ------------------------------------------
unc_f  <- file.path(INFREP_DIR, sprintf("%s_og_uncertainty.tsv", STUDY))
pair_f <- file.path(INFREP_DIR, sprintf("%s_readsteal_pairs.tsv", STUDY))
for (f in c(unc_f, pair_f))
  if (!file.exists(f))
    stop(basename(f), " missing\n  run  ./run.sh infreps ", STUDY, "  first",
         call. = FALSE)
unc <- fread(unc_f)
members <- merge(members, unc[, .(Orthogroup = og_2sp, infrv_member_mean,
                                  infrv_collapsed, infrv_gain, boot_r_within)],
                 by = "Orthogroup", all.x = TRUE)
pp <- fread(pair_f)
steal <- pp[, .(n_pairs = .N,
                n_pairs_determined = sum(n_samples_usable == 0L),
                n_pairs_stealing = sum(verdict == "read_stealing"),
                frac_stealing = sum(verdict == "read_stealing") /
                  pmax(sum(n_samples_usable > 0L), 1L)), by = .(Orthogroup = og_2sp)]
members <- merge(members, steal, by = "Orthogroup", all.x = TRUE)

# --- the classification ------------------------------------------------------
# Separability: a group is INSEPARABLE when most of its measurable pairs are
# trading reads. A group whose pairs are all `determined` had no contention at
# all, which is the strongest evidence of separability, not missing data.
members[, separability := fifelse(
  is.na(n_pairs) | n < 2L, NA_character_,
  fifelse(n_pairs_determined == n_pairs, "separable",
  fifelse(frac_stealing >= 0.5, "inseparable", "separable")))]
members[, uniformity := fifelse(
  is.na(expr_concordance) | n < 2L, NA_character_,
  fifelse(expr_concordance >= null_cut, "uniform", "divergent"))]
members[, cell := fifelse(is.na(separability) | is.na(uniformity), NA_character_,
                          paste(uniformity, separability, sep = "/"))]

out <- members[, .(og_2sp = Orthogroup, n_members = n,
                   expr_concordance, null_mean, null_cut,
                   n_tested, n_r_mean, n_r_sd, sign_concordance, n_responsive,
                   n_nodes = if ("n_nodes" %in% names(members)) n_nodes else NA_integer_,
                   degree_median = if ("degree_median" %in% names(members)) degree_median else NA_real_,
                   degree_fold_range = if ("degree_fold_range" %in% names(members)) degree_fold_range else NA_real_,
                   n_distinct_modules = if ("n_distinct_modules" %in% names(members)) n_distinct_modules else NA_integer_,
                   infrv_member_mean, infrv_collapsed, infrv_gain,
                   boot_r_within, n_pairs, n_pairs_determined, n_pairs_stealing,
                   frac_stealing, separability, uniformity, cell)]
cw <- fread(CROSSWALK, select = c("og_2sp", "og_3sp", "ambiguous", "sorghum_anchor",
                                  sprintf("class_%s", STUDY)))
out <- merge(out, cw, by = "og_2sp", all.x = TRUE)
# Round for reading. Full double precision on a concordance is noise, and these
# tables are meant to be opened and sorted by hand.
for (cl in c("expr_concordance", "null_mean", "null_cut", "n_r_mean", "n_r_sd",
             "sign_concordance", "degree_median", "degree_fold_range",
             "infrv_member_mean", "infrv_collapsed", "infrv_gain",
             "boot_r_within", "frac_stealing"))
  if (cl %in% names(out)) out[, (cl) := round(get(cl), 4)]
write_tsv(out[order(-n_members, og_2sp)],
          file.path(OUT_DIR, sprintf("og_uniformity_%s.tsv", STUDY)))

# --- the 2x2 -----------------------------------------------------------------
tab <- out[!is.na(cell), .N, by = .(uniformity, separability)][order(uniformity, separability)]
tab[, pct := round(100 * N / sum(N), 2)]
write_tsv(tab, file.path(OUT_DIR, sprintf("og_classification_%s.tsv", STUDY)))
banner("the 2x2")
for (i in seq_len(nrow(tab)))
  say(sprintf("  %-9s / %-12s %7s groups (%5.2f%%)",
              tab$uniformity[i], tab$separability[i], fmt_n(tab$N[i]), tab$pct[i]))
disc <- out[cell == "divergent/separable"]
say(sprintf("\nthe discovery set -- divergent AND separable: %s groups",
            fmt_n(nrow(disc))))
if (nrow(disc))
  say(sprintf("  median members %d | median degree fold-range %.1f | %s with a responsive member",
              as.integer(median(disc$n_members)),
              median(disc$degree_fold_range, na.rm = TRUE),
              fmt_n(disc[n_responsive > 0, .N])))

# --- the global ICCs ---------------------------------------------------------
# These are between-vs-within decompositions over ALL groups, so one number each.
icc_of <- function(vals, grp, label) {
  r <- icc_oneway(vals, grp)
  if (is.null(r)) return(data.table(metric = label, value = "NA"))
  say(sprintf("  ICC %-22s %.4f  (%s groups, %s values)",
              label, r$icc, fmt_n(r$n_groups), fmt_n(r$n_values)))
  data.table(metric = paste0("icc_", label), value = sprintf("%.6f", r$icc))
}
banner("ICC: how much of each quantity sits BETWEEN orthogroups")
md <- merge(og, tr, by.x = "Gene", by.y = "gene")
iccs <- rbindlist(list(icc_of(md$r_partial, md$Orthogroup, "nitrogen_r_partial")))
if (nzchar(NODES) && file.exists(NODES)) {
  nd <- merge(og, fread(NODES, select = c("gene", "degree"))[
    , .(gene = strip_version(gene), degree)], by.x = "Gene", by.y = "gene")
  iccs <- rbindlist(list(iccs, icc_of(log10(nd$degree + 1), nd$Orthogroup, "log_degree")))
}

mc <- out[n_members >= 2]
summ <- rbindlist(list(
  data.table(metric = c("study", "n_genes", "n_orthogroups", "n_multicopy",
                        "n_classified", "uniform_quantile", "steal_threshold",
                        "expr_concordance_median_multicopy",
                        "null_mean_n2", "null_cut_n2",
                        "pct_uniform", "pct_separable",
                        "n_divergent_separable", "n_divergent_inseparable",
                        "infrv_gain_median", "pct_infrv_gain_above_1"),
             value = c(STUDY, fmt_n(nrow(vst)), fmt_n(nrow(out)), fmt_n(nrow(mc)),
                       fmt_n(out[!is.na(cell), .N]),
                       sprintf("%.2f", UNIFORM_Q), sprintf("%.2f", STEAL_THR),
                       sprintf("%.4f", median(mc$expr_concordance, na.rm = TRUE)),
                       sprintf("%.4f", null_q[n == 2, null_mean]),
                       sprintf("%.4f", null_q[n == 2, null_cut]),
                       sprintf("%.2f", 100 * mean(out$uniformity == "uniform", na.rm = TRUE)),
                       sprintf("%.2f", 100 * mean(out$separability == "separable", na.rm = TRUE)),
                       fmt_n(out[cell == "divergent/separable", .N]),
                       fmt_n(out[cell == "divergent/inseparable", .N]),
                       sprintf("%.4f", median(mc$infrv_gain, na.rm = TRUE)),
                       sprintf("%.2f", 100 * mean(mc$infrv_gain > 1, na.rm = TRUE)))),
  iccs))
write_tsv(summ, file.path(OUT_DIR, sprintf("og_uniformity_%s_summary.tsv", STUDY)))
write_tsv(null_q, file.path(OUT_DIR, sprintf("og_uniformity_%s_null.tsv", STUDY)))

say("done")
