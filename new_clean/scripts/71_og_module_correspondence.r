#!/usr/bin/env Rscript
# =============================================================================
# 71_og_module_correspondence.r -- do the two species put the same orthogroups in
# the same module?
#
# WHY THIS COULD NOT BE ASKED BEFORE. Comparing two clusterings needs them to be
# clusterings OF THE SAME THINGS. At gene level they are not: sugarcane's modules
# partition 101,990 sugarcane genes and purple's partition 170,135 purple genes,
# and the only bridge is a many-to-many ortholog map. ARI and NMI have no meaning
# across two different vertex sets. 68 put both networks on one shared vertex set,
# so the question becomes well posed for the first time.
#
# WHAT IS COMPARED, AND WHAT IS LEFT OUT. Only orthogroups that are a node in
# BOTH networks -- i.e. carry at least one edge in each species. An orthogroup
# isolated in one species has no module there, and scoring it as "its own cluster"
# would manufacture agreement: both species would be recorded as placing it alone,
# which is true and uninformative, and at this scale the singletons would dominate
# every index. The count used is reported, and so is the count discarded.
#
# TWO NULLS, BECAUSE THE INDICES ANSWER DIFFERENT WORRIES.
#   label    permute one partition's module labels over the same nodes. Destroys
#            the correspondence while keeping both module-size distributions. This
#            is the standard null for ARI/NMI.
#   degree   relabel nodes WITHIN degree strata, as 69 does. Answers the harder
#            question: is the agreement more than two graphs of similar degree
#            structure would produce anyway?
#
# ARI AND NMI DISAGREE ON PURPOSE. ARI is chance-corrected and punishes splitting;
# NMI is not chance-corrected and rises with the number of clusters. Reporting one
# alone would be choosing an answer. Both are reported against both nulls.
#
# RUN: through run.sh -> ./run.sh ogmodules
# =============================================================================

suppressMessages({
  library(data.table)
  library(igraph)
})

source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE),
                                                 value = TRUE)[1])), "lib", "common.R"))

STUDIES    <- env_list("CLEAN_STUDIES", c("sugarcane", "purple"))
OUT_DIR    <- env_req("CLEAN_OUT_DIR")
NULL_REPS  <- as.integer(env_num("CLEAN_NULL_REPS", 100))
DEG_STRATA <- as.integer(env_num("CLEAN_DEGREE_STRATA", 20))
MIN_SIZE   <- as.integer(env_num("CLEAN_MIN_MODULE_SIZE", 2))
SEED       <- as.integer(env_num("CLEAN_SEED", 1))

memb_f <- setNames(vapply(STUDIES, function(s)
  env_req(sprintf("CLEAN_MEMBERSHIP_%s", toupper(s))), ""), STUDIES)

banner("module correspondence on a shared vertex set")
ensure_dir(OUT_DIR)
set.seed(SEED)

# --- load both partitions ----------------------------------------------------
M <- list()
for (s in STUDIES) {
  if (!file.exists(memb_f[[s]]))
    stop(basename(memb_f[[s]]), " missing\n  run the MCL chain for ", s,
         " in this tree first", call. = FALSE)
  d <- fread(memb_f[[s]], select = c("gene", "module_name", "degree"))
  setnames(d, "gene", "og")
  M[[s]] <- d
  say(sprintf("%-10s %s nodes | %s named modules | %s unassigned",
              s, fmt_n(nrow(d)), fmt_n(d[module_name != "Unassigned",
                                         uniqueN(module_name)]),
              fmt_n(d[module_name == "Unassigned", .N])))
}

A <- STUDIES[1]; B <- STUDIES[2]
shared <- intersect(M[[A]]$og, M[[B]]$og)
say(sprintf("\nnodes in both networks: %s", fmt_n(length(shared))))
for (s in STUDIES)
  say(sprintf("  %-10s has %s nodes the other does not",
              s, fmt_n(nrow(M[[s]]) - length(shared))))

dt <- merge(M[[A]][og %chin% shared, .(og, a = module_name, deg_a = degree)],
            M[[B]][og %chin% shared, .(og, b = module_name, deg_b = degree)],
            by = "og")

# Unassigned is mcl's label for a module below MCL_MIN_MODULE_SIZE, not a module.
# Keeping it would merge every such node into one giant pseudo-cluster per species
# and both indices would then be measuring that artefact.
n_before <- nrow(dt)
dt <- dt[a != "Unassigned" & b != "Unassigned"]
say(sprintf("dropping %s nodes unassigned in either species -> %s compared",
            fmt_n(n_before - nrow(dt)), fmt_n(nrow(dt))))
if (nrow(dt) < 100)
  stop("only ", nrow(dt), " nodes are in a named module in both species -- ",
       "nothing to compare", call. = FALSE)

ga <- as.integer(factor(dt$a))
gb <- as.integer(factor(dt$b))
say(sprintf("modules among the compared nodes: %s (%s) | %s (%s)",
            fmt_n(uniqueN(ga)), A, fmt_n(uniqueN(gb)), B))

ari <- igraph::compare(ga, gb, method = "adjusted.rand")
nmi <- igraph::compare(ga, gb, method = "nmi")
say(sprintf("\nobserved   ARI %.6f | NMI %.6f", ari, nmi))

# --- nulls -------------------------------------------------------------------
run_null <- function(kind) {
  out <- matrix(NA_real_, NULL_REPS, 2,
                dimnames = list(NULL, c("ari", "nmi")))
  if (kind == "degree") {
    ord <- order(dt$deg_b)
    strata <- split(ord, cut(seq_along(ord), DEG_STRATA, labels = FALSE))
  }
  for (i in seq_len(NULL_REPS)) {
    if (kind == "label") {
      perm <- sample(seq_along(gb))
    } else {
      perm <- seq_along(gb)
      for (st in strata) perm[st] <- st[sample(length(st))]
    }
    out[i, 1] <- igraph::compare(ga, gb[perm], method = "adjusted.rand")
    out[i, 2] <- igraph::compare(ga, gb[perm], method = "nmi")
  }
  out
}

res <- list()
for (kind in c("label", "degree")) {
  say(sprintf("null (%s): %d replicates", kind, NULL_REPS))
  nl <- run_null(kind)
  for (ix in c(ari = 1L, nmi = 2L)) {
    stat <- c(ari, nmi)[ix]
    v <- nl[, ix]
    nm <- mean(v); sd_ <- sd(v)
    z <- if (is.finite(sd_) && sd_ > 0) (stat - nm) / sd_ else NA_real_
    p <- (sum(v >= stat) + 1) / (NULL_REPS + 1)
    # A FOLD IS MEANINGLESS FOR ARI. It is already chance-corrected, so its null
    # mean sits at 0 by construction (measured -5e-05 here) and dividing by that
    # produces a number like -1133 that reads as enormous and means nothing. The
    # excess over the null is the interpretable quantity for ARI; NMI is not
    # chance-corrected, so for it the fold is meaningful and both are carried.
    idx <- c("ARI", "NMI")[ix]
    fold <- if (idx == "NMI" && is.finite(nm) && abs(nm) > 1e-6) stat / nm else NA_real_
    res[[length(res) + 1]] <- data.table(
      index = idx, null = kind, observed = stat,
      null_mean = nm, null_sd = sd_, excess = stat - nm,
      fold = fold, z = z, p_empirical = p)
    say(sprintf("  %-3s observed %.6f | null %.6f (sd %.2g) | excess %+.6f | %s z %6.1f | p %.4g",
                idx, stat, nm, sd_, stat - nm,
                if (is.na(fold)) "            " else sprintf("fold %5.2f |", fold), z, p))
  }
}
summ <- rbindlist(res)

# --- where the agreement actually is -----------------------------------------
# One index over 38,000 nodes says nothing about WHICH modules correspond. For
# every module of A, the B module it overlaps most, with the overlap as a
# fraction of A's module -- one denominator, so the column is rankable.
best <- dt[, .N, by = .(a, b)][order(a, -N)]
sz_a <- dt[, .(size_a = .N), by = a]
sz_b <- dt[, .(size_b = .N), by = b]
best <- best[, .SD[1], by = a]
best <- merge(merge(best, sz_a, by = "a"), sz_b, by = "b")
best[, frac_of_a := N / size_a]
setorder(best, -size_a)
# SIZE MATTERS HERE. The median compared module holds a handful of nodes, and for
# a 3-node module a 2-node overlap already scores 67%. An unstratified "median
# best-match fraction" therefore says more about module size than about
# correspondence, so it is reported per size band.
best[, size_band := cut(size_a, c(0, 3, 5, 10, 25, Inf),
                        labels = c("2-3", "4-5", "6-10", "11-25", ">25"))]
setnames(best, c("a", "b", "N"), c(paste0("module_", A), paste0("module_", B),
                                   "n_shared"))
write_tsv(best, file.path(OUT_DIR, "og_module_best_match.tsv"))
say(sprintf("\nbest-match table: %s modules of %s", fmt_n(nrow(best)), A))
say(sprintf("  median overlap with its best partner: %.1f%% of the %s module",
            100 * median(best$frac_of_a), A))
say(sprintf("  modules whose best partner holds >= 50%% of them: %s (%.1f%%)",
            fmt_n(best[frac_of_a >= 0.5, .N]),
            100 * best[frac_of_a >= 0.5, .N] / nrow(best)))
band <- best[, .(n_modules = .N, median_frac = round(median(frac_of_a), 4),
                 pct_ge_50 = round(100 * mean(frac_of_a >= 0.5), 1)), by = size_band]
setorder(band, size_band)
say("  by module size, because a 2-of-3 overlap already scores 67%:")
for (i in seq_len(nrow(band)))
  with(band[i], say(sprintf("    %-6s %6s modules | median overlap %5.1f%% | %5.1f%% at >= 50%%",
                            as.character(size_band), fmt_n(n_modules),
                            100 * median_frac, pct_ge_50)))
write_tsv(band, file.path(OUT_DIR, "og_module_best_match_by_size.tsv"))

write_tsv(summ, file.path(OUT_DIR, "og_module_correspondence.tsv"))
write_tsv(data.table(
  metric = c("nodes_in_both_networks", "nodes_compared",
             "dropped_unassigned", paste0("modules_", c(A, B)),
             "null_reps", "degree_strata", "ari", "nmi",
             "median_best_match_frac", "pct_modules_best_match_ge_50"),
  value = c(fmt_n(length(shared)), fmt_n(nrow(dt)), fmt_n(n_before - nrow(dt)),
            fmt_n(uniqueN(ga)), fmt_n(uniqueN(gb)), NULL_REPS, DEG_STRATA,
            sprintf("%.6f", ari), sprintf("%.6f", nmi),
            sprintf("%.4f", median(best$frac_of_a)),
            sprintf("%.2f", 100 * best[frac_of_a >= 0.5, .N] / nrow(best)))),
  file.path(OUT_DIR, "og_module_correspondence_summary.tsv"))
say("done")
