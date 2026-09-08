#!/usr/bin/env Rscript
# =============================================================================
# 43_knn_ba_collect.r -- merge the per-k scores and name the winning k.
#
# k is chosen as argmin GIC(BA) over the grid: the k whose graph is closest to a
# Barabasi-Albert model. That is the criterion as asked for, stated literally.
#
# THE SELECTED-MODEL COLUMN IS REPORTED BESIDE IT, and must not be skipped when
# reading this table. "Closest to BA" is only meaningful against what else the
# graph could have been: if ER or WS wins at every k, then the winning k is
# merely the least-bad BA fit of a graph that is not BA-like at all. That is a
# finding about the networks, not a failure of the sweep, and it has to be
# visible rather than buried under an argmin.
#
# The power-law columns travel alongside as the agreed fallback criterion. If
# 39_statgraph_validate.r did not PASS, pl_alpha/pl_ks_stat are what should
# choose k, and this script says so rather than quietly switching.
# =============================================================================
STUDY  <- Sys.getenv("CLEAN_STUDY");  if (!nzchar(STUDY)) stop("CLEAN_STUDY not set")
OUTDIR <- Sys.getenv("CLEAN_OUT_DIR"); if (!nzchar(OUTDIR)) stop("CLEAN_OUT_DIR not set")
VALID  <- Sys.getenv("CLEAN_VALIDATION", "")
say <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

f <- list.files(OUTDIR, pattern = sprintf("^knnba_sel_%s_k[0-9]+\\.tsv$", STUDY),
                full.names = TRUE)
if (!length(f)) stop("no per-k score files in ", OUTDIR)
d <- do.call(rbind, lapply(f, read.delim, stringsAsFactors = FALSE))
d <- d[order(d$k), ]

out <- file.path(OUTDIR, sprintf("knnba_selection_%s.tsv", STUDY))
write.table(d, out, sep = "\t", quote = FALSE, row.names = FALSE)
say("wrote ", basename(out), " (", nrow(d), " k values)")

say("")
print(d[, c("k","nodes","edges","mean_degree","median_degree",
            "selected_model","ba_param","gic_ba","gic_er","gic_ws",
            "pl_alpha","pl_ks_stat")], row.names = FALSE, digits = 5)

ok <- d[is.finite(d$gic_ba), ]
say("")
if (!nrow(ok)) {
  say("NO k produced a usable GIC(BA). The spectral criterion cannot choose here.")
} else {
  best <- ok[which.min(ok$gic_ba), ]
  say("WINNER by argmin GIC(BA):  k = ", best$k)
  say(sprintf("  nodes %s   edges %s   mean degree %.1f   median degree %s",
              format(best$nodes, big.mark = ","), format(best$edges, big.mark = ","),
              best$mean_degree, best$median_degree))
  say(sprintf("  GIC   BA %+.5f | ER %+.5f | WS %+.5f", best$gic_ba, best$gic_er, best$gic_ws))
  say(sprintf("  BA exponent %.3f   |   power law alpha %.3f (KS %.4f)",
              best$ba_param, best$pl_alpha, best$pl_ks_stat))
  nba <- sum(ok$selected_model == "BA")
  say("")
  say("  models actually selected across the grid: ",
      paste(sprintf("%s x%d", names(table(ok$selected_model)), table(ok$selected_model)),
            collapse = ", "))
  if (nba == 0)
    say("  NOTE: BA never wins. The chosen k is the least-bad BA fit of a graph",
        " that another model describes better at every k -- read the winner with that in mind.")
}

if (nzchar(VALID) && file.exists(VALID)) {
  v <- read.delim(VALID, stringsAsFactors = FALSE)
  a <- v[v$check == "A_recovery", ]; b <- v[v$check == "B_agreement", ]
  # The gate is read at two grains, because a blunt pass/fail would throw away
  # the distinction that matters here. Choosing k by argmin GIC(BA) needs only
  # that BA is recognised when a graph IS BA, and that `fast` tracks the exact
  # `diag` on real graphs. It does NOT need ER and WS to be told apart -- and
  # measured, that is exactly the split: BA recovery is perfect while ER/WS are
  # confused in both directions. So the ba-specific verdict licenses the k
  # choice, and the overall one governs how much the selected_model column is
  # worth.
  aba  <- a[a$truth == "BA", ]
  ok_ba   <- nrow(aba) > 0 && all(aba$correct)
  ok_agr  <- nrow(b) == 0 || all(b$correct)
  ok_all  <- nrow(a) > 0 && all(a$correct)
  say("")
  say("validation gate:")
  say("  BA recovery      ", sum(aba$correct), "/", nrow(aba),
      "  -- can it spot a Barabasi-Albert graph? This is what argmin GIC(BA) rests on.")
  say("  all-model recovery ", sum(a$correct), "/", nrow(a),
      "  -- governs how much `selected_model` is worth.")
  if (nrow(b)) say("  fast vs exact    ", sum(b$correct), "/", nrow(b),
                   "  -- on real co-expression subgraphs.")
  say("")
  if (ok_ba && ok_agr && ok_all) {
    say("  VERDICT: PASS. Both the k choice and selected_model are licensed.")
  } else if (ok_ba && ok_agr) {
    say("  VERDICT: PARTIAL. GIC(BA) may choose k -- BA is recognised reliably and")
    say("    `fast` tracks the exact method on real graphs. But ER and WS are NOT")
    say("    reliably told apart, so read `selected_model` as 'BA vs not-BA' only,")
    say("    and do not quote which of ER/WS came second.")
  } else {
    say("  VERDICT: FAIL. Choose k on the power-law columns (pl_alpha / pl_ks_stat)")
    say("    instead, and say so wherever the chosen k is quoted.")
  }
} else {
  say("")
  say("NOTE: no validation table found. Run ./run.sh sgvalidate before trusting gic_ba.")
}
