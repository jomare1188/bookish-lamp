#!/usr/bin/env Rscript
# =============================================================================
# 42_knn_model_selection.r -- score one k-NN-reduced graph against ER / WS / BA.
#
# ONE k PER PROCESS, deliberately. statGraph's own numCores path builds and tears
# down a PSOCK cluster for every spectral density; across a few dozen GIC calls
# that leaks connections and deadlocks (measured: hangs at numCores=32, fine at
# 4-16). So statGraph runs single-threaded here and the parallelism sits outside,
# across k values -- which is also the better split, since the k jobs are
# independent and there are 256 cores.
#
# THE PARAMETER GRIDS ARE GIVEN EXPLICITLY, and that is not a detail.
# statGraph's default ER grid is seq(0, 1, 0.01). At p = 0 the model graph is
# empty, its spectral density is undefined, and the whole call dies with
# "missing value where TRUE/FALSE needed". Worse, at p = 1e-06 it returns a
# spuriously NEGATIVE GIC that beats every real model -- so a graph generated as
# BA gets called ER. Verified on statGraph 1.0.6. Every range below excludes the
# degenerate end, and the ER range is centred on the observed mean degree rather
# than spanning [0,1], which is both correct and far cheaper.
#
# TWO CRITERIA ARE ALWAYS WRITTEN, not one.
#   gic_ba  -- the spectral criterion the user asked for. Its `fast` method is a
#              cavity/locally-tree-like approximation, and co-expression graphs
#              are dense and clustered, so 39_statgraph_validate.r has to pass
#              before this column is allowed to choose k.
#   pl_alpha, pl_ks_p -- a direct Clauset-Shalizi-Newman power-law fit on the
#              degree sequence. Cheap, assumption-free about topology, and the
#              agreed fallback. Computing both always means the fallback needs no
#              re-run if the gate fails, and the two can be compared.
#
# RUN: through run.sh  ->  RESULTS=... ./run.sh knnselect purple [k]
# =============================================================================
suppressPackageStartupMessages({library(igraph)})
LIB <- Sys.getenv("CLEAN_RLIB", ""); if (nzchar(LIB)) .libPaths(c(LIB, .libPaths()))
suppressPackageStartupMessages({library(statGraph)})

STUDY  <- Sys.getenv("CLEAN_STUDY");   if (!nzchar(STUDY)) stop("CLEAN_STUDY not set")
GRID   <- Sys.getenv("CLEAN_GRID");    if (!nzchar(GRID))  stop("CLEAN_GRID not set")
OUTDIR <- Sys.getenv("CLEAN_OUT_DIR"); if (!nzchar(OUTDIR)) stop("CLEAN_OUT_DIR not set")
ONE_K  <- Sys.getenv("CLEAN_K", "")
NPAR   <- as.integer(Sys.getenv("CLEAN_NPARAM", "20"))
set.seed(as.integer(Sys.getenv("CLEAN_SEED", "1188")))
say <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

g <- read.delim(GRID, stringsAsFactors = FALSE)
if (nzchar(ONE_K)) g <- g[g$k == as.integer(ONE_K), , drop = FALSE]
if (!nrow(g)) stop("no rows in ", GRID, if (nzchar(ONE_K)) paste0(" for k=", ONE_K) else "")

# Ranges that exclude the degenerate ends; ER is centred on the observed density.
param_grids <- function(n, meandeg) list(
  ER = seq(max(1e-6, meandeg / (n - 1) / 4), min(0.999, 4 * meandeg / (n - 1)), length.out = NPAR),
  WS = seq(0.01, 1.00, length.out = NPAR),
  BA = seq(0.10, 3.00, length.out = NPAR))

for (i in seq_len(nrow(g))) {
  k  <- g$k[i]
  el <- g$edge_list[i]
  outf <- file.path(OUTDIR, sprintf("knnba_sel_%s_k%s.tsv", STUDY, k))
  if (file.exists(outf) && file.size(outf) > 0) { say("k=", k, " already scored"); next }
  if (!file.exists(el)) { say("k=", k, " edge list missing: ", el); next }

  e <- data.table::fread(el, header = FALSE, select = 1:2, col.names = c("a", "b"))
  G <- simplify(graph_from_data_frame(e, directed = FALSE))
  G <- delete_vertices(G, which(degree(G) == 0))
  n <- vcount(G); md <- mean(degree(G))
  say(sprintf("k=%-4s  %s nodes, %s edges, mean degree %.1f", k,
              format(n, big.mark = ","), format(ecount(G), big.mark = ","), md))

  # --- the spectral criterion ---
  t0 <- Sys.time()
  p <- param_grids(n, md)
  r <- try(graph.model.selection(G, models = c("ER", "WS", "BA"),
                                 parameters = list(p$ER, p$WS, p$BA),
                                 method = "fast"), silent = TRUE)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(r, "try-error")) {
    sel <- "ERROR"; gic <- setNames(rep(NA_real_, 3), c("ER","WS","BA"))
    par <- gic; say("  model selection FAILED: ", sub("\n.*", "", as.character(r)))
  } else {
    sel <- r$model[1]
    gic <- setNames(as.numeric(r$estimates[, "GIC"]),   rownames(r$estimates))
    par <- setNames(as.numeric(r$estimates[, "param"]), rownames(r$estimates))
    say(sprintf("  selected %-3s | GIC  BA %+.5f  ER %+.5f  WS %+.5f | BA param %.3f | %.0fs",
                sel, gic["BA"], gic["ER"], gic["WS"], par["BA"], secs))
  }

  # --- the fallback criterion, always computed ---
  d <- degree(G); d <- d[d > 0]
  f <- try(fit_power_law(d, implementation = "plfit"), silent = TRUE)
  alpha <- if (inherits(f, "try-error")) NA_real_ else f$alpha
  xmin  <- if (inherits(f, "try-error")) NA_real_ else f$xmin
  ksp   <- if (inherits(f, "try-error") || is.null(f$KS.p)) NA_real_ else f$KS.p
  kss   <- if (inherits(f, "try-error") || is.null(f$KS.stat)) NA_real_ else f$KS.stat
  say(sprintf("  power law: alpha %.3f  xmin %.0f  KS %.4f", alpha, xmin, kss))

  write.table(data.frame(
      study = STUDY, k = k, mode = g$mode[i], nodes = n, edges = ecount(G),
      mean_degree = round(md, 2), median_degree = median(degree(G)),
      selected_model = sel,
      ba_param = unname(par["BA"]), gic_ba = unname(gic["BA"]),
      gic_er = unname(gic["ER"]), gic_ws = unname(gic["WS"]),
      pl_alpha = alpha, pl_xmin = xmin, pl_ks_stat = kss, pl_ks_p = ksp,
      secs = round(secs, 1), stringsAsFactors = FALSE),
    outf, sep = "\t", quote = FALSE, row.names = FALSE)
  say("  wrote ", basename(outf))
}
