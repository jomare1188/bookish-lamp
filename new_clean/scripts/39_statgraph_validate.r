#!/usr/bin/env Rscript
# =============================================================================
# 39_statgraph_validate.r -- can statGraph's criterion actually choose k here?
#
# The plan is to pick the k-NN parameter per network by whichever k makes the
# graph closest to a Barabasi-Albert model, using statGraph::graph.model.selection.
# Before that criterion is allowed to choose anything, it has to earn it. This
# script is the gate.
#
# WHY A GATE IS NEEDED, specifically:
#
#   1. The default path cannot run. method="diag" computes a FULL dense
#      eigendecomposition, and not once: 502 grid points x 50 simulated graphs =
#      25,100 of them. Measured here, eigen() is cleanly cubic (exponent 3.10);
#      at n = 170,736 one call is ~34.6 h and 233 GB (466 GB peak, i.e. OOM), so
#      the default is ~99 years. Only method="fast" is available to us.
#
#   2. method="fast" derives the spectral density from the DEGREE DISTRIBUTION
#      via a cavity/locally-tree-like approximation. Co-expression networks are
#      dense and highly clustered -- precisely the regime where that
#      approximation is weakest. Whether it still ranks ER/WS/BA correctly on
#      our graphs is an empirical question.
#
#   3. statGraph's default ER grid starts at p = 0, which is an empty graph with
#      no spectral density, and the whole call dies with "missing value where
#      TRUE/FALSE needed". Every model selection here must pass explicit,
#      non-degenerate parameter ranges. (Verified on statGraph 1.0.6.)
#
# THE THREE CHECKS
#   A  recovery   -- on synthetic graphs of KNOWN type, does `fast` name the
#                    right model? Run at several n, including our real size.
#   B  agreement  -- on real subgraphs small enough for `diag`, does `fast`
#                    agree with the exact method?
#   C  cost       -- how long does one spectral density take on a REAL reduced
#                    graph, which has thousands of unique degrees rather than
#                    the ~25 a sparse synthetic one has?
#
# A FAIL is not a crash. It is a reported verdict that switches the sweep to the
# fallback criterion (igraph::fit_power_law, Clauset-Shalizi-Newman), which is
# what the user chose. It must never be substituted silently.
#
# RUN: through run.sh  ->  ./run.sh sgvalidate
# =============================================================================
suppressPackageStartupMessages({library(igraph)})
LIB <- Sys.getenv("CLEAN_RLIB", "")
if (nzchar(LIB)) .libPaths(c(LIB, .libPaths()))
suppressPackageStartupMessages({library(statGraph)})

OUT      <- Sys.getenv("CLEAN_OUT", "statgraph_validation.tsv")
NCORES   <- as.integer(Sys.getenv("CLEAN_CORES", "16"))
BIG_N    <- as.integer(Sys.getenv("CLEAN_BIG_N", "170736"))
REAL_EDG <- Sys.getenv("CLEAN_REAL_EDGES", "")   # optional: a real edge list for B/C
set.seed(as.integer(Sys.getenv("CLEAN_SEED", "1188")))

say <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
rows <- list()
add  <- function(...) rows[[length(rows) + 1L]] <<- data.frame(..., stringsAsFactors = FALSE)

# --- parameter ranges that exclude the degenerate ends ------------------------
# ER at p = 0 and WS at p = 0 are empty / fully-regular graphs whose spectral
# density is undefined or singular. statGraph's defaults walk straight into the
# first of those, so the ranges are always given explicitly here.
PARAMS <- function(n, m) list(
  ER = seq(max(1e-6, m / (n - 1) / 4), min(1, 4 * m / (n - 1)), length.out = 25),
  WS = seq(0.01, 1.00, length.out = 25),
  BA = seq(0.10, 3.00, length.out = 25))

# SINGLE-THREADED ON PURPOSE. statGraph's numCores path creates and destroys a
# PSOCK cluster per spectral density; across dozens of GIC calls it leaks
# connections and deadlocks -- measured: hangs at numCores=32 with the main
# process at 0% CPU, fine at 4-16. Parallelism belongs outside, across k.
select_model <- function(G, method = "fast", ncores = 1) {
  n <- vcount(G); m <- mean(degree(G))
  p <- PARAMS(n, m)
  t0 <- Sys.time()
  r <- try(graph.model.selection(G, models = c("ER", "WS", "BA"),
                                 parameters = list(p$ER, p$WS, p$BA),
                                 method = method), silent = TRUE)
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(r, "try-error"))
    return(list(model = "ERROR", gic = NA_real_, param = NA_real_, secs = el,
                est = NULL, msg = as.character(r)))
  e <- r$estimates
  gic <- setNames(as.numeric(e[, "GIC"]), rownames(e))
  par <- setNames(as.numeric(e[, "param"]), rownames(e))
  list(model = r$model[1], gic = gic, param = par, secs = el, est = e, msg = "")
}

# =============================================================================
say("CHECK A -- does `fast` recover a KNOWN model?")
# If it cannot name the generator of a graph it was handed, it cannot be trusted
# to say which k makes our graph most BA-like.
for (n in as.integer(strsplit(Sys.getenv("CLEAN_A_N","1000,10000"),",")[[1]])) {
  truth <- list(
    BA = sample_pa(n, power = 1, m = 4, directed = FALSE),
    ER = sample_gnp(n, 8 / (n - 1)),
    WS = sample_smallworld(1, n, 4, 0.05))
  for (nm in names(truth)) {
    G <- simplify(truth[[nm]])
    r <- select_model(G)
    ok <- identical(r$model, nm)
    say(sprintf("  n=%-7d truth=%-3s -> selected=%-6s %s  (%.1fs)",
                n, nm, r$model, if (ok) "OK" else "MISMATCH", r$secs))
    add(check = "A_recovery", n = n, truth = nm, selected = r$model,
        correct = ok, secs = round(r$secs, 2),
        gic_ba = unname(r$gic["BA"]), gic_er = unname(r$gic["ER"]),
        gic_ws = unname(r$gic["WS"]), ba_param = unname(r$param["BA"]))
  }
}

# =============================================================================
say("")
say("CHECK B -- does `fast` agree with the exact `diag` on REAL subgraphs?")
if (!nzchar(REAL_EDG) || !file.exists(REAL_EDG)) {
  say("  no real edge list supplied (CLEAN_REAL_EDGES) -- skipped")
} else {
  el <- data.table::fread(REAL_EDG, header = FALSE, select = 1:2,
                          col.names = c("a", "b"))
  Gfull <- simplify(graph_from_data_frame(el, directed = FALSE))
  say("  real graph: ", vcount(Gfull), " nodes, ", ecount(Gfull), " edges")
  for (ns in c(600L, 1200L)) {
    # SNOWBALL, not uniform node sampling. A uniform sample of 600 nodes from a
    # 98,000-node graph induces almost no edges -- measured: 133 nodes and 78
    # edges, a fragment with none of the density or clustering that makes this
    # check worth running. Growing outward from a seed keeps a real
    # neighbourhood, which is exactly the regime where the fast method's
    # locally-tree-like assumption is under test.
    seed <- sample(V(Gfull)[degree(Gfull) > 0], 1)
    frontier <- as.integer(seed); keep <- frontier
    while (length(keep) < ns && length(frontier)) {
      nb <- unique(unlist(adjacent_vertices(Gfull, frontier)))
      nb <- setdiff(nb, keep)
      if (!length(nb)) break
      take <- head(sample(nb), ns - length(keep))
      keep <- c(keep, take); frontier <- take
    }
    Gs <- simplify(induced_subgraph(Gfull, keep))
    Gs <- delete_vertices(Gs, which(degree(Gs) == 0))
    if (vcount(Gs) < 50) { say("  n=", ns, " subgraph too sparse -- skipped"); next }
    rf <- select_model(Gs, "fast")
    rd <- select_model(Gs, "diag", ncores = 1)
    agree <- identical(rf$model, rd$model)
    say(sprintf("  sub n=%-5d (%d edges)  fast=%-6s diag=%-6s %s  (%.0fs / %.0fs)",
                vcount(Gs), ecount(Gs), rf$model, rd$model,
                if (agree) "AGREE" else "DISAGREE", rf$secs, rd$secs))
    add(check = "B_agreement", n = vcount(Gs), truth = "real_subgraph",
        selected = paste0("fast:", rf$model, "|diag:", rd$model),
        correct = agree, secs = round(rf$secs + rd$secs, 2),
        gic_ba = unname(rf$gic["BA"]), gic_er = unname(rf$gic["ER"]),
        gic_ws = unname(rf$gic["WS"]), ba_param = unname(rf$param["BA"]))
  }
}

# =============================================================================
say("")
say("CHECK C -- what does one spectral density cost on a REALISTIC degree profile?")
# The published 0.96 s benchmark used sparse ER graphs with ~25 unique degrees.
# A k-NN-reduced co-expression graph has thousands, and the fast method's cost
# scales with that, not with n. This is the number that sets the sweep's budget.
Gc <- simplify(sample_pa(BIG_N, power = 1, m = 20, directed = FALSE))
say("  probe graph: ", vcount(Gc), " nodes, ", ecount(Gc), " edges, ",
    length(unique(degree(Gc))), " unique degrees")
for (nc in c(1L)) {
  t0 <- Sys.time()
  invisible(graph.spectral.density(Gc, method = "fast", numCores = nc))
  s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  say(sprintf("  numCores=%-4d %.1f s per spectral density", nc, s))
  add(check = "C_cost", n = vcount(Gc), truth = paste0("numCores=", nc),
      selected = NA, correct = NA, secs = round(s, 2),
      gic_ba = NA, gic_er = NA, gic_ws = NA, ba_param = NA)
}

# =============================================================================
res <- do.call(rbind, rows)
write.table(res, OUT, sep = "\t", quote = FALSE, row.names = FALSE)
say("")
say("wrote ", OUT)

a <- res[res$check == "A_recovery", ]
b <- res[res$check == "B_agreement", ]
passA <- nrow(a) > 0 && all(a$correct)
passB <- nrow(b) == 0 || all(b$correct)
say("")
say("CHECK A recovery : ", sum(a$correct, na.rm = TRUE), "/", nrow(a),
    " known models named correctly")
if (nrow(b)) say("CHECK B agreement: ", sum(b$correct), "/", nrow(b),
                 " fast/diag agreements on real subgraphs")
say("")
if (passA && passB) {
  say("VERDICT: PASS -- graph.model.selection(method='fast') may choose k.")
} else {
  say("VERDICT: FAIL -- the spectral criterion does not reliably identify the")
  say("  model class at this scale on these graphs. The sweep must use the")
  say("  fallback (igraph::fit_power_law, Clauset-Shalizi-Newman) and say so.")
  say("  Failing rows:")
  bad <- rbind(a[!a$correct, ], b[!b$correct, ])
  print(bad[, c("check", "n", "truth", "selected")], row.names = FALSE)
}
