#!/usr/bin/env Rscript
# =============================================================================
# 51_mcl_membership_from_cls.r — the clustering contract, from a partition that
# already exists.
#
# WHY NOT 05_mcl_clustering.r. 05 re-runs mcl from the 9-column edge table. For
# these networks that table was deleted (70 GB whose only consumer was mcxload),
# and re-running would spend hours reproducing a partition that is already on
# disk: the inflation ladder wrote `cls.knone.I15` (sugarcane) and
# `cls.knone.I35` (purple), and those two cells ARE the chosen settings — each
# species' modularity optimum, the argument figure 12 makes.
#
# So this adopts the existing partition instead of recomputing it, the same move
# 27_sbm_membership.r makes for an SBM fit. The contract is copied from 27
# exactly, because every stage downstream takes these two files as given:
#
#   <prefix>_membership.tsv      gene, module_name, strength, degree
#   <prefix>_module_summary.tsv  module, n_genes, modularity_Q, inflation,
#                                min_module_size, shown_in_plot
#
# `module_name` is Module_%03d ranked LARGEST FIRST, or the literal "Unassigned"
# for genes in blocks below MIN_SIZE — both conventions matter downstream:
# 22_fig_topology.r separates real modules with grepl("^Module_", module) and
# 14_module_eigengene.r drops "Unassigned" by name.
#
# THE GUARD THAT MATTERS. An mcl cell tag strips the decimal point, so `-I 1.2`
# and `-I 12` collide on `cls.knone.I12`; the sweep now writes the actual
# inflation into a `.inflation` sidecar beside every cell for exactly this
# reason. This script REFUSES to adopt a partition whose sidecar disagrees with
# the inflation requested. Adopting the wrong cell would be invisible — the
# module table would look perfectly normal and every downstream number would be
# about a different clustering.
#
# strength and degree are whole-graph values, as they are in 05 (:121-122) and
# 27, and are read from network_<study>_node_metrics.tsv rather than recomputed.
# Only one consumer uses `strength` at all: 16_module_heatmaps.r, to take the top
# N genes of a module too large to draw whole.
#
# Modularity is READ from mcl_sweep_<study>.tsv rather than recomputed: it was
# already measured there by `clm info` against this exact matrix and partition,
# and recomputing it in igraph would need the edge table back.
#
# RUN: through run.sh  ->  ./run.sh membership sugarcane
# =============================================================================
suppressMessages({ library(data.table) })
source(file.path(dirname(sub("--file=", "",
       grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "lib", "common.R"))

STUDY     <- env_req("CLEAN_STUDY")
CLS       <- env_req("CLEAN_CLS")           # mcl cluster file
TAB       <- env_req("CLEAN_TAB")           # index -> gene
NODES     <- env_req("CLEAN_NODE_METRICS")
SWEEP     <- env_opt("CLEAN_SWEEP_TSV", "")
PREFIX    <- env_req("CLEAN_PREFIX")
INFLATION <- env_num("CLEAN_INFLATION", NA_real_)
MIN_SIZE  <- as.integer(env_num("CLEAN_MIN_MODULE_SIZE", 2))
MIN_PLOT  <- as.integer(env_num("CLEAN_MIN_MODULE_SIZE_PLOT", 10))
setDTthreads(as.integer(env_num("CLEAN_CORES", 32)))

banner(sprintf("adopting MCL partition: %s  (inflation %s)", STUDY, INFLATION))

for (f in c(CLS, TAB, NODES))
  if (!file.exists(f)) stop("missing input: ", f, call. = FALSE)

# --- the guard: this cell must be the inflation we asked for -----------------
side <- paste0(CLS, ".inflation")
if (file.exists(side)) {
  got <- as.numeric(trimws(readLines(side, warn = FALSE)[1]))
  if (!is.na(INFLATION) && !isTRUE(all.equal(got, INFLATION))) {
    stop(sprintf(paste("%s was produced with -I %s, not -I %s.",
                       "\n  Cell tags strip the decimal point, so -I 1.2 and -I 12 share a file.",
                       "\n  Refusing to adopt a partition from a different inflation."),
                 basename(CLS), got, INFLATION), call. = FALSE)
  }
  say("inflation sidecar checks out: -I ", got)
} else {
  warning("no .inflation sidecar beside ", basename(CLS),
          " -- cannot verify which inflation produced it", call. = FALSE)
}

# --- the partition ------------------------------------------------------------
# mcl cluster records run from the cluster index to a '$' and MAY WRAP across
# lines, so tokens are accumulated rather than lines parsed.
say("reading ", basename(CLS))
ln <- readLines(CLS, warn = FALSE)
i0 <- grep("^begin", ln)
if (!length(i0)) stop("no 'begin' in ", basename(CLS), call. = FALSE)
i1 <- grep("^\\)", ln)
i1 <- i1[i1 > i0][1]
tok <- unlist(strsplit(paste(ln[(i0 + 1):(i1 - 1)], collapse = " "), "[ \t]+"))
tok <- tok[nzchar(tok)]
ends <- which(tok == "$")
starts <- c(1L, head(ends, -1L) + 1L)
parts <- lapply(seq_along(ends), function(i) {
  m <- tok[(starts[i] + 1L):(ends[i] - 1L)]      # drop the leading cluster index
  if (!length(m)) integer(0) else as.integer(m)
})
say("  ", fmt_n(length(parts)), " clusters, ",
    fmt_n(sum(lengths(parts))), " node slots")

tab <- fread(TAB, header = FALSE, col.names = c("idx", "gene"))
say("  ", fmt_n(nrow(tab)), " genes in ", basename(TAB))
if (sum(lengths(parts)) != nrow(tab))
  stop("partition covers ", sum(lengths(parts)), " nodes but the tab has ",
       nrow(tab), " -- these are not the same graph", call. = FALSE)

# --- name the modules, largest first ------------------------------------------
ord <- order(lengths(parts), decreasing = TRUE)
parts <- parts[ord]
sizes <- lengths(parts)
keep  <- sizes >= MIN_SIZE
name  <- rep("Unassigned", length(parts))
name[keep] <- sprintf("Module_%03d", seq_len(sum(keep)))

memb <- data.table(idx = unlist(parts),
                   module_name = rep(name, sizes))
setkey(tab, idx)
memb <- tab[memb, on = "idx"]
if (anyNA(memb$gene)) stop("a partition index is absent from the tab file", call. = FALSE)

nm <- fread(NODES, select = c("gene", "degree", "strength"))
memb <- merge(memb[, .(gene, module_name)], nm, by = "gene", all.x = TRUE, sort = FALSE)
if (anyNA(memb$degree))
  stop(sum(is.na(memb$degree)), " genes have no node metrics -- ",
       basename(NODES), " is from a different network", call. = FALSE)

if (nrow(memb) != nrow(tab))
  stop("membership has ", nrow(memb), " rows for ", nrow(tab), " genes", call. = FALSE)
setcolorder(memb, c("gene", "module_name", "strength", "degree"))
write_tsv(memb, paste0(PREFIX, "_membership.tsv"))

# --- module summary, in the same contract -------------------------------------
Q <- NA_real_
if (nzchar(SWEEP) && file.exists(SWEEP)) {
  sw <- fread(SWEEP)
  hit <- sw[as.numeric(inflation) == INFLATION &
            (!"resource" %in% names(sw) | resource == "default")]
  if (nrow(hit) == 1) {
    Q <- as.numeric(hit$modularity)
    if (as.numeric(hit$n_clusters) != length(parts))
      stop("mcl_sweep says ", hit$n_clusters, " clusters at -I ", INFLATION,
           " but this file has ", length(parts), call. = FALSE)
    if (as.numeric(hit$largest) != max(sizes))
      stop("mcl_sweep says largest ", hit$largest, " but this file has ", max(sizes),
           call. = FALSE)
    say("cross-checked against ", basename(SWEEP), ": clusters and largest agree")
  } else say("no unique row for -I ", INFLATION, " in ", basename(SWEEP), "; Q = NA")
}

summ <- data.table(module = name[keep], n_genes = sizes[keep],
                   modularity_Q = Q, inflation = INFLATION,
                   min_module_size = MIN_SIZE,
                   shown_in_plot = sizes[keep] >= MIN_PLOT)
write_tsv(summ, paste0(PREFIX, "_module_summary.tsv"))

n_un <- sum(sizes[!keep])
say("")
say(STUDY, " at -I ", INFLATION, ":")
say("  clusters              ", fmt_n(length(parts)))
say("  named modules (>= ", MIN_SIZE, ")  ", fmt_n(sum(keep)))
say("  largest module        ", fmt_n(max(sizes)), sprintf(" (%.2f%%)", 100 * max(sizes) / nrow(tab)))
say("  median named size     ", median(sizes[keep]))
say("  Unassigned genes      ", fmt_n(n_un), sprintf(" (%.2f%%)", 100 * n_un / nrow(tab)))
say("  modularity Q          ", ifelse(is.na(Q), "NA", format(Q)))
say("  genes accounted for   ", fmt_n(nrow(memb)), " of ", fmt_n(nrow(tab)))
say("done: ", STUDY)
