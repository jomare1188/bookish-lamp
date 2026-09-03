# =============================================================================
# 27_sbm_membership.r — a stochastic block model fit, in the MCL contract
#
# The pipeline reads its clustering through exactly two files, and every stage
# downstream of them — eigengenes, module-trait Spearman, TF enrichment, module
# GO, the heatmaps and the topology figure — takes them as given. This script
# writes those two files from a graph-tool nested SBM fit, so an alternative
# clustering can be carried through the identical downstream analysis and the
# two compared on their biology rather than on their partition shape.
#
# THE CONTRACT, and it is honoured exactly:
#   <prefix>_membership.tsv      gene, module_name, strength, degree
#   <prefix>_module_summary.tsv  module, n_genes, modularity_Q,
#                                modularity_Q_weighted, inflation,
#                                min_module_size, shown_in_plot
# `module_name` is Module_%03d ranked LARGEST FIRST, or the literal "Unassigned"
# for genes in blocks below MIN_SIZE -- both conventions copied from
# 05_mcl_clustering.r, because 22_fig_topology.r separates real modules from the
# leftover pool with `grepl("^Module_", module)` and 14_module_eigengene.r drops
# "Unassigned" by name.
#
# WHY NO EDGE TABLE IS READ FOR strength AND degree. MCL's own `strength` and
# `degree` columns are whole-graph igraph values, not intramodular ones
# (05_mcl_clustering.r:121-122) -- which makes them identical to the columns
# already sitting in network_<study>_node_metrics.tsv. So they are read from
# there. Only one consumer uses `strength` at all: 16_module_heatmaps.r, to take
# the top N genes of a module too large to draw whole.
#
# MODULARITY IS COMPUTED, and it is the one slow step (CLEAN_SBM_COMPUTE_Q=0
# skips it and writes NA). An SBM does not optimise modularity -- it minimises a
# description length -- so a low Q here is not a defect, it is the two methods
# optimising different things. It is computed anyway because it is the number
# 05_mcl_clustering.r reports for MCL, and a comparison needs both sides measured
# the same way.
#
# THE LEVEL. A nested fit gives a hierarchy; CLEAN_SBM_LEVEL picks which level is
# "the modules". Level 0 is the finest and the only one at module granularity.
# The level's block ids are taken from the b_l<N> column of the node-block table,
# which is a direct gene -> block map and needs no parsing.
#
# RUN: through run.sh  ->  ./run.sh sbmclust sugarcane
# =============================================================================

suppressMessages({ library(data.table); library(igraph) })

source(file.path(dirname(sub("--file=", "",
       grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "lib", "common.R"))

STUDY      <- env_req("CLEAN_STUDY")
NODE_BLOCKS<- env_req("CLEAN_SBM_NODE_BLOCKS")
NODES      <- env_req("CLEAN_NODE_METRICS")
EDGES      <- env_req("CLEAN_EDGES")
PREFIX     <- env_req("CLEAN_PREFIX")
LEVEL      <- as.integer(env_num("CLEAN_SBM_LEVEL", 0))
MIN_SIZE   <- as.integer(env_num("CLEAN_MIN_MODULE_SIZE", 2))
MIN_PLOT   <- as.integer(env_num("CLEAN_MIN_MODULE_SIZE_PLOT", 10))
COMPUTE_Q  <- env_flag("CLEAN_SBM_COMPUTE_Q", TRUE)
setDTthreads(as.integer(env_num("CLEAN_CORES", 100)))

banner(paste0("SBM clustering -> MCL contract: ", STUDY, "  (level ", LEVEL, ")"))

for (f in c(NODE_BLOCKS, NODES))
  if (!file.exists(f)) stop("missing ", basename(f), call. = FALSE)

# --- the partition -----------------------------------------------------------
lev_col <- sprintf("b_l%d", LEVEL)
hdr <- names(fread(NODE_BLOCKS, nrows = 0L))
if (!lev_col %in% hdr)
  stop("level ", LEVEL, " is not in ", basename(NODE_BLOCKS), "\n",
       "  available: ", paste(grep("^b_l", hdr, value = TRUE), collapse = ", "),
       call. = FALSE)
nb <- fread(NODE_BLOCKS, select = c("name", lev_col))
setnames(nb, c("gene", "block"))
say("SBM nodes at level ", LEVEL, ": ", fmt_n(nrow(nb)),
    " in ", fmt_n(uniqueN(nb$block)), " blocks")

# --- strength and degree, from the network the SBM was fitted on -------------
nm <- fread(NODES, select = c("gene", "degree", "strength"))
say("network nodes: ", fmt_n(nrow(nm)))

# The two must describe the SAME graph. If they do not, every module downstream
# is a different object from the one the fit produced, and the comparison this
# script exists to enable would be meaningless -- so it is fatal, not a warning.
only_sbm <- setdiff(nb$gene, nm$gene)
only_net <- setdiff(nm$gene, nb$gene)
if (length(only_sbm) || length(only_net))
  stop("the SBM fit and the network do not describe the same gene set:\n",
       "  only in the SBM fit: ", fmt_n(length(only_sbm)),
       if (length(only_sbm)) paste0(" e.g. ", paste(head(only_sbm, 3), collapse = ", ")) else "",
       "\n  only in the network: ", fmt_n(length(only_net)),
       if (length(only_net)) paste0(" e.g. ", paste(head(only_net, 3), collapse = ", ")) else "",
       "\n  the fit must be run on this study's network_", STUDY, "_edges.tsv",
       call. = FALSE)
say("gene sets identical (", fmt_n(nrow(nb)), " genes) — the fit is on this network")

# --- name the blocks, largest first, exactly as MCL does ---------------------
mem <- merge(nb, nm, by = "gene")
sizes      <- mem[, .N, by = block][order(-N)]
large      <- sizes[N >= MIN_SIZE, block]
rank_map   <- data.table(block = large,
                         module_name = sprintf("Module_%03d", seq_along(large)))
mem <- merge(mem, rank_map, by = "block", all.x = TRUE)
mem[is.na(module_name), module_name := "Unassigned"]
mem[, block := NULL]
setorder(mem, module_name, -strength)
say(sprintf("blocks >= %d genes: %s named | %s genes unassigned",
            MIN_SIZE, fmt_n(length(large)), fmt_n(mem[module_name == "Unassigned", .N])))

# --- modularity of the SBM partition, on the same graph ----------------------
q <- NA_real_
qw <- NA_real_
if (COMPUTE_Q) {
  assert_network_schema(EDGES)
  say("reading edges for modularity — the slow step (CLEAN_SBM_COMPUTE_Q=0 skips it)")
  e <- fread(EDGES, header = TRUE, select = c("gene1", "gene2", "weight"))
  say("  ", fmt_n(nrow(e)), " edges")
  g <- graph_from_data_frame(e, directed = FALSE)
  rm(e); invisible(gc())
  g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE,
                edge.attr.comb = list(weight = "max"))
  # Modularity needs an integer community per vertex in the graph's own order.
  # "Unassigned" genes each become their own community, which is what MCL's
  # partition does for them too -- they were never merged into anything.
  key <- mem[match(V(g)$name, gene)]
  comm <- ifelse(key$module_name == "Unassigned", NA_character_, key$module_name)
  comm <- ifelse(is.na(comm), paste0("solo_", V(g)$name), comm)
  # BOTH, and this is a correction. This stage used to compute only the WEIGHTED
  # modularity while 05_mcl_clustering.r computed only the UNWEIGHTED one --
  # igraph does not pick up the `weight` attribute when `weights` is NULL -- so
  # the headline "MCL 0.1005 vs SBM 0.0081" was a comparison between two
  # different statistics, asserted as like-for-like in four places in the docs.
  cm <- as.integer(factor(comm))
  q  <- modularity(g, membership = cm)
  qw <- modularity(g, membership = cm, weights = E(g)$weight)
  say(sprintf("modularity Q of the SBM partition: %.4f (unweighted), %.4f (weighted)",
              q, qw))
  rm(g); invisible(gc())
} else {
  say("NOTE: CLEAN_SBM_COMPUTE_Q=0 — modularity_Q written as NA")
}

# --- write the contract ------------------------------------------------------
write_tsv(mem[, .(gene, module_name, strength, degree)],
          paste0(PREFIX, "_membership.tsv"))

summ <- mem[, .N, by = module_name][order(-N)]
setnames(summ, c("module", "n_genes"))
summ[, modularity_Q    := if (is.na(q)) NA_real_ else round(q, 4)]
summ[, modularity_Q_weighted := if (is.na(qw)) NA_real_ else round(qw, 4)]
# `inflation` is an MCL parameter with no SBM equivalent; the column is kept so
# the schema matches and carries the level instead, which is the analogous knob.
summ[, inflation       := NA_real_]
summ[, sbm_level       := LEVEL]
summ[, min_module_size := MIN_SIZE]
summ[, shown_in_plot   := n_genes >= MIN_PLOT & module != "Unassigned"]
write_tsv(summ[, .(module, n_genes, modularity_Q, modularity_Q_weighted, inflation, sbm_level,
                   min_module_size, shown_in_plot)],
          paste0(PREFIX, "_module_summary.tsv"))

named <- summ[module != "Unassigned"]
say("")
say(sprintf("modules (>= %d genes): %s | median %s | largest %s | shown in plot (>= %d): %s",
            MIN_SIZE, fmt_n(nrow(named)), fmt_n(median(named$n_genes)),
            fmt_n(max(named$n_genes)), MIN_PLOT, fmt_n(sum(summ$shown_in_plot))))
say("done: ", STUDY)
