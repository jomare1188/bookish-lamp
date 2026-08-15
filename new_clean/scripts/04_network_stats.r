# =============================================================================
# 04_network_stats.r — node and global topology of one network
#
# Reads the merged network (03_merge_layers.py output) and writes:
#
#   <prefix>_node_metrics.tsv    gene, degree, strength, transitivity
#   <prefix>_global_metrics.tsv  nodes, edges, density, components, ...
#   <prefix>_degree_distribution.{pdf,png}
#   <prefix>_strength_vs_degree.{pdf,png}
#   <prefix>_transitivity_vs_degree.{pdf,png}    (only if transitivity is on)
#
# node_metrics is the most load-bearing file in the pipeline: it is the
# BACKGROUND UNIVERSE for the GO enrichment and for the degree-matched
# permutation nulls in the MYB61 and Module-20 readouts. A stale copy silently
# biases those tests rather than failing, which is why it is regenerated from
# the network rather than carried over.
#
# NOTE ON WHAT THIS NO LONGER DOES: the old general_stats.r both FILTERED the
# raw edge list (|r| in [0.8, 0.9999], padj <= 0.05) and computed topology, and
# wrote a filtered edge file as a side effect. The filtering now happens in the
# engine, so the network arrives already filtered and weighted and this stage
# only measures it. It never writes an edge file, so there is exactly one edge
# table per study and no way for a "filtered" and an "augmented" copy to drift
# apart -- which is precisely what went wrong in the old tree.
#
# RUN: through run.sh  ->  ./run.sh stats sugarcane
# =============================================================================

suppressMessages({
  library(data.table)
  library(igraph)
  library(ggplot2)
  library(hexbin)
  library(scales)
})

source(file.path(dirname(sub("--file=", "",
       grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "lib", "common.R"))

STUDY      <- env_req("CLEAN_STUDY")
EDGE_FILE  <- env_req("CLEAN_EDGES")
PREFIX     <- env_req("CLEAN_PREFIX")
# Local transitivity is O(sum deg^2) and took ~15 h on purple's 710 M edges.
# Nothing in the pipeline TESTS it -- it appears only as a reported column in
# the TF and MYB61 readouts -- so it is off by default and the column is NA.
TRANSITIVITY <- env_flag("CLEAN_TRANSITIVITY", FALSE)
setDTthreads(as.integer(env_num("CLEAN_CORES", 100)))

banner(paste("network topology:", STUDY))
assert_network_schema(EDGE_FILE)

# --- read --------------------------------------------------------------------
# Only three columns are needed, and restricting them is what keeps purple's
# 65 GB edge file inside RAM.
say("reading ", basename(EDGE_FILE))
edges <- fread(EDGE_FILE, header = TRUE,
               select = c("gene1", "gene2", "weight"))
say("  ", fmt_n(nrow(edges)), " edges")
if (nrow(edges) == 0L) stop("the network is empty", call. = FALSE)

# --- provenance: how many edges came from which layer? -----------------------
# This is the number the whole augmentation exists to produce, so report it
# here where it is cheap rather than making someone awk for it later.
src <- fread(EDGE_FILE, header = TRUE, select = "source")[, .N, by = source]
setorder(src, -N)
say("edge provenance:")
for (i in seq_len(nrow(src)))
  say(sprintf("  %-8s %15s  %6.2f%%", src$source[i], fmt_n(src$N[i]),
              100 * src$N[i] / nrow(edges)))
rm(src)

# --- graph -------------------------------------------------------------------
g <- graph_from_data_frame(edges, directed = FALSE)
rm(edges); invisible(gc())
say("nodes ", fmt_n(vcount(g)), " | edges ", fmt_n(ecount(g)))

# --- global metrics ----------------------------------------------------------
comp <- components(g)
global_metrics <- data.table(
  Metric = c("Study", "Nodes", "Edges", "Edge_density",
             "Global_transitivity", "Mean_normalised_weight",
             "N_connected_components", "Giant_component_size",
             "Edge_file"),
  Value  = c(STUDY, vcount(g), ecount(g),
             round(edge_density(g), 8),
             if (TRANSITIVITY) round(transitivity(g, "global"), 6) else NA,
             round(mean(E(g)$weight), 6),
             comp$no, max(comp$csize),
             basename(EDGE_FILE))
)
write_tsv(global_metrics, paste0(PREFIX, "_global_metrics.tsv"))
print(global_metrics, row.names = FALSE)

# --- node metrics ------------------------------------------------------------
say("computing degree and strength")
node_stats <- data.table(
  gene         = V(g)$name,
  degree       = degree(g),
  strength     = round(strength(g), 6),
  transitivity = NA_real_
)
if (TRANSITIVITY) {
  say("computing local transitivity (this is the slow one)")
  node_stats[, transitivity := transitivity(g, type = "local")]
  node_stats[is.nan(transitivity), transitivity := NA_real_]
} else {
  say("local transitivity SKIPPED (CLEAN_TRANSITIVITY=0); column is NA")
}
write_tsv(node_stats, paste0(PREFIX, "_node_metrics.tsv"))

# --- plots -------------------------------------------------------------------
pub_theme <- theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text  = element_text(color = "black"))

deg_dist <- node_stats[, .N, by = degree][order(degree)]
p1 <- ggplot(deg_dist, aes(x = degree, y = N)) +
  geom_point(alpha = 0.6, color = "#2c3e50", size = 1.5) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  labs(title = paste("Degree Distribution —", STUDY),
       x = "Degree (k)", y = "Frequency P(k)") + pub_theme
ggsave(paste0(PREFIX, "_degree_distribution.pdf"), p1, width = 6, height = 5)
ggsave(paste0(PREFIX, "_degree_distribution.png"), p1, width = 6, height = 5, dpi = 300)

str_dist <- node_stats[, .(mean_strength = mean(strength)), by = degree][order(degree)]
p2 <- ggplot(str_dist, aes(x = degree, y = mean_strength)) +
  geom_point(alpha = 0.6, color = "#8e44ad", size = 1.5) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  labs(title = paste("Mean Strength vs Degree —", STUDY),
       x = "Degree (k)", y = "Mean Strength s(k)") + pub_theme
ggsave(paste0(PREFIX, "_strength_vs_degree.pdf"), p2, width = 6, height = 5)
ggsave(paste0(PREFIX, "_strength_vs_degree.png"), p2, width = 6, height = 5, dpi = 300)

if (TRANSITIVITY && node_stats[!is.na(transitivity), .N] > 0L) {
  p3 <- ggplot(node_stats[!is.na(transitivity)], aes(x = degree, y = transitivity)) +
    geom_hex(bins = 75) +
    scale_fill_viridis_c(trans = "log10", name = "Node Count") +
    scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
    labs(title = paste("Transitivity vs Degree —", STUDY),
         x = "Degree (k)", y = "Local Transitivity C(k)") + pub_theme
  ggsave(paste0(PREFIX, "_transitivity_vs_degree.pdf"), p3, width = 7, height = 5)
  ggsave(paste0(PREFIX, "_transitivity_vs_degree.png"), p3, width = 7, height = 5, dpi = 300)
}

say("done: ", STUDY)
