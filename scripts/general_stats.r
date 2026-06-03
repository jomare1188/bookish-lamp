library(data.table)
library(igraph)
library(ggplot2)
library(hexbin)
library(scales)


# ==============================================================================
# CONFIGURATION
# ==============================================================================
PEARSON_MIN   <- 0.8
PEARSON_MAX   <- 0.9999
PADJ_FILTER   <- 0.05
WEIGHT_METHOD <- "minmax"



# ==============================================================================
# INPUT FILES — one entry per study
# ==============================================================================
files_to_process <- list(
#  purple = list(
#    edge_file     = "/home/genomics/jorge/files/purple/new/edgelist_purple_pearson.tsv",
#    output_prefix = file.path("/home/genomics/jorge/files/purple/new/", "network_purple")
  
#   sugarcane = list(
#     edge_file     = "/home/genomics/jorge/files/sugarcane/edgelist_sugarcane_pearson.tsv",
#     output_prefix = file.path("/home/genomics/jorge/files/sugarcane/", "network_sugarcane")
   )
)

# ==============================================================================
# HELPER: normalise weights
# ==============================================================================
normalise_weights <- function(abs_r, method) {
  switch(method,
   minmax = {
  lo <- min(abs_r, na.rm = TRUE); hi <- max(abs_r, na.rm = TRUE)
  if (hi == lo) return(rep(1, length(abs_r)))
  normalized <- (abs_r - lo) / (hi - lo)
  0.01 + (normalized * 0.99)
  },
    rank = {
      r <- rank(abs_r, ties.method = "average", na.last = "keep")
      r / max(r, na.rm = TRUE)
    },
    raw  = abs_r,
    stop("Unknown WEIGHT_METHOD: '", method, "'.")
  )
}

# ==============================================================================
# MAIN ANALYSIS
# ==============================================================================
analyze_network <- function(edge_file, output_prefix,
                            pearson_min = PEARSON_MIN,
                            pearson_max = PEARSON_MAX,
                            padj_filter = PADJ_FILTER,
                            weight_method = WEIGHT_METHOD) {

  grp <- basename(output_prefix)
  cat("\n", strrep("=", 60), "\n", sep = "")
  cat("Processing: ", grp, "\n", sep = "")
  cat(strrep("=", 60), "\n", sep = "")

  # ── Read ────────────────────────────────────────────────────────────
  edges <- fread(edge_file, header = TRUE,
                 col.names = c("gene1", "gene2", "pearson", "pval", "padj"))
  cat(sprintf("  Raw edges: %s\n", format(nrow(edges), big.mark = ",")))

  # ── Filter ──────────────────────────────────────────────────────────
  abs_r <- abs(edges$pearson)
  keep  <- abs_r >= pearson_min & abs_r <= pearson_max
  if (!is.null(padj_filter)) keep <- keep & edges$padj <= padj_filter
  edges <- edges[keep]
  cat(sprintf("  After filters: %s edges\n", format(nrow(edges), big.mark = ",")))

  if (nrow(edges) == 0L) { cat("  No edges — skipping.\n"); return(invisible(NULL)) }

  # ── Normalise ───────────────────────────────────────────────────────
  edges[, weight := normalise_weights(abs(pearson), weight_method)]

  out_edges <- paste0(output_prefix, "_filtered_edges.tsv")
  fwrite(edges, file = out_edges, sep = "\t", quote = FALSE)
  cat("  Saved filtered edges: ", basename(out_edges), "\n", sep = "")

  # ── Build graph ─────────────────────────────────────────────────────
  g <- graph_from_data_frame(edges[, .(gene1, gene2, weight)], directed = FALSE)
  rm(edges); gc()
  cat(sprintf("  Nodes: %s | Edges: %s\n",
              format(vcount(g), big.mark = ","),
              format(ecount(g), big.mark = ",")))

  # ── Global metrics ──────────────────────────────────────────────────
  comp <- components(g)
  global_metrics <- data.table(
    Metric = c("Nodes", "Edges", "Edge_density",
               "Global_transitivity", "Mean_normalised_weight",
               "N_connected_components", "Giant_component_size",
               "Pearson_min_filter", "Pearson_max_filter",
               "Padj_filter", "Weight_method"),
    Value  = c(vcount(g), ecount(g),
               round(edge_density(g),             6),
               round(transitivity(g, "global"),   6),
               round(mean(E(g)$weight),            6),
               comp$no, max(comp$csize),
               pearson_min, pearson_max,
               ifelse(is.null(padj_filter), "none", padj_filter),
               weight_method)
  )
  out_global <- paste0(output_prefix, "_global_metrics.tsv")
  fwrite(global_metrics, file = out_global, sep = "\t", quote = FALSE)
  cat("  Saved global metrics: ", basename(out_global), "\n", sep = "")
  print(global_metrics, row.names = FALSE)

  # ── Node metrics ────────────────────────────────────────────────────
  node_stats <- data.table(
    gene         = V(g)$name,
    degree       = degree(g),
    strength     = round(strength(g), 6),
    transitivity = transitivity(g, type = "local")
  )
  node_stats[is.nan(transitivity), transitivity := NA]

  out_nodes <- paste0(output_prefix, "_node_metrics.tsv")
  fwrite(node_stats, file = out_nodes, sep = "\t", quote = FALSE)
  cat("  Saved node metrics: ", basename(out_nodes), "\n", sep = "")

  # ── Plots ───────────────────────────────────────────────────────────
  pub_theme <- theme_classic(base_size = 14) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          axis.text  = element_text(color = "black"))

  deg_dist <- node_stats[, .N, by = degree][order(degree)]

  p1 <- ggplot(deg_dist, aes(x = degree, y = N)) +
    geom_point(alpha = 0.6, color = "#2c3e50", size = 1.5) +
    scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
    labs(title = paste("Degree Distribution —", grp),
         x = "Degree (k)", y = "Frequency P(k)") +
    pub_theme

  ggsave(paste0(output_prefix, "_degree_distribution.pdf"), p1, width = 6, height = 5)
  ggsave(paste0(output_prefix, "_degree_distribution.png"), p1, width = 6, height = 5, dpi = 300)

  str_dist <- node_stats[, .(mean_strength = mean(strength)), by = degree][order(degree)]

  p2 <- ggplot(str_dist, aes(x = degree, y = mean_strength)) +
    geom_point(alpha = 0.6, color = "#8e44ad", size = 1.5) +
    scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
    labs(title = paste("Mean Strength vs Degree —", grp),
         x = "Degree (k)", y = "Mean Strength s(k)") +
    pub_theme

  ggsave(paste0(output_prefix, "_strength_vs_degree.pdf"), p2, width = 6, height = 5)
  ggsave(paste0(output_prefix, "_strength_vs_degree.png"), p2, width = 6, height = 5, dpi = 300)

  p3 <- ggplot(node_stats[!is.na(transitivity)], aes(x = degree, y = transitivity)) +
    geom_hex(bins = 75) +
    scale_fill_viridis_c(trans = "log10", name = "Node Count") +
    scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
    labs(title = paste("Transitivity vs Degree —", grp),
         x = "Degree (k)", y = "Local Transitivity C(k)") +
    pub_theme

  ggsave(paste0(output_prefix, "_transitivity_vs_degree.pdf"), p3, width = 7, height = 5)
  ggsave(paste0(output_prefix, "_transitivity_vs_degree.png"), p3, width = 7, height = 5, dpi = 300)

  cat("Analysis complete for: ", grp, "\n", sep = "")
  rm(g, node_stats, deg_dist, str_dist, p1, p2, p3); gc()
}

# ==============================================================================
# ENTRY POINT
# ==============================================================================
for (grp in names(files_to_process)) {
  f <- files_to_process[[grp]]
  analyze_network(edge_file = f$edge_file, output_prefix = f$output_prefix)
}

message("\nAll networks processed.")
