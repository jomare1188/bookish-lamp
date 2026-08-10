library(data.table)
library(igraph)
library(ggplot2)
library(hexbin)
library(scales)

# ==============================================================================
# CONFIGURATION — adjust these parameters as needed
# ==============================================================================
PEARSON_MIN   <- 0.8     # minimum |r| to keep an edge
PEARSON_MAX   <- 0.9999  # maximum |r| — excludes rounding artefacts / homoeologs
PADJ_FILTER   <- 0.05    # maximum adjusted p-value (set NULL to skip)

# Normalisation method for edge weights loaded into igraph:
#   "minmax"  — rescales |r| to [0, 1] within the filtered edge set
#   "rank"    — replaces |r| with its rank / max_rank (robust to outliers)
#   "raw"     — uses the raw Pearson r value as-is (no normalisation)
WEIGHT_METHOD <- "minmax"

# Output directory
OUT_DIR <- "/home/genomics/jorge/files/purple/"

# ==============================================================================
# INPUT FILES
# ==============================================================================
files_to_process <- list(
  responsive = list(
    edge_file     = "/home/genomics/jorge/files/purple/edgelist_purple_responsive_pearson.tsv",
    output_prefix = file.path(OUT_DIR, "network_purple_responsive")
  ),
  non_responsive = list(
    edge_file     = "/home/genomics/jorge/files/purple/edgelist_purple_non_responsive_pearson.tsv",
    output_prefix = file.path(OUT_DIR, "network_purple_non_responsive")
  )
)

# ==============================================================================
# HELPER: normalise a numeric vector of absolute correlations
# ==============================================================================
normalise_weights <- function(abs_r, method) {
  switch(method,
    minmax = {
      lo <- min(abs_r, na.rm = TRUE)
      hi <- max(abs_r, na.rm = TRUE)
      if (hi == lo) return(rep(1, length(abs_r)))
      (abs_r - lo) / (hi - lo)
    },
    rank = {
      r <- rank(abs_r, ties.method = "average", na.last = "keep")
      r / max(r, na.rm = TRUE)
    },
    raw = abs_r,
    stop("Unknown WEIGHT_METHOD: '", method,
         "'. Choose 'minmax', 'rank', or 'raw'.")
  )
}

# ==============================================================================
# MAIN ANALYSIS FUNCTION
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
  cat(sprintf("Filters:  |r| in [%.4f, %.4f]", pearson_min, pearson_max))
  if (!is.null(padj_filter)) cat(sprintf(" | padj <= %.4f", padj_filter))
  cat(sprintf("\nWeights:  %s normalisation of |r|\n", weight_method))

  # ── 1. Read edge list ───────────────────────────────────────────────
  cat("Reading edge list...\n")
  edges <- fread(edge_file, header = TRUE,
                 col.names = c("gene1", "gene2", "pearson", "pval", "padj"))
  cat(sprintf("  Raw edges: %s\n", format(nrow(edges), big.mark = ",")))

  # ── 2. Filter ───────────────────────────────────────────────────────
  abs_r <- abs(edges$pearson)
  keep  <- abs_r >= pearson_min & abs_r <= pearson_max

  if (!is.null(padj_filter)) {
    keep <- keep & edges$padj <= padj_filter
  }

  edges <- edges[keep]
  cat(sprintf("  After filters: %s edges\n", format(nrow(edges), big.mark = ",")))

  if (nrow(edges) == 0L) {
    cat("  No edges passed filters — skipping.\n")
    return(invisible(NULL))
  }

  # ── 3. Normalise weights ────────────────────────────────────────────
  edges[, weight := normalise_weights(abs(pearson), weight_method)]

  cat(sprintf("  Weight range after normalisation: [%.6f, %.6f]\n",
              min(edges$weight), max(edges$weight)))

  # ── 4. Save filtered + normalised edge list ─────────────────────────
  out_edges <- paste0(output_prefix, "_filtered_edges.tsv")
  fwrite(edges, file = out_edges, sep = "\t", quote = FALSE)
  cat("  Saved filtered edge list: ", basename(out_edges), "\n", sep = "")

  # ── 5. Build weighted igraph ────────────────────────────────────────
  cat("Building weighted igraph object...\n")
  g <- graph_from_data_frame(
    edges[, .(gene1, gene2, weight)],
    directed = FALSE
  )
  rm(edges); gc()

  cat(sprintf("  Nodes: %s | Edges: %s\n",
              format(vcount(g), big.mark = ","),
              format(ecount(g), big.mark = ",")))

  # ── 6. Global metrics ───────────────────────────────────────────────
  cat("Calculating global metrics...\n")

  global_trans   <- transitivity(g, type = "global")
  mean_w         <- mean(E(g)$weight)
  density        <- edge_density(g)
  n_components   <- components(g)$no
  giant_size     <- max(components(g)$csize)

  global_metrics <- data.table(
    Metric = c("Nodes", "Edges", "Edge_density",
               "Global_transitivity", "Mean_normalised_weight",
               "N_connected_components", "Giant_component_size",
               "Pearson_min_filter", "Pearson_max_filter",
               "Padj_filter", "Weight_method"),
    Value  = c(vcount(g), ecount(g),
               round(density, 6),
               round(global_trans, 6),
               round(mean_w, 6),
               n_components, giant_size,
               pearson_min, pearson_max,
               ifelse(is.null(padj_filter), "none", padj_filter),
               weight_method)
  )

  out_global <- paste0(output_prefix, "_global_metrics.tsv")
  fwrite(global_metrics, file = out_global, sep = "\t", quote = FALSE)
  cat("  Saved global metrics: ", basename(out_global), "\n", sep = "")
  print(global_metrics, row.names = FALSE)

  # ── 7. Node-level metrics ───────────────────────────────────────────
  cat("Calculating node-level metrics...\n")

  node_degree  <- degree(g, mode = "all")
  node_str     <- strength(g, mode = "all")       # weighted degree
  node_trans   <- transitivity(g, type = "local")

  node_stats <- data.table(
    gene         = V(g)$name,
    degree       = node_degree,
    strength     = round(node_str, 6),            # sum of normalised weights
    transitivity = node_trans
  )
  node_stats[is.nan(transitivity), transitivity := NA]

  out_nodes <- paste0(output_prefix, "_node_metrics.tsv")
  fwrite(node_stats, file = out_nodes, sep = "\t", quote = FALSE)
  cat("  Saved node metrics: ", basename(out_nodes), "\n", sep = "")

  # ── 8. Plots ────────────────────────────────────────────────────────
  cat("Generating plots...\n")

  pub_theme <- theme_classic(base_size = 14) +
    theme(plot.title   = element_text(face = "bold", hjust = 0.5),
          axis.text    = element_text(color = "black"))

  # Plot A: Degree distribution (log-log)
  deg_dist <- node_stats[, .N, by = degree][order(degree)]

  p1 <- ggplot(deg_dist, aes(x = degree, y = N)) +
    geom_point(alpha = 0.6, color = "#2c3e50", size = 1.5) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    labs(title = paste("Degree Distribution —", grp),
         x = "Degree (k)", y = "Frequency P(k)") +
    pub_theme

  ggsave(paste0(output_prefix, "_degree_distribution.pdf"),
         plot = p1, width = 6, height = 5, device = "pdf")
  ggsave(paste0(output_prefix, "_degree_distribution.png"),
         plot = p1, width = 6, height = 5, dpi = 300)

  # Plot B: Strength distribution (log-log) — weighted analogue of degree
  str_dist <- node_stats[, .(mean_strength = mean(strength)), by = degree][order(degree)]

  p2 <- ggplot(str_dist, aes(x = degree, y = mean_strength)) +
    geom_point(alpha = 0.6, color = "#8e44ad", size = 1.5) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    labs(title = paste("Mean Strength vs Degree —", grp),
         x = "Degree (k)", y = "Mean Strength s(k)") +
    pub_theme

  ggsave(paste0(output_prefix, "_strength_vs_degree.pdf"),
         plot = p2, width = 6, height = 5, device = "pdf")
  ggsave(paste0(output_prefix, "_strength_vs_degree.png"),
         plot = p2, width = 6, height = 5, dpi = 300)

  # Plot C: Transitivity vs degree (hexbin to avoid overplotting)
  p3 <- ggplot(node_stats[!is.na(transitivity)],
               aes(x = degree, y = transitivity)) +
    geom_hex(bins = 75) +
    scale_fill_viridis_c(trans = "log10", name = "Node Count") +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    labs(title = paste("Transitivity vs Degree —", grp),
         x = "Degree (k)", y = "Local Transitivity C(k)") +
    pub_theme

  ggsave(paste0(output_prefix, "_transitivity_vs_degree.pdf"),
         plot = p3, width = 7, height = 5, device = "pdf")
  ggsave(paste0(output_prefix, "_transitivity_vs_degree.png"),
         plot = p3, width = 7, height = 5, dpi = 300)

  cat("Analysis complete for: ", grp, "\n", sep = "")

  rm(g, node_stats, deg_dist, str_dist, p1, p2, p3); gc()
}

# ==============================================================================
# ENTRY POINT
# ==============================================================================
for (grp in names(files_to_process)) {
  f <- files_to_process[[grp]]
  analyze_network(
    edge_file     = f$edge_file,
    output_prefix = f$output_prefix
  )
}

message("\nAll networks processed.")#
