library("ggplot2")
library("dplyr")
library("patchwork")
library("scales")

# ============================================================================
# PARAMETERS
# ============================================================================

file_sugarcane <- "/home/genomics/jorge/files/sugarcane/network_sugarcane_node_metrics.tsv"
file_purple    <- "/home/genomics/jorge/files/purple/new/network_purple_node_metrics.tsv"

panel_title_a  <- "A. Sugarcane network degree distribution"
panel_title_b  <- "B. S. officinarum and S. robustum network degree distribution"

output_dir     <- "/home/genomics/jorge/files/degree_distribution_panel"
output_file    <- file.path(output_dir, "degree_distribution_panel.png")

# Plot parameters
metric          <- "degree"
bin_width       <- NULL  # auto-detect via Sturges rule
fill_color_a    <- "black"  # viridis dark
fill_color_b    <- "black"  # viridis light
text_color      <- "#000000"

# ============================================================================
# LOAD AND PREPARE DATA
# ============================================================================

cat("Loading network node metrics...\n\n")

# Load files
data_sugarcane <- read.table(file_sugarcane, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
data_purple    <- read.table(file_purple, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cat(sprintf("Sugarcane: %d nodes\n", nrow(data_sugarcane)))
cat(sprintf("  %s — min: %.2f, max: %.2f, mean: %.2f\n",
            metric,
            min(data_sugarcane[[metric]]),
            max(data_sugarcane[[metric]]),
            mean(data_sugarcane[[metric]])))

cat(sprintf("Purple: %d nodes\n", nrow(data_purple)))
cat(sprintf("  %s — min: %.2f, max: %.2f, mean: %.2f\n",
            metric,
            min(data_purple[[metric]]),
            max(data_purple[[metric]]),
            mean(data_purple[[metric]])))

# ============================================================================
# BUILD PLOTS
# ============================================================================

build_degree_scatter <- function(data, metric_col, title, point_color) {

  # Count frequency of each degree value
  deg_dist <- data %>%
    group_by(.data[[metric_col]]) %>%
    summarise(N = n(), .groups = "drop") %>%
    arrange(.data[[metric_col]])

  colnames(deg_dist)[1] <- "degree"

  p <- ggplot(deg_dist, aes(x = degree, y = N)) +
    geom_point(alpha = 0.7, color = point_color, size = 2.5) +
    scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
    xlab("Degree (k)") +
    ylab("Frequency P(k)") +
    ggtitle(title) +
    theme_bw(base_size = 14) +
    theme(
      plot.title        = element_text(face = "bold", size = 16,
                                       hjust = 0.0, margin = margin(b = 8)),
      axis.title        = element_text(size = 13, face = "bold"),
      axis.text         = element_text(size = 11),
      panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      panel.grid.major  = element_line(colour = "grey90", linewidth = 0.3),
      panel.grid.minor  = element_blank(),
      plot.background   = element_rect(fill = "white", colour = NA)
    )

  p
}

# Build individual scatter plots
p1 <- build_degree_scatter(data_sugarcane, metric, panel_title_a, "black")
p2 <- build_degree_scatter(data_purple, metric, panel_title_b, "black")

# Assemble panel
cat("\nAssembling panel...\n")

panel <- p1 + p2 +
  plot_layout(ncol = 2, widths = c(1, 1))

# Save
cat("Saving panel...\n\n")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(
  filename = output_file,
  plot     = panel,
  width    = 40,
  height   = 15,
  units    = "cm",
  dpi      = 320,
  bg       = "white"
)

cat(sprintf("✓ Panel saved: %s\n\n", output_file))
cat("=================================================================\n")
