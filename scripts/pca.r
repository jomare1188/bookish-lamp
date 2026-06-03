library("ggplot2")
library("dplyr")
library("DESeq2")
library("tximport")
library("tidyverse")
library("viridis")
library("patchwork")

# ============================================================================
# PARAMETERS
# ============================================================================

# --- Study 1: Sugarcane ---
base_dir_s1        <- "/dados02/jorge/comparative_saccharum"
metadata_file_s1   <- file.path(base_dir_s1, "samplesheet.csv")
salmon_subdir_s1   <- "run1/salmon"
panel_title_s1     <- "Sugarcane"

# --- Study 2: S. officinarum and S. robustum ---
base_dir_s2        <- "/dados02/jorge/comparative_saccharum/china"
metadata_file_s2   <- file.path(base_dir_s2, "samplesheet_china.csv")
salmon_subdir_s2   <- "run2_onlyL/salmon"
panel_title_s2     <- expression(italic("S. officinarum") ~ "and" ~ italic("S. robustum"))

# --- Shared settings ---
grouping_col       <- "genotype"   # column used for colouring PCA
cv_threshold       <- 0.15         # filter OUT genes with CV < 15%
pca_top_genes      <- 99990000        # top variable genes for PCA

# --- Output ---
output_dir         <- "/dados02/jorge/comparative_saccharum/pca_panel"
output_file        <- file.path(output_dir, "pca_panel_saccharum.png")

# ============================================================================
# HELPERS
# ============================================================================

load_study <- function(base_dir, metadata_file, salmon_subdir, grouping_col) {

  cat(sprintf("\n── Loading: %s\n", metadata_file))

  metadata <- read.csv(metadata_file, header = TRUE, stringsAsFactors = FALSE)

  if (!grouping_col %in% colnames(metadata))
    stop(sprintf("Column '%s' not found in %s.\nAvailable columns: %s",
                 grouping_col, metadata_file,
                 paste(colnames(metadata), collapse = ", ")))

  # Determine the sample-name column (first column, or 'sample')
  sample_col <- if ("sample" %in% colnames(metadata)) "sample" else colnames(metadata)[1]
  cat(sprintf("  Using '%s' as sample identifier column\n", sample_col))

  sample_files <- file.path(base_dir, salmon_subdir, metadata[[sample_col]], "quant.sf")
  names(sample_files) <- metadata[[sample_col]]

  missing <- sample_files[!file.exists(sample_files)]
  if (length(missing) > 0) {
    cat("  ⚠ Missing quant.sf files:\n")
    for (f in missing) cat(sprintf("    %s\n", f))
    stop("Aborting: one or more quant.sf files not found.")
  }

  # Auto-detect tx2gene — common naming conventions
  tx2gene_candidates <- c(
    file.path(base_dir, salmon_subdir, "salmon.merged.tx2gene.tsv"),
    file.path(base_dir, salmon_subdir, "tx2gene.tsv"),
    file.path(base_dir, "tx2gene.tsv")
  )
  tx2gene_file <- tx2gene_candidates[file.exists(tx2gene_candidates)][1]
  if (is.na(tx2gene_file))
    stop(sprintf("No tx2gene file found. Tried:\n%s",
                 paste(tx2gene_candidates, collapse = "\n")))
  cat(sprintf("  tx2gene: %s\n", tx2gene_file))

  tx2gene <- read.table(tx2gene_file, header = TRUE)

  txi <- tximport(
    files            = sample_files,
    type             = "salmon",
    tx2gene          = tx2gene,
    ignoreTxVersion  = FALSE,
    ignoreAfterBar   = TRUE
  )

  coldata <- metadata
  rownames(coldata) <- coldata[[sample_col]]
  coldata[[grouping_col]] <- as.factor(coldata[[grouping_col]])

  stopifnot(all(colnames(txi$counts) == rownames(coldata)))

  cat(sprintf("  ✓ Loaded  |  Samples: %d  |  Genes: %d\n",
              ncol(txi$counts), nrow(txi$counts)))
  cat(sprintf("  Genotype groups: %s\n",
              paste(levels(coldata[[grouping_col]]), collapse = ", ")))

  list(txi = txi, coldata = coldata)
}


filter_by_cv <- function(counts_mat, cv_threshold) {
  # Coefficient of Variation = sd / mean  (per gene, across samples)
  gene_means <- rowMeans(counts_mat)
  gene_sds   <- apply(counts_mat, 1, sd)

  # Avoid division by zero: genes with mean == 0 get CV = 0
  cv <- ifelse(gene_means == 0, 0, gene_sds / gene_means)

  keep <- cv >= cv_threshold
  cat(sprintf("  CV filter (>= %.0f%%): kept %d / %d genes (removed %d)\n",
              cv_threshold * 100, sum(keep), length(keep), sum(!keep)))
  counts_mat[keep, , drop = FALSE]
}


build_pca_plot <- function(study_data, grouping_col, title,
                           ntop = 10000, palette_colors = NULL) {

  txi     <- study_data$txi
  coldata <- study_data$coldata

  # --- Build minimal DESeq2 object for VST ---
  formula_str <- as.formula(paste("~", grouping_col))

  dds <- DESeqDataSetFromTximport(
    txi     = txi,
    colData = coldata,
    design  = formula_str
  )

  # --- CV-based gene filtering on raw counts ---
  cat("  Applying CV filter...\n")
  raw_counts_filtered <- filter_by_cv(counts(dds), cv_threshold)

  # Subset the dds to CV-passing genes
  dds_filtered <- dds[rownames(raw_counts_filtered), ]

  # Estimate size factors (needed before VST)
  dds_filtered <- estimateSizeFactors(dds_filtered)

  # --- VST ---
  vst <- varianceStabilizingTransformation(dds_filtered, blind = TRUE)

  # --- PCA ---
  pca_data    <- plotPCA(vst, intgroup = grouping_col,
                         returnData = TRUE, ntop = ntop)
  percent_var <- round(100 * attr(pca_data, "percentVar"))

  n_groups <- nlevels(coldata[[grouping_col]])
  if (is.null(palette_colors)) {
    palette_colors <- viridis(n_groups, option = "D")
  }

  p <- ggplot(pca_data, aes(x = PC1, y = PC2,
                             color = .data[[grouping_col]])) +
    geom_point(size = 3.5, alpha = 0.9) +
    xlab(paste0("PC1: ", percent_var[1], "% variance")) +
    ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    ggtitle(title) +
    scale_colour_viridis_d(option = "D", name = "Genotype") +
    theme_bw(base_size = 16) +
    theme(
      plot.title        = element_text(face = "bold", size = 18*1.5,
                                       hjust = 0.5, margin = margin(b = 10)),
      legend.title      = element_text(face = "bold", size = 14*1.5),
      legend.text       = element_text(size = 12*1.5),
      legend.position   = "right",
      panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      panel.grid.major  = element_line(colour = "grey90", linewidth = 0.4),
      panel.grid.minor  = element_blank(),
      axis.title        = element_text(size = 14*1.5),
      axis.text         = element_text(size = 12*1.5),
      plot.background   = element_rect(fill = "white", colour = NA)
    )

  p
}

# ============================================================================
# MAIN
# ============================================================================

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("=================================================================\n")
cat("PCA PANEL: SACCHARUM STUDIES\n")
cat("=================================================================\n")

# Load both studies
cat("\n[1/4] Loading Study 1 (Sugarcane)...\n")
study1 <- load_study(base_dir_s1, metadata_file_s1, salmon_subdir_s1, grouping_col)

cat("\n[2/4] Loading Study 2 (S. officinarum / S. robustum)...\n")
study2 <- load_study(base_dir_s2, metadata_file_s2, salmon_subdir_s2, grouping_col)

# Build individual PCA plots
cat("\n[3/4] Building PCA plots...\n")

cat("\n  ── Study 1 PCA:\n")
p1 <- build_pca_plot(study1, grouping_col, panel_title_s1, ntop = pca_top_genes)

cat("\n  ── Study 2 PCA:\n")
p2 <- build_pca_plot(study2, grouping_col, panel_title_s2, ntop = pca_top_genes)

# Assemble panel
cat("\n[4/4] Assembling and saving panel...\n")

panel <- p1 + p2 +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag         = element_text(size = 20, face = "bold"),
      plot.background  = element_rect(fill = "white", colour = NA)
    )
  ) +
  plot_layout(ncol = 2, widths = c(1, 1))

ggsave(
  filename = output_file,
  plot     = panel,
  width    = 40,        # cm — wide format for poster
  height   = 18,
  units    = "cm",
  dpi      = 320,
  bg       = "white"
)

cat(sprintf("\n✓ Panel saved: %s\n", output_file))
cat("\n=================================================================\n")
cat("DONE\n")
cat("=================================================================\n")
