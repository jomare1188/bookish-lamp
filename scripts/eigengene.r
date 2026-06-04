library(data.table)
library(DESeq2)

# ==============================================================================
# CONFIGURATION
# ==============================================================================
Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1)

# ==============================================================================
# INPUT FILES
# ==============================================================================
DDS_PATH <- "/dados02/jorge/comparative_saccharum/run1/salmon/deseq2_qc/deseq2.dds.RData"
# /dados02/jorge/comparative_saccharum/china/run2_onlyL/salmon/deseq2_qc/deseq2.dds.RData
# /dados02/jorge/comparative_saccharum/run1/salmon/deseq2_qc/deseq2.dds.RData
networks <- list(
  sugarcane = list(
    membership_file = "/home/genomics/jorge/files/sugarcane/mcl_sugarcane_membership.tsv",
    col_filter      = function(dds) seq_len(ncol(dds)),
    out_dir         = "/home/genomics/jorge/files/sugarcane/"
  )
#  purple = list(
#    membership_file = "/home/genomics/jorge/files/purple/new/mcl_purple_membership.tsv",
#    col_filter      = function(dds) seq_len(ncol(dds)),
#    out_dir         = "/home/genomics/jorge/files/purple/new/"
#  )
)

# ==============================================================================
# LOAD VST DATA ONCE
# ==============================================================================
message("Loading DESeq2 object...")
load(DDS_PATH)   # loads 'dds'
all_vst <- as.matrix(assay(dds, "vst"))
message("VST matrix loaded: ", paste(dim(all_vst), collapse = " x "))

# ==============================================================================
# HELPER: compute eigengene (PC1) for a single module
#
# - Scales genes (unit variance) before PCA so high-variance genes don't
#   dominate — consistent with standard WGCNA eigengene definition.
# - Sign convention: flip PC1 so it correlates positively with the mean
#   expression of the module (avoids arbitrary sign ambiguity).
# - Returns both the PC1 scores (length = n_samples) and the variance
#   explained by PC1.
# ==============================================================================
compute_module_eigengene <- function(vst_sub) {

# Remove zero variance genes 
  vst_sub <- vst_sub[apply(vst_sub, 1, var, na.rm = TRUE) > 0, ]
  
  ###
  pca     <- prcomp(t(vst_sub), center = TRUE, scale. = TRUE)
  pc1     <- pca$x[, 1]
  var_pct <- summary(pca)$importance[2, 1] * 100   # % variance explained

  # Orient PC1 toward mean expression
  if (cor(pc1, colMeans(vst_sub)) < 0) pc1 <- -pc1

  list(pc1 = pc1, var_pct = var_pct)
}

# ==============================================================================
# MAIN: compute eigengenes for one network
# ==============================================================================
compute_network_eigengenes <- function(network_name, membership_file,
                                       col_filter, out_dir) {

  cat("\n", strrep("=", 60), "\n", sep = "")
  cat("Eigengenes: ", network_name, "\n", sep = "")
  cat(strrep("=", 60), "\n", sep = "")

  # ── 1. Sample selection (mirrors pearson_cor.r col_filter) ──────────
  col_idx    <- col_filter(dds)
  vst        <- all_vst[, col_idx, drop = FALSE]
  sample_ids <- colnames(vst)
  vst_genes  <- rownames(vst)
  cat(sprintf("  Samples: %d | Genes in VST: %s\n",
              length(sample_ids), format(length(vst_genes), big.mark = ",")))

  # ── 2. Load module membership ────────────────────────────────────────
  membership <- fread(membership_file, header = TRUE)
  setnames(membership, tolower(names(membership)))   # defensive normalisation
  modules    <- sort(unique(membership[module_name != "Unassigned", module_name]))
  cat(sprintf("  Modules (excl. Unassigned): %d\n", length(modules)))

  # ── 3. Eigengene per module ──────────────────────────────────────────
  eigen_list   <- vector("list", length(modules))
  var_pct_list <- numeric(length(modules))
  gene_n_list  <- integer(length(modules))
  names(eigen_list) <- modules

  for (i in seq_along(modules)) {
    mod       <- modules[i]
    mod_genes <- membership[module_name == mod, gene]
    mod_genes <- intersect(mod_genes, vst_genes)

    if (length(mod_genes) < 2L) {
      warning("Module ", mod, " has < 2 genes in VST matrix — skipped.")
      next
    }

    res               <- compute_module_eigengene(vst[mod_genes, , drop = FALSE])
    eigen_list[[i]]   <- res$pc1
    var_pct_list[i]   <- res$var_pct
    gene_n_list[i]    <- length(mod_genes)

    cat(sprintf("  %-15s  %5d genes  PC1 var: %5.1f%%\n",
                mod, length(mod_genes), res$var_pct))
  }

  # Drop any modules that were skipped (< 2 genes)
  keep        <- !vapply(eigen_list, is.null, logical(1L))
  eigen_list  <- eigen_list[keep]
  var_pct_list <- var_pct_list[keep]
  gene_n_list  <- gene_n_list[keep]

  # ── 4. Assemble samples × modules matrix ────────────────────────────
  eigen_mat        <- do.call(cbind, eigen_list)   # samples × modules
  rownames(eigen_mat) <- sample_ids

  eigen_dt <- data.table(sample_id = sample_ids, as.data.table(eigen_mat))

  # ── 5. Variance-explained summary ────────────────────────────────────
  var_dt <- data.table(
    module        = names(eigen_list),
    n_genes       = gene_n_list,
    pc1_var_pct   = round(var_pct_list, 2)
  )

  # ── 6. Save outputs ──────────────────────────────────────────────────
  out_eigen <- file.path(out_dir, paste0("eigengenes_", network_name, ".tsv"))
  out_var   <- file.path(out_dir, paste0("eigengenes_", network_name, "_pc1_variance.tsv"))

  fwrite(eigen_dt, file = out_eigen, sep = "\t", quote = FALSE)
  fwrite(var_dt,   file = out_var,   sep = "\t", quote = FALSE)

  cat(sprintf("\n  Saved: %s  [%d samples × %d modules]\n",
              basename(out_eigen), nrow(eigen_mat), ncol(eigen_mat)))
  cat(sprintf("  Saved: %s\n", basename(out_var)))

  invisible(list(eigengenes = eigen_dt, variance = var_dt))
}

# ==============================================================================
# ENTRY POINT
# ==============================================================================
for (net in names(networks)) {
  n <- networks[[net]]
  compute_network_eigengenes(
    network_name    = net,
    membership_file = n$membership_file,
    col_filter      = n$col_filter,
    out_dir         = n$out_dir
  )
}

message("\nAll eigengenes computed.")
