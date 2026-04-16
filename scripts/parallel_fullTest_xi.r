library(XICOR)
library(DESeq2)
library(Matrix)
library(data.table)
library(matrixStats)
library(parallel)

# --- CONFIGURATION ---
#n_genes_test <- 1000
n_cores    <- 250
chunk_size <- 50
options(cores = n_cores)
closeAllConnections()

# --- DATA LOADING ---
load("deseq2.dds.RData")
all_vst_full <- Matrix::as.matrix(assay(dds, "vst"))

# --- PROCESSING FUNCTION ---
process_group_linux <- function(group_name, genotype_label, n_genes_test = NULL) {
  message("\n>>> Starting Group: ", group_name)

  # 1. Subset by genotype name and CV Filter
  col_names <- colnames(dds)[dds$Group1 == genotype_label]
  vst <- all_vst_full[, col_names, drop = FALSE]
  cv  <- (rowSds(vst) / rowMeans(vst)) * 100
  vst <- vst[which(cv >= 15), ]

  # --- TEST MODE: slice to first n genes ---
  if (!is.null(n_genes_test)) {
    n_genes_test <- min(n_genes_test, nrow(vst))
    vst <- vst[seq_len(n_genes_test), ]
    message("*** TEST MODE: using ", n_genes_test, " genes only ***")
  }

  genes   <- rownames(vst)
  n_genes <- nrow(vst)
  vst_t   <- t(vst)  # samples x genes

  out_file_xi <- paste0("/home/genomics/jorge/files/matrix_sugarcane_", group_name, "_xicor.tsv")
  out_file_pv <- paste0("/home/genomics/jorge/files/matrix_sugarcane_", group_name, "_xipvals.tsv")

  # 2. Initialize Files with Headers (one row, genes as columns)
  fwrite(as.data.table(t(c("GeneID", genes))), file = out_file_xi, sep = "\t", col.names = FALSE)
  fwrite(as.data.table(t(c("GeneID", genes))), file = out_file_pv, sep = "\t", col.names = FALSE)

  # 3. Chunked Processing
  main_indices <- seq_len(n_genes)
  block_size   <- chunk_size * n_cores
  blocks       <- split(main_indices, ceiling(seq_along(main_indices) / block_size))

  message("Total Genes: ", n_genes,
          " | Block size: ", block_size,
          " | Total Blocks: ", length(blocks))

  for (b in seq_along(blocks)) {
    curr_block <- blocks[[b]]
    message(sprintf("[%s] Block %d/%d — genes %d to %d",
                    group_name, b, length(blocks),
                    min(curr_block), max(curr_block)))

    results <- mclapply(curr_block, function(i) {
      tryCatch({
        x_vals <- vst_t[, i]
        both <- lapply(seq_len(n_genes), function(j) {
          obj <- XICOR::xicor(x_vals, vst_t[, j], pvalue = TRUE, ties = TRUE)
          c(obj$xi, obj$pval)
        })
        mat <- do.call(rbind, both)   # n_genes x 2 matrix
        list(xi_row = round(mat[, 1], 4),
             pv_row = signif(mat[, 2], 4))
      }, error = function(e) {
        structure(conditionMessage(e), class = "worker-error")
      })
    }, mc.cores = n_cores)

    # 4. Error checking
    failed <- vapply(results, inherits, logical(1L), what = "worker-error")
    if (any(failed)) {
      err_msgs <- unique(unlist(results[failed]))
      stop("Worker error(s) in block ", b, ":\n", paste(err_msgs, collapse = "\n"))
    }

    # 5. Write Xi values
    block_xi <- as.data.table(do.call(rbind, lapply(results, `[[`, "xi_row")))
    setnames(block_xi, genes)
    block_xi[, GeneID := genes[curr_block]]
    setcolorder(block_xi, "GeneID")
    fwrite(block_xi, file = out_file_xi, sep = "\t", append = TRUE, col.names = FALSE)

    # 6. Write P-values
    block_pv <- as.data.table(do.call(rbind, lapply(results, `[[`, "pv_row")))
    setnames(block_pv, genes)
    block_pv[, GeneID := genes[curr_block]]
    setcolorder(block_pv, "GeneID")
    fwrite(block_pv, file = out_file_pv, sep = "\t", append = TRUE, col.names = FALSE)

    progress <- (max(curr_block) / n_genes) * 100
    message(sprintf("[%s] Progress: %.1f%% (%d/%d genes)",
                    group_name, progress, max(curr_block), n_genes))
    rm(results, block_xi, block_pv)
    gc()
  }
  message(">>> Done: ", group_name)
}

# --- MAIN EXECUTION ---
genotype_groups <- list(
  responsive     = "R",
  non_responsive = "NR"
)

# --- TEST RUN (100 genes, one group) ---
# Uncomment to test before full run:
 process_group_linux("responsive", "R", n_genes_test = 1000)

# --- FULL RUN ---
#for (name in names(genotype_groups)) {
#  process_group_linux(name, genotype_groups[[name]])
#}

message("\nAll groups processed.")
