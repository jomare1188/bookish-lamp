library(DESeq2)
library(Matrix)
library(data.table)
library(matrixStats)
library(parallel)


# CONFIGURE THIS
# For sugarcane

in_dir="/dados02/jorge/comparative_saccharum/run1/salmon/deseq2_qc"
out_dir="/home/genomics/jorge/files/sugarcane/"
label="sugarcane"

# For the other study

#in_dir="/dados02/jorge/comparative_saccharum/china/run2_onlyL/salmon/deseq2_qc"
#out_dir="/home/genomics/jorge/files/purple/"
#label="purple"





setwd(in_dir)

Sys.setenv(OMP_NUM_THREADS      = 1,
           OPENBLAS_NUM_THREADS = 1,
           MKL_NUM_THREADS      = 1,
           BLAS_NUM_THREADS     = 1)

# ==============================================================================
# CONFIGURATION
# ==============================================================================
MIN_CV     <- 15
N_CORES    <- 100
CHUNK_ROWS <- 500

# ==============================================================================
# INPUT / OUTPUT
# ==============================================================================
# One entry per study — combines all genotypes within each study
study_groups <- list(
#  purple = list(
#    label      = "purple",                          # Ta Quang Kiet et al. 2025
#    col_filter = function(dds) which(dds$Group1 == "L"),  # all leaf samples, both genotypes
#    out_dir    = "/home/genomics/jorge/files/purple/"
 
  # Add the Muñoz-Perez study here analogously when ready:
   sugarcane = list(
     label      = label,
     col_filter = function(dds) seq_len(ncol(dds)),   # all samples
     out_dir    = out_dir
   )
)

# ==============================================================================
# LOAD DATA ONCE
# ==============================================================================
message("Loading data...")
load("deseq2.dds.RData")
all_vst <- Matrix::as.matrix(assay(dds, "vst"))
all_raw <- counts(dds, normalized = FALSE)
message("VST matrix: ", paste(dim(all_vst), collapse = " x "))

# ==============================================================================
# WORKER
# ==============================================================================
worker_chunk <- function(chunk_idx, row_start, row_end) {
  tryCatch({
    chunk_genes <- genes[row_start:row_end]
    vst_chunk   <- vst[row_start:row_end, , drop = FALSE]
    df          <- n_samp - 2L

    cor_block  <- cor(t(vst_chunk), vst_t)

    t_stat     <- cor_block * sqrt(df / (1 - cor_block^2 + 1e-15))
    pval_block <- 2 * pt(-abs(t_stat), df = df)
    rm(t_stat, vst_chunk); gc()

    tmp_cor  <- file.path(tmp_dir, sprintf("chunk_%05d_cor.tsv",  chunk_idx))
    tmp_pval <- file.path(tmp_dir, sprintf("chunk_%05d_pval.tsv", chunk_idx))

    cor_dt  <- data.table(GeneID = chunk_genes, as.data.table(round(cor_block,  4)))
    pval_dt <- data.table(GeneID = chunk_genes, as.data.table(signif(pval_block, 4)))
    setnames(cor_dt,  c("GeneID", genes))
    setnames(pval_dt, c("GeneID", genes))

    fwrite(cor_dt,  file = tmp_cor,  sep = "\t", quote = FALSE, col.names = FALSE)
    fwrite(pval_dt, file = tmp_pval, sep = "\t", quote = FALSE, col.names = FALSE)

    rm(cor_block, pval_block, cor_dt, pval_dt); gc()

    list(ok = TRUE, cor = tmp_cor, pval = tmp_pval, chunk = chunk_idx)

  }, error = function(e) {
    list(ok = FALSE, chunk = chunk_idx, msg = conditionMessage(e))
  })
}

# ==============================================================================
# MAIN LOOP — one iteration per study
# ==============================================================================
for (study_name in names(study_groups)) {
  study  <- study_groups[[study_name]]
  out_dir <- study$out_dir

  message("\n", strrep("=", 60))
  message("Processing study: ", study_name)
  message(strrep("=", 60))

  # ── 1. Sample selection (all genotypes in this study) ─────────────
  col_idx <- study$col_filter(dds)
  message("Samples selected: ", length(col_idx))

  # ── 2. CV filtering ───────────────────────────────────────────────
  raw_sub <- all_raw[, col_idx, drop = FALSE]
  m_raw   <- rowMeans(raw_sub)
  sd_raw  <- rowSds(raw_sub)
  cv_raw  <- (sd_raw / (m_raw + 1e-6)) * 100
  keep    <- which(cv_raw >= MIN_CV)
  rm(raw_sub, m_raw, sd_raw, cv_raw); gc()

  vst     <- all_vst[keep, col_idx, drop = FALSE]
  genes   <- rownames(vst)
  n_genes <- length(genes)
  n_samp  <- ncol(vst)
  vst_t   <- t(vst)

  n_chunks <- ceiling(n_genes / CHUNK_ROWS)

  message("Genes after CV filter:   ", format(n_genes, big.mark = ","))
  message("Samples: ", n_samp, " | df: ", n_samp - 2L)
  message("Chunks: ", n_chunks, " | Workers: ", N_CORES)
  message(sprintf("RAM per parallel batch:  ~%.0f GB",
                  (N_CORES * CHUNK_ROWS * n_genes * 8 * 2) / 1e9))

  # ── 3. Temp directory ─────────────────────────────────────────────
  tmp_dir <- file.path(out_dir, paste0(".tmp_", study_name))
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

  # ── 4. Chunk index ────────────────────────────────────────────────
  chunk_starts <- seq(1L, n_genes, by = CHUNK_ROWS)
  chunk_ends   <- pmin(chunk_starts + CHUNK_ROWS - 1L, n_genes)
  chunk_ids    <- seq_along(chunk_starts)

  # ── 5. Parallel computation ───────────────────────────────────────
  message("Launching parallel computation...")
  results <- mclapply(
    chunk_ids,
    function(ci) worker_chunk(ci, chunk_starts[ci], chunk_ends[ci]),
    mc.cores           = N_CORES,
    mc.preschedule     = FALSE,
    mc.allow.recursive = FALSE
  )

  # ── 6. Check errors ───────────────────────────────────────────────
  failed <- which(!vapply(results, `[[`, logical(1L), "ok"))
  if (length(failed) > 0L) {
    msgs <- vapply(results[failed], `[[`, character(1L), "msg")
    stop("Errors in chunks ", paste(failed, collapse = ", "), ":\n",
         paste(msgs, collapse = "\n"))
  }
  message("All ", n_chunks, " chunks computed successfully.")

  # ── 7. Concatenate in order ───────────────────────────────────────
  # Output: one matrix per study (not per genotype)
  out_cor  <- file.path(out_dir, paste0("matrix_", study_name, "_pearson.tsv"))
  out_pval <- file.path(out_dir, paste0("matrix_", study_name, "_pvalues.tsv"))

  header_line <- paste(c("GeneID", genes), collapse = "\t")
  writeLines(header_line, con = out_cor)
  writeLines(header_line, con = out_pval)

  results <- results[order(vapply(results, `[[`, integer(1L), "chunk"))]

  for (res in results) {
    file.append(out_cor,  res$cor)
    file.append(out_pval, res$pval)
    file.remove(res$cor, res$pval)
  }

  unlink(tmp_dir, recursive = TRUE)

  message("Saved: ", out_cor)
  message("Saved: ", out_pval)

  rm(vst, vst_t, keep, results); gc()
}

message("\nAll studies processed.")
