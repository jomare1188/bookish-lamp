library(DESeq2)
library(Matrix)
library(data.table)
library(matrixStats)
library(parallel)

setwd("/dados02/jorge/comparative_saccharum/china/run1/salmon/deseq2_qc")

# BLAS threading off — let mclapply own all cores
Sys.setenv(OMP_NUM_THREADS      = 1,
           OPENBLAS_NUM_THREADS = 1,
           MKL_NUM_THREADS      = 1,
           BLAS_NUM_THREADS     = 1)

# ==============================================================================
# CONFIGURATION
# ==============================================================================
MIN_CV     <- 15    # CV threshold (%)
N_CORES    <- 100   # parallel workers (mclapply forks — no socket limit)
CHUNK_ROWS <- 500   # genes per chunk per worker
                    # RAM per batch = N_CORES × CHUNK_ROWS × n_genes × 8 × 2
                    # e.g. 256 workers × 500 rows × 161k genes × 16 B ≈ 330 GB
                    # Lower CHUNK_ROWS if RAM is tight; raise if headroom allows

# ==============================================================================
# INPUT / OUTPUT
# ==============================================================================
genotype_groups <- list(
  responsive     = "51N03",
  non_responsive = "TA0Z"
)
out_dir <- "/home/genomics/jorge/files/purple/"

# ==============================================================================
# LOAD DATA ONCE  (shared across all forked workers at zero copy cost)
# ==============================================================================
message("Loading data...")
load("deseq2.dds.RData")
all_vst <- Matrix::as.matrix(assay(dds, "vst"))
all_raw <- counts(dds, normalized = FALSE)
message("VST matrix: ", paste(dim(all_vst), collapse = " x "))

# ==============================================================================
# WORKER: compute Pearson r + p-values for one chunk, write to temp files
#
# Each worker:
#   1. Slices its assigned rows from vst (inherited via fork)
#   2. Computes cor(chunk, full_matrix)
#   3. Computes p-values in-place and frees t_stat immediately
#   4. Writes cor and pval blocks to numbered temp files
#   5. Returns the temp file paths for the main process to concatenate
#
# Arguments passed explicitly (everything else inherited from parent env):
#   chunk_idx   — chunk number (for ordering and temp file naming)
#   row_start   — first gene index in this chunk (1-based)
#   row_end     — last  gene index in this chunk (1-based)
#   vst         — full genes × samples VST matrix (inherited via fork)
#   vst_t       — transposed: samples × genes     (inherited via fork)
#   genes       — gene name vector                (inherited via fork)
#   n_samp      — number of samples               (inherited via fork)
#   tmp_dir     — directory for temp files        (inherited via fork)
# ==============================================================================
worker_chunk <- function(chunk_idx, row_start, row_end) {
  tryCatch({
    chunk_genes <- genes[row_start:row_end]
    vst_chunk   <- vst[row_start:row_end, , drop = FALSE]
    df          <- n_samp - 2L

    # Pearson r: chunk_genes × all_genes
    cor_block  <- cor(t(vst_chunk), vst_t)

    # P-values — computed and freed before writing to minimise peak RAM
    t_stat     <- cor_block * sqrt(df / (1 - cor_block^2 + 1e-15))
    pval_block <- 2 * pt(-abs(t_stat), df = df)
    rm(t_stat, vst_chunk); gc()

    # Write to numbered temp files so main process can cat in correct order
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
# MAIN LOOP
# ==============================================================================
for (group_name in names(genotype_groups)) {
  genotype_label <- genotype_groups[[group_name]]

  message("\n", strrep("=", 60))
  message("Processing group: ", group_name, " (", genotype_label, ")")
  message(strrep("=", 60))

  # ── 1. Sample selection ───────────────────────────────────────────

col_idx <- which(dds$Group1 == "L" & dds$Group2 == genotype_label)
#  col_idx <- which(dds$Group1 == genotype_label) # for sugarcane only 
  message("Samples: ", length(col_idx))

  # ── 2. CV filtering ───────────────────────────────────────────────
  raw_sub <- all_raw[, col_idx, drop = FALSE]
  m_raw   <- rowMeans(raw_sub)
  sd_raw  <- rowSds(raw_sub)
  cv_raw  <- (sd_raw / (m_raw + 1e-6)) * 100
  keep    <- which(cv_raw >= MIN_CV)
  rm(raw_sub, m_raw, sd_raw, cv_raw); gc()

  # These are inherited by all forked workers at zero copy cost
  vst     <- all_vst[keep, col_idx, drop = FALSE]
  genes   <- rownames(vst)
  n_genes <- length(genes)
  n_samp  <- ncol(vst)
  vst_t   <- t(vst)   # samples × genes — built once, shared by all workers

  n_chunks <- ceiling(n_genes / CHUNK_ROWS)

  message("Genes after CV filter:   ", format(n_genes,   big.mark = ","))
  message("Samples: ", n_samp, " | df: ", n_samp - 2L)
  message("Chunk size: ", CHUNK_ROWS, " | Chunks: ", n_chunks,
          " | Workers: ", N_CORES)
  message(sprintf("RAM per parallel batch:  ~%.0f GB",
                  (N_CORES * CHUNK_ROWS * n_genes * 8 * 2) / 1e9))
  message(sprintf("Full matrix (reference): ~%.0f GB",
                  (n_genes^2 * 8 * 2) / 1e9))

  # ── 3. Temp directory for chunk files ─────────────────────────────
  tmp_dir <- file.path(out_dir, paste0(".tmp_", group_name))
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

  # ── 4. Build chunk index table ────────────────────────────────────
  chunk_starts <- seq(1L, n_genes, by = CHUNK_ROWS)
  chunk_ends   <- pmin(chunk_starts + CHUNK_ROWS - 1L, n_genes)
  chunk_ids    <- seq_along(chunk_starts)

  # ── 5. Process all chunks in parallel ────────────────────────────
  # mclapply forks — workers inherit vst, vst_t, genes, n_samp, tmp_dir
  # from the parent process (copy-on-write, no explicit export needed)
  message("Launching parallel computation...")

  results <- mclapply(
    chunk_ids,
    function(ci) worker_chunk(ci, chunk_starts[ci], chunk_ends[ci]),
    mc.cores             = N_CORES,
    mc.preschedule       = FALSE,   # dynamic scheduling — better load balance
    mc.allow.recursive   = FALSE
  )

  # ── 6. Check for worker errors ────────────────────────────────────
  failed <- which(!vapply(results, `[[`, logical(1L), "ok"))
  if (length(failed) > 0L) {
    msgs <- vapply(results[failed], `[[`, character(1L), "msg")
    stop("Errors in chunks ", paste(failed, collapse = ", "), ":\n",
         paste(msgs, collapse = "\n"))
  }
  message("All ", n_chunks, " chunks computed successfully.")

  # ── 7. Concatenate temp files in order → final output files ───────
  out_cor  <- file.path(out_dir,
               paste0("matrix_purple_", group_name, "_pearson.tsv"))
  out_pval <- file.path(out_dir,
               paste0("matrix_purple_", group_name, "_pvalues.tsv"))

  # Write header lines
  header_line <- paste(c("GeneID", genes), collapse = "\t")
  writeLines(header_line, con = out_cor)
  writeLines(header_line, con = out_pval)

  message("Concatenating ", n_chunks, " chunk files in order...")

  # Sort results by chunk index to guarantee row order
  results <- results[order(vapply(results, `[[`, integer(1L), "chunk"))]

  for (res in results) {
    # Append cor chunk
    con <- file(out_cor, open = "a")
    file.append(out_cor, res$cor)
    close(con)
    # Append pval chunk
    con <- file(out_pval, open = "a")
    file.append(out_pval, res$pval)
    close(con)
    # Remove temp files as we go to free disk space
    file.remove(res$cor, res$pval)
  }

  unlink(tmp_dir, recursive = TRUE)

  message("Saved: ", out_cor)
  message("Saved: ", out_pval)

  rm(vst, vst_t, keep, results); gc()
}

message("\nAll groups processed.")
