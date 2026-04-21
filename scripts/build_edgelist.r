library(data.table)
library(parallel)

# ==============================================================================
# CONFIGURATION
# ==============================================================================
PEARSON_THRESHOLD <- 0.7
PVAL_RAW_MAX      <- 0.05
PADJ_MAX          <- 0.05
N_CORES           <- 100
ROWS_PER_CORE     <- 50

Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1)

# ==============================================================================
# INPUT FILES
# ==============================================================================
input_groups <- list(
  purple = list(
    cor_file = "/home/genomics/jorge/files/sugarcane/matrix_sugarcane_pearson.tsv",
    pv_file  = "/home/genomics/jorge/files/sugarcane/matrix_sugarcane_pvalues.tsv",
    out_file = "/home/genomics/jorge/files/sugarcane/edgelist_sugarcane_pearson.tsv"
  )
)

# ==============================================================================
# HELPER: detect file structure
# ==============================================================================
detect_file_structure <- function(file) {
  if (!file.exists(file)) stop("File not found: ", file)
  header_line <- readLines(file, n = 1L, warn = FALSE)
  fields      <- strsplit(header_line, "\t", fixed = TRUE)[[1L]]
  if (fields[1L] != "GeneID") {
    stop("Expected 'GeneID' as first field in: ", file,
         "\n  Got: '", fields[1L], "'")
  }
  col_genes <- fields[-1L]
  n_genes   <- length(col_genes)
  message("  ", format(n_genes, big.mark = ","), " column genes detected")
  list(col_genes = col_genes, n_genes = n_genes, skip_n = 1L)
}

# ==============================================================================
# HELPER: process one block of lines → data.table or NULL
# ==============================================================================
process_block <- function(cor_lines, pv_lines, col_genes, row_offset,
                          cor_thr, pval_max) {
  n_genes <- length(col_genes)

  results <- mclapply(seq_along(cor_lines), function(local_i) {
    global_i  <- row_offset + local_i - 1L
    j_indices <- (global_i + 1L):n_genes
    if (length(j_indices) == 0L) return(NULL)

    cor_parts <- strsplit(cor_lines[local_i], "\t", fixed = TRUE)[[1L]]
    pv_parts  <- strsplit(pv_lines[local_i],  "\t", fixed = TRUE)[[1L]]

    col_pos  <- j_indices + 1L
    cor_vals <- as.numeric(cor_parts[col_pos])
    pv_vals  <- as.numeric(pv_parts[col_pos])

    keep <- which(abs(cor_vals) >= cor_thr & pv_vals <= pval_max)
    if (length(keep) == 0L) return(NULL)

    data.table(
      gene1   = cor_parts[1L],
      gene2   = col_genes[j_indices[keep]],
      pearson = round(cor_vals[keep], 4),
      pval    = pv_vals[keep]
    )
  }, mc.cores = N_CORES)

  results <- Filter(Negate(is.null), results)
  if (length(results) == 0L) return(NULL)
  rbindlist(results)
}

# ==============================================================================
# MAIN: process one study
#
# Three-pass strategy to avoid the 2^31 rbindlist row limit:
#   Pass 1 — stream blocks → append candidate edges directly to temp TSV on disk
#   Pass 2 — read pval column only → compute BH padj → identify surviving indices
#   Pass 3 — re-read temp file in chunks → write only surviving rows to final file
# ==============================================================================
process_group <- function(group_name, cor_file, pv_file, out_file) {

  message("\n", strrep("=", 60))
  message("Processing study: ", group_name)
  message(strrep("=", 60))

  if (!file.exists(cor_file)) stop("File not found: ", cor_file)
  if (!file.exists(pv_file))  stop("File not found: ", pv_file)

  message("Detecting structure: ", basename(cor_file))
  cor_info <- detect_file_structure(cor_file)
  message("Detecting structure: ", basename(pv_file))
  pv_info  <- detect_file_structure(pv_file)

  if (cor_info$n_genes != pv_info$n_genes) {
    stop("Gene count mismatch: cor=", cor_info$n_genes, " pval=", pv_info$n_genes)
  }
  if (!identical(cor_info$col_genes, pv_info$col_genes)) {
    stop("Gene order differs between correlation and p-value files.")
  }

  col_genes  <- cor_info$col_genes
  n_genes    <- cor_info$n_genes
  block_size <- ROWS_PER_CORE * N_CORES
  n_blocks   <- ceiling(n_genes / block_size)

  message("\nGenes:                ", format(n_genes,           big.mark = ","))
  message("Block size:           ", format(block_size,         big.mark = ","),
          " rows | Total blocks: ", n_blocks)
  message("Upper-triangle pairs: ", format(choose(n_genes, 2), big.mark = ","))

  tmp_candidates <- paste0(out_file, ".tmp_candidates.tsv")

  # ── PASS 1: stream blocks → disk ──────────────────────────────────
  message("\n[Pass 1] Streaming candidate edges to disk...")
  message("  Temp file: ", tmp_candidates)

  # Write header to temp file
  writeLines("gene1\tgene2\tpearson\tpval", con = tmp_candidates)

  cor_con <- file(cor_file, open = "r")
  pv_con  <- file(pv_file,  open = "r")

  # Skip header lines
  readLines(cor_con, n = cor_info$skip_n, warn = FALSE)
  readLines(pv_con,  n = pv_info$skip_n,  warn = FALSE)

  total_candidates <- 0    # numeric (double) — avoids 32-bit overflow at ~2.1B rows
  block_count      <- 0L
  row_start        <- 1L

  repeat {
    cor_lines <- readLines(cor_con, n = block_size, warn = FALSE)
    pv_lines  <- readLines(pv_con,  n = block_size, warn = FALSE)

    if (length(cor_lines) == 0L) break

    if (length(cor_lines) != length(pv_lines)) {
      close(cor_con); close(pv_con)
      stop("Line count mismatch between cor and pval files at block ",
           block_count + 1L)
    }

    block_count <- block_count + 1L

    block_edges <- process_block(
      cor_lines, pv_lines, col_genes,
      row_offset = row_start,
      cor_thr    = PEARSON_THRESHOLD,
      pval_max   = PVAL_RAW_MAX
    )

    if (!is.null(block_edges) && nrow(block_edges) > 0L) {
      fwrite(block_edges, file = tmp_candidates,
             sep = "\t", quote = FALSE, append = TRUE, col.names = FALSE)
      total_candidates <- total_candidates + nrow(block_edges)
    }

    row_start <- row_start + length(cor_lines)

    message(sprintf(
      "  [%s] Block %d/%d — %.1f%% (%s/%s genes) | candidates on disk: %s",
      group_name,
      block_count, n_blocks,
      (min(row_start - 1L, n_genes) / n_genes) * 100,
      format(min(row_start - 1L, n_genes), big.mark = ","),
      format(n_genes,          big.mark = ","),
      format(total_candidates, big.mark = ",")
    ))

    rm(cor_lines, pv_lines, block_edges)
    gc()
  }

  close(cor_con)
  close(pv_con)

  message("Total candidate edges before FDR: ",
          format(total_candidates, big.mark = ","))

  if (total_candidates == 0) {
    message("No edges passed pre-filters for study: ", group_name)
    if (file.exists(tmp_candidates)) file.remove(tmp_candidates)
    return(invisible(NULL))
  }

  # ── PASS 2: read pval column only → BH correction ──────────────────
  message("\n[Pass 2] Reading p-values for BH correction...")
  message("  Loading pval column from temp file...")

  pvals <- fread(tmp_candidates, select = "pval", header = TRUE)$pval

  message("  Computing BH adjusted p-values over ",
          format(length(pvals), big.mark = ","), " candidates...")

  padj     <- p.adjust(pvals, method = "BH")
  keep_idx <- which(padj <= PADJ_MAX)
  padj_kept <- padj[keep_idx]

  rm(pvals, padj)
  gc()

  message("  Edges surviving FDR (padj <= ", PADJ_MAX, "): ",
          format(length(keep_idx), big.mark = ","))

  if (length(keep_idx) == 0L) {
    message("No edges survived FDR correction for study: ", group_name)
    if (file.exists(tmp_candidates)) file.remove(tmp_candidates)
    return(invisible(NULL))
  }

  # ── PASS 3: extract surviving rows → final output ──────────────────
  # Read the full candidate file, subset by index, attach padj, sort, write.
  # If this is still too large for RAM, see chunked alternative below.
  message("\n[Pass 3] Writing final edge list...")
  message("  Reading surviving rows from temp file...")

  candidates  <- fread(tmp_candidates, header = TRUE)
  final_edges <- candidates[keep_idx]
  final_edges[, padj := padj_kept]

  rm(candidates, padj_kept, keep_idx)
  gc()

  setorder(final_edges, padj)
  fwrite(final_edges, file = out_file, sep = "\t", quote = FALSE)

  message("  Saved: ", out_file)
  message("  Final edges: ", format(nrow(final_edges), big.mark = ","))

  rm(final_edges)
  gc()

  if (file.exists(tmp_candidates)) file.remove(tmp_candidates)
  message("  Temp file removed.")
}

# ==============================================================================
# ENTRY POINT
# ==============================================================================
for (grp in names(input_groups)) {
  g <- input_groups[[grp]]
  process_group(
    group_name = grp,
    cor_file   = g$cor_file,
    pv_file    = g$pv_file,
    out_file   = g$out_file
  )
}

message("\nAll studies processed.")
