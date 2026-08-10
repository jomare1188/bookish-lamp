library(data.table)
library(parallel)

# ==============================================================================
# CONFIGURATION — adjust these parameters as needed
# ==============================================================================
PEARSON_THRESHOLD <- 0.6    # minimum ABSOLUTE Pearson r to keep an edge
                             # catches both positive (>0.6) and negative (<-0.6)
PVAL_RAW_MAX      <- 0.05   # maximum RAW p-value (pre-FDR) — pre-filter
PADJ_MAX          <- 0.05   # maximum FDR-adjusted p-value in final edge list
N_CORES           <- 100    # cores available
ROWS_PER_CORE     <- 50   # rows per worker per block; tune for RAM/speed

Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1)

# ==============================================================================
# INPUT FILES
# ==============================================================================
input_groups <- list(
  responsive = list(
    cor_file = "/home/genomics/jorge/files/purple/matrix_purple_responsive_pearson.tsv",
    pv_file  = "/home/genomics/jorge/files/purple/matrix_purple_responsive_pvalues.tsv",
    out_file = "/home/genomics/jorge/files/purple/edgelist_purple_responsive_pearson.tsv"
  ),
  non_responsive = list(
    cor_file = "/home/genomics/jorge/files/purple/matrix_purple_non_responsive_pearson.tsv",
    pv_file  = "/home/genomics/jorge/files/purple/matrix_purple_non_responsive_pvalues.tsv",
    out_file = "/home/genomics/jorge/files/purple/edgelist_purple_non_responsive_pearson.tsv"
  )
)

# ==============================================================================
# HELPER: detect file structure for the Pearson matrix format
#
# The Pearson script writes a proper single tab-separated header line:
#   Line 1     : GeneID \t gene1 \t gene2 \t ... \t geneN
#   Lines 2..  : data rows: gene_id \t val1 \t val2 \t ...
#
# Returns: list(col_genes, n_genes, skip_n = 1)
# ==============================================================================
detect_file_structure <- function(file) {
  if (!file.exists(file)) stop("File not found: ", file)

  # Read first line — this is the full tab-separated header
  header_line <- readLines(file, n = 1L, warn = FALSE)
  fields      <- strsplit(header_line, "\t", fixed = TRUE)[[1L]]

  # First field is "GeneID" label, rest are gene names
  if (fields[1L] != "GeneID") {
    stop("Expected 'GeneID' as first field in header of: ", file,
         "\n  Got: '", fields[1L], "' — check file format.")
  }

  col_genes <- fields[-1L]
  n_genes   <- length(col_genes)

  message("  Header: single tab-separated line | ",
          format(n_genes, big.mark = ","), " column genes detected")

  list(col_genes = col_genes, n_genes = n_genes, skip_n = 1L)
}

# ==============================================================================
# HELPER: process one block of raw text lines
#
# cor_lines / pv_lines : character vectors — one raw TSV line per element
# col_genes            : column gene names from header (in order)
# row_offset           : global 1-based index of first row in this block
# cor_thr              : absolute Pearson threshold
# pval_max             : raw p-value upper limit
#
# Upper triangle (j > global_i) removes self-loops and duplicate edges.
# strsplit is called once per line; only needed field positions extracted.
# ==============================================================================
process_block <- function(cor_lines, pv_lines, col_genes, row_offset,
                          cor_thr, pval_max) {
  n_genes <- length(col_genes)

  results <- mclapply(seq_along(cor_lines), function(local_i) {
    global_i <- row_offset + local_i - 1L

    # Upper triangle only
    j_indices <- (global_i + 1L):n_genes
    if (length(j_indices) == 0L) return(NULL)

    # Field 1 = gene name, fields 2..N+1 = values
    # Gene j → field position j+1
    cor_parts <- strsplit(cor_lines[local_i], "\t", fixed = TRUE)[[1L]]
    pv_parts  <- strsplit(pv_lines[local_i],  "\t", fixed = TRUE)[[1L]]

    col_pos   <- j_indices + 1L
    cor_vals  <- as.numeric(cor_parts[col_pos])
    pv_vals   <- as.numeric(pv_parts[col_pos])

    # Filter: absolute Pearson >= threshold AND raw pval <= max
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
# MAIN: process one group end-to-end
# ==============================================================================
process_group <- function(group_name, cor_file, pv_file, out_file) {
  message("\n", strrep("=", 60))
  message("Processing group: ", group_name)
  message(strrep("=", 60))

  for (f in c(cor_file, pv_file)) {
    if (!file.exists(f)) stop("File not found: ", f)
  }

  # ── Detect structure ──────────────────────────────────────────────
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

  message("\nGenes:                 ", format(n_genes,            big.mark = ","))
  message("Block size:            ", format(block_size,          big.mark = ","),
          " rows | Total blocks: ", n_blocks)
  message("Upper-triangle pairs:  ", format(choose(n_genes, 2), big.mark = ","))

  # ── Open persistent connections — skip header once, read sequentially ────
  # This avoids fread's large-skip segfault on wide matrices
  cor_con <- file(cor_file, open = "r")
  pv_con  <- file(pv_file,  open = "r")
  on.exit({ close(cor_con); close(pv_con) }, add = TRUE)

  # Skip exactly 1 header line in each file
  readLines(cor_con, n = cor_info$skip_n, warn = FALSE)
  readLines(pv_con,  n = pv_info$skip_n,  warn = FALSE)

  # ── Collect candidate edges across all blocks ─────────────────────
  all_edges   <- vector("list", n_blocks)
  block_count <- 0L
  row_start   <- 1L

  repeat {
    cor_lines <- readLines(cor_con, n = block_size, warn = FALSE)
    pv_lines  <- readLines(pv_con,  n = block_size, warn = FALSE)

    if (length(cor_lines) == 0L) break
    if (length(cor_lines) != length(pv_lines)) {
      stop("Line count mismatch between correlation and p-value files at block ",
           block_count + 1L)
    }

    block_count <- block_count + 1L

    block_edges <- process_block(
      cor_lines, pv_lines, col_genes,
      row_offset = row_start,
      cor_thr    = PEARSON_THRESHOLD,
      pval_max   = PVAL_RAW_MAX
    )

    all_edges[[block_count]] <- block_edges
    row_start <- row_start + length(cor_lines)

    edges_so_far <- sum(vapply(all_edges,
                               function(x) if (is.null(x)) 0L else nrow(x),
                               integer(1L)))
    message(sprintf(
      "[%s] Block %d/%d — %.1f%% complete (%s/%s genes) | candidate edges: %s",
      group_name,
      block_count, n_blocks,
      (min(row_start - 1L, n_genes) / n_genes) * 100,
      format(min(row_start - 1L, n_genes), big.mark = ","),
      format(n_genes,      big.mark = ","),
      format(edges_so_far, big.mark = ",")
    ))

    rm(cor_lines, pv_lines, block_edges); gc()
  }

  # ── Combine candidate edges ───────────────────────────────────────
  message("\nCombining candidate edges...")
  edges <- rbindlist(Filter(Negate(is.null), all_edges))
  rm(all_edges); gc()
  message("Total candidate edges before FDR: ", format(nrow(edges), big.mark = ","))

  if (nrow(edges) == 0L) {
    message("No edges passed pre-filters for group: ", group_name)
    return(invisible(NULL))
  }

  # ── FDR correction (BH) across ALL candidate edges ────────────────
  message("Applying FDR correction (BH method)...")
  edges[, padj := p.adjust(pval, method = "BH")]

  edges <- edges[padj <= PADJ_MAX]
  message("Edges after FDR filter (padj <= ", PADJ_MAX, "): ",
          format(nrow(edges), big.mark = ","))

  if (nrow(edges) == 0L) {
    message("No edges survived FDR correction for group: ", group_name)
    return(invisible(NULL))
  }

  # ── Sort by adjusted p-value and write ───────────────────────────
  setorder(edges, padj)
  fwrite(edges, file = out_file, sep = "\t", quote = FALSE)
  message("Saved: ", out_file)
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

message("\nAll groups processed.")
