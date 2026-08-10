# ============================================================================
# 00_export_vst.r — export the VST matrix that the MI engine will consume
#
# WHY A SEPARATE STEP: the GPU engine is Python, and the expression matrix
# lives inside a DESeq2 `dds` R object. Rather than round-tripping through a
# 200 GB TSV, this writes a flat little binary triple:
#
#     <prefix>.f32        row-major float32, n_genes x n_samples
#     <prefix>.genes.txt  one gene id per line, in row order
#     <prefix>.meta.json  dimensions + sample names + provenance
#
# For sugarcane that is 170,790 x 48 x 4 B = 33 MB. numpy reads it in one
# fromfile() call with no parsing.
#
# CRITICAL — COMPARABILITY: the gene filter here is character-for-character the
# same as scripts/pearson_cor.r (coefficient of variation on RAW counts,
# CV >= MIN_CV, MIN_CV = 15). If the MI network were built on a different gene
# set, none of the conservation / module comparisons against the Pearson
# networks would be valid. The script therefore VERIFIES its gene list against
# the header of the existing Pearson matrix and aborts on any mismatch.
#
# RUN (see run_all.sh; both studies):
#   MI_DDS=... MI_PREFIX=... MI_PEARSON=... \
#     /home/genomics/miniconda3/envs/cor_env/bin/Rscript 00_export_vst.r
# ============================================================================

suppressMessages({
  library(DESeq2)
  library(matrixStats)
})

dds_file    <- Sys.getenv("MI_DDS")
out_prefix  <- Sys.getenv("MI_PREFIX")
pearson_ref <- Sys.getenv("MI_PEARSON")           # optional; "" disables the check
label       <- Sys.getenv("MI_LABEL", "unnamed")
min_cv      <- as.numeric(Sys.getenv("MIN_CV", "15"))

# Optional library subset, as "colDataColumn=value" (e.g. "Group1=L").
#
# This is not cosmetic. The purple Pearson network was built from china/run1 --
# 36 libraries, leaf AND root -- narrowed to the 18 leaf libraries with
# Group1 == "L", and the CV filter was then applied to that subset. Running the
# same filter on china/run2_onlyL instead gives a gene set that differs by 68
# genes, because the two quantifications do not agree on which genes are flat.
# Since MI is meant to be compared edge-for-edge against those Pearson networks,
# the library selection has to be reproduced exactly, not approximated.
col_filter  <- Sys.getenv("MI_COL_FILTER")

stopifnot(nzchar(dds_file), nzchar(out_prefix))
if (!file.exists(dds_file)) stop("dds not found: ", dds_file)
dir.create(dirname(out_prefix), showWarnings = FALSE, recursive = TRUE)

message("=== ", label, " ===")
message("Loading ", dds_file)
load(dds_file)                                     # provides `dds`

vst_all <- as.matrix(assay(dds, "vst"))
raw_all <- counts(dds, normalized = FALSE)
message("VST matrix: ", paste(dim(vst_all), collapse = " x "))

# --- library subset ---------------------------------------------------------
if (nzchar(col_filter)) {
  parts <- strsplit(col_filter, "=", fixed = TRUE)[[1L]]
  if (length(parts) != 2L) stop("MI_COL_FILTER must look like 'Group1=L'")
  fld <- parts[1L]; val <- parts[2L]
  if (!fld %in% names(colData(dds)))
    stop("colData has no column '", fld, "'; available: ",
         paste(names(colData(dds)), collapse = ", "))
  col_idx <- which(as.character(colData(dds)[[fld]]) == val)
  if (length(col_idx) == 0L)
    stop("MI_COL_FILTER '", col_filter, "' selected 0 libraries; ",
         fld, " values are: ",
         paste(unique(as.character(colData(dds)[[fld]])), collapse = ", "))
  message("Library subset ", col_filter, ": ", length(col_idx), " of ", ncol(vst_all))
  vst_all <- vst_all[, col_idx, drop = FALSE]
  raw_all <- raw_all[, col_idx, drop = FALSE]
}

# --- gene filter: identical to pearson_cor.r --------------------------------
m_raw  <- rowMeans(raw_all)
sd_raw <- rowSds(raw_all)
cv_raw <- (sd_raw / (m_raw + 1e-6)) * 100
keep   <- which(cv_raw >= min_cv)

vst   <- vst_all[keep, , drop = FALSE]
genes <- rownames(vst)
message(sprintf("Genes after CV >= %g filter: %s of %s",
                min_cv, format(length(genes), big.mark = ","),
                format(nrow(vst_all), big.mark = ",")))
message("Samples: ", ncol(vst))

# --- verify against the Pearson run -----------------------------------------
# Reads only the first line of the (200 GB) matrix.
if (nzchar(pearson_ref) && file.exists(pearson_ref)) {
  hdr <- strsplit(readLines(pearson_ref, n = 1L, warn = FALSE), "\t", fixed = TRUE)[[1L]]
  ref <- hdr[-1L]
  if (!identical(ref, genes)) {
    message("  Pearson matrix genes: ", format(length(ref),   big.mark = ","))
    message("  this export genes:    ", format(length(genes), big.mark = ","))
    message("  in Pearson only: ", length(setdiff(ref, genes)),
            " | in export only: ", length(setdiff(genes, ref)))
    stop("Gene set does not match the Pearson matrix — the MI network would ",
         "not be comparable. Fix the filter before continuing.")
  }
  message("Gene set matches ", basename(pearson_ref), " exactly.")
} else {
  message("NOTE: no Pearson matrix given — skipping the comparability check.")
}

# --- write -------------------------------------------------------------------
# t() because R is column-major and we want row-major (gene-major) on disk.
con <- file(paste0(out_prefix, ".f32"), "wb")
writeBin(as.vector(t(vst)), con, size = 4L)
close(con)

writeLines(genes, paste0(out_prefix, ".genes.txt"))

json_esc <- function(x) paste0('"', gsub('"', '\\\\"', x), '"')
writeLines(paste0(
  '{\n',
  '  "label": ',      json_esc(label),                                  ',\n',
  '  "n_genes": ',    length(genes),                                    ',\n',
  '  "n_samples": ',  ncol(vst),                                        ',\n',
  '  "min_cv": ',     min_cv,                                           ',\n',
  '  "source_dds": ', json_esc(dds_file),                               ',\n',
  '  "created": ',    json_esc(format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ',\n',
  '  "samples": [',   paste(json_esc(colnames(vst)), collapse = ", "),  ']\n',
  '}'), paste0(out_prefix, ".meta.json"))

message("Wrote ", out_prefix, ".{f32,genes.txt,meta.json}  (",
        round(file.size(paste0(out_prefix, ".f32")) / 1e6, 1), " MB)\n")
