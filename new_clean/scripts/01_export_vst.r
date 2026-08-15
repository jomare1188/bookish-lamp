# =============================================================================
# 01_export_vst.r — export the VST matrix that both network layers consume
#
# WHY A SEPARATE STEP: the GPU engine is Python and the expression matrix lives
# inside a DESeq2 `dds` R object. Rather than round-tripping through a 200 GB
# TSV, this writes a flat binary triple:
#
#     <prefix>.f32        row-major float32, n_genes x n_samples
#     <prefix>.genes.txt  one gene id per line, in row order
#     <prefix>.meta.json  dimensions + sample names + provenance
#
# For sugarcane that is 170,790 x 48 x 4 B = 33 MB, read by numpy in a single
# fromfile() call with no parsing.
#
# WHY THIS IS THE ONLY GENE FILTER IN THE PIPELINE: both the Pearson layer and
# the MI layer read this one file. In the old tree the two layers were filtered
# independently and 00_export_vst.r had to VERIFY its gene list against the
# header of the 200 GB Pearson matrix to prove they agreed. That check is gone
# because the situation it guarded against can no longer arise -- there is one
# matrix, so the gene sets are identical by construction.
#
# GENE IDS ARE STRIPPED HERE, ONCE. See the note above strip_version() below.
#
# RUN (via run.sh; env-driven):
#   CLEAN_DDS=... CLEAN_PREFIX=... CLEAN_LABEL=... MIN_CV=15 \
#     /home/genomics/miniconda3/envs/cor_env/bin/Rscript 01_export_vst.r
# =============================================================================

suppressMessages({
  library(DESeq2)
  library(matrixStats)
})

dds_file   <- Sys.getenv("CLEAN_DDS")
out_prefix <- Sys.getenv("CLEAN_PREFIX")
label      <- Sys.getenv("CLEAN_LABEL", "unnamed")
min_cv     <- as.numeric(Sys.getenv("MIN_CV", "15"))
strip_ver  <- Sys.getenv("CLEAN_STRIP_VERSION", "0") == "1"

# Optional library subset, as "colDataColumn=value" (e.g. "Group1=L"). Both
# studies currently use the whole object -- purple's run2_onlyL is already
# exactly the 18 leaf libraries -- but the hook stays because narrowing a dds is
# a normal thing to want and doing it by hand is how the old tree drifted.
col_filter <- Sys.getenv("CLEAN_COL_FILTER")

stopifnot(nzchar(dds_file), nzchar(out_prefix))
if (!file.exists(dds_file)) stop("dds not found: ", dds_file)
dir.create(dirname(out_prefix), showWarnings = FALSE, recursive = TRUE)

message("=== ", label, " ===")
message("Loading ", dds_file)
load(dds_file)                                     # provides `dds`

vst_all <- as.matrix(assay(dds, "vst"))
raw_all <- counts(dds, normalized = FALSE)
message("VST matrix: ", paste(dim(vst_all), collapse = " x "))

# --- library subset ----------------------------------------------------------
if (nzchar(col_filter)) {
  parts <- strsplit(col_filter, "=", fixed = TRUE)[[1L]]
  if (length(parts) != 2L) stop("CLEAN_COL_FILTER must look like 'Group1=L'")
  fld <- parts[1L]; val <- parts[2L]
  if (!fld %in% names(colData(dds)))
    stop("colData has no column '", fld, "'; available: ",
         paste(names(colData(dds)), collapse = ", "))
  col_idx <- which(as.character(colData(dds)[[fld]]) == val)
  if (length(col_idx) == 0L)
    stop("CLEAN_COL_FILTER '", col_filter, "' selected 0 libraries; ",
         fld, " values are: ",
         paste(unique(as.character(colData(dds)[[fld]])), collapse = ", "))
  message("Library subset ", col_filter, ": ", length(col_idx), " of ", ncol(vst_all))
  vst_all <- vst_all[, col_idx, drop = FALSE]
  raw_all <- raw_all[, col_idx, drop = FALSE]
}

# --- gene filter: coefficient of variation on RAW counts ---------------------
# Character-for-character what the original pearson_cor.r did. Applied AFTER any
# library subset, so the CV reflects the libraries actually correlated.
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

# --- gene ids: strip the version suffix, ONCE, here --------------------------
# The R570 annotation carries ids like `Sh_236E10_t000010.v2.1` while the
# proteomes, orthogroups and TF tables use the bare form. In the old tree the
# strip happened in different places at different times -- node_metrics from
# Jun 2 had the suffix, the edge list from Jun 25 did not -- and every readout
# carried a defensive sub() to cope. Those merges use all.x = TRUE, so a
# mismatch produced silent NAs rather than an error.
#
# Doing it here means every file in new_clean/ carries the bare id and the
# defensive strips downstream become harmless no-ops. The collision assert is
# the part that matters: if two distinct versioned ids ever collapse to the same
# bare id, silently keeping one would corrupt the network.
if (strip_ver) {
  bare <- sub("\\.v[0-9]+(\\.[0-9]+)*$", "", genes)
  n_changed <- sum(bare != genes)
  if (anyDuplicated(bare)) {
    dup <- unique(bare[duplicated(bare)])
    stop("stripping the version suffix collapses ", length(dup),
         " distinct gene ids, e.g. ",
         paste(head(genes[bare %in% dup], 6L), collapse = ", "),
         ". Refusing to continue -- fix the id convention first.")
  }
  message("Stripped version suffix from ", format(n_changed, big.mark = ","),
          " of ", format(length(genes), big.mark = ","), " gene ids")
  genes <- bare
} else {
  message("Gene ids used as-is (no version suffix to strip)")
}

# --- write -------------------------------------------------------------------
# t() because R is column-major and the engine wants row-major (gene-major).
con <- file(paste0(out_prefix, ".f32"), "wb")
writeBin(as.vector(t(vst)), con, size = 4L)
close(con)

writeLines(genes, paste0(out_prefix, ".genes.txt"))

json_esc <- function(x) paste0('"', gsub('"', '\\\\"', x), '"')
writeLines(paste0(
  '{\n',
  '  "label": ',        json_esc(label),                                  ',\n',
  '  "n_genes": ',      length(genes),                                    ',\n',
  '  "n_samples": ',    ncol(vst),                                        ',\n',
  '  "min_cv": ',       min_cv,                                           ',\n',
  '  "col_filter": ',   json_esc(col_filter),                             ',\n',
  '  "stripped_version": ', if (strip_ver) "true" else "false",           ',\n',
  '  "source_dds": ',   json_esc(dds_file),                               ',\n',
  '  "created": ',      json_esc(format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ',\n',
  '  "samples": [',     paste(json_esc(colnames(vst)), collapse = ", "),  ']\n',
  '}'), paste0(out_prefix, ".meta.json"))

message("Wrote ", out_prefix, ".{f32,genes.txt,meta.json}  (",
        round(file.size(paste0(out_prefix, ".f32")) / 1e6, 1), " MB)\n")
