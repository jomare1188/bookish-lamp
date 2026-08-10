library(data.table)
library(DESeq2)

setDTthreads(100)

# ==============================================================================
# WHY THIS SCRIPT EXISTS
# ==============================================================================
# module_trait_cor.r correlates trait values against module *eigengenes*
# (PC1 of all genes in a module). That's a deliberate summary, but it can
# dilute a real nitrogen signal if it's carried by a subset of genes within
# an otherwise heterogeneous module.
#
# This script skips the eigengene step entirely and correlates each gene's
# own VST expression directly against trait, restricted to the gene set
# that came out of network_conservation_join.r (genes sitting on at least
# one conserved cross-species edge). The correlation itself uses the same
# vectorised analytic t-test as pearson_cor.r (fast even across tens of
# thousands of genes), and BH correction / rounding follows the same
# convention as module_trait_cor.r.
# ==============================================================================

# ==============================================================================
# CONFIGURATION
# ==============================================================================
NETWORK_NAME <- "purple"

# Restrict the correlation universe to genes on conserved cross-species
# edges (output of network_conservation_join.r). Set to NULL to correlate
# every gene in the VST matrix instead.
# GENE_FILTER_FILE <- "/dados04/jorge/comparative_saccharum/files/network_conservation/conserved_genes_purple_FULL.txt"
# "/dados04/jorge/comparative_saccharum/files/network_conservation/conserved_genes_sugarcane_FULL.txt"
# "/dados04/jorge/comparative_saccharum/files/network_conservation/conserved_genes_purple_FULL.txt"
PEARSON_THR <- 0.6     # |r| threshold for gene selection
PADJ_THR    <- 0.05    # FDR threshold

#TRAIT_ENCODING <- list(
#  genotype  = c("RB975375" = 1L, "RB937570" = 0L),   # responsive=1, non-responsive=0
#  treatment = c("High Nitrogen" = 1L, "Low Nitrogen" = 0L)
#)

# ── INPUT FILES -- SUGARCANE ─────────────────────────────────────────────────
#DDS_PATH      <- "/dados04/jorge/comparative_saccharum/run1/salmon/deseq2_qc/deseq2.dds.RData"
#METADATA_FILE <- "/dados04/jorge/comparative_saccharum/samplesheet.csv"
#OUT_DIR       <- "/dados04/jorge/comparative_saccharum/files/sugarcane/"

# ── INPUT FILES -- PURPLE (swap in, comment the sugarcane block above) ──────
NETWORK_NAME     <- "purple"
 GENE_FILTER_FILE <- "/dados04/jorge/comparative_saccharum/files/network_conservation/conserved_genes_purple_FULL.txt"
 DDS_PATH         <- "/dados04/jorge/comparative_saccharum/china/run2_onlyL/salmon/deseq2_qc/deseq2.dds.RData"
 METADATA_FILE    <- "/dados04/jorge/comparative_saccharum/china/samplesheet_china.csv"
 OUT_DIR          <- "/dados04/jorge/comparative_saccharum/files/purple/new/"
 TRAIT_ENCODING <- list(
   genotype  = c("51NG3" = 1L, "TAGZ" = 0L),
   treatment = c("0N" = 0L, "2N" = 2L, "6N" = 6L)
 )

# ==============================================================================
# 1. LOAD VST EXPRESSION MATRIX
# ==============================================================================
message("Loading DESeq2 object...")
load(DDS_PATH)   # loads 'dds'
all_vst    <- as.matrix(assay(dds, "vst"))
sample_ids <- colnames(all_vst)
message(sprintf("  VST matrix: %s genes x %s samples",
                format(nrow(all_vst), big.mark = ","), ncol(all_vst)))

# ==============================================================================
# 2. RESTRICT TO CONSERVED-EDGE GENES (if requested)
# ==============================================================================
if (!is.null(GENE_FILTER_FILE) && nzchar(GENE_FILTER_FILE)) {
  message("Restricting to conserved-edge gene set: ", basename(GENE_FILTER_FILE))
  rownames(all_vst) <- sub("\\.v[0-9]+\\.[0-9]+$", "", rownames(all_vst))   # strip .v2.1 etc.
  filter_genes <- readLines(GENE_FILTER_FILE)
  keep_genes   <- intersect(filter_genes, rownames(all_vst))
  message(sprintf("  %s / %s conserved genes found in VST matrix",
                  format(length(keep_genes), big.mark = ","),
                  format(length(filter_genes), big.mark = ",")))
  vst <- all_vst[keep_genes, , drop = FALSE]
} else {
  message("No gene filter supplied -- using all genes in VST matrix.")
  rownames(all_vst) <- sub("\\.v[0-9]+\\.[0-9]+$", "", rownames(all_vst))   # strip .v2.1 etc.
  vst <- all_vst
}
rm(all_vst); gc()
genes   <- rownames(vst)
n_genes <- length(genes)
message(sprintf("  Genes to correlate: %s", format(n_genes, big.mark = ",")))

# ==============================================================================
# 3. LOAD METADATA AND ENCODE TRAITS
# ==============================================================================
meta <- fread(METADATA_FILE, header = TRUE)
setnames(meta, tolower(names(meta)))

for (tr in names(TRAIT_ENCODING)) {
  enc <- TRAIT_ENCODING[[tr]]
  meta[, (paste0(tr, "_num")) := enc[get(tr)]]
}
trait_num_cols <- paste0(names(TRAIT_ENCODING), "_num")

# Match metadata rows to VST sample columns, preserving VST's sample order
meta_matched <- meta[match(sample_ids, sample)]
if (any(is.na(meta_matched$sample)))
  stop("Some VST samples were not found in metadata -- check sample naming.")

# ==============================================================================
# 4. VECTORISED PER-GENE CORRELATION  (same analytic t-test as pearson_cor.r)
# ==============================================================================
message("\nComputing gene-level correlations...")

vst_t <- t(vst)   # samples x genes -- reused for every trait

cor_results <- rbindlist(lapply(names(TRAIT_ENCODING), function(tr) {

  trait_vec <- meta_matched[[paste0(tr, "_num")]]
  ok        <- !is.na(trait_vec)
  n_ok      <- sum(ok)

  if (n_ok < length(trait_vec))
    warning(sprintf("%s: %d samples with unrecognised trait value excluded.",
                    tr, length(trait_vec) - n_ok))

  df <- n_ok - 2L
  r  <- as.numeric(cor(vst_t[ok, , drop = FALSE], trait_vec[ok]))

  t_stat <- r * sqrt(df / (1 - r^2 + 1e-15))
  pval   <- 2 * pt(-abs(t_stat), df = df)

  data.table(
    gene    = genes,
    trait   = tr,
    pearson = round(r, 4),
    pval    = pval
  )
}))

# BH correction within each trait independently
cor_results[, padj := signif(p.adjust(pval, method = "BH"), 4), by = trait]
setorder(cor_results, trait, padj)
setcolorder(cor_results, c("gene", "trait", "pearson", "pval", "padj"))

out_cor <- file.path(OUT_DIR, sprintf("gene_trait_correlations_%s.tsv", NETWORK_NAME))
fwrite(cor_results, file = out_cor, sep = "\t", quote = FALSE)
cat("Saved:", basename(out_cor), "\n")

for (tr in names(TRAIT_ENCODING)) {
  sig <- cor_results[trait == tr & padj <= PADJ_THR, .N]
  cat(sprintf("  %-12s  %d genes with padj <= %.2f\n", tr, sig, PADJ_THR))
}

# ==============================================================================
# 5. SELECT NITROGEN-RESPONSIVE GENES
# ==============================================================================
message("\nSelecting genes...")

selected <- cor_results[
  trait        == "treatment" &
  abs(pearson) >= PEARSON_THR &
  padj         <= PADJ_THR
][order(-abs(pearson))]

cat(sprintf("  Genes passing filters (|r|>=%.1f, padj<=%.2f): %d\n",
            PEARSON_THR, PADJ_THR, nrow(selected)))

if (nrow(selected) > 0L) print(selected, row.names = FALSE)

out_sel <- file.path(OUT_DIR, sprintf("selected_genes_treatment_%s.tsv", NETWORK_NAME))
fwrite(selected, file = out_sel, sep = "\t", quote = FALSE)
cat("Saved:", basename(out_sel), "\n")

message("\nDone.")
message("\nNote: these genes are already restricted to the conserved-edge gene")
message("set, so 'selected_genes_treatment_*.tsv' is, by construction, the")
message("intersection of 'nitrogen-responsive' and 'cross-species conserved'.")
message("To get edges (not just nodes) for visualisation, re-filter")
message("conserved_edges_*_FULL.tsv to rows where both gene1 and gene2 are in")
message("this selected gene list.")
