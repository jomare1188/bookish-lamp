library(data.table)

# ==============================================================================
# CONFIGURATION
# ==============================================================================
NETWORK_NAME <- "purple"
PEARSON_THR  <- 0.6       # |r| threshold for module selection
PADJ_THR     <- 0.05      # FDR threshold
MIN_GENES    <- 2       # minimum module size for selection

# Trait encoding — numeric values assigned to each factor level
# (Pearson on binary = point-biserial correlation, mathematically equivalent)




#TRAIT_ENCODING <- list(
#  genotype  = c("RB975375" = 1L, "RB937570" = 0L),   # responsive=1, non-responsive=0
#  treatment = c("High Nitrogen" = 1L, "Low Nitrogen" = 0L)
#)


TRAIT_ENCODING <- list(
  genotype  = c("51NG3" = 1L, "TAGZ" = 0L),   # responsive=1, non-responsive=0
  treatment = c("0N" = 0L, "2N" = 2L, "6N" = 6L )
)

# ==============================================================================
# INPUT FILES PURPLE
# ==============================================================================
EIGENGENE_FILE  <- "/home/genomics/jorge/files/purple/new/eigengenes_purple.tsv"
METADATA_FILE   <- "/dados02/jorge/comparative_saccharum/china/samplesheet_china.csv"
SUMMARY_FILE    <- "/home/genomics/jorge/files/purple/new/mcl_purple_module_summary.tsv"
MEMBERSHIP_FILE <- "/home/genomics/jorge/files/purple/new/mcl_purple_membership.tsv"
EDGE_FILE       <- "/home/genomics/jorge/files/purple/new/network_purple_filtered_edges.tsv"
OUT_DIR         <- "/home/genomics/jorge/files/purple/new/"

# ==============================================================================
# INPUT FILES SUGARCANE
# ==============================================================================

#EIGENGENE_FILE  <- "/home/genomics/jorge/files/sugarcane/eigengenes_sugarcane.tsv"
#METADATA_FILE   <- "/dados02/jorge/comparative_saccharum/samplesheet.csv"
#SUMMARY_FILE    <- "/home/genomics/jorge/files/sugarcane/mcl_sugarcane_module_summary.tsv"
#MEMBERSHIP_FILE <- "/home/genomics/jorge/files/sugarcane/mcl_sugarcane_membership.tsv"
#EDGE_FILE       <- "/home/genomics/jorge/files/sugarcane/network_sugarcane_filtered_edges.tsv"
#OUT_DIR         <- "/home/genomics/jorge/files/sugarcane/"


# ==============================================================================
# LOAD DATA
# ==============================================================================
message("Loading inputs...")

eigengenes <- fread(EIGENGENE_FILE, header = TRUE)
modules    <- setdiff(names(eigengenes), "sample_id")
cat(sprintf("  Eigengenes: %d modules × %d samples\n", length(modules), nrow(eigengenes)))

meta <- fread(METADATA_FILE, header = TRUE)
setnames(meta, tolower(names(meta)))
cat(sprintf("  Metadata: %d samples\n", nrow(meta)))

summary_dt <- fread(SUMMARY_FILE, header = TRUE)
setnames(summary_dt, tolower(names(summary_dt)))
n_genes_dt <- summary_dt[module != "Unassigned", .(module, n_genes)]

# ==============================================================================
# ENCODE TRAITS
# ==============================================================================
for (tr in names(TRAIT_ENCODING)) {
  enc <- TRAIT_ENCODING[[tr]]
  meta[, (paste0(tr, "_num")) := enc[get(tr)]]
}

trait_num_cols <- paste0(names(TRAIT_ENCODING), "_num")
join_cols      <- c("sample", trait_num_cols)

# Validate all samples have known trait values
for (col in trait_num_cols) {
  n_na <- meta[is.na(get(col)), .N]
  if (n_na > 0)
    warning(sprintf("%d samples with unrecognised %s — will be excluded from that correlation.", n_na, col))
}

# ── Join metadata onto eigengene rows ──────────────────────────────────────
joined <- merge(eigengenes,
                meta[, ..join_cols],
                by.x = "sample_id", by.y = "sample",
                sort = FALSE)

if (nrow(joined) != nrow(eigengenes))
  warning(sprintf("Sample ID mismatch: %d eigengene rows, %d after join — check sample names.",
                  nrow(eigengenes), nrow(joined)))

cat(sprintf("  Samples matched: %d\n", nrow(joined)))

# ==============================================================================
# MODULE × TRAIT CORRELATIONS
# ==============================================================================
message("Computing Pearson correlations...")

cor_results <- rbindlist(lapply(names(TRAIT_ENCODING), function(tr) {
  trait_vec <- joined[[paste0(tr, "_num")]]

  rbindlist(lapply(modules, function(mod) {
    eg_vec <- joined[[mod]]
    ct     <- cor.test(eg_vec, trait_vec, method = "pearson")
    data.table(
      module  = mod,
      trait   = tr,
      pearson = round(unname(ct$estimate), 4),
      pval    = ct$p.value
    )
  }))
}))

# BH correction within each trait independently
cor_results[, padj := round(p.adjust(pval, method = "BH"), 6), by = trait]

# Attach n_genes
cor_results <- merge(cor_results, n_genes_dt, by = "module", all.x = TRUE)
setorder(cor_results, trait, padj)
setcolorder(cor_results, c("module", "trait", "pearson", "pval", "padj", "n_genes"))

out_cor <- file.path(OUT_DIR, sprintf("module_trait_correlations_%s.tsv", NETWORK_NAME))
fwrite(cor_results, file = out_cor, sep = "\t", quote = FALSE)
cat("Saved:", basename(out_cor), "\n")

# Print a quick summary per trait
for (tr in names(TRAIT_ENCODING)) {
  sig <- cor_results[trait == tr & padj <= PADJ_THR, .N]
  cat(sprintf("  %-12s  %d modules with padj <= %.2f\n", tr, sig, PADJ_THR))
}

# ==============================================================================
# SELECT MODULES — genotype, |r| >= PEARSON_THR, padj <= PADJ_THR, n_genes > MIN_GENES
# ==============================================================================
message("\nSelecting modules...")

selected <- cor_results[
  trait        == "treatment" &   ## THISSSS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  abs(pearson) >= PEARSON_THR &
  padj         <= PADJ_THR &
  n_genes      >  MIN_GENES
][order(-abs(pearson))]

cat(sprintf("  Modules passing filters (|r|>=%.1f, padj<=%.2f, n_genes>%d): %d\n",
            PEARSON_THR, PADJ_THR, MIN_GENES, nrow(selected)))

if (nrow(selected) > 0L) {
  print(selected[, .(module, pearson, padj, n_genes)], row.names = FALSE)
}

out_sel <- file.path(OUT_DIR, sprintf("selected_modules_genotype_%s.tsv", NETWORK_NAME))
fwrite(selected, file = out_sel, sep = "\t", quote = FALSE)
cat("Saved:", basename(out_sel), "\n")

if (nrow(selected) == 0L) {
  message("No modules passed the filters — subnetwork step skipped.")
  quit(save = "no")
}

# ==============================================================================
# SUBNETWORK — all edges whose both endpoints belong to the selected modules
# ==============================================================================
message("\nBuilding subnetwork...")

sel_modules  <- selected$module

membership <- fread(MEMBERSHIP_FILE, header = TRUE)
setnames(membership, tolower(names(membership)))
sel_genes_dt <- membership[module_name %in% sel_modules, .(gene, module_name)]
sel_genes    <- sel_genes_dt$gene

cat(sprintf("  Genes in selected modules: %s\n", format(length(sel_genes), big.mark = ",")))

# Index genes for fast lookup
sel_set <- data.table(gene = sel_genes, key = "gene")

cat("  Reading edge list...\n")
edges <- fread(EDGE_FILE, header = TRUE)
cat(sprintf("  Total edges before subsetting: %s\n", format(nrow(edges), big.mark = ",")))

edges_sub <- edges[gene1 %in% sel_genes & gene2 %in% sel_genes]
rm(edges); gc()

cat(sprintf("  Subnetwork edges: %s\n", format(nrow(edges_sub), big.mark = ",")))

# Annotate endpoints with their module (useful for visualisation in Cytoscape etc.)
gene_to_mod <- setNames(sel_genes_dt$module_name, sel_genes_dt$gene)
edges_sub[, module_gene1 := gene_to_mod[gene1]]
edges_sub[, module_gene2 := gene_to_mod[gene2]]

out_sub <- file.path(OUT_DIR, sprintf("subnetwork_selected_modules_%s.tsv", NETWORK_NAME))
fwrite(edges_sub, file = out_sub, sep = "\t", quote = FALSE)
cat("Saved:", basename(out_sub), "\n")

message("\nDone.")
