library(data.table)
library(Matrix)
library(cogeqc)
setDTthreads(100)
# ==============================================================================
# FULL-NETWORK EDGE CONSERVATION (matrix approach)
#
# Same logic as the subnetwork script, but run over the entire filtered
# co-expression networks instead of only the nitrogen-correlated modules.
#
# WARNING — SCALE: the triple product O %*% A_b %*% t(O) materializes an
# A-gene x A-gene sparse matrix and the intermediate O %*% A_b can densify
# heavily when B has a high mean degree (purple). Expect large memory use.
# Run sugarcane->purple first (fewer outer edges) and monitor RAM. If it OOMs,
# the matrix approach is not viable at this scale and a chunked join is needed.
# ==============================================================================

# ==============================================================================
# CONFIGURATION
# ==============================================================================
ORTHOFINDER_TSV     <- "/dados04/jorge/comparative_saccharum/files/fix_orthofinder/proteins/OrthoFinder/Results_Jun04_2/Orthogroups/Orthogroups.tsv"
EDGE_FILE_SUGARCANE <- "/dados04/jorge/comparative_saccharum/files/sugarcane/network_sugarcane_filtered_edges.tsv"
EDGE_FILE_PURPLE    <- "/dados04/jorge/comparative_saccharum/files/purple/new/network_purple_filtered_edges.tsv"
OUT_DIR             <- "/dados04/jorge/comparative_saccharum/files/network_conservation/"

# OrthoFinder species labels (protein filename without extension)
SPECIES_SUGARCANE <- "sugarcane_one_transcript"
SPECIES_PURPLE    <- "one_transcript_purple_proteins"

# Strip the .p<N> isoform suffix OrthoFinder appends to protein gene IDs?
STRIP_ISOFORM <- FALSE

# ==============================================================================
# LOAD ORTHOGROUPS
# ==============================================================================
message("Loading OrthoFinder results...")
og <- as.data.table(read_orthogroups(ORTHOFINDER_TSV))

if (STRIP_ISOFORM) og[, Gene := sub("\\.p[0-9]+$", "", Gene)]

# Strip R570 v2.1 assembly suffix (e.g. "GeneID.v2.1" -> "GeneID")
og[, Gene := sub("\\.v[0-9]+\\.[0-9]+$", "", Gene)]

sc_og <- og[Species == SPECIES_SUGARCANE, .(Orthogroup, sugarcane_gene = Gene)]
pu_og <- og[Species == SPECIES_PURPLE,    .(Orthogroup, purple_gene    = Gene)]
ortholog_pairs <- merge(sc_og, pu_og, by = "Orthogroup", allow.cartesian = TRUE)

message(sprintf("Ortholog pairs (many-to-many included): %s",
                format(nrow(ortholog_pairs), big.mark = ",")))

# ==============================================================================
# LOAD FULL EDGE LISTS  (only gene1 / gene2 / weight needed)
# ==============================================================================
message("Loading full edge lists...")
edges_sc <- fread(EDGE_FILE_SUGARCANE, select = c("gene1", "gene2", "weight"))
edges_pu <- fread(EDGE_FILE_PURPLE,    select = c("gene1", "gene2", "weight"))

# Strip .v2.1 from the sugarcane network IDs (match the ortholog table)
edges_sc[, gene1 := sub("\\.v[0-9]+\\.[0-9]+$", "", gene1)]
edges_sc[, gene2 := sub("\\.v[0-9]+\\.[0-9]+$", "", gene2)]

message(sprintf("  Sugarcane: %s edges | %s unique genes",
                format(nrow(edges_sc), big.mark = ","),
                format(length(unique(c(edges_sc$gene1, edges_sc$gene2))), big.mark = ",")))
message(sprintf("  Purple:    %s edges | %s unique genes",
                format(nrow(edges_pu), big.mark = ","),
                format(length(unique(c(edges_pu$gene1, edges_pu$gene2))), big.mark = ",")))

# Gene ID format sanity check
cat("\n  Sample IDs — sugarcane network:   ", head(edges_sc$gene1, 3), "\n")
cat("  Sample IDs — purple network:      ", head(edges_pu$gene1, 3), "\n")
cat("  Sample IDs — sugarcane orthologs: ", head(ortholog_pairs$sugarcane_gene, 3), "\n")
cat("  Sample IDs — purple orthologs:    ", head(ortholog_pairs$purple_gene, 3), "\n\n")

# ==============================================================================
# CORE: map edges from network A -> network B through ortholog matrix
#
#   O         — sparse ortholog matrix: |genes_a| x |genes_b|
#   A_a       — symmetric adjacency of network A (unweighted)
#   A_b       — symmetric adjacency of network B (unweighted)
#   mapped    — (O %*% A_b %*% t(O)) > 0  -> B edges projected into A-gene space
#   conserved — mapped & (A_a > 0)        -> A edges whose orthologs connect in B
# ==============================================================================
map_conserved_edges <- function(edges_a, edges_b, pairs,
                                label_a = "A", label_b = "B") {

  genes_a <- unique(c(edges_a$gene1, edges_a$gene2))
  genes_b <- unique(c(edges_b$gene1, edges_b$gene2))
  n_a     <- length(genes_a);  n_b <- length(genes_b)
  idx_a   <- setNames(seq_len(n_a), genes_a)
  idx_b   <- setNames(seq_len(n_b), genes_b)

  # Per-side overlap diagnostic — distinguishes ID mismatch from a real result
  a_in_og <- sum(genes_a %in% pairs$gene_a)
  b_in_og <- sum(genes_b %in% pairs$gene_b)
  message(sprintf("  [%s->%s] %s genes in ortholog table: %s / %s",
                  label_a, label_b, label_a,
                  format(a_in_og, big.mark = ","), format(n_a, big.mark = ",")))
  message(sprintf("  [%s->%s] %s genes in ortholog table: %s / %s",
                  label_a, label_b, label_b,
                  format(b_in_og, big.mark = ","), format(n_b, big.mark = ",")))

  pairs_filt <- pairs[gene_a %in% genes_a & gene_b %in% genes_b]
  message(sprintf("  [%s->%s] Ortholog pairs overlapping both networks: %s",
                  label_a, label_b, format(nrow(pairs_filt), big.mark = ",")))
  if (nrow(pairs_filt) == 0L) {
    if (a_in_og == 0L || b_in_og == 0L)
      stop(sprintf("Gene ID mismatch: %s side has 0 overlap with the ortholog table.",
                   if (a_in_og == 0L) label_a else label_b))
    stop("Both sides map to the ortholog table individually, but no ortholog PAIR ",
         "connects the two networks — genuine biological result, not an ID problem.")
  }

  message(sprintf("  [%s->%s] Building matrices (n_a=%s, n_b=%s)...",
                  label_a, label_b,
                  format(n_a, big.mark = ","), format(n_b, big.mark = ",")))

  # O: ortholog matrix, rows = A genes, cols = B genes
  O <- sparseMatrix(
    i    = idx_a[pairs_filt$gene_a],
    j    = idx_b[pairs_filt$gene_b],
    x    = 1,
    dims = c(n_a, n_b)
  )

  # A_a: symmetric, unweighted adjacency of network A
  i_a <- idx_a[edges_a$gene1];  j_a <- idx_a[edges_a$gene2]
  A_a <- sparseMatrix(
    i    = c(i_a, j_a),
    j    = c(j_a, i_a),
    x    = 1,
    dims = c(n_a, n_a)
  )

  # A_b: symmetric, unweighted adjacency of network B
  i_b <- idx_b[edges_b$gene1];  j_b <- idx_b[edges_b$gene2]
  A_b <- sparseMatrix(
    i    = c(i_b, j_b),
    j    = c(j_b, i_b),
    x    = 1,
    dims = c(n_b, n_b)
  )

  message(sprintf("  [%s->%s] Step 1: O %%*%% A_b ...", label_a, label_b))
  M <- O %*% A_b
  gc()

  message(sprintf("  [%s->%s] Step 2: (O A_b) %%*%% t(O) ...", label_a, label_b))
  mapped <- (M %*% t(O)) > 0
  rm(M); gc()

  message(sprintf("  [%s->%s] Step 3: AND with A_a ...", label_a, label_b))
  conserved_matrix <- mapped & (A_a > 0)
  rm(mapped, A_a, A_b, O); gc()

  out <- copy(edges_a)
  out[, conserved := as.logical(conserved_matrix[cbind(i_a, j_a)])]
  rm(conserved_matrix); gc()
  out
}

# ==============================================================================
# BIDIRECTIONAL CONSERVATION
# ==============================================================================
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

message("\n[1/2] Sugarcane -> Purple (run first: fewer outer edges)")
pairs_sc2pu <- ortholog_pairs[, .(gene_a = sugarcane_gene, gene_b = purple_gene)]
res_sc      <- map_conserved_edges(edges_sc, edges_pu, pairs_sc2pu, "sugarcane", "purple")
n_con_sc    <- sum(res_sc$conserved)
jac_sc      <- signif(n_con_sc / nrow(res_sc), 4)
message(sprintf("  Conserved: %s / %s edges | Jaccard: %s",
                format(n_con_sc,     big.mark = ","),
                format(nrow(res_sc), big.mark = ","), jac_sc))
fwrite(res_sc, file.path(OUT_DIR, "conserved_edges_full_sugarcane_to_purple.tsv"),
       sep = "\t", quote = FALSE)
rm(res_sc); gc()

message("\n[2/2] Purple -> Sugarcane")
pairs_pu2sc <- ortholog_pairs[, .(gene_a = purple_gene, gene_b = sugarcane_gene)]
res_pu      <- map_conserved_edges(edges_pu, edges_sc, pairs_pu2sc, "purple", "sugarcane")
n_con_pu    <- sum(res_pu$conserved)
jac_pu      <- signif(n_con_pu / nrow(res_pu), 4)
message(sprintf("  Conserved: %s / %s edges | Jaccard: %s",
                format(n_con_pu,     big.mark = ","),
                format(nrow(res_pu), big.mark = ","), jac_pu))
fwrite(res_pu, file.path(OUT_DIR, "conserved_edges_full_purple_to_sugarcane.tsv"),
       sep = "\t", quote = FALSE)

# ==============================================================================
# SUMMARY
# ==============================================================================
summary_dt <- data.table(
  direction       = c("sugarcane_to_purple", "purple_to_sugarcane"),
  total_edges     = c(nrow(edges_sc), nrow(edges_pu)),
  conserved_edges = c(n_con_sc, n_con_pu),
  jaccard_index   = c(jac_sc, jac_pu)
)
fwrite(summary_dt, file.path(OUT_DIR, "conservation_summary_full.tsv"),
       sep = "\t", quote = FALSE)

message("\nSummary:")
print(summary_dt)
message("\nOutputs saved to: ", OUT_DIR)
