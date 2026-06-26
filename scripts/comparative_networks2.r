library(data.table)
library(Matrix)
library(cogeqc)
setDTthreads(100) 
# ==============================================================================
# CONFIGURATION
# ==============================================================================
ORTHOFINDER_TSV  <- "/dados02/jorge/comparative_saccharum/files/fix_orthofinder/proteins/OrthoFinder/Results_Jun04_2/Orthogroups/Orthogroups.tsv"
SUBNET_SUGARCANE <- "/dados02/jorge/comparative_saccharum/files/sugarcane/subnetwork_selected_modules_sugarcane.tsv"
SUBNET_PURPLE    <- "/dados02/jorge/comparative_saccharum/files/purple/new/subnetwork_selected_modules_purple.tsv"
OUT_DIR          <- "/dados02/jorge/comparative_saccharum/files/network_conservation/"

# ── Species detection ──────────────────────────────────────────────────────────
# Run once with DETECT_SPECIES = TRUE to print species names, then set to FALSE
# and fill in the two strings below (protein filename without extension).
DETECT_SPECIES    <- FALSE
SPECIES_SUGARCANE <- "sugarcane"   # update after detection
SPECIES_PURPLE    <- "purple"       # update after detection

# Strip the .p<N> isoform suffix OrthoFinder appends to protein gene IDs?
# Set FALSE if your gene IDs already match between the network and OrthoFinder.
STRIP_ISOFORM <- FALSE

# ==============================================================================
# LOAD ORTHOGROUPS
# ==============================================================================
message("Loading OrthoFinder results...")
og <- as.data.table(read_orthogroups(ORTHOFINDER_TSV))

if (DETECT_SPECIES) {
  cat("\nSpecies names found in the OrthoFinder TSV:\n")
  print(unique(og$Species))
  stop(paste(
    "\n→ Set DETECT_SPECIES = FALSE and update SPECIES_SUGARCANE / SPECIES_PURPLE",
    "with the names printed above, then re-run."
  ))
}

if (STRIP_ISOFORM) og[, Gene := sub("\\.p[0-9]+$", "", Gene)]

og[Species == SPECIES_SUGARCANE, Species := "sugarcane"]
og[Species == SPECIES_PURPLE,    Species := "purple"]

# Build ortholog pairs (allow many-to-many from polyploidy / gene families)
sc_og <- og[Species == "sugarcane_one_transcript", .(Orthogroup, sugarcane_gene = Gene)]
pu_og <- og[Species == "one_transcript_purple_proteins",    .(Orthogroup, purple_gene    = Gene)]
ortholog_pairs <- merge(sc_og, pu_og, by = "Orthogroup", allow.cartesian = TRUE)

message(sprintf("Ortholog pairs (many-to-many included): %s",
                format(nrow(ortholog_pairs), big.mark = ",")))

# ==============================================================================
# LOAD SUBNETWORKS  (only gene1 / gene2 / weight needed)
# ==============================================================================
message("Loading subnetworks...")
edges_sc <- fread(SUBNET_SUGARCANE, select = c("gene1", "gene2", "weight"))
edges_pu <- fread(SUBNET_PURPLE,    select = c("gene1", "gene2", "weight"))

message(sprintf("  Sugarcane: %s edges | %s unique genes",
                format(nrow(edges_sc), big.mark = ","),
                format(length(unique(c(edges_sc$gene1, edges_sc$gene2))), big.mark = ",")))
message(sprintf("  Purple:    %s edges | %s unique genes",
                format(nrow(edges_pu), big.mark = ","),
                format(length(unique(c(edges_pu$gene1, edges_pu$gene2))), big.mark = ",")))

# ==============================================================================
# CORE: map edges from network A → network B through ortholog matrix
#
# pairs must have columns gene_a (source species) and gene_b (target species).
# Returns edges_a with an extra logical column `conserved`.
#
# Logic (from shining-octopus/shared_edges.r):
#   O          — sparse ortholog matrix (|genes_a| × |genes_b|)
#   A_b        — symmetric adjacency of network B (unweighted, for mapping)
#   mapped     — O %*% A_b %*% t(O)  → A-indexed matrix
#   An edge (i,j) in A is conserved when mapped[i,j] > 0, meaning at least one
#   pair of B-orthologs of gene_i and gene_j are connected in network B.
# ==============================================================================
map_conserved_edges <- function(edges_a, edges_b, pairs,
                                label_a = "A", label_b = "B") {

  # Gene universe: network nodes + any gene appearing in the ortholog table
  genes_a <- unique(c(edges_a$gene1, edges_a$gene2, pairs$gene_a))
  genes_b <- unique(c(edges_b$gene1, edges_b$gene2, pairs$gene_b))
  n_a     <- length(genes_a);  n_b <- length(genes_b)
  idx_a   <- setNames(seq_len(n_a), genes_a)
  idx_b   <- setNames(seq_len(n_b), genes_b)

  # Ortholog matrix O (rows = A, cols = B)
  pairs_filt <- pairs[gene_a %in% genes_a & gene_b %in% genes_b]
  message(sprintf("  [%s→%s] Ortholog pairs overlapping both networks: %s",
                  label_a, label_b, format(nrow(pairs_filt), big.mark = ",")))
  if (nrow(pairs_filt) == 0L)
    stop("No ortholog pairs overlap with both gene sets — check gene ID format.")

  O <- sparseMatrix(
    i    = idx_a[pairs_filt$gene_a],
    j    = idx_b[pairs_filt$gene_b],
    x    = 1L,
    dims = c(n_a, n_b)
  )

  # Symmetric (undirected) adjacency of B, unweighted
  i_b <- idx_b[edges_b$gene1];  j_b <- idx_b[edges_b$gene2]
  A_b <- sparseMatrix(
    i    = c(i_b, j_b),
    j    = c(j_b, i_b),
    x    = 1L,
    dims = c(n_b, n_b)
  )

  # Map B edges into A space, classify each edge in A
  mapped    <- (O %*% A_b %*% t(O)) > 0L
  i_a       <- idx_a[edges_a$gene1];  j_a <- idx_a[edges_a$gene2]
  conserved <- as.logical(mapped[cbind(i_a, j_a)])

  out <- copy(edges_a)
  out[, conserved := conserved]
  out
}

# ==============================================================================
# BIDIRECTIONAL CONSERVATION
# ==============================================================================
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

message("\n[1/2] Sugarcane → Purple")
pairs_sc2pu <- ortholog_pairs[, .(gene_a = sugarcane_gene, gene_b = purple_gene)]
res_sc      <- map_conserved_edges(edges_sc, edges_pu, pairs_sc2pu, "sugarcane", "purple")
n_con_sc    <- sum(res_sc$conserved)
jac_sc      <- n_con_sc / nrow(res_sc)
message(sprintf("  Conserved: %s / %s edges | Jaccard: %.4f",
                format(n_con_sc,     big.mark = ","),
                format(nrow(res_sc), big.mark = ","),
                jac_sc))

message("\n[2/2] Purple → Sugarcane")
pairs_pu2sc <- ortholog_pairs[, .(gene_a = purple_gene, gene_b = sugarcane_gene)]
res_pu      <- map_conserved_edges(edges_pu, edges_sc, pairs_pu2sc, "purple", "sugarcane")
n_con_pu    <- sum(res_pu$conserved)
jac_pu      <- n_con_pu / nrow(res_pu)
message(sprintf("  Conserved: %s / %s edges | Jaccard: %.4f",
                format(n_con_pu,     big.mark = ","),
                format(nrow(res_pu), big.mark = ","),
                jac_pu))

# ==============================================================================
# SAVE
# ==============================================================================
fwrite(res_sc,
       file.path(OUT_DIR, "conserved_edges_sugarcane_to_purple.tsv"),
       sep = "\t", quote = FALSE)

fwrite(res_pu,
       file.path(OUT_DIR, "conserved_edges_purple_to_sugarcane.tsv"),
       sep = "\t", quote = FALSE)

summary_dt <- data.table(
  direction       = c("sugarcane_to_purple", "purple_to_sugarcane"),
  total_edges     = c(nrow(res_sc),  nrow(res_pu)),
  conserved_edges = c(n_con_sc,      n_con_pu),
  jaccard_index   = round(c(jac_sc, jac_pu), 4)
)
fwrite(summary_dt,
       file.path(OUT_DIR, "conservation_summary.tsv"),
       sep = "\t", quote = FALSE)

message("\nSummary:")
print(summary_dt)
message("\nOutputs saved to: ", OUT_DIR)
