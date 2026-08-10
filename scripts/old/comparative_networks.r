library(data.table)
library(Matrix)
library(cogeqc)

# ==============================================================================
# CONFIGURATION
# ==============================================================================
ORTHOFINDER_TSV  <- "/home/genomics/jorge/files/fix_orthofinder/proteins/OrthoFinder/Results_Jun04_2/Orthogroups/Orthogroups.tsv"
SUBNET_SUGARCANE <- "/home/genomics/jorge/files/sugarcane/subnetwork_selected_modules_sugarcane.tsv"
SUBNET_PURPLE    <- "/home/genomics/jorge/files/purple/new/subnetwork_selected_modules_purple.tsv"
OUT_DIR          <- "/home/genomics/jorge/files/network_conservation/"

# ── Species detection ──────────────────────────────────────────────────────────
# Run once with DETECT_SPECIES = TRUE to print species names, then set to FALSE
# and fill in the two strings below (protein filename without extension).
DETECT_SPECIES    <- FALSE
SPECIES_SUGARCANE <- "sugarcane_one_transcript"   # update after detection
SPECIES_PURPLE    <- "one_transcript_purple_proteins"       # update after detection

# Strip the .p<N> isoform suffix OrthoFinder appends to protein gene IDs?
# Set FALSE if your gene IDs already match between the network and OrthoFinder.
STRIP_ISOFORM <- TRUE

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

#if (STRIP_ISOFORM) og[, Gene := sub("\\v2.1", "", Gene)]
if (STRIP_ISOFORM)  og[, Gene := sub("\\.v[0-9]+\\.[0-9]+$", "", Gene)]

og[Species == SPECIES_SUGARCANE, Species := "sugarcane_one_transcript"]
og[Species == SPECIES_PURPLE,    Species := "one_transcript_purple_proteins"]

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

# Gene ID format sanity check — compare a few IDs from each source
cat("\n  Sample IDs — sugarcane network:   ", head(edges_sc$gene1, 3), "\n")
cat("  Sample IDs — purple network:      ", head(edges_pu$gene1, 3), "\n")
cat("  Sample IDs — sugarcane orthologs: ", head(ortholog_pairs$sugarcane_gene, 3), "\n")
cat("  Sample IDs — purple orthologs:    ", head(ortholog_pairs$purple_gene, 3), "\n")
cat("  If these formats differ, adjust STRIP_ISOFORM or add a manual gsub() below.\n\n")

# ==============================================================================
# CORE: map edges from network A → network B through ortholog matrix
#
# pairs must have columns gene_a (source species) and gene_b (target species).
# Returns edges_a with an extra logical column `conserved`.
#
# Logic mirrors simple_code.r / shared_edges.r exactly:
#   O         — sparse ortholog matrix: |genes_a| × |genes_b|
#   A_a       — symmetric adjacency of network A (unweighted)
#   A_b       — symmetric adjacency of network B (unweighted)
#   mapped    — (O %*% A_b %*% t(O)) > 0  → projects B edges into A-gene space;
#               entry [i,j] is TRUE iff at least one ortholog pair of gene_i and
#               gene_j is connected in B
#   conserved — mapped & (A_a > 0)  → edges that exist in A AND whose orthologs
#               are connected in B
# ==============================================================================
map_conserved_edges <- function(edges_a, edges_b, pairs,
                                label_a = "A", label_b = "B") {

  # Gene universes: NETWORK GENES ONLY — not inflated with ortholog-table genes
  genes_a <- unique(c(edges_a$gene1, edges_a$gene2))
  genes_b <- unique(c(edges_b$gene1, edges_b$gene2))
  n_a     <- length(genes_a);  n_b <- length(genes_b)
  idx_a   <- setNames(seq_len(n_a), genes_a)
  idx_b   <- setNames(seq_len(n_b), genes_b)

  # Filter ortholog pairs to genes present in BOTH networks
  pairs_filt <- pairs[gene_a %in% genes_a & gene_b %in% genes_b]
  message(sprintf("  [%s→%s] Ortholog pairs overlapping both networks: %s",
                  label_a, label_b, format(nrow(pairs_filt), big.mark = ",")))
  if (nrow(pairs_filt) == 0L)
    stop("No ortholog pairs overlap with both gene sets — check gene ID formatting.")

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

  # Project B edges into A-gene space, then AND with A's own adjacency
  mapped            <- (O %*% A_b %*% t(O)) > 0
  conserved_matrix  <- mapped & (A_a > 0)

  # Annotate each edge in edges_a
  out <- copy(edges_a)
  out[, conserved := as.logical(conserved_matrix[cbind(i_a, j_a)])]
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
