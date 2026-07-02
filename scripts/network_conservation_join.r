library(data.table)
library(cogeqc)

setDTthreads(100)

# ==============================================================================
# WHY THIS SCRIPT EXISTS
# ==============================================================================
# comparative_networks.r / comparative_networks2.r detect conserved edges via
# sparse matrix multiplication: mapped <- (O %*% A_b %*% t(O)) > 0. That works
# fine at module-subnetwork scale, but it materializes a sparse product over
# the *entire* gene-pair space of network B. At full-network scale -- and
# especially for the purple network (681M edges, mean degree ~8,000) -- that
# product is far too large to be worth computing densely.
#
# This script instead does the same logical check (does at least one ortholog
# pair of gene_i/gene_j form an edge in network B?) as a chunked data.table
# join: for each edge in network A, project gene1/gene2 through orthology,
# then look up the resulting candidate pairs in a keyed adjacency table of B.
# Cost scales with actual ortholog fan-out, not with B's edge density.
#
# Chunking follows the same streaming pattern as build_edgelist.r (read a
# block, process, append to disk, discard) to keep memory bounded regardless
# of total network size.
# ==============================================================================

# ==============================================================================
# CONFIGURATION
# ==============================================================================
ORTHOFINDER_TSV <- "/dados02/jorge/comparative_saccharum/files/fix_orthofinder/proteins/OrthoFinder/Results_Jun04_2/Orthogroups/Orthogroups.tsv"

EDGE_FILE_SUGARCANE <- "/dados02/jorge/comparative_saccharum/files/sugarcane/network_sugarcane_filtered_edges.tsv"
EDGE_FILE_PURPLE    <- "/dados02/jorge/comparative_saccharum/files/purple/new/network_purple_filtered_edges.tsv"

OUT_DIR <- "/dados02/jorge/comparative_saccharum/files/network_conservation/"

SPECIES_SUGARCANE <- "sugarcane_one_transcript"
SPECIES_PURPLE    <- "one_transcript_purple_proteins"

# "sugarcane_to_purple": walk sugarcane edges (75M, the cheaper outer loop),
#   project through orthology, check existence in the purple network.
# "purple_to_sugarcane": walk purple edges (681M) instead -- ~9x more outer
#   edges, run this direction second / separately once the first completes.
DIRECTION <- "purple_to_sugarcane"

CHUNK_SIZE <- 2e6   # edges per chunk -- lower this first if you see memory pressure

# ==============================================================================
# 1. LOAD ORTHOLOGS  (same ID cleaning as comparative_networks.r)
# ==============================================================================
message("Loading OrthoFinder results...")
og <- as.data.table(read_orthogroups(ORTHOFINDER_TSV))
og[, Gene := sub("\\.p[0-9]+$", "", Gene)]
og[, Gene := sub("\\.v[0-9]+\\.[0-9]+$", "", Gene)]

og[Species == SPECIES_SUGARCANE, Species := "sugarcane"]
og[Species == SPECIES_PURPLE,    Species := "purple"]

sc_og <- og[Species == "sugarcane", .(Orthogroup, sugarcane_gene = Gene)]
pu_og <- og[Species == "purple",    .(Orthogroup, purple_gene    = Gene)]
ortholog_pairs <- merge(sc_og, pu_og, by = "Orthogroup", allow.cartesian = TRUE)
message(sprintf("  Ortholog pairs (many-to-many): %s", format(nrow(ortholog_pairs), big.mark = ",")))

# Fan-out diagnostic -- tells you roughly how expensive the join will be,
# and is a quick early signal of ID-format mismatches (fan-out of ~0 would
# mean almost nothing is joining)
fanout_sc <- ortholog_pairs[, .N, by = sugarcane_gene][, mean(N)]
fanout_pu <- ortholog_pairs[, .N, by = purple_gene][, mean(N)]
message(sprintf("  Mean orthologs per sugarcane gene: %.2f", fanout_sc))
message(sprintf("  Mean orthologs per purple gene:    %.2f", fanout_pu))

# ==============================================================================
# 2. SET UP DIRECTION  (A = network walked in chunks, B = lookup target)
# ==============================================================================
if (DIRECTION == "sugarcane_to_purple") {
  edge_file_a <- EDGE_FILE_SUGARCANE
  edge_file_b <- EDGE_FILE_PURPLE
  pairs_ab    <- ortholog_pairs[, .(gene_a = sugarcane_gene, gene_b = purple_gene)]
  label_a <- "sugarcane"; label_b <- "purple"
} else if (DIRECTION == "purple_to_sugarcane") {
  edge_file_a <- EDGE_FILE_PURPLE
  edge_file_b <- EDGE_FILE_SUGARCANE
  pairs_ab    <- ortholog_pairs[, .(gene_a = purple_gene, gene_b = sugarcane_gene)]
  label_a <- "purple"; label_b <- "sugarcane"
} else {
  stop("DIRECTION must be 'sugarcane_to_purple' or 'purple_to_sugarcane'")
}
setkey(pairs_ab, gene_a)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
out_edges   <- file.path(OUT_DIR, sprintf("conserved_edges_%s_to_%s_FULL.tsv", label_a, label_b))
out_genes   <- file.path(OUT_DIR, sprintf("conserved_genes_%s_FULL.txt", label_a))
out_summary <- file.path(OUT_DIR, sprintf("conservation_summary_%s_to_%s_FULL.tsv", label_a, label_b))

# ==============================================================================
# 3. BUILD TARGET ADJACENCY (network B, both directions, keyed for fast join)
# ==============================================================================
message("\nLoading target network (", label_b, ") into adjacency lookup...")
edges_b <- fread(edge_file_b, header = TRUE, select = c("gene1", "gene2"))
message(sprintf("  %s edges loaded", format(nrow(edges_b), big.mark = ",")))

adj_b <- rbindlist(list(
  edges_b[, .(g1 = gene1, g2 = gene2)],
  edges_b[, .(g1 = gene2, g2 = gene1)]
))
rm(edges_b); gc()
setkey(adj_b, g1, g2)
message(sprintf("  Adjacency table (both directions): %s rows", format(nrow(adj_b), big.mark = ",")))

# ==============================================================================
# 4. STREAM NETWORK A IN CHUNKS, PROJECT THROUGH ORTHOLOGS, CHECK AGAINST B
# ==============================================================================
n_total <- tryCatch({
  as.integer(system(sprintf("wc -l < %s", shQuote(edge_file_a)), intern = TRUE)) - 1L
}, error = function(e) NA_integer_)
message(sprintf("\nStreaming network %s (%s edges)...", label_a,
                ifelse(is.na(n_total), "unknown total", format(n_total, big.mark = ","))))

con <- file(edge_file_a, open = "r")
readLines(con, n = 1L, warn = FALSE)   # discard header

writeLines("gene1\tgene2\tconserved", out_edges)
conserved_gene_set <- character(0)
total_conserved    <- 0
total_processed     <- 0
chunk_i             <- 0L
checked_id_format   <- FALSE

repeat {
  lines <- readLines(con, n = CHUNK_SIZE, warn = FALSE)
  if (length(lines) == 0L) break
  chunk_i <- chunk_i + 1L

  chunk <- fread(text = lines, header = FALSE, select = 1:2,
                 col.names = c("gene1", "gene2"))
  rm(lines)
  chunk[, edge_id := .I]

  if (!checked_id_format) {
    cat("\n  Sample IDs -- network ", label_a, ": ", head(chunk$gene1, 3), "\n", sep = "")
    cat("  Sample IDs -- ortholog table (gene_a): ", head(pairs_ab$gene_a, 3), "\n")
    cat("  If these formats clearly differ, stop and fix ID stripping before continuing.\n\n")
    checked_id_format <- TRUE
  }

  # Project gene1 -> orthologs, gene2 -> orthologs (fan-out)
  map1 <- merge(chunk[, .(edge_id, gene1)], pairs_ab,
                by.x = "gene1", by.y = "gene_a", allow.cartesian = TRUE)
  setnames(map1, "gene_b", "ortho1")

  map2 <- merge(chunk[, .(edge_id, gene2)], pairs_ab,
                by.x = "gene2", by.y = "gene_a", allow.cartesian = TRUE)
  setnames(map2, "gene_b", "ortho2")

  # All candidate ortholog-edge combinations per original edge (cross join
  # within each edge_id group -- this is what gives the k1 x k2 fan-out)
  cand <- merge(map1[, .(edge_id, ortho1)], map2[, .(edge_id, ortho2)],
                by = "edge_id", allow.cartesian = TRUE)
  rm(map1, map2)

  # Which candidates actually exist as edges in network B?
  matched       <- merge(cand, adj_b, by.x = c("ortho1", "ortho2"), by.y = c("g1", "g2"))
  conserved_ids <- unique(matched$edge_id)
  rm(cand, matched)

  chunk[, conserved := edge_id %in% conserved_ids]
  total_conserved <- total_conserved + length(conserved_ids)
  total_processed <- total_processed + nrow(chunk)

  fwrite(chunk[, .(gene1, gene2, conserved)], file = out_edges,
         sep = "\t", quote = FALSE, append = TRUE, col.names = FALSE)

  if (any(chunk$conserved)) {
    conserved_gene_set <- union(conserved_gene_set,
                                 unique(c(chunk[conserved == TRUE, gene1],
                                          chunk[conserved == TRUE, gene2])))
  }

  pct <- if (is.na(n_total)) NA else 100 * total_processed / n_total
  message(sprintf("  [%s] Chunk %d -- %s edges processed%s | conserved so far: %s | genes flagged: %s",
                  label_a, chunk_i,
                  format(total_processed, big.mark = ","),
                  ifelse(is.na(pct), "", sprintf(" (%.1f%%)", pct)),
                  format(total_conserved, big.mark = ","),
                  format(length(conserved_gene_set), big.mark = ",")))

  rm(chunk); gc()
}
close(con)

writeLines(conserved_gene_set, out_genes)

jac <- if (total_processed > 0) total_conserved / total_processed else NA
summary_dt <- data.table(
  direction       = sprintf("%s_to_%s", label_a, label_b),
  total_edges     = total_processed,
  conserved_edges = total_conserved,
  conserved_genes = length(conserved_gene_set),
  jaccard_index   = round(jac, 4)
)
fwrite(summary_dt, file = out_summary, sep = "\t", quote = FALSE)

message("\nDone.")
message(sprintf("  Total edges processed: %s", format(total_processed, big.mark = ",")))
message(sprintf("  Conserved edges:       %s (%.4f fraction)", format(total_conserved, big.mark = ","), jac))
message(sprintf("  Genes in conserved edges (%s): %s", label_a, format(length(conserved_gene_set), big.mark = ",")))
message("  Saved: ", out_edges)
message("  Saved: ", out_genes)
message("  Saved: ", out_summary)
