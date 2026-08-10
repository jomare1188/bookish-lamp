library(data.table)
library(cogeqc)

setDTthreads(100)

# ==============================================================================
# TREATMENT-CORRELATED CONSERVED EDGES -- STRICT vs LOOSE
# ==============================================================================
# Combines the conserved-edge table with the per-network treatment gene lists
# (selected_genes_treatment_<network>.tsv). An edge is treatment-correlated if
# both its genes correlate with nitrogen in BOTH networks. Two readings:
#
#   LOOSE : both A-endpoints are treatment-correlated in network A, and each
#           endpoint has at least one ortholog that is treatment-correlated in
#           network B. (Does NOT require those orthologs to be connected in B.)
#
#   STRICT: same, but additionally the treatment-correlated ortholog pair must
#           itself form an edge in network B -- i.e. the SAME connection is
#           treatment-correlated on both sides, not two unrelated ortholog
#           mappings. STRICT is always a subset of LOOSE.
#
# Both are computed in one pass from a candidate-pair table carrying an
# is_b_edge flag; the strict/loose split is just a filter on that flag.
# ==============================================================================

# ==============================================================================
# CONFIGURATION
# ==============================================================================
DIRECTION <- "purple_to_sugarcane"   # or "purple_to_sugarcane"

BASE      <- "/dados04/jorge/comparative_saccharum/files/"
CONS_DIR  <- paste0(BASE, "network_conservation/")

ORTHOFINDER_TSV <- "/dados04/jorge/comparative_saccharum/files/fix_orthofinder/proteins/OrthoFinder/Results_Jun04_2/Orthogroups/Orthogroups.tsv"

SPECIES_SUGARCANE <- "sugarcane_one_transcript"
SPECIES_PURPLE    <- "one_transcript_purple_proteins"

EDGE_FILE_SUGARCANE <- paste0(BASE, "sugarcane/network_sugarcane_filtered_edges.tsv")
EDGE_FILE_PURPLE    <- paste0(BASE, "purple/new/network_purple_filtered_edges.tsv")

SELECTED_SUGARCANE  <- paste0(BASE, "sugarcane/selected_genes_treatment_sugarcane.tsv")
SELECTED_PURPLE     <- paste0(BASE, "purple/new/selected_genes_treatment_purple.tsv")

# ── DIRECTION SWAP ───────────────────────────────────────────────────────────
if (DIRECTION == "sugarcane_to_purple") {
  CONSERVED_EDGES <- paste0(CONS_DIR, "conserved_edges_full_sugarcane_to_purple.tsv")
  EDGE_FILE_B     <- EDGE_FILE_PURPLE
  SELECTED_A      <- SELECTED_SUGARCANE
  SELECTED_B      <- SELECTED_PURPLE
  label_a <- "sugarcane"; label_b <- "purple"
} else if (DIRECTION == "purple_to_sugarcane") {
  CONSERVED_EDGES <- paste0(CONS_DIR, "conserved_edges_full_purple_to_sugarcane.tsv")
  EDGE_FILE_B     <- EDGE_FILE_SUGARCANE
  SELECTED_A      <- SELECTED_PURPLE
  SELECTED_B      <- SELECTED_SUGARCANE
  label_a <- "purple"; label_b <- "sugarcane"
} else stop("DIRECTION must be 'sugarcane_to_purple' or 'purple_to_sugarcane'")

OUT_DIR <- CONS_DIR

strip_id <- function(x) sub("\\.v[0-9]+\\.[0-9]+$", "", x)

# ==============================================================================
# 1. CONSERVED EDGES (network A), conserved == TRUE
# ==============================================================================
message(sprintf("[%s -> %s]", label_a, label_b))
ce <- fread(CONSERVED_EDGES, header = TRUE)
ce[, conserved := as.logical(conserved)]
ce[, gene1 := strip_id(gene1)]
ce[, gene2 := strip_id(gene2)]
ce <- ce[conserved == TRUE]
ce[, edge_id := .I]
n_conserved <- nrow(ce)
message(sprintf("  Conserved edges: %s", format(n_conserved, big.mark = ",")))

# ==============================================================================
# 2. TREATMENT-CORRELATED GENE SETS
# ==============================================================================
sel_a <- strip_id(fread(SELECTED_A, header = TRUE)$gene)
sel_b <- strip_id(fread(SELECTED_B, header = TRUE)$gene)
message(sprintf("  Treatment genes -- %s: %s | %s: %s",
                label_a, format(length(sel_a), big.mark = ","),
                label_b, format(length(sel_b), big.mark = ",")))

# ==============================================================================
# 3. ORTHOLOG PAIRS (A -> B)
# ==============================================================================
og <- as.data.table(read_orthogroups(ORTHOFINDER_TSV))
og[, Gene := strip_id(Gene)]
og[Species == SPECIES_SUGARCANE, Species := "sugarcane"]
og[Species == SPECIES_PURPLE,    Species := "purple"]
sc_og <- og[Species == "sugarcane", .(Orthogroup, sugarcane_gene = Gene)]
pu_og <- og[Species == "purple",    .(Orthogroup, purple_gene    = Gene)]
op    <- merge(sc_og, pu_og, by = "Orthogroup", allow.cartesian = TRUE)

pairs_ab <- if (DIRECTION == "sugarcane_to_purple")
  unique(op[, .(gene_a = sugarcane_gene, gene_b = purple_gene)]) else
  unique(op[, .(gene_a = purple_gene, gene_b = sugarcane_gene)])

# ==============================================================================
# 4. NETWORK B ADJACENCY (both directions)
# ==============================================================================
eb <- fread(EDGE_FILE_B, header = TRUE, select = c("gene1", "gene2"))
eb[, gene1 := strip_id(gene1)]; eb[, gene2 := strip_id(gene2)]
adj_b <- unique(rbindlist(list(
  eb[, .(g1 = gene1, g2 = gene2)],
  eb[, .(g1 = gene2, g2 = gene1)]
)))
adj_b[, is_b_edge := TRUE]
rm(eb); gc()
setkey(adj_b, g1, g2)

# ==============================================================================
# 5. FILTER STAGES
# ==============================================================================
# 5a. both A-endpoints correlated in network A
ce_a <- ce[gene1 %in% sel_a & gene2 %in% sel_a]
n_both_a <- nrow(ce_a)
message(sprintf("\n  [stage 1] conserved edges w/ both %s endpoints correlated: %s / %s",
                label_a, format(n_both_a, big.mark = ","), format(n_conserved, big.mark = ",")))

if (n_both_a == 0L) {
  message("  -> nothing to cross-check on the ", label_b, " side. Stopping.")
  quit(save = "no")
}

# 5b. project endpoints to B orthologs restricted to B treatment list
pairs_sel <- pairs_ab[gene_b %in% sel_b]
map1 <- merge(ce_a[, .(edge_id, gene1)], pairs_sel,
              by.x = "gene1", by.y = "gene_a", allow.cartesian = TRUE)[, .(edge_id, o1 = gene_b)]
map2 <- merge(ce_a[, .(edge_id, gene2)], pairs_sel,
              by.x = "gene2", by.y = "gene_a", allow.cartesian = TRUE)[, .(edge_id, o2 = gene_b)]

# candidate treatment-correlated ortholog pairs (inner join => both endpoints
# have >=1 correlated ortholog in B). This set defines LOOSE.
cand <- merge(map1, map2, by = "edge_id", allow.cartesian = TRUE)

# flag which candidate pairs are real B edges. Any TRUE for an edge => STRICT.
cand <- merge(cand, adj_b, by.x = c("o1", "o2"), by.y = c("g1", "g2"), all.x = TRUE)
cand[is.na(is_b_edge), is_b_edge := FALSE]

loose_ids  <- unique(cand$edge_id)
strict_ids <- unique(cand[is_b_edge == TRUE, edge_id])

message(sprintf("  [stage 2 - LOOSE ] + each endpoint has a correlated %s ortholog: %s",
                label_b, format(length(loose_ids), big.mark = ",")))
message(sprintf("  [stage 3 - STRICT] + those orthologs form a %s edge:          %s",
                label_b, format(length(strict_ids), big.mark = ",")))

# ==============================================================================
# 6. OUTPUT
# ==============================================================================
gcol <- paste0(label_a, c("_gene1", "_gene2"))
ocol <- paste0(label_b, c("_gene1", "_gene2"))

# (a) full candidate-pair table with is_b_edge flag (both definitions in one file)
detail <- merge(cand, ce_a, by = "edge_id")
setnames(detail, c("gene1", "gene2", "o1", "o2"), c(gcol, ocol))
keep <- c(gcol, ocol, "is_b_edge", if ("weight" %in% names(detail)) "weight")
detail <- unique(detail[, ..keep])
fwrite(detail, file.path(OUT_DIR,
       sprintf("treatment_conserved_pairs_ALL_%s_to_%s.tsv", label_a, label_b)),
       sep = "\t", quote = FALSE)

# (b) LOOSE network-A edge list
ecols <- c(gcol, if ("weight" %in% names(ce_a)) "weight")
loose_edges <- unique(merge(data.table(edge_id = loose_ids), ce_a,
                            by = "edge_id")[, setNames(.SD, c(gcol, if ("weight" %in% names(ce_a)) "weight")),
                            .SDcols = c("gene1", "gene2", if ("weight" %in% names(ce_a)) "weight")])
fwrite(loose_edges, file.path(OUT_DIR,
       sprintf("treatment_conserved_edges_LOOSE_%s.tsv", label_a)),
       sep = "\t", quote = FALSE)

# (c) STRICT network-A edge list
strict_edges <- unique(merge(data.table(edge_id = strict_ids), ce_a,
                             by = "edge_id")[, setNames(.SD, c(gcol, if ("weight" %in% names(ce_a)) "weight")),
                             .SDcols = c("gene1", "gene2", if ("weight" %in% names(ce_a)) "weight")])
fwrite(strict_edges, file.path(OUT_DIR,
       sprintf("treatment_conserved_edges_STRICT_%s.tsv", label_a)),
       sep = "\t", quote = FALSE)

# (d) STRICT cross-species pairs (figure-ready: A edge <-> the B edge)
strict_cols   <- c(gcol, ocol, if ("weight" %in% names(detail)) "weight")
strict_detail <- detail[is_b_edge == TRUE, ..strict_cols]
fwrite(strict_detail, file.path(OUT_DIR,
       sprintf("treatment_conserved_edge_pairs_STRICT_%s_to_%s.tsv", label_a, label_b)),
       sep = "\t", quote = FALSE)

# ==============================================================================
# 7. SUMMARY TABLE
# ==============================================================================
summary_dt <- data.table(
  direction              = sprintf("%s_to_%s", label_a, label_b),
  conserved_edges        = n_conserved,
  both_A_endpoints_corr  = n_both_a,
  loose_edges            = length(loose_ids),
  strict_edges           = length(strict_ids),
  loose_minus_strict     = length(loose_ids) - length(strict_ids)
)
fwrite(summary_dt, file.path(OUT_DIR,
       sprintf("treatment_conserved_summary_%s_to_%s.tsv", label_a, label_b)),
       sep = "\t", quote = FALSE)

message("\nSummary:")
print(summary_dt)
message("\nSaved LOOSE edges, STRICT edges, STRICT pairs, full flagged pair table, and summary to:")
message("  ", OUT_DIR)
