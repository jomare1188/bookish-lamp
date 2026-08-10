library(data.table)
library(cogeqc)

setDTthreads(100)

# ==============================================================================
# NODE-LEVEL CONSERVATION OF THE NITROGEN RESPONSE
# ==============================================================================
# The edge-level analysis asked whether conserved co-expression LINKS carry the
# nitrogen signal on both sides (answer: ~0 -- the networks rewire heavily).
# This script asks the complementary, node-level question:
#
#   For a gene that correlates with nitrogen in one network, does it have an
#   ortholog in the other network, and is that ortholog ALSO nitrogen-correlated?
#
# This measures conservation of the RESPONSE at the gene, independent of whether
# the surrounding network topology is conserved. The two can diverge: a gene's
# nitrogen response can be conserved across species even when its co-expression
# neighbourhood is not.
#
# Every correlated gene is classified into one of:
#   - no_ortholog              : no ortholog of it exists in the other network
#   - ortholog_not_correlated  : ortholog(s) exist, but none are nitrogen-correlated
#   - conserved_correlated     : >=1 ortholog exists AND is nitrogen-correlated  (the hit)
#
# so a zero (or small) result can be attributed to orthology coverage vs. a real
# lack of shared response, not left ambiguous.
# ==============================================================================

# ==============================================================================
# CONFIGURATION
# ==============================================================================
BASE <- "/dados04/jorge/comparative_saccharum/files/"
OUT_DIR <- paste0(BASE, "network_conservation/")

ORTHOFINDER_TSV <- "/dados04/jorge/comparative_saccharum/files/fix_orthofinder/proteins/OrthoFinder/Results_Jun04_2/Orthogroups/Orthogroups.tsv"

SPECIES_SUGARCANE <- "sugarcane_one_transcript"
SPECIES_PURPLE    <- "one_transcript_purple_proteins"

SELECTED_SUGARCANE <- paste0(BASE, "sugarcane/selected_genes_treatment_sugarcane.tsv")
SELECTED_PURPLE    <- paste0(BASE, "purple/new/selected_genes_treatment_purple.tsv")

strip_id <- function(x) sub("\\.v[0-9]+\\.[0-9]+$", "", x)

# ==============================================================================
# 1. ORTHOLOG PAIRS (sugarcane <-> purple), keeping orthogroup
# ==============================================================================
message("Loading orthologs...")
og <- as.data.table(read_orthogroups(ORTHOFINDER_TSV))
og[, Gene := strip_id(Gene)]
og[Species == SPECIES_SUGARCANE, Species := "sugarcane"]
og[Species == SPECIES_PURPLE,    Species := "purple"]

sc_og <- og[Species == "sugarcane", .(Orthogroup, sugarcane_gene = Gene)]
pu_og <- og[Species == "purple",    .(Orthogroup, purple_gene    = Gene)]
pairs <- unique(merge(sc_og, pu_og, by = "Orthogroup", allow.cartesian = TRUE))
message(sprintf("  Ortholog pairs (many-to-many): %s", format(nrow(pairs), big.mark = ",")))

# ==============================================================================
# 2. TREATMENT-CORRELATED GENES (with their Pearson r)
# ==============================================================================
sc <- fread(SELECTED_SUGARCANE, header = TRUE)[, .(sugarcane_gene = strip_id(gene), sc_r = pearson)]
pu <- fread(SELECTED_PURPLE,    header = TRUE)[, .(purple_gene    = strip_id(gene), pu_r = pearson)]
sc <- unique(sc); pu <- unique(pu)
message(sprintf("  Correlated genes -- sugarcane: %s | purple: %s",
                format(nrow(sc), big.mark = ","), format(nrow(pu), big.mark = ",")))

# ==============================================================================
# 3. COVERAGE DIAGNOSTICS  (distinguish 'no ortholog' from 'not correlated')
# ==============================================================================
sc_has_any_ortho <- sc$sugarcane_gene %in% pairs$sugarcane_gene
pu_has_any_ortho <- pu$purple_gene    %in% pairs$purple_gene
message(sprintf("\n  [diag] correlated sugarcane genes with ANY purple ortholog: %d / %d",
                sum(sc_has_any_ortho), nrow(sc)))
message(sprintf("  [diag] correlated purple genes with ANY sugarcane ortholog:  %d / %d",
                sum(pu_has_any_ortho), nrow(pu)))

# ==============================================================================
# 4. CORE INTERSECTION: ortholog pairs correlated on BOTH sides
# ==============================================================================
corr_pairs <- pairs[sugarcane_gene %in% sc$sugarcane_gene &
                    purple_gene    %in% pu$purple_gene]
corr_pairs <- merge(corr_pairs, sc, by = "sugarcane_gene")
corr_pairs <- merge(corr_pairs, pu, by = "purple_gene")

# direction concordance (positive r => expression rises with nitrogen on both
# sides; note purple uses graded 0N/2N/6N, sugarcane High/Low -- both oriented
# so that higher = more nitrogen, so signs are comparable)
corr_pairs[, concordant := sign(sc_r) == sign(pu_r)]
setcolorder(corr_pairs,
            c("Orthogroup", "sugarcane_gene", "sc_r", "purple_gene", "pu_r", "concordant"))
setorder(corr_pairs, -concordant, Orthogroup)

fwrite(corr_pairs, file.path(OUT_DIR, "conserved_correlated_ortholog_pairs.tsv"),
       sep = "\t", quote = FALSE)

# ==============================================================================
# 5. PER-GENE CLASSIFICATION (each side)
# ==============================================================================
classify <- function(sel, self_col, other_col, pairs, corr_pairs) {
  genes_with_ortho <- unique(pairs[[self_col]])
  genes_conserved  <- unique(corr_pairs[[self_col]])
  n_ortho    <- pairs[, .N, by = self_col]        # orthologs per gene
  n_corr_ort <- corr_pairs[, .(n_corr = uniqueN(get(other_col))), by = self_col]

  out <- copy(sel)
  out <- merge(out, n_ortho,    by = self_col, all.x = TRUE)
  out <- merge(out, n_corr_ort, by = self_col, all.x = TRUE)
  setnames(out, "N", "n_orthologs")
  out[is.na(n_orthologs), n_orthologs := 0L]
  out[is.na(n_corr),      n_corr      := 0L]
  out[, status := fifelse(n_orthologs == 0L, "no_ortholog",
                  fifelse(n_corr == 0L, "ortholog_not_correlated",
                          "conserved_correlated"))]
  out[]
}

sc_status <- classify(sc, "sugarcane_gene", "purple_gene",    pairs, corr_pairs)
pu_status <- classify(pu, "purple_gene",    "sugarcane_gene", pairs, corr_pairs)

fwrite(sc_status, file.path(OUT_DIR, "sugarcane_correlated_conservation_status.tsv"),
       sep = "\t", quote = FALSE)
fwrite(pu_status, file.path(OUT_DIR, "purple_correlated_conservation_status.tsv"),
       sep = "\t", quote = FALSE)

# ==============================================================================
# 6. SUMMARY
# ==============================================================================
sc_tab <- sc_status[, .N, by = status]
pu_tab <- pu_status[, .N, by = status]

message("\n=== sugarcane correlated genes ===")
print(sc_tab)
message("\n=== purple correlated genes ===")
print(pu_tab)

summary_dt <- data.table(
  metric = c("correlated_sugarcane_genes",
             "correlated_purple_genes",
             "sugarcane_genes_conserved_correlated",
             "purple_genes_conserved_correlated",
             "conserved_correlated_ortholog_pairs",
             "conserved_correlated_orthogroups",
             "pairs_sign_concordant",
             "pairs_sign_discordant"),
  value  = c(nrow(sc),
             nrow(pu),
             sc_status[status == "conserved_correlated", .N],
             pu_status[status == "conserved_correlated", .N],
             nrow(corr_pairs),
             uniqueN(corr_pairs$Orthogroup),
             corr_pairs[concordant == TRUE,  .N],
             corr_pairs[concordant == FALSE, .N])
)
fwrite(summary_dt, file.path(OUT_DIR, "conserved_correlated_summary.tsv"),
       sep = "\t", quote = FALSE)

message("\nSummary:")
print(summary_dt)
message("\nSaved to: ", OUT_DIR)
message("  conserved_correlated_ortholog_pairs.tsv        (the hits, with both r values + concordance)")
message("  sugarcane_correlated_conservation_status.tsv   (per-gene: no_ortholog / ortholog_not_correlated / conserved_correlated)")
message("  purple_correlated_conservation_status.tsv")
message("  conserved_correlated_summary.tsv")
