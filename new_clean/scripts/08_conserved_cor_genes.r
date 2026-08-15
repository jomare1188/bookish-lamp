# =============================================================================
# 08_conserved_cor_genes.r — node-level conservation of the nitrogen response
#
# 06_conservation_join.r asks whether an EDGE is conserved. This asks whether a
# gene's RESPONSE is: given the trait-correlated genes on each side (from
# 07_gene_trait_cor.r), which ortholog pairs are correlated in BOTH species, and
# do they move in the same direction?
#
# Writes, into results/conservation/:
#   conserved_correlated_ortholog_pairs.tsv    pair-level, with both r's
#   <study>_correlated_conservation_status.tsv per-gene classification
#   conserved_correlated_summary.tsv
#
# WHY THE PER-GENE STATUS COLUMN MATTERS: a gene can fail to be "conserved
# correlated" for two completely different reasons -- it has no ortholog at all,
# or it has one that simply is not correlated. Collapsing those would make a
# small result look like an absence of shared response when it is really an
# absence of orthology coverage. The status column keeps them apart:
#     no_ortholog | ortholog_not_correlated | conserved_correlated
#
# SIGN CONCORDANCE: both studies are encoded so that higher = more nitrogen
# (sugarcane High/Low, purple the 0/2/6 mM gradient), so the signs of the two
# r's are directly comparable.
#
# RUN: through run.sh  ->  ./run.sh conscor
# =============================================================================

suppressMessages(library(data.table))

source(file.path(dirname(sub("--file=", "",
       grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "lib", "common.R"))

RESULTS      <- env_req("CLEAN_RESULTS")
OUT_DIR      <- ensure_dir(env_req("CLEAN_OUT_DIR"))
ORTHOGROUPS  <- env_req("CLEAN_ORTHOGROUPS")
SELECT_TRAIT <- env_opt("CLEAN_SELECT_TRAIT", "treatment")
SPECIES_SUGARCANE <- env_opt("CLEAN_OG_SPECIES_SUGARCANE", "sugarcane_one_transcript")
SPECIES_PURPLE    <- env_opt("CLEAN_OG_SPECIES_PURPLE", "one_transcript_purple_proteins")
setDTthreads(as.integer(env_num("CLEAN_CORES", 100)))

banner("conserved correlated genes")

sel_file <- function(study)
  file.path(RESULTS, study, sprintf("selected_genes_%s_%s.tsv", SELECT_TRAIT, study))

for (s in c("sugarcane", "purple"))
  if (!file.exists(sel_file(s)))
    stop("missing ", basename(sel_file(s)), "\n  run  ./run.sh trait ", s,
         "  first", call. = FALSE)

# --- orthologs ---------------------------------------------------------------
say("loading orthologs")
og <- as.data.table(read_orthogroups_tsv(ORTHOGROUPS))
og[, Gene := sub("\\.p[0-9]+$", "", Gene)]
og[, Gene := strip_version(Gene)]
og[Species == SPECIES_SUGARCANE, Species := "sugarcane"]
og[Species == SPECIES_PURPLE,    Species := "purple"]

pairs <- unique(merge(
  og[Species == "sugarcane", .(Orthogroup, sugarcane_gene = Gene)],
  og[Species == "purple",    .(Orthogroup, purple_gene    = Gene)],
  by = "Orthogroup", allow.cartesian = TRUE))
say("  ortholog pairs (many-to-many): ", fmt_n(nrow(pairs)))

# --- trait-correlated genes on each side -------------------------------------
sc <- unique(fread(sel_file("sugarcane"))[, .(sugarcane_gene = strip_version(gene), sc_r = pearson)])
pu <- unique(fread(sel_file("purple"))[,    .(purple_gene    = strip_version(gene), pu_r = pearson)])
say("  correlated genes: sugarcane ", fmt_n(nrow(sc)), " | purple ", fmt_n(nrow(pu)))

say("  [coverage] correlated sugarcane genes with ANY purple ortholog: ",
    fmt_n(sum(sc$sugarcane_gene %in% pairs$sugarcane_gene)), " / ", fmt_n(nrow(sc)))
say("  [coverage] correlated purple genes with ANY sugarcane ortholog: ",
    fmt_n(sum(pu$purple_gene %in% pairs$purple_gene)), " / ", fmt_n(nrow(pu)))

# --- pairs correlated on BOTH sides ------------------------------------------
corr_pairs <- pairs[sugarcane_gene %in% sc$sugarcane_gene &
                    purple_gene    %in% pu$purple_gene]
corr_pairs <- merge(corr_pairs, sc, by = "sugarcane_gene")
corr_pairs <- merge(corr_pairs, pu, by = "purple_gene")
corr_pairs[, concordant := sign(sc_r) == sign(pu_r)]
setcolorder(corr_pairs, c("Orthogroup", "sugarcane_gene", "sc_r",
                          "purple_gene", "pu_r", "concordant"))
setorder(corr_pairs, -concordant, Orthogroup)
write_tsv(corr_pairs, file.path(OUT_DIR, "conserved_correlated_ortholog_pairs.tsv"))

# --- per-gene classification -------------------------------------------------
classify <- function(sel, self_col, other_col) {
  n_ortho    <- pairs[, .N, by = c(self_col)]
  n_corr_ort <- corr_pairs[, .(n_corr = uniqueN(get(other_col))), by = c(self_col)]
  out <- merge(copy(sel), n_ortho,    by = self_col, all.x = TRUE)
  out <- merge(out,       n_corr_ort, by = self_col, all.x = TRUE)
  setnames(out, "N", "n_orthologs")
  out[is.na(n_orthologs), n_orthologs := 0L]
  out[is.na(n_corr),      n_corr      := 0L]
  out[, status := fifelse(n_orthologs == 0L, "no_ortholog",
                  fifelse(n_corr == 0L, "ortholog_not_correlated",
                          "conserved_correlated"))]
  out[]
}

sc_status <- classify(sc, "sugarcane_gene", "purple_gene")
pu_status <- classify(pu, "purple_gene",    "sugarcane_gene")
write_tsv(sc_status, file.path(OUT_DIR, "sugarcane_correlated_conservation_status.tsv"))
write_tsv(pu_status, file.path(OUT_DIR, "purple_correlated_conservation_status.tsv"))

banner("sugarcane correlated genes"); print(sc_status[, .N, by = status], row.names = FALSE)
banner("purple correlated genes");    print(pu_status[, .N, by = status], row.names = FALSE)

summary_dt <- data.table(
  metric = c("correlated_sugarcane_genes", "correlated_purple_genes",
             "sugarcane_genes_conserved_correlated",
             "purple_genes_conserved_correlated",
             "conserved_correlated_ortholog_pairs",
             "conserved_correlated_orthogroups",
             "pairs_sign_concordant", "pairs_sign_discordant"),
  value  = c(nrow(sc), nrow(pu),
             sc_status[status == "conserved_correlated", .N],
             pu_status[status == "conserved_correlated", .N],
             nrow(corr_pairs), uniqueN(corr_pairs$Orthogroup),
             corr_pairs[concordant == TRUE,  .N],
             corr_pairs[concordant == FALSE, .N]))
write_tsv(summary_dt, file.path(OUT_DIR, "conserved_correlated_summary.tsv"))
print(summary_dt, row.names = FALSE)
say("done")
