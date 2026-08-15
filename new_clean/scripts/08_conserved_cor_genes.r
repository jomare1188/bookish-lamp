# =============================================================================
# 08_conserved_cor_genes.r — conservation of the nitrogen response
#
# 06_conservation_join.r asks whether an EDGE is conserved. This asks whether the
# nitrogen RESPONSE on it is, at two levels:
#
#   NODE level  which ortholog pairs are trait-responsive in BOTH species, and do
#               they move in the same direction?
#   EDGE level  which CONSERVED EDGES join two genes that are responsive in both
#               species -- and which layer (pearson / mi / both) found them?
#
# The edge level is the one the MI layer was added for. A gene pair related by a
# saturating or threshold nitrogen response is invisible to Pearson, so a
# Pearson-only search can return nothing while a real signal is present. Running
# the same funnel under all three selection rules is what distinguishes "there is
# no shared response" from "Pearson cannot see the shared response".
#
# SELECTION RULE -- `CLEAN_SELECTION`, one of:
#   pearson   |r| >= TRAIT_R_THR and padj <= TRAIT_PADJ_THR   (the historical rule)
#   mi        gene-trait mutual information, padj <= TRAIT_PADJ_THR
#   union     either                                          (default)
#
# Both statistics are read from 12_gene_trait_mi.py's output, which computes them
# on the SAME samples and the same VST, so the two selections are on identical
# footing and switching between them cannot introduce a provenance difference.
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
# Pearson r's are directly comparable. MI is unsigned, so a pair selected only by
# MI has no direction to compare and is reported as NA rather than concordant.
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
SELECTION    <- env_opt("CLEAN_SELECTION", "union")
# DIRECTED test: correct the confirmation species' p-values over only the
# orthologs of the discovery species' responsive genes, instead of genome-wide.
DIRECTED     <- env_flag("CLEAN_DIRECTED", FALSE)
DISCOVERY    <- env_opt("CLEAN_DISCOVERY", "sugarcane")
R_THR        <- env_num("CLEAN_TRAIT_R_THR", 0.6)
PADJ_THR     <- env_num("CLEAN_TRAIT_PADJ_THR", 0.05)
SPECIES_SUGARCANE <- env_opt("CLEAN_OG_SPECIES_SUGARCANE", "sugarcane_one_transcript")
SPECIES_PURPLE    <- env_opt("CLEAN_OG_SPECIES_PURPLE", "one_transcript_purple_proteins")
setDTthreads(as.integer(env_num("CLEAN_CORES", 100)))

if (!SELECTION %in% c("pearson", "mi", "union"))
  stop("CLEAN_SELECTION must be pearson, mi or union (got '", SELECTION, "')",
       call. = FALSE)

TAG <- if (DIRECTED) "_directed" else ""
banner(paste0("conserved nitrogen response  [selection: ", SELECTION,
              if (DIRECTED) paste0(", DIRECTED from ", DISCOVERY) else "", "]"))

# --- inputs ------------------------------------------------------------------
trait_file <- function(s) file.path(RESULTS, s, sprintf("gene_trait_mi_%s.tsv", s))
cons_file  <- function(s) file.path(OUT_DIR, sprintf("conserved_genes_%s_FULL.txt", s))
for (s in c("sugarcane", "purple")) {
  if (!file.exists(trait_file(s)))
    stop("missing ", basename(trait_file(s)),
         "\n  run  ./run.sh traitmi ", s, "  first -- it carries BOTH the Pearson",
         "\n  and the MI statistic, computed on the same samples.", call. = FALSE)
  if (!file.exists(cons_file(s)))
    stop("missing ", basename(cons_file(s)), "\n  run  ./run.sh conserve  first",
         call. = FALSE)
}

# --- orthologs ---------------------------------------------------------------
say("loading orthologs")
og <- as.data.table(read_orthogroups_tsv(ORTHOGROUPS))
og[, Gene := strip_version(sub("\\.p[0-9]+$", "", Gene))]
og[Species == SPECIES_SUGARCANE, Species := "sugarcane"]
og[Species == SPECIES_PURPLE,    Species := "purple"]
pairs <- unique(merge(
  og[Species == "sugarcane", .(Orthogroup, sugarcane_gene = Gene)],
  og[Species == "purple",    .(Orthogroup, purple_gene    = Gene)],
  by = "Orthogroup", allow.cartesian = TRUE))
say("  ortholog pairs (many-to-many): ", fmt_n(nrow(pairs)))

# --- responsive genes, restricted to the conserved-edge gene set -------------
responsive <- function(study, restrict = NULL) {
  keep <- strip_version(readLines(cons_file(study), warn = FALSE))
  t <- fread(trait_file(study))
  t[, gene := strip_version(gene)]
  t <- t[gene %chin% keep]
  if (!is.null(restrict)) {
    # The genome-wide BH denominator is the wrong burden for a comparative
    # question. "Does THIS gene's ortholog also respond?" is a test over the
    # candidate orthologs, not over the transcriptome, and correcting over
    # 44,118 genes when the hypothesis names a few thousand throws away most of
    # the power for nothing. Recompute BH over the restricted set.
    n_before <- nrow(t)
    t <- t[gene %chin% restrict]
    t[, padj         := p.adjust(pval,         method = "BH")]
    t[, pearson_padj := p.adjust(pearson_pval, method = "BH")]
    say(sprintf("  %-9s DIRECTED: %s candidate orthologs of %s's responsive genes",
                study, fmt_n(nrow(t)), DISCOVERY))
    say(sprintf("            BH denominator %s -> %s", fmt_n(n_before), fmt_n(nrow(t))))
  }
  t[, by_pearson := !is.na(pearson_padj) & pearson_padj <= PADJ_THR &
                    abs(pearson) >= R_THR]
  t[, by_mi      := !is.na(padj) & padj <= PADJ_THR]
  t[, selected := switch(SELECTION,
                         pearson = by_pearson,
                         mi      = by_mi,
                         union   = by_pearson | by_mi)]
  say(sprintf("  %-9s conserved genes %s | pearson %s | mi %s | union %s -> selected %s",
              study, fmt_n(nrow(t)), fmt_n(sum(t$by_pearson)), fmt_n(sum(t$by_mi)),
              fmt_n(sum(t$by_pearson | t$by_mi)), fmt_n(sum(t$selected))))
  t[selected == TRUE]
}
say("selecting responsive genes")
if (DIRECTED) {
  CONFIRM <- setdiff(c("sugarcane", "purple"), DISCOVERY)
  say("DIRECTED design: discover in ", DISCOVERY, ", confirm in ", CONFIRM)
  disc <- responsive(DISCOVERY)
  cand <- if (DISCOVERY == "sugarcane")
            pairs[sugarcane_gene %chin% disc$gene, unique(purple_gene)]
          else
            pairs[purple_gene %chin% disc$gene, unique(sugarcane_gene)]
  conf <- responsive(CONFIRM, restrict = cand)
  sc_t <- if (DISCOVERY == "sugarcane") disc else conf
  pu_t <- if (DISCOVERY == "sugarcane") conf else disc
} else {
  sc_t <- responsive("sugarcane")
  pu_t <- responsive("purple")
}
sc <- unique(sc_t[, .(sugarcane_gene = gene, sc_r = pearson, sc_mi = mi,
                      sc_by_pearson = by_pearson, sc_by_mi = by_mi)])
pu <- unique(pu_t[, .(purple_gene = gene, pu_r = pearson, pu_mi = mi,
                      pu_by_pearson = by_pearson, pu_by_mi = by_mi)])

# =============================================================================
# NODE level
# =============================================================================
corr_pairs <- pairs[sugarcane_gene %chin% sc$sugarcane_gene &
                    purple_gene    %chin% pu$purple_gene]
corr_pairs <- merge(corr_pairs, sc, by = "sugarcane_gene")
corr_pairs <- merge(corr_pairs, pu, by = "purple_gene")
# Direction is only comparable when BOTH sides were called by Pearson; MI is
# unsigned, so a pair that only MI selected has no direction to agree about.
corr_pairs[, both_by_pearson := sc_by_pearson & pu_by_pearson]
corr_pairs[, concordant := fifelse(both_by_pearson, sign(sc_r) == sign(pu_r), NA)]
setcolorder(corr_pairs, c("Orthogroup", "sugarcane_gene", "sc_r", "sc_mi",
                          "purple_gene", "pu_r", "pu_mi", "concordant"))
setorder(corr_pairs, -concordant, Orthogroup)
write_tsv(corr_pairs, file.path(OUT_DIR,
          sprintf("conserved_correlated_ortholog_pairs_%s%s.tsv", SELECTION, TAG)))

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
write_tsv(sc_status, file.path(OUT_DIR,
          sprintf("sugarcane_correlated_conservation_status_%s%s.tsv", SELECTION, TAG)))
write_tsv(pu_status, file.path(OUT_DIR,
          sprintf("purple_correlated_conservation_status_%s%s.tsv", SELECTION, TAG)))

banner("node level")
print(sc_status[, .N, by = status], row.names = FALSE)
print(pu_status[, .N, by = status], row.names = FALSE)
say("conserved correlated ortholog pairs: ", fmt_n(nrow(corr_pairs)))

# =============================================================================
# EDGE level -- the question the MI layer was added for
# =============================================================================
# A conserved edge counts only if BOTH its genes are responsive in BOTH species.
# The funnel is reported at every step, because where it collapses is the result:
# an empty answer caused by purple having 30 responsive genes is a power finding,
# not a biological one.
banner("edge level")

# genes responsive in this species AND having a responsive ortholog in the other
sc_both <- intersect(sc$sugarcane_gene, pairs[purple_gene %chin% pu$purple_gene, unique(sugarcane_gene)])
say("sugarcane genes responsive here AND with a responsive purple ortholog: ",
    fmt_n(length(sc_both)))

edge_file <- file.path(OUT_DIR, "conserved_edges_sugarcane_to_purple_FULL.tsv")
if (!file.exists(edge_file)) {
  say("NOTE: ", basename(edge_file), " not found — skipping the edge level.")
} else {
  hdr <- strsplit(readLines(edge_file, n = 1L), "\t", fixed = TRUE)[[1L]]
  stopifnot(all(c("gene1", "gene2", "source", "conserved") %in% hdr))

  count_edges <- function(genes, label) {
    if (length(genes) == 0L) {
      say(sprintf("  %-42s %s", label, "0  (no qualifying genes)"))
      return(data.table(source = character(), n = integer()))
    }
    gf <- tempfile(); writeLines(genes, gf); on.exit(unlink(gf), add = TRUE)
    out <- system(sprintf(
      "awk -F'\\t' -v G=%s 'BEGIN{while((getline l<G)>0) r[l]=1} NR>1 && $4==\"TRUE\" && ($1 in r) && ($2 in r){c[$3]++} END{for(k in c) print k\"\\t\"c[k]}' %s",
      shQuote(gf), shQuote(edge_file)), intern = TRUE)
    unlink(gf)
    if (!length(out)) { say(sprintf("  %-42s 0", label)); return(data.table(source = character(), n = integer())) }
    d <- fread(text = paste(out, collapse = "\n"), header = FALSE,
               col.names = c("source", "n"))
    say(sprintf("  %-42s %s   (%s)", label, fmt_n(sum(d$n)),
                paste(sprintf("%s %s", d$source, fmt_n(d$n)), collapse = " | ")))
    d
  }

  say("conserved edges whose BOTH endpoints are...")
  e_sc   <- count_edges(sc$sugarcane_gene, "responsive in sugarcane")
  e_both <- count_edges(sc_both,           "responsive in BOTH species")

  edge_summary <- rbind(
    data.table(level = "both endpoints responsive in sugarcane",
               selection = SELECTION, e_sc),
    data.table(level = "both endpoints responsive in both species",
               selection = SELECTION, e_both), fill = TRUE)
  if (nrow(edge_summary))
    write_tsv(edge_summary, file.path(OUT_DIR,
              sprintf("conserved_correlated_edges_%s%s.tsv", SELECTION, TAG)))

  if (sum(e_both$n) == 0L) {
    say("")
    say("No conserved edge joins two genes responsive in both species.")
    say("An edge needs TWO such genes and there are ", length(sc_both),
        ", so this is arithmetic, not a weak signal: the funnel closes at the")
    say("purple selection (", fmt_n(nrow(pu)), " responsive genes at n = 18), not at the edges.")
  }
}

# =============================================================================
summary_dt <- data.table(
  metric = c("selection_rule", "design",
             "correlated_sugarcane_genes", "correlated_purple_genes",
             "sugarcane_genes_conserved_correlated",
             "purple_genes_conserved_correlated",
             "conserved_correlated_ortholog_pairs",
             "conserved_correlated_orthogroups",
             "pairs_sign_concordant", "pairs_sign_discordant",
             "pairs_direction_not_comparable",
             "genes_responsive_both_species"),
  value  = c(SELECTION, if (DIRECTED) paste0("directed from ", DISCOVERY) else "genome-wide",
             nrow(sc), nrow(pu),
             sc_status[status == "conserved_correlated", .N],
             pu_status[status == "conserved_correlated", .N],
             nrow(corr_pairs), uniqueN(corr_pairs$Orthogroup),
             corr_pairs[concordant %in% TRUE,  .N],
             corr_pairs[concordant %in% FALSE, .N],
             corr_pairs[is.na(concordant), .N],
             length(sc_both)))
write_tsv(summary_dt, file.path(OUT_DIR,
          sprintf("conserved_correlated_summary_%s%s.tsv", SELECTION, TAG)))
print(summary_dt, row.names = FALSE)
say("done  [selection: ", SELECTION, "]")
