#!/usr/bin/env Rscript
# =============================================================================
# 61_conserved_blocked_nodes.r — which ortholog genes are nitrogen-correlated in
# BOTH species. NODE level, blocked rule only, on the Pearson-only networks.
#
# WHY A NEW SCRIPT AND NOT A BRANCH IN 08. 08_conserved_cor_genes.r answers the
# same question on the MERGED (Pearson + MI) network, under `pearson|mi|union`
# selection, and carries an EDGE level as well. None of that survives the move to
# Pearson-only graphs: there is no MI layer to union, and the edge level needs
# conserved_edges_*_FULL.tsv, which describes the merged graph. 08 stays as the
# record of that analysis; this is the answer on the current one.
#
# THE UNIVERSE. Ortholog pairs whose BOTH sides are nodes of their species'
# Pearson-only network -- the same gene sets the modules are built on, and the same
# sets 31 corrected over. Measured before this script was written:
#
#     sugarcane   101,990 nodes,  65,469 with a purple ortholog
#     purple      170,135 nodes, 109,641 with a sugarcane ortholog
#     43,390 orthogroups have nodes on both sides -> 114,681 candidate PAIRS
#
# FOUR STATUSES, NOT TWO. "Not conserved" has causes that must not be collapsed,
# because collapsing them makes an ORTHOLOGY or NETWORK coverage gap look like a
# biological absence. 08 separated two; a third matters now that the gene universe
# is the network rather than the transcriptome -- only 59.7% of sugarcane's VST
# genes are network nodes, so a gene can have a perfectly good ortholog that simply
# is not in the other species' graph:
#
#     no_ortholog            no ortholog in the orthogroups at all
#     ortholog_not_a_node    has one, but it is not in the other network
#     ortholog_not_correlated in the network, and not nitrogen-correlated
#     conserved_correlated   the answer
#
# PURPLE GETS A SECOND TIER, SUGARCANE CANNOT. Purple has THREE nitrogen levels, so
# a gene moved the same way by deficiency and excess is invisible to every monotone
# test. A purple gene therefore counts as nitrogen-correlated if it is a monotone
# hit (31) OR a quadratic-contrast hit (30) -- two test families, each BH-corrected
# within itself over the SAME genes, then unioned. This script ASSERTS that the two
# gene sets are identical, because a union across two denominators is not a union.
#
# THE ASYMMETRY IS FORCED AND MUST NOT BE READ PAST. Sugarcane's design has TWO
# nitrogen levels and cannot express curvature at all. So a pair admitted by the
# quadratic tier says "purple responds non-monotonically and its ortholog responds
# monotonically in sugarcane" -- never "the two species share a non-monotone
# response". That second question is unanswerable with these two designs.
#
# DIRECTION IS READ OFF THE STATISTIC THAT DID THE SELECTING. Both monotone sides
# carry a signed partial correlation, so sign concordance is defined. A purple side
# admitted ONLY by the quadratic tier has no comparable sign -- a trough at the
# control is neither up nor down -- and is reported as not comparable rather than
# being counted as agreement or disagreement.
#
# RUN: through run.sh  ->  ./run.sh consblocked [0|1]
# =============================================================================
suppressPackageStartupMessages({library(data.table)})
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE),
                                                 value = TRUE)[1])), "lib", "common.R"))

RESULTS      <- env_req("CLEAN_RESULTS")
OUT_DIR      <- ensure_dir(env_req("CLEAN_OUT_DIR"))
ORTHOGROUPS  <- env_req("CLEAN_ORTHOGROUPS")
NODES        <- c(sugarcane = env_req("CLEAN_NODES_SUGARCANE"),
                  purple    = env_req("CLEAN_NODES_PURPLE"))
SELECT_TRAIT <- env_opt("CLEAN_SELECT_TRAIT", "treatment")
DIRECTED     <- env_flag("CLEAN_DIRECTED", FALSE)
DISCOVERY    <- env_opt("CLEAN_DISCOVERY", "sugarcane")
USHAPE_TIER  <- env_flag("CLEAN_USHAPE_TIER", TRUE)
R_THR        <- env_num("CLEAN_TRAIT_R_THR", 0.6)
PADJ_THR     <- env_num("CLEAN_TRAIT_PADJ_THR", 0.05)
NULL_REPS    <- as.integer(env_num("CLEAN_NULL_REPS", 1000))
SEED         <- as.integer(env_num("CLEAN_SEED", 1188))
SPECIES_SUGARCANE <- env_opt("CLEAN_OG_SPECIES_SUGARCANE", "sugarcane_one_transcript")
SPECIES_PURPLE    <- env_opt("CLEAN_OG_SPECIES_PURPLE", "one_transcript_purple_proteins")
setDTthreads(as.integer(env_num("CLEAN_CORES", 100)))
set.seed(SEED)

TAG <- if (DIRECTED) "_directed" else ""
banner(paste0("conserved nitrogen response, NODE level, BLOCKED rule",
              if (DIRECTED) paste0("  [DIRECTED from ", DISCOVERY, "]") else "  [genome-wide]"))

blocked_file <- function(s) file.path(RESULTS, s, sprintf("gene_trait_blocked_%s.tsv", s))
ushape_file  <- function(s) file.path(RESULTS, s, sprintf("gene_trait_ushape_%s.tsv", s))
for (s in c("sugarcane", "purple"))
  if (!file.exists(blocked_file(s)))
    stop("missing ", basename(blocked_file(s)),
         "\n  run  ./run.sh traitblocked ", s, "  first", call. = FALSE)
if (USHAPE_TIER && !file.exists(ushape_file("purple")))
  stop("missing ", basename(ushape_file("purple")),
       "\n  run  ./run.sh ushape purple  first, or set TRAIT_USHAPE_TIER=0",
       call. = FALSE)
say("rule: blocked padj <= ", PADJ_THR, " and |partial r| >= ", R_THR,
    if (USHAPE_TIER) "  (purple also: quadratic-contrast padj <= threshold)" else "")

# --- the two node universes ---------------------------------------------------
nodes <- lapply(NODES, read_gene_universe)
for (s in names(nodes)) say("  ", s, " network nodes: ", fmt_n(length(nodes[[s]])))

# --- orthologs, restricted to pairs the analysis could possibly see ------------
say("")
say("loading orthologs")
og <- as.data.table(read_orthogroups_tsv(ORTHOGROUPS))
og[, Gene := strip_version(sub("\\.p[0-9]+$", "", Gene))]
og[Species == SPECIES_SUGARCANE, Species := "sugarcane"]
og[Species == SPECIES_PURPLE,    Species := "purple"]
pairs_all <- unique(merge(
  og[Species == "sugarcane", .(Orthogroup, sugarcane_gene = Gene)],
  og[Species == "purple",    .(Orthogroup, purple_gene    = Gene)],
  by = "Orthogroup", allow.cartesian = TRUE))
say("  ortholog pairs, unrestricted (many-to-many): ", fmt_n(nrow(pairs_all)))

# A pair only exists for this analysis if BOTH genes are nodes of their network:
# a gene absent from the graph was never clustered and never tested, so a pair
# through it could not have been found however well conserved it is.
pairs <- pairs_all[sugarcane_gene %chin% nodes$sugarcane &
                   purple_gene    %chin% nodes$purple]
say("  pairs with BOTH sides in their network: ", fmt_n(nrow(pairs)),
    "  over ", fmt_n(uniqueN(pairs$Orthogroup)), " orthogroups")
say("  sugarcane nodes with a purple ortholog IN the purple network: ",
    fmt_n(uniqueN(pairs$sugarcane_gene)))
say("  purple nodes with a sugarcane ortholog IN the sugarcane network: ",
    fmt_n(uniqueN(pairs$purple_gene)))

# --- selection -----------------------------------------------------------------
# `restrict` is the DIRECTED design: re-correct the confirmation species over the
# candidate orthologs of the discovery species' responsive genes, not over its whole
# network. The genome-wide denominator is the wrong burden for "does THIS gene's
# ortholog also respond" -- that hypothesis names a few thousand genes, not 170,135.
responsive <- function(study, restrict = NULL) {
  b <- fread(blocked_file(study),
             select = c("gene", "r_partial", "pval", "padj", "r_marginal",
                        "p_marginal", "padj_marginal"))
  b[, gene := strip_version(gene)]
  if (!all(b$gene %chin% nodes[[study]]))
    stop(basename(blocked_file(study)), " holds genes that are not network nodes",
         " -- it was not run on ", basename(NODES[[study]]), call. = FALSE)

  if (!is.null(restrict)) {
    n_before <- nrow(b)
    b <- b[gene %chin% restrict]
    b[, padj          := p.adjust(pval,       method = "BH")]
    b[, padj_marginal := p.adjust(p_marginal, method = "BH")]
    say(sprintf("  %-9s DIRECTED: %s candidate orthologs of %s's responsive genes",
                study, fmt_n(nrow(b)), DISCOVERY))
    say(sprintf("            BH denominator %s -> %s", fmt_n(n_before), fmt_n(nrow(b))))
  }
  b[, by_monotone := !is.na(padj) & padj <= PADJ_THR & abs(r_partial) >= R_THR]
  b[, by_marginal := !is.na(padj_marginal) & padj_marginal <= PADJ_THR &
                     abs(r_marginal) >= R_THR]

  # --- the non-monotone tier, purple only -------------------------------------
  b[, `:=`(by_quadratic = FALSE, u_est = NA_real_, u_padj = NA_real_,
           u_direction = NA_character_)]
  if (USHAPE_TIER && file.exists(ushape_file(study))) {
    u <- fread(ushape_file(study), select = c("gene", "u_est", "u_pval", "u_padj",
                                              "direction"))
    u[, gene := strip_version(gene)]
    setnames(u, "direction", "u_direction")
    # THE TWO FAMILIES MUST SHARE A DENOMINATOR. Each is BH-corrected within itself
    # and the selections are unioned, which is only a union if both corrected over
    # the same genes. 30 takes CLEAN_GENE_FILTER for exactly this reason.
    if (is.null(restrict) && !setequal(u$gene, b$gene))
      stop("the quadratic family (", fmt_n(nrow(u)), " genes) and the monotone ",
           "family (", fmt_n(nrow(b)), " genes) were corrected over DIFFERENT gene ",
           "sets.\n  Re-run  ./run.sh ushape ", study,
           "  so it uses the network node list.", call. = FALSE)
    if (!is.null(restrict)) {
      u <- u[gene %chin% restrict]
      u[, u_padj := p.adjust(u_pval, method = "BH")]
    }
    b[, c("u_est", "u_padj", "u_direction") := NULL]
    b <- merge(b, u[, .(gene, u_est, u_padj, u_direction)], by = "gene", all.x = TRUE)
    b[, by_quadratic := !is.na(u_padj) & u_padj <= PADJ_THR]
    say(sprintf("  %-9s quadratic tier: %s tested, %s selected (trough %s, peak %s)",
                study, fmt_n(sum(!is.na(b$u_padj))), fmt_n(sum(b$by_quadratic)),
                fmt_n(b[by_quadratic & u_direction == "trough_at_control", .N]),
                fmt_n(b[by_quadratic & u_direction == "peak_at_control", .N])))
  } else if (USHAPE_TIER) {
    say(sprintf("  %-9s has no quadratic table -- a two-level design cannot express",
                study), " curvature; monotone only")
  }

  b[, selected := by_monotone | by_quadratic]
  b[, tier := fifelse(by_monotone & by_quadratic, "both",
              fifelse(by_monotone, "monotone",
              fifelse(by_quadratic, "quadratic", "none")))]
  say(sprintf("  %-9s tested %s | monotone %s | quadratic %s | SELECTED %s   (marginal rule: %s)",
              study, fmt_n(nrow(b)), fmt_n(sum(b$by_monotone)),
              fmt_n(sum(b$by_quadratic)), fmt_n(sum(b$selected)),
              fmt_n(sum(b$by_marginal))))
  out <- b[selected == TRUE]
  attr(out, "tested") <- b$gene
  attr(out, "n_marginal") <- sum(b$by_marginal)
  out
}

say("")
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
TESTED_SC <- attr(sc_t, "tested"); TESTED_PU <- attr(pu_t, "tested")
MARG_SC   <- attr(sc_t, "n_marginal"); MARG_PU <- attr(pu_t, "n_marginal")

sc <- unique(sc_t[, .(sugarcane_gene = gene, sc_r = r_partial, sc_padj = padj,
                      sc_tier = tier)])
pu <- unique(pu_t[, .(purple_gene = gene, pu_r = r_partial, pu_padj = padj,
                      pu_tier = tier, pu_u_est = u_est, pu_u_padj = u_padj,
                      pu_u_direction = u_direction)])

# =============================================================================
# NODE LEVEL — the ortholog pairs correlated in BOTH species
# =============================================================================
corr_pairs <- pairs[sugarcane_gene %chin% sc$sugarcane_gene &
                    purple_gene    %chin% pu$purple_gene]
corr_pairs <- merge(corr_pairs, sc, by = "sugarcane_gene")
corr_pairs <- merge(corr_pairs, pu, by = "purple_gene")

# Direction: defined only where both sides were selected by the MONOTONE test. A
# purple side admitted only by the quadratic contrast responds at both extremes,
# which is neither "up" nor "down", so there is nothing for a sign to agree with.
corr_pairs[, direction_comparable := sc_tier != "quadratic" & pu_tier != "quadratic"]
corr_pairs[, concordant := fifelse(direction_comparable,
                                   sign(sc_r) == sign(pu_r), NA)]
corr_pairs[, pair_tier := fifelse(pu_tier == "quadratic", "purple_non_monotone",
                                  "monotone_both")]
setcolorder(corr_pairs, c("Orthogroup", "sugarcane_gene", "sc_r", "sc_padj", "sc_tier",
                          "purple_gene", "pu_r", "pu_padj", "pu_tier",
                          "concordant", "pair_tier"))
setorder(corr_pairs, -concordant, Orthogroup)
write_tsv(corr_pairs, file.path(OUT_DIR,
          sprintf("conserved_correlated_pairs_blocked_nodes%s.tsv", TAG)))

# --- per-gene status, four causes kept apart ----------------------------------
classify <- function(sel, self, other) {
  # every ortholog this gene has, before the network restriction
  n_any  <- pairs_all[get(self) %chin% sel[[self]], .N, by = c(self)]
  # ... of which how many are nodes of the other network
  n_node <- pairs[,      .N, by = c(self)]
  n_corr <- corr_pairs[, .(n_corr = uniqueN(get(other))), by = c(self)]
  out <- merge(copy(sel), n_any,  by = self, all.x = TRUE)
  setnames(out, "N", "n_orthologs")
  out <- merge(out, n_node, by = self, all.x = TRUE)
  setnames(out, "N", "n_orthologs_in_network")
  out <- merge(out, n_corr, by = self, all.x = TRUE)
  for (cn in c("n_orthologs", "n_orthologs_in_network", "n_corr"))
    out[is.na(get(cn)), (cn) := 0L]
  out[, status := fifelse(n_orthologs == 0L, "no_ortholog",
                  fifelse(n_orthologs_in_network == 0L, "ortholog_not_a_node",
                  fifelse(n_corr == 0L, "ortholog_not_correlated",
                          "conserved_correlated")))]
  out[]
}
sc_status <- classify(sc, "sugarcane_gene", "purple_gene")
pu_status <- classify(pu, "purple_gene",    "sugarcane_gene")
write_tsv(sc_status, file.path(OUT_DIR,
          sprintf("sugarcane_status_blocked_nodes%s.tsv", TAG)))
write_tsv(pu_status, file.path(OUT_DIR,
          sprintf("purple_status_blocked_nodes%s.tsv", TAG)))

banner("node level")
say("sugarcane responsive genes, by what stopped them:")
print(sc_status[, .N, by = status][order(-N)], row.names = FALSE)
say("purple responsive genes, by what stopped them:")
print(pu_status[, .N, by = status][order(-N)], row.names = FALSE)
say("")
say("conserved correlated ortholog PAIRS: ", fmt_n(nrow(corr_pairs)),
    "  over ", fmt_n(uniqueN(corr_pairs$Orthogroup)), " orthogroups")
say("  of which purple's side is non-monotone only: ",
    fmt_n(corr_pairs[pair_tier == "purple_non_monotone", .N]))
sc_both <- sc_status[status == "conserved_correlated", sugarcane_gene]
pu_both <- pu_status[status == "conserved_correlated", purple_gene]
say("  sugarcane genes with a responsive purple ortholog: ", fmt_n(length(sc_both)))
say("  purple genes with a responsive sugarcane ortholog: ", fmt_n(length(pu_both)))

# =============================================================================
# HOW MANY WOULD CHANCE GIVE? — computed, not asserted
# =============================================================================
# THE NULL is 13_conservation_null.r's, applied to genes instead of edges: permute
# the purple column of the ortholog table. That destroys which purple gene each
# sugarcane gene maps to while preserving exactly
#   - each sugarcane gene's number of orthologs (fan-out)
#   - which sugarcane genes have any ortholog at all (coverage)
#   - the multiset of purple genes used
# so the only thing removed is the assignment. A product-of-margins calculation
# would ignore fan-out entirely, which is why the prose expectations this project
# used to quote (README ~4.4, results.md ~5.5 for the SAME quantity) disagreed with
# each other and with the data.
#
# THE UNIVERSE MATTERS and getting it wrong inflates the fold badly. The null runs
# over pairs whose BOTH sides were actually TESTED. A gene never tested could not
# have been selected, so putting it in the null's universe credits the observation
# with an impossibility. Under DIRECTED this is what keeps the comparison honest:
# purple was only testable on the candidate orthologs.
say("")
univ <- pairs[sugarcane_gene %chin% TESTED_SC & purple_gene %chin% TESTED_PU]
say("ortholog-shuffle null, ", fmt_n(NULL_REPS), " replicates, seed ", SEED)
say("  universe: ", fmt_n(nrow(univ)), " of ", fmt_n(nrow(pairs)),
    " network pairs have both sides in the tested pools",
    if (DIRECTED) "  (DIRECTED: purple restricted to the candidates)" else "")
if (nrow(univ) > nrow(pairs))
  stop("the null universe is larger than the pair universe -- impossible", call. = FALSE)

null_stats <- if (nrow(sc) && nrow(pu) && nrow(univ)) {
  pu_col <- univ$purple_gene
  sc_in  <- univ$sugarcane_gene %chin% sc$sugarcane_gene
  pu_sel <- pu$purple_gene
  vapply(seq_len(NULL_REPS), function(i) {
    hit <- sc_in & (sample(pu_col) %chin% pu_sel)
    c(pairs = sum(hit), genes = uniqueN(univ$sugarcane_gene[hit]))
  }, numeric(2))
} else matrix(0, nrow = 2, ncol = NULL_REPS,
              dimnames = list(c("pairs", "genes"), NULL))

emp_p <- function(obs, nullv) (sum(nullv >= obs) + 1) / (length(nullv) + 1)
exp_pairs <- mean(null_stats["pairs", ]); sd_pairs <- sd(null_stats["pairs", ])
exp_genes <- mean(null_stats["genes", ]); sd_genes <- sd(null_stats["genes", ])
obs_pairs <- nrow(corr_pairs); obs_genes <- length(sc_both)
say(sprintf("  conserved correlated pairs : observed %s   null %.1f +/- %.1f   fold %.2f   p %.4g",
            fmt_n(obs_pairs), exp_pairs, sd_pairs,
            if (exp_pairs > 0) obs_pairs / exp_pairs else NA_real_,
            emp_p(obs_pairs, null_stats["pairs", ])))
say(sprintf("  genes responsive in both   : observed %s   null %.1f +/- %.1f   fold %.2f   p %.4g",
            fmt_n(obs_genes), exp_genes, sd_genes,
            if (exp_genes > 0) obs_genes / exp_genes else NA_real_,
            emp_p(obs_genes, null_stats["genes", ])))

# --- does direction agree more often than a coin? ------------------------------
# TWO TESTS, because the pairs are NOT independent trials. Orthology here is
# many-to-many -- 392 pairs over 270 orthogroups in the genome-wide run -- so a
# binomial over pairs counts one duplicated locus several times, the same
# pseudoreplication 31's plant-level control exists to answer at gene level. The
# ORTHOGROUP-level test is the conservative one: each orthogroup votes once, by the
# majority sign agreement of its comparable pairs, and a tied orthogroup abstains.
conc <- corr_pairs[!is.na(concordant)]
bt <- if (nrow(conc)) binom.test(sum(conc$concordant), nrow(conc), p = 0.5) else NULL
if (nrow(conc))
  say(sprintf("  sign concordance, PAIRS    : %s of %s (%.1f%%)  binom p %.4g",
              fmt_n(sum(conc$concordant)), fmt_n(nrow(conc)),
              100 * mean(conc$concordant), bt$p.value))

og_vote <- if (nrow(conc)) {
  v <- conc[, .(agree = sum(concordant), disagree = sum(!concordant)), by = Orthogroup]
  v[, vote := fifelse(agree > disagree, TRUE, fifelse(agree < disagree, FALSE, NA))]
  v
} else data.table(Orthogroup = character(), vote = logical())
og_ok <- og_vote[!is.na(vote)]
btg <- if (nrow(og_ok)) binom.test(sum(og_ok$vote), nrow(og_ok), p = 0.5) else NULL
if (nrow(og_ok))
  say(sprintf("  sign concordance, ORTHOGROUPS: %s of %s (%.1f%%)  binom p %.4g   [%s tied, abstaining]",
              fmt_n(sum(og_ok$vote)), fmt_n(nrow(og_ok)),
              100 * mean(og_ok$vote), btg$p.value, fmt_n(og_vote[is.na(vote), .N])))
say(sprintf("  direction not comparable    : %s pairs (purple selected only by the ",
            fmt_n(corr_pairs[is.na(concordant), .N])), "quadratic contrast)")

# =============================================================================
summary_dt <- data.table(
  metric = c("rule", "design", "level", "ushape_tier",
             "sugarcane_network_nodes", "purple_network_nodes",
             "ortholog_pairs_unrestricted", "ortholog_pairs_both_in_network",
             "orthogroups_both_in_network",
             "responsive_sugarcane", "responsive_purple",
             "responsive_sugarcane_marginal_rule", "responsive_purple_marginal_rule",
             "purple_selected_quadratic_only",
             "sugarcane_conserved_correlated", "purple_conserved_correlated",
             "conserved_correlated_pairs", "conserved_correlated_orthogroups",
             "pairs_purple_non_monotone",
             "pairs_sign_concordant", "pairs_sign_discordant",
             "pairs_direction_not_comparable", "sign_concordance_binom_p",
             "orthogroups_voting", "orthogroups_concordant", "orthogroups_tied",
             "sign_concordance_orthogroup_binom_p",
             "null_reps", "null_universe_pairs",
             "expected_pairs", "expected_pairs_sd", "fold_over_null_pairs", "p_emp_pairs",
             "expected_genes_both", "expected_genes_both_sd",
             "fold_over_null_genes", "p_emp_genes"),
  value  = c("blocked", if (DIRECTED) paste0("directed from ", DISCOVERY) else "genome-wide",
             "nodes", as.integer(USHAPE_TIER),
             length(nodes$sugarcane), length(nodes$purple),
             nrow(pairs_all), nrow(pairs), uniqueN(pairs$Orthogroup),
             nrow(sc), nrow(pu), MARG_SC, MARG_PU,
             pu[pu_tier == "quadratic", .N],
             length(sc_both), length(pu_both),
             obs_pairs, uniqueN(corr_pairs$Orthogroup),
             corr_pairs[pair_tier == "purple_non_monotone", .N],
             corr_pairs[concordant %in% TRUE, .N],
             corr_pairs[concordant %in% FALSE, .N],
             corr_pairs[is.na(concordant), .N],
             if (!is.null(bt)) signif(bt$p.value, 4) else NA,
             nrow(og_ok), sum(og_ok$vote), og_vote[is.na(vote), .N],
             if (!is.null(btg)) signif(btg$p.value, 4) else NA,
             NULL_REPS, nrow(univ),
             round(exp_pairs, 3), round(sd_pairs, 3),
             if (exp_pairs > 0) round(obs_pairs / exp_pairs, 3) else NA,
             signif(emp_p(obs_pairs, null_stats["pairs", ]), 4),
             round(exp_genes, 3), round(sd_genes, 3),
             if (exp_genes > 0) round(obs_genes / exp_genes, 3) else NA,
             signif(emp_p(obs_genes, null_stats["genes", ]), 4)))
write_tsv(summary_dt, file.path(OUT_DIR,
          sprintf("conserved_correlated_summary_blocked_nodes%s.tsv", TAG)))
print(summary_dt, row.names = FALSE)
say("done  [blocked, nodes", if (DIRECTED) ", directed" else "", "]")
