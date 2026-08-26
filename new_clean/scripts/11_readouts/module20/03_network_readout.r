# ============================================================================
# 03_network_readout.r — do Module-20 genes keep their characteristics in OUR
# two networks?
#
# Muñoz define Module 20 on their own de novo network, where it is ~75%
# MYB/MYB-related and high-betweenness. Asking "is it MYB-rich and central?"
# against that definition is circular. This script asks the non-circular
# versions:
#
#   1. Do our own TF calls independently agree these genes are MYB?
#   2. Are they central IN OUR NETWORKS — and specifically, do they stand out
#      against the MYB/MYB-related background rather than against all nodes?
#      (Comparing MYBs to all nodes just re-measures "MYBs are well connected".)
#   3. Are they on cross-species conserved edges more than degree-matched nodes?
#   4. Do they still cluster together in our own MCL modules — i.e. does our
#      pipeline re-derive their module from the same reads?
#   5. Does any of it hold in purple, the other study's network?
#
# Degree is reported as a PERCENTILE within each network: the two networks
# differ ~9x in density, and Muñoz's own degrees (2-14) come from a third,
# much sparser network, so raw degrees are not comparable across any of them.
#
# Copies of one transcript are NOT independent (polyploid haplotypes), so every
# test is run at TRANSCRIPT level (one value per Module-20 member, the median of
# its copies) as well as copy level, and the transcript-level result is the one
# to trust.
#
# RUN: /home/genomics/miniconda3/envs/r_env/bin/Rscript 03_network_readout.r
# ============================================================================

suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
set.seed(1)

base   <- "/dados04/jorge/comparative_saccharum"
outdir <- "/dados04/jorge/comparative_saccharum/new_clean/results/readouts/module20"
# Cached inputs from steps 01-02b (sequence-only, network-independent) live in
# the original GET_TFS tree; only 03-05 output goes under results/.
cache  <- file.path(base, "GET_TFS/new/results/module20")

COL    <- c(sugarcane = "#1b7837", purple = "#762a83")
N_PERM <- 20000

cfg <- list(
  sugarcane = list(
    map   = file.path(cache, "module20_map_sugarcane.tsv"),
    nodes = file.path(base, "new_clean/results/sugarcane/network_sugarcane_node_metrics.tsv"),
    strip = "\\.v[0-9.]+$",
    # Overridable, but MCL by default and deliberately so. This readout asks
    # which of OUR modules hold the Module-20 genes, which is a question about a
    # specific clustering -- and run.sh invokes this stage with no environment,
    # so an unqualified path here would silently follow whatever CLUSTERING
    # happened to be set. Point CLEAN_M20_MODS_SUGARCANE at another membership
    # file to re-run the readout against it.
    mods  = Sys.getenv("CLEAN_M20_MODS_SUGARCANE",
              file.path(base, "new_clean/results/sugarcane/mcl_sugarcane_membership.tsv")),
    tf    = "/dados04/jorge/comparative_saccharum/new_clean/results/readouts/get_tfs/sugarcane/TF_in_network.tsv",
    cons  = file.path(base, "new_clean/results/conservation/conserved_genes_sugarcane_FULL.txt")),
  purple = list(
    map   = file.path(cache, "module20_map_purple.tsv"),
    nodes = file.path(base, "new_clean/results/purple/network_purple_node_metrics.tsv"),
    strip = NULL,
    mods  = Sys.getenv("CLEAN_M20_MODS_PURPLE",
              file.path(base, "new_clean/results/purple/mcl_purple_membership.tsv")),
    tf    = "/dados04/jorge/comparative_saccharum/new_clean/results/readouts/get_tfs/purple/TF_in_network.tsv",
    cons  = file.path(base, "new_clean/results/conservation/conserved_genes_purple_FULL.txt"))
)

members <- fread(file.path(cache, "module20_members.tsv"))
MYB_FAMS <- c("MYB", "MYB-related")

per_copy <- list(); bg_all <- list(); tests <- list(); modinfo <- list()

for (sp in names(cfg)) {
  cc <- cfg[[sp]]

  nodes <- fread(cc$nodes)
  if (!is.null(cc$strip)) nodes[, gene := sub(cc$strip, "", gene)]
  nodes[, deg_pct := 100 * frank(degree, ties.method = "average") / .N]
  cons <- fread(cc$cons, header = FALSE)$V1
  nodes[, on_cons := gene %in% cons]

  tf <- fread(cc$tf); setnames(tf, c("gene", "Family", "Type"))
  tf_g <- tf[, .(our_family = paste(sort(unique(Family)), collapse = ";")), by = gene]
  myb_genes <- unique(tf[Family %in% MYB_FAMS, gene])

  mods <- fread(cc$mods)[, .(gene, module_name)]
  if (!is.null(cc$strip)) mods[, gene := sub(cc$strip, "", gene)]

  m <- fread(cc$map)[verdict == "reciprocal"]
  # One row per GENE. A gene can be hit by several Module-20 transcripts (they
  # include two isoforms of one locus, and several members are copies of the
  # same orthologue); keep its best-scoring transcript and record the rest.
  setorder(m, -bitscore)
  m <- m[, .(transcript = transcript[1], at_anchor = at_anchor[1],
             n_transcripts = uniqueN(transcript),
             bitscore = max(bitscore), pident = max(pident)), by = gene]
  m <- merge(m, members[, .(transcript, family, munoz_degree, munoz_betweenness)],
             by = "transcript", all.x = TRUE)
  # locus = transcript without its isoform suffix (_i2/_i4 are one locus)
  m[, locus := sub("_i[0-9]+$", "", transcript)]
  m <- merge(m, nodes[, .(gene, degree, deg_pct, strength, transitivity, on_cons)],
             by = "gene", all.x = TRUE)
  m <- merge(m, tf_g,  by = "gene", all.x = TRUE)
  m <- merge(m, mods,  by = "gene", all.x = TRUE)
  m[, in_network := !is.na(degree)]
  m[is.na(our_family), our_family := "not_a_TF"]
  # our_family may hold several families for one gene (multi-isoform), e.g.
  # "MYB;MYB-related" — split before testing membership.
  m[, our_is_MYB := vapply(strsplit(our_family, ";", fixed = TRUE),
                           function(v) any(v %in% MYB_FAMS), logical(1))]
  m[, species := sp]
  per_copy[[sp]] <- m

  inet <- m[in_network == TRUE]

  # ---- backgrounds ---------------------------------------------------------
  mk <- function(g, lab) { d <- nodes[gene %in% g]
                           if (!nrow(d)) return(NULL)
                           data.table(species = sp, group = lab, gene = d$gene,
                                      degree = d$degree, deg_pct = d$deg_pct,
                                      on_cons = d$on_cons) }
  bg_all[[sp]] <- rbindlist(list(
    mk(nodes$gene,  "All nodes"),
    mk(tf_g$gene,   "All TFs"),
    mk(myb_genes,   "All MYB/MYB-rel"),
    mk(inet$gene,   "Module-20 orthologs")), use.names = TRUE)

  # ---- levels of independence ---------------------------------------------
  # gene   : every mapped haplotype copy (NOT independent - polyploid copies)
  # locus  : Module-20 member, isoforms merged
  # anchor : Arabidopsis orthologue group - the most conservative unit, because
  #          7 of the 8 MYB members resolve to the SAME Arabidopsis gene
  #          (AT5G59780 / AtMYB59), i.e. they are one orthologous group, not 7.
  lo <- inet[, .(deg_pct = median(deg_pct), n = .N), by = .(locus, family)]
  an <- inet[, .(deg_pct = median(deg_pct), n = .N), by = at_anchor]

  myb_bg <- nodes[gene %in% setdiff(myb_genes, inet$gene), deg_pct]
  w_copy <- suppressWarnings(wilcox.test(inet$deg_pct, myb_bg))
  w_lo   <- suppressWarnings(wilcox.test(lo$deg_pct,   myb_bg))
  w_an   <- suppressWarnings(wilcox.test(an$deg_pct,   myb_bg))

  # ---- test 2: conserved edges vs a DEGREE-MATCHED null --------------------
  pool <- nodes[order(degree), .(degree, on_cons)]
  idx  <- sapply(inet$degree, function(d) which.min(abs(pool$degree - d)))
  obs  <- sum(inet$on_cons)
  draws <- replicate(N_PERM, sum(sapply(idx, function(i) {
             lo <- max(1, i - 250); hi <- min(nrow(pool), i + 250)
             pool$on_cons[sample(lo:hi, 1)] })))

  # ---- test 3: is Module 20 recovered as a unit in our own modules? --------
  # A raw "fewer modules than random" count is CONFOUNDED: haplotype copies of
  # one gene co-express and land in the same module regardless of Module 20.
  # The meaningful quantity is co-membership across DIFFERENT Module-20 loci.
  obs_mod <- uniqueN(inet$module_name)
  mn <- nodes[gene %in% myb_genes & gene %in% mods$gene, gene]
  nullmod <- replicate(N_PERM, uniqueN(mods[gene %in% sample(mn, nrow(inet)), module_name]))

  msize   <- mods[, .N, by = module_name]
  cross   <- inet[, .(n_loci = uniqueN(locus), n_genes = .N), by = module_name][n_loci > 1]
  cross   <- merge(cross, msize, by = "module_name", all.x = TRUE)
  setnames(cross, "N", "module_size")
  modinfo[[sp]] <- merge(inet[, .(n_genes = .N, n_loci = uniqueN(locus)), by = module_name],
                         msize, by = "module_name", all.x = TRUE)[order(-n_genes)]
  modinfo[[sp]][, species := sp]; setnames(modinfo[[sp]], "N", "module_size")

  tests[[sp]] <- data.table(
    species = sp,
    n_loci_mapped        = uniqueN(inet$locus),
    n_anchors            = uniqueN(inet$at_anchor),
    n_genes_in_network   = nrow(inet),
    our_MYB_call_pct     = round(100 * mean(inet$our_is_MYB), 1),
    median_deg_pct       = round(median(inet$deg_pct), 1),
    median_deg_pct_MYBbg = round(median(myb_bg), 1),
    p_vs_MYB_gene        = signif(w_copy$p.value, 3),
    p_vs_MYB_locus       = signif(w_lo$p.value, 3),
    p_vs_MYB_anchor      = signif(w_an$p.value, 3),
    on_cons_obs          = obs,
    on_cons_exp_degmatch = round(mean(draws), 2),
    p_cons_degmatch      = round(mean(draws >= obs), 4),
    modules_observed     = obs_mod,
    modules_exp_random   = round(mean(nullmod), 1),
    p_modules_fewer      = round(mean(nullmod <= obs_mod), 4),
    modules_with_2plus_loci = nrow(cross),
    largest_cross_module    = if (nrow(cross)) cross[which.max(n_genes), module_name] else NA_character_,
    largest_cross_module_size = if (nrow(cross)) cross[which.max(n_genes), module_size] else NA_integer_)

  cat(sprintf("[%s] %d loci / %d At-anchors -> %d genes in network | our MYB call %.0f%% | median deg pct %.1f (MYB bg %.1f)\n",
              sp, uniqueN(inet$locus), uniqueN(inet$at_anchor), nrow(inet),
              100 * mean(inet$our_is_MYB), median(inet$deg_pct), median(myb_bg)))
}

per_copy <- rbindlist(per_copy); bg_all <- rbindlist(bg_all); tests <- rbindlist(tests)
fwrite(per_copy, file.path(outdir, "module20_network_table.tsv"), sep = "\t")
fwrite(tests,    file.path(outdir, "module20_tests.tsv"), sep = "\t")
fwrite(rbindlist(modinfo), file.path(outdir, "module20_module_membership.tsv"), sep = "\t")

summ <- bg_all[, .(n = .N, median_deg_pct = round(median(deg_pct), 1),
                   median_degree = round(median(degree), 1),
                   pct_on_conserved = round(100 * mean(on_cons), 1)),
               by = .(species, group)]
summ[, group := factor(group, levels = c("All nodes","All TFs","All MYB/MYB-rel","Module-20 orthologs"))]
setorder(summ, species, group)
fwrite(summ, file.path(outdir, "module20_vs_background.tsv"), sep = "\t")

cat("\n=== Module-20 orthologs vs backgrounds ===\n"); print(summ, row.names = FALSE)
cat("\n=== tests ===\n"); print(tests, row.names = FALSE)

# ---------------------------------------------------------------------------
# FIGURE
# ---------------------------------------------------------------------------
sp_lab <- function(x) factor(x, levels = c("sugarcane", "purple"))
thm <- theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank(),
        legend.position = "bottom")

ecdf_d <- bg_all[group == "All nodes"][order(species, degree),
                 .(degree, pct = 100 * seq_len(.N) / .N), by = species]
ecdf_d <- ecdf_d[, .SD[seq(1, .N, length.out = min(.N, 4000))], by = species]
pts <- per_copy[in_network == TRUE]

pA <- ggplot() +
  geom_step(data = ecdf_d, aes(degree, pct, colour = sp_lab(species)), linewidth = 0.7) +
  geom_point(data = pts, aes(degree, deg_pct, fill = sp_lab(species)),
             shape = 21, size = 2.6, stroke = 0.5, colour = "white", alpha = 0.9) +
  scale_x_log10() + scale_colour_manual(values = COL, guide = "none") +
  scale_fill_manual(values = COL, name = NULL) +
  labs(title = "A  Where Module-20 orthologs sit in each network",
       subtitle = "degree ECDF (line) with every mapped copy marked (points)",
       x = "degree (log10)", y = "percentile of nodes below") + thm

pB <- ggplot(bg_all, aes(group, deg_pct, fill = sp_lab(species))) +
  geom_boxplot(outlier.size = 0.4, alpha = 0.85, position = position_dodge(0.85)) +
  scale_fill_manual(values = COL, name = NULL) + guides(fill = "none") +
  labs(title = "B  Degree percentile vs nested backgrounds",
       subtitle = "the honest comparison is against MYB/MYB-rel, not all nodes",
       x = NULL, y = "degree percentile") +
  thm + theme(axis.text.x = element_text(angle = 20, hjust = 1))

pC_d <- merge(summ[group == "Module-20 orthologs", .(species, observed = pct_on_conserved)],
              tests[, .(species, exp_pct = 100 * on_cons_exp_degmatch / n_genes_in_network,
                        p = p_cons_degmatch)], by = "species")
pC_m <- melt(pC_d, id.vars = c("species", "p"), variable.name = "what", value.name = "pct")
pC_m[, what := factor(what, labels = c("observed", "degree-matched\nexpectation"))]
pC <- ggplot(pC_m, aes(what, pct, fill = sp_lab(species))) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.0f%%", pct)), position = position_dodge(0.8),
            vjust = -0.35, size = 3.1, colour = "grey25") +
  scale_fill_manual(values = COL, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) + guides(fill = "none") +
  labs(title = "C  On cross-species conserved edges",
       subtitle = sprintf("permutation p: sugarcane %.3f, purple %.3f",
                          pC_d[species == "sugarcane", p], pC_d[species == "purple", p]),
       x = NULL, y = "% of mapped copies") + thm

fig <- pA / (pB | pC) + plot_layout(heights = c(1, 1), guides = "collect") &
       theme(legend.position = "bottom")
ggsave(file.path(outdir, "module20_network_overview.png"), fig,
       width = 30, height = 26, units = "cm", dpi = 300)
ggsave(file.path(outdir, "module20_network_overview.pdf"), fig,
       width = 30, height = 26, units = "cm")

cat(sprintf("\nOutputs in %s\n", outdir))
