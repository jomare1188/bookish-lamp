# ============================================================================
# 07_network_readout.r — what the MYB61 copies actually do in each network (H1)
#
# For every confirmed MYB61 orthologue (steps 04 + 06), report its position in
# its own co-expression network and place it against three backgrounds:
# all nodes, all TFs, and all MYB-family TFs. The question H1 asks is whether
# the conserved core is the MYB regulatory layer — so the readout is:
#   is MYB61 a hub?  is it nitrogen-responsive?  is it on conserved edges?
#
# Betweenness is NOT computed: exact betweenness on 75M-681M-edge networks is
# infeasible. Degree / strength / transitivity are the centrality proxies, and
# degree is additionally expressed as a PERCENTILE within each network so the
# two networks (very different edge counts) can be compared at all.
#
# Palette: ColorBrewer PRGn endpoints, as used by 05_tf_characterization.r.
# (The dataviz validator could not be executed here — no Node runtime on this
# box — so the palette is justified by provenance: PRGn is a published
# colorblind-safe diverging scheme and these are its two poles.)
#
# RUN:  /home/genomics/miniconda3/envs/r_env/bin/Rscript 07_network_readout.r
# ============================================================================

suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })

base   <- "/dados04/jorge/comparative_saccharum"
# cached: steps 01-06 output (candidates, phylogeny, work/)
myb61  <- file.path(base, "GET_TFS/new/results/myb61")
outdir <- "/dados04/jorge/comparative_saccharum/new_clean/results/readouts/myb61"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
R_THR    <- 0.6
PADJ_THR <- 0.05
COL      <- c(sugarcane = "#1b7837", purple = "#762a83")

cfg <- list(
  sugarcane = list(
    cand      = file.path(myb61, "sugarcane_MYB61_candidates.tsv"),
    nodes     = file.path(base, "new_clean/results/sugarcane/network_sugarcane_node_metrics.tsv"),
    node_strip= "\\.v[0-9.]+$",
    trait     = file.path(base, "new_clean/results/sugarcane/gene_trait_correlations_sugarcane.tsv"),
    conserved = file.path(base, "new_clean/results/conservation/conserved_genes_sugarcane_FULL.txt"),
    status    = file.path(base, "new_clean/results/conservation/sugarcane_correlated_conservation_status.tsv"),
    tf        = "/dados04/jorge/comparative_saccharum/new_clean/results/readouts/get_tfs/sugarcane/TF_in_network.tsv"),
  purple = list(
    cand      = file.path(myb61, "purple_MYB61_candidates.tsv"),
    nodes     = file.path(base, "new_clean/results/purple/network_purple_node_metrics.tsv"),
    node_strip= NULL,
    trait     = file.path(base, "new_clean/results/purple/gene_trait_correlations_purple.tsv"),
    conserved = file.path(base, "new_clean/results/conservation/conserved_genes_purple_FULL.txt"),
    status    = file.path(base, "new_clean/results/conservation/purple_correlated_conservation_status.tsv"),
    tf        = "/dados04/jorge/comparative_saccharum/new_clean/results/readouts/get_tfs/purple/TF_in_network.tsv")
)

# copies rescued by the tree (step 06) count as MYB61 even though the strict
# domain rules downgraded them — the gene model is broken, the orthology is not
rescued <- character(0)
rf <- file.path(myb61, "phylogeny", "MYB61_recovered_by_tree.ids")
if (file.exists(rf)) rescued <- sub("\\.[0-9]+\\.p[0-9]*$", "", readLines(rf))

per_gene <- list(); bg <- list()

for (sp in names(cfg)) {
  cc <- cfg[[sp]]

  cand <- fread(cc$cand)
  cand[, gene_id := gene]
  cand[, call := fifelse(verdict == "MYB61_ortholog", "strict (RBH + MYB rule)",
                  fifelse(gene_id %in% rescued, "rescued by tree", NA_character_))]
  m <- cand[!is.na(call), .(gene = gene_id, clade = best_anchor, bitscore, qcov, call)]
  m[, clade := fifelse(grepl("3003G136600", clade), "chr3-type co-ortholog",
                                                    "chr9-type co-ortholog")]

  nodes <- fread(cc$nodes)
  if (!is.null(cc$node_strip)) nodes[, gene := sub(cc$node_strip, "", gene)]
  nodes[, deg_pct := 100 * frank(degree, ties.method = "average") / .N]

  trait <- fread(cc$trait)
  trt <- trait[trait == "treatment",
               .(gene, treatment_r = suppressWarnings(as.numeric(pearson)),
                      treatment_padj = suppressWarnings(as.numeric(padj)))]

  cons   <- fread(cc$conserved, header = FALSE)$V1
  status <- fread(cc$status); setnames(status, 1:5,
             c("gene","status_r","n_orthologs","n_corr","conservation_status"))
  tfin <- fread(cc$tf); setnames(tfin, c("gene","Family","Type"))

  m <- merge(m, nodes, by = "gene", all.x = TRUE)
  m <- merge(m, trt,   by = "gene", all.x = TRUE)
  m <- merge(m, status[, .(gene, conservation_status)], by = "gene", all.x = TRUE)
  m[, in_network       := !is.na(degree)]
  m[, on_conserved_edge:= gene %in% cons]
  m[, N_correlated     := !is.na(treatment_r) & abs(treatment_r) >= R_THR &
                          !is.na(treatment_padj) & treatment_padj <= PADJ_THR]
  m[, species := sp]
  setcolorder(m, c("species","gene","call","clade","in_network","degree","deg_pct",
                   "strength","transitivity","treatment_r","treatment_padj",
                   "N_correlated","on_conserved_edge","conservation_status",
                   "bitscore","qcov"))
  per_gene[[sp]] <- m

  # backgrounds for the comparison panels
  mk <- function(genes, lab) {
    d <- nodes[gene %in% genes]
    if (!nrow(d)) return(NULL)
    data.table(species = sp, group = lab, gene = d$gene,
               degree = d$degree, deg_pct = d$deg_pct,
               on_conserved_edge = d$gene %in% cons)
  }
  bg[[sp]] <- rbindlist(list(
    mk(nodes$gene,                          "All nodes"),
    mk(tfin$gene,                           "All TFs"),
    mk(tfin[Family %in% c("MYB","MYB-related"), gene], "All MYB/MYB-rel"),
    mk(m[in_network == TRUE, gene],         "MYB61")), use.names = TRUE)

  cat(sprintf("[%s] MYB61 copies=%d  in network=%d  median degree pct=%.1f  N-corr=%d  on conserved edge=%d\n",
      sp, nrow(m), sum(m$in_network), median(m$deg_pct, na.rm = TRUE),
      sum(m$N_correlated), sum(m$on_conserved_edge)))
}

per_gene <- rbindlist(per_gene)
bg       <- rbindlist(bg)

fwrite(per_gene, file.path(outdir, "MYB61_network_table.tsv"), sep = "\t")

# ---------------------------------------------------------------------------
# summary table: MYB61 vs the three backgrounds
# ---------------------------------------------------------------------------
summ <- bg[, .(n = .N,
               median_degree      = round(median(degree), 1),
               median_degree_pct  = round(median(deg_pct), 1),
               pct_on_conserved   = round(100 * mean(on_conserved_edge), 1)),
           by = .(species, group)]
summ[, group := factor(group, levels = c("All nodes","All TFs","All MYB/MYB-rel","MYB61"))]
setorder(summ, species, group)
fwrite(summ, file.path(outdir, "MYB61_vs_background_summary.tsv"), sep = "\t")
cat("\n=== MYB61 vs background ===\n"); print(summ, row.names = FALSE)

# --- is MYB61 enriched on conserved edges relative to its own network? ------
# one-sided exact binomial against that network's all-node conserved-edge rate
enr <- rbindlist(lapply(names(cfg), function(sp) {
  m   <- per_gene[species == sp & in_network == TRUE]
  p0  <- summ[species == sp & group == "All nodes", pct_on_conserved] / 100
  k   <- sum(m$on_conserved_edge); n <- nrow(m)
  bt  <- binom.test(k, n, p = p0, alternative = "greater")
  data.table(species = sp, n_in_network = n, n_on_conserved = k,
             pct_MYB61 = round(100 * k / n, 1), pct_background = round(100 * p0, 1),
             fold = round((k / n) / p0, 2), binom_p = signif(bt$p.value, 4),
             max_abs_r_with_N = round(max(abs(m$treatment_r), na.rm = TRUE), 3),
             n_N_correlated = sum(m$N_correlated))
}))

# The "on conserved edge" flag is BINARY and STRONGLY DEGREE-DEPENDENT (a gene
# is flagged if >=1 of its edges is conserved, so a hub has far more chances).
# The all-node null above is therefore not sufficient on its own: a gene set with
# an unusual degree profile could look enriched for that reason alone. Repeat the
# test against a DEGREE-MATCHED null — resample, for each MYB61 copy, a random
# node of near-identical degree (+/-250 rank positions) and count how often the
# resampled set reaches the observed number.
set.seed(1)
N_PERM <- 20000
deg_match <- rbindlist(lapply(names(cfg), function(sp) {
  cc <- cfg[[sp]]
  nodes <- fread(cc$nodes)
  if (!is.null(cc$node_strip)) nodes[, gene := sub(cc$node_strip, "", gene)]
  cons <- fread(cc$conserved, header = FALSE)$V1
  nodes[, oc := gene %in% cons]
  pool <- nodes[order(degree), .(degree, oc)]
  m    <- per_gene[species == sp & in_network == TRUE]
  obs  <- sum(m$on_conserved_edge)
  idx  <- sapply(m$degree, function(d) which.min(abs(pool$degree - d)))
  draws <- replicate(N_PERM, sum(sapply(idx, function(i) {
             lo <- max(1, i - 250); hi <- min(nrow(pool), i + 250)
             pool$oc[sample(lo:hi, 1)] })))
  data.table(species = sp, n = nrow(m), observed = obs,
             expected_degree_matched = round(mean(draws), 2),
             fold = round(obs / mean(draws), 2),
             p_degree_matched = round(mean(draws >= obs), 4))
}))
enr <- merge(enr, deg_match[, .(species, expected_degree_matched, fold_dm = fold,
                                p_degree_matched)], by = "species")
fwrite(enr, file.path(outdir, "MYB61_conserved_edge_enrichment.tsv"), sep = "\t")
cat("\n=== conserved-edge enrichment (exact binomial vs own network) ===\n")
print(enr, row.names = FALSE)

# ---------------------------------------------------------------------------
# FIGURE
# ---------------------------------------------------------------------------
sp_lab <- function(x) factor(x, levels = c("sugarcane","purple"))
thm <- theme_bw(base_size = 12) +
       theme(plot.title = element_text(face = "bold"),
             panel.grid.minor = element_blank(),
             legend.position = "bottom")

## A — is MYB61 a hub? each copy on the network's own degree ECDF ------------
ecdf_d <- bg[group == "All nodes", .(degree = sort(unique(degree))), by = species]
ecdf_d <- bg[group == "All nodes"][order(species, degree),
             .(degree, pct = 100 * seq_len(.N) / .N), by = species]
ecdf_d <- ecdf_d[, .SD[seq(1, .N, length.out = min(.N, 4000))], by = species]
pts    <- per_gene[in_network == TRUE]

pA <- ggplot() +
  geom_step(data = ecdf_d, aes(degree, pct, colour = sp_lab(species)), linewidth = 0.7) +
  geom_point(data = pts, aes(degree, deg_pct, fill = sp_lab(species)),
             shape = 21, size = 3.2, stroke = 0.6, colour = "white") +
  scale_x_log10() +
  scale_colour_manual(values = COL, guide = "none") +
  scale_fill_manual(values = COL, name = NULL) +
  labs(title = "A  Is MYB61 a hub?",
       subtitle = "network degree ECDF (line) with every MYB61 copy marked (points)",
       x = "degree (log10)", y = "percentile of nodes below") + thm

## B — nitrogen response of each copy ---------------------------------------
pB_d <- pts[!is.na(treatment_r)][order(treatment_r)]
pB_d[, gene := factor(gene, levels = gene)]
pB <- ggplot(pB_d, aes(treatment_r, gene, fill = sp_lab(species))) +
  annotate("rect", xmin = -R_THR, xmax = R_THR, ymin = -Inf, ymax = Inf,
           fill = "grey85", alpha = 0.5) +
  geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.3) +
  geom_point(shape = 21, size = 3, stroke = 0.5, colour = "white") +
  scale_fill_manual(values = COL, name = NULL) +
  coord_cartesian(xlim = c(-1, 1)) +
  guides(fill = "none") +
  labs(title = "B  Nitrogen response per copy",
       subtitle = sprintf("grey band = |r| < %.1f; %d of %d in-network copies scored",
                          R_THR, nrow(pB_d), nrow(pts)),
       x = "correlation with nitrogen treatment", y = NULL) +
  thm + theme(axis.text.y = element_text(size = 6.5))

## C — conservation, MYB61 vs backgrounds -----------------------------------
pC <- ggplot(summ, aes(group, pct_on_conserved, fill = sp_lab(species))) +
  geom_col(position = position_dodge(0.8), width = 0.72) +
  geom_text(aes(label = sprintf("%.0f%%", pct_on_conserved)),
            position = position_dodge(0.8), vjust = -0.35, size = 3.1, colour = "grey25") +
  scale_fill_manual(values = COL, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  guides(fill = "none") +
  labs(title = "C  On conserved edges",
       subtitle = "MYB61 vs three backgrounds, same network",
       x = NULL, y = "% of genes") +
  thm + theme(axis.text.x = element_text(angle = 20, hjust = 1))

fig <- pA / (pB | pC) + plot_layout(heights = c(1, 1.1), guides = "collect") &
       theme(legend.position = "bottom")

ggsave(file.path(outdir, "MYB61_network_overview.png"), fig,
       width = 30, height = 26, units = "cm", dpi = 300)
ggsave(file.path(outdir, "MYB61_network_overview.pdf"), fig,
       width = 30, height = 26, units = "cm")

cat(sprintf("\nOutputs in %s\n", outdir))
