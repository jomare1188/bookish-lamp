# ============================================================================
# 05_tf_characterization.r — network characteristics of the TFs (H1, general view)
#
# For every TF that is a node of each co-expression network, assemble:
#   family, type, degree, strength, transitivity (centrality proxies),
#   nitrogen (treatment) correlation, whether it sits on a cross-species
#   conserved edge, and its ortholog-conservation status.
# Then summarise per family and draw a MYB / MYB-related comparative figure.
#
# NB: betweenness is NOT computed — exact betweenness on 75M–681M-edge networks
#     is infeasible; degree / strength / transitivity are the centrality proxies.
#
# RUN:  conda activate r_env ; Rscript 05_tf_characterization.r
# ============================================================================

suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })

base   <- "/dados04/jorge/comparative_saccharum"
outdir <- file.path(base, "GET_TFS/new/results/characterization")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

MYB_FAMS <- c("MYB", "MYB-related")
R_THR    <- 0.6      # |Pearson| for "nitrogen-correlated"
PADJ_THR <- 0.05
COL      <- c(sugarcane = "#1b7837", purple = "#762a83")

cfg <- list(
  sugarcane = list(
    tf        = file.path(base, "GET_TFS/new/results/sugarcane/TF_in_network.tsv"),
    nodes     = file.path(base, "files/sugarcane/network_sugarcane_node_metrics.tsv"),
    node_strip= "\\.v[0-9.]+$",   # node ids carry a .v2.1 suffix; TF/trait ids do not
    trait     = file.path(base, "files/sugarcane/gene_trait_correlations_sugarcane.tsv"),
    conserved = file.path(base, "files/network_conservation/conserved_genes_sugarcane_FULL.txt"),
    status    = file.path(base, "files/network_conservation/sugarcane_correlated_conservation_status.tsv")),
  purple = list(
    tf        = file.path(base, "GET_TFS/new/results/purple/TF_in_network.tsv"),
    nodes     = file.path(base, "files/purple/new/network_purple_node_metrics.tsv"),
    node_strip= NULL,             # node ids already match TF/trait ids
    trait     = file.path(base, "files/purple/new/gene_trait_correlations_purple.tsv"),
    conserved = file.path(base, "files/network_conservation/conserved_genes_purple_FULL.txt"),
    status    = file.path(base, "files/network_conservation/purple_correlated_conservation_status.tsv"))
)

# ---------------------------------------------------------------------------
# per-species assembly
# ---------------------------------------------------------------------------
node_bg  <- list()   # all-node degree background, for the distribution panel
per_gene <- list()   # one row per TF gene
tf_long  <- list()   # one row per TF gene x family (for family-level counts)

for (sp in names(cfg)) {
  c <- cfg[[sp]]
  nodes <- fread(c$nodes)                                  # gene degree strength transitivity
  if (!is.null(c$node_strip)) nodes[, gene := sub(c$node_strip, "", gene)]
  tf    <- fread(c$tf); setnames(tf, c("gene","Family","Type"))
  trait <- fread(c$trait)
  trt   <- trait[trait == "treatment",
                 .(gene, treatment_r = suppressWarnings(as.numeric(pearson)),
                        treatment_padj = suppressWarnings(as.numeric(padj)))]
  cons   <- fread(c$conserved, header = FALSE)$V1
  status <- fread(c$status); setnames(status, 1:5,
                 c("gene","status_r","n_orthologs","n_corr","conservation_status"))

  # collapse isoform-driven multi-family calls to one row per gene
  g <- tf[, .(family  = paste(sort(unique(Family)), collapse = ";"),
              type    = paste(sort(unique(Type)),   collapse = ";"),
              is_MYB  = any(Family %in% MYB_FAMS)), by = gene]

  g <- merge(g, nodes, by = "gene", all.x = TRUE)
  g <- merge(g, trt,   by = "gene", all.x = TRUE)
  g[, on_conserved_edge := gene %in% cons]
  g[, N_correlated := !is.na(treatment_r) & abs(treatment_r) >= R_THR &
                      !is.na(treatment_padj) & treatment_padj <= PADJ_THR]
  g <- merge(g, status[, .(gene, conservation_status)], by = "gene", all.x = TRUE)
  g[, species := sp]
  setcolorder(g, c("species","gene","family","type","is_MYB","degree","strength",
                   "transitivity","treatment_r","treatment_padj","N_correlated",
                   "on_conserved_edge","conservation_status"))
  per_gene[[sp]] <- g

  # exploded gene x family, carrying the per-gene flags (for counts by family)
  l <- merge(unique(tf[, .(gene, Family, Type)]),
             g[, .(gene, degree, strength, transitivity, treatment_r,
                   treatment_padj, N_correlated, on_conserved_edge)],
             by = "gene", all.x = TRUE)
  l[, species := sp]
  tf_long[[sp]] <- l

  nb <- nodes[, .(degree, strength, transitivity)]; nb[, species := sp]
  node_bg[[sp]] <- nb

  cat(sprintf("[%s] TF genes=%d  MYB/MYB-related=%d  on-conserved-edge=%d  N-correlated=%d\n",
              sp, nrow(g), sum(g$is_MYB), sum(g$on_conserved_edge), sum(g$N_correlated)))
}

per_gene <- rbindlist(per_gene)
tf_long  <- rbindlist(tf_long)
node_bg  <- rbindlist(node_bg)

# ---------------------------------------------------------------------------
# TABLES
# ---------------------------------------------------------------------------
fwrite(per_gene[species=="sugarcane"], file.path(outdir,"TF_network_characteristics_sugarcane.tsv"), sep="\t")
fwrite(per_gene[species=="purple"],    file.path(outdir,"TF_network_characteristics_purple.tsv"),    sep="\t")
fwrite(per_gene,                        file.path(outdir,"TF_network_characteristics_ALL.tsv"),        sep="\t")

fam_summary <- tf_long[, .(
  n_genes            = uniqueN(gene),
  median_degree      = round(median(degree, na.rm=TRUE), 1),
  mean_degree        = round(mean(degree,   na.rm=TRUE), 1),
  median_strength    = round(median(strength, na.rm=TRUE), 2),
  median_transitivity= round(median(transitivity, na.rm=TRUE), 3),
  pct_on_conserved   = round(100*mean(on_conserved_edge), 1),
  pct_N_correlated   = round(100*mean(N_correlated), 1)
), by = .(Family, Type, species)][order(species, -n_genes)]
fwrite(fam_summary, file.path(outdir,"TF_family_summary.tsv"), sep="\t")

# ---------------------------------------------------------------------------
# FIGURE — MYB / MYB-related, comparative
# ---------------------------------------------------------------------------
sp_lab <- function(x) factor(x, levels=c("sugarcane","purple"))

## Panel A — family landscape (top families, MYB highlighted) --------------
famN <- tf_long[, .(n = uniqueN(gene)), by=.(species, Family)]
top  <- famN[, .(tot=sum(n)), by=Family][order(-tot)][1:15, Family]
pA_d <- famN[Family %in% top]
pA_d[, Family := factor(Family, levels=rev(top))]
pA_d[, hl := ifelse(Family %in% MYB_FAMS, "MYB / MYB-related", "other TF family")]
pA <- ggplot(pA_d, aes(Family, n, fill=sp_lab(species), alpha=hl)) +
  geom_col(position=position_dodge(width=0.8), width=0.75) +
  coord_flip() +
  scale_fill_manual(values=COL, name=NULL) +
  scale_alpha_manual(values=c("MYB / MYB-related"=1, "other TF family"=0.45), name=NULL) +
  labs(title="A  TF family landscape in each network",
       subtitle="TF genes that are network nodes; MYB / MYB-related highlighted",
       x=NULL, y="TF genes in network") +
  theme_bw(base_size=12) + theme(plot.title=element_text(face="bold"))

## Panel B — are MYBs hubs? degree vs background ---------------------------
mk <- function(dt, grp) dt[, .(degree, species)][, group := grp]
pB_d <- rbindlist(list(
  mk(node_bg,                    "All network nodes"),
  mk(per_gene,                   "All TFs"),
  mk(per_gene[is_MYB==TRUE],     "MYB / MYB-related")))
pB_d[, group := factor(group, levels=c("All network nodes","All TFs","MYB / MYB-related"))]
pB <- ggplot(pB_d, aes(group, degree, fill=sp_lab(species))) +
  geom_violin(scale="width", alpha=0.55, colour=NA, position=position_dodge(0.85)) +
  geom_boxplot(width=0.14, outlier.shape=NA, position=position_dodge(0.85), alpha=0.9) +
  scale_y_log10() + scale_fill_manual(values=COL, name=NULL) +
  labs(title="B  Are MYBs central? degree vs background",
       subtitle="log10 degree; MYBs vs all TFs vs all nodes", x=NULL, y="degree (log10)") +
  theme_bw(base_size=12) +
  theme(plot.title=element_text(face="bold"), axis.text.x=element_text(angle=20, hjust=1))

## Panel C — conservation & nitrogen response of MYBs ----------------------
grp_stat <- function(d, g) d[, .(species,
                                 pct_on_conserved = 100*mean(on_conserved_edge),
                                 pct_N_correlated = 100*mean(N_correlated)),
                             by=species][, group := g][]
pC_d <- rbindlist(list(
  grp_stat(per_gene,                                   "All TFs"),
  grp_stat(tf_long[Family=="MYB"],                     "MYB"),
  grp_stat(tf_long[Family=="MYB-related"],             "MYB-related")))
pC_d <- unique(pC_d)
pC_m <- melt(pC_d, id.vars=c("species","group"),
             measure.vars=c("pct_on_conserved","pct_N_correlated"),
             variable.name="metric", value.name="pct")
pC_m[, metric := factor(metric, labels=c("on conserved edge","N-correlated (|r|≥0.6)"))]
pC_m[, group := factor(group, levels=c("All TFs","MYB","MYB-related"))]
pC <- ggplot(pC_m, aes(group, pct, fill=sp_lab(species))) +
  geom_col(position=position_dodge(0.8), width=0.72) +
  facet_wrap(~metric, scales="free_y") +
  scale_fill_manual(values=COL, name=NULL) +
  labs(title="C  Conservation & nitrogen response",
       subtitle="% of genes in each group", x=NULL, y="percent of genes") +
  theme_bw(base_size=12) +
  theme(plot.title=element_text(face="bold"), axis.text.x=element_text(angle=20, hjust=1))

fig <- pA / (pB | pC) + plot_layout(heights=c(1.05,1), guides="collect") &
  theme(legend.position="bottom")

ggsave(file.path(outdir,"MYB_TF_network_overview.png"), fig, width=30, height=26, units="cm", dpi=300)
ggsave(file.path(outdir,"MYB_TF_network_overview.pdf"), fig, width=30, height=26, units="cm")

# ---------------------------------------------------------------------------
# console summary of the MYB comparison
# ---------------------------------------------------------------------------
cat("\n=== MYB / MYB-related summary ===\n")
print(fam_summary[Family %in% MYB_FAMS][order(Family, species)], row.names=FALSE)
cat(sprintf("\nOutputs in %s\n", outdir))
