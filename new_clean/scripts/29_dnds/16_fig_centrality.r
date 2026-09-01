#!/usr/bin/env Rscript
# ============================================================================
# 16_fig_centrality.r -- local clustering as a second description of position
#
# Panel A carries the whole argument's precondition: if clustering is just
# degree upside down, nothing after it is independent evidence. It is drawn
# first and largest for that reason, not as decoration.
# ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(scico); library(scales); library(svglite)
})
env <- function(k, d = NA) { v <- Sys.getenv(k); if (nzchar(v)) v else d }
OUTDIR <- env("OUTDIR"); PREFIX <- file.path(OUTDIR, "fig_centrality")

theme_f <- theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(face = "bold", size = 7.5),
        legend.key.size = unit(3.5, "mm"),
        legend.title = element_text(face = "bold", size = 7.5),
        plot.margin = margin(2, 4, 2, 2),
        plot.subtitle = element_text(size = 6.5, colour = "grey35"),
        plot.tag = element_text(face = "bold", size = 10))
PAL2 <- scico(3, palette = "batlow")[1:2]

d   <- fread(file.path(OUTDIR, "dnds_gene_table.tsv"))
cen <- fread(file.path(OUTDIR, "centrality_sugarcane.tsv"))
setnames(cen, "gene", "sugarcane_gene")
d <- merge(d, cen[, .(sugarcane_gene, clustering)], by = "sugarcane_gene", all.x = TRUE)
a <- d[in_network == TRUE & ds_ok_sc == TRUE & omega_outlier == FALSE &
         !is.na(sc_frac_conserved) & !is.na(sc_mean_log_tpm) &
         sc_sb_omega > 0 & !is.na(clustering)]
tests <- fread(file.path(OUTDIR, "centrality_tests.tsv"))
rho_cd <- tests[statistic == "spearman_clustering_vs_degree", value][1]

# --- A: is clustering just degree upside down? ------------------------------
ck <- fread(file.path(OUTDIR, "centrality_ck_curve.tsv"))
pA <- ggplot(a, aes(sc_degree, clustering)) +
  geom_hex(bins = 45) +
  scale_fill_scico(palette = "devon", direction = -1, trans = "log10", name = "genes") +
  scale_x_log10(labels = label_number(accuracy = 1)) +
  geom_line(data = ck, aes(median_degree, median_clustering), inherit.aes = FALSE,
            colour = "grey15", linewidth = 0.6) +
  geom_point(data = ck, aes(median_degree, median_clustering), inherit.aes = FALSE,
             colour = "grey15", size = 1.2) +
  labs(x = "degree (log)", y = "local clustering coefficient",
       title = "C(k) rises here, so clustering is not degree relabelled",
       subtitle = sprintf("black = median per degree decile; spearman = %+.3f (redundancy screen passes)",
                          rho_cd)) +
  theme_f

# --- B: clustering against constraint ---------------------------------------
pB <- ggplot(a, aes(clustering, sc_sb_omega)) +
  geom_hex(bins = 45) +
  scale_fill_scico(palette = "devon", direction = -1, trans = "log10", name = "genes") +
  scale_y_log10(labels = label_number(accuracy = 0.01)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), colour = "grey15",
              linewidth = 0.5, fill = "grey70") +
  labs(x = "local clustering coefficient", y = expression(omega~"(sugarcane vs sorghum)"),
       title = "Clustering against constraint",
       subtitle = "marginal view -- confounded with degree, see A") +
  theme_f

# --- C: the honest summary --------------------------------------------------
# All three network terms on one axis, in the only currency that compares them.
cmp <- tests[test == "4 partial" & statistic %in% c("clustering", "ldeg", "sc_frac_conserved")]
cmp[, term := factor(fcase(statistic == "clustering", "clustering",
                           statistic == "ldeg", "degree",
                           default = "edge conservation"),
                     levels = c("edge conservation", "degree", "clustering"))]
cmp[, rd := factor(fifelse(readout == "omega", "omega", "constraint score"),
                   levels = c("omega", "constraint score"))]
pC <- ggplot(cmp, aes(term, dR2, fill = rd)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62) +
  scale_fill_manual(values = setNames(PAL2, c("omega", "constraint score")), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = expression(Delta*R^2~"when the term is dropped"),
       title = "How much each network term explains",
       subtitle = "the only currency in which the three are comparable") +
  theme_f + theme(legend.position = "bottom", axis.text.x = element_text(size = 6.4))

# --- D: is the sign stable once degree is fixed? ----------------------------
st <- rbindlist(lapply(c("omega", "constraint"), function(r) {
  f <- file.path(OUTDIR, sprintf("centrality_by_degree_decile_%s.tsv", r))
  if (!file.exists(f)) return(NULL)
  x <- fread(f); x[, readout := fifelse(r == "omega", "omega", "constraint score")]; x
}), fill = TRUE)
pD <- ggplot(st[!is.na(rho)], aes(factor(deg_decile), rho, colour = readout)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.35) +
  geom_point(size = 1.8, position = position_dodge(width = 0.4)) +
  scale_colour_manual(values = setNames(PAL2, c("omega", "constraint score")), name = NULL) +
  labs(x = "degree decile (1 = least connected)",
       y = expression(rho~"(clustering, constraint)"),
       title = "Within degree deciles",
       subtitle = "degree held roughly fixed; consistency of sign is the evidence") +
  theme_f + theme(legend.position = "bottom")

# --- E: do the two readouts agree? ------------------------------------------
pE <- ggplot(a[!is.na(constraint_score)], aes(log(sc_sb_omega), constraint_score)) +
  geom_hex(bins = 45) +
  scale_fill_scico(palette = "devon", direction = -1, trans = "log10", name = "genes") +
  labs(x = expression(log~omega), y = "constraint score (residual log dN)",
       title = "The two readouts",
       subtitle = "both oriented low = constrained; they agree where omega is undistorted") +
  theme_f

fig <- (pA | pB) / (pC | pD) / (pE | plot_spacer()) +
  plot_annotation(tag_levels = "A") + plot_layout(heights = c(1, 1, 0.85))
W <- 22; H <- 22
ggsave(paste0(PREFIX, ".png"), fig, width = W, height = H, units = "cm", dpi = 300)
ggsave(paste0(PREFIX, ".pdf"), fig, width = W, height = H, units = "cm")
ggsave(paste0(PREFIX, ".svg"), fig, width = W, height = H, units = "cm")
cat(sprintf("== wrote %s.{png,pdf,svg}  (%d genes)\n", PREFIX, nrow(a)))
