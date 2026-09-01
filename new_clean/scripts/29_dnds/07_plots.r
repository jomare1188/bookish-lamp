#!/usr/bin/env Rscript
# ============================================================================
# 07_plots.r -- the dN/dS figure
#
# Panel choices follow from what each piece of data has to say:
#   A  omega vs conservation over ~25k genes  -> hexbin (a scatter at this n is
#      a solid blob; density is the honest encoding), single-hue sequential
#   B  omega by conservation bin              -> median + bootstrap CI, because
#      the claim is about a shift in central tendency, not about spread
#   C  rho within each degree decile          -> dot + zero reference line; the
#      question is only "is the sign consistent", so that is all it shows
#   D  the Saccharum pair, counts summed      -> point-range; per-gene omega
#      does not exist here, so no distribution can be drawn
#   E  hubs, the positive control             -> two boxes, thin marks
#
# House style throughout: scico palettes (perceptually uniform and CVD-safe by
# construction), theme_bw(8), PNG + PDF + SVG at 22 cm, exactly as figures
# 20-28 in this pipeline are built.
# ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(scico); library(scales); library(svglite)
})

env <- function(k, d = NA) { v <- Sys.getenv(k); if (nzchar(v)) v else d }
OUTDIR <- env("OUTDIR")
PREFIX <- file.path(OUTDIR, "fig_dnds")

theme_f <- theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(face = "bold", size = 7.5),
        legend.key.size = unit(3.5, "mm"),
        legend.title = element_text(face = "bold", size = 7.5),
        plot.margin = margin(2, 4, 2, 2),
        plot.tag = element_text(face = "bold", size = 10))

d <- fread(file.path(OUTDIR, "dnds_gene_table.tsv"))
DS_MIN <- as.numeric(env("DNDS_DS_MIN")); DS_MAX <- as.numeric(env("DNDS_DS_MAX"))
OM <- as.numeric(env("DNDS_OMEGA_FLAG"))
a <- d[in_network == TRUE & ds_ok_sc == TRUE & omega_outlier == FALSE &
         !is.na(sc_frac_conserved) & sc_sb_omega > 0]

PAL2 <- scico(3, palette = "batlow")[1:2]

# Fixed, interpretable conservation bins, shared by panels B and D. Quantile
# cut-points collapse here because a large share of genes sit at exactly zero
# conserved edges. Declared once so the two panels cannot drift apart, and the
# levels are re-imposed after fread, which returns them as plain character and
# would otherwise sort "0" to the END alphabetically.
CONS_BREAKS <- c(-Inf, 0, 0.05, 0.1, 0.2, Inf)
CONS_LABS   <- c("0", "(0,0.05]", "(0.05,0.1]", "(0.1,0.2]", ">0.2")

# --- A: the relationship, at its real density -------------------------------
pA <- ggplot(a, aes(sc_frac_conserved, sc_sb_omega)) +
  geom_hex(bins = 45) +
  scale_fill_scico(palette = "devon", direction = -1, trans = "log10",
                   name = "genes", labels = label_number(accuracy = 1)) +
  scale_y_log10(labels = label_number(accuracy = 0.01)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), colour = "grey15",
              linewidth = 0.5, se = TRUE, fill = "grey70") +
  labs(x = "fraction of edges conserved (sugarcane → purple)",
       y = expression(omega~"(sugarcane vs sorghum)"),
       title = "Constraint against neighbourhood conservation") +
  theme_f

# --- B: the effect size, which is what should be read off the figure --------
bq <- a[, .(bin = cut(sc_frac_conserved, breaks = CONS_BREAKS, labels = CONS_LABS),
            omega = sc_sb_omega)]
bs <- bq[, {
  set.seed(1)
  bm <- replicate(2000, median(sample(omega, .N, replace = TRUE)))
  .(n = .N, med = median(omega), lo = quantile(bm, .025), hi = quantile(bm, .975))
}, by = bin][!is.na(bin)][order(bin)]
pB <- ggplot(bs, aes(bin, med)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.12, linewidth = 0.4,
                colour = PAL2[1]) +
  geom_point(size = 1.9, colour = PAL2[1]) +
  labs(x = "fraction of edges conserved", y = expression("median"~omega),
       title = "Long-term constraint by conservation",
       subtitle = "bars: 95% bootstrap CI of the median; note the y-range") +
  theme_f + theme(axis.text.x = element_text(size = 6.2),
                  plot.subtitle = element_text(size = 6.5, colour = "grey35"))

# --- C: is the sign consistent once degree is held roughly fixed? -----------
st <- fread(file.path(OUTDIR, "dnds_by_degree_decile.tsv"))
pC <- ggplot(st, aes(factor(deg_decile), rho)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.35) +
  geom_point(size = 1.9, colour = PAL2[2]) +
  labs(x = "degree decile (1 = least connected)",
       y = expression(rho~"("*omega*", conservation)"),
       title = "Within degree deciles",
       subtitle = "degree and conservation covary; this holds degree roughly fixed") +
  theme_f + theme(plot.subtitle = element_text(size = 6.5, colour = "grey35"))

# --- D: the Saccharum pair, the only honest form for it ---------------------
# The pooled trend across conservation bins looks like a clean decline -- and it
# is an EXPRESSION trend. Median log2 TPM runs 0.45 -> 2.16 across these bins, so
# they are not otherwise comparable. Drawing the pooled line alone would be the
# misleading version of this panel; it is drawn in grey BEHIND the same ratio
# re-formed within expression tertiles, which is where the decline disappears.
pb <- fread(file.path(OUTDIR, "dnds_saccharum_pair_binned.tsv"))
pb <- pb[!is.na(omega)]
pb[, cons_bin := factor(cons_bin, levels = CONS_LABS)]
ps <- fread(file.path(OUTDIR, "dnds_saccharum_pair_stratified.tsv"))
ps <- ps[stratum == "expression" & !is.na(omega)]
ps[, cons_bin := factor(cons_bin, levels = CONS_LABS)]
ps[, level := factor(level, levels = c("low", "mid", "high"))]
PAL3 <- scico(5, palette = "batlow")[c(1, 3, 4)]

pD <- ggplot(ps, aes(cons_bin, omega, colour = level, group = level)) +
  geom_line(data = pb, aes(cons_bin, omega), inherit.aes = FALSE, group = 1,
            colour = "grey60", linewidth = 1.1) +
  geom_point(data = pb, aes(cons_bin, omega), inherit.aes = FALSE,
             colour = "grey60", size = 1.6) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 1.7) +
  annotate("text", x = 1.15, y = pb[1, omega], label = "pooled", hjust = 0,
           vjust = -0.9, size = 2.3, colour = "grey40") +
  scale_colour_manual(values = setNames(PAL3, c("low", "mid", "high")),
                      name = "expression\ntertile") +
  labs(x = "fraction of edges conserved",
       y = expression(omega == (Sigma*N[d]/Sigma*N)/(Sigma*S[d]/Sigma*S)),
       title = "R570 vs LA purple: the pooled decline is expression",
       subtitle = "pooled (grey) falls; within expression tertiles it does not") +
  theme_f + theme(axis.text.x = element_text(size = 6.2),
                  plot.subtitle = element_text(size = 6.5, colour = "grey35"))

# --- E: the positive control ------------------------------------------------
a[, hub := factor(fifelse(sc_is_hub, "hub", "non-hub"), levels = c("non-hub", "hub"))]
pE <- ggplot(a, aes(hub, sc_sb_omega, fill = hub)) +
  geom_boxplot(outlier.size = 0.25, outlier.alpha = 0.25, linewidth = 0.3,
               width = 0.55) +
  scale_fill_manual(values = setNames(PAL2, c("non-hub", "hub")), guide = "none") +
  scale_y_log10(labels = label_number(accuracy = 0.01)) +
  labs(x = NULL, y = expression(omega),
       title = "Positive control: hubs",
       subtitle = "hubs are expected to be more constrained") +
  theme_f + theme(plot.subtitle = element_text(size = 6.5, colour = "grey35"))

fig <- (pA | pB) / (pC | pD) / (pE | plot_spacer()) +
  plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(1, 1, 0.85))

W <- 22; H <- 22
ggsave(paste0(PREFIX, ".png"), fig, width = W, height = H, units = "cm", dpi = 300)
ggsave(paste0(PREFIX, ".pdf"), fig, width = W, height = H, units = "cm")
ggsave(paste0(PREFIX, ".svg"), fig, width = W, height = H, units = "cm")
cat(sprintf("== wrote %s.{png,pdf,svg}  (%d genes plotted)\n", PREFIX, nrow(a)))
