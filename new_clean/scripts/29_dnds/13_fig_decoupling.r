#!/usr/bin/env Rscript
# ============================================================================
# 13_fig_decoupling.r -- the polyploid figure
#
# The panels are ordered so the figure argues rather than just displays: the
# claim, then the scale it should be read against, then the two controls that
# could destroy it, then the contrasts. Panel D repeats the device 07_plots.r
# panel D uses -- pooled line in grey behind the stratified ones -- because the
# same confounder (expression) is the same threat here.
# ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(scico); library(scales); library(svglite)
})
env <- function(k, d = NA) { v <- Sys.getenv(k); if (nzchar(v)) v else d }
OUTDIR <- env("OUTDIR"); PREFIX <- file.path(OUTDIR, "fig_decoupling")
FLOOR <- as.numeric(env("MIN_FRAC_UNIQUE")); NEAR <- as.numeric(env("PID_NEAR_IDENTICAL"))

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
PAL3 <- scico(5, palette = "batlow")[c(1, 3, 4)]
ID_LABS <- c("<0.95", "0.95-0.98", "0.98-0.99", "0.99-0.999", ">=0.999")

d <- fread(file.path(OUTDIR, "decoupling_pairs_sugarcane.tsv"))
h <- d[pair_class == "homeolog" & is.finite(net_div)]
h[, id_bin := cut(pid_cds, c(-Inf, 0.95, 0.98, 0.99, 0.999, Inf), labels = ID_LABS)]

# --- A: the claim, at its real density --------------------------------------
pA <- ggplot(h, aes(pid_cds, 2^net_div)) +
  geom_hex(bins = 45) +
  scale_fill_scico(palette = "devon", direction = -1, trans = "log10", name = "pairs") +
  scale_y_log10(labels = label_number(accuracy = 1)) +
  geom_vline(xintercept = NEAR, linetype = "22", colour = "grey30", linewidth = 0.4) +
  annotate("text", x = NEAR, y = Inf, label = "near-identical ", hjust = 1, vjust = 1.6,
           size = 2.3, colour = "grey30") +
  labs(x = "CDS identity between the two copies",
       y = "fold difference in network degree",
       title = "Sequence identity barely predicts network position",
       subtitle = expression(paste("homeologous copy pairs, sugarcane; Spearman ", rho, " = 0.048"))) +
  theme_f

# --- B: what does 'far apart' even mean? ------------------------------------
tests <- fread(file.path(OUTDIR, "decoupling_tests.tsv"))
gv <- function(st, stat) tests[study == st & statistic == stat, value][1]
nullcmp <- data.table(
  what = factor(c("homeologous copies", "random pairs,\nexpression-matched"),
                levels = c("homeologous copies", "random pairs,\nexpression-matched")),
  fold = c(gv("sugarcane", "median_fold_real_homeolog_pairs"),
           gv("sugarcane", "median_fold_expression_matched_random_pairs")))
pB <- ggplot(nullcmp, aes(what, fold, fill = what)) +
  geom_col(width = 0.55) +
  scale_fill_manual(values = setNames(c(PAL2[1], "grey65"), levels(nullcmp$what)),
                    guide = "none") +
  geom_text(aes(label = sprintf("%.1fx", fold)), vjust = -0.4, size = 2.6) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = NULL, y = "median fold difference in degree",
       title = "Copies are MORE similar than chance",
       subtitle = "unrelated genes matched on expression decile diverge further") +
  theme_f + theme(axis.text.x = element_text(size = 6.2))

# --- C: do they even land in the same module? -------------------------------
mod <- h[!is.na(diff_module), .(pct = 100 * mean(diff_module), n = .N), by = id_bin][order(id_bin)]
nullmod <- gv("sugarcane", "pct_diff_module_random_pairs")
pC <- ggplot(mod[!is.na(id_bin)], aes(id_bin, pct)) +
  geom_col(fill = PAL2[2], width = 0.6) +
  geom_hline(yintercept = nullmod, colour = "grey35", linetype = "22", linewidth = 0.45) +
  annotate("text", x = 0.6, y = nullmod, label = "random pairs", hjust = 0, vjust = -0.5,
           size = 2.3, colour = "grey35") +
  scale_y_continuous(limits = c(0, max(100, nullmod * 1.05))) +
  labs(x = "CDS identity between copies", y = "% of pairs in DIFFERENT modules",
       title = "Different module -- but less often than chance",
       subtitle = "MCL modules; dashed line is expression-matched random pairs") +
  theme_f + theme(axis.text.x = element_text(size = 6.2))

# --- D: the expression control ----------------------------------------------
st <- fread(file.path(OUTDIR, "decoupling_by_identity_bin_sugarcane.tsv"))
st[, id_bin := factor(id_bin, levels = ID_LABS)]
st[, expr_stratum := factor(expr_stratum, levels = c("low", "mid", "high"))]
pooled <- h[!is.na(id_bin), .(median_fold = 2^median(net_div)), by = id_bin][order(id_bin)]
pD <- ggplot(st[!is.na(id_bin)], aes(id_bin, median_fold, colour = expr_stratum,
                                     group = expr_stratum)) +
  geom_line(data = pooled, aes(id_bin, median_fold), inherit.aes = FALSE, group = 1,
            colour = "grey60", linewidth = 1.1) +
  geom_point(data = pooled, aes(id_bin, median_fold), inherit.aes = FALSE,
             colour = "grey60", size = 1.6) +
  geom_line(linewidth = 0.55) + geom_point(size = 1.7) +
  scale_colour_manual(values = setNames(PAL3, c("low", "mid", "high")),
                      name = "expression\ndifference") +
  labs(x = "CDS identity between copies", y = "median fold difference in degree",
       title = "It is expression, not sequence",
       subtitle = "within each stratum the identity trend is flat; the strata differ 2x") +
  theme_f + theme(axis.text.x = element_text(size = 6.2))

# --- E: the control that can kill the result --------------------------------
mc <- fread(file.path(OUTDIR, "decoupling_mapping_control_sugarcane.tsv"))
mc[, grp := factor(fifelse(distinguishable, "salmon can separate them",
                           "k-mer ambiguous"),
                   levels = c("salmon can separate them", "k-mer ambiguous"))]
mc[, idl := factor(fifelse(near_identical, sprintf("CDS id >= %.2f", NEAR),
                           sprintf("CDS id < %.2f", NEAR)))]
pE <- ggplot(mc, aes(idl, median_fold, fill = grp)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62) +
  scale_fill_manual(values = setNames(c(PAL2[1], "grey65"), levels(mc$grp)), name = NULL) +
  labs(x = NULL, y = "median fold difference in degree",
       title = "It is partly a mapping artefact",
       subtitle = sprintf("copies salmon CAN separate diverge less (floor: %.2f unique 31-mers)",
                          FLOOR)) +
  theme_f + theme(axis.text.x = element_text(size = 6.2), legend.position = "bottom")

# --- F: contrasts ------------------------------------------------------------
cls <- rbindlist(lapply(c("sugarcane", "purple"), function(s) {
  f <- file.path(OUTDIR, sprintf("decoupling_by_class_%s.tsv", s))
  if (!file.exists(f)) return(NULL)
  x <- fread(f); x[, study := s]; x
}), fill = TRUE)
pF <- ggplot(cls[pair_class %in% c("homeolog", "tandem", "dispersed")],
             aes(pair_class, median_fold, fill = study)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62) +
  scale_fill_manual(values = setNames(PAL2, c("purple", "sugarcane")), name = NULL) +
  labs(x = NULL, y = "median fold difference in degree",
       title = "By duplicate type, both species",
       subtitle = "homeolog = polyploid copy; tandem = same haplotype; dispersed = other chromosome") +
  theme_f + theme(axis.text.x = element_text(size = 6.2), legend.position = "bottom")

fig <- (pA | pB) / (pC | pD) / (pE | pF) +
  plot_annotation(tag_levels = "A") + plot_layout(heights = c(1, 1, 1))
W <- 22; H <- 24
ggsave(paste0(PREFIX, ".png"), fig, width = W, height = H, units = "cm", dpi = 300)
ggsave(paste0(PREFIX, ".pdf"), fig, width = W, height = H, units = "cm")
ggsave(paste0(PREFIX, ".svg"), fig, width = W, height = H, units = "cm")
cat(sprintf("== wrote %s.{png,pdf,svg}  (%d homeolog pairs)\n", PREFIX, nrow(h)))
