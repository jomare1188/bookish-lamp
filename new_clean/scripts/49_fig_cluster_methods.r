#!/usr/bin/env Rscript
# ============================================================================
# 49_fig_cluster_methods.r -- MCL vs Leiden vs Louvain, all on one graph.
#
# Panel A is the figure. eff against af, one point per setting: every method
# traces a curve through the same trade-off space, so "which method" and "which
# granularity" become one comparison instead of two incomparable ones. The
# one-cluster baseline is drawn as a labelled point because it is the reductio
# for mass fraction -- it scores mf = 1.0 and is not a clustering.
#
# Panel B is the giant module against each method's own granularity knob
# (inflation for MCL, gamma for Leiden CPM). The x axes are not commensurable
# and are not pretended to be: the panel is faceted by method, and what is being
# compared is the SHAPE of each curve and where it lands, not x for x.
#
# Panel C is the check that keeps A honest. eff says how well a partition
# explains the graph; it cannot say whether the modules mean anything. PFAM
# Sorensen-Dice homogeneity above a SIZE-MATCHED null is the independent axis,
# and a method has to do well on both to be adopted.
#
# RUN: ./run.sh figclustermethods
# ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(scico); library(scales); library(svglite)
})
env <- function(k, d = NA) { v <- Sys.getenv(k); if (nzchar(v)) v else d }
RESULTS <- env("CLEAN_RESULTS"); OUTDIR <- env("CLEAN_OUT_DIR")
STUDIES <- strsplit(trimws(env("CLEAN_STUDIES", "sugarcane purple")), "[ ,]+")[[1]]
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
PREFIX <- file.path(OUTDIR, "figure12_cluster_methods")

theme_f <- theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(face = "bold", size = 7.5),
        legend.key.size = unit(3.5, "mm"),
        legend.title = element_text(face = "bold", size = 7.5),
        plot.margin = margin(2, 4, 2, 2),
        plot.subtitle = element_text(size = 6.5, colour = "grey35"),
        plot.tag = element_text(face = "bold", size = 10))

grab <- function(f) if (file.exists(f)) fread(f) else NULL

cmp <- rbindlist(lapply(STUDIES, function(s)
  grab(file.path(RESULTS, s, sprintf("cluster_methods_%s.tsv", s)))), fill = TRUE)
if (!length(cmp) || !nrow(cmp)) stop("no cluster_methods_*.tsv -- run ./run.sh clustercompare <study>")

hom <- rbindlist(lapply(STUDIES, function(s)
  grab(file.path(RESULTS, s, sprintf("cluster_homogeneity_%s.tsv", s)))), fill = TRUE)

# Numeric where it matters; param stays character because it means different
# things per method and must never be silently pooled across them.
for (v in c("eff", "mf", "af", "mod", "largest_pct", "n_clusters"))
  if (v %in% names(cmp)) cmp[[v]] <- suppressWarnings(as.numeric(cmp[[v]]))
cmp[, param_num := suppressWarnings(as.numeric(param))]
cmp[, method := factor(method,
      levels = c("mcl", "leiden_cpm", "leiden_mod", "louvain", "fastgreedy", "baseline"),
      labels = c("MCL", "Leiden CPM", "Leiden modularity", "Louvain",
                 "fast-greedy", "baseline (1 cluster)"))]
PAL <- c("MCL" = "#B2182B", "Leiden CPM" = "#2166AC", "Leiden modularity" = "#4393C3",
         "Louvain" = "#5AAE61", "fast-greedy" = "#8073AC",
         "baseline (1 cluster)" = "grey45")

is_base <- cmp$method == "baseline (1 cluster)"

pA <- ggplot(cmp[!is_base], aes(af, eff, colour = method)) +
  geom_path(aes(group = interaction(study, method)), linewidth = 0.35, alpha = 0.8) +
  geom_point(size = 1.1) +
  geom_point(data = cmp[is_base], shape = 4, size = 2, stroke = 0.7) +
  geom_text(data = cmp[is_base], aes(label = "1 cluster: mf = 1.0"),
            hjust = 1.08, size = 2, show.legend = FALSE) +
  facet_wrap(~ study, scales = "free") +
  scale_colour_manual(values = PAL, name = NULL) +
  scale_x_continuous(trans = "sqrt", labels = label_number(accuracy = 0.01)) +
  labs(tag = "A", x = "area fraction  (giant-module statistic, sqrt scale)",
       y = "efficiency  (primary criterion)",
       subtitle = "every point scored by clm info against the SAME unpruned graph") +
  theme_f + theme(legend.position = "bottom")

pB <- ggplot(cmp[!is_base & !is.na(param_num) & param_num > 0],
             aes(param_num, largest_pct, colour = method)) +
  geom_line(aes(group = interaction(study, method)), linewidth = 0.35) +
  geom_point(size = 1) +
  facet_grid(study ~ method, scales = "free_x") +
  scale_colour_manual(values = PAL, guide = "none") +
  scale_x_log10() +
  scale_y_continuous(labels = label_percent(scale = 1, accuracy = 1)) +
  labs(tag = "B", x = "granularity knob (inflation for MCL, gamma for Leiden CPM) - log scale",
       y = "largest module, % of nodes",
       subtitle = "the x axes are not commensurable; compare the shapes, not the values") +
  theme_f

pC <- if (!is.null(hom) && nrow(hom) && "H_excess" %in% names(hom)) {
  h <- copy(hom)
  h[, key := sub("^cls\\.", "", clustering)]
  c2 <- copy(cmp)
  c2[, key := paste0(as.character(method), ".", param)]
  # Join on the manifest's own method/param when 37 recorded them; otherwise the
  # panel degrades to whatever matched rather than inventing a mapping.
  if ("method" %in% names(h)) {
    h2 <- merge(h[, .(study, method, param, H_excess)],
                cmp[, .(study, method_lab = method, param, eff)],
                by.x = c("study", "param"), by.y = c("study", "param"),
                allow.cartesian = TRUE)
    ggplot(h2, aes(eff, H_excess, colour = method_lab)) +
      geom_point(size = 1.1) +
      facet_wrap(~ study, scales = "free") +
      scale_colour_manual(values = PAL, name = NULL) +
      labs(tag = "C", x = "efficiency", y = "PFAM homogeneity above a size-matched null",
           subtitle = "the independent axis: does the partition mean anything biologically") +
      theme_f + theme(legend.position = "bottom")
  } else ggplot() + labs(tag = "C", title = "homogeneity table has no method column") + theme_f
} else {
  ggplot() + labs(tag = "C",
    title = "homogeneity not computed yet (./run.sh clusterhomog <study>)") + theme_f
}

fig <- pA / pB / pC + plot_layout(heights = c(1, 1, 0.9))
W <- 18; H <- 21
ggsave(paste0(PREFIX, ".png"), fig, width = W, height = H, units = "cm", dpi = 300)
ggsave(paste0(PREFIX, ".pdf"), fig, width = W, height = H, units = "cm")
ggsave(paste0(PREFIX, ".svg"), fig, width = W, height = H, units = "cm")
cat("wrote ", PREFIX, ".{png,pdf,svg}\n", sep = "")
