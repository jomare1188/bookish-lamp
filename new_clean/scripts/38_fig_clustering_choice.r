#!/usr/bin/env Rscript
# ============================================================================
# 38_fig_clustering_choice.r -- choosing MCL's granularity, with the evidence
#
# Panel A carries the claim the whole exercise rests on: the giant module is a
# node-degree artefact, not an inflation setting. The k = none line is the
# control -- if raising inflation alone flattened it, the k lines would be
# unnecessary and the diagnosis would be wrong.
# ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(scico); library(scales); library(svglite)
})
env <- function(k, d = NA) { v <- Sys.getenv(k); if (nzchar(v)) v else d }
RESULTS <- env("CLEAN_RESULTS"); OUTDIR <- env("CLEAN_OUT_DIR")
STUDIES <- strsplit(trimws(env("CLEAN_STUDIES", "sugarcane purple")), "[ ,]+")[[1]]
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
PREFIX <- file.path(OUTDIR, "figure10_clustering_choice")

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
sweep <- rbindlist(lapply(STUDIES, function(s)
  grab(file.path(RESULTS, s, sprintf("mcl_sweep_%s.tsv", s)))), fill = TRUE)
if (!nrow(sweep)) stop("no sweep results -- run ./run.sh mclsweep <study> first", call. = FALSE)
homog <- rbindlist(lapply(STUDIES, function(s)
  grab(file.path(RESULTS, s, sprintf("cluster_homogeneity_%s.tsv", s)))), fill = TRUE)

sweep[, knn_f := factor(knn, levels = c("none", sort(unique(suppressWarnings(
        as.numeric(knn[knn != "none"]))))))]
KL <- levels(sweep$knn_f)
PAL <- setNames(c("grey25", scico(max(2, length(KL) - 1), palette = "batlow")), KL)

# --- A: does the giant module survive? ---------------------------------------
pA <- ggplot(sweep, aes(inflation, largest_pct, colour = knn_f, group = knn_f)) +
  geom_hline(yintercept = 10, linetype = "22", colour = "grey60", linewidth = 0.35) +
  geom_line(linewidth = 0.55) + geom_point(size = 1.4) +
  facet_wrap(~ study, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = PAL, name = "k-NN") +
  labs(x = "inflation (-I)", y = "largest module, % of network",
       title = "Does the giant module survive?",
       subtitle = paste0("grey = no reduction, the pipeline's current setting. ",
                         "If inflation alone were the handle,\nthat line would fall on its own.")) +
  theme_f

# --- B: MCL's own criteria ----------------------------------------------------
# Area fraction IS the giant-module statistic (sum of squared cluster sizes over
# N^2); efficiency balances captured mass against that footprint.
mb <- melt(sweep[, .(study, knn_f, inflation, efficiency, area_fraction, mass_fraction)],
           id.vars = c("study", "knn_f", "inflation"), variable.name = "metric")
mb[, metric := factor(metric, levels = c("efficiency", "mass_fraction", "area_fraction"),
                      labels = c("efficiency", "mass fraction", "area fraction"))]
pB <- ggplot(mb, aes(inflation, value, colour = knn_f, group = knn_f)) +
  geom_line(linewidth = 0.5) + geom_point(size = 1.1) +
  facet_grid(metric ~ study, scales = "free_y") +
  scale_colour_manual(values = PAL, name = "k-NN") +
  labs(x = "inflation (-I)", y = NULL,
       title = "clm info: the author's own criteria",
       subtitle = "area fraction is the giant-module statistic; lower is finer-grained") +
  theme_f

# --- C: where does the partition stop moving? ---------------------------------
dist_rows <- rbindlist(lapply(STUDIES, function(s) {
  d <- file.path(env("CLEAN_WORK_DIR", "/dados04/jorge/tmp/mcl_work"),
                 sprintf("sweep_%s", s))
  fs <- list.files(d, pattern = "^dist\\.k.*\\.txt$", full.names = TRUE)
  rbindlist(lapply(fs, function(f) {
    ln <- readLines(f, warn = FALSE); ln <- ln[grepl("^d=", ln)]
    if (!length(ln)) return(NULL)
    data.table(study = s, knn = sub("^dist\\.k(.*)\\.txt$", "\\1", basename(f)),
               step = seq_along(ln),
               pct_differing = 100 * as.numeric(sub(".*^d=([0-9]+).*", "\\1",
                                 sub("\t.*", "", ln))) /
                               as.numeric(sub(".*nn=([0-9]+).*", "\\1", ln)),
               from = sub(".*n1=.*cls\\.k[^.]*\\.I", "", sub("\tn2=.*", "", ln)),
               to   = sub(".*cls\\.k[^.]*\\.I", "", ln))
  }), fill = TRUE)
}), fill = TRUE)
pC <- if (nrow(dist_rows)) {
  dist_rows[, knn_f := factor(knn, levels = KL)]
  ggplot(dist_rows, aes(step, pct_differing, colour = knn_f, group = knn_f)) +
    geom_line(linewidth = 0.5) + geom_point(size = 1.2) +
    facet_wrap(~ study, nrow = 1) +
    scale_colour_manual(values = PAL, name = "k-NN") +
    scale_x_continuous(breaks = seq_len(max(dist_rows$step))) +
    labs(x = "step along the inflation ladder", y = "% of nodes that move",
         title = "clm dist: where the partition settles",
         subtitle = "a flattening curve means further inflation is buying nothing") +
    theme_f
} else {
  ggplot() + labs(title = "clm dist: no chain output found") + theme_f
}

# --- D: are the smaller modules meaningful? -----------------------------------
# H alone rises trivially as modules shrink, so only the excess over a
# size-matched random partition is evidence.
pD <- if (!is.null(homog) && nrow(homog) && "H_excess" %in% names(homog)) {
  h <- homog[!is.na(inflation)]
  h[, knn_f := factor(knn, levels = KL)]
  ggplot(h, aes(inflation, H_excess, colour = knn_f, group = knn_f)) +
    geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.35) +
    geom_line(linewidth = 0.55) + geom_point(size = 1.4) +
    facet_wrap(~ study, nrow = 1, scales = "free_y") +
    scale_colour_manual(values = PAL, name = "k-NN") +
    labs(x = "inflation (-I)",
         y = "PFAM homogeneity above a size-matched null",
         title = "Are the smaller modules actually more coherent?",
         subtitle = "zero = no better than a random partition of the same module sizes") +
    theme_f
} else {
  ggplot() + labs(title = "homogeneity not computed yet (./run.sh clusterhomog <study>)") + theme_f
}

fig <- (pA / pB / pC / pD) +
  plot_annotation(tag_levels = "A") + plot_layout(heights = c(1, 1.7, 0.9, 1))
W <- 22; H <- 28
ggsave(paste0(PREFIX, ".png"), fig, width = W, height = H, units = "cm", dpi = 300)
ggsave(paste0(PREFIX, ".pdf"), fig, width = W, height = H, units = "cm")
ggsave(paste0(PREFIX, ".svg"), fig, width = W, height = H, units = "cm")
cat(sprintf("== wrote %s.{png,pdf,svg}  (%d sweep cells)\n", PREFIX, nrow(sweep)))
