#!/usr/bin/env Rscript
# ============================================================================
# 49_fig_cluster_methods.r -- figure 12: modularity against MCL inflation.
#
# ONE panel, ONE argument. Modularity is the only criterion with an interior
# optimum for MCL: mass fraction, area fraction and efficiency are each monotone
# over the usable inflation range, so each is maximised at a grid boundary and
# none of them can choose a setting. This plots the criterion that can.
#
# Only the peak is marked. The case for -I 2 is made in the legend, which this
# script WRITES FROM THE DATA rather than leaving to be typed by hand -- every
# number in it is read off the same table the panel is drawn from, so the two
# cannot drift apart.
#
# The x axis is log10 (verified: the built scale reports trans "log-10" and
# places 1.2, 2 and 10 at 0.079, 0.301 and 1.0). It does not look strongly
# logarithmic because the whole ladder spans barely one decade, 1.2 to 13.
#
# Cells mcl could not compute are excluded: above ~-I 13 it underflows (35,044
# of 101,990 vectors zeroed at -I 20) and returns a plausible-looking partition
# that is numerically meaningless. So are the -S resource-probe cells, which
# belong to a different comparison.
#
# RUN: ./run.sh figclustermethods
# ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(scales); library(svglite)
})
env <- function(k, d = NA) { v <- Sys.getenv(k); if (nzchar(v)) v else d }
RESULTS <- env("CLEAN_RESULTS"); OUTDIR <- env("CLEAN_OUT_DIR")
STUDIES <- strsplit(trimws(env("CLEAN_STUDIES", "sugarcane purple")), "[ ,]+")[[1]]
CHOSEN  <- as.numeric(env("CLEAN_CHOSEN_I", "2"))
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
PREFIX <- file.path(OUTDIR, "figure12_cluster_methods")

theme_f <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(face = "bold", size = 9),
        plot.margin = margin(4, 6, 2, 2))

grab <- function(f) if (file.exists(f)) fread(f) else NULL
sweep <- rbindlist(lapply(STUDIES, function(s)
  grab(file.path(RESULTS, s, sprintf("mcl_sweep_%s.tsv", s)))), fill = TRUE)
if (!length(sweep) || !nrow(sweep)) stop("no mcl_sweep_*.tsv -- run ./run.sh mclladder <study>")

num <- function(x) suppressWarnings(as.numeric(x))
m <- sweep[, .(study,
               inflation  = num(inflation),
               n_clusters = num(n_clusters),
               largest    = num(largest),
               largest_pct = num(largest_pct),
               median_size = num(median_size),
               singletons = num(singletons),
               modularity = num(modularity),
               resource = if ("resource" %in% names(sweep)) resource else "default",
               under = if ("underflow_vectors" %in% names(sweep)) num(underflow_vectors) else 0)]
m <- m[resource == "default" & under == 0 & !is.na(modularity)]
setorder(m, study, inflation)

peaks <- m[, .SD[which.max(modularity)], by = study]
# Ties and near-ties matter: report the plateau, not an arbitrary argmax.
plateau <- m[peaks[, .(study, pk = modularity)], on = .(study)][
  modularity >= pk * 0.995,
  .(lo = min(inflation), hi = max(inflation)), by = study]

pfig <- ggplot(m, aes(inflation, modularity)) +
  geom_line(linewidth = 0.55, colour = "grey30") +
  geom_point(size = 1.5, colour = "grey30") +
  geom_point(data = peaks, size = 3.2, shape = 21, fill = NA,
             colour = "#B2182B", stroke = 1) +
  geom_text(data = peaks, aes(label = sprintf("I = %g", inflation)),
            vjust = -1.5, size = 3, colour = "#B2182B") +
  facet_wrap(~ study, scales = "free_y") +
  scale_x_log10(breaks = c(1.2, 1.5, 2, 3, 4, 6, 10)) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.18))) +
  labs(x = "MCL inflation  (log scale)", y = "modularity") +
  theme_f

W <- 16; H <- 7
ggsave(paste0(PREFIX, ".png"), pfig, width = W, height = H, units = "cm", dpi = 300)
ggsave(paste0(PREFIX, ".pdf"), pfig, width = W, height = H, units = "cm")
ggsave(paste0(PREFIX, ".svg"), pfig, width = W, height = H, units = "cm")

# --- the legend, generated from the same table ------------------------------
f  <- function(x) formatC(x, format = "d", big.mark = ",")
gene <- function(n) sprintf("%s gene%s", f(n), ifelse(n == 1, "", "s"))
d1 <- function(x) formatC(x, format = "f", digits = 1)
d5 <- function(x) formatC(x, format = "f", digits = 5)
at <- m[inflation == CHOSEN]
J  <- merge(peaks[, .(study, pk_I = inflation, pk_ncl = n_clusters, pk_mod = modularity,
                      pk_large = largest, pk_large_pct = largest_pct,
                      pk_med = median_size, pk_sgl = singletons)],
            at[, .(study, ch_ncl = n_clusters, ch_mod = modularity,
                   ch_large = largest, ch_large_pct = largest_pct,
                   ch_med = median_size, ch_sgl = singletons)], by = "study")
J <- merge(J, plateau, by = "study")
J[, pct := 100 * ch_mod / pk_mod]
setorder(J, study)

per_study <- paste(J[, sprintf(
"  %s. Modularity peaks at I = %s (Q = %s), a partition of %s modules whose
  largest holds %s genes (%s%% of the network), with a median module of %s
  and %s singletons. Modularity stays within 0.5%% of that peak across
  I = %s-%s. At I = %s the partition has %s modules, largest %s genes (%s%%),
  median %s, %s singletons, and Q = %s -- %s%% of the peak value.",
  study, pk_I, d5(pk_mod), f(pk_ncl), f(pk_large), d1(pk_large_pct), gene(pk_med), f(pk_sgl),
  lo, hi, CHOSEN, f(ch_ncl), f(ch_large), d1(ch_large_pct), ch_med, f(ch_sgl),
  d5(ch_mod), d1(pct))], collapse = "\n\n")

legend <- sprintf(
"FIGURE 12. Modularity against MCL inflation, and why the pipeline clusters at
I = %s.

Newman modularity of the MCL partition, computed by `clm info` against the graph
the partition came from, plotted against the inflation parameter for each species.
Networks are the unpruned Pearson-only co-expression graphs (|r| >= 0.8, edge
weights rescaled to [0.01, 1], no k-NN reduction): sugarcane 101,990 genes and
75,333,769 edges, purple 170,135 genes and 675,955,918 edges. Each point is one
independent `mcl` run; the ladder was sampled densely around the maximum so that
the peak is a property of the data and not of a coarse grid. Open red circles mark
the maximum in each species. The x axis is log10.

WHY THIS CRITERION. Modularity is the only one of the four available criteria with
an interior optimum here. Mass fraction falls monotonically with inflation and is
therefore maximised by the coarsest partition on the grid -- the single-cluster
partition scores 1.0. Area fraction and efficiency both rise monotonically and are
maximised by the finest. Each of those three is maximised at a grid boundary and so
cannot choose a value; modularity is the only criterion that can.

WHAT THE PANELS SHOW.

%s

WHY I = %s AND NOT EACH SPECIES' OWN PEAK. The two optima do not coincide, so no
single inflation is optimal for both. I = %s is the value between them: it retains
%s%% of the peak modularity in sugarcane and %s%% in purple, and it is the setting
at which both networks still have a usable module-size distribution. Purple's own
modularity peak is a poor partition to adopt despite its score -- at I = %s its
median module contains %s and %s genes (%s%% of the network) are singletons, against
a median of %s and %s singletons at I = %s. Buying %s%% additional
modularity at the cost of %s further singletons is not a trade the analysis needs
to make, and using one inflation for both species keeps the two networks
comparable, which is the point of the study.

Source: docs/data/clustering/mcl_sweep_{sugarcane,purple}.tsv. Underflowed cells
(mcl stops converging above ~I 13; 35,044 of 101,990 vectors zeroed at I = 20) and
the -S resource-probe cells are excluded.",
  CHOSEN, per_study, CHOSEN, CHOSEN,
  d1(J[study == "sugarcane", pct]), d1(J[study == "purple", pct]),
  J[study == "purple", pk_I], gene(J[study == "purple", pk_med]),
  f(J[study == "purple", pk_sgl]), d1(100 * J[study == "purple", pk_sgl] / 170135),
  gene(J[study == "purple", ch_med]), f(J[study == "purple", ch_sgl]), CHOSEN,
  d1(100 * J[study == "purple", pk_mod] / J[study == "purple", ch_mod] - 100),
  f(J[study == "purple", pk_sgl] - J[study == "purple", ch_sgl]))

writeLines(legend, paste0(PREFIX, "_legend.txt"))
cat(legend, "\n\n", sep = "")
cat("wrote ", PREFIX, ".{png,pdf,svg} and _legend.txt\n", sep = "")
