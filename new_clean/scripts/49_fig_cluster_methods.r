#!/usr/bin/env Rscript
# ============================================================================
# 49_fig_cluster_methods.r -- why -I 2, in two panels and nothing else.
#
# This figure makes ONE argument. Everything that does not serve it was cut:
# the eff/af trade-off curves, the giant-module-vs-parameter panel and the
# homogeneity scatter are all real results, but they belong to the write-up,
# not to the figure that has to justify a parameter choice.
#
# PANEL A -- the decision. Modularity against inflation, per species. This is
# the ONLY criterion with an interior optimum for MCL: mass fraction, area
# fraction and efficiency are each monotone over the usable inflation range, so
# each is maximised at a grid boundary and none of them can choose a setting.
#
# The peaks are at -I 1.7 (sugarcane) and -I 3-4 (purple), so -I 2 is NOT the
# optimum in either species -- it is the value that sits BETWEEN them, within
# ~2% of both. That is the real argument for using one setting across two
# networks, and the panel is drawn to show exactly that rather than to imply a
# peak at 2 which is not there.
#
# PANEL B -- the same criterion, every method, on ONE x axis. Parameter values
# are not commensurable across methods (inflation is not gamma), but the number
# of clusters is, so that is the axis. It shows three things at once:
#   * MCL's interior peak is real, and MCL's MAXIMUM modularity exceeds Leiden
#     CPM's maximum in both species (0.0821 vs 0.0803 sugarcane, 0.1529 vs
#     0.1426 purple) -- reached at a COARSER granularity than the range where
#     the two methods overlap.
#   * Which method is higher INSIDE the overlap is species-dependent, and is not
#     claimed either way here: at ~19,300 clusters sugarcane's Leiden CPM is
#     above MCL (0.0799 vs 0.0664), while purple's MCL is above Leiden CPM
#     throughout. The panel shows both curves and lets them speak.
#   * Louvain and Leiden-modularity sit far left with much higher modularity --
#     50 to 1,107 clusters for 102k-170k genes. That corner is the resolution
#     limit, and the write-up shows those partitions carry no PFAM signal above
#     a size-matched null. High modularity there is a warning, not a result.
#
# Cells mcl could not compute are excluded: above ~-I 13 it underflows (35,044
# of 101,990 vectors zeroed at -I 20) and returns a plausible-looking partition
# that is numerically meaningless.
#
# RUN: ./run.sh figclustermethods
# ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(scales); library(svglite)
})
env <- function(k, d = NA) { v <- Sys.getenv(k); if (nzchar(v)) v else d }
RESULTS <- env("CLEAN_RESULTS"); OUTDIR <- env("CLEAN_OUT_DIR")
STUDIES <- strsplit(trimws(env("CLEAN_STUDIES", "sugarcane purple")), "[ ,]+")[[1]]
CHOSEN  <- as.numeric(env("CLEAN_CHOSEN_I", "2"))
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
PREFIX <- file.path(OUTDIR, "figure12_cluster_methods")

theme_f <- theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(face = "bold", size = 8),
        legend.key.size = unit(3.5, "mm"),
        plot.margin = margin(2, 4, 2, 2),
        plot.subtitle = element_text(size = 6.8, colour = "grey35"),
        plot.tag = element_text(face = "bold", size = 10))

grab <- function(f) if (file.exists(f)) fread(f) else NULL

sweep <- rbindlist(lapply(STUDIES, function(s)
  grab(file.path(RESULTS, s, sprintf("mcl_sweep_%s.tsv", s)))), fill = TRUE)
cmp <- rbindlist(lapply(STUDIES, function(s)
  grab(file.path(RESULTS, s, sprintf("cluster_methods_%s.tsv", s)))), fill = TRUE)
if (!length(sweep) || !nrow(sweep)) stop("no mcl_sweep_*.tsv -- run ./run.sh mclladder <study>")

# --- panel A: MCL at default resources, numerically valid cells only --------
m <- sweep[, .(study, inflation = as.numeric(inflation),
               n_clusters = as.numeric(n_clusters),
               modularity = as.numeric(modularity),
               resource = if ("resource" %in% names(sweep)) resource else "default",
               under = if ("underflow_vectors" %in% names(sweep))
                         as.numeric(underflow_vectors) else 0)]
m <- m[resource == "default" & under == 0 & !is.na(modularity)]
setorder(m, study, inflation)

peaks <- m[, .SD[which.max(modularity)], by = study]
at2   <- m[inflation == CHOSEN]
lab   <- merge(peaks[, .(study, peak_I = inflation, peak_mod = modularity)],
               at2[, .(study, mod2 = modularity)], by = "study")
lab[, pct := 100 * mod2 / peak_mod]
# Ties matter: purple's -I 3 and -I 4 give identical modularity to five decimals
# on genuinely different partitions (22,545 vs 35,561 clusters). Report the
# plateau rather than an arbitrary argmax.
tied <- m[lab, on = .(study), .(study, inflation, modularity, peak_mod)][
           abs(modularity - peak_mod) < 1e-5]
plateau <- tied[, .(peak_lab = paste(sort(unique(inflation)), collapse = "-")), by = study]
lab <- merge(lab, plateau, by = "study")
pcts <- lab[, sprintf("%.1f", pct)]
keep_txt <- if (uniqueN(pcts) == 1) {
  sprintf("I = %g keeps %s%% of the peak in both species.", CHOSEN, pcts[1])
} else {
  paste0(sprintf("I = %g keeps ", CHOSEN),
         paste(sprintf("%s%% (%s)", pcts, lab$study), collapse = ", "), " of the peak.")
}

pA <- ggplot(m, aes(inflation, modularity)) +
  geom_vline(xintercept = CHOSEN, colour = "#B2182B", linewidth = 0.4, linetype = "22") +
  geom_line(linewidth = 0.5, colour = "grey30") +
  geom_point(size = 1.3, colour = "grey30") +
  geom_point(data = peaks, size = 2.8, shape = 21, fill = NA,
             colour = "#2166AC", stroke = 0.9) +
  geom_point(data = at2, size = 2.2, colour = "#B2182B") +
  geom_text(data = lab, aes(x = peak_I, y = peak_mod,
                            label = sprintf("peak I = %s", peak_lab)),
            vjust = -1.3, size = 2.4, colour = "#2166AC") +
  facet_wrap(~ study, scales = "free_y") +
  scale_x_log10(breaks = c(1.2, 1.5, 2, 3, 4, 6, 10)) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.20))) +
  labs(tag = "A", x = "MCL inflation  (log scale)", y = "modularity",
       subtitle = paste0("the only criterion with an interior optimum for MCL. mass ",
                         "fraction, area fraction and efficiency are all monotone over ",
                         "this range,\nso each is maximised at a grid boundary and none of ",
                         "them can choose a setting.  ", keep_txt)) +
  theme_f

# --- panel B: every method, x = number of clusters --------------------------
b <- cmp[, .(study, method, param = as.character(param),
             resource = if ("resource" %in% names(cmp)) as.character(resource) else "default",
             n_clusters = as.numeric(n_clusters), modularity = as.numeric(mod))]
b <- b[method != "baseline" & !is.na(modularity) & n_clusters > 1]
# DEFAULT RESOURCES ONLY. The -S probe cells are in this table too, and one of
# them is sugarcane's highest-modularity mcl cell (0.08239 at -S 4000 against
# 0.08212 at the default). Letting it into a comparison BETWEEN METHODS would
# compare mcl-with-more-memory against everything else at stock settings.
b <- b[resource %in% c("", "default", NA)]
# drop the underflowed mcl cells exactly as panel A does
if ("underflow_vectors" %in% names(sweep)) {
  badI <- unique(sweep[as.numeric(underflow_vectors) > 0,
                       .(study, method = "mcl", param = as.character(inflation))])
  if (nrow(badI)) b <- b[!badI, on = .(study, method, param)]
}
b[, method := factor(method,
      levels = c("mcl", "leiden_cpm", "leiden_mod", "louvain"),
      labels = c("MCL", "Leiden CPM", "Leiden modularity", "Louvain"))]
b <- b[!is.na(method)]
setorder(b, study, method, n_clusters)
lines  <- b[method %in% c("MCL", "Leiden CPM")]
points <- b[method %in% c("Leiden modularity", "Louvain")]
PAL <- c("MCL" = "#B2182B", "Leiden CPM" = "#2166AC",
         "Leiden modularity" = "#762A83", "Louvain" = "#1B7837")

pB <- ggplot(mapping = aes(n_clusters, modularity, colour = method)) +
  geom_line(data = lines, linewidth = 0.5) +
  geom_point(data = lines, size = 1.2) +
  geom_point(data = points, size = 2.8, shape = 17) +
  facet_wrap(~ study, scales = "free") +
  scale_colour_manual(values = PAL, name = NULL) +
  scale_x_log10(labels = label_number(big.mark = ",", accuracy = 1)) +
  labs(tag = "B",
       x = "number of clusters  (log scale) - the only axis the methods share",
       y = "modularity",
       subtitle = paste("MCL's maximum exceeds Leiden CPM's in both species, at a coarser",
                        "granularity than where the two overlap.\nLouvain and",
                        "Leiden-modularity win this criterion outright with 50-1,107 clusters",
                        "for 102k-170k genes -\nthat corner is the resolution limit, and those",
                        "partitions carry no PFAM signal above a size-matched null.")) +
  theme_f + theme(legend.position = "bottom")

fig <- pA / pB
W <- 18; H <- 13
ggsave(paste0(PREFIX, ".png"), fig, width = W, height = H, units = "cm", dpi = 300)
ggsave(paste0(PREFIX, ".pdf"), fig, width = W, height = H, units = "cm")
ggsave(paste0(PREFIX, ".svg"), fig, width = W, height = H, units = "cm")

cat("\npanel A -- modularity peak vs the chosen setting\n")
print(lab[, .(study, peak_I, peak_mod = round(peak_mod, 5),
              chosen = CHOSEN, mod_at_chosen = round(mod2, 5),
              pct_of_peak = round(pct, 1))], row.names = FALSE)
cat("\npanel B -- best modularity per method\n")
print(b[, .SD[which.max(modularity)], by = .(study, method)][
        order(study, -modularity), .(study, method, param, n_clusters,
                                     modularity = round(modularity, 5))],
      row.names = FALSE)
cat("\nwrote ", PREFIX, ".{png,pdf,svg}\n", sep = "")
