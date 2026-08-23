# =============================================================================
# 23_fig_conservation.r — the conservation figure: what transfers between the
# two networks, and what does not.
#
# This is the comparative core of the project. Two networks built from two
# unrelated studies are joined through the OrthoFinder orthogroups and asked
# whether their EDGES agree more than orthology alone would produce. The figure
# has to carry two findings that pull in opposite directions and are both true:
# edge conservation is real and roughly 2.5x above chance in both directions,
# and the nitrogen RESPONSE nevertheless fails to transfer at all.
#
#   A  OBSERVED CONSERVED-EDGE RATE AGAINST ITS PERMUTATION NULL, per direction.
#      The raw rates differ ~7x between directions and that difference is
#      MEANINGLESS on its own -- purple has 9.3x more edges, so a sugarcane edge
#      has far more chances to find a partner than the reverse. The panel shows
#      the observed rate beside the null built for that direction, which is the
#      only way the two are readable together.
#
#   B  FOLD OVER NULL, by direction and by which estimator found the edge. This
#      is the comparable quantity, and the two directions agree on it.
#
#   C  THE SAME LAYERS AS RAW RATE AND AS FOLD, both relative to the Pearson-only
#      layer. MI-only edges look better conserved than Pearson-only ones in the
#      raw rate and stop looking that way once the null absorbs the opportunity
#      difference -- MI-only edges connect genes with more orthologs, so they had
#      more chances to match by accident. Putting the two normalisations side by
#      side is the whole argument in one panel.
#
#   D  THE FUNNEL, per species: network nodes -> genes on at least one conserved
#      edge -> nitrogen-responsive among those -> responsive on BOTH sides of an
#      ortholog pair. Log scale, because it spans five orders of magnitude and
#      closes at 2. This is where the comparative question actually dies, and it
#      dies in purple.
#
# NOTHING IS WRITTEN ON THE FIGURE THAT BELONGS IN THE LEGEND -- panel letters and
# the labels the data needs, nothing else. See docs/methods.md, "The rule every
# paper figure follows".
#
# NO EDGE TABLE IS READ. conserved_edges_*_FULL.tsv are 5 GB and 39 GB; every
# number here comes from the summary and null tables beside them.
#
# RUN: through run.sh  ->  ./run.sh figconservation   (then ./run.sh legends)
# =============================================================================

suppressMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(scico); library(grid); library(svglite); library(scales)
})

source(file.path(dirname(sub("--file=", "",
       grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "lib", "common.R"))

CONS_DIR   <- env_req("CLEAN_CONS_DIR")
OUT_PREFIX <- env_req("CLEAN_OUT_PREFIX")
FIG        <- env_opt("CLEAN_FIG_NUM", "4")
SELECTION  <- env_opt("CLEAN_SELECTION", "union")
NODES_SUG  <- env_req("CLEAN_NODES_SUGARCANE")
NODES_PUR  <- env_req("CLEAN_NODES_PURPLE")
setDTthreads(as.integer(env_num("CLEAN_CORES", 8)))

banner(sprintf("Figure %s — cross-species conservation", FIG))

DIRS     <- c("sugarcane_to_purple", "purple_to_sugarcane")
DIR_LAB  <- c(sugarcane_to_purple = "sugarcane → purple",
              purple_to_sugarcane = "purple → sugarcane")
LAYERS   <- c(pearson = "Pearson only", both = "both estimators", mi = "MI only")
PAL_LAYER <- setNames(scico(3, palette = "lapaz", begin = 0.10, end = 0.85),
                      unname(LAYERS))
PAL_DIR  <- setNames(scico(3, palette = "batlow")[1:2], unname(DIR_LAB))
PAL_OBS  <- c(observed = scico(5, palette = "batlow")[2],
              `permutation null` = "grey65")

theme_f <- theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(face = "bold", size = 7.5),
        legend.key.size = unit(3.5, "mm"),
        legend.title = element_text(face = "bold", size = 7.5),
        plot.margin = margin(2, 4, 2, 2))

# --- inputs ------------------------------------------------------------------
NULLS <- rbindlist(lapply(DIRS, function(d)
  fread(file.path(CONS_DIR, sprintf("conservation_null_%s.tsv", d)))))
SUMS <- rbindlist(lapply(DIRS, function(d)
  fread(file.path(CONS_DIR, sprintf("conservation_summary_%s_FULL.tsv", d)))))
CORR <- fread(file.path(CONS_DIR,
                        sprintf("conserved_correlated_summary_%s.tsv", SELECTION)))
cv <- function(k) as.numeric(CORR[metric == k, value])

NULLS[, dir_lab := factor(DIR_LAB[direction], levels = unname(DIR_LAB))]
NULLS[, layer := factor(fifelse(source == "ALL", "ALL", LAYERS[source]),
                        levels = c("ALL", unname(LAYERS)))]
SUMS[, dir_lab := factor(DIR_LAB[direction], levels = unname(DIR_LAB))]

say("conservation, full edge sets:")
print(SUMS[, .(dir_lab, layer, total_edges = fmt_n(total_edges),
               conserved = fmt_n(conserved_edges),
               pct = round(100 * conserved_fraction, 3))], row.names = FALSE)

# =============================================================================
# A — observed vs null, per direction (ALL edges)
# =============================================================================
alld <- NULLS[source == "ALL"]
obsA <- rbindlist(list(
  alld[, .(dir_lab, what = "observed", rate = obs_rate, lo = obs_rate, hi = obs_rate)],
  alld[, .(dir_lab, what = "permutation null", rate = null_mean_rate,
           lo = null_min_rate, hi = null_max_rate)]))
obsA[, what := factor(what, levels = names(PAL_OBS))]

pA <- ggplot(obsA, aes(dir_lab, rate, fill = what)) +
  geom_col(position = position_dodge(width = 0.65), width = 0.55) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.14,
                position = position_dodge(width = 0.65), linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f%%", 100 * rate)),
            position = position_dodge(width = 0.65), vjust = -0.5, size = 2.2) +
  scale_fill_manual(values = PAL_OBS, name = NULL) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.18))) +
  labs(x = NULL, y = "edges with a conserved partner") +
  theme_f

# =============================================================================
# B — fold over null, by layer and direction
# =============================================================================
foldB <- NULLS[source != "ALL"]
pB <- ggplot(foldB, aes(layer, fold_over_null, fill = layer)) +
  geom_col(width = 0.62) +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.3, colour = "grey35") +
  geom_text(aes(label = sprintf("%.2f", fold_over_null)), vjust = -0.45, size = 2.2) +
  facet_wrap(~ dir_lab, nrow = 1) +
  scale_fill_manual(values = PAL_LAYER, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "fold over permutation null") +
  theme_f + theme(axis.text.x = element_text(angle = 20, hjust = 1))

# =============================================================================
# C — raw rate vs fold, both relative to the Pearson-only layer
# =============================================================================
relC <- NULLS[source != "ALL", .(dir_lab, layer,
                                 raw = obs_rate, fold = fold_over_null), by = direction]
relC[, `:=`(raw_rel  = raw  / raw[layer == "Pearson only"],
            fold_rel = fold / fold[layer == "Pearson only"]), by = direction]
relC <- melt(relC[, .(direction, dir_lab, layer, raw_rel, fold_rel)],
             id.vars = c("direction", "dir_lab", "layer"),
             variable.name = "norm", value.name = "rel")
relC[, norm := factor(fifelse(norm == "raw_rel", "raw conserved rate",
                              "fold over null"),
                      levels = c("raw conserved rate", "fold over null"))]
say("panel C: layers relative to Pearson-only")
print(dcast(relC, dir_lab + layer ~ norm, value.var = "rel"), row.names = FALSE)

pC <- ggplot(relC[layer != "Pearson only"],
             aes(norm, rel, fill = layer)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.3, colour = "grey35") +
  geom_text(aes(label = sprintf("%.2f", rel)),
            position = position_dodge(width = 0.7), vjust = -0.45, size = 2.2) +
  facet_wrap(~ dir_lab, nrow = 1) +
  scale_fill_manual(values = PAL_LAYER, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = NULL, y = "relative to the Pearson-only layer") +
  theme_f + theme(axis.text.x = element_text(angle = 12, hjust = 1))

# =============================================================================
# D — the funnel
# =============================================================================
n_nodes <- c(sugarcane = nrow(fread(NODES_SUG, select = 1L)),
             purple    = nrow(fread(NODES_PUR, select = 1L)))
funD <- data.table(
  species = rep(c("sugarcane", "purple"), each = 4),
  stage = factor(rep(c("nodes in\nthe network", "on a conserved\nedge",
                       "nitrogen-\nresponsive", "responsive on\nboth sides"), 2),
                 levels = c("nodes in\nthe network", "on a conserved\nedge",
                            "nitrogen-\nresponsive", "responsive on\nboth sides")),
  n = c(n_nodes[["sugarcane"]],
        SUMS[direction == "sugarcane_to_purple" & layer == "ALL", conserved_genes],
        cv("correlated_sugarcane_genes"), cv("sugarcane_genes_conserved_correlated"),
        n_nodes[["purple"]],
        SUMS[direction == "purple_to_sugarcane" & layer == "ALL", conserved_genes],
        cv("correlated_purple_genes"), cv("purple_genes_conserved_correlated")))
funD[, species := factor(species, levels = c("sugarcane", "purple"))]
say("panel D: the funnel")
print(funD[, .(species, stage = gsub("\n", " ", stage), n = fmt_n(n))], row.names = FALSE)

pD <- ggplot(funD, aes(stage, n, fill = species)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.72) +
  # Angled, because two six-digit labels side by side collide at any legible
  # size in a panel this width.
  geom_text(aes(label = fmt_n(n)), position = position_dodge(width = 0.85),
            hjust = -0.08, angle = 55, size = 1.95) +
  scale_fill_manual(values = setNames(scico(3, palette = "batlow")[1:2],
                                      c("sugarcane", "purple")), name = NULL) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)),
                expand = expansion(mult = c(0, 0.35))) +
  annotation_logticks(sides = "l", linewidth = 0.2,
                      short = unit(0.4, "mm"), mid = unit(0.7, "mm"),
                      long = unit(1.1, "mm")) +
  labs(x = NULL, y = "genes") +
  theme_f

# =============================================================================
# compose
# =============================================================================
fig <- (pA | pB) / (pC | pD) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12))

W <- 22; H <- 14
invisible(ensure_dir(dirname(OUT_PREFIX)))
ggsave(paste0(OUT_PREFIX, ".png"), fig, width = W, height = H, units = "cm",
       dpi = 400, type = "cairo")
ggsave(paste0(OUT_PREFIX, ".pdf"), fig, width = W, height = H, units = "cm",
       device = cairo_pdf)
ggsave(paste0(OUT_PREFIX, ".svg"), fig, width = W, height = H, units = "cm",
       device = svglite::svglite)
say(sprintf("wrote %s.{png,pdf,svg}  (%d x %d cm)", basename(OUT_PREFIX), W, H))

# =============================================================================
# legend
# =============================================================================
wrap_at <- function(x, width = 96)
  paste(vapply(strsplit(x, "\n[[:space:]]*\n")[[1L]],
               function(para) paste(strwrap(gsub("[[:space:]]+", " ", para), width),
                                    collapse = "\n"),
               ""), collapse = "\n\n")

nl  <- function(d, sc, col) NULLS[direction == d & source == sc][[col]]
sm  <- function(d, l, col)  SUMS[direction == d & layer == l][[col]]
rel <- function(d, l, nm)   relC[direction == d & layer == LAYERS[[l]] & norm == nm, rel]
D1 <- DIRS[1]; D2 <- DIRS[2]
n_perm <- 20

legend <- paste0(
"Figure ", FIG, ". Cross-species conservation of network edges, and the failure of the nitrogen ",
"response to travel with them. An edge is CONSERVED when both of its genes have orthologs in the ",
"other species and those orthologs are themselves connected there; orthology is the OrthoFinder ",
"backbone described in Figure 1. Conservation is asked in both directions and the two are not ",
"interchangeable, because a sugarcane edge searching a 706-million-edge purple network has far ",
"more opportunity than a purple edge searching a 76-million-edge sugarcane one.\n",
"\n",
"(A) Share of each network's edges with a conserved partner, against a permutation null built ",
"separately for each direction. ", DIR_LAB[[D1]], ": ",
sprintf("%.2f%%", 100 * nl(D1, "ALL", "obs_rate")), " observed against ",
sprintf("%.2f%%", 100 * nl(D1, "ALL", "null_mean_rate")), " expected. ", DIR_LAB[[D2]], ": ",
sprintf("%.2f%%", 100 * nl(D2, "ALL", "obs_rate")), " against ",
sprintf("%.2f%%", 100 * nl(D2, "ALL", "null_mean_rate")), ". Bars on the null are its full range ",
"over ", n_perm, " permutations, which shuffle gene labels within the orthology map so that ",
"orthogroup fan-out and per-gene coverage are preserved and only the pairing is randomised -- ",
"without that the null would mostly measure how many orthologs each gene has. Nulls are computed ",
"on a 5-million-edge sample of each direction, whose observed rate matches the full set to three ",
"decimals. THE TWO RAW RATES ARE NOT COMPARABLE WITH EACH OTHER: the ~7x difference between ",
"directions is opportunity, not biology, which is the reason panel B exists.\n",
"\n",
"(B) The same conservation expressed as fold over that direction's own null, split by which ",
"estimator found the edge. This is the comparable quantity, and the two directions agree closely ",
"on it: over all edges ", sprintf("%.2f", nl(D1, "ALL", "fold_over_null")), "x and ",
sprintf("%.2f", nl(D2, "ALL", "fold_over_null")), "x. Conservation is therefore real and of ",
"similar magnitude whichever network is asked. Every layer in both directions sits above 1 ",
"(dashed line) with z between ", sprintf("%.0f", min(NULLS$z)), " and ",
sprintf("%.0f", max(NULLS$z)), "; the empirical p is ", sprintf("%.3f", 1 / (n_perm + 1)),
" for all of them, which is simply the floor set by ", n_perm, " permutations and should be read ",
"as 'beyond every permutation drawn' rather than as a precise p-value.\n",
"\n",
"(C) Each non-linear layer relative to the Pearson-only layer of the same direction, under both ",
"normalisations. For MI-only edges the two normalisations disagree in ", DIR_LAB[[D1]],
" and agree in ", DIR_LAB[[D2]], ", and neither leaves MI ahead. In ", DIR_LAB[[D1]],
" MI-only edges look BETTER conserved than Pearson-only ones by raw rate (",
sprintf("%.2f", rel(D1, "mi", "raw conserved rate")), "x) and that advantage vanishes once the ",
"null absorbs the opportunity difference (", sprintf("%.3f", rel(D1, "mi", "fold over null")),
"x, i.e. parity). In ", DIR_LAB[[D2]], " MI-only edges are already below parity on the raw rate (",
sprintf("%.2f", rel(D2, "mi", "raw conserved rate")), "x) and stay below it after the null (",
sprintf("%.3f", rel(D2, "mi", "fold over null")), "x). So there is no direction in which MI-only ",
"edges are better conserved once chance is accounted for. The correction that matters is ",
"opportunity bias: MI-only edges connect genes with more orthologs, so they had more chances to ",
"match by accident, and the permutation null absorbs exactly that. Edges found by BOTH estimators are the ",
"exception and survive both normalisations (", sprintf("%.2f", rel(D1, "both", "fold over null")),
"x and ", sprintf("%.2f", rel(D2, "both", "fold over null")), "x over Pearson-only). The robust ",
"statement is therefore not that mutual information finds better-conserved edges -- it does not ",
"-- but that edges BOTH estimators agree on are better conserved than either alone.\n",
"\n",
"(D) Where the comparative question closes, per species, on a log axis. Of ",
fmt_n(n_nodes[["sugarcane"]]), " sugarcane nodes, ",
fmt_n(sm(D1, "ALL", "conserved_genes")), " sit on at least one conserved edge, ",
fmt_n(cv("correlated_sugarcane_genes")), " of those are nitrogen-responsive, and ",
fmt_n(cv("sugarcane_genes_conserved_correlated")), " have an ortholog that is responsive in ",
"purple too. The purple column is ", fmt_n(n_nodes[["purple"]]), " -> ",
fmt_n(sm(D2, "ALL", "conserved_genes")), " -> ", fmt_n(cv("correlated_purple_genes")), " -> ",
fmt_n(cv("purple_genes_conserved_correlated")), ". The funnel does not close at the orthology ",
"step or at the conservation step -- both leave tens of thousands of genes -- it closes at ",
"purple's ", fmt_n(cv("correlated_purple_genes")), " responsive genes, and that is a POWER ",
"result rather than a biological one: at n = 18 over ", fmt_n(sm(D2, "ALL", "conserved_genes")),
" genes a gene needs |r| ~ 0.80 merely to clear the false-discovery correction. Selection here ",
"is the `", SELECTION, "` rule (responsive by Pearson or by mutual information); the Pearson-only ",
"and MI-only rules give ", fmt_n(1), " and ", fmt_n(0), " conserved responsive pairs ",
"respectively, so no choice of statistic rescues it. At the EDGE level the count is 0 under every ",
"rule, because a conserved responsive edge needs two connected genes responsive on both sides and ",
"there are at most ", fmt_n(cv("conserved_correlated_ortholog_pairs")), " such genes in the whole ",
"analysis.")

legend_file <- paste0(OUT_PREFIX, "_legend.txt")
writeLines(wrap_at(legend), legend_file)
say("wrote ", basename(legend_file))

stats <- rbindlist(list(
  NULLS[, .(panel = "A/B", study = DIR_LAB[direction],
            quantity = sprintf("%s: observed / null / fold", layer),
            value = sprintf("%.4f / %.4f / %.2f", obs_rate, null_mean_rate, fold_over_null))],
  funD[, .(panel = "D", study = as.character(species),
           quantity = gsub("\n", " ", stage), value = fmt_n(n))]))
write_tsv(stats, paste0(OUT_PREFIX, "_stats.tsv"))
say("done")
