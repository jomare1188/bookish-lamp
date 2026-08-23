# =============================================================================
# 22_fig_topology.r — the topology figure: what the two networks are shaped like.
#
# The dataset figure establishes that the re-quantification worked and the
# reproduction figure that each source study's own finding survives it. This one
# describes the objects everything after it is computed on. It exists mainly to
# stop a reader importing intuitions from sparse biological networks, because
# these are not that: mean degree is ~1,500 and ~8,300, edge density 1.4% and
# 4.8%, and purple is a single connected component. Conclusions drawn about hubs,
# modularity or centrality have to be read against that.
#
#   A  DEGREE DISTRIBUTION, as the complementary CDF P(K >= k) on log-log axes.
#      Not the raw frequency histogram the per-study diagnostic plots plot: at
#      this size the upper tail of a frequency plot is one node per bin and reads
#      as noise. The CCDF is monotone by construction, so the shape of the tail
#      is actually legible. Evaluated on a log-spaced grid rather than at every
#      unique degree, which keeps the vector output small without changing the
#      curve.
#
#   B  MEAN WEIGHT PER EDGE AGAINST DEGREE. Not node strength, which is the SUM
#      of a node's edge weights and therefore rises with degree by construction --
#      plotting that answers nothing. Dividing by degree gives the average weight
#      of the edges a node actually carries, which is what decides whether a
#      high-degree gene is well connected or merely connected to many weak
#      partners. Degrees are binned on the same log grid and averaged.
#
#   C  EDGE COMPOSITION BY LAYER -- Pearson only, both estimators, MI only. This
#      is the project's central methodological question in one panel, and it is
#      also a topology statement: the union is what the modules and every
#      centrality number are computed on. Shown as a proportion so the two
#      networks are comparable despite a 9.3x difference in edge count, with the
#      absolute counts printed.
#
#   D  MCL MODULE SIZE, again as a CCDF on log-log axes. Both networks are
#      clustered by one giant module and a long tail of small ones, and a CCDF
#      shows that without the binning choices a histogram would need.
#
# NOTHING IS WRITTEN ON THE FIGURE THAT BELONGS IN THE LEGEND -- panel letters and
# the labels the data needs, nothing else. Global metrics (nodes, edges, density,
# components, modularity) are quoted in the generated legend, not printed on the
# panels. See docs/methods.md, "The rule every paper figure follows".
#
# NO EDGE TABLE IS READ. The two edge files are 76 M and 706 M rows; everything
# here comes from the node-metric tables, the MCL module summaries and the merge
# step's summary.json, all of which are small.
#
# RUN: through run.sh  ->  ./run.sh figtopology   (then ./run.sh legends)
# =============================================================================

suppressMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(scico); library(grid); library(svglite); library(scales)
})

source(file.path(dirname(sub("--file=", "",
       grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "lib", "common.R"))

STUDIES    <- env_list("CLEAN_STUDIES", c("sugarcane", "purple"))
OUT_PREFIX <- env_req("CLEAN_OUT_PREFIX")
FIG        <- env_opt("CLEAN_FIG_NUM", "3")
NGRID      <- as.integer(env_num("CLEAN_TOPO_GRID", 300))
setDTthreads(as.integer(env_num("CLEAN_CORES", 8)))

cfg <- lapply(STUDIES, function(st) list(
  nodes  = env_req(sprintf("CLEAN_NODES_%s", toupper(st))),
  global = env_req(sprintf("CLEAN_GLOBAL_%s", toupper(st))),
  mcl    = env_req(sprintf("CLEAN_MCL_%s", toupper(st))),
  layers = env_req(sprintf("CLEAN_LAYERS_%s", toupper(st)))))
names(cfg) <- STUDIES

banner(sprintf("Figure %s — network topology", FIG))

# Study colours are the same two as the dataset figure's gene funnel, so a reader
# carries one association across the paper.
PAL_STUDY <- setNames(scico(3, palette = "batlow")[1:2], STUDIES)
LAYERS    <- c("Pearson only", "both estimators", "MI only")
PAL_LAYER <- setNames(scico(3, palette = "lapaz", begin = 0.10, end = 0.85), LAYERS)

theme_f <- theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(face = "bold", size = 7.5),
        legend.key.size = unit(3.5, "mm"),
        legend.title = element_text(face = "bold", size = 7.5),
        plot.margin = margin(2, 4, 2, 2))
log_lab <- function() trans_format("log10", math_format(10^.x))

# A log-spaced grid, deduplicated to integers. Evaluating a CCDF at every unique
# degree would put ~41,000 points in the SVG for purple alone and draw the same
# curve.
grid_of <- function(x, n = NGRID)
  unique(round(10^seq(log10(max(1, min(x))), log10(max(x)), length.out = n)))

# --- inputs ------------------------------------------------------------------
NODES <- rbindlist(lapply(STUDIES, function(st) {
  d <- fread(cfg[[st]]$nodes, select = c("gene", "degree", "strength"))
  # `strength` from igraph is the SUM of a node's edge weights, so it tracks
  # degree trivially. The per-edge mean is the quantity panel B needs.
  d[, w_mean := strength / degree]
  d[, study := st][]
}))
NODES[, study := factor(study, levels = STUDIES)]

GLOB <- rbindlist(lapply(STUDIES, function(st) {
  g <- fread(cfg[[st]]$global)
  data.table(study = st, metric = g$Metric, value = g$Value)
}))
gv <- function(st, m) GLOB[study == st & metric == m, value][1]

MCL <- rbindlist(lapply(STUDIES, function(st) {
  d <- fread(cfg[[st]]$mcl, select = c("module", "n_genes", "modularity_Q"))
  d[, study := st][]
}))
MCL[, study := factor(study, levels = STUDIES)]

LAY <- rbindlist(lapply(STUDIES, function(st) {
  txt <- paste(readLines(cfg[[st]]$layers, warn = FALSE), collapse = " ")
  num <- function(k) as.numeric(sub(sprintf('.*"%s"[[:space:]]*:[[:space:]]*([0-9.eE+-]+).*', k),
                                    "\\1", txt))
  data.table(study = st, layer = factor(LAYERS, levels = LAYERS),
             n = c(num("n_pearson_only"), num("n_both"), num("n_mi_only")))
}))
LAY[, study := factor(study, levels = STUDIES)]
LAY[, frac := n / sum(n), by = study]

say("nodes / edges / density / components:")
for (st in STUDIES)
  say(sprintf("  %-9s %s nodes | %s edges | density %s | %s components (giant %s)",
              st, fmt_n(as.integer(gv(st, "Nodes"))), fmt_n(as.numeric(gv(st, "Edges"))),
              gv(st, "Edge_density"), gv(st, "N_connected_components"),
              fmt_n(as.integer(gv(st, "Giant_component_size")))))

# =============================================================================
# A — degree CCDF
# =============================================================================
ccdf <- function(x, grid) vapply(grid, function(k) mean(x >= k), 0)
degA <- rbindlist(lapply(STUDIES, function(st) {
  d <- NODES[study == st, degree]
  g <- grid_of(d)
  data.table(study = st, k = g, p = ccdf(d, g))
}))
degA[, study := factor(study, levels = STUDIES)]
degA <- degA[p > 0]

deg_sum <- NODES[, .(nodes = .N, mean_deg = mean(degree), med_deg = median(degree),
                     max_deg = max(degree)), by = study]
say("panel A: degree")
print(deg_sum, row.names = FALSE)

pA <- ggplot(degA, aes(k, p, colour = study)) +
  geom_line(linewidth = 0.7) +
  scale_x_log10(labels = log_lab()) + scale_y_log10(labels = log_lab()) +
  scale_colour_manual(values = PAL_STUDY, name = NULL) +
  annotation_logticks(sides = "bl", linewidth = 0.2,
                      short = unit(0.4, "mm"), mid = unit(0.7, "mm"),
                      long = unit(1.1, "mm")) +
  labs(x = "degree  k", y = expression(P(K >= k))) +
  theme_f

# =============================================================================
# B — mean strength vs degree
# =============================================================================
# A coarser grid than panel A uses. The CCDF is a cumulative quantity and stays
# smooth at 300 points; a per-bin MEAN does not -- the mid-range bins of a
# 300-point log grid are narrow enough to hold a handful of nodes each, and the
# mean spikes in them. Roughly 60 bins per decade of degree is enough to show
# the shape without drawing sampling noise.
strB <- rbindlist(lapply(STUDIES, function(st) {
  d <- NODES[study == st, .(degree, w_mean)]
  g <- grid_of(d$degree, n = 60)
  d[, bin := g[findInterval(degree, g)]]
  d[, .(w = mean(w_mean), n = .N), by = bin][order(bin)][, study := st][]
}))
strB[, study := factor(study, levels = STUDIES)]
say("panel B: mean weight per edge vs degree")
for (st in STUDIES) {
  z <- NODES[study == st]
  say(sprintf("  %-9s rho(degree, per-edge weight) = %+.4f | rho(degree, strength) = %+.4f | per-edge weight %.3f -> %.3f (low to high decile)",
              st, z[, cor(degree, w_mean, method = "spearman")],
              z[, cor(degree, strength, method = "spearman")],
              z[degree <= quantile(degree, 0.1), mean(w_mean)],
              z[degree >= quantile(degree, 0.9), mean(w_mean)]))
}

# Decile profile: the summary the legend quotes. The two extreme deciles alone
# hide the shape -- the curve is flat across most of the range and moves only at
# the top -- so the whole profile is computed and its flat span reported.
DEC <- rbindlist(lapply(STUDIES, function(st) {
  z <- NODES[study == st, .(degree, w_mean)]
  z[, dec := cut(rank(degree, ties.method = "first"), 10, labels = FALSE)]
  z[, .(median_degree = as.integer(median(degree)), w = mean(w_mean)),
    by = dec][order(dec)][, study := st][]
}))
for (st in STUDIES) {
  z <- DEC[study == st]
  say(sprintf("  %-9s deciles 1-9 span %.3f-%.3f (min at decile %d), top decile %.3f",
              st, min(z[dec < 10, w]), max(z[dec < 10, w]),
              z[dec < 10][which.min(w), dec], z[dec == 10, w]))
}

# n >= 20, not 5: the mid-range bins of the log grid hold few nodes and the
# per-edge mean spikes in them, which is sampling noise rather than structure.
pB <- ggplot(strB[n >= 20], aes(bin, w, colour = study)) +
  geom_line(linewidth = 0.7) +
  scale_x_log10(labels = log_lab()) +
  scale_colour_manual(values = PAL_STUDY, name = NULL, guide = "none") +
  annotation_logticks(sides = "b", linewidth = 0.2,
                      short = unit(0.4, "mm"), mid = unit(0.7, "mm"),
                      long = unit(1.1, "mm")) +
  labs(x = "degree  k", y = "mean weight per edge") +
  theme_f

# =============================================================================
# C — edge composition by layer
# =============================================================================
say("panel C: edge composition")
print(LAY[, .(study, layer, n = fmt_n(n), pct = round(100 * frac, 2))], row.names = FALSE)

pC <- ggplot(LAY, aes(study, frac, fill = layer)) +
  geom_col(width = 0.6, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(frac > 0.03, sprintf("%.1f%%", 100 * frac), "")),
            position = position_stack(vjust = 0.5), size = 2.3,
            colour = c("white", "white", "grey15")[as.integer(LAY$layer)]) +
  scale_fill_manual(values = PAL_LAYER, name = NULL) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "share of edges") +
  theme_f

# =============================================================================
# D — MCL module size CCDF
# =============================================================================
modD <- rbindlist(lapply(STUDIES, function(st) {
  z <- MCL[study == st, n_genes]
  g <- grid_of(z)
  data.table(study = st, size = g, p = ccdf(z, g))
}))
modD[, study := factor(study, levels = STUDIES)]
modD <- modD[p > 0]

mod_sum <- MCL[, .(modules = .N, med = median(n_genes), max = max(n_genes),
                   Q = modularity_Q[1]), by = study]
say("panel D: modules")
print(mod_sum, row.names = FALSE)

pD <- ggplot(modD, aes(size, p, colour = study)) +
  geom_line(linewidth = 0.7) +
  scale_x_log10(labels = log_lab()) + scale_y_log10(labels = log_lab()) +
  scale_colour_manual(values = PAL_STUDY, name = NULL, guide = "none") +
  annotation_logticks(sides = "bl", linewidth = 0.2,
                      short = unit(0.4, "mm"), mid = unit(0.7, "mm"),
                      long = unit(1.1, "mm")) +
  labs(x = "module size  s  (genes)", y = expression(P(S >= s))) +
  theme_f

# =============================================================================
# compose
# =============================================================================
fig <- (pA | pB) / (pC | pD) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12))

W <- 20; H <- 13
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

S1 <- STUDIES[1]; S2 <- STUDIES[2]
d  <- function(st, col) deg_sum[study == st][[col]]
mo <- function(st, col) mod_sum[study == st][[col]]
ly <- function(st, l, col) LAY[study == st & layer == l][[col]]
rho_str <- vapply(STUDIES, function(st)
  NODES[study == st, cor(degree, strength, method = "spearman")], 0)
rho_w   <- vapply(STUDIES, function(st)
  NODES[study == st, cor(degree, w_mean, method = "spearman")], 0)
dmin   <- function(st) min(DEC[study == st & dec < 10, w])
dmax   <- function(st) max(DEC[study == st & dec < 10, w])
dfirst <- function(st) DEC[study == st & dec == 1, w]
dtop   <- function(st) DEC[study == st & dec == 10, w]

legend <- paste0(
"Figure ", FIG, ". Topology of the two co-expression networks. Each network is the union of a ",
"linear (Pearson) and a non-linear (Kraskov-Stogbauer-Grassberger mutual information) layer, ",
"thresholded at the same per-edge false-positive rate so the union is legitimate, built from the ",
"variance-stabilised matrix of its own study. ", S1, ": ",
fmt_n(as.integer(gv(S1, "Nodes"))), " nodes, ", fmt_n(as.numeric(gv(S1, "Edges"))), " edges, ",
"edge density ", sprintf("%.3f", as.numeric(gv(S1, "Edge_density"))), ", ",
gv(S1, "N_connected_components"), " connected components with a giant component of ",
fmt_n(as.integer(gv(S1, "Giant_component_size"))), ". ", S2, ": ",
fmt_n(as.integer(gv(S2, "Nodes"))), " nodes, ", fmt_n(as.numeric(gv(S2, "Edges"))), " edges, ",
"density ", sprintf("%.3f", as.numeric(gv(S2, "Edge_density"))), ", and a SINGLE connected ",
"component. Both are dense thresholded correlation networks, not sparse interaction graphs: ",
"mean degree is ", fmt_n(round(d(S1, "mean_deg"))), " and ", fmt_n(round(d(S2, "mean_deg"))),
". Every statement made later about hubs, modules or centrality has to be read against that, and ",
"in particular against the fact that ", S2, " has no separable components at all.\n",
"\n",
"(A) Degree distribution as the complementary cumulative distribution, P(K >= k), on log-log ",
"axes. The CCDF rather than a frequency histogram because at this size the upper tail of a ",
"histogram holds one node per bin and reads as noise, whereas the CCDF is monotone by ",
"construction. ", S1, " has a median degree of ", fmt_n(round(d(S1, "med_deg"))), " and a maximum ",
"of ", fmt_n(d(S1, "max_deg")), "; ", S2, " a median of ", fmt_n(round(d(S2, "med_deg"))),
" and a maximum of ", fmt_n(d(S2, "max_deg")), ". Neither curve is a straight line, so neither ",
"network is scale-free in the usual sense; both bend down at the tail, which is what a hard ",
"correlation threshold does to a dense graph. The two curves are offset rather than differently ",
"shaped -- ", S2, " is the same kind of object, shifted to higher degree, and the shift is a ",
"consequence of n = 18 rather than of biology: at 18 libraries the correlation needed to clear ",
"the threshold is far easier to reach by chance.\n",
"\n",
"(B) Mean weight PER EDGE of a node against its degree, degrees binned on the same log-spaced ",
"grid and bins of fewer than twenty nodes dropped, since the sparser mid-range bins spike on ",
"sampling noise. The quantity is deliberately not node strength: ",
"strength is the SUM of a node's edge weights, so it rises with degree by construction (Spearman ",
"rho = ", sprintf("%.3f", rho_str[[S1]]), " and ", sprintf("%.3f", rho_str[[S2]]),
") and plotting it answers nothing. Dividing by degree asks the question that matters -- is a ",
"high-degree gene well connected, or merely connected to many weak partners? The answer is not a ",
"trend but a step, and the step is at the very top. Across the first NINE degree deciles per-edge ",
"weight is essentially flat -- ", sprintf("%.3f", dmin(S1)), " to ", sprintf("%.3f", dmax(S1)),
" in ", S1, ", and drifting mildly DOWN from ", sprintf("%.3f", dfirst(S2)), " to ",
sprintf("%.3f", dmin(S2)), " in ", S2, " -- and then rises sharply in the tenth, to ",
sprintf("%.3f", dtop(S1)), " and ", sprintf("%.3f", dtop(S2)), ". The positive rank correlation ",
"over all nodes (Spearman rho = ", sprintf("%+.3f", rho_w[[S1]]), " and ",
sprintf("%+.3f", rho_w[[S2]]), ") is produced almost entirely by that top decile and should not ",
"be read as a gradient. For most of the network degree carries no information about edge quality; ",
"only the very highest-degree genes are also connected by stronger-than-average edges. Degree and ",
"per-edge weight are therefore reported separately throughout and never treated as one axis.\n",
"\n",
"(C) Composition of each network by which estimator found the edge, as a share of that network's ",
"edges so the two are comparable despite a ", sprintf("%.1f", as.numeric(gv(S2, "Edges")) /
as.numeric(gv(S1, "Edges"))), "x difference in edge count. In ", S1, ", ",
fmt_n(ly(S1, "Pearson only", "n")), " edges (", sprintf("%.1f%%", 100 * ly(S1, "Pearson only", "frac")),
") are found by Pearson alone, ", fmt_n(ly(S1, "both estimators", "n")), " (",
sprintf("%.1f%%", 100 * ly(S1, "both estimators", "frac")), ") by both, and ",
fmt_n(ly(S1, "MI only", "n")), " (", sprintf("%.1f%%", 100 * ly(S1, "MI only", "frac")),
") by mutual information alone; in ", S2, " the same three are ",
fmt_n(ly(S2, "Pearson only", "n")), " (", sprintf("%.1f%%", 100 * ly(S2, "Pearson only", "frac")),
"), ", fmt_n(ly(S2, "both estimators", "n")), " (",
sprintf("%.1f%%", 100 * ly(S2, "both estimators", "frac")), ") and ",
fmt_n(ly(S2, "MI only", "n")), " (", sprintf("%.1f%%", 100 * ly(S2, "MI only", "frac")), "). ",
S2, "'s larger non-linear share is NOT more non-linear biology: the two layers are matched on ",
"p-value, and at n = 18 the p implied by |r| = 0.8 is 6.7e-05 against 9.0e-12 at n = 48, so its ",
"MI floor is far softer in absolute terms. The comparison to make between the two studies is of ",
"the `both` class, which is the one an estimator cannot enter by accident.\n",
"\n",
"(D) MCL module size as a complementary cumulative distribution on log-log axes. ", S1,
" resolves into ", fmt_n(mo(S1, "modules")), " modules (median ", mo(S1, "med"), " genes, largest ",
fmt_n(mo(S1, "max")), ", modularity Q = ", mo(S1, "Q"), ") and ", S2, " into ",
fmt_n(mo(S2, "modules")), " (median ", mo(S2, "med"), ", largest ", fmt_n(mo(S2, "max")),
", Q = ", mo(S2, "Q"), "). Both partitions have the same shape -- one very large module and a ",
"long tail of small ones -- and both modularity values are low, which is what clustering a dense ",
"graph produces. ", S2, "'s largest module holds ",
sprintf("%.0f%%", 100 * mo(S2, "max") / as.integer(gv(S2, "Nodes"))), " of its nodes against ",
sprintf("%.0f%%", 100 * mo(S1, "max") / as.integer(gv(S1, "Nodes"))), " for ", S1,
", so ", S2, "'s module structure in particular should be treated as a working partition rather ",
"than as biology, and any per-module claim from it read with that in mind.")

legend_file <- paste0(OUT_PREFIX, "_legend.txt")
writeLines(wrap_at(legend), legend_file)
say("wrote ", basename(legend_file))

stats <- rbindlist(list(
  deg_sum[, .(panel = "A", study = as.character(study), quantity = "nodes / mean / median / max degree",
              value = sprintf("%s / %s / %s / %s", fmt_n(nodes), fmt_n(round(mean_deg)),
                              fmt_n(round(med_deg)), fmt_n(max_deg)))],
  data.table(panel = "B", study = STUDIES,
             quantity = "per-edge weight: deciles 1-9 span, then top decile",
             value = sprintf("%.3f-%.3f, then %.3f",
                             vapply(STUDIES, dmin, 0), vapply(STUDIES, dmax, 0),
                             vapply(STUDIES, dtop, 0))),
  LAY[, .(panel = "C", study = as.character(study), quantity = as.character(layer),
          value = sprintf("%s (%.2f%%)", fmt_n(n), 100 * frac))],
  mod_sum[, .(panel = "D", study = as.character(study),
              quantity = "modules / median / largest / Q",
              value = sprintf("%s / %s / %s / %s", fmt_n(modules), med, fmt_n(max), Q))]))
write_tsv(stats, paste0(OUT_PREFIX, "_stats.tsv"))
say("done")
