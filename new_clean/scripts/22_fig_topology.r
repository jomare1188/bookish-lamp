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
#   B  MCL MODULE SIZE, again as a CCDF on log-log axes. Both networks are
#      clustered by one giant module and a long tail of small ones, and a CCDF
#      shows that without the binning choices a histogram would need.
#
# TWO PANELS DROPPED, 2026-09-12. The mean-weight-against-degree panel was cut as
# uninformative. The edge-composition panel -- Pearson only / both / MI only --
# was cut because it no longer describes anything: the analysis runs on the
# Pearson-only graphs, so that panel is a single 100% bar by construction. The MI
# layer's contribution is a methodological result reported in docs/results.md, not
# a property of the network being described here.
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

NTERMS     <- as.integer(env_num("CLEAN_GO_NTERMS", 7))

cfg <- lapply(STUDIES, function(st) list(
  nodes  = env_req(sprintf("CLEAN_NODES_%s", toupper(st))),
  global = env_req(sprintf("CLEAN_GLOBAL_%s", toupper(st))),
  mcl    = env_req(sprintf("CLEAN_MCL_%s", toupper(st))),
  degreego = env_opt(sprintf("CLEAN_DEGREEGO_%s", toupper(st)), "")))
names(cfg) <- STUDIES

banner(sprintf("Figure %s — network topology", FIG))

# Study colours are the same two as the dataset figure's gene funnel, so a reader
# carries one association across the paper.
PAL_STUDY <- setNames(scico(3, palette = "batlow")[1:2], STUDIES)

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

# The summary table carries an `Unassigned` PSEUDO-MODULE alongside the real
# Module_NNN rows -- the genes MCL left in groups below MCL_MIN_MODULE_SIZE,
# pooled into one row. It is not a module: counting it inflates the module count
# by one, and in purple it would enter the size distribution as a 14,594-gene
# "module", second only to the real largest. Dropped from both, and reported
# separately in the legend, where it is genuinely informative.
MCL <- rbindlist(lapply(STUDIES, function(st) {
  d <- fread(cfg[[st]]$mcl, select = c("module", "n_genes", "modularity_Q"))
  d[, study := st][]
}))
UNASSIGNED <- MCL[!grepl("^Module_", module), .(study, unassigned_genes = n_genes)]
MCL <- MCL[grepl("^Module_", module)]
MCL[, study := factor(study, levels = STUDIES)]
for (st in STUDIES)
  say(sprintf("  %-9s %s named modules | %s genes unassigned (below the min module size)",
              st, fmt_n(MCL[study == st, .N]),
              fmt_n(UNASSIGNED[study == st, unassigned_genes])))

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
# =============================================================================
# B — MCL module size CCDF
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
say("panel B: modules")
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
# C — what the hubs are for, and what the periphery is for
# =============================================================================
# Panels A and B describe the degree distribution; this asks what sits at its two
# ends. 63_degree_go.r does the testing (weight01 + KS on the FULL degree ranking,
# not an over-representation test on a decile) and writes the table; nothing is
# computed here.
#
# BOTH THE p AND THE EFFECT SIZE ARE DRAWN. At n = 64,178 a KS test calls terms
# significant on small shifts -- sugarcane's "carbohydrate metabolic process" clears
# p = 2.6e-08 with its median degree at 1.04x the background, i.e. no shift at all.
# Point size is the median degree of the term's genes over the background median, so
# a term that is significant but flat is visibly small and cannot be mistaken for a
# strong one.
DG <- rbindlist(lapply(STUDIES, function(st) {
  f <- cfg[[st]]$degreego
  if (!nzchar(f) || !file.exists(f)) {
    say("NOTE: no degree-GO table for ", st, " -- run ./run.sh degreego ", st)
    return(NULL)
  }
  fread(f)[, study := st]
}), fill = TRUE)

pC <- NULL
if (length(DG) && nrow(DG)) {
  DG[, study := factor(study, levels = STUDIES)]
  setorder(DG, study, direction, pvalue)
  topC <- DG[, head(.SD, NTERMS), by = .(study, direction)]
  # Opposed on one axis -- hubs to the right, periphery to the left -- the same
  # idiom 25_fig_module_go.r uses for its two response directions, so a reader meets
  # it twice in the same form rather than learning two layouts.
  topC[, mlp := -log10(pmax(pvalue, 1e-300))]
  topC[, x := ifelse(direction == "hub", mlp, -mlp)]
  topC[, dir_lab := factor(ifelse(direction == "hub", "hubs", "periphery"),
                           levels = c("periphery", "hubs"))]
  wrap_term <- function(z, w = 40)
    vapply(z, function(t) paste(strwrap(t, w), collapse = "\n"), "")
  topC[, Term_w := wrap_term(Term)]
  topC[, key := paste(study, direction, Term, sep = "|")]
  setorder(topC, study, -x)
  topC[, key := factor(key, levels = rev(unique(key)))]
  say("panel C: degree-ranked GO, top ", NTERMS, " per direction per study")
  print(topC[, .(study, direction, Term = substr(Term, 1, 42),
                 p = signif(pvalue, 3), ratio = degree_ratio)], row.names = FALSE)

  PAL_DIR3 <- setNames(rev(scico(3, palette = "managua", begin = 0.08, end = 0.92))[c(3, 1)],
                       c("periphery", "hubs"))
  pC <- ggplot(topC, aes(x, key, fill = dir_lab, size = degree_ratio)) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey45") +
    geom_segment(aes(x = 0, xend = x, y = key, yend = key, colour = dir_lab),
                 linewidth = 0.35, show.legend = FALSE) +
    geom_point(shape = 21, colour = "grey25", stroke = 0.25) +
    facet_wrap(~ study, nrow = 1, scales = "free") +
    scale_fill_manual(values = PAL_DIR3, name = NULL) +
    scale_colour_manual(values = PAL_DIR3) +
    scale_size_continuous(range = c(1.1, 4.4),
                          name = "median degree\n/ background",
                          breaks = c(0.25, 1, 4, 9)) +
    scale_y_discrete(labels = setNames(topC$Term_w, topC$key)) +
    scale_x_continuous(labels = function(z) abs(z),
                       expand = expansion(mult = c(0.10, 0.10))) +
    labs(x = expression(-log[10](p)~", topGO weight01 + KS"), y = NULL) +
    theme_f +
    theme(axis.text.y = element_text(size = 5.4, lineheight = 0.92),
          panel.grid.major.y = element_line(linewidth = 0.2),
          legend.position = "right",
          legend.title = element_text(face = "bold", size = 6.4),
          legend.text = element_text(size = 6.2))
}

# =============================================================================
# compose
# =============================================================================
fig <- if (is.null(pC)) {
  (pA | pD) + plot_annotation(tag_levels = "A")
} else {
  (pA | pD) / pC + plot_layout(heights = c(0.62, 1)) +
    plot_annotation(tag_levels = "A")
}
fig <- fig & theme(plot.tag = element_text(face = "bold", size = 12))

W <- 20; H <- if (is.null(pC)) 7.5 else 18
invisible(ensure_dir(dirname(OUT_PREFIX)))
ggsave(paste0(OUT_PREFIX, ".png"), fig, width = W, height = H, units = "cm",
       dpi = 400, type = "cairo")
ggsave(paste0(OUT_PREFIX, ".pdf"), fig, width = W, height = H, units = "cm",
       device = cairo_pdf)
ggsave(paste0(OUT_PREFIX, ".svg"), fig, width = W, height = H, units = "cm",
       device = svglite::svglite)
say(sprintf("wrote %s.{png,pdf,svg}  (%g x %g cm)", basename(OUT_PREFIX), W, H))

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
un <- function(st) UNASSIGNED[study == st, unassigned_genes]
rho_str <- vapply(STUDIES, function(st)
  NODES[study == st, cor(degree, strength, method = "spearman")], 0)
dmax   <- function(st) max(DEC[study == st & dec < 10, w])
dfirst <- function(st) DEC[study == st & dec == 1, w]
dtop   <- function(st) DEC[study == st & dec == 10, w]

legend <- paste0(
"Figure ", FIG, ". Topology of the two co-expression networks. Each is the PEARSON-ONLY graph ",
"at |r| >= 0.8, with edge weights rescaled to [0.01, 1] and no degree reduction applied, built ",
"from the variance-stabilised matrix of its own study. The mutual-information layer is excluded: ",
"it contributed 1.14% of sugarcane's edges and 4.20% of purple's, and dropping it is a ",
"methodological result reported in the text rather than a property of the graphs described here. ",
S1, ": ",
fmt_n(as.integer(gv(S1, "Nodes"))), " nodes, ", fmt_n(as.numeric(gv(S1, "Edges"))), " edges, ",
"edge density ", sprintf("%.3f", as.numeric(gv(S1, "Edge_density"))), ", ",
gv(S1, "N_connected_components"), " connected components with a giant component of ",
fmt_n(as.integer(gv(S1, "Giant_component_size"))), ". ", S2, ": ",
fmt_n(as.integer(gv(S2, "Nodes"))), " nodes, ", fmt_n(as.numeric(gv(S2, "Edges"))), " edges, ",
"density ", sprintf("%.3f", as.numeric(gv(S2, "Edge_density"))), ", ",
gv(S2, "N_connected_components"), " connected components with a giant component of ",
fmt_n(as.integer(gv(S2, "Giant_component_size"))), " (",
sprintf("%.2f%%", 100 * as.numeric(gv(S2, "Giant_component_size")) / as.numeric(gv(S2, "Nodes"))),
" of its nodes). Both are dense thresholded correlation networks, not sparse interaction graphs: ",
"mean degree is ", fmt_n(round(d(S1, "mean_deg"))), " and ", fmt_n(round(d(S2, "mean_deg"))),
". Every statement made later about hubs, modules or centrality has to be read against that, and ",
"in particular against the fact that each network is dominated by one component that holds ",
"essentially every gene.\n",
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
"(B) MCL module size as a complementary cumulative distribution on log-log axes. ", S1,
" resolves into ", fmt_n(mo(S1, "modules")), " modules (median ", mo(S1, "med"), " genes, largest ",
fmt_n(mo(S1, "max")), ", modularity Q = ", mo(S1, "Q"), ") and ", S2, " into ",
fmt_n(mo(S2, "modules")), " (median ", mo(S2, "med"), ", largest ", fmt_n(mo(S2, "max")),
", Q = ", mo(S2, "Q"), "). These counts are of NAMED modules only: the clustering also leaves ",
fmt_n(un(S1)), " and ", fmt_n(un(S2)), " genes in groups below the minimum module size, pooled ",
"as `Unassigned` and excluded here -- a pooled leftover is not a module, and in ", S2,
" including it would put a spurious ", fmt_n(un(S2)), "-gene point into the size distribution, ",
"second only to the real largest module. Both partitions have the same shape -- one very large ",
"module and a ",
"long tail of small ones -- and both modularity values are low, which is what clustering a dense ",
"graph produces. ", S2, "'s largest module holds ",
sprintf("%.0f%%", 100 * mo(S2, "max") / as.integer(gv(S2, "Nodes"))), " of its nodes against ",
sprintf("%.0f%%", 100 * mo(S1, "max") / as.integer(gv(S1, "Nodes"))), " for ", S1,
", so ", S2, "'s module structure in particular should be treated as a working partition rather ",
"than as biology, and any per-module claim from it read with that in mind.",
if (is.null(pC)) "" else paste0(
"\n\n",
"(C) WHAT SITS AT THE TWO ENDS OF THAT DEGREE DISTRIBUTION. Every GO-annotated network node is ",
"ranked by degree and each GO term tested for concentration at one end, by a Kolmogorov-Smirnov ",
"statistic under topGO's weight01 algorithm; terms enriched among HUBS extend to the right and ",
"terms enriched at the PERIPHERY to the left, the top ", NTERMS, " each. THE RANKING IS USED ",
"WHOLE, with no cut: degree spans four orders of magnitude (", S1, "'s 10th percentile is 2, its ",
"median ", fmt_n(round(median(as.numeric(fread(cfg[[S1]]$nodes, select = 'degree')[[1]])))),
", its 90th percentile in the thousands), so any decile boundary is arbitrary and discards the ",
"middle of the data. weight01 is kept rather than moving to a GSEA implementation because it ",
"decorrelates the GO DAG -- without it a parent and its children score on the same genes and the ",
"table fills with near-duplicates.\n",
"\n",
"POINT SIZE IS THE EFFECT, and it is drawn because the p-value alone is misleading at this n. ",
"Over ", fmt_n(nrow(DG[study == S1 & direction == 'hub'])), " tested terms a KS test reaches ",
"significance on small shifts: ", S1, "'s `carbohydrate metabolic process` clears p = 2.6e-08 ",
"with its genes' median degree at 1.04x the background, which is no shift at all, while ",
"`trehalose biosynthetic process` sits at 0.15x. Size is that ratio -- the median degree of the ",
"term's genes over the background median -- so a term that is significant but flat is visibly ",
"small.\n",
"\n",
"THE TWO SPECIES AGREE. Hubs in both are photosynthesis and light harvesting (median degree ",
"5.7x and 9.8x the background), translation, intracellular protein transport, mRNA splicing and ",
"protein deubiquitination. The periphery in both is regulatory and metabolic: RNA modification, ",
"protein phosphorylation, hydrogen peroxide catabolism, carbohydrate metabolism. So the densely ",
"connected core of a co-expression network is the housekeeping machinery, and signalling and ",
"regulation sit at its edge -- in two species whose networks were built from unrelated ",
"experiments.\n",
"\n",
"READ THIS WITH FIGURE 4 PANEL B, which is largely the same contrast reached another way. Genes ",
"on a conserved edge are hubs: median degree 228 against 18 in ", S1, " and 6,916 against 376 in ",
S2, ", and 63.3% of ", S1, "'s top-degree decile sits on a conserved edge against 7.1% of the ",
"bottom. Part of that is arithmetic -- one conserved edge out of 6,677 is near-certain at a 10.4% ",
"per-edge rate and unlikely out of 2 -- so panel B's housekeeping-versus-regulation split and this ",
"one are entangled, and neither is independent evidence for the other.\n",
"\n",
"Selection is on the RAW weight01 p. weight01 deliberately makes a term's score depend on its ",
"neighbours', so its p-values are not an exchangeable family and BH's assumptions do not hold on ",
"them; a BH column is carried in the source table but does not select. One caveat measured rather ",
"than assumed: GO coverage is mildly degree-dependent, 59.7% in ", S1, "'s bottom degree decile ",
"against 66.8% in its top (62.6% and 65.4% in ", S2, "), so the periphery is slightly the less ",
"well annotated end. Source: 63_degree_go.r -> results/<study>/degree_go_<study>.tsv."))

legend_file <- paste0(OUT_PREFIX, "_legend.txt")
writeLines(wrap_at(legend), legend_file)
say("wrote ", basename(legend_file))

stats <- rbindlist(list(
  deg_sum[, .(panel = "A", study = as.character(study), quantity = "nodes / mean / median / max degree",
              value = sprintf("%s / %s / %s / %s", fmt_n(nodes), fmt_n(round(mean_deg)),
                              fmt_n(round(med_deg)), fmt_n(max_deg)))],
  mod_sum[, .(panel = "B", study = as.character(study),
              quantity = "modules / median / largest / Q",
              value = sprintf("%s / %s / %s / %s", fmt_n(modules), med, fmt_n(max), Q))]))
write_tsv(stats, paste0(OUT_PREFIX, "_stats.tsv"))
say("done")
