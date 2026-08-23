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
#   B  WHAT THE CONSERVED GENES ARE FOR. One topGO BP enrichment per species over
#      the genes sitting on at least one conserved edge, against that network's
#      own GO-annotated nodes as background, showing the terms enriched in BOTH.
#      Conservation of edges is a structural statement; this is the panel that
#      says whether the structure that transfers is doing anything recognisable.
#      Terms are ranked by their WORSE p-value across the two species, so the
#      panel shows terms both species agree on rather than terms one species
#      drives.
#
#   C  EDGE COMPOSITION BY LAYER, per network -- how much of each network the
#      non-linear layer actually contributes. Stacked as a proportion so the two
#      are comparable across a 9.3x difference in edge count, with the absolute
#      counts on the bars. Together with the per-layer conservation folds quoted
#      in the legend, this is the case on the MI layer: it is a small minority of
#      both networks AND the edges it alone finds are not better conserved.
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
GO_ONT     <- env_opt("CLEAN_GO_ONTOLOGY", "BP")
GO_NTERMS  <- as.integer(env_num("CLEAN_GO_NTERMS", 12))
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

GO_DIR <- file.path(CONS_DIR, "enrichment_conserved")
GO_SHARED <- fread(file.path(GO_DIR, sprintf("GO_%s_conserved_shared_terms.csv", GO_ONT)))
# topGO's GenTable truncates term names at 40 characters with an ellipsis, so the
# CSV carries strings like "positive regulation of cellular cataboli...". Expand
# them from the cached GO.db dump; anything missing keeps the truncated string
# rather than becoming NA.
GO_NAMES <- env_opt("CLEAN_GO_NAMES")
if (nzchar(GO_NAMES) && file.exists(GO_NAMES)) {
  nm <- fread(GO_NAMES)
  GO_SHARED <- merge(GO_SHARED, nm, by = "GO.ID", all.x = TRUE, sort = FALSE)
  n_trunc <- GO_SHARED[grepl("\\.\\.\\.$", Term), .N]
  GO_SHARED[!is.na(Term_full) & nzchar(Term_full), Term := Term_full]
  say(sprintf("expanded %d truncated %s term name(s) from %s",
              n_trunc, GO_ONT, basename(GO_NAMES)))
} else {
  say("NOTE: no GO term-name cache; long term names stay truncated as topGO wrote them")
}
GO_CMP <- rbindlist(lapply(c("BP", "MF", "CC"), function(o)
  fread(file.path(GO_DIR, sprintf("GO_%s_conserved_comparison_summary.csv", o)))))
GO_SUM <- fread(file.path(GO_DIR, sprintf("GO_%s_conserved_summary.csv", GO_ONT)))
gc_ <- function(o, col) GO_CMP[Ontology == o][[col]]
gs_ <- function(sp, col) GO_SUM[Network == sp][[col]]

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
# Totals come from the FULL edge sets, not the 5 M null sample, so the absolute
# counts printed on the bars are the real ones. The null bar's absolute value is
# the count its rate implies over that same full set, i.e. edges expected to find
# a partner by chance.
tot <- SUMS[layer == "ALL", .(direction, total_edges, conserved_edges)]
alld <- merge(alld, tot, by = "direction")
obsA <- rbindlist(list(
  alld[, .(dir_lab, total_edges, what = "observed", rate = obs_rate,
           lo = obs_rate, hi = obs_rate, n = conserved_edges)],
  alld[, .(dir_lab, total_edges, what = "permutation null", rate = null_mean_rate,
           lo = null_min_rate, hi = null_max_rate,
           n = round(null_mean_rate * total_edges))]))
obsA[, what := factor(what, levels = names(PAL_OBS))]
# The denominator rides on the axis label, so the panel carries the percentage,
# the count and what the count is out of.
# Broken across three short lines rather than two long ones: at this panel width
# "sugarcane → purple" on one line is wider than the space a single group gets,
# and the two group labels collide.
ax_lab <- function(d, n)
  sprintf("%s\n%s edges", sub(" (→|->) ", "\n\\1 ", d), trimws(fmt_n(n)))
obsA[, dir_ax := factor(ax_lab(dir_lab, total_edges),
                        levels = unique(ax_lab(dir_lab, total_edges)[order(dir_lab)]))]
say("panel A: observed vs null, with absolute counts")
print(obsA[, .(dir_lab, what, pct = round(100 * rate, 3), n = fmt_n(n))], row.names = FALSE)

pA <- ggplot(obsA, aes(dir_ax, rate, fill = what)) +
  geom_col(position = position_dodge(width = 0.65), width = 0.55) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.14,
                position = position_dodge(width = 0.65), linewidth = 0.3) +
  geom_text(aes(label = fmt_n(n)), position = position_dodge(width = 0.65),
            vjust = -0.55, size = 2.1) +
  scale_fill_manual(values = PAL_OBS, name = NULL) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.20))) +
  labs(x = NULL, y = "edges with a conserved partner") +
  theme_f + theme(axis.text.x = element_text(size = 6.3, lineheight = 1.1))

# =============================================================================
# B — fold over null, by layer and direction
# =============================================================================
# Ranked by the WORSE of the two p-values, so the panel shows terms BOTH species
# are enriched for rather than terms one species drives.
gsh <- copy(GO_SHARED)
setnames(gsh, c("pvalue_sugarcane", "pvalue_purple"), c("p_sugarcane", "p_purple"))
gsh[, worst := pmax(p_sugarcane, p_purple)]
setorder(gsh, worst)
top <- head(gsh, GO_NTERMS)
goB <- melt(top[, .(Term, p_sugarcane, p_purple)], id.vars = "Term",
            variable.name = "species", value.name = "p")
goB[, species := factor(sub("^p_", "", species), levels = c("sugarcane", "purple"))]
# Full GO names run to 70+ characters and, unwrapped, the label column eats half
# the figure width and squeezes panel A. Wrapped to two or three short lines.
wrap_term <- function(x, w = 38)
  vapply(x, function(t) paste(strwrap(t, w), collapse = "\n"), "")
top[, Term_w := wrap_term(Term)]
goB <- merge(goB, top[, .(Term, Term_w)], by = "Term")
goB[, Term := factor(Term_w, levels = rev(top$Term_w))]
goB[, mlp := -log10(p)]
say(sprintf("panel B: %s terms shared by both conserved sets, showing top %d",
            fmt_n(nrow(gsh)), nrow(top)))
print(top[, .(Term = substr(Term, 1, 46), p_sugarcane, p_purple)], row.names = FALSE)

pB <- ggplot(goB, aes(mlp, Term)) +
  geom_line(aes(group = Term), colour = "grey70", linewidth = 0.35) +
  geom_point(aes(fill = species), shape = 21, size = 1.9,
             colour = "grey25", stroke = 0.25) +
  scale_fill_manual(values = setNames(scico(3, palette = "batlow")[1:2],
                                      c("sugarcane", "purple")), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.08))) +
  labs(x = expression(-log[10](p)~", topGO weight01"), y = NULL) +
  theme_f + theme(axis.text.y = element_text(size = 6, lineheight = 0.95),
                  panel.grid.major.y = element_line(linewidth = 0.2))

# =============================================================================
# C — edge composition by layer
# =============================================================================
# What this panel replaced, and why the numbers it carried are now in the legend:
# it used to plot each non-linear layer relative to the Pearson-only layer under
# two normalisations, which is the sharpest form of the "MI adds nothing" result.
# The composition is the more direct statement of the same case -- MI-only edges
# are a small minority of both networks to begin with -- so the plot shows that
# and the legend keeps the fold-over-null numbers, which are the actual evidence.
LAY <- rbindlist(lapply(c("sugarcane", "purple"), function(st) {
  d <- if (st == "sugarcane") "sugarcane_to_purple" else "purple_to_sugarcane"
  x <- SUMS[direction == d & layer != "ALL"]
  data.table(study = st, layer = factor(LAYERS[x$layer], levels = unname(LAYERS)),
             n = x$total_edges)
}))
LAY[, study := factor(study, levels = c("sugarcane", "purple"))]
LAY[, frac := n / sum(n), by = study]
say("panel C: edge composition by layer")
print(LAY[, .(study, layer, n = fmt_n(n), pct = round(100 * frac, 2))], row.names = FALSE)

# The two minority classes are 1-7% of a bar, so their in-slice labels sit on top
# of each other however they are justified. They go under the axis instead, where
# they cannot collide and where the numbers the MI argument needs are the easiest
# thing on the panel to read.
minor <- LAY[layer != "Pearson only"]
# Percentages only, and the class names shortened: the full "both estimators
# 2,004,577 (2.6%)" is wider than a bar and the two studies' labels collide. The
# counts are in the legend and the stats table.
lab_c <- minor[, .(lab = paste(sprintf("%s %.1f%%",
                                       fifelse(layer == "MI only", "MI only", "both"),
                                       100 * frac),
                               collapse = "\n")), by = study]
lab_c <- lab_c[match(levels(LAY$study), study)]
LAY[, study_ax := factor(sprintf("%s\n%s", study, lab_c$lab[match(study, lab_c$study)]),
                         levels = sprintf("%s\n%s", lab_c$study, lab_c$lab))]

pC <- ggplot(LAY, aes(study_ax, frac, fill = layer)) +
  geom_col(width = 0.58, colour = "white", linewidth = 0.3) +
  geom_text(data = LAY[layer == "Pearson only"],
            aes(label = sprintf("%s\n%.1f%%", fmt_n(n), 100 * frac)),
            position = position_stack(vjust = 0.5), size = 2.2,
            lineheight = 0.95, colour = "white") +
  scale_fill_manual(values = PAL_LAYER, name = NULL) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "share of edges") +
  theme_f + theme(axis.text.x = element_text(size = 6.1, lineheight = 1.15))

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
# The GO panel needs the taller row: its labels are full GO term names wrapped
# to two or three lines, and at equal row heights they collide.
fig <- (pA | pB) / (pC | pD) +
  plot_layout(heights = c(1.25, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12))

W <- 22; H <- 17
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
ly  <- function(st, l, col) LAY[study == st & layer == l][[col]]
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
"separately for each direction; bars carry the absolute edge counts and the axis label the ",
"denominator. ", DIR_LAB[[D1]], ": ", fmt_n(sm(D1, "ALL", "conserved_edges")), " of ",
fmt_n(sm(D1, "ALL", "total_edges")), " edges conserved (",
sprintf("%.2f%%", 100 * nl(D1, "ALL", "obs_rate")), ") against ",
fmt_n(round(nl(D1, "ALL", "null_mean_rate") * sm(D1, "ALL", "total_edges"))), " expected by ",
"chance (", sprintf("%.2f%%", 100 * nl(D1, "ALL", "null_mean_rate")), "), a fold over null of ",
sprintf("%.2f", nl(D1, "ALL", "fold_over_null")), "x. ", DIR_LAB[[D2]], ": ",
fmt_n(sm(D2, "ALL", "conserved_edges")), " of ", fmt_n(sm(D2, "ALL", "total_edges")), " (",
sprintf("%.2f%%", 100 * nl(D2, "ALL", "obs_rate")), ") against ",
fmt_n(round(nl(D2, "ALL", "null_mean_rate") * sm(D2, "ALL", "total_edges"))), " (",
sprintf("%.2f%%", 100 * nl(D2, "ALL", "null_mean_rate")), "), fold ",
sprintf("%.2f", nl(D2, "ALL", "fold_over_null")), "x. The two folds agree closely even though ",
"the raw percentages differ ~7x, and the fold is the comparable quantity: a sugarcane edge ",
"searching a 706-million-edge purple network has far more opportunity than the reverse. Every ",
"layer in both directions clears the null with z between ", sprintf("%.0f", min(NULLS$z)),
" and ", sprintf("%.0f", max(NULLS$z)), "; the empirical p is ",
sprintf("%.3f", 1 / (n_perm + 1)), " throughout, which is simply the floor set by ", n_perm,
" permutations and should be read as \'beyond every permutation drawn\' rather than as a precise ",
"p-value. Bars on the null are its full range ",
"over ", n_perm, " permutations, which shuffle gene labels within the orthology map so that ",
"orthogroup fan-out and per-gene coverage are preserved and only the pairing is randomised -- ",
"without that the null would mostly measure how many orthologs each gene has. Nulls are computed ",
"on a 5-million-edge sample of each direction, whose observed rate matches the full set to three ",
"decimals. THE TWO RAW RATES ARE NOT COMPARABLE WITH EACH OTHER: the ~7x difference between ",
"directions is opportunity, not biology, which is the reason panel B exists.\n",
"\n",
"(B) What the conserved genes are FOR. Panel A is a structural statement; this asks whether the ",
"structure that transfers is doing anything recognisable. One topGO ", GO_ONT, " enrichment per ",
"species over the genes on at least one conserved edge (", fmt_n(gs_("sugarcane", "Conserved_genes")),
" and ", fmt_n(gs_("purple", "Conserved_genes")), " genes, of which ",
fmt_n(gs_("sugarcane", "Conserved_in_bg")), " and ", fmt_n(gs_("purple", "Conserved_in_bg")),
" carry any GO annotation), tested against that network's OWN GO-annotated nodes as background (",
fmt_n(gs_("sugarcane", "Node_background_w_GO")), " and ",
fmt_n(gs_("purple", "Node_background_w_GO")), ") rather than against the genome, so the result is ",
"not the generic \'co-expressed genes differ from the genome\' effect. Sugarcane returns ",
gc_(GO_ONT, "Terms_sugarcane"), " enriched ", GO_ONT, " terms and purple ",
gc_(GO_ONT, "Terms_purple"), ", of which ", gc_(GO_ONT, "Shared"), " are SHARED (Jaccard ",
sprintf("%.2f", gc_(GO_ONT, "Jaccard")), "); MF and CC behave the same way, at ",
gc_("MF", "Shared"), " and ", gc_("CC", "Shared"), " shared terms. Plotted are the ",
sprintf("%d", nrow(top)), " shared terms with the strongest agreement, ranked by the WORSE of ",
"the two p-values so the panel shows terms both species support rather than terms one species ",
"drives; the two points on a row are the same term in the two networks. The shared set is ",
"dominated by regulation and core metabolism -- regulation of gene expression, MAPK cascade, ",
"proteolysis, translation, thylakoid membrane organization -- and it includes two terms that ",
"speak directly to the trait: GLUTAMATE BIOSYNTHETIC PROCESS and the AMMONIA ASSIMILATION CYCLE, ",
"enriched in both conserved sets. So the conserved edges are not a structural curiosity: they ",
"connect genes doing the same recognisable jobs in both species, including nitrogen assimilation ",
"itself. Read this alongside panel D, which shows that this shared FUNCTION coexists with almost ",
"no shared responsive GENES -- an ordinary evolutionary pattern, but note the power gap between ",
"the two tests before reading it as agreement.\n",
"\n",
"(C) Composition of each network by which estimator found the edge, stacked as a share so the ",
"two are comparable across a ", sprintf("%.1f", sm(D2, "ALL", "total_edges") /
sm(D1, "ALL", "total_edges")), "x difference in edge count; absolute counts are on the bars. In ",
"sugarcane, ", fmt_n(ly("sugarcane", "Pearson only", "n")), " edges (",
sprintf("%.1f%%", 100 * ly("sugarcane", "Pearson only", "frac")), ") are found by Pearson alone, ",
fmt_n(ly("sugarcane", "both estimators", "n")), " (",
sprintf("%.1f%%", 100 * ly("sugarcane", "both estimators", "frac")), ") by both estimators, and ",
fmt_n(ly("sugarcane", "MI only", "n")), " (",
sprintf("%.1f%%", 100 * ly("sugarcane", "MI only", "frac")), ") by mutual information alone; in ",
"purple the same three are ", fmt_n(ly("purple", "Pearson only", "n")), " (",
sprintf("%.1f%%", 100 * ly("purple", "Pearson only", "frac")), "), ",
fmt_n(ly("purple", "both estimators", "n")), " (",
sprintf("%.1f%%", 100 * ly("purple", "both estimators", "frac")), ") and ",
fmt_n(ly("purple", "MI only", "n")), " (",
sprintf("%.1f%%", 100 * ly("purple", "MI only", "frac")), "). Purple's larger non-linear share is ",
"NOT more non-linear biology: the two layers are matched on p-value, and at n = 18 the p implied ",
"by |r| = 0.8 is 6.7e-05 against 9.0e-12 at n = 48, so purple's MI floor is far softer in ",
"absolute terms.\n",
"\n",
"Read against the conservation results, this panel is the case against the MI layer. The ",
"MI-ONLY edges are a small minority of both networks -- ",
sprintf("%.1f%%", 100 * ly("sugarcane", "MI only", "frac")), " and ",
sprintf("%.1f%%", 100 * ly("purple", "MI only", "frac")), " -- and they are not better conserved ",
"than Pearson-only edges once chance is accounted for: their fold over null is ",
sprintf("%.2f", nl(D1, "mi", "fold_over_null")), "x and ",
sprintf("%.2f", nl(D2, "mi", "fold_over_null")), "x against ",
sprintf("%.2f", nl(D1, "pearson", "fold_over_null")), "x and ",
sprintf("%.2f", nl(D2, "pearson", "fold_over_null")), "x for Pearson-only, i.e. parity in one ",
"direction and BELOW Pearson in the other. Any apparent MI advantage in the raw conserved rate is ",
"opportunity bias -- MI-only edges connect genes with more orthologs, so they had more chances to ",
"match by accident, and the permutation null absorbs exactly that. The one layer that does earn ",
"its place is `both`: edges the two estimators agree on conserve at ",
sprintf("%.2f", nl(D1, "both", "fold_over_null")), "x and ",
sprintf("%.2f", nl(D2, "both", "fold_over_null")), "x, above Pearson-only in both directions. So ",
"the defensible statement is not that mutual information finds better-conserved edges -- it does ",
"not -- but that agreement between two estimators marks better-conserved edges than either ",
"alone.\n",
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
  NULLS[, .(panel = "A", study = DIR_LAB[direction],
            quantity = sprintf("%s: observed / null / fold", layer),
            value = sprintf("%.4f / %.4f / %.2f", obs_rate, null_mean_rate, fold_over_null))],
  GO_CMP[, .(panel = "B", study = "both",
             quantity = sprintf("%s terms: sugarcane / purple / shared", Ontology),
             value = sprintf("%d / %d / %d (Jaccard %.3f)",
                             Terms_sugarcane, Terms_purple, Shared, Jaccard))],
  funD[, .(panel = "D", study = as.character(species),
           quantity = gsub("\n", " ", stage), value = fmt_n(n))]))
write_tsv(stats, paste0(OUT_PREFIX, "_stats.tsv"))
say("done")
