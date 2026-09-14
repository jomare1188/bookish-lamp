# =============================================================================
# 23_fig_conservation.r — the conservation figure: what transfers between the
# two networks, and what does not.
#
# This is the comparative core of the project. Two networks built from two
# unrelated studies are joined through the OrthoFinder orthogroups and asked
# whether their EDGES agree more than orthology alone would produce. The figure
# has to carry findings that pull in opposite directions and are all true: edge
# conservation is real and ~1.7x above chance in both directions, it rises with
# edge strength, and the nitrogen RESPONSE still barely travels with it.
#
# REBUILT 2026-09-13 ON THE PEARSON-ONLY GRAPHS. Every number now comes from
# 62_conserved_edges_pearson.py and 61_conserved_blocked_nodes.r. The old panel C
# (edge composition by layer) is gone: there is one layer, so it was a 100% bar.
# What replaced it is the finding that took its place in the analysis too.
#
#   A  OBSERVED CONSERVED-EDGE RATE AGAINST ITS PERMUTATION NULL, per direction.
#      The raw rates differ ~7x between directions and that difference is
#      MEANINGLESS on its own -- purple has 9x more edges, so a sugarcane edge has
#      far more chances to find a partner than the reverse. The panel shows the
#      observed rate beside the null built for that direction, which is the only
#      way the two are readable together.
#
#   B  WHAT THE CONSERVED GENES ARE FOR. One topGO BP enrichment per species over
#      the genes on at least one conserved edge, against that network's own
#      GO-annotated nodes as background, showing the terms enriched in BOTH.
#      Conservation of edges is a structural statement; this is the panel that
#      says whether the structure that transfers is doing anything recognisable.
#      Terms are ranked by their WORSE p-value across the two species, so the
#      panel shows terms both species agree on rather than terms one species
#      drives.
#
#   C  CONSERVATION AGAINST EDGE STRENGTH. Rank-based weight deciles, one line per
#      direction, on free y scales because the two rates differ ~7x for the reason
#      panel A exists. It is monotone in BOTH species and -- the load-bearing part
#      -- so is the fold over null, which is why the fold is drawn as the second
#      series rather than left in a table. If strong edges merely joined
#      better-annotated genes the rate would rise and the fold would not.
#
#   D  THE FUNNEL, per species: network nodes -> genes on at least one conserved
#      edge -> nitrogen-responsive -> responsive on BOTH sides of an ortholog
#      pair. Log scale, because it spans four orders of magnitude. This is where
#      the comparative question narrows hardest.
#
# NOTHING IS WRITTEN ON THE FIGURE THAT BELONGS IN THE LEGEND -- panel letters and
# the labels the data needs, nothing else. See docs/methods.md, "The rule every
# paper figure follows".
#
# NO EDGE TABLE IS READ. The conserved-edge files are 472 MB and 517 MB; every
# number here comes from the summary tables beside them.
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
# ONE table per direction now: 62 writes the observed rate, the null and the
# decile breakdown together, so the figure cannot pair a rate with the wrong null.
SUMS <- rbindlist(lapply(DIRS, function(d)
  fread(file.path(CONS_DIR, sprintf("conservation_summary_%s_pearson.tsv", d)))))
STAT <- rbindlist(lapply(c("sugarcane", "purple"), function(sp)
  fread(file.path(CONS_DIR, sprintf("%s_status_blocked_nodes.tsv", sp)))[
    , .(species = sp, status)]))

# BOTH SIDES OF THE PARTITION. `conserved` is the genes on at least one conserved
# edge; `nonconserved` is the exact complement -- network nodes with NO conserved
# edge (37,867 + 64,123 = 101,990 in sugarcane). Same background for both, so the
# two are a contrast rather than two unrelated tests.
SETS <- c(conserved = "on a conserved edge", nonconserved = "no conserved edge")
godir <- function(tag) file.path(CONS_DIR, sprintf("enrichment_%s", tag))
gof <- function(tag, pat) file.path(godir(tag), sprintf(pat, GO_ONT, tag))
GO_DIR <- godir("conserved")
GO_SHARED <- rbindlist(lapply(names(SETS), function(tg)
  fread(gof(tg, "GO_%s_%s_shared_terms.csv"))[, set := tg]))
# topGO's GenTable truncates term names at 40 characters with an ellipsis, so the
# CSV carries strings like "positive regulation of cellular cataboli...". Expand
# them from the cached GO.db dump; anything missing keeps the truncated string
# rather than becoming NA.
GO_NAMES <- env_opt("CLEAN_GO_NAMES")
if (nzchar(GO_NAMES) && file.exists(GO_NAMES)) {
  nm <- fread(GO_NAMES)
  GO_SHARED <- merge(GO_SHARED, nm, by = "GO.ID", all.x = TRUE, sort = FALSE, allow.cartesian = TRUE)
  n_trunc <- GO_SHARED[grepl("\\.\\.\\.$", Term), .N]
  GO_SHARED[!is.na(Term_full) & nzchar(Term_full), Term := Term_full]
  say(sprintf("expanded %d truncated %s term name(s) from %s",
              n_trunc, GO_ONT, basename(GO_NAMES)))
} else {
  say("NOTE: no GO term-name cache; long term names stay truncated as topGO wrote them")
}
GO_CMP <- rbindlist(lapply(names(SETS), function(tg)
  rbindlist(lapply(c("BP", "MF", "CC"), function(o)
    fread(file.path(godir(tg), sprintf("GO_%s_%s_comparison_summary.csv", o, tg)))[
      , set := tg]))))
GO_SUM <- fread(gof("conserved", "GO_%s_%s_summary.csv"))
gc_ <- function(o, col, tg = "conserved") GO_CMP[Ontology == o & set == tg][[col]]
gs_ <- function(sp, col) GO_SUM[Network == sp][[col]]

SUMS[, dir_lab := factor(DIR_LAB[direction], levels = unname(DIR_LAB))]
ALLROW <- SUMS[stratum == "all" & weight_bin == "ALL"]
sm <- function(d, col) ALLROW[direction == d][[col]]

say("conservation, Pearson-only graphs:")
print(ALLROW[, .(dir_lab, edges = fmt_n(edges), conserved = fmt_n(conserved),
                 pct = round(100 * conserved_fraction, 3),
                 fold = fold_over_null, p = p_emp)], row.names = FALSE)

# =============================================================================
# A — observed vs null, per direction (ALL edges)
# =============================================================================
# The observed rate comes from the FULL edge set; the null from whatever set 62
# could afford for that stratum (`null_basis`). For `all` that is a 5M sample, so
# the null bar is the count its rate implies over the full set -- edges expected to
# find a partner by chance. Pairing a full-pass rate with a sampled null would
# confound the two, which is why 62 carries the basis and this panel prints it.
obsA <- rbindlist(list(
  ALLROW[, .(dir_lab, edges, what = "observed", rate = conserved_fraction,
             n = conserved)],
  ALLROW[, .(dir_lab, edges, what = "permutation null", rate = null_mean,
             n = round(null_mean * edges))]))
obsA[, what := factor(what, levels = names(PAL_OBS))]
# The denominator rides on the axis label, so the panel carries the percentage,
# the count and what the count is out of. Broken across three short lines rather
# than two long ones: at this panel width "sugarcane → purple" on one line is
# wider than the space a single group gets, and the two group labels collide.
ax_lab <- function(d, n)
  sprintf("%s\n%s edges", sub(" (→|->) ", "\n\\1 ", d), trimws(fmt_n(n)))
obsA[, dir_ax := factor(ax_lab(dir_lab, edges),
                        levels = unique(ax_lab(dir_lab, edges)[order(dir_lab)]))]
say("panel A: observed vs null, with absolute counts")
print(obsA[, .(dir_lab, what, pct = round(100 * rate, 3), n = fmt_n(n))], row.names = FALSE)

pA <- ggplot(obsA, aes(dir_ax, rate, fill = what)) +
  geom_col(position = position_dodge(width = 0.65), width = 0.55) +
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
# WHAT THE PANEL CONTRASTS. Top terms of the genes ON a conserved edge against top
# terms of the exact complement -- genes with NO conserved edge. The two sets
# partition the network and share a background, so the comparison is of the same
# genes split one way, not of two unrelated tests.
#
# Each set gets its OWN top terms rather than the two being scored on one term list:
# the finding is that the lists differ, and forcing a common list would hide exactly
# that. Within each set, terms are ranked by the WORSE of the two species' p-values,
# so these are terms BOTH species agree on rather than terms one species drives.
gsh <- copy(GO_SHARED)
setnames(gsh, c("pvalue_sugarcane", "pvalue_purple"), c("p_sugarcane", "p_purple"))
gsh[, worst := pmax(p_sugarcane, p_purple)]
setorder(gsh, set, worst)
NPER <- max(4L, GO_NTERMS %/% 2L)
top <- gsh[, head(.SD, NPER), by = set]
top[, set_lab := factor(SETS[set], levels = unname(SETS))]

# Full GO names run to 70+ characters and, unwrapped, the label column eats half the
# figure width. Wrapped to two short lines. The set tag is carried in the factor
# level because two sets can share a term name and would otherwise collapse onto one
# row of the shared y axis.
wrap_term <- function(x, w = 34)
  vapply(x, function(t) paste(strwrap(t, w), collapse = "\n"), "")
top[, Term_w := wrap_term(Term)]
top[, key := paste(set, Term, sep = "|")]
setorder(top, set_lab, -worst)
top[, key := factor(key, levels = key)]

goB <- melt(top[, .(key, set_lab, Term_w, p_sugarcane, p_purple)],
            id.vars = c("key", "set_lab", "Term_w"),
            variable.name = "species", value.name = "p")
goB[, species := factor(sub("^p_", "", species), levels = c("sugarcane", "purple"))]
goB[, mlp := -log10(p)]
say(sprintf("panel B: %d shared %s terms (conserved) and %d (non-conserved); top %d each",
            gsh[set == "conserved", .N], GO_ONT, gsh[set == "nonconserved", .N], NPER))
print(top[, .(set, Term = substr(Term, 1, 44), p_sugarcane, p_purple)], row.names = FALSE)

pB <- ggplot(goB, aes(mlp, key)) +
  geom_line(aes(group = key), colour = "grey70", linewidth = 0.35) +
  geom_point(aes(fill = species), shape = 21, size = 1.9,
             colour = "grey25", stroke = 0.25) +
  facet_grid(set_lab ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_y_discrete(labels = setNames(top$Term_w, top$key)) +
  scale_fill_manual(values = setNames(scico(3, palette = "batlow")[1:2],
                                      c("sugarcane", "purple")), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.08))) +
  labs(x = expression(-log[10](p)~", topGO weight01"), y = NULL) +
  theme_f +
  theme(axis.text.y = element_text(size = 5.7, lineheight = 0.92),
        panel.grid.major.y = element_line(linewidth = 0.2),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 90, size = 6.3, face = "bold"),
        panel.spacing.y = unit(1.6, "mm"))

# =============================================================================
# C — conservation against edge strength
# =============================================================================
# WHAT THIS REPLACED. The old panel C was edge composition by layer, the sharpest
# form of the "MI adds nothing" argument. The networks are Pearson-only now, so
# that panel is one 100% bar. This is what took its place in the analysis as well:
# the question a single weighted layer can still be asked is whether the STRENGTH
# of co-expression predicts whether an edge survives across species.
#
# BOTH SERIES ARE DRAWN. The rate alone would be consistent with strong edges
# simply joining better-annotated genes; the fold over null rising with it is what
# rules that out, so it is a panel series rather than a number in the legend.
DEC <- SUMS[stratum == "all" & weight_bin != "ALL"]
DEC[, decile := as.integer(sub("^D", "", weight_bin))]
setorder(DEC, direction, decile)
say("panel C: conservation by weight decile")
print(DEC[, .(dir_lab, decile, edges = fmt_n(edges),
              pct = round(100 * conserved_fraction, 3), fold = fold_over_null)],
      row.names = FALSE)

# WHAT THIS PANEL SAYS, in one sentence: the stronger the co-expression, the likelier
# the edge exists in the other species too -- in both directions.
#
# ONE PLOT, NOT TWO. An earlier version drew the fold over null beside the rate. The
# fold is still computed, still in the summary table and still quoted in the legend --
# it is what rules out the rate rising merely because strongly co-expressed genes have
# more or better orthologs -- but it is a second statistic about the same curve, and
# the panel reads better making one point than two.
#
# THE AXIS IS LOGARITHMIC, and that is not only to fit both directions on it. Sugarcane
# runs at ~10% and purple at ~1.5% -- the 7x that panel A exists to explain -- so a
# linear axis flattens purple onto the floor. On a log axis both fit AND slopes read as
# RELATIVE change, which is the comparison that means something here: purple rises 84%
# across the deciles against sugarcane's 16%, so the effect is the LARGER one in purple
# despite the lower absolute rate. A linear panel would show the opposite.
say("panel C: conservation by weight decile")
print(DEC[, .(dir_lab, decile, edges = fmt_n(edges),
              pct = round(100 * conserved_fraction, 3), fold = fold_over_null)],
      row.names = FALSE)

pC <- ggplot(DEC, aes(decile, conserved_fraction, colour = dir_lab)) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 1.4) +
  scale_colour_manual(values = PAL_DIR, name = NULL) +
  scale_y_log10(labels = function(x) sprintf("%g%%", 100 * x),
                breaks = c(0.0125, 0.025, 0.05, 0.10),
                expand = expansion(mult = c(0.10, 0.10))) +
  scale_x_continuous(breaks = c(1, 4, 7, 10),
                     labels = c("D1\nweakest", "D4", "D7", "D10\nstrongest"),
                     expand = expansion(mult = 0.06)) +
  labs(x = "edge weight decile",
       y = "edges with a conserved partner\n(log scale)") +
  theme_f +
  theme(axis.title.y = element_text(size = 6.5, lineheight = 1.05),
        axis.title.x = element_text(size = 6.7),
        axis.text.x = element_text(size = 6.0, lineheight = 1.0),
        legend.position = "bottom",
        legend.box.margin = margin(-7, 0, 0, 0),
        legend.text = element_text(size = 6.5))

# =============================================================================
# D — the funnel
# =============================================================================
n_nodes <- c(sugarcane = nrow(fread(NODES_SUG, select = 1L)),
             purple    = nrow(fread(NODES_PUR, select = 1L)))
# Stages 3 and 4 come from 61_conserved_blocked_nodes.r's status tables, which ARE
# the responsive genes (one row each) with `conserved_correlated` marking the ones
# whose ortholog responds too. Reading them rather than a separate summary is what
# keeps the funnel and the node-level analysis from drifting apart.
# The conserved-gene count lives in the gene list 62 writes, not in the summary
# table -- read the file rather than a number that would have to be kept in sync.
n_cons <- vapply(c("sugarcane", "purple"), function(sp)
  length(readLines(file.path(CONS_DIR, sprintf("conserved_genes_%s_pearson.txt", sp)),
                   warn = FALSE)), integer(1))
n_resp <- STAT[, .N, by = species]
n_both <- STAT[status == "conserved_correlated", .N, by = species]
gv <- function(d, sp) d[species == sp, N]
funD <- data.table(
  species = rep(c("sugarcane", "purple"), each = 4),
  stage = factor(rep(c("nodes in\nthe network", "on a conserved\nedge",
                       "nitrogen-\nresponsive", "responsive on\nboth sides"), 2),
                 levels = c("nodes in\nthe network", "on a conserved\nedge",
                            "nitrogen-\nresponsive", "responsive on\nboth sides")),
  n = c(n_nodes[["sugarcane"]], n_cons[["sugarcane"]],
        gv(n_resp, "sugarcane"), gv(n_both, "sugarcane"),
        n_nodes[["purple"]], n_cons[["purple"]],
        gv(n_resp, "purple"), gv(n_both, "purple")))
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

D1 <- DIRS[1]; D2 <- DIRS[2]
st <- function(d, stratum, col) SUMS[direction == d & stratum == ..stratum &
                                     weight_bin == "ALL"][[col]]
dc <- function(d, dec, col) DEC[direction == d & decile == dec][[col]]
pc <- function(x) sprintf("%.2f%%", 100 * x)
gsum <- function(sp, col) gs_(sp, col)
fn <- function(sp) funD[species == sp, n]

legend <- paste0(
"Figure ", FIG, ". Cross-species conservation of network edges on the unpruned Pearson-only ",
"graphs, and how little of the nitrogen response travels with it. An edge is CONSERVED when both ",
"of its genes have orthologs in the other species and those orthologs are themselves joined by an ",
"edge there. Orthology is many-to-many and any ortholog pair counts, so this is the most ",
"permissive definition available; it is scored against a null that removes only the orthology ",
"ASSIGNMENT.\n",
"\n",
"(A) Observed conserved-edge rate against its own permutation null, per direction; bar labels are ",
"absolute edge counts. In ", fmt_n(ALLROW[direction == D1, edges]), " sugarcane edges, ",
fmt_n(ALLROW[direction == D1, conserved]), " (", pc(ALLROW[direction == D1, conserved_fraction]),
") have a purple counterpart, ", ALLROW[direction == D1, fold_over_null],
"x the null rate; in ", fmt_n(ALLROW[direction == D2, edges]), " purple edges, ",
fmt_n(ALLROW[direction == D2, conserved]), " (", pc(ALLROW[direction == D2, conserved_fraction]),
"), ", ALLROW[direction == D2, fold_over_null], "x. THE TWO DIRECTIONS ARE NOT COMPARABLE TO EACH ",
"OTHER and the panel is not an invitation to compare them: purple's edge density is 3.2x ",
"sugarcane's, so a random ortholog pair is far likelier to land on an edge there, which is exactly ",
"what each direction's own null absorbs. The null permutes the target column of the ortholog ",
"table, preserving every gene's ortholog count, which genes have any ortholog at all, and the ",
"multiset of targets -- so the only thing destroyed is which target each source maps to. It is ",
"computed over ortholog pairs whose BOTH sides are network nodes (114,681 pairs over 43,390 ",
"orthogroups, the same universe the node-level analysis uses), which is a stricter null than the ",
"merged-network version this replaces: that one drew from all 275,047 pairs including genes that ",
"could never have matched, and reported a correspondingly larger fold.\n",
"\n",
"(B) WHAT TRANSFERS IS HOUSEKEEPING; WHAT DOES NOT IS REGULATION. One topGO enrichment (",
GO_ONT, ", weight01) per species over the genes ON at least one conserved edge (top) and over the ",
"EXACT COMPLEMENT -- network nodes with no conserved edge (bottom). The two sets partition the ",
"node set (sugarcane ", fmt_n(fn("sugarcane")[2]), " + ",
fmt_n(fn("sugarcane")[1] - fn("sugarcane")[2]), " = ", fmt_n(fn("sugarcane")[1]),
") and are tested against the SAME background, that network's own GO-annotated nodes, so this is ",
"one gene set split two ways rather than two unrelated tests. Points are the two species' ",
"p-values for the same term, joined by a line; within each set terms are ranked by the WORSE of ",
"the two and the top ", NPER, " are drawn, so these are terms BOTH species agree on rather than ",
"terms one species drives. Each set gets its own term list because the finding IS that the lists ",
"differ -- conserved edges are enriched for translation, protein folding and transport, ",
"photosynthesis and splicing, the core cellular machinery; the complement is enriched for ",
"regulation of transcription, auxin and ethylene signalling, ubiquitination and the cell cycle, ",
"i.e. the regulatory layer. ", fmt_n(gc_(GO_ONT, "Shared", "conserved")), " ", GO_ONT,
" terms are shared by both species in the conserved set (Jaccard ",
sprintf("%.3f", gc_(GO_ONT, "Jaccard", "conserved")), ") and ",
fmt_n(gc_(GO_ONT, "Shared", "nonconserved")), " in the complement (Jaccard ",
sprintf("%.3f", gc_(GO_ONT, "Jaccard", "nonconserved")), "). ",
"TWO THINGS THE SPLIT IS ENTANGLED WITH, both measured. FIRST, ORTHOLOGY: a gene with no ",
"ortholog at all cannot have a conserved edge, so it is in the bottom set by construction -- ",
"2,775 sugarcane genes and 1,217 purple ones among the nitrogen-responsive alone -- which makes ",
"the regulatory signal partly a statement about which gene families have clean one-to-one ",
"orthology across these two genomes. SECOND, DEGREE, and this one is partly arithmetic: genes on ",
"a conserved edge are hubs (median degree 228 against 18 in sugarcane, 6,916 against 376 in ",
"purple; 63.3% of sugarcane's top-degree decile against 7.1% of its bottom), because one ",
"conserved edge out of 6,677 is near-certain at a 10.4% per-edge rate and unlikely out of 2. ",
"Ranking the SAME genes by degree instead of by conservation reproduces almost this same term ",
"split -- figure 3 panel C -- so the two are one contrast seen twice and neither is independent ",
"evidence for the other. GO comes from the adopted ",
"full-InterProScan annotation, which reaches 62.9% and 64.4% of network nodes. topGO reports ",
"anything below 1e-30 AS 1e-30, so a term at the right edge is censored rather than exactly that ",
"significant, and two points there coincide because both are at the floor -- not because the two ",
"species agree to that precision.\n",
"\n",
"(C) THE STRONGER THE CO-EXPRESSION, THE LIKELIER THE EDGE EXISTS IN THE OTHER SPECIES -- in ",
"both directions. Edges are split into ten rank-based weight deciles, so each bin holds a tenth ",
"of that network's edges, D1 the weakest and D10 the strongest, and the panel plots the share of ",
"each decile that has a conserved partner: sugarcane from ", pc(dc(D1, 1, "conserved_fraction")),
" to ", pc(dc(D1, 10, "conserved_fraction")), " and purple from ",
pc(dc(D2, 1, "conserved_fraction")), " to ", pc(dc(D2, 10, "conserved_fraction")),
", monotonically in both. THE AXIS IS LOGARITHMIC and that is not only to fit both curves on it: ",
"the two directions differ ~7x in absolute rate for the reason panel A gives, so a linear axis ",
"lays purple flat on the floor, and on a log axis slopes read as RELATIVE change, which is the ",
"comparison that means anything here. Read that way the result inverts -- purple's conservation ",
"rises ", sprintf("%.0f%%", 100 * (dc(D2, 10, "conserved_fraction") / dc(D2, 1, "conserved_fraction") - 1)),
" across the deciles against sugarcane's ",
sprintf("%.0f%%", 100 * (dc(D1, 10, "conserved_fraction") / dc(D1, 1, "conserved_fraction") - 1)),
", so the effect is the LARGER one in purple despite its far lower absolute rate. ",
"A RISING RATE ALONE WOULD NOT SETTLE THIS, because it is also what you would see if strongly ",
"co-expressed genes simply had more or better orthologs. What rules that out is the fold over ",
"each decile's own permutation null, which rises with it: ", dc(D1, 1, "fold_over_null"), " to ",
dc(D1, 10, "fold_over_null"), " in sugarcane and ", dc(D2, 1, "fold_over_null"), " to ",
dc(D2, 10, "fold_over_null"), " in purple, never below 1. Those folds are in ",
basename(paste0(OUT_PREFIX, "_stats.tsv")), " and in the summary tables rather than on the panel, ",
"which makes one point. Weight is the per-study rescaling of |r| to [0.01, 1], so a decile does ",
"NOT stand for the same |r| in both species; the boundaries are in the summary table. What this ",
"does not settle is WHY -- stronger selective constraint on tight co-expression, and a noisier ",
"weakest decile sitting nearest the |r| >= 0.8 threshold, predict the same shape.\n",
"\n",
"(D) The funnel, per species, log scale. Sugarcane: ", fmt_n(fn("sugarcane")[1]),
" network nodes -> ", fmt_n(fn("sugarcane")[2]), " on at least one conserved edge -> ",
fmt_n(fn("sugarcane")[3]), " nitrogen-responsive under the blocked test -> ",
fmt_n(fn("sugarcane")[4]), " whose ortholog is responsive in purple too. Purple: ",
fmt_n(fn("purple")[1]), " -> ", fmt_n(fn("purple")[2]), " -> ", fmt_n(fn("purple")[3]), " -> ",
fmt_n(fn("purple")[4]), ". The funnel does NOT close at orthology or at edge conservation -- ",
"37% and 25% of nodes sit on a conserved edge -- but at the last step, and the last step is small ",
"in both species.\n",
"\n",
"WHAT THE NITROGEN-RESPONSIVE EDGES DO, since it is the question the figure is built around and ",
"the raw numbers mislead. Edges joining two genes responsive in BOTH species are conserved at ",
pc(SUMS[direction == D1 & stratum == "resp_both" & weight_bin == "ALL", conserved_fraction]),
" in sugarcane (", fmt_n(SUMS[direction == D1 & stratum == "resp_both" & weight_bin == "ALL", edges]),
" edges) and ",
pc(SUMS[direction == D2 & stratum == "resp_both" & weight_bin == "ALL", conserved_fraction]),
" in purple (", fmt_n(SUMS[direction == D2 & stratum == "resp_both" & weight_bin == "ALL", edges]),
" edges), against backgrounds of ", pc(ALLROW[direction == D1, conserved_fraction]), " and ",
pc(ALLROW[direction == D2, conserved_fraction]), " -- a 3-4x raw enrichment. Against the ",
"ortholog-shuffle null MOST OF THAT IS ORTHOLOGY, because genes responsive in both species are by ",
"construction genes with good orthologs: sugarcane keeps a modest real excess (",
SUMS[direction == D1 & stratum == "resp_both" & weight_bin == "ALL", fold_over_null], "x, p = ",
SUMS[direction == D1 & stratum == "resp_both" & weight_bin == "ALL", p_emp], ") and purple's ",
"disappears entirely (",
SUMS[direction == D2 & stratum == "resp_both" & weight_bin == "ALL", fold_over_null], "x, p = ",
SUMS[direction == D2 & stratum == "resp_both" & weight_bin == "ALL", p_emp],
"). A shared nitrogen response therefore predicts edge conservation in sugarcane and not in ",
"purple, and the raw rates should not be quoted without their nulls.\n",
"\n",
"Sources: 62_conserved_edges_pearson.py (edges and nulls, from the mcxdump edge stream), ",
"61_conserved_blocked_nodes.r (who is responsive), 09_go_enrichment.r (panel B). The edge join was ",
"verified against an independent route -- gene names resolved through a fresh Orthogroups.tsv ",
"parse and membership confirmed by a direct pass over the raw edge dump -- confirming 4,000 of ",
"4,000 conserved edges and rejecting 4,000 of 4,000 non-conserved ones.")

legend_file <- paste0(OUT_PREFIX, "_legend.txt")
writeLines(wrap_at(legend), legend_file)
say("wrote ", basename(legend_file))

stats <- rbindlist(list(
  SUMS[weight_bin == "ALL", .(panel = "A", study = DIR_LAB[direction],
       quantity = sprintf("%s: observed / null / fold / p  [null on %s]",
                          stratum, null_basis),
       value = sprintf("%.5f / %.5f / %.2f / %.4g",
                       conserved_fraction, null_mean, fold_over_null, p_emp))],
  DEC[, .(panel = "C", study = DIR_LAB[direction],
          quantity = sprintf("decile %d (w %s-%s): rate / fold",
                             decile, weight_lo, weight_hi),
          value = sprintf("%.5f / %.2f", conserved_fraction, fold_over_null))],
  GO_CMP[, .(panel = "B", study = "both",
             quantity = sprintf("%s terms: sugarcane / purple / shared", Ontology),
             value = sprintf("%d / %d / %d (Jaccard %.3f)",
                             Terms_sugarcane, Terms_purple, Shared, Jaccard))],
  funD[, .(panel = "D", study = as.character(species),
           quantity = gsub("\n", " ", stage), value = fmt_n(n))]))
write_tsv(stats, paste0(OUT_PREFIX, "_stats.tsv"))
say("done")
