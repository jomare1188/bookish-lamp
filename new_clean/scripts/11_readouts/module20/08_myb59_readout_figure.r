#!/usr/bin/env Rscript
# =============================================================================
# 08_myb59_readout_figure.r -- the AtMYB59 copies as network objects, in one page.
#
# THE READOUT IS A FOUR-PART NEGATIVE, and the figure is built to show that rather
# than to rescue a positive:
#
#   A  they are NOT hubs. Munoz define Module 20 as high-betweenness on their own
#      network; in ours the copies are mid-degree and the purple focus gene sits
#      BELOW purple's median.
#   B  no sugarcane copy reproduces the purple gene's neighbourhood. Nothing clears
#      BH across the nine, and the 1:1 ortholog is specifically indistinguishable
#      from a degree-matched purple gene.
#   C  the copies DO share neighbourhoods -- with each other, in sugarcane, up to
#      83% of the smaller one. So the ~1-2% cross-species figure in B is not a
#      limitation of the measure.
#   D  and the neighbourhoods are functionally different: the purple gene's is
#      stress and specialised metabolism, the sugarcane ortholog's is core
#      housekeeping.
#
# Panels C and B share an axis definition on purpose -- the same overlap
# coefficient -- because the contrast between them IS the result.
#
# Nothing is written on the figure that belongs in the legend, and every number in
# the legend is read from the same tables the panels are drawn from.
#
# RUN: ./run_all.sh 08
# =============================================================================

suppressMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(scico); library(grid); library(scales)
})

BASE <- "/dados04/jorge/comparative_saccharum"
RES  <- file.path(BASE, "new_clean/results")
OUT  <- file.path(RES, "readouts/module20")
PREFIX <- file.path(OUT, "myb59_copies_readout")
FOCUS <- "Soffic.09G0001580-9H"
say <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
short <- function(g) sub("^SoffiXsponR570\\.", "", g)

PAL <- setNames(scico(3, palette = "batlow")[1:2], c("sugarcane", "purple"))
theme_f <- theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(face = "bold", size = 7.2),
        legend.key.size = unit(3.5, "mm"),
        plot.margin = margin(3, 5, 3, 3))

OV  <- fread(file.path(OUT, "myb59_neighbour_overlap.tsv"))
WIN <- fread(file.path(OUT, "myb59_within_sugarcane_overlap.tsv"))
GO  <- fread(file.path(OUT, "myb59_neighbour_go.tsv"))
TAB <- fread(file.path(OUT, "module20_network_table.tsv"))[at_anchor == "AT5G59780"]

# =============================================================================
# A — degree percentile: they are not hubs
# =============================================================================
degs <- rbindlist(lapply(c("sugarcane", "purple"), function(sp)
  fread(file.path(RES, sp, sprintf("network_%s_node_metrics.tsv", sp)),
        select = c("gene", "degree"))[, study := sp]))
tgt <- unique(TAB[, .(gene, study = species)])
A <- merge(tgt, degs, by = c("gene", "study"))
A[, pct := {
  d <- degs[study == .BY$study, degree]
  100 * vapply(degree, function(x) mean(d < x), 0)
}, by = study]
A[, study := factor(study, levels = c("sugarcane", "purple"))]
A[, lab := short(gene)]
setorder(A, study, pct)
A[, lab := factor(lab, levels = unique(lab))]
say("panel A: degree percentile")
print(A[, .(study, lab, degree, pct = round(pct, 1))], row.names = FALSE)

pA <- ggplot(A, aes(pct, lab, colour = study)) +
  geom_vline(xintercept = 50, linetype = "22", linewidth = 0.35, colour = "grey45") +
  geom_segment(aes(x = 50, xend = pct, y = lab, yend = lab), linewidth = 0.4) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = PAL, name = NULL) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100),
                     labels = function(x) paste0(x, "%")) +
  labs(x = "degree percentile in its own network", y = NULL) +
  theme_f + theme(axis.text.y = element_text(size = 5.6),
                  legend.position = "bottom",
                  legend.box.margin = margin(-7, 0, 0, 0))

# =============================================================================
# B — cross-species neighbourhood overlap, against a degree-matched null
# =============================================================================
B <- copy(OV)
B[, lab := short(sugarcane_gene)]
setorder(B, frac_of_focus)
B[, lab := factor(lab, levels = lab)]
B[, orth := fifelse(is_ortholog_of_focus %in% c(TRUE, "TRUE"), "1:1 ortholog", "other copy")]
say("panel B: cross-species overlap")
print(B[, .(lab, shared, shared_expected, frac = round(100 * frac_of_focus, 2),
            p_adj)], row.names = FALSE)

pB <- ggplot(B, aes(100 * frac_of_focus, lab)) +
  geom_segment(aes(x = 100 * null_mean_frac, xend = 100 * frac_of_focus,
                   y = lab, yend = lab), linewidth = 0.4, colour = "grey70") +
  geom_point(aes(x = 100 * null_mean_frac), size = 1.5, shape = 4,
             colour = "grey45", stroke = 0.5) +
  geom_point(aes(shape = orth), size = 2, colour = unname(PAL["sugarcane"])) +
  # nudged past the null marker, which sits left of the observed point on the one
  # row where the observed value is zero
  geom_text(aes(label = sprintf("%d shared", shared),
                x = pmax(100 * frac_of_focus, 100 * null_mean_frac)),
            hjust = -0.28, size = 2.0, colour = "grey30") +
  scale_shape_manual(values = c(`1:1 ortholog` = 17, `other copy` = 16), name = NULL) +
  scale_x_continuous(limits = c(0, 2.9), labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0.02, 0.10))) +
  labs(x = sprintf("%% of %s's %s neighbours recovered\n(x = degree-matched null)",
                   "the purple gene", format(B$focus_neighbours[1], big.mark = ",")),
       y = NULL) +
  theme_f + theme(axis.text.y = element_text(size = 5.6),
                  axis.title.x = element_text(size = 6.2, lineheight = 1.05),
                  legend.position = "bottom",
                  legend.box.margin = margin(-7, 0, 0, 0))

# =============================================================================
# C — the same measure WITHIN sugarcane: the copies do share neighbourhoods
# =============================================================================
# Same statistic as B. If the cross-species figure were a limitation of the measure
# rather than a result, this panel would be low too. It is not.
C <- rbindlist(list(
  WIN[, .(comparison = "within sugarcane\n(copy vs copy)", v = 100 * overlap_coef)],
  OV[, .(comparison = "across species\n(copy vs purple gene)", v = 100 * frac_of_focus)]))
C[, comparison := factor(comparison, levels = unique(comparison))]
say("panel C: within-sugarcane vs cross-species overlap")
print(C[, .(median = round(median(v), 2), max = round(max(v), 2), n = .N),
        by = comparison], row.names = FALSE)

pC <- ggplot(C, aes(comparison, v, colour = comparison)) +
  geom_boxplot(outlier.shape = NA, width = 0.45, linewidth = 0.35, colour = "grey50") +
  geom_jitter(width = 0.12, height = 0, size = 1.3, alpha = 0.8) +
  scale_colour_manual(values = c("grey25", unname(PAL["purple"])), guide = "none") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(x = NULL, y = "shared neighbours\n(% of the smaller set)") +
  theme_f + theme(axis.text.x = element_text(size = 6.0, lineheight = 1.05),
                  axis.title.y = element_text(size = 6.4, lineheight = 1.05))

# =============================================================================
# D — the neighbourhoods do different things
# =============================================================================
ORTH <- OV[is_ortholog_of_focus %in% c(TRUE, "TRUE"), sugarcane_gene][1]
NTERM <- 7
pick <- function(g, sp) {
  d <- GO[gene == g & ontology %in% c("BP", "MF")][order(pvalue)][seq_len(min(NTERM, .N))]
  if (!nrow(d)) return(NULL)
  d[, .(who = sp, Term, ontology, mlp = -log10(pmax(pvalue, 1e-300)))]
}
# Short strip labels: the full gene ids are wider than a half-panel facet strip and
# get clipped. They are named in the legend and in panel B instead.
D <- rbindlist(list(pick(FOCUS, "purple gene"),
                    pick(ORTH, "sugarcane 1:1 ortholog")),
               fill = TRUE)
D[, who := factor(who, levels = unique(who))]
D[, Term_w := vapply(Term, function(t) paste(strwrap(t, 34), collapse = "\n"), "")]
D[, key := paste(who, Term, sep = "|")]
setorder(D, who, mlp)
D[, key := factor(key, levels = key)]
say("panel D: top neighbourhood terms")
print(D[, .(who = substr(as.character(who), 1, 22), ontology,
            Term = substr(Term, 1, 40), mlp = round(mlp, 2))], row.names = FALSE)

pD <- ggplot(D, aes(mlp, key, fill = ontology)) +
  geom_col(width = 0.62) +
  facet_wrap(~ who, ncol = 2, scales = "free") +
  scale_fill_manual(values = setNames(scico(3, palette = "lapaz",
                                            begin = 0.15, end = 0.75)[1:2],
                                      c("BP", "MF")), name = NULL) +
  scale_y_discrete(labels = setNames(D$Term_w, D$key)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = expression(-log[10](p)~", topGO weight01"), y = NULL) +
  theme_f + theme(axis.text.y = element_text(size = 5.2, lineheight = 0.92),
                  legend.position = "bottom",
                  legend.box.margin = margin(-7, 0, 0, 0))

# =============================================================================
fig <- (pA | pB) / (pC | pD) +
  plot_layout(heights = c(1, 1.15), widths = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12))

W <- 22; H <- 17
dir.create(dirname(PREFIX), showWarnings = FALSE, recursive = TRUE)
ggsave(paste0(PREFIX, ".png"), fig, width = W, height = H, units = "cm", dpi = 400)
ggsave(paste0(PREFIX, ".pdf"), fig, width = W, height = H, units = "cm")
say("wrote ", basename(PREFIX), ".{png,pdf}  (", W, " x ", H, " cm)")
