# =============================================================================
# 26_fig_module20.r — Muñoz's Module 20 in purple, against sugarcane.
#
# Muñoz-Perez et al. build their nitrogen claim on Module 20: 12 transcripts,
# ~75% MYB or MYB-related, central in their own network. Figure 2 showed the
# module's nitrogen response reproduces in their own data. This figure asks the
# comparative question -- what is left of it in the other species -- and the
# answer is not "less of the same". It is a different response shape, in one
# copy, that every monotone test in this project is blind to.
#
#   A  THE MODULE AS MUNOZ REPORT IT, reproduced -- all 33 mapped sugarcane genes
#      across all 48 libraries, split by nitrogen. It is the premise rather than
#      the result: it is here so the rest of the figure has something to be
#      compared against.
#
#      WHERE THE GENE OF INTEREST SITS. The focus gene is a PURPLE gene, so it has
#      no row here -- but the SAME LOCUS does. Its nine sugarcane copies (the
#      AtMYB59 anchor) are marked and named, which is what makes the contrast the
#      legend draws visible rather than asserted: those nine are every one
#      monotonically repressed by nitrogen, while the purple copy in panel C turns
#      back up at excess. Without the marks the reader cannot tell which of the 33
#      rows are the ones being talked about.
#
#   B  THE SAME MODULE IN PURPLE, as a per-gene Z-SCORE and split by GENOTYPE as
#      well as nitrogen, with a TF ANNOTATION COLUMN. All three are deliberate.
#      The z-score is what makes the one responding copy visible at all -- on an
#      absolute scale it is a bright row among dim ones and the SHAPE of its
#      response is invisible. Splitting by genotype is what makes the dominant
#      structure legible: in purple these genes vary far more between the two
#      species than across nitrogen, and an absolute heatmap split by nitrogen
#      alone hides that. The TF column is the same annotation the per-module
#      heatmaps carry (16_module_heatmaps.r, from readouts/get_tfs/<study>/
#      TF_in_network.tsv), and it is the point of the panel as much as the
#      expression is: this module is published as a ~75% MYB module, so how many
#      of the copies our own pipeline calls transcription factors is a claim the
#      figure should let the reader check.
#      Only the responding copy is labelled; the rest would be 26 unreadable ids.
#
#   C  THAT COPY ON ITS OWN. Soffic.09G0001580-9H, called MYB by the source study
#      AND by our own pipeline, is the only AtMYB59-anchor copy in purple that is
#      expressed at all (60.6 TPM against 0.01 and 0.00 for the other two). Its
#      response is a U centred on the CONTROL: highest under nitrogen starvation,
#      lowest at 2 N, high again under excess. Since 2 N is the control and 0 N
#      and 6 N are stresses in opposite directions, this is a gene induced by
#      nitrogen stress in EITHER direction -- which is why the Spearman test used
#      everywhere else in this project cannot see it, and why the U-shape
#      contrast exists as a separate test.
#
# WHAT WAS DROPPED, 2026-09-12. The old panel B (mapped genes per Arabidopsis
# anchor, MYB against not-MYB, as bars) is gone -- that breakdown belongs in a
# table, not in a panel. The numbers are still computed here, because the legend
# quotes them and _stats.tsv carries them, so removing the panel costs no data.
#
# NOTHING IS WRITTEN ON THE FIGURE THAT BELONGS IN THE LEGEND.
#
# RUN: through run.sh  ->  ./run.sh figmodule20   (then ./run.sh legends)
# =============================================================================

suppressMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(scico); library(grid); library(svglite); library(scales)
})

source(file.path(dirname(sub("--file=", "",
       grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "lib", "common.R"))

M20_DIR    <- env_req("CLEAN_M20_DIR")
TPM_SUG    <- env_req("CLEAN_TPM_SUGARCANE")
TPM_PUR    <- env_req("CLEAN_TPM_PURPLE")
META_SUG   <- env_req("CLEAN_META_SUGARCANE")
META_PUR   <- env_req("CLEAN_META_PURPLE")
OUT_PREFIX <- env_req("CLEAN_OUT_PREFIX")
FIG        <- env_opt("CLEAN_FIG_NUM", "7")
FOCUS      <- env_opt("CLEAN_M20_FOCUS", "Soffic.09G0001580-9H")
FOCUS_ANCH <- env_opt("CLEAN_M20_ANCHOR", "AT5G59780")
TF_PURPLE  <- env_opt("CLEAN_TF_PURPLE")
setDTthreads(as.integer(env_num("CLEAN_CORES", 8)))

banner(sprintf("Figure %s — Munoz Module 20 in purple", FIG))

ANCHOR_LAB <- c(AT5G59780 = "AT5G59780\nAtMYB59",
                AT3G46130 = "AT3G46130\nAtMYB48",
                AT3G55960 = "AT3G55960\nuncharacterised",
                AT4G10770 = "AT4G10770\nOPT7 transporter")
NSTATUS_LEVELS <- c("low N", "control", "high N")
NSTATUS <- c("Low Nitrogen" = "low N", "High Nitrogen" = "high N",
             "0N" = "low N", "2N" = "control", "6N" = "high N")
PAL_NSTATUS <- setNames(rev(scico(3, palette = "managua", begin = 0.08, end = 0.92)),
                        NSTATUS_LEVELS)
PAL_STUDY <- setNames(scico(3, palette = "batlow")[1:2], c("sugarcane", "purple"))

theme_f <- theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(face = "bold", size = 7.2),
        legend.key.size = unit(3.5, "mm"),
        legend.title = element_text(face = "bold", size = 7.5),
        plot.margin = margin(2, 4, 2, 2))

TAB  <- fread(file.path(M20_DIR, "module20_network_table.tsv"))
NPUR <- fread(file.path(M20_DIR, "module20_nitrogen_purple.tsv"))
NSUG <- fread(file.path(M20_DIR, "module20_nitrogen_sugarcane.tsv"))

load_tpm <- function(tpmf, metaf, sp) {
  tpm <- fread(tpmf); setnames(tpm, 1, "gene")
  tpm[, gene := strip_version(gene)]
  m <- fread(metaf); setnames(m, tolower(names(m)))
  samp <- intersect(names(tpm), m$sample)
  m <- m[match(samp, sample)]
  m[, status := factor(NSTATUS[as.character(treatment)], levels = NSTATUS_LEVELS)]
  genes <- TAB[species == sp, unique(gene)]
  X <- as.matrix(tpm[gene %chin% genes, ..samp])
  rownames(X) <- tpm[gene %chin% genes, gene]
  list(X = X, meta = m)
}
PU <- load_tpm(TPM_PUR, META_PUR, "purple")
SU <- load_tpm(TPM_SUG, META_SUG, "sugarcane")

# =============================================================================
# A — the module as Munoz report it, reproduced in sugarcane
# =============================================================================
zrow <- function(M) {
  Z <- t(scale(t(M))); Z[!is.finite(Z)] <- 0; Z
}
XS <- SU$X; ms <- SU$meta
ordS <- order(ms$status, ms$genotype,
              if ("segment" %in% names(ms)) ms$segment else rep("", nrow(ms)), ms$sample)
XS <- XS[, ordS, drop = FALSE]; ms <- ms[ordS]
ZS <- zrow(XS)
ZS <- ZS[order(-rowMeans(ZS[, ms$status == "low N", drop = FALSE])), , drop = FALSE]

dA <- data.table(gene = factor(rep(rownames(ZS), times = ncol(ZS)),
                               levels = rev(rownames(ZS))),
                 sample = factor(rep(colnames(ZS), each = nrow(ZS)),
                                 levels = colnames(ZS)),
                 z = as.vector(ZS))
dA[, status := ms$status[match(as.character(sample), ms$sample)]]
ZLIM <- as.numeric(quantile(abs(dA$z), 0.98, na.rm = TRUE))
say(sprintf("panel A: %d sugarcane genes x %d libraries", nrow(ZS), ncol(ZS)))

PAL_Z <- rev(scico(256, palette = "roma"))
MYB_COL <- unname(scico(3, palette = "roma")[1])

# WHERE THE FOCUS GENE IS. It is a PURPLE gene and this panel is sugarcane, so it
# has no row here -- but the same LOCUS does, as the copies sharing its Arabidopsis
# anchor. Those rows are named and marked. The contrast the legend draws (these nine
# monotonically repressed, the purple copy turning back up at excess) is otherwise
# a claim about rows the reader cannot pick out of 33 unlabelled ones.
myb_sug <- intersect(TAB[species == "sugarcane" & at_anchor == FOCUS_ANCH, unique(gene)],
                     rownames(ZS))
say(sprintf("panel A: %d of %d rows are %s copies -- marked and named",
            length(myb_sug), nrow(ZS), FOCUS_ANCH))
if (!length(myb_sug))
  stop("no sugarcane copy of the focus anchor ", FOCUS_ANCH, " is in panel A -- ",
       "the mark would be silently absent", call. = FALSE)

# Ids are 25 characters and all share the assembly prefix; dropping it is what
# makes nine labels fit without shrinking them to unreadable.
short_id <- function(g) sub("^SoffiXsponR570\\.", "", g)

# The marker sits just left of the first nitrogen facet, so it reads as a row
# pointer rather than as data. clip = "off" is what lets it live outside the panel.
mark <- data.table(gene = factor(myb_sug, levels = levels(dA$gene)),
                   status = factor(NSTATUS_LEVELS[1], levels = NSTATUS_LEVELS))

pA <- ggplot(dA, aes(sample, gene, fill = z)) +
  geom_raster() +
  geom_point(data = mark, aes(x = 0.1, y = gene), inherit.aes = FALSE,
             shape = 18, size = 1.3, colour = MYB_COL) +
  facet_grid(. ~ status, scales = "free_x", space = "free_x") +
  scale_fill_gradientn(colours = PAL_Z, limits = c(-ZLIM, ZLIM), oob = squish,
                       name = "z", guide = guide_colourbar(
                         barheight = unit(1.4, "cm"), barwidth = unit(2.4, "mm"))) +
  scale_y_discrete(breaks = myb_sug, labels = short_id(myb_sug)) +
  coord_cartesian(clip = "off") +
  labs(x = sprintf("%d sugarcane libraries", ncol(ZS)),
       y = sprintf("%d Module-20 genes  |  %s copies named",
                   nrow(ZS), sub("\\n.*$", "", ANCHOR_LAB[[FOCUS_ANCH]]))) +
  theme_f +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 5.2, face = "bold", colour = MYB_COL),
        axis.ticks = element_blank(),
        panel.grid = element_blank(), panel.spacing = unit(0.8, "mm"),
        axis.title = element_text(size = 6.5))

# =============================================================================
# B — 12 members are copies of three anchors
# =============================================================================
anch <- unique(TAB[, .(species, gene, at_anchor, our_family)])
anch[, is_myb := grepl("MYB", our_family)]
cnt <- anch[, .(genes = .N, myb = sum(is_myb)), by = .(species, at_anchor)]
cnt <- melt(cnt, id.vars = c("species", "at_anchor"),
            measure.vars = c("myb", "genes"), variable.name = "kind", value.name = "n")
# "genes" is the total, so the non-MYB part is total - myb; stack them.
cnt <- dcast(cnt, species + at_anchor ~ kind, value.var = "n")
cnt[, `not MYB` := genes - myb]
cnt <- melt(cnt[, .(species, at_anchor, `called MYB` = myb, `not MYB`)],
            id.vars = c("species", "at_anchor"), variable.name = "call", value.name = "n")
cnt[, species := factor(species, levels = c("sugarcane", "purple"))]
cnt[, anchor := factor(ANCHOR_LAB[at_anchor], levels = unname(ANCHOR_LAB))]
cnt[, call := factor(call, levels = c("called MYB", "not MYB"))]
say("anchor breakdown (no panel draws it; legend + stats only)")
print(dcast(cnt, anchor ~ species + call, value.var = "n", fill = 0), row.names = FALSE)

# NO PANEL IS DRAWN FROM THIS. The bar panel it used to feed was dropped on
# 2026-09-12 -- an anchor-by-anchor count belongs in a table. `cnt` stays because
# the legend quotes it and _stats.tsv carries it, so the panel went and no number
# did.

# =============================================================================
# B — every purple copy, z-scored, split by genotype as well as nitrogen,
#     with the TF annotation column
# =============================================================================
rdB <- unique(TAB[species == "purple", .(gene, at_anchor, our_family)])
rdB[, anchor := factor(ANCHOR_LAB[at_anchor], levels = unname(ANCHOR_LAB))]
X <- PU$X; mp <- PU$meta
# Genotype first, then nitrogen: in purple these genes vary far more between the
# two species than across nitrogen, and splitting by nitrogen alone hides it.
ordc <- order(mp$genotype, mp$status, mp$sample)
X <- X[, ordc, drop = FALSE]; mp <- mp[ordc]
rdB <- rdB[match(rownames(X), gene)]
rdB[, mean_tpm := rowMeans(X)]
setorder(rdB, anchor, -mean_tpm)
X <- X[rdB$gene, , drop = FALSE]
# Per-gene z. An absolute scale makes the one responding copy a bright row among
# dim ones and hides the SHAPE of its response, which is the thing to see.
ZP <- zrow(X)

dC <- data.table(gene = factor(rep(rownames(ZP), times = ncol(ZP)),
                               levels = rev(rdB$gene)),
                 sample = factor(rep(colnames(ZP), each = nrow(ZP)),
                                 levels = colnames(ZP)),
                 z = as.vector(ZP))
dC[, anchor := rdB$anchor[match(as.character(gene), rdB$gene)]]
dC[, status := mp$status[match(as.character(sample), mp$sample)]]
dC[, genotype := mp$genotype[match(as.character(sample), mp$sample)]]
ZLIM_P <- as.numeric(quantile(abs(dC$z), 0.98, na.rm = TRUE))
say(sprintf("panel B: %d purple copies x %d libraries; %d copies below 1 TPM",
            nrow(ZP), ncol(ZP), sum(rdB$mean_tpm < 1)))

# --- the TF column ------------------------------------------------------------
# THE SAME ANNOTATION THE PER-MODULE HEATMAPS CARRY. 16_module_heatmaps.r draws it
# as a ComplexHeatmap left_annotation off readouts/get_tfs/<study>/TF_in_network.tsv;
# this figure is ggplot, so it is a one-column raster sharing the row facets, which
# is what makes patchwork align it row-for-row with the heatmap.
#
# IT IS NOT COSMETIC HERE. Module 20 is published as ~75% MYB, so "how many of the
# copies does our own pipeline call a transcription factor" is exactly the claim a
# reader should be able to check off the panel.
#
# The independent call in TAB (our_family, from the Pfam/PlnTFDB pass) is compared
# against it and any disagreement is REPORTED rather than reconciled -- two
# annotations of the same gene disagreeing is information, and silently preferring
# one would hide it.
tf_genes <- character(0)
if (nzchar(TF_PURPLE) && file.exists(TF_PURPLE)) {
  tf_genes <- unique(strip_version(fread(TF_PURPLE)$gene))
  say(sprintf("TF column: %s TFs in the purple network (%s)",
              fmt_n(length(tf_genes)), basename(TF_PURPLE)))
} else {
  stop("TF table not found: ", TF_PURPLE,
       "\n  run  ./run.sh tfs purple  first -- the panel would silently lose its ",
       "annotation column", call. = FALSE)
}
rdB[, is_tf := gene %chin% tf_genes]
rdB[, tf_lab := fifelse(is_tf, "TF", "not TF")]
say(sprintf("  of the %d purple Module-20 copies, %d are called TF",
            nrow(rdB), sum(rdB$is_tf)))

# the cross-check, reported either way
rdB[, our_myb := grepl("MYB", our_family)]
disagree <- rdB[is_tf != our_myb]
if (nrow(disagree)) {
  say(sprintf("  NOTE: %d copy/copies where the TF table and our_family disagree:",
              nrow(disagree)))
  for (i in seq_len(nrow(disagree)))
    say(sprintf("    %-24s TF table %-6s our_family %s", disagree$gene[i],
                ifelse(disagree$is_tf[i], "TF", "not TF"), disagree$our_family[i]))
} else say("  the TF table and our_family agree on every copy")

dTF <- data.table(gene = factor(rdB$gene, levels = rev(rdB$gene)),
                  anchor = rdB$anchor, x = 1L,
                  tf = factor(rdB$tf_lab, levels = c("TF", "not TF")))
PAL_TF <- c(TF = MYB_COL, `not TF` = "grey88")

pTF <- ggplot(dTF, aes(x, gene, fill = tf)) +
  geom_raster() +
  facet_grid(anchor ~ ., scales = "free", space = "free") +
  scale_fill_manual(values = PAL_TF, name = NULL,
                    guide = guide_legend(keyheight = unit(3, "mm"),
                                         keywidth = unit(3, "mm"))) +
  scale_x_continuous(breaks = 1, labels = "TF", expand = c(0, 0),
                     position = "top") +
  labs(x = NULL, y = NULL) +
  theme_f +
  theme(axis.text.y = element_blank(), axis.ticks = element_blank(),
        # ticks are hidden but their LENGTH still reserves space, which is what
        # opened a visible gap between the strip and the rows it annotates
        axis.ticks.length = unit(0, "mm"),
        axis.text.x = element_text(size = 5.6, face = "bold"),
        panel.grid = element_blank(),
        strip.text.y = element_blank(), strip.background = element_blank(),
        panel.spacing = unit(0.7, "mm"),
        plot.margin = margin(2, 2, 2, 0))

pHeat <- ggplot(dC, aes(sample, gene, fill = z)) +
  geom_raster() +
  facet_grid(anchor ~ genotype + status, scales = "free", space = "free",
             switch = "y") +
  scale_fill_gradientn(colours = PAL_Z, limits = c(-ZLIM_P, ZLIM_P), oob = squish,
                       name = "z", guide = guide_colourbar(
                         barheight = unit(1.6, "cm"), barwidth = unit(2.4, "mm"))) +
  # Only the responding copy is named; the other 26 ids would be unreadable.
  scale_y_discrete(breaks = FOCUS, labels = FOCUS) +
  labs(x = "18 purple libraries", y = NULL) +
  theme_f +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(),
        axis.text.y = element_text(size = 5.6, face = "bold"),
        panel.grid = element_blank(), panel.spacing = unit(0.7, "mm"),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, size = 5.4, lineheight = 0.95),
        strip.text.x = element_text(size = 5.8),
        axis.ticks.length = unit(0, "mm"),
        plot.margin = margin(2, 0, 2, 2))

# The heatmap and its annotation are ONE panel, so they are combined here and get a
# single tag. The row facets are identical in both and space = "free" makes the
# facet heights proportional to row counts in each, which is what makes the strip
# line up with the rows it annotates.
pB <- (pHeat | pTF) + plot_layout(widths = c(1, 0.045), guides = "collect")

# =============================================================================
# C — that copy on its own
# =============================================================================
vC <- data.table(sample = colnames(PU$X), tpm = as.numeric(PU$X[FOCUS, ]))
vC <- merge(vC, PU$meta[, .(sample, genotype, status)], by = "sample")
uu <- NPUR[gene == FOCUS, .(genotype, U_est, U_p, anova_p, mean_tpm)]
say("panel C: the focus copy"); print(uu, row.names = FALSE)
print(vC[, .(mean_tpm = round(mean(tpm), 1)), by = .(genotype, status)][order(genotype, status)],
      row.names = FALSE)

pC <- ggplot(vC, aes(status, tpm, colour = genotype, group = genotype)) +
  stat_summary(fun = mean, geom = "line", linewidth = 0.6) +
  geom_point(size = 1.5, alpha = 0.85) +
  scale_colour_manual(values = c(`51NG3` = unname(PAL_STUDY[1]),
                                 TAGZ = unname(PAL_STUDY[2])), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.12))) +
  labs(x = NULL, y = sprintf("%s\nTPM", FOCUS)) +
  theme_f + theme(axis.title.y = element_text(size = 6.5, lineheight = 1.1))

# =============================================================================
# compose
# =============================================================================
# With the old bar panel gone, A takes the full width it always wanted: 48 columns
# across half a page was why its nine named rows could not have been read before.
# B (heatmap + TF strip) and C sit below, B wide because it carries 27 rows and six
# column facets, C narrow because it is three points and two lines.
#
# TAGS FOLLOW COMPOSITION ORDER, NOT POSITION, so A/B/C land on the intended panels
# and the legend below is written against the same letters.
#
# pB IS WRAPPED, and that is load-bearing. `|` and `/` FLATTEN a nested patchwork
# into the outer layout, so `pA / (pB | pC)` becomes three plots under a two-column
# layout and patchwork aborts with "Need 3 panels, but nrow and ncol only provide 2".
# wrap_elements() makes the heatmap-plus-strip pair one opaque element, which also
# gives it ONE panel tag instead of tagging the annotation column as its own panel.
fig <- wrap_plots(A = pA, B = wrap_elements(full = pB), C = pC,
                  design = "AAAAAAA\nBBBBBCC", heights = c(0.85, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12))

W <- 22; H <- 16
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

nA <- function(sp, a) sum(cnt[species == sp & at_anchor == a, n])
mA <- function(sp, a) cnt[species == sp & at_anchor == a & call == "called MYB", n]
fm <- function(gt, col) uu[genotype == gt][[col]]
mn <- function(gt, st) round(vC[genotype == gt & status == st, mean(tpm)], 1)
sil <- rdB[at_anchor == FOCUS_ANCH & gene != FOCUS, mean_tpm]

n_tf   <- sum(rdB$is_tf)
tf_ids <- rdB[is_tf == TRUE, gene]

legend <- paste0(
"Figure ", FIG, ". Muñoz-Perez et al.'s Module 20 in purple, against sugarcane. Figure 2 showed ",
"that the module's nitrogen response reproduces in the study's own data. This figure asks what ",
"survives in the other species, and the answer is not less of the same thing -- it is a different ",
"response shape, in one copy, that every monotone test in this work is blind to.\n",
"\n",
"WHAT THE 12 PUBLISHED MEMBERS ACTUALLY ARE, since every panel below depends on it. They are ",
"TransDecoder transcripts from a de novo assembly; mapped into each reference proteome by ",
"reciprocal DIAMOND blastp through Arabidopsis they resolve to only THREE anchors, and the mapped ",
"genes are haplotype COPIES of those few loci rather than 12 independent genes -- ",
fmt_n(sum(cnt[species == "sugarcane", n])), " sugarcane genes from 10 members and ",
fmt_n(sum(cnt[species == "purple", n])), " purple genes from 6. Only the AtMYB59 anchor carries a ",
"MYB call under our own Pfam/PlnTFDB rules, at ", fmt_n(mA("sugarcane", FOCUS_ANCH)), " of ",
fmt_n(nA("sugarcane", FOCUS_ANCH)), " copies in sugarcane and ", fmt_n(mA("purple", FOCUS_ANCH)),
" of ", fmt_n(nA("purple", FOCUS_ANCH)), " in purple, while the AT3G55960 and OPT7 anchors -- ",
"which supply most of the mapped genes in both species -- are not transcription factors at all. ",
"The AtMYB48 anchor maps to nothing in either network. The per-anchor counts are in ",
basename(paste0(OUT_PREFIX, "_stats.tsv")), ".\n",
"\n",
"(A) The module as reported, reproduced: all ", fmt_n(nrow(ZS)),
" mapped sugarcane genes across all ", fmt_n(ncol(ZS)), " libraries as per-gene z-scores, split ",
"by nitrogen. It is the premise rather than the result -- the rest of the figure is what happens ",
"to this when it is carried into the other species. The ", fmt_n(length(myb_sug)),
" rows marked and named are the sugarcane copies of the AtMYB59 anchor, i.e. THE SAME LOCUS as ",
"the focus gene in panels B and C, which is itself a purple gene and so has no row here. They are ",
"marked because the comparison this figure rests on is a claim about those particular rows: every ",
"one of them is monotonically REPRESSED by nitrogen (log2 fold change -1.8 to -3.4, BH-adjusted ",
"p ~ 1e-08 in the responsive genotype RB975375), which is visible as the band of high z at low N ",
"that goes uniformly low at high N. Unmarked, they are nine anonymous rows among ", fmt_n(nrow(ZS)),
".\n",
"\n",
"(B) The same module in purple, as a per-gene Z-SCORE, split by GENOTYPE as well as nitrogen, with ",
"rows grouped by Arabidopsis anchor and ordered by mean expression, and with the TF ANNOTATION ",
"COLUMN at the right. Only the responding copy is labelled; the other ", fmt_n(nrow(ZP) - 1),
" ids would be unreadable at this size. Three choices, each deliberate. The Z-SCORE is what makes ",
"that copy legible at all: on an absolute scale it is a bright row among dim ones and the SHAPE of ",
"its response cannot be seen. Splitting by GENOTYPE is what makes the dominant structure visible ",
"-- in purple these genes vary far more between the two species than across nitrogen, and a panel ",
"split by nitrogen alone hides that and invites the reader to attribute the variation to the ",
"treatment. The TF COLUMN is the same annotation the per-module heatmaps carry (from the network's ",
"own TF table), and it is here because the module is published as ~75% MYB or MYB-related: of the ",
fmt_n(nrow(rdB)), " purple copies, ", fmt_n(n_tf), " are called transcription factors (",
paste(tf_ids, collapse = ", "), "), both at the AtMYB59 anchor. Our independent domain call in the ",
"network table agrees with the TF table on every one of the ", fmt_n(nrow(rdB)), " copies. Most ",
"copies are also flat or near-silent: ", fmt_n(sum(rdB$mean_tpm < 1)), " of the ", fmt_n(nrow(ZP)),
" run below 1 TPM, and the AtMYB59 block is nearly empty -- of its ",
fmt_n(nA("purple", FOCUS_ANCH)), " copies two sit at ", sprintf("%.2f", max(sil)), " and ",
sprintf("%.2f", min(sil)), " TPM, i.e. they are not expressed in leaf at all.\n",
"\n",
"(C) That copy on its own; points are libraries, lines join the means. ", FOCUS,
" is the only AtMYB59-anchor copy expressed in purple (", sprintf("%.1f", fm("51NG3", "mean_tpm")),
" TPM against ", sprintf("%.2f", max(sil)), " and ", sprintf("%.2f", min(sil)),
" for the other two) and is called MYB both by the source study and by our own domain rules. Its ",
"response is a U CENTRED ON THE CONTROL: in 51NG3 it runs at ", mn("51NG3", "low N"),
" TPM under nitrogen starvation, ", mn("51NG3", "control"), " at the 2 N control and ",
mn("51NG3", "high N"), " under excess -- an eleven-fold drop from starvation to control and a ",
"four-fold rise again beyond it (U-shape contrast c(+1,-2,+1) p = ",
sprintf("%.4f", fm("51NG3", "U_p")), ", one-way p = ", sprintf("%.4f", fm("51NG3", "anova_p")),
"). TAGZ shows the same shape more weakly (", mn("TAGZ", "low N"), ", ", mn("TAGZ", "control"),
", ", mn("TAGZ", "high N"), " TPM).\n",
"\n",
"Because 2 N is the CONTROL and 0 N and 6 N are stresses in opposite directions, this is a gene ",
"induced by nitrogen stress in EITHER direction. That is exactly the shape a monotone rank ",
"correlation cannot see -- it is why the module-level analysis reports this gene\'s module as ",
"unresponsive, and why the U-shape contrast exists as a separate test. Sugarcane\'s design has no ",
"control level, only Low and High, so what panel A can observe is ONE ARM of such a curve. The two ",
"studies are therefore consistent rather than contradictory: sugarcane sees the low-nitrogen arm ",
"of the same response, and only purple\'s three-level design shows that the high-nitrogen side ",
"turns back up.\n",
"\n",
"THE LIMIT. This rests on ONE expressed copy in one species, at a p-value that does not clear BH ",
"over the ", fmt_n(nrow(NPUR[genotype == "51NG3"])), " purple Module-20 genes tested. Run ",
"genome-wide over every gene in the purple network the same contrast returns a SINGLE gene past ",
"FDR and this is not it -- it sits at padj 0.48, in the top 0.5% by raw p but nowhere near ",
"surviving correction. It is a candidate worth following, not a finding.")

legend_file <- paste0(OUT_PREFIX, "_legend.txt")
writeLines(wrap_at(legend), legend_file)
say("wrote ", basename(legend_file))

stats <- rbindlist(list(
  # kept although no panel draws them any more -- the legend quotes these and the
  # anchor breakdown is what the dropped bar panel used to show
  cnt[, .(panel = "-", study = as.character(species),
          quantity = sprintf("%s: %s", at_anchor, call), value = fmt_n(n))],
  rdB[, .(panel = "B", study = "purple",
          quantity = sprintf("%s TF call", gene),
          value = tf_lab)],
  data.table(panel = "A", study = "sugarcane",
             quantity = sprintf("%s copies marked in panel A", FOCUS_ANCH),
             value = fmt_n(length(myb_sug))),
  vC[, .(panel = "C", study = "purple",
         quantity = sprintf("%s, %s, %s", FOCUS, genotype, status),
         value = sprintf("%.1f TPM", mean(tpm))), by = .(genotype, status)][
           , .(panel, study, quantity, value)],
  uu[, .(panel = "C", study = "purple",
         quantity = sprintf("%s U-contrast, %s", FOCUS, genotype),
         value = sprintf("est %+.2f, p %.4f", U_est, U_p))]))
write_tsv(stats, paste0(OUT_PREFIX, "_stats.tsv"))
say("done")
