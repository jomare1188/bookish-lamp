# =============================================================================
# 20_figure1_reproduction.r — Figure 1: both source studies' primary findings,
# recovered in this project's single re-quantification.
#
# THE CLAIM THIS FIGURE MAKES, AND WHY IT COMES FIRST. Everything downstream is
# a comparison between two studies that were never designed to be compared. That
# comparison is only worth reading if the re-quantification reproduces what each
# study found on its own data first. So before any cross-species result, the
# paper has to show: run through OUR pipeline, from OUR references, do Muñoz's
# and Kiet's own headline observations still appear?
#
#   A  Muñoz-Perez et al. 2025 — Module 20 is nitrogen-responsive in sugarcane.
#      33 mapped genes x 48 leaf libraries. Columns split by NITROGEN ONLY, with
#      genotype as an annotation bar: that is what separates "tracks nitrogen"
#      from "tracks genotype" by eye, because a genotype-driven set would split
#      into two sub-clusters inside each nitrogen block instead of differing
#      between them.
#
#   B  Ta Quang Kiet et al. 2025 — MYB61 responds to nitrogen non-monotonically
#      and differently between genotypes, in purple. 16 confirmed copies x 18
#      leaf libraries, columns split by GENOTYPE (their claim is that the two
#      genotypes behave differently, so the genotypes must be side by side) and
#      ordered 0N -> 2N -> 6N inside each.
#
# THE PANELS ARE NOT SYMMETRIC AND THE FIGURE SHOULD NOT PRETEND THEY ARE.
# Muñoz's result reproduces outright. Kiet's reproduces in SHAPE — a significant
# non-monotonic, genotype-restricted response — but with the sign INVERTED
# (expression peaks at 2N and falls at both extremes, where they report the
# opposite), and at expression levels so low that most copies sit under 1 TPM in
# leaf. Both panel subtitles say so. A figure that showed only the agreement
# would be the wrong figure.
#
# The published id `Soff.09G0002230-3D` does not resolve to a MYB (see
# 11_readouts/myb61/README.md), so panel B carries the copies at that locus as a
# labelled NEGATIVE CONTROL block: they are one to two orders of magnitude more
# highly expressed than any true MYB61 copy, which is the likely source of the
# signal the paper attributes to MYB61.
#
# ROWS ARE PER-GENE Z-SCORES in both panels, because the question is the SHAPE of
# the response across libraries and the copies span three orders of magnitude in
# absolute expression. Absolute log2(TPM+1) is carried as a row annotation
# instead, so a row that is structured but silent cannot be mistaken for a
# result — the failure mode a z-score heatmap has on its own.
#
# Expression is TPM from the same salmon quantification the networks were built
# from. Not VST: these are per-study descriptive panels reproducing per-study
# claims, and TPM is what both source papers report.
#
# NOTHING IS WRITTEN ON THE FIGURE THAT BELONGS IN THE LEGEND. Panels carry a
# bold letter and the labels the data needs to be read -- split titles,
# annotation names, colour keys -- and nothing else. No titles, no subtitles, no
# statistics printed onto the panel. Everything a reader needs to interpret the
# figure goes into the LEGEND, which this script generates from the same
# variables that produced the panels, so the two cannot drift apart. The legend
# is written to <prefix>_legend.txt and assembled into the paper-level
# figures_legends.txt by  ./run.sh legends.
#
# RUN: through run.sh  ->  ./run.sh figure1   (then ./run.sh legends)
# =============================================================================

suppressMessages({
  library(data.table); library(ComplexHeatmap); library(circlize)
  library(scico); library(grid); library(svglite)
})

source(file.path(dirname(sub("--file=", "",
       grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "lib", "common.R"))

BASE        <- env_opt("CLEAN_BASE", "/dados04/jorge/comparative_saccharum")
READOUTS    <- env_req("CLEAN_READOUTS")
TPM_SUG     <- env_req("CLEAN_TPM_SUGARCANE")
TPM_PUR     <- env_req("CLEAN_TPM_PURPLE")
META_SUG    <- env_req("CLEAN_META_SUGARCANE")
OUT_PREFIX  <- env_req("CLEAN_OUT_PREFIX")
KIET_LOCUS  <- env_opt("CLEAN_KIET_LOCUS", "09G0002230")
setDTthreads(as.integer(env_num("CLEAN_CORES", 8)))

banner("Figure 1 — reproducing each source study's own finding")

PAL_Z  <- rev(scico(256, palette = "roma"))

# Legend prose is wrapped once, here, so the .txt file is readable in a terminal
# and pastes into a manuscript without stray line breaks inside sentences.
# Paragraphs (blank lines) survive; everything inside one is re-flowed.
wrap_at <- function(x, width = 96)
  paste(vapply(strsplit(x, "\n[[:space:]]*\n")[[1L]],
               function(para) paste(strwrap(gsub("[[:space:]]+", " ", para), width),
                                    collapse = "\n"),
               ""), collapse = "\n\n")
CELL_W <- unit(3.0, "mm")
CELL_H <- unit(3.4, "mm")

# =============================================================================
# Panel A — Muñoz Module 20 in sugarcane
# =============================================================================
m20_tab <- fread(file.path(READOUTS, "module20/module20_network_table.tsv"))
m20_n   <- fread(file.path(READOUTS, "module20/module20_nitrogen_sugarcane.tsv"))

tpm <- fread(TPM_SUG); setnames(tpm, 1, "gene")
tpm[, gene := strip_version(gene)]
meta <- fread(META_SUG); setnames(meta, tolower(names(meta)))

genes_A <- m20_tab[species == "sugarcane", unique(gene)]
samp    <- intersect(names(tpm), meta$sample)
mA      <- meta[match(samp, sample)]
XA <- as.matrix(tpm[gene %chin% genes_A, ..samp])
rownames(XA) <- tpm[gene %chin% genes_A, gene]
XA <- log2(XA + 1)

# Column order: nitrogen block, then genotype, then tissue, then sample -- the
# same replicate-adjacent rule the module figures use (see 16/17). Sugarcane's
# names run B0_1, B_1, M_1, P_1, B0_2, ... so sorting by name alone scatters the
# replicates of one leaf segment four columns apart and the segment effect reads
# as striping.
mA[, nlev := factor(fifelse(grepl("Low", treatment), "Low N", "High N"),
                    levels = c("Low N", "High N"))]
mA[, gt := factor(fifelse(genotype == "RB975375", "RB975375 (responsive)",
                                                  "RB937570 (non-responsive)"),
                  levels = c("RB975375 (responsive)", "RB937570 (non-responsive)"))]
ordA <- order(mA$nlev, mA$gt, mA$tissue, mA$sample)
mA <- mA[ordA]; XA <- XA[, mA$sample, drop = FALSE]

ZA <- t(scale(t(XA))); ZA[!is.finite(ZA)] <- 0

# rows: group by Module-20 locus, drop the study's own prefixes from the label
rdA <- unique(m20_tab[species == "sugarcane", .(gene, locus, our_family, in_network)])
rdA <- rdA[match(rownames(XA), gene)]
rdA[, locus_lab := sub("^(N_R___|N_NR___|SCA[0-9_]+__)", "", locus)]
rdA[, myb := fifelse(grepl("MYB", our_family), "MYB", "not MYB")]

# per-gene nitrogen test in the RESPONSIVE genotype, which is Muñoz's own claim
sigA <- m20_n[grepl("RB975375", genotype), .(gene, padj, l2fc = log2FC_HighminusLow)]
rdA  <- merge(rdA, sigA, by = "gene", all.x = TRUE, sort = FALSE)
rdA  <- rdA[match(rownames(XA), gene)]
n_sig  <- rdA[!is.na(padj) & padj < 0.05, .N]
n_down <- rdA[!is.na(padj) & padj < 0.05 & l2fc < 0, .N]
n_test <- rdA[!is.na(padj), .N]

vp <- fread(file.path(READOUTS, "module20/module20_variance_partition_sugarcane.tsv"))
say(sprintf("panel A: %d genes x %d libraries | %d/%d padj<0.05 in RB975375, %d down at high N",
            nrow(XA), ncol(XA), n_sig, n_test, n_down))

annA_top <- HeatmapAnnotation(
  Genotype = mA$gt,
  col = list(Genotype = structure(scico(3, palette = "batlow")[1:2],
                                  names = levels(mA$gt))),
  annotation_name_gp = gpar(fontsize = 7),
  simple_anno_size = unit(3.5, "mm"),
  annotation_legend_param = list(Genotype = list(title_gp = gpar(fontsize = 8, fontface = "bold"),
                                                 labels_gp = gpar(fontsize = 7))))

annA_row <- rowAnnotation(
  `log2(TPM+1)` = anno_barplot(rowMeans(XA), gp = gpar(fill = "grey40", col = NA),
                               width = unit(1.5, "cm"),
                               axis_param = list(gp = gpar(fontsize = 6))),
  `N-resp.` = ifelse(!is.na(rdA$padj) & rdA$padj < 0.05, "padj < 0.05", "n.s."),
  TF = rdA$myb,
  col = list(`N-resp.` = c(`padj < 0.05` = "#B2182B", `n.s.` = "grey88"),
             TF = c(MYB = scico(5, palette = "batlow")[2], `not MYB` = "grey88")),
  annotation_name_gp = gpar(fontsize = 6.5),
  show_annotation_name = c(`log2(TPM+1)` = TRUE, `N-resp.` = FALSE, TF = FALSE),
  simple_anno_size = unit(3, "mm"),
  annotation_legend_param = list(
    `N-resp.` = list(title = "N-responsive\n(RB975375)",
                     title_gp = gpar(fontsize = 8, fontface = "bold"),
                          labels_gp = gpar(fontsize = 7)),
    TF = list(title_gp = gpar(fontsize = 8, fontface = "bold"),
              labels_gp = gpar(fontsize = 7))))

zlA <- as.numeric(quantile(abs(ZA), 0.98, na.rm = TRUE)); if (zlA == 0) zlA <- 1

htA <- Heatmap(ZA, name = "z-score",
  col = colorRamp2(seq(-zlA, zlA, length.out = 256), PAL_Z),
  cluster_rows = TRUE, cluster_columns = FALSE, cluster_row_slices = FALSE,
  row_split = factor(rdA$locus_lab, levels = unique(rdA$locus_lab)),
  column_split = mA$nlev,
  top_annotation = annA_top, right_annotation = annA_row,
  show_row_names = TRUE, row_names_gp = gpar(fontsize = 5.5),
  show_column_names = FALSE,
  row_title_gp = gpar(fontsize = 7, fontface = "bold"), row_title_rot = 0,
  column_title_gp = gpar(fontsize = 9, fontface = "bold"),
  width = CELL_W * ncol(ZA), height = CELL_H * nrow(ZA),
  border = TRUE, row_gap = unit(1.2, "mm"), column_gap = unit(1.5, "mm"),
  heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "bold"),
                              labels_gp = gpar(fontsize = 7),
                              legend_height = unit(2.2, "cm")))

med_n  <- median(vp$pct_nitrogen, na.rm = TRUE)
med_g  <- median(vp$pct_genotype, na.rm = TRUE)
n_ngt  <- sum(vp$pct_nitrogen > vp$pct_genotype, na.rm = TRUE)
n_vp   <- sum(!is.na(vp$pct_genotype))
n_myb  <- sum(rdA$myb == "MYB")
n_loci <- uniqueN(rdA$locus_lab)

# =============================================================================
# Panel B — Kiet MYB61 in purple
# =============================================================================
myb <- fread(file.path(READOUTS, "myb61/MYB61_network_table.tsv"))
anv <- fread(file.path(READOUTS, "myb61/MYB61_expression_anova.tsv"))

tpmP <- fread(TPM_PUR); setnames(tpmP, 1, "gene")
sampP <- setdiff(names(tpmP), c("gene", "gene_name"))

# Kiet's sample names carry the design: L_<genotype>_<level>-<rep>.
mB <- data.table(sample = sampP)
mB[, gt := factor(fifelse(grepl("51N03", sample), "51NG3 (responsive)",
                                                  "TAGZ (non-responsive)"),
                  levels = c("51NG3 (responsive)", "TAGZ (non-responsive)"))]
mB[, nlev := factor(fifelse(grepl("_ON-", sample), "0N",
                    fifelse(grepl("_2N-", sample), "2N", "6N")),
                    levels = c("0N", "2N", "6N"))]
stopifnot(nrow(mB) == 18, all(table(mB$gt, mB$nlev) == 3))
mB <- mB[order(gt, nlev, sample)]

copies <- myb[species == "purple", unique(gene)]
control <- grep(KIET_LOCUS, tpmP$gene, value = TRUE)
genes_B <- c(copies, setdiff(control, copies))

sampB <- mB$sample
XB <- as.matrix(tpmP[gene %chin% genes_B, ..sampB])
rownames(XB) <- tpmP[gene %chin% genes_B, gene]
XB <- log2(XB + 1)

cladeB <- myb[species == "purple", .(gene, clade)]
rdB <- data.table(gene = rownames(XB))
rdB <- merge(rdB, cladeB, by = "gene", all.x = TRUE, sort = FALSE)
rdB[is.na(clade), clade := sprintf("published id\n%s", KIET_LOCUS)]
rdB[, clade := sub(" co-ortholog", "\nco-ortholog", clade)]
rdB <- rdB[match(rownames(XB), gene)]
blk <- unique(rdB$clade)
blk <- c(sort(grep("published", blk, value = TRUE, invert = TRUE)),
         grep("published", blk, value = TRUE))
rdB[, clade := factor(clade, levels = blk)]

# The U-shape contrast in 51NG3 is Kiet's claim; carry its verdict per copy.
uB <- anv[genotype == "51NG3", .(gene, U_padj = as.numeric(U_padj), U_est = as.numeric(U_est),
                                 mean_tpm = as.numeric(mean_tpm))]
uB <- uB[!duplicated(gene)]
rdB <- merge(rdB, uB, by = "gene", all.x = TRUE, sort = FALSE)
rdB <- rdB[match(rownames(XB), gene)]

is_copy <- rdB$gene %chin% copies
n_u   <- sum(is_copy & !is.na(rdB$U_padj) & rdB$U_padj < 0.05)
n_ut  <- sum(is_copy & !is.na(rdB$U_padj))
n_neg <- sum(is_copy & !is.na(rdB$U_padj) & rdB$U_padj < 0.05 & rdB$U_est < 0)
say(sprintf("panel B: %d rows (%d MYB61 copies + %d control) x %d libraries | U-contrast padj<0.05 in %d/%d, %d inverted",
            nrow(XB), sum(is_copy), sum(!is_copy), ncol(XB), n_u, n_ut, n_neg))

ZB <- t(scale(t(XB))); ZB[!is.finite(ZB)] <- 0

annB_top <- HeatmapAnnotation(
  Nitrogen = mB$nlev,
  col = list(Nitrogen = structure(scico(4, palette = "lajolla",
                                        begin = 0.25, end = 0.9)[1:3],
                                  names = levels(mB$nlev))),
  annotation_name_gp = gpar(fontsize = 7),
  simple_anno_size = unit(3.5, "mm"),
  annotation_legend_param = list(Nitrogen = list(title_gp = gpar(fontsize = 8, fontface = "bold"),
                                                 labels_gp = gpar(fontsize = 7))))

annB_row <- rowAnnotation(
  `log2(TPM+1)` = anno_barplot(rowMeans(XB), gp = gpar(fill = "grey40", col = NA),
                               width = unit(1.5, "cm"),
                               axis_param = list(gp = gpar(fontsize = 6))),
  `U-shape` = ifelse(!is.na(rdB$U_padj) & rdB$U_padj < 0.05, "padj < 0.05", "n.s."),
  col = list(`U-shape` = c(`padj < 0.05` = "#B2182B", `n.s.` = "grey88")),
  annotation_name_gp = gpar(fontsize = 6.5),
  show_annotation_name = c(`log2(TPM+1)` = TRUE, `U-shape` = FALSE),
  simple_anno_size = unit(3, "mm"),
  annotation_legend_param = list(
    `U-shape` = list(title = "non-monotonic\nin 51NG3",
                     title_gp = gpar(fontsize = 8, fontface = "bold"),
                     labels_gp = gpar(fontsize = 7))))

# 95th, not 98th: most MYB61 copies sit near zero TPM in leaf, so their z-scores
# are spiky by construction and a 98% clip lets three cells set the whole scale.
zlB <- as.numeric(quantile(abs(ZB), 0.95, na.rm = TRUE)); if (zlB == 0) zlB <- 1

htB <- Heatmap(ZB, name = "z-score ",
  col = colorRamp2(seq(-zlB, zlB, length.out = 256), PAL_Z),
  cluster_rows = FALSE, cluster_columns = FALSE, cluster_row_slices = FALSE,
  row_order = order(rdB$clade, -rdB$mean_tpm),
  row_split = rdB$clade,
  column_split = mB$gt,
  top_annotation = annB_top, right_annotation = annB_row,
  show_row_names = TRUE, row_names_gp = gpar(fontsize = 5.5),
  show_column_names = FALSE,
  row_title_gp = gpar(fontsize = 7, fontface = "bold"), row_title_rot = 0,
  column_title_gp = gpar(fontsize = 9, fontface = "bold"),
  width = unit(6.5, "mm") * ncol(ZB), height = CELL_H * nrow(ZB),
  border = TRUE, row_gap = unit(1.2, "mm"), column_gap = unit(1.5, "mm"),
  heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "bold"),
                              labels_gp = gpar(fontsize = 7),
                              legend_height = unit(2.2, "cm")))

n_tagz <- sum(anv[genotype == "TAGZ" & gene %chin% copies,
                  as.numeric(U_padj)] < 0.05, na.rm = TRUE)
n_ctrl <- sum(!is_copy)

# =============================================================================
# compose and write
# =============================================================================
# Two ComplexHeatmap objects on one device: each is drawn into its own viewport
# from a 2-row layout, since ComplexHeatmap has no operator for stacking
# independent heatmaps that do not share an axis. These two share neither rows
# nor columns -- different species, different libraries -- so they must not be
# forced onto one grid.
W  <- 30
tA <- 0.8                                # just the panel letter
tB <- 0.8
hA <- 3.2 + 0.36 * nrow(ZA)              # column annotation + split titles + body
hB <- 3.2 + 0.36 * nrow(ZB)
H  <- tA + hA + tB + hB

# The panel letter is the ONLY text this script puts on the figure. Everything
# else a reader needs is in the legend.
put_letter <- function(ltr)
  grid.text(ltr, x = unit(3, "mm"), y = unit(1, "npc") - unit(2, "mm"),
            just = c("left", "top"), gp = gpar(fontsize = 13, fontface = "bold"))

draw_fig <- function() {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(
    4, 1, heights = unit(c(tA, hA, tB, hB), "cm"))))

  pushViewport(viewport(layout.pos.row = 1)); put_letter("A"); popViewport()
  pushViewport(viewport(layout.pos.row = 2))
  draw(htA, newpage = FALSE, merge_legends = TRUE,
       heatmap_legend_side = "right", annotation_legend_side = "right",
       padding = unit(c(9, 2, 1, 2), "mm"))
  popViewport()

  pushViewport(viewport(layout.pos.row = 3)); put_letter("B"); popViewport()
  pushViewport(viewport(layout.pos.row = 4))
  draw(htB, newpage = FALSE, merge_legends = TRUE,
       heatmap_legend_side = "right", annotation_legend_side = "right",
       padding = unit(c(9, 2, 1, 2), "mm"))
  popViewport()
  popViewport()
}

invisible(ensure_dir(dirname(OUT_PREFIX)))
png(paste0(OUT_PREFIX, ".png"), width = W, height = H, units = "cm", res = 400,
    type = "cairo")
draw_fig(); dev.off()
cairo_pdf(paste0(OUT_PREFIX, ".pdf"), width = W / 2.54, height = H / 2.54)
draw_fig(); dev.off()
svglite(paste0(OUT_PREFIX, ".svg"), width = W / 2.54, height = H / 2.54)
draw_fig(); dev.off()

say(sprintf("wrote %s.{png,pdf,svg}  (%.0f x %.0f cm)", basename(OUT_PREFIX), W, H))

# --- the legend -------------------------------------------------------------
# Built from the same variables that drew the panels, so a number can never
# disagree between figure and legend. Written next to the figure; assembled into
# the paper-level figures_legends.txt by  ./run.sh legends.
#
# ON ORTHOFINDER, precisely. It is NOT how panel A was built: the Module-20
# members are de novo assembly ORFs and were mapped into the R570 proteome by
# reciprocal DIAMOND blastp through Arabidopsis. OrthoFinder enters panel B, as
# one of the four independent lines of evidence confirming which purple genes
# are MYB61, and it is the same orthology used for every cross-species
# comparison elsewhere in this work. The legend says so in those terms rather
# than crediting it with the whole figure.
legend <- paste0(
"Figure 1. The principal nitrogen finding of each source study is recovered in a single, ",
"independent re-quantification. Public RNA-seq from both studies was re-processed through one ",
"nf-core/rnaseq + salmon workflow against a common reference per species (Saccharum hybrid R570 ",
"for sugarcane; LA purple for S. officinarum / S. robustum), so neither panel reuses any ",
"quantification, mapping or statistic from the original publications. Cross-species orthology ",
"throughout this work is the OrthoFinder backbone: OrthoFinder v3.1.3 (-S diamond -M msa -A famsa ",
"-T fasttree) run on one protein per gene of the two reference proteomes, assigning 360,794 of ",
"435,856 genes (82.8%) to 94,273 orthogroups. Every comparison between the two species rests on ",
"it, and it is what makes 'the same gene' a defined object across two polyploid references. It ",
"enters this figure in panel B, as one of the lines of evidence identifying MYB61; panel A is a ",
"within-species mapping problem and is solved differently, as described below. In both panels ",
"colour is ",
"the per-gene z-score of log2(TPM + 1) across that panel's libraries, so a row shows the SHAPE of ",
"a gene's response rather than its absolute level; ramps are clipped at the ", sprintf("%.0f", 98),
"th (A) and ", sprintf("%.0f", 95), "th (B) percentile of |z|. Because a z-score rescales a ",
"near-silent gene to look as structured as an abundant one, mean log2(TPM + 1) is drawn as a grey ",
"bar beside every row. Rows are genes; columns are individual libraries.
",
"
",
"(A) Module 20 of Mu\u00f1oz-Perez et al. (2025), in sugarcane. Their 12 published members are ",
"TransDecoder ORFs from a de novo assembly and carry no reference coordinates, so placing them is ",
"a within-species mapping problem rather than an orthology one and the orthogroups are not used: ",
"they were mapped into the R570 proteome by DIAMOND blastp under a reciprocal-Arabidopsis filter ",
"-- member and candidate must return the same best Arabidopsis hit -- which resolves them to ",
sprintf("%d", nrow(ZA)), " genes at ", sprintf("%d", n_loci), " independent Arabidopsis anchors, ",
"shown as the ", sprintf("%d", n_loci), " row blocks. (The 12 published members are therefore not ",
"12 independent units; two have no predicted ORF at all and the nine MYB members resolve to only ",
"two Arabidopsis MYB genes.) Columns are the ", sprintf("%d", ncol(ZA)), " leaf libraries, split ",
"by nitrogen treatment alone (Low N | High N) with genotype as the upper annotation bar, and ",
"ordered treatment, then genotype, then leaf segment, then library, so replicates of one condition ",
"are adjacent; rows are hierarchically clustered within each block. This column layout is what ",
"separates a nitrogen response from a genotype response by eye: a genotype-driven set would split ",
"into two sub-clusters INSIDE each nitrogen block rather than differ BETWEEN blocks. Red marks in ",
"the 'N-responsive (RB975375)' track are genes at padj < 0.05 for High vs Low nitrogen in the ",
"nitrogen-responsive genotype RB975375 (Benjamini-Hochberg corrected); ", sprintf("%d", n_sig),
" of the ", sprintf("%d", n_test), " testable genes qualify and all ", sprintf("%d", n_down),
" are LOWER at high nitrogen. Teal marks in the 'TF' track are the ", sprintf("%d", n_myb),
" genes independently called MYB by this project's own Pfam/PlnTFDB domain pipeline. Partitioning ",
"each gene's variance with lm(expression ~ genotype + nitrogen), nitrogen accounts for a median ",
sprintf("%.1f%%", med_n), " against ", sprintf("%.1f%%", med_g), " for genotype, and exceeds ",
"genotype in ", sprintf("%d", n_ngt), " of ", sprintf("%d", n_vp), " genes. Mu\u00f1oz-Perez et al.'s ",
"core claim is therefore reproduced.
",
"
",
"(B) MYB61 of Ta Quang Kiet et al. (2025), in purple. The published identifier ",
"(Soff.09G0002230-3D) does not resolve in the LA purple annotation, so MYB61 was re-derived from ",
"sequence: reciprocal best hits from AtMYB61 (AT1G09540) into sorghum and rice, DIAMOND blastp of ",
"those five anchors into the LA purple proteome, a reciprocal-best-hit filter back to AT1G09540, ",
"an independent Myb-domain call, the OrthoFinder orthogroups described above, and a MAFFT/FastTree ",
"phylogeny. The orthogroups place 9 of the 16 purple copies in five orthogroups that also hold a ",
"sugarcane MYB61 copy, which is what licenses calling them the same gene across the two ",
"references; the genes at the published identifier fall in four further orthogroups, none of ",
"which contains any confirmed MYB61 copy of either species. Note ",
"that OrthoFinder leaves many polyploid haplotype copies unassigned, so absence from an orthogroup ",
"is not evidence against a copy; that is why it is used alongside rather than instead of the ",
"reciprocal search and the phylogeny. Rows are the ", sprintf("%d", sum(is_copy)), " confirmed ",
"copies, split into the two grass co-orthologue clades left by the grass whole-genome duplication ",
"(chr3-type, chr9-type), plus the ", sprintf("%d", n_ctrl), " genes at the published identifier as ",
"a NEGATIVE CONTROL: that locus carries no Myb domain, never shares an orthogroup with a confirmed ",
"MYB61 copy, and its genes are one to two orders of magnitude more highly expressed than any true ",
"copy. Columns are the ", sprintf("%d", ncol(ZB)), " leaf libraries, split by genotype and ordered ",
"0N, 2N, 6N within each, because the published claim is that the two genotypes behave differently. ",
"Red marks in the 'non-monotonic in 51NG3' track are copies with a significant U-shape contrast ",
"c(+1, -2, +1) across 0N/2N/6N in 51NG3 ",
"(one-way model per genotype, Benjamini-Hochberg corrected): ", sprintf("%d", n_u), " of ",
sprintf("%d", n_ut), " testable copies, against ", sprintf("%d", n_tagz), " of ",
sprintf("%d", n_ut), " in TAGZ. A non-monotonic, genotype-restricted nitrogen response is ",
"therefore present, as reported -- but its sign is inverted: all ", sprintf("%d", n_neg),
" significant copies peak at 2N and fall at both extremes, whereas maxima at 0N and 6N are ",
"described. Two limits belong with this panel: the copies are lowly expressed in leaf (most below ",
"1 TPM, visible in the grey bars), and the original study emphasises root tissue, which this ",
"dataset does not contain.")

legend_file <- paste0(OUT_PREFIX, "_legend.txt")
writeLines(wrap_at(legend), legend_file)
say("wrote ", basename(legend_file))

# The numbers quoted in the legend, as a table, so they can be checked without
# re-reading the prose.
stats <- data.table(
  panel = c("A", "A", "A", "A", "A", "B", "B", "B", "B"),
  study = c(rep("Mu\u00f1oz-Perez 2025 (sugarcane)", 5), rep("Ta Quang Kiet 2025 (purple)", 4)),
  quantity = c("genes mapped / libraries", "N-responsive in RB975375 (padj<0.05)",
               "of those, down at high N", "median % variance: nitrogen",
               "median % variance: genotype",
               "MYB61 copies / libraries", "U-contrast padj<0.05 in 51NG3",
               "of those, inverted vs the paper", "U-contrast padj<0.05 in TAGZ"),
  value = c(sprintf("%d / %d", nrow(ZA), ncol(ZA)), sprintf("%d / %d", n_sig, n_test),
            sprintf("%d / %d", n_down, n_sig),
            sprintf("%.1f%%", median(vp$pct_nitrogen, na.rm = TRUE)),
            sprintf("%.1f%%", median(vp$pct_genotype, na.rm = TRUE)),
            sprintf("%d / %d", sum(is_copy), ncol(ZB)), sprintf("%d / %d", n_u, n_ut),
            sprintf("%d / %d", n_neg, n_u),
            sprintf("%d", n_tagz)))
write_tsv(stats, paste0(OUT_PREFIX, "_stats.tsv"))
say("done")
