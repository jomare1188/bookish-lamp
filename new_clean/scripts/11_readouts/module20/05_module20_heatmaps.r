# ============================================================================
# 05_module20_heatmaps.r — publication heatmaps of every Module-20 candidate
# across all libraries of BOTH studies.
#
# One figure per species (the sample sets are different, so columns cannot be
# shared), each with two panels:
#   A  log2(TPM + 1)  — sequential scico palette: which copies are ACTUALLY
#      expressed. Without this panel a z-score heatmap makes near-zero copies
#      look as structured as real ones.
#   B  per-gene z-score — diverging scico `roma`: the SHAPE of the response.
#
# Layout: columns split by genotype then nitrogen level (so the design stays
# readable and conditions are directly comparable), rows split by Module-20
# locus; hierarchical clustering runs INSIDE each block. Genotype and nitrogen
# are also drawn as labelled annotation bars.
#
# Genes shown: every gene that passed the reciprocal-Arabidopsis filter in
# step 02, whether or not it is a node of that species' network — the
# `in_network` state is drawn as a row annotation instead of used as a filter.
#
# RUN: /home/genomics/miniconda3/envs/r_env/bin/Rscript 05_module20_heatmaps.r
# ============================================================================

suppressMessages({
  library(data.table); library(ComplexHeatmap); library(circlize)
  library(scico); library(grid)
})

base   <- "/dados04/jorge/comparative_saccharum"
outdir <- "/dados04/jorge/comparative_saccharum/new_clean/results/readouts/module20"
# Cached inputs from steps 01-02b (sequence-only, network-independent) live in
# the original GET_TFS tree; only 03-05 output goes under results/.
cache  <- file.path(base, "GET_TFS/new/results/module20")


# Column grouping. "nitrogen" splits columns by N level ONLY, pooling both
# genotypes inside each block, with genotype kept as an annotation bar. That is
# the layout that answers "does Module 20 track nitrogen or genotype?": if it
# tracks nitrogen the blocks differ; if it tracks genotype the samples separate
# by the Genotype bar WITHIN each block (column clustering runs inside blocks,
# so a genotype-driven set visibly splits into two sub-clusters).
# Use "genotype_nitrogen" for the previous genotype x N layout.
SPLIT_BY <- "nitrogen"

CELL_W <- unit(3.2, "mm")     # small cells: these heatmaps are wide
CELL_H <- unit(2.6, "mm")
PAL_ABS <- scico(256, palette = "batlow")            # sequential: expression
PAL_Z   <- scico(256, palette = "roma", direction = -1)  # diverging: z-score

tab <- fread(file.path(outdir, "module20_network_table.tsv"))
mem <- fread(file.path(cache, "module20_members.tsv"))

# ---------------------------------------------------------------------------
# expression, with genotype / nitrogen parsed from the sample sheets
# ---------------------------------------------------------------------------
load_expr <- function(tpmf, metaf, sp) {
  tpm <- fread(tpmf); setnames(tpm, 1, "gene")
  tpm[, gene := sub("\\.v[0-9.]+$", "", gene)]
  meta <- fread(metaf)
  if (sp == "sugarcane") {
    # `segment`, NOT the sheet's `tissue` column. The study sampled FOUR leaf
    # segments -- base0, base, mid, tip, 12 libraries each, carried in the
    # library names as B0/B/M/P -- and `tissue` collapses base0 and base into a
    # single "Leaf Base" of 24, calls mid "Leaf" and tip "Leaf Apex". Using it
    # here fitted the covariate on three levels instead of four.
    if (!"segment" %in% names(meta))
      stop("the sugarcane sample sheet has no `segment` column; `tissue` collapses\n",
           "  base0 and base and must not be used as the leaf-segment factor",
           call. = FALSE)
    m <- meta[, .(sample, genotype, treatment, segment)]
    m[, segment := factor(segment, levels = c("base0", "base", "mid", "tip"))]
    m[, nlev := factor(fifelse(grepl("Low", treatment), "LowN", "HighN"),
                       levels = c("LowN", "HighN"))]
    m[, genotype := fifelse(genotype == "RB975375", "RB975375\nresponsive",
                                                    "RB937570\nnon-responsive")]
  } else {
    m <- meta[tissue == "leaf", .(sample, genotype, treatment)]
    m[, segment := "leaf"]
    m[, nlev := factor(treatment, levels = c("0N", "2N", "6N"))]
  }
  m[, genotype := factor(genotype, levels = unique(genotype))]
  genes <- tab[species == sp, unique(gene)]
  samp  <- intersect(names(tpm), m$sample)
  X <- as.matrix(tpm[gene %in% genes, ..samp])
  rownames(X) <- tpm[gene %in% genes, gene]
  list(X = log2(X[, m$sample, drop = FALSE] + 1), meta = m)
}

cfg <- list(
  sugarcane = list(tpm = file.path(base, "run1/salmon/salmon.merged.gene_tpm.tsv"),
                   meta = file.path(base, "samplesheet.csv"),
                   ttl = "Muñoz (sugarcane / R570) — 48 leaf libraries"),
  purple    = list(tpm = file.path(base, "china/run2_onlyL/salmon/salmon.merged.gene_tpm.tsv"),
                   meta = file.path(base, "china/samplesheet_china.csv"),
                   ttl = "Kiet (purple / LA purple) — 18 leaf libraries")
)

# ---------------------------------------------------------------------------
# one species -> a two-panel figure
# ---------------------------------------------------------------------------
draw_species <- function(sp) {
  cc <- cfg[[sp]]
  e  <- load_expr(cc$tpm, cc$meta, sp)
  X  <- e$X; m <- e$meta

  # row metadata, in the same order as X
  rd <- unique(tab[species == sp, .(gene, locus, family, in_network, our_family)])
  rd <- rd[match(rownames(X), gene)]
  rd[, locus_lab := sub("^(N_R___|N_NR___|SCA[0-9_]+__)", "", locus)]
  rd[, mybcall := fifelse(grepl("MYB", our_family), "MYB", "not MYB")]

  # z-score per gene; genes flat across all libraries have no defined z
  Z <- t(scale(t(X)))
  Z[!is.finite(Z)] <- 0

  # ---- variance partition: does this gene set track NITROGEN or GENOTYPE? --
  # The column grouping above makes this visible; this makes it quantitative.
  # Per gene, lm(expression ~ genotype + nitrogen) and the share of the total
  # sum of squares each factor takes.
  vp <- rbindlist(lapply(rownames(X), function(g) {
    y <- X[g, ]
    if (sd(y) == 0) return(data.table(gene = g, pct_genotype = NA_real_, pct_nitrogen = NA_real_))
    a  <- anova(lm(y ~ m$genotype + m$nlev)); ss <- a[["Sum Sq"]]
    data.table(gene = g, pct_genotype = 100 * ss[1] / sum(ss),
               pct_nitrogen = 100 * ss[2] / sum(ss))
  }))
  vp[, species := sp]
  fwrite(vp, file.path(outdir, sprintf("module20_variance_partition_%s.tsv", sp)), sep = "\t")
  med_g <- median(vp$pct_genotype, na.rm = TRUE); med_n <- median(vp$pct_nitrogen, na.rm = TRUE)
  n_gt  <- sum(vp$pct_nitrogen > vp$pct_genotype, na.rm = TRUE)
  vp_lab <- sprintf("median variance explained: nitrogen %.1f%% vs genotype %.1f%%  |  nitrogen > genotype in %d/%d genes",
                    med_n, med_g, n_gt, sum(!is.na(vp$pct_genotype)))
  cat(sprintf("[%s] %s\n", sp, vp_lab))

  # ---- annotations --------------------------------------------------------
  gt_lev <- levels(m$genotype); n_lev <- levels(m$nlev)
  gt_col <- structure(scico(length(gt_lev) + 1, palette = "batlow")[seq_len(length(gt_lev))],
                      names = gt_lev)
  n_col  <- structure(scico(length(n_lev) + 1, palette = "lajolla",
                            begin = 0.25, end = 0.9)[seq_len(length(n_lev))],
                      names = n_lev)
  top_ann <- HeatmapAnnotation(
    Genotype = m$genotype, Nitrogen = m$nlev,
    col = list(Genotype = gt_col, Nitrogen = n_col),
    annotation_name_gp = gpar(fontsize = 8),
    simple_anno_size = unit(3.5, "mm"),
    annotation_legend_param = list(
      Genotype = list(title_gp = gpar(fontsize = 8, fontface = "bold"),
                      labels_gp = gpar(fontsize = 7)),
      Nitrogen = list(title_gp = gpar(fontsize = 8, fontface = "bold"),
                      labels_gp = gpar(fontsize = 7))))

  left_ann <- rowAnnotation(
    `TF call` = rd$mybcall,
    `in network` = fifelse(rd$in_network, "yes", "no"),
    col = list(`TF call` = c(MYB = scico(3, palette = "roma")[1], `not MYB` = "grey85"),
               `in network` = c(yes = "grey25", no = "grey85")),
    annotation_name_gp = gpar(fontsize = 7),
    simple_anno_size = unit(3, "mm"),
    annotation_legend_param = list(
      `TF call` = list(title_gp = gpar(fontsize = 8, fontface = "bold"),
                       labels_gp = gpar(fontsize = 7)),
      `in network` = list(title_gp = gpar(fontsize = 8, fontface = "bold"),
                          labels_gp = gpar(fontsize = 7))))

  common <- list(
    cluster_rows = TRUE, cluster_columns = TRUE,      # clustering inside blocks
    cluster_row_slices = FALSE, cluster_column_slices = FALSE,  # keep design order
    row_split = rd$locus_lab,
    column_split = if (SPLIT_BY == "nitrogen") m$nlev else m[, .(genotype, nlev)],
    row_title_gp = gpar(fontsize = 7), row_title_rot = 0,
    column_title_gp = gpar(fontsize = 8, fontface = "bold"),
    show_row_names = TRUE, row_names_gp = gpar(fontsize = 5),
    show_column_names = TRUE, column_names_gp = gpar(fontsize = 5),
    width = CELL_W * ncol(X), height = CELL_H * nrow(X),
    border = TRUE, row_gap = unit(0.8, "mm"), column_gap = unit(0.8, "mm"))

  # A few copies are far more expressed than the rest; without clipping they
  # flatten the whole panel. colorRamp2 clamps out-of-range values to the ends.
  amax <- as.numeric(quantile(X, 0.99, na.rm = TRUE))
  hA <- do.call(Heatmap, c(list(
    matrix = X, name = "log2(TPM+1)\n(99% clip)",
    col = colorRamp2(seq(0, amax, length.out = 256), PAL_ABS),
    top_annotation = top_ann, left_annotation = left_ann,
    column_title = paste0("A  absolute expression — ", cc$ttl),
    column_title_side = "top",
    heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "bold"),
                                labels_gp = gpar(fontsize = 7),
                                legend_height = unit(2.5, "cm"))), common))

  zl <- as.numeric(quantile(abs(Z), 0.98, na.rm = TRUE))
  hB <- do.call(Heatmap, c(list(
    matrix = Z, name = "z-score\n(98% clip)",
    col = colorRamp2(seq(-zl, zl, length.out = 256), PAL_Z),
    top_annotation = top_ann, left_annotation = left_ann,
    column_title = paste0("B  per-gene z-score — ",
      if (SPLIT_BY == "nitrogen")
        "columns grouped by NITROGEN only; genotype shown as a bar\n"
      else "shape of the nitrogen response\n", vp_lab),
    column_title_side = "top",
    heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "bold"),
                                labels_gp = gpar(fontsize = 7),
                                legend_height = unit(2.5, "cm"))), common))

  w <- 12 + 0.32 * ncol(X)
  h <- 4 + 0.30 * nrow(X)
  for (dev in c("png", "pdf")) {
    f <- file.path(outdir, sprintf("module20_heatmap_%s.%s", sp, dev))
    # cairo, not the base devices: these titles carry an em dash and "Muñoz",
    # and base pdf() transliterates UTF-8 -- it was silently writing "Munoz" and
    # a hyphen, with a warning per draw.
    if (dev == "png") png(f, width = w, height = 2 * h, units = "cm", res = 300,
                          type = "cairo")
    else              cairo_pdf(f, width = w / 2.54, height = 2 * h / 2.54)
    pushViewport(viewport(layout = grid.layout(2, 1)))
    pushViewport(viewport(layout.pos.row = 1))
    draw(hA, heatmap_legend_side = "right", annotation_legend_side = "right",
         merge_legend = TRUE, newpage = FALSE)
    popViewport()
    pushViewport(viewport(layout.pos.row = 2))
    draw(hB, heatmap_legend_side = "right", annotation_legend_side = "right",
         merge_legend = TRUE, newpage = FALSE)
    popViewport(2)
    dev.off()
  }
  cat(sprintf("[%s] %d genes x %d libraries -> module20_heatmap_%s.{png,pdf}\n",
              sp, nrow(X), ncol(X), sp))
  invisible(NULL)
}

for (sp in names(cfg)) draw_species(sp)
cat(sprintf("\nOutputs in %s\n", outdir))
