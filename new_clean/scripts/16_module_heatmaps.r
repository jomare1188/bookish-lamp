# =============================================================================
# 16_module_heatmaps.r — expression heatmaps for nitrogen-responsive modules
#
# Draws the member genes of selected modules across the study's samples, so a
# module called "responsive" can be looked at rather than taken on faith. This is
# the step where an `mi_only` module either visibly does something a correlation
# could not have caught -- a non-monotone dose response, a dispersion change --
# or reveals itself as an artifact.
#
# PER-GENE Z-SCORE, one panel. The pattern across samples is what a module
# heatmap is read for, and an absolute-VST panel alongside it mostly duplicated
# the row labels while doubling the figure width.
#
# The trade-off that costs: a z-score rescales a gene varying by 0.01 VST units
# to look exactly as structured as one varying by 5. Nothing in the figure now
# distinguishes them. Absolute expression per gene remains available in the VST
# export if a specific gene's magnitude matters.
#
# The colour ramp is PERCENTILE-CLIPPED at the 98th percentile of |z|; a few
# extreme cells otherwise flatten the whole panel. colorRamp2 clamps
# out-of-range values to the endpoints.
#
# Rendering is adapted from 11_readouts/module20/05_module20_heatmaps.r, whose
# plotting half is good. Three changes: expression comes from read_vst() (the
# matrix the network, the modules and the eigengenes were all computed from,
# rather than TPM); the sample annotation is built from the config trait spec
# instead of a hardcoded per-species if/else; and the gene set is a module rather
# than Module-20's tables.
#
# LARGE MODULES ARE SUBSET. Sugarcane's biggest module holds 19,604 genes and
# purple's 47,887 -- no heatmap renders that, and no reader reads it. Those are
# cut to the top CLEAN_HEATMAP_MAX_GENES by intramodular strength, the `strength`
# column MCL already wrote into the membership table, and the subsetting is
# stated in the panel subtitle.
#
# RUN: through run.sh  ->  ./run.sh moduleheatmap sugarcane
# =============================================================================

suppressMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(scico)
  library(grid)
})

source(file.path(dirname(sub("--file=", "",
       grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "lib", "common.R"))

STUDY      <- env_req("CLEAN_STUDY")
VST_PREFIX <- env_req("CLEAN_VST_PREFIX")
PROFILE    <- env_req("CLEAN_MODULE_PROFILE")
MEMBERSHIP <- env_req("CLEAN_MEMBERSHIP")
META_FILE  <- env_req("CLEAN_META")
TRAIT_SPEC <- env_req("CLEAN_TRAITS")
TF_FILE    <- env_opt("CLEAN_TF_FILE")
OUT_DIR    <- ensure_dir(env_req("CLEAN_OUT_DIR"))
TOP_N      <- as.integer(env_num("CLEAN_HEATMAP_TOP_N", 20))
MAX_GENES  <- as.integer(env_num("CLEAN_HEATMAP_MAX_GENES", 100))
ONLY       <- env_opt("CLEAN_HEATMAP_MODULES")   # optional explicit module list
setDTthreads(as.integer(env_num("CLEAN_CORES", 100)))

CELL_W <- unit(3.2, "mm"); CELL_H <- unit(2.6, "mm")

banner(paste("module heatmaps:", STUDY))

# --- trait spec, reusing 07's parser ----------------------------------------
parse_traits <- function(spec) {
  out <- list()
  for (block in strsplit(spec, ";", fixed = TRUE)[[1L]]) {
    block <- trimws(block); if (!nzchar(block)) next
    nm <- sub(":.*$", "", block); body <- sub("^[^:]*:", "", block)
    kv <- strsplit(strsplit(body, ",", fixed = TRUE)[[1L]], "=", fixed = TRUE)
    v <- as.numeric(vapply(kv, `[`, "", 2L))
    names(v) <- trimws(vapply(kv, `[`, "", 1L))
    out[[trimws(nm)]] <- v
  }
  out
}
TRAITS <- parse_traits(TRAIT_SPEC)

# --- data --------------------------------------------------------------------
vst <- read_vst(VST_PREFIX)
prof <- fread(PROFILE)
mem <- fread(MEMBERSHIP, select = c("gene", "module_name", "strength"))
mem[, gene := strip_version(gene)]
setnames(mem, "module_name", "module")

meta <- fread(META_FILE); setnames(meta, tolower(names(meta)))
m <- meta[match(colnames(vst), sample)]
if (anyNA(m$sample)) stop("VST samples missing from ", basename(META_FILE), call. = FALSE)
for (tr in names(TRAITS))
  if (!tr %in% names(m)) stop("trait column '", tr, "' not in the samplesheet",
                              call. = FALSE)

# order samples by the design so the annotation reads left-to-right
ord <- order(m[[names(TRAITS)[1]]], m[[names(TRAITS)[2]]], m$sample)
m <- m[ord]; vst <- vst[, m$sample, drop = FALSE]

tf_genes <- character(0)
if (nzchar(TF_FILE) && file.exists(TF_FILE))
  tf_genes <- unique(strip_version(fread(TF_FILE)$gene))

# --- which modules -----------------------------------------------------------
if (nzchar(ONLY)) {
  sel <- prof[module %chin% strsplit(trimws(ONLY), "[[:space:],]+")[[1L]]]
} else {
  sel <- prof[finding != "neither"]
  setorder(sel, padj, pearson_padj)
  # take the top N of EACH response class, so the non-linear ones are not
  # crowded out by the far more numerous linear ones
  sel <- sel[, head(.SD, TOP_N), by = finding]
}
say("modules to draw: ", nrow(sel), "  (",
    paste(sprintf("%s %d", sel[, .N, by = finding]$finding,
                  sel[, .N, by = finding]$N), collapse = " | "), ")")
if (!nrow(sel)) { say("nothing to draw"); quit(save = "no", status = 0) }

# --- annotation, built once --------------------------------------------------
t1 <- names(TRAITS)[1]; t2 <- names(TRAITS)[2]
f1 <- factor(m[[t1]], levels = names(TRAITS[[t1]]))
f2 <- factor(m[[t2]], levels = names(TRAITS[[t2]]))
c1 <- structure(scico(nlevels(f1) + 1, palette = "batlow")[seq_len(nlevels(f1))],
                names = levels(f1))
c2 <- structure(scico(nlevels(f2) + 1, palette = "lajolla",
                      begin = 0.25, end = 0.9)[seq_len(nlevels(f2))],
                names = levels(f2))
ann_list <- list(f1, f2); names(ann_list) <- c(t1, t2)
top_ann <- do.call(HeatmapAnnotation, c(
  ann_list,
  list(col = setNames(list(c1, c2), c(t1, t2)),
       annotation_name_gp = gpar(fontsize = 8),
       simple_anno_size = unit(3.5, "mm"))))

PAL_Z   <- rev(scico(256, palette = "roma"))

# --- draw --------------------------------------------------------------------
for (i in seq_len(nrow(sel))) {
  mod <- sel$module[i]
  g <- mem[module == mod][order(-strength)]
  n_all <- nrow(g)
  subset_note <- ""
  if (n_all > MAX_GENES) {
    g <- head(g, MAX_GENES)
    subset_note <- sprintf("  |  showing top %d of %d genes by intramodular strength",
                           MAX_GENES, n_all)
  }
  genes <- intersect(g$gene, rownames(vst))
  if (length(genes) < 2L) { say("  ", mod, ": <2 genes in the VST, skipped"); next }

  X <- vst[genes, , drop = FALSE]
  Z <- t(scale(t(X))); Z[!is.finite(Z)] <- 0

  zl <- as.numeric(quantile(abs(Z), 0.98, na.rm = TRUE)); if (zl == 0) zl <- 1
  col_z <- colorRamp2(seq(-zl, zl, length.out = 256), PAL_Z)

  left_ann <- NULL
  if (length(tf_genes)) {
    left_ann <- rowAnnotation(
      TF = fifelse(genes %chin% tf_genes, "TF", "-"),
      col = list(TF = c(TF = scico(3, palette = "roma")[1], `-` = "grey88")),
      annotation_name_gp = gpar(fontsize = 7),
      simple_anno_size = unit(3, "mm"))
  }

  r <- sel[i]
  sub <- sprintf("%s  |  %s genes, PC1 %.0f%% var  |  r = %.2f (padj %.1e), MI = %.2f (padj %.1e)%s",
                 r$finding, fmt_n(r$n_genes), r$pc1_var_pct, r$pearson,
                 r$pearson_padj, r$mi, r$padj, subset_note)

  show_names <- length(genes) <= 60

  ht <- Heatmap(Z, col = col_z, name = "z",
    cluster_rows = TRUE, cluster_columns = FALSE,   # keep the design column order
    cluster_row_slices = FALSE,
    column_split = f2,                              # split by the nitrogen axis
    top_annotation = top_ann, left_annotation = left_ann,
    show_row_names = show_names, row_names_gp = gpar(fontsize = 5),
    show_column_names = TRUE, column_names_gp = gpar(fontsize = 5),
    width = CELL_W * ncol(X), height = CELL_H * nrow(X),
    border = TRUE, row_gap = unit(0.8, "mm"), column_gap = unit(0.8, "mm"),
    column_title_gp = gpar(fontsize = 8, fontface = "bold"))

  # Device size derived from the ACTUAL body size, not a guess: body is
  # CELL_W x CELL_H per cell, everything else is fixed overhead in cm.
  body_w_cm <- 0.32 * ncol(X)
  body_h_cm <- 0.26 * nrow(X)
  w <- body_w_cm + (if (show_names) 5 else 1.5) + 6.5
  h <- body_h_cm + 8
  h <- max(h, 9)                                  # a 3-gene module still needs a title

  for (dev in c("png", "pdf")) {
    f <- file.path(OUT_DIR, sprintf("module_%s_%s.%s", mod, STUDY, dev))
    if (dev == "png") png(f, width = w, height = h, units = "cm", res = 300)
    else pdf(f, width = w / 2.54, height = h / 2.54)
    draw(ht,
         column_title = paste0(mod, " - ", STUDY, "\n", sub),
         column_title_gp = gpar(fontsize = 10, fontface = "bold"),
         merge_legends = TRUE, heatmap_legend_side = "right")
    dev.off()
  }

  say(sprintf("  %-14s %s  %d genes drawn  (%.0f x %.0f cm)",
              mod, r$finding, length(genes), w, h))
}

say("wrote ", nrow(sel), " module heatmaps to ", OUT_DIR)
say("done: ", STUDY)
