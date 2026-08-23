# =============================================================================
# 17_module_summary.r — one figure per study: every responsive module at once
#
# The per-module heatmaps (16) show the genes inside one module. This shows the
# MODULES themselves: each row is a module's eigengene across the study's
# samples, so the whole responsive set can be read in one picture -- how many
# distinct nitrogen response shapes there are, and which modules carry enough of
# their members' variance to be worth believing.
#
# ROWS are eigengenes, already z-scored by 14_module_eigengene.r, split by the
# SIGN of the module's Spearman rho: modules that rise with nitrogen above,
# modules that fall below. That replaces the old split by response class, which
# is gone with the MI and Pearson layers -- and it is the more useful split
# anyway, since the two halves are different biology rather than two ways of
# detecting the same thing. Clustering runs inside each half, not across them.
#
# ROW ANNOTATIONS carry what a reader needs to judge each row:
#   PC1 var %   how much of the module's gene-level variance the eigengene
#               actually represents. A module at 45% is a much weaker summary
#               than one at 95%, and the row should not be read the same way.
#   n genes     module size, log10 -- they span 3 to ~20,000
#   TF          whether the module is TF-enriched (hypergeometric, BH)
#
# The two studies get SEPARATE figures. Their sample sets differ (48 vs 18
# libraries), so the columns cannot be shared and forcing them onto one axis
# would be a lie about the design.
#
# RUN: through run.sh  ->  ./run.sh modulesummary sugarcane
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
EIG_PREFIX <- env_req("CLEAN_EIGENGENE_PREFIX")
PROFILE    <- env_req("CLEAN_MODULE_PROFILE")
META_FILE  <- env_req("CLEAN_META")
TRAIT_SPEC <- env_req("CLEAN_TRAITS")
OUT_PREFIX <- env_req("CLEAN_OUT_PREFIX")
MAX_MODULES <- as.integer(env_num("CLEAN_SUMMARY_MAX_MODULES", 250))
DIRS <- c("positive", "negative")
setDTthreads(as.integer(env_num("CLEAN_CORES", 100)))

banner(paste("module summary figure:", STUDY))

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
E <- read_vst(EIG_PREFIX)                       # modules x samples, z-scored
prof <- fread(PROFILE)
say("modules with an eigengene: ", fmt_n(nrow(E)))

meta <- fread(META_FILE); setnames(meta, tolower(names(meta)))
m <- meta[match(colnames(E), sample)]
if (anyNA(m$sample)) stop("eigengene samples missing from ", basename(META_FILE),
                          call. = FALSE)
t1 <- names(TRAITS)[1]; t2 <- names(TRAITS)[2]
# Column order groups REPLICATES of the same condition next to each other.
# Sorting by sample name alone interleaves the leaf positions -- sugarcane's
# names run B0_1, B_1, M_1, P_1, B0_2, ... so the three replicates of one segment
# land four columns apart, and the segment effect reads as vertical striping that
# obscures the treatment pattern the figure is for. Ordering by the split trait,
# then the other trait, then any CLEAN_HEATMAP_GROUP_BY columns (`segment` -- NOT
# the sheet's `tissue` column, which collapses base0 and base; see config.sh), then
# the sample puts replicates adjacent.
grp_cols <- env_list("CLEAN_HEATMAP_GROUP_BY", character(0))
grp_cols <- intersect(grp_cols, names(m))
ord_cols <- c(t2, t1, grp_cols, "sample")
ord <- do.call(order, lapply(ord_cols, function(cc) m[[cc]]))
m <- m[ord]; E <- E[, m$sample, drop = FALSE]

# --- which modules -----------------------------------------------------------
sel <- prof[responsive == TRUE & module %chin% rownames(E)]
if (!nrow(sel)) { say("no responsive modules — nothing to draw"); quit(save = "no") }
sel[, direction := factor(direction, levels = DIRS)]
sel <- sel[order(direction, padj, -abs(rho))]
n_all <- nrow(sel)
if (n_all > MAX_MODULES) {
  per <- max(1L, floor(MAX_MODULES / uniqueN(sel$direction)))
  sel <- sel[, head(.SD, per), by = direction]
  say(sprintf("showing %s of %s responsive modules (top %d per direction by significance)",
              fmt_n(nrow(sel)), fmt_n(n_all), per))
} else {
  say("showing all ", fmt_n(n_all), " responsive modules")
}
say("  ", paste(sprintf("%s %s", sel[, .N, by = direction]$direction,
                        fmt_n(sel[, .N, by = direction]$N)), collapse = " | "))

M <- E[sel$module, , drop = FALSE]

# --- annotations -------------------------------------------------------------
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
       simple_anno_size = unit(4, "mm"))))

# PC1 variance is the "should I believe this row" annotation, so it gets a
# barplot on a fixed 0-100 axis rather than a colour anyone has to decode.
var_col <- colorRamp2(c(40, 70, 100), scico(3, palette = "bamako"))
right_ann <- rowAnnotation(
  `PC1 var %` = anno_barplot(sel$pc1_var_pct, ylim = c(0, 100),
                             gp = gpar(fill = "grey35", col = NA),
                             width = unit(2.2, "cm"),
                             axis_param = list(gp = gpar(fontsize = 6))),
  `genes (log10)` = anno_barplot(log10(sel$n_genes),
                             gp = gpar(fill = "grey60", col = NA),
                             width = unit(1.6, "cm"),
                             axis_param = list(gp = gpar(fontsize = 6))),
  TF = ifelse(sel$tf_enriched %in% TRUE, "enriched", "-"),
  col = list(TF = c(enriched = scico(3, palette = "roma")[1], `-` = "grey88")),
  annotation_name_gp = gpar(fontsize = 7),
  simple_anno_size = unit(3, "mm"))

zl <- as.numeric(quantile(abs(M), 0.98, na.rm = TRUE)); if (zl == 0) zl <- 1
col_z <- colorRamp2(seq(-zl, zl, length.out = 256), rev(scico(256, palette = "roma")))

# Sample labels are dropped on this figure -- with 48 columns they add ~4 cm of
# rotated text for names nobody reads at this scale, and the design is already
# carried by the column annotation bars and the treatment split. That frees the
# cells to shrink.
CELL_W <- unit(2, "mm")
CELL_H <- unit(if (nrow(M) > 120) 1 else if (nrow(M) > 40) 1.8 else 3.5, "mm")
cell_h_cm <- as.numeric(gsub("mm", "", format(CELL_H))) / 10
cell_w_cm <- 0.2

ht <- Heatmap(M, col = col_z, name = "eigengene z",
  cluster_rows = TRUE, cluster_columns = FALSE,
  cluster_row_slices = FALSE,
  row_split = factor(as.character(sel$direction), levels = DIRS),
  column_split = f2,
  top_annotation = top_ann, right_annotation = right_ann,
  show_row_names = nrow(M) <= 60, row_names_gp = gpar(fontsize = 5),
  show_column_names = FALSE,
  row_title_gp = gpar(fontsize = 9, fontface = "bold"), row_title_rot = 0,
  column_title_gp = gpar(fontsize = 9, fontface = "bold"),
  width = CELL_W * ncol(M), height = CELL_H * nrow(M),
  border = TRUE, row_gap = unit(2, "mm"), column_gap = unit(1, "mm"))

# Kept short and wrapped: the figure is ~23 cm wide and a single long title line
# is clipped at the edges rather than shrunk to fit.
drawn <- if (nrow(sel) < n_all)
  sprintf("  |  %s drawn, top per direction", fmt_n(nrow(sel))) else ""
sub <- sprintf("%s of %s modules responsive%s  |  median PC1 var %.0f%%\nSpearman |rho| >= %s, padj <= %s",
               fmt_n(n_all), fmt_n(nrow(prof)), drawn, median(sel$pc1_var_pct),
               env_opt("CLEAN_MODULE_R_THR", "0.6"),
               env_opt("CLEAN_MODULE_PADJ_THR", "0.05"))

w <- cell_w_cm * ncol(M) + (if (nrow(M) <= 60) 5 else 1.5) + 12
h <- cell_h_cm * nrow(M) + 7

invisible(ensure_dir(dirname(OUT_PREFIX)))
for (dev in c("png", "pdf")) {
  f <- paste0(OUT_PREFIX, ".", dev)
  if (dev == "png") png(f, width = w, height = h, units = "cm", res = 300)
  else pdf(f, width = w / 2.54, height = h / 2.54)
  draw(ht, column_title = paste0("Nitrogen-responsive modules - ", STUDY, "\n", sub),
       column_title_gp = gpar(fontsize = 9, fontface = "bold"),
       merge_legends = TRUE, heatmap_legend_side = "right")
  dev.off()
}
say(sprintf("wrote %s.{png,pdf}  (%.0f x %.0f cm, %d modules x %d samples)",
            basename(OUT_PREFIX), w, h, nrow(M), ncol(M)))
say("done: ", STUDY)
