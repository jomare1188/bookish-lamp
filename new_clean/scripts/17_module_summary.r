# =============================================================================
# 17_module_summary.r — one figure per study: every responsive module at once
#
# The per-module heatmaps (16) show the genes inside one module. This shows the
# MODULES themselves: each row is a module's eigengene across the study's
# samples, so the whole responsive set can be read in one picture -- how many
# distinct nitrogen response shapes there are, whether the linear and non-linear
# classes look different, and which modules carry enough of their members'
# variance to be worth believing.
#
# ROWS are eigengenes, already z-scored by 14_module_eigengene.r, split by
# response class (`both` / `mi_only` / `pearson_only`). Clustering runs inside
# each class, not across them, so the classes stay visually separated.
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
ord <- order(m[[t1]], m[[t2]], m$sample)
m <- m[ord]; E <- E[, m$sample, drop = FALSE]

# --- which modules -----------------------------------------------------------
sel <- prof[finding != "neither" & module %chin% rownames(E)]
if (!nrow(sel)) { say("no responsive modules — nothing to draw"); quit(save = "no") }
setorder(sel, finding, padj, pearson_padj)
n_all <- nrow(sel)
if (n_all > MAX_MODULES) {
  per <- max(1L, floor(MAX_MODULES / uniqueN(sel$finding)))
  sel <- sel[, head(.SD, per), by = finding]
  say(sprintf("showing %s of %s responsive modules (top %d per class by significance)",
              fmt_n(nrow(sel)), fmt_n(n_all), per))
} else {
  say("showing all ", fmt_n(n_all), " responsive modules")
}
say("  ", paste(sprintf("%s %s", sel[, .N, by = finding]$finding,
                        fmt_n(sel[, .N, by = finding]$N)), collapse = " | "))

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

CELL_W <- unit(4, "mm")
CELL_H <- unit(if (nrow(M) > 120) 1.2 else if (nrow(M) > 40) 2.2 else 4, "mm")
cell_h_cm <- as.numeric(gsub("mm", "", format(CELL_H))) / 10

ht <- Heatmap(M, col = col_z, name = "eigengene z",
  cluster_rows = TRUE, cluster_columns = FALSE,
  cluster_row_slices = FALSE,
  row_split = factor(sel$finding, levels = c("both", "mi_only", "pearson_only")),
  column_split = f2,
  top_annotation = top_ann, right_annotation = right_ann,
  show_row_names = nrow(M) <= 60, row_names_gp = gpar(fontsize = 5),
  show_column_names = TRUE, column_names_gp = gpar(fontsize = 6),
  row_title_gp = gpar(fontsize = 9, fontface = "bold"), row_title_rot = 0,
  column_title_gp = gpar(fontsize = 9, fontface = "bold"),
  width = CELL_W * ncol(M), height = CELL_H * nrow(M),
  border = TRUE, row_gap = unit(2, "mm"), column_gap = unit(1, "mm"))

sub <- sprintf("%s responsive modules of %s tested  |  median PC1 var %.0f%%  |  |r| >= %s, MI >= calibrated floor",
               fmt_n(nrow(sel)), fmt_n(nrow(prof)), median(sel$pc1_var_pct),
               env_opt("CLEAN_MODULE_R_THR", "0.6"))

w <- 0.4 * ncol(M) + (if (nrow(M) <= 60) 5 else 1.5) + 12
h <- cell_h_cm * nrow(M) + 10

invisible(ensure_dir(dirname(OUT_PREFIX)))
for (dev in c("png", "pdf")) {
  f <- paste0(OUT_PREFIX, ".", dev)
  if (dev == "png") png(f, width = w, height = h, units = "cm", res = 300)
  else pdf(f, width = w / 2.54, height = h / 2.54)
  draw(ht, column_title = paste0("Nitrogen-responsive modules - ", STUDY, "\n", sub),
       column_title_gp = gpar(fontsize = 11, fontface = "bold"),
       merge_legends = TRUE, heatmap_legend_side = "right")
  dev.off()
}
say(sprintf("wrote %s.{png,pdf}  (%.0f x %.0f cm, %d modules x %d samples)",
            basename(OUT_PREFIX), w, h, nrow(M), ncol(M)))
say("done: ", STUDY)
