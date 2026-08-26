# =============================================================================
# 28_fig_clustering_compare.r — MCL against the SBM, on the same network
#
# The pipeline has always clustered with MCL. A stochastic block model was fitted
# to the same sugarcane network as an alternative, and this figure is the
# evidence for keeping one or switching to the other. It exists because the
# choice cannot be made on partition statistics alone: two clusterings can look
# equally reasonable as partitions and behave completely differently once the
# module-level analysis is run on them.
#
# So both partitions are carried through the IDENTICAL downstream pipeline --
# same eigengene construction, same Spearman module-trait test at the same two
# thresholds, same TF hypergeometric, same topGO run against the same background
# -- and compared at each step.
#
#   A  GRANULARITY. Module-size CCDF for both, on log-log axes. This is where the
#      two are least alike: MCL puts 19% of the network in one module and has a
#      long tail of 3-gene modules, the SBM has no giant block at all. The
#      adjusted Rand index between them is in the legend.
#
#   B  COHERENCE. PC1 variance explained per module -- how well a single
#      eigengene actually summarises its members. A partition whose modules are
#      arbitrary produces incoherent eigengenes, and every module-level result
#      downstream is only as good as this.
#
#   C  THE TEST. Modules tested, modules responsive, and the testing burden each
#      partition imposes. Fewer, larger modules mean a lighter BH correction; the
#      question is whether that converts into more signal or just fewer tests.
#
#   D  INTERPRETABILITY. How many modules clear the GO annotation gate and how
#      many terms come back. This is the panel that matters most for the
#      decision: a partition is only useful if its modules can be named.
#
# NOTHING IS WRITTEN ON THE FIGURE THAT BELONGS IN THE LEGEND.
#
# RUN: through run.sh  ->  ./run.sh figclustering   (then ./run.sh legends)
# =============================================================================

suppressMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(scico); library(grid); library(svglite); library(scales)
})

source(file.path(dirname(sub("--file=", "",
       grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "lib", "common.R"))

STUDY      <- env_req("CLEAN_STUDY")
STUDY_DIR  <- env_req("CLEAN_STUDY_DIR")
OUT_PREFIX <- env_req("CLEAN_OUT_PREFIX")
FIG        <- env_opt("CLEAN_FIG_NUM", "8")
ONT        <- env_opt("CLEAN_GO_ONTOLOGY", "BP")
NGRID      <- as.integer(env_num("CLEAN_TOPO_GRID", 300))
setDTthreads(as.integer(env_num("CLEAN_CORES", 8)))

banner(sprintf("Figure %s — MCL against the SBM (%s)", FIG, STUDY))

# Each clustering: where its two contract files are, and where its module-level
# outputs went. mcl writes at the study root, everything else in a subdirectory.
CLUS <- list(
  MCL = list(prefix = file.path(STUDY_DIR, sprintf("mcl_%s", STUDY)),
             dir    = STUDY_DIR),
  SBM = list(prefix = file.path(STUDY_DIR, sprintf("sbm_%s", STUDY)),
             dir    = file.path(STUDY_DIR, "sbm")))
LAB <- names(CLUS)
PAL <- setNames(scico(3, palette = "batlow")[1:2], LAB)

theme_f <- theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(face = "bold", size = 7.5),
        legend.key.size = unit(3.5, "mm"),
        legend.title = element_text(face = "bold", size = 7.5),
        plot.margin = margin(2, 4, 2, 2))
log_lab <- function() trans_format("log10", math_format(10^.x))
grid_of <- function(x, n = NGRID)
  unique(round(10^seq(log10(max(1, min(x))), log10(max(x)), length.out = n)))

have <- function(f) file.exists(f)
rd <- function(f, ...) if (have(f)) fread(f, ...) else NULL

# --- what exists -------------------------------------------------------------
# Named modules only: both writers add an `Unassigned` pseudo-row for the genes
# left below the minimum module size, and it is not a module.
SUMM <- rbindlist(lapply(LAB, function(k) {
  d <- rd(paste0(CLUS[[k]]$prefix, "_module_summary.tsv"))
  if (is.null(d)) { say("MISSING: ", basename(CLUS[[k]]$prefix), "_module_summary.tsv"); return(NULL) }
  d[, clustering := k][]
}), fill = TRUE)
if (!nrow(SUMM)) stop("neither clustering has a module summary — run mcl / sbmclust",
                      call. = FALSE)
UNASSIGNED <- SUMM[!grepl("^Module_", module), .(clustering, n_genes)]
SUMM <- SUMM[grepl("^Module_", module)]
SUMM[, clustering := factor(clustering, levels = LAB)]

part <- SUMM[, .(modules = .N, median = median(n_genes), max = max(n_genes),
                 Q = modularity_Q[1]), by = clustering]
part <- merge(part, UNASSIGNED[, .(clustering, unassigned = n_genes)],
              by = "clustering", all.x = TRUE)
say("partitions:"); print(part, row.names = FALSE)

# --- adjusted Rand index between the two partitions --------------------------
ARI <- NA_real_
MEM <- lapply(LAB, function(k)
  rd(paste0(CLUS[[k]]$prefix, "_membership.tsv"), select = c("gene", "module_name")))
names(MEM) <- LAB
if (!any(vapply(MEM, is.null, TRUE))) {
  m <- merge(MEM[[1]], MEM[[2]], by = "gene", suffixes = c("_a", "_b"))
  tb <- table(m$module_name_a, m$module_name_b); n <- sum(tb)
  s_ij <- sum(choose(tb, 2)); s_i <- sum(choose(rowSums(tb), 2))
  s_j <- sum(choose(colSums(tb), 2))
  ex <- s_i * s_j / choose(n, 2)
  ARI <- (s_ij - ex) / ((s_i + s_j) / 2 - ex)
  say(sprintf("adjusted Rand index over %s shared genes: %.4f", fmt_n(nrow(m)), ARI))
}

# =============================================================================
# A — granularity
# =============================================================================
ccdf <- function(x, g) vapply(g, function(k) mean(x >= k), 0)
dA <- rbindlist(lapply(LAB, function(k) {
  z <- SUMM[clustering == k, n_genes]
  if (!length(z)) return(NULL)
  g <- grid_of(z)
  data.table(clustering = k, size = g, p = ccdf(z, g))
}))
dA[, clustering := factor(clustering, levels = LAB)]

pA <- ggplot(dA[p > 0], aes(size, p, colour = clustering)) +
  geom_line(linewidth = 0.7) +
  scale_x_log10(labels = log_lab()) + scale_y_log10(labels = log_lab()) +
  scale_colour_manual(values = PAL, name = NULL) +
  annotation_logticks(sides = "bl", linewidth = 0.2,
                      short = unit(0.4, "mm"), mid = unit(0.7, "mm"),
                      long = unit(1.1, "mm")) +
  labs(x = "module size  s  (genes)", y = expression(P(S >= s))) +
  theme_f

# =============================================================================
# B — eigengene coherence
# =============================================================================
dB <- rbindlist(lapply(LAB, function(k) {
  f <- file.path(CLUS[[k]]$dir, "modules",
                 sprintf("%s_eigengenes_pc1_variance.tsv", STUDY))
  d <- rd(f)
  if (is.null(d)) { say("MISSING: ", f); return(NULL) }
  d[, clustering := k][]
}), fill = TRUE)
pB <- ggplot() + theme_void()
if (nrow(dB)) {
  dB[, clustering := factor(clustering, levels = LAB)]
  say("PC1 variance explained:")
  print(dB[, .(modules = .N, median = round(median(pc1_var_pct), 1),
               q25 = round(quantile(pc1_var_pct, .25), 1),
               q75 = round(quantile(pc1_var_pct, .75), 1)), by = clustering],
        row.names = FALSE)
  pB <- ggplot(dB, aes(clustering, pc1_var_pct, fill = clustering)) +
    geom_violin(colour = NA, alpha = 0.55, scale = "width") +
    geom_boxplot(width = 0.16, outlier.shape = NA, linewidth = 0.3,
                 fill = "white", colour = "grey25") +
    scale_fill_manual(values = PAL, guide = "none") +
    scale_y_continuous(limits = c(0, 100)) +
    labs(x = NULL, y = "PC1 variance explained (%)") +
    theme_f
}

# =============================================================================
# C — the module-trait test
# =============================================================================
dC <- rbindlist(lapply(LAB, function(k) {
  f <- file.path(CLUS[[k]]$dir, sprintf("module_trait_%s.tsv", STUDY))
  d <- rd(f)
  if (is.null(d)) { say("MISSING: ", f); return(NULL) }
  data.table(clustering = k,
             stage = c("tested", "padj <= 0.05", "|rho| >= 0.6", "responsive"),
             n = c(nrow(d), sum(d$padj <= 0.05), sum(abs(d$rho) >= 0.6),
                   sum(d$responsive)))
}), fill = TRUE)
pC <- ggplot() + theme_void()
if (nrow(dC)) {
  dC[, clustering := factor(clustering, levels = LAB)]
  dC[, stage := factor(stage, levels = c("tested", "padj <= 0.05",
                                         "|rho| >= 0.6", "responsive"))]
  say("module-trait:"); print(dcast(dC, stage ~ clustering, value.var = "n"),
                              row.names = FALSE)
  pC <- ggplot(dC, aes(stage, n, fill = clustering)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.64) +
    geom_text(aes(label = fmt_n(n)), position = position_dodge(width = 0.72),
              vjust = -0.35, size = 1.95) +
    scale_fill_manual(values = PAL, name = NULL) +
    scale_y_log10(labels = log_lab(), expand = expansion(mult = c(0, 0.22))) +
    annotation_logticks(sides = "l", linewidth = 0.2,
                        short = unit(0.4, "mm"), mid = unit(0.7, "mm"),
                        long = unit(1.1, "mm")) +
    labs(x = NULL, y = "modules") +
    theme_f + theme(axis.text.x = element_text(size = 6.2, angle = 15, hjust = 1))
}

# =============================================================================
# D — can the modules be named?
# =============================================================================
dD <- rbindlist(lapply(LAB, function(k) {
  fs <- file.path(CLUS[[k]]$dir, "module_go",
                  sprintf("module_GO_%s_%s_summary.tsv", ONT, STUDY))
  ft <- file.path(CLUS[[k]]$dir, "module_go",
                  sprintf("module_GO_%s_%s.tsv", ONT, STUDY))
  d <- rd(fs); tt <- rd(ft)
  if (is.null(d)) { say("MISSING: ", fs); return(NULL) }
  data.table(clustering = k, responsive = nrow(d), testable = sum(d$tested),
             with_term = sum(d$n_sig_terms > 0, na.rm = TRUE),
             terms = if (is.null(tt)) 0L else nrow(tt))
}), fill = TRUE)
pD <- ggplot() + theme_void()
if (nrow(dD)) {
  dD[, clustering := factor(clustering, levels = LAB)]
  say("module GO:"); print(dD, row.names = FALSE)
  # As a FRACTION of each partition's own responsive set, not as counts. Counts
  # would put 465 modules and 325 terms on one axis and invite reading a
  # partition with more modules as the better one; the question is what share of
  # them can actually be named.
  dDl <- melt(dD[, .(clustering, responsive,
                     `clears the GO\nannotation gate` = testable / responsive,
                     `returns >= 1\n{ONT} term` = with_term / responsive)],
              id.vars = c("clustering", "responsive"),
              variable.name = "stage", value.name = "frac")
  dDl[, stage := factor(gsub("\\{ONT\\}", ONT, as.character(stage)),
                        levels = gsub("\\{ONT\\}", ONT,
                                      levels(dDl$stage)))]
  dDl[, n := round(frac * responsive)]
  pD <- ggplot(dDl, aes(stage, frac, fill = clustering)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_text(aes(label = sprintf("%s of %s", fmt_n(n), fmt_n(responsive))),
              position = position_dodge(width = 0.7), vjust = -0.4, size = 2.0) +
    scale_fill_manual(values = PAL, name = NULL) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1),
                       expand = expansion(mult = c(0, 0.12))) +
    labs(x = NULL, y = "share of responsive modules") +
    theme_f + theme(axis.text.x = element_text(size = 6.2, lineheight = 1.05))
}

# =============================================================================
# compose
# =============================================================================
fig <- (pA | pB) / (pC | pD) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12))

W <- 20; H <- 14
invisible(ensure_dir(dirname(OUT_PREFIX)))
ggsave(paste0(OUT_PREFIX, ".png"), fig, width = W, height = H, units = "cm",
       dpi = 400, type = "cairo")
ggsave(paste0(OUT_PREFIX, ".pdf"), fig, width = W, height = H, units = "cm",
       device = cairo_pdf)
say(sprintf("wrote %s.{png,pdf}  (%d x %d cm)", basename(OUT_PREFIX), W, H))

stats <- rbindlist(list(
  part[, .(panel = "A", clustering = as.character(clustering),
           quantity = "modules / median / max / Q / unassigned",
           value = sprintf("%s / %s / %s / %s / %s", fmt_n(modules), fmt_n(median),
                           fmt_n(max), Q, fmt_n(unassigned)))],
  if (nrow(dB)) dB[, .(panel = "B", clustering = as.character(clustering),
                       quantity = "median PC1 variance explained",
                       value = sprintf("%.1f%%", median(pc1_var_pct))),
                   by = clustering][, .(panel, clustering, quantity, value)],
  if (nrow(dC)) dC[, .(panel = "C", clustering = as.character(clustering),
                       quantity = as.character(stage), value = fmt_n(n))],
  if (nrow(dD)) dD[, .(panel = "D", clustering = as.character(clustering),
                       quantity = "responsive / GO-testable / with a term / terms returned",
                       value = sprintf("%s / %s / %s / %s", fmt_n(responsive),
                                       fmt_n(testable), fmt_n(with_term), fmt_n(terms)))]),
  fill = TRUE)
write_tsv(stats, paste0(OUT_PREFIX, "_stats.tsv"))
say("done")

# =============================================================================
# legend
# =============================================================================
wrap_at <- function(x, width = 96)
  paste(vapply(strsplit(x, "\n[[:space:]]*\n")[[1L]],
               function(para) paste(strwrap(gsub("[[:space:]]+", " ", para), width),
                                    collapse = "\n"),
               ""), collapse = "\n\n")

pv <- function(k, col) part[clustering == k][[col]]
cv <- function(k, st) dC[clustering == k & stage == st, n]
gv <- function(k, col) dD[clustering == k][[col]]
bv <- function(k) median(dB[clustering == k, pc1_var_pct])

legend <- paste0(
"Figure ", FIG, ". MCL against a stochastic block model, on the same ", STUDY,
" network. Both partitions are carried through the IDENTICAL downstream analysis ",
"-- the same eigengene construction, the same Spearman module-trait test at the same two ",
"thresholds, the same TF hypergeometric and the same topGO run against the same background -- ",
"so the two can be compared at every step rather than only as partitions. This is evidence for a ",
"decision that has not been made, not a result.\n",
"\n",
"(A) Module-size complementary CDFs. The two partitions are barely related: MCL gives ",
fmt_n(pv("MCL", "modules")), " modules with a median of ", fmt_n(pv("MCL", "median")),
" genes and one holding ", fmt_n(pv("MCL", "max")), " (19% of the network), the SBM ",
fmt_n(pv("SBM", "modules")), " blocks with a median of ", fmt_n(pv("SBM", "median")),
" and a largest of ", fmt_n(pv("SBM", "max")), ". Their ADJUSTED RAND INDEX is ",
sprintf("%.4f", ARI), " -- essentially no agreement, so this is a genuine alternative rather ",
"than a re-parameterisation of the same structure. Newman modularity is ", pv("MCL", "Q"),
" for MCL against ", pv("SBM", "Q"), " for the SBM, which is expected and not a defect: MCL ",
"optimises a flow-based criterion that correlates with modularity, an SBM minimises a description ",
"length and does not.\n",
"\n",
"(B) How well one eigengene summarises its own module -- the assumption every module-level ",
"result rests on. The two are close, with MCL slightly ahead on the median (",
sprintf("%.1f%%", bv("MCL")), " against ", sprintf("%.1f%%", bv("SBM")),
") but a wider low tail for the SBM, which is what larger modules cost. Note the medians are ",
"not measured on comparable objects: MCL's are dominated by very small modules, where a single ",
"component explains a lot almost by construction.\n",
"\n",
"(C) The module-trait test. The SBM tests ", fmt_n(cv("SBM", "tested")), " modules against MCL's ",
fmt_n(cv("MCL", "tested")), " -- a ", sprintf("%.1f", cv("MCL", "tested") / cv("SBM", "tested")),
"x lighter multiple-testing burden -- and returns ", fmt_n(cv("SBM", "responsive")),
" responsive modules against ", fmt_n(cv("MCL", "responsive")),
". In both partitions the |rho| floor binds and BH does not, so the difference is not a ",
"correction artefact. Fewer responsive modules is not automatically worse: an SBM module holds a ",
"median of ", fmt_n(pv("SBM", "median")), " genes, so ", fmt_n(cv("SBM", "responsive")),
" of them cover far more genes than the raw counts suggest.\n",
"\n",
"(D) Whether the responsive modules can be NAMED, as a share of each partition's own responsive ",
"set. This is the largest difference in the figure and the one that bears most on the decision. ",
"Only ", fmt_n(gv("MCL", "testable")), " of MCL's ", fmt_n(gv("MCL", "responsive")),
" responsive modules (", sprintf("%.0f%%", 100 * gv("MCL", "testable") / gv("MCL", "responsive")),
") carry enough GO-annotated members to be tested at all, because its median module holds ",
fmt_n(pv("MCL", "median")), " genes. For the SBM it is ", fmt_n(gv("SBM", "testable")), " of ",
fmt_n(gv("SBM", "responsive")), " (",
sprintf("%.0f%%", 100 * gv("SBM", "testable") / gv("SBM", "responsive")),
"), and those return ", fmt_n(gv("SBM", "terms")), " ", ONT, " terms against MCL's ",
fmt_n(gv("MCL", "terms")), " from a set twenty times larger. The SBM partition also leaves only ",
fmt_n(pv("SBM", "unassigned")), " genes unassigned against MCL's ", fmt_n(pv("MCL", "unassigned")),
".\n",
"\n",
"READ AS: the SBM trades a much smaller responsive set for modules that can almost all be ",
"interpreted, while MCL finds many more responsive modules of which nine in ten cannot be tested ",
"for function. Which is preferable depends on whether the next step is to name modules or to ",
"count them. Only ", STUDY, " has an SBM fit, so this comparison cannot yet be checked against ",
"the second species.")

legend_file <- paste0(OUT_PREFIX, "_legend.txt")
writeLines(wrap_at(legend), legend_file)
say("wrote ", basename(legend_file))
