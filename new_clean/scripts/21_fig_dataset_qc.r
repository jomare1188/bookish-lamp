# =============================================================================
# 21_fig_dataset_qc.r — the dataset figure: what was reused, and what came back
# out of re-quantifying it.
#
# This is the paper's opening figure. Every later result rests on two claims that
# a reader has no way to check from the text alone: that the two reused studies
# have designs that can be compared at all, and that re-quantifying somebody
# else's libraries through one pipeline actually worked. This figure is the
# evidence for both, and it is the figure that has to earn the reader's trust in
# the 1.2 TB of network output that follows.
#
#   A  THE DESIGN. Libraries per study x genotype x nitrogen level (x leaf
#      segment, in sugarcane). This is where the sample-size asymmetry that runs
#      through every downstream result becomes visible: 48 against 18. It is also
#      where the two designs stop matching -- sugarcane contrasts two nitrogen
#      levels across four leaf segments, purple runs a three-point dose in one
#      tissue -- which is why the trait is coded ordinally and tested by rank.
#
#   B  SEQUENCING DEPTH AND MAPPING RATE, per library, from the salmon logs of
#      each nf-core/rnaseq run. Plotted against each other rather than as two
#      bar charts: the failure mode worth seeing is a library that is both
#      shallow and poorly mapped, and that is a position on this plane, not a
#      value on either axis.
#
#   C  THE GENE FUNNEL. Annotated genes -> quantified -> surviving the CV filter
#      -> nodes in the network. A reader who only sees "103,336 nodes" cannot
#      tell whether that is most of the annotation or a tenth of it. Shown as
#      counts and as the fraction retained, per study, because the two references
#      differ enough that the same filter keeps very different proportions.
#
#   D  PCA of the VST matrix per study, on the 2,000 most variable genes. The
#      question is whether the design is visible in the data at all before any
#      network is built: if nitrogen does not separate here, nothing downstream
#      can be believed. It is also where sugarcane's leaf-segment effect shows
#      itself as the largest axis of variation, which is the reason the module
#      heatmaps order columns by segment.
#
# NOTHING IS WRITTEN ON THE FIGURE THAT BELONGS IN THE LEGEND -- panel letters and
# the labels the data needs, nothing else. The legend is generated from the same
# variables that drew the panels and written to <prefix>_legend.txt; see
# docs/methods.md, "The rule every paper figure follows".
#
# All four panels are ggplot2 assembled with patchwork, unlike the reproduction
# figure's ComplexHeatmap panels -- these are scatter/bar/tile plots and a
# grid-viewport assembly would buy nothing here.
#
# RUN: through run.sh  ->  ./run.sh figdataset   (then ./run.sh legends)
# =============================================================================

suppressMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(scico); library(grid); library(svglite)
})

source(file.path(dirname(sub("--file=", "",
       grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "lib", "common.R"))

STUDIES    <- env_list("CLEAN_STUDIES", c("sugarcane", "purple"))
OUT_PREFIX <- env_req("CLEAN_OUT_PREFIX")
FIG        <- env_opt("CLEAN_FIG_NUM", "1")
N_PCA      <- as.integer(env_num("CLEAN_PCA_NTOP", 2000))
setDTthreads(as.integer(env_num("CLEAN_CORES", 8)))

# per-study inputs, all passed in by run.sh so nothing is hardcoded here
cfg <- lapply(STUDIES, function(st) list(
  label   = env_req(sprintf("CLEAN_LABEL_%s", toupper(st))),
  meta    = env_req(sprintf("CLEAN_META_%s", toupper(st))),
  salmon  = env_req(sprintf("CLEAN_SALMONQC_%s", toupper(st))),
  tpm     = env_req(sprintf("CLEAN_TPM_%s", toupper(st))),
  vst     = env_req(sprintf("CLEAN_VST_%s", toupper(st))),
  nodes   = env_req(sprintf("CLEAN_NODES_%s", toupper(st)))))
names(cfg) <- STUDIES

banner(sprintf("Figure %s — dataset, QC and quantification", FIG))

STUDY_LAB <- vapply(cfg, `[[`, "", "label")
# base -> tip, the order they sit on the leaf; any other order makes the design
# panel and the heatmap column blocks arbitrary.
SEGMENT_LEVELS <- c("base0", "base", "mid", "tip")
PAL_N   <- scico(4, palette = "lajolla", begin = 0.25, end = 0.9)
theme_f <- theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(face = "bold", size = 7.5),
        legend.key.size = unit(3.5, "mm"),
        legend.title = element_text(face = "bold", size = 7.5),
        plot.margin = margin(2, 4, 2, 2))

# --- metadata, harmonised across two differently-shaped sample sheets ---------
read_meta <- function(st) {
  m <- fread(cfg[[st]]$meta); setnames(m, tolower(names(m)))
  m[, study := STUDY_LAB[[st]]]
  # Sugarcane names its nitrogen levels, purple gives a dose. Ordering is what
  # matters downstream, so both become an ordered factor with the dose as label.
  if (st == "sugarcane") {
    m[, nlev := factor(fifelse(grepl("Low", treatment), "Low N", "High N"),
                       levels = c("Low N", "High N"))]
    # THE `tissue` COLUMN OF THE SUGARCANE SHEET IS NOT THE DESIGN. The study
    # sampled FOUR leaf segments -- base0, base, mid, tip, 12 libraries each,
    # carried in the library names as B0/B/M/P -- but `tissue` collapses base0
    # and base both into "Leaf Base" (24), calls mid "Leaf" and tip "Leaf Apex".
    # Use the `segment` column, which is the four real levels.
    if (!"segment" %in% names(m))
      stop("the sugarcane sample sheet has no `segment` column; `tissue` collapses\n",
           "  base0 and base and must not be used as the leaf-segment factor",
           call. = FALSE)
    m[, seg := factor(segment, levels = SEGMENT_LEVELS)]
  } else {
    m[, nlev := factor(treatment, levels = c("0N", "2N", "6N"))]
    m[, seg := factor("leaf")]
  }
  m[, .(sample, genotype, treatment, study, nlev, seg)]
}
META <- rbindlist(lapply(STUDIES, read_meta), fill = TRUE)
META[, study := factor(study, levels = unname(STUDY_LAB))]

# =============================================================================
# A — the design
# =============================================================================
des <- META[, .(n = .N), by = .(study, genotype, nlev, seg)]
say("panel A: design cells")
print(META[, .N, by = .(study, genotype, nlev)][order(study, genotype, nlev)],
      row.names = FALSE)

# facet_GRID would draw a cell for every study x genotype pair including the six
# that do not exist (no purple genotype is in the sugarcane study and vice
# versa); facet_wrap makes only the facets the data actually has. Fill carries
# the study rather than the count -- with every cell at 3 or 6 a continuous
# colour bar is decoration, and the count is printed anyway.
des[, panel := factor(sprintf("%s  |  %s", study, genotype),
                      levels = unique(sprintf("%s  |  %s", study, genotype)[order(study, genotype)]))]
pA <- ggplot(des, aes(nlev, seg, fill = study)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = n), size = 2.7, colour = "white", fontface = "bold") +
  facet_wrap(~ panel, scales = "free", nrow = 2) +
  scale_fill_manual(values = scico(3, palette = "batlow")[1:2], guide = "none") +
  labs(x = "nitrogen supplied", y = "leaf segment") +
  theme_f + theme(panel.grid = element_blank())

# =============================================================================
# B — depth vs mapping rate
# =============================================================================
qc <- rbindlist(lapply(STUDIES, function(st) {
  q <- fread(cfg[[st]]$salmon,
             select = c("Sample", "num_processed", "num_mapped", "percent_mapped"))
  setnames(q, "Sample", "sample")
  q[, study := STUDY_LAB[[st]]]
  q
}))
qc <- merge(qc, META[, .(sample, nlev, genotype, study)], by = c("sample", "study"))
qc[, study := factor(study, levels = unname(STUDY_LAB))]
qc[, m_proc := num_processed / 1e6]

qc_sum <- qc[, .(libs = .N,
                 proc_med = median(m_proc), proc_min = min(m_proc), proc_max = max(m_proc),
                 map_med = median(percent_mapped), map_min = min(percent_mapped),
                 map_max = max(percent_mapped)), by = study]
say("panel B: per-library depth and mapping")
print(qc_sum, row.names = FALSE)

# the worst library in each study, named, because panel B exists to show it
worst <- qc[, .SD[which.min(percent_mapped)], by = study]
say("  lowest mapping rate per study: ",
    paste(sprintf("%s %s %.1f%% (%.1f M fragments)", worst$study, worst$sample,
                  worst$percent_mapped, worst$m_proc), collapse = " | "))

pB <- ggplot(qc, aes(m_proc, percent_mapped, fill = nlev, shape = study)) +
  geom_point(size = 1.9, colour = "grey25", stroke = 0.25, alpha = 0.95) +
  scale_shape_manual(values = c(21, 24), name = NULL) +
  scale_fill_manual(values = setNames(PAL_N[c(1, 3, 2, 4)][seq_len(nlevels(META$nlev))],
                                      levels(META$nlev)),
                    name = "nitrogen",
                    guide = guide_legend(override.aes = list(shape = 21))) +
  scale_x_continuous(name = "fragments processed (millions)") +
  scale_y_continuous(name = "mapped to the reference (%)",
                     limits = c(min(30, floor(min(qc$percent_mapped))), 90)) +
  theme_f

# =============================================================================
# C — the gene funnel
# =============================================================================
fun <- rbindlist(lapply(STUDIES, function(st) {
  tpm  <- fread(cfg[[st]]$tpm, select = 1L)
  ann  <- nrow(tpm)
  meta_f <- paste0(cfg[[st]]$vst, ".meta.json")
  txt  <- paste(readLines(meta_f, warn = FALSE), collapse = " ")
  kept <- as.integer(sub('.*"n_genes"[[:space:]]*:[[:space:]]*([0-9]+).*', "\\1", txt))
  cv   <- as.numeric(sub('.*"min_cv"[[:space:]]*:[[:space:]]*([0-9.]+).*', "\\1", txt))
  nod  <- nrow(fread(cfg[[st]]$nodes, select = 1L))
  data.table(study = STUDY_LAB[[st]], min_cv = cv,
             stage = factor(c("annotated &\nquantified",
                              sprintf("CV >= %g%%\nretained", cv),
                              "in the\nnetwork"),
                            levels = c("annotated &\nquantified",
                                       sprintf("CV >= %g%%\nretained", cv),
                                       "in the\nnetwork")),
             n = c(ann, kept, nod))
}))
fun[, study := factor(study, levels = unname(STUDY_LAB))]
fun[, frac := n / n[1], by = study]
say("panel C: gene funnel")
print(fun[, .(study, stage = gsub("\n", " ", stage), n, pct = round(100 * frac, 1))],
      row.names = FALSE)

pC <- ggplot(fun, aes(stage, n / 1000, fill = study)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.66) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * frac)),
            position = position_dodge(width = 0.72), vjust = -0.45, size = 2.3) +
  scale_fill_manual(values = scico(3, palette = "batlow")[1:2], name = NULL) +
  scale_y_continuous(name = "genes (thousands)",
                     expand = expansion(mult = c(0, 0.14))) +
  labs(x = NULL) + theme_f

# =============================================================================
# D — PCA per study
# =============================================================================
pca <- rbindlist(lapply(STUDIES, function(st) {
  V <- read_vst(cfg[[st]]$vst)
  # rowVars without a matrixStats dependency: E[x^2] - E[x]^2 is numerically
  # fine here because VST values are O(10) and the matrix is float32 anyway.
  mu <- rowMeans(V)
  v  <- rowMeans(V * V) - mu * mu
  keep <- head(order(v, decreasing = TRUE), min(N_PCA, nrow(V)))
  p <- prcomp(t(V[keep, , drop = FALSE]), center = TRUE, scale. = FALSE)
  ve <- 100 * p$sdev^2 / sum(p$sdev^2)
  d <- data.table(sample = rownames(p$x), PC1 = p$x[, 1], PC2 = p$x[, 2],
                  study = STUDY_LAB[[st]], ve1 = ve[1], ve2 = ve[2])
  # What does each axis TRACK? Assuming is how a legend acquires a false claim,
  # so it is measured: R^2 of each component on each design factor.
  m <- fread(cfg[[st]]$meta); setnames(m, tolower(names(m)))
  m <- m[match(d$sample, sample)]
  # `segment`, not `tissue` -- see read_meta(). On the collapsed 3-level `tissue`
  # column this R^2 came out at 0.450; on the real four segments it is higher,
  # and it is the number the legend quotes.
  for (fac in c("genotype", "treatment", "segment")) {
    f <- factor(m[[fac]])
    if (nlevels(f) < 2) next
    r1 <- summary(lm(d$PC1 ~ f))$r.squared; r2 <- summary(lm(d$PC2 ~ f))$r.squared
    say(sprintf("  %-9s %-10s R2 on PC1 = %.3f, PC2 = %.3f", STUDY_LAB[[st]], fac, r1, r2))
    d[[paste0("r2_pc1_", fac)]] <- r1
    d[[paste0("r2_pc2_", fac)]] <- r2
  }
  d
}), fill = TRUE)      # purple has one tissue, so it has no tissue R^2 columns
pca <- merge(pca, META[, .(sample, nlev, genotype, seg, study)],
             by = c("sample", "study"))
pca[, study := factor(study, levels = unname(STUDY_LAB))]
pca[, facet := sprintf("%s\nPC1 %.0f%% | PC2 %.0f%%", study, ve1, ve2)]
pca[, facet := factor(facet, levels = unique(facet[order(study)]))]
say("panel D: PCA variance explained")
print(unique(pca[, .(study, PC1 = round(ve1, 1), PC2 = round(ve2, 1))]),
      row.names = FALSE)

pD <- ggplot(pca, aes(PC1, PC2, fill = nlev, shape = genotype)) +
  geom_point(size = 2.1, colour = "grey25", stroke = 0.25) +
  facet_wrap(~ facet, scales = "free", nrow = 1) +
  scale_shape_manual(values = c(21, 24, 22, 23), name = "genotype") +
  scale_fill_manual(values = setNames(PAL_N[c(1, 3, 2, 4)][seq_len(nlevels(META$nlev))],
                                      levels(META$nlev)),
                    name = "nitrogen",
                    guide = guide_legend(override.aes = list(shape = 21))) +
  theme_f

# =============================================================================
# compose
# =============================================================================
fig <- (pA | pB) / (pC | pD) +
  plot_layout(heights = c(1, 1), widths = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12))

W <- 24; H <- 15
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
# Generated from the variables that drew the panels, so a number cannot disagree
# between figure and legend. Assembled by ./run.sh legends.
wrap_at <- function(x, width = 96)
  paste(vapply(strsplit(x, "\n[[:space:]]*\n")[[1L]],
               function(para) paste(strwrap(gsub("[[:space:]]+", " ", para), width),
                                    collapse = "\n"),
               ""), collapse = "\n\n")

g <- function(st, col) qc_sum[study == STUDY_LAB[[st]]][[col]]
f <- function(st, i)   fun[study == STUDY_LAB[[st]]][i]
pv <- function(st, col) unique(pca[study == STUDY_LAB[[st]]])[[col]][1]
S1 <- STUDIES[1]; S2 <- STUDIES[2]

legend <- paste0(
"Figure ", FIG, ". The two reused datasets, and the result of re-quantifying them through a ",
"single pipeline. Both studies are public RNA-seq that was downloaded and re-processed here from ",
"raw reads with one nf-core/rnaseq + salmon workflow, each against its own reference genome ",
"(Saccharum hybrid R570 for ", STUDY_LAB[[S1]], "; LA purple for ", STUDY_LAB[[S2]], "); no ",
"count matrix, mapping or quality metric is taken from the original publications. Cross-species ",
"orthology throughout this work comes from OrthoFinder v3.1.3 (-S diamond -M msa -A famsa ",
"-T fasttree) run on one protein per gene of the two reference proteomes, which assigns 360,794 ",
"of 435,856 genes (82.8%) to 94,273 orthogroups; it is the backbone every comparison between the ",
"two species rests on, and it is what makes 'the same gene' a defined object across two polyploid ",
"references.
",
"
",
"(A) Library counts per study, genotype, nitrogen level and leaf segment; the number in each cell ",
"is the number of sequenced libraries, and the segments run base to tip in the order sampled. The two designs are deliberately shown on the same axes ",
"because their mismatch sets what the rest of the work can ask. ", STUDY_LAB[[S1]], " contrasts ",
"two nitrogen levels across four leaf segments -- base0, base, mid and tip, sampled along the ",
"developmental gradient of the leaf -- in two genotypes of contrasting nitrogen-use ",
"efficiency (", sprintf("%d", META[study == STUDY_LAB[[S1]], .N]), " libraries); ",
STUDY_LAB[[S2]], " runs a three-point nitrogen DOSE (0, 2, 6 mM) in one tissue in two species ",
"(", sprintf("%d", META[study == STUDY_LAB[[S2]], .N]), " libraries). The ",
sprintf("%d", META[study == STUDY_LAB[[S1]], .N]), " against ",
sprintf("%d", META[study == STUDY_LAB[[S2]], .N]), " asymmetry is the single most important fact ",
"about this project and nearly every underpowered result traces back to it. That the second ",
"design is an ordered dose rather than a two-level contrast is why the trait is coded ordinally ",
"and tested by rank correlation rather than by Pearson.
",
"
",
"(B) Sequencing depth against mapping rate, one point per library, from the salmon logs of each ",
"nf-core/rnaseq run; fill gives the nitrogen level and shape the study. Depth and mapping rate ",
"are plotted against each other rather than as two separate distributions because the failure ",
"worth seeing is a library that is both shallow and poorly mapped, which is a position on this ",
"plane and not a value on either axis. ", STUDY_LAB[[S1]], ": median ",
sprintf("%.1f", g(S1, "proc_med")), " M fragments (range ", sprintf("%.1f", g(S1, "proc_min")),
"-", sprintf("%.1f", g(S1, "proc_max")), ") and ", sprintf("%.1f%%", g(S1, "map_med")),
" mapped (", sprintf("%.1f", g(S1, "map_min")), "-", sprintf("%.1f", g(S1, "map_max")), "). ",
STUDY_LAB[[S2]], ": median ", sprintf("%.1f", g(S2, "proc_med")), " M fragments (",
sprintf("%.1f", g(S2, "proc_min")), "-", sprintf("%.1f", g(S2, "proc_max")), ") and ",
sprintf("%.1f%%", g(S2, "map_med")), " mapped (", sprintf("%.1f", g(S2, "map_min")), "-",
sprintf("%.1f", g(S2, "map_max")), "). Depth is more variable in ", STUDY_LAB[[S1]],
" and mapping rate lower at its worst (", sprintf("%.1f%%", min(qc[study == STUDY_LAB[[S1]],
percent_mapped])), " in one library), which is expected for a hybrid genotype read against a ",
"single haplotype-resolved reference; no library was excluded on these metrics.
",
"
",
"(C) Genes surviving each stage from annotation to network, per study, with the percentage of the ",
"annotated set retained above each bar. Of ", sprintf("%s", fmt_n(f(S1, 1)$n)), " annotated and ",
"quantified genes in ", STUDY_LAB[[S1]], ", ", sprintf("%s", fmt_n(f(S1, 2)$n)), " (",
sprintf("%.0f%%", 100 * f(S1, 2)$frac), ") clear the coefficient-of-variation filter used to ",
"build the network and ", sprintf("%s", fmt_n(f(S1, 3)$n)), " (",
sprintf("%.0f%%", 100 * f(S1, 3)$frac), ") end up as nodes; for ", STUDY_LAB[[S2]],
" the same three numbers are ", sprintf("%s", fmt_n(f(S2, 1)$n)), ", ",
sprintf("%s", fmt_n(f(S2, 2)$n)), " (", sprintf("%.0f%%", 100 * f(S2, 2)$frac), ") and ",
sprintf("%s", fmt_n(f(S2, 3)$n)), " (", sprintf("%.0f%%", 100 * f(S2, 3)$frac), "). The CV filter ",
"removes a similar share from both, but the step from filtered genes to network nodes does not: a ",
"gene becomes a node only by carrying at least one edge above the significance threshold, and at ",
"n = ", sprintf("%d", META[study == STUDY_LAB[[S2]], .N]), " that threshold is far softer, which ",
"is why ", STUDY_LAB[[S2]], " retains nearly all of its filtered genes as nodes and ",
STUDY_LAB[[S1]], " does not. This panel is the context for every node and edge count reported ",
"later.
",
"
",
"(D) Principal components of the variance-stabilised expression matrix of each study, computed ",
"independently per study on the ", fmt_n(N_PCA), " most variable genes, with the percentage of ",
"total variance carried by each axis in the panel strip. Fill gives nitrogen level, shape gives ",
"genotype. The purpose is to establish, before any network is built, that the experimental design ",
"is visible in the data at all, and to say plainly what the dominant structure is. It is genotype, ",
"in both studies and almost exactly: regressing each component on each design factor, PC1 tracks ",
"genotype at R2 = ", sprintf("%.3f", pv(S1, "r2_pc1_genotype")), " in ", STUDY_LAB[[S1]],
" (", sprintf("%.0f%%", pv(S1, "ve1")), " of total variance) and R2 = ",
sprintf("%.3f", pv(S2, "r2_pc1_genotype")), " in ", STUDY_LAB[[S2]], " (",
sprintf("%.0f%%", pv(S2, "ve2") + pv(S2, "ve1") - pv(S2, "ve2")), " of total variance). ",
"Nitrogen loads on neither first component (R2 = ", sprintf("%.3f", pv(S1, "r2_pc1_treatment")),
" and ", sprintf("%.3f", pv(S2, "r2_pc1_treatment")), "). It appears on PC2, and only clearly in ",
STUDY_LAB[[S2]], ": there PC2 (", sprintf("%.0f%%", pv(S2, "ve2")), ") tracks the nitrogen dose at ",
"R2 = ", sprintf("%.3f", pv(S2, "r2_pc2_treatment")), ", whereas in ", STUDY_LAB[[S1]],
" PC2 (", sprintf("%.0f%%", pv(S1, "ve2")), ") is mostly leaf segment (R2 = ",
sprintf("%.3f", pv(S1, "r2_pc2_segment")), ") with nitrogen a distant second (R2 = ",
sprintf("%.3f", pv(S1, "r2_pc2_treatment")), "). Two things follow, and they shape everything ",
"after this figure. First, the nitrogen response is real but it is nowhere near the dominant ",
"structure in either dataset, so it has to be found AGAINST genotype and tissue rather than read ",
"off a leading component -- which is why the analyses are run within study and the trait is ",
"tested per gene and per module rather than by unsupervised decomposition. Second, because leaf ",
"segment is the second-largest source of variation in ", STUDY_LAB[[S1]], ", every heatmap in ",
"this work orders columns so that replicates of one segment sit adjacent; sorted by library name ",
"they interleave and the segment effect reads as vertical striping across the panel.")

legend_file <- paste0(OUT_PREFIX, "_legend.txt")
writeLines(wrap_at(legend), legend_file)
say("wrote ", basename(legend_file))

stats <- rbindlist(list(
  qc_sum[, .(panel = "B", study = as.character(study),
             quantity = "libraries / median M fragments / median % mapped",
             value = sprintf("%d / %.1f / %.1f", libs, proc_med, map_med))],
  fun[, .(panel = "C", study = as.character(study),
          quantity = gsub("\n", " ", stage),
          value = sprintf("%s (%.1f%%)", fmt_n(n), 100 * frac))],
  unique(pca[, .(panel = "D", study = as.character(study),
                 quantity = "PC1 / PC2 % variance",
                 value = sprintf("%.1f / %.1f", ve1, ve2))])))
write_tsv(stats, paste0(OUT_PREFIX, "_stats.tsv"))
say("done")
