# =============================================================================
# 24_fig_modules.r — the module-level figure: the nitrogen response at the
# resolution where it can actually be tested.
#
# The gene level dies on a testing burden of 39,226 / 44,118 genes. One eigengene
# per MCL module cuts that to ~6,500 tests per study, and a module signal can
# survive where a single gene cannot. This figure shows what that buys, how the
# responsive set was selected, that it is not noise, and what the responding
# modules look like.
#
#   A  SELECTION. Every module as |rho| against its BH-adjusted p, with both
#      thresholds drawn. The panel exists to show WHICH THRESHOLD BINDS, and it
#      is a different one in each study: at n = 48 an |rho| of 0.6 already implies
#      p ~ 6e-06 so BH never binds for sugarcane and the effect-size floor is the
#      only active constraint; at n = 18 it is the reverse.
#
#   B  WHY SPEARMAN. The same eigengenes, the SAME blocked model, scored by
#      Pearson and by Spearman at identical thresholds -- one variable changes. Purple's responsive set more than doubles, which is
#      exactly where the ordinal argument predicts it -- purple is the study with
#      three nitrogen levels, and Pearson reads 0/2/6 mM as arithmetic the design
#      never claimed.
#
#   C, D  EVERY RESPONSIVE MODULE'S EIGENGENE, one panel per species. Rows are
#      modules, columns are libraries, colour is the z-scored eigengene. Rows are
#      split by the sign of rho and columns by nitrogen status, so both splits
#      label themselves and no separate annotation strip is needed. A mean
#      profile would show the same two shapes in four points; this shows whether
#      the set is coherent -- whether every module in a block really moves the
#      same way, or whether the block average is carried by a minority.
#
#      Drawn with geom_raster rather than ComplexHeatmap: these have to compose
#      with two ggplots in one patchwork, and a heatmap object would have to be
#      grabbed into a grob to get there. Rows are ordered by rho within each
#      block, which is deterministic and interpretable, rather than clustered.
#
# NOTHING IS WRITTEN ON THE FIGURE THAT BELONGS IN THE LEGEND -- panel letters and
# the labels the data needs, nothing else. The permutation null is a pair of
# numbers and lives in the legend rather than costing a panel.
#
# RUN: through run.sh  ->  ./run.sh figmodules   (then ./run.sh legends)
# =============================================================================

suppressMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(scico); library(grid); library(svglite); library(scales)
})

source(file.path(dirname(sub("--file=", "",
       grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "lib", "common.R"))

STUDIES    <- env_list("CLEAN_STUDIES", c("sugarcane", "purple"))
OUT_PREFIX <- env_req("CLEAN_OUT_PREFIX")
FIG        <- env_opt("CLEAN_FIG_NUM", "5")
R_THR      <- env_num("CLEAN_MODULE_R_THR", 0.6)
PADJ_THR   <- env_num("CLEAN_MODULE_PADJ_THR", 0.05)
SELECT_TRAIT <- env_opt("CLEAN_SELECT_TRAIT", "treatment")
setDTthreads(as.integer(env_num("CLEAN_CORES", 8)))

cfg <- lapply(STUDIES, function(st) list(
  profile = env_req(sprintf("CLEAN_PROFILE_%s", toupper(st))),
  null    = env_req(sprintf("CLEAN_NULL_%s", toupper(st))),
  eigen   = env_req(sprintf("CLEAN_EIGEN_%s", toupper(st))),
  meta    = env_req(sprintf("CLEAN_META_%s", toupper(st))),
  traits  = env_req(sprintf("CLEAN_TRAITS_%s", toupper(st)))))
names(cfg) <- STUDIES

banner(sprintf("Figure %s — module-level nitrogen response", FIG))

DIRS <- c("positive", "negative")
PAL_STUDY <- setNames(scico(3, palette = "batlow")[1:2], STUDIES)
PAL_DIR   <- setNames(rev(scico(3, palette = "managua", begin = 0.08, end = 0.92))[c(1, 3)],
                      DIRS)
# The shared nitrogen-status scale of the dataset figure, so one colour means the
# same thing across the paper.
NSTATUS_LEVELS <- c("low N", "control", "high N")
NSTATUS <- c("Low Nitrogen" = "low N", "High Nitrogen" = "high N",
             "0N" = "low N", "2N" = "control", "6N" = "high N")
PAL_NSTATUS <- setNames(rev(scico(3, palette = "managua", begin = 0.08, end = 0.92)),
                        NSTATUS_LEVELS)

theme_f <- theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(face = "bold", size = 7.5),
        legend.key.size = unit(3.5, "mm"),
        legend.title = element_text(face = "bold", size = 7.5),
        plot.margin = margin(2, 4, 2, 2))

PROF <- rbindlist(lapply(STUDIES, function(st)
  fread(cfg[[st]]$profile)[, study := st]), fill = TRUE)
PROF[, study := factor(study, levels = STUDIES)]
NUL <- rbindlist(lapply(STUDIES, function(st)
  fread(cfg[[st]]$null)[, study := st]))

resp <- PROF[, .(tested = .N, responsive = sum(responsive),
                 pos = sum(direction == "positive"), neg = sum(direction == "negative"),
                 med_rho = median(abs(rho[responsive])),
                 n_padj = sum(padj <= PADJ_THR), n_rho = sum(abs(rho) >= R_THR)),
             by = study]
say("module response:"); print(resp, row.names = FALSE)

# =============================================================================
# A — selection: |rho| against padj, with both thresholds
# =============================================================================
# -log10(padj) rather than padj on a log axis: significance should increase
# upwards, which is what a reader expects of anything volcano-shaped.
pA <- ggplot(PROF, aes(abs(rho), -log10(pmax(padj, 1e-16)))) +
  geom_point(aes(colour = responsive), size = 0.22, alpha = 0.35, shape = 16) +
  geom_vline(xintercept = R_THR, linetype = 2, linewidth = 0.3, colour = "grey30") +
  geom_hline(yintercept = -log10(PADJ_THR), linetype = 2, linewidth = 0.3,
             colour = "grey30") +
  facet_wrap(~ study, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = c(`FALSE` = "grey78", `TRUE` = unname(PAL_STUDY[1])),
                      labels = c("not responsive", "responsive"), name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 1.6, alpha = 1))) +
  labs(x = expression("|Spearman "*rho*"|  with nitrogen"),
       y = expression(-log[10]*"(BH-adjusted p)")) +
  theme_f

# =============================================================================
# B — Spearman against Pearson on the same eigengenes
# =============================================================================
parse_traits <- function(spec) {
  out <- list()
  for (block in strsplit(spec, ";", fixed = TRUE)[[1L]]) {
    block <- trimws(block); if (!nzchar(block)) next
    nm <- sub(":.*$", "", block); body <- sub("^[^:]*:", "", block)
    kv <- strsplit(strsplit(body, ",", fixed = TRUE)[[1L]], "=", fixed = TRUE)
    v <- as.numeric(vapply(kv, `[`, "", 2L)); names(v) <- trimws(vapply(kv, `[`, "", 1L))
    out[[trimws(nm)]] <- v
  }
  out
}

eig_and_trait <- function(st) {
  E <- read_vst(cfg[[st]]$eigen)
  m <- fread(cfg[[st]]$meta); setnames(m, tolower(names(m)))
  m <- m[match(colnames(E), sample)]
  enc <- parse_traits(cfg[[st]]$traits)[[SELECT_TRAIT]]
  tv <- unname(enc[as.character(m[[SELECT_TRAIT]])])
  ok <- !is.na(tv)
  list(E = E[, ok, drop = FALSE], t = tv[ok],
       status = factor(NSTATUS[as.character(m[[SELECT_TRAIT]])[ok]],
                       levels = NSTATUS_LEVELS))
}
ET <- lapply(STUDIES, eig_and_trait); names(ET) <- STUDIES

cmpB <- rbindlist(lapply(STUDIES, function(st) {
  E <- ET[[st]]$E; tv <- ET[[st]]$t
  n <- length(tv); df <- n - 2L
  pe <- as.numeric(cor(t(E), tv))
  pp <- p.adjust(2 * pt(-abs(pe * sqrt(df / (1 - pe^2 + 1e-15))), df), "BH")
  sp <- PROF[study == st]; setkey(sp, module); sp <- sp[rownames(E)]
  # HOLD THE MODEL FIXED, VARY ONLY THE STATISTIC. Since 53 the responsive call
  # is a BLOCKED partial correlation, so pairing it with the marginal Pearson
  # computed above would change the statistic AND the model in one bar chart --
  # and this panel's entire claim is "why Spearman". 53 writes both blocked
  # columns, so when they are present the comparison is blocked-vs-blocked.
  if (all(c("rho_blocked_pearson", "padj_blocked_pearson") %in% names(sp))) {
    n_pear <- sum(sp$padj_blocked_pearson <= PADJ_THR &
                  abs(sp$rho_blocked_pearson) >= R_THR, na.rm = TRUE)
  } else n_pear <- sum(pp <= PADJ_THR & abs(pe) >= R_THR)
  data.table(study = st, Pearson = n_pear, Spearman = sum(sp$responsive))
}))
cmpB <- melt(cmpB, id.vars = "study", variable.name = "statistic", value.name = "n")
cmpB[, study := factor(study, levels = STUDIES)]
say("panel B: responsive modules by statistic"); print(cmpB, row.names = FALSE)

pB <- ggplot(cmpB, aes(study, n, fill = statistic)) +
  geom_col(position = position_dodge(width = 0.65), width = 0.55) +
  geom_text(aes(label = n), position = position_dodge(width = 0.65),
            vjust = -0.4, size = 2.3) +
  scale_fill_manual(values = c(Pearson = "grey70",
                               Spearman = unname(PAL_STUDY[1])), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = NULL, y = "nitrogen-responsive modules") +
  theme_f

# =============================================================================
# C, D — every responsive module's eigengene, per species
# =============================================================================
PAL_Z <- rev(scico(256, palette = "roma"))

heat_data <- function(st) {
  E <- ET[[st]]$E; sta <- ET[[st]]$status
  m <- fread(cfg[[st]]$meta); setnames(m, tolower(names(m)))
  m <- m[match(colnames(E), sample)]
  sel <- PROF[study == st & responsive == TRUE]
  # rho descending inside "rises", ascending inside "falls", so the strongest
  # module of each block sits at the top of it.
  sel[, direction := factor(direction, levels = DIRS)]
  sel <- sel[order(direction, -abs(rho))]
  M <- E[sel$module, , drop = FALSE]

  # Columns: nitrogen status blocks, then genotype, then leaf segment (sugarcane
  # only), then library, so replicates of one condition sit together.
  segcol <- if ("segment" %in% names(m)) m$segment else rep("", nrow(m))
  ord <- order(sta, m$genotype, segcol, m$sample)
  M <- M[, ord, drop = FALSE]

  d <- data.table(module = factor(rep(rownames(M), times = ncol(M)),
                                  levels = rev(sel$module)),
                  sample = factor(rep(colnames(M), each = nrow(M)),
                                  levels = colnames(M)),
                  z = as.vector(M))
  d[, direction := factor(c(positive = "rises with N", negative = "falls with N")[
      as.character(sel$direction[match(as.character(module), sel$module)])],
      levels = c("rises with N", "falls with N"))]
  d[, status := sta[ord][match(as.character(sample), colnames(M))]]
  d[, study := st]
  d[]
}
HD <- lapply(STUDIES, heat_data); names(HD) <- STUDIES

# ONE colour scale across both species, so a cell in C means what a cell in D
# means. Clipped at the 98th percentile of |z| over both responsive sets
# together; the colourbar is drawn once, on C.
ZL <- as.numeric(quantile(abs(rbindlist(HD)$z), 0.98, na.rm = TRUE))
if (ZL == 0) ZL <- 1

heat_panel <- function(st, show_legend) {
  d <- HD[[st]]
  n_mod <- uniqueN(d$module)
  ggplot(d, aes(sample, module, fill = z)) +
    geom_raster() +
    facet_grid(direction ~ status, scales = "free", space = "free", switch = "y") +
    scale_fill_gradientn(colours = PAL_Z, limits = c(-ZL, ZL), oob = squish,
                         name = "eigengene z",
                         guide = if (show_legend)
                           guide_colourbar(barheight = unit(2.2, "cm"),
                                           barwidth = unit(2.6, "mm")) else "none") +
    labs(x = sprintf("%s libraries", ncol(ET[[st]]$E)),
         y = sprintf("%s  —  %s responsive modules", st, fmt_n(n_mod))) +
    theme_f +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(), panel.spacing = unit(0.8, "mm"),
          strip.placement = "outside",
          strip.text.y.left = element_text(angle = 90, size = 6.5),
          strip.text.x = element_text(size = 6.5),
          axis.title = element_text(size = 7))
}
pC <- heat_panel(STUDIES[1], TRUE)
pD <- heat_panel(STUDIES[2], FALSE)

# The coherence the panels are for, as a number: how much of each module's
# eigengene variance the nitrogen split accounts for, per direction.
coh <- rbindlist(lapply(STUDIES, function(st) {
  E <- ET[[st]]$E; sta <- ET[[st]]$status
  sel <- PROF[study == st & responsive == TRUE]
  r2 <- apply(E[sel$module, , drop = FALSE], 1,
              function(v) summary(lm(v ~ sta))$r.squared)
  data.table(study = st, direction = sel$direction, r2 = r2)[
    , .(modules = .N, median_r2 = median(r2)), by = .(study, direction)]
}))
say("panels C/D: variance of each eigengene explained by nitrogen status")
print(coh, row.names = FALSE)

# TF enrichment keeps its numbers for the legend even though it lost its panel.
tfD <- rbindlist(lapply(STUDIES, function(st) {
  d <- PROF[study == st]
  grp <- list(`not responsive` = d[responsive == FALSE],
              responsive       = d[responsive == TRUE],
              `rises with N`   = d[direction == "positive"],
              `falls with N`   = d[direction == "negative"])
  base <- grp[["not responsive"]]
  rbindlist(lapply(names(grp), function(g) {
    r <- grp[[g]]
    if (!nrow(r)) return(NULL)
    tab <- matrix(c(sum(r$tf_enriched), nrow(r) - sum(r$tf_enriched),
                    sum(base$tf_enriched), nrow(base) - sum(base$tf_enriched)), nrow = 2)
    ft <- if (g == "not responsive") NULL else fisher.test(tab)
    data.table(study = st, group = g, modules = nrow(r),
               pct = 100 * mean(r$tf_enriched),
               OR = if (is.null(ft)) NA_real_ else unname(ft$estimate),
               p  = if (is.null(ft)) NA_real_ else ft$p.value)
  }))
}))
say("TF enrichment (legend only, no longer a panel)")
print(tfD, row.names = FALSE)

# =============================================================================
# compose
# =============================================================================
fig <- (pA | pB) / (pC | pD) +
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

rs <- function(st, col) resp[study == st][[col]]
cb <- function(st, stat) cmpB[study == st & statistic == stat, n]
tf <- function(st, g, col) tfD[study == st & group == g][[col]]
nn <- function(st, col) NUL[study == st][[col]]
S1 <- STUDIES[1]; S2 <- STUDIES[2]
n_perm <- NUL[study == S1, .N]
obs_ge <- function(st) NUL[study == st & n_responsive >= rs(st, "responsive"), .N]

legend <- paste0(
"Figure ", FIG, ". The nitrogen response at module level. Each MCL module of >= 3 genes is ",
"summarised by an eigengene -- PC1 of its members' variance-stabilised expression, genes scaled ",
"before the decomposition per the WGCNA convention, oriented to mean module expression and ",
"z-scored -- giving ", fmt_n(rs(S1, "tested")), " and ", fmt_n(rs(S2, "tested")),
" tests per study instead of the 39,226 and 44,118 the gene level carries. Association with ",
"nitrogen is SPEARMAN'S RHO ALONE, at padj <= ", sprintf("%.2f", PADJ_THR),
" (Benjamini-Hochberg over every module in the study) and |rho| >= ", sprintf("%.1f", R_THR),
", the same two thresholds the gene level uses. Neither Pearson nor mutual information is used ",
"here: the trait is ordinal, and MI at this resolution proved to be an omnibus test firing on ",
"quirks in one or two libraries.\n",
"\n",
"(A) Every module as |rho| against its adjusted p, with both thresholds dashed; coloured points ",
"clear both. The panel is here to show WHICH THRESHOLD BINDS, and it is a different one in each ",
"study. In ", S1, " (n = 48) an |rho| of ", sprintf("%.1f", R_THR), " already implies p ~ 6e-06, ",
"so BH never binds: ", fmt_n(rs(S1, "n_padj")), " modules clear padj and ",
fmt_n(rs(S1, "n_rho")), " clear the effect-size floor, and the floor alone defines the ",
fmt_n(rs(S1, "responsive")), " responsive modules -- the corrected and uncorrected readings are ",
"the same set. In ", S2, " (n = 18) it is the reverse: ", fmt_n(rs(S2, "n_rho")),
" modules reach |rho| >= ", sprintf("%.1f", R_THR), " and only ", fmt_n(rs(S2, "n_padj")),
" survive BH, so ", S2, "'s uncorrected count would be ", fmt_n(rs(S2, "n_rho")),
" and should be labelled as such wherever it is used.\n",
"\n",
"(B) The same eigengenes under the same blocked model, scored by Pearson and by Spearman at ",
"identical thresholds -- only the statistic changes. ", S1, ": ",
fmt_n(cb(S1, "Pearson")), " modules against ", fmt_n(cb(S1, "Spearman")), ". ", S2, ": ",
fmt_n(cb(S2, "Pearson")), " against ", fmt_n(cb(S2, "Spearman")), " -- MORE THAN DOUBLE, and ",
S2, " is exactly where the ordinal argument predicts the gain, because it is the study with ",
"three nitrogen levels and Pearson reads 0/2/6 mM as arithmetic the design never claimed. ", S1,
"'s trait is two-level, where Spearman is the rank-biserial correlation, so the gain there is ",
"smaller and is mostly robustness to the outliers a PC1 can carry.\n",
"\n",
"(C, D) Every responsive module\'s eigengene, one panel per species: rows are modules, columns ",
"are libraries, colour is the z-scored eigengene on a single scale shared by both panels and ",
"clipped at the 98th percentile of |z| (the colourbar is drawn once, on C). Rows are split by the ",
"sign of rho and ordered by |rho| within each block; columns are split by nitrogen status and ",
"ordered by genotype, then leaf segment, then library, so replicates of one condition sit ",
"together. Both splits label themselves, so no annotation strip is needed. A mean profile would ",
"show the same two shapes in four points; what these panels add is whether the SET is coherent -- ",
"whether every module in a block really moves the same way or whether the block is carried by a ",
"minority. It is coherent: the nitrogen split accounts for a median ",
sprintf("%.2f", coh[study == S1 & direction == "positive", median_r2]), " and ",
sprintf("%.2f", coh[study == S1 & direction == "negative", median_r2]),
" of each eigengene\'s variance in ", S1, "\'s two blocks, and ",
sprintf("%.2f", coh[study == S2 & direction == "positive", median_r2]), " and ",
sprintf("%.2f", coh[study == S2 & direction == "negative", median_r2]), " in ", S2, "\'s.\n",
"\n",
"The two panels are not on the same footing and should not be read as if they were. ", S1,
"\'s ", fmt_n(rs(S1, "responsive")), " modules come from 48 libraries in a two-level contrast, so ",
"the blocks are simply warm-on-one-side and cool-on-the-other. ", S2, "\'s ",
fmt_n(rs(S2, "responsive")), " come from 18 libraries across three levels that are NOT a dose ",
"series -- 2 N is the CONTROL and 0 N and 6 N are stresses in opposite directions -- so its ",
"control column sits near zero between two opposite extremes, which is the design showing through ",
"rather than a weak response. A monotone rank test on that design asks whether a module tracks ",
"nitrogen SUPPLY, not whether it responds to nitrogen STRESS: a module moved the same way by both ",
"stresses is invisible to it. The cost was measured rather than assumed -- a U-shape contrast ",
"c(+1,-2,+1) over the three levels finds ONE significant ", S2, " module, which Spearman misses, ",
"against Spearman\'s ", fmt_n(rs(S2, "responsive")), ". At n = 18 the U-shape test has almost no ",
"power, the same wall every other ", S2, " result hits, but the ", fmt_n(rs(S2, "responsive")),
" should be described as tracking nitrogen supply rather than as stress responders.\n",
"\n",
"NOT SHOWN AS A PANEL, but the biological payoff of the module level: transcription-factor ",
"enrichment. A hypergeometric test per module against that network\'s own node universe (BH across ",
"modules), then Fisher against the non-responsive modules, puts ", S1, "\'s responsive modules at ",
sprintf("%.2f%%", tf(S1, "responsive", "pct")), " TF-enriched against ",
sprintf("%.2f%%", tf(S1, "not responsive", "pct")), " (odds ratio ",
sprintf("%.2f", tf(S1, "responsive", "OR")), ", p = ", sprintf("%.4f", tf(S1, "responsive", "p")),
"), and the signal sits in the modules that FALL with nitrogen (OR ",
sprintf("%.2f", tf(S1, "falls with N", "OR")), ", p = ",
sprintf("%.4f", tf(S1, "falls with N", "p")), ") rather than those that rise (OR ",
sprintf("%.2f", tf(S1, "rises with N", "OR")), ", p = ",
sprintf("%.3f", tf(S1, "rises with N", "p")), "). Do not over-read the direction split -- it is 9 ",
"enriched modules against 6 -- but the pooled result is solid, and it is one the previous ",
"three-statistic call had buried: split across pearson_only / mi_only / both it read as OR 1.96, ",
"p = 0.06, with the strongest class at 3 of 41 modules. In ", S2, " nothing is enriched (0 of ",
fmt_n(tf(S2, "responsive", "modules")), " responsive modules), and its TF-enrichment rate is an ",
"order of magnitude below ", S1, "\'s everywhere in the network.\n",
"\n",
"NOT SHOWN, because it is two numbers rather than a panel: the responsive sets are not noise. ",
"Permuting the trait labels against the same eigengenes ", fmt_n(n_perm), " times gives a median ",
"of 0 responsive modules in both studies; ", fmt_n(obs_ge(S1)), " of ", fmt_n(n_perm),
" permutations reach ", S1, "'s ", fmt_n(rs(S1, "responsive")), ", and ", fmt_n(obs_ge(S2)),
" of ", fmt_n(n_perm), " reach ", S2, "'s ", fmt_n(rs(S2, "responsive")), " (empirical p = ",
sprintf("%.3f", (obs_ge(S2) + 1) / (n_perm + 1)), "). ", S2, "'s null has a long right tail, ",
"however -- its worst shuffle produced ", fmt_n(max(nn(S2, "n_responsive"))),
" responsive modules, far above the observed count -- because its network is a single dense ",
"component whose eigengenes are strongly correlated, so one lucky label assignment lights up many ",
"at once and BH, being adaptive, then loosens for all of them. ", S2, "'s ",
fmt_n(rs(S2, "responsive")), " modules are a real signal and are NOT ", fmt_n(rs(S2, "responsive")),
" independent findings.")

legend_file <- paste0(OUT_PREFIX, "_legend.txt")
writeLines(wrap_at(legend), legend_file)
say("wrote ", basename(legend_file))

stats <- rbindlist(list(
  resp[, .(panel = "A", study = as.character(study),
           quantity = "tested / padj<=thr / |rho|>=thr / responsive",
           value = sprintf("%s / %s / %s / %s", fmt_n(tested), fmt_n(n_padj),
                           fmt_n(n_rho), fmt_n(responsive)))],
  cmpB[, .(panel = "B", study = as.character(study),
           quantity = sprintf("responsive by %s", statistic), value = fmt_n(n))],
  coh[, .(panel = "C/D", study = as.character(study),
          quantity = sprintf("median R2 on nitrogen, %s modules", direction),
          value = sprintf("%.3f (n = %d)", median_r2, modules))],
  tfD[!is.na(OR), .(panel = "legend", study = as.character(study),
                    quantity = sprintf("TF enrichment: %s", group),
                    value = sprintf("%.2f%% OR %.2f p %.4g", pct, OR, p))]))
write_tsv(stats, paste0(OUT_PREFIX, "_stats.tsv"))
say("done")
