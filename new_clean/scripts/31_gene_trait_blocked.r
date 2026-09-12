#!/usr/bin/env Rscript
# =============================================================================
# 31_gene_trait_blocked.r — the gene-trait test with the design in the model
#
# WHY THIS EXISTS. 07_gene_trait_cor.r:135 is a marginal `cor(expr, trait)`. It
# has no genotype term and no segment term. This project's own QC log records
#
#     purple    genotype R2 on PC1 = 0.999      sugarcane genotype R2 = 0.998
#     sugarcane segment  R2 on PC2 = 0.802
#
# so the single largest source of variance in either matrix is sitting in the
# RESIDUAL of every gene-level nitrogen test. That inflates the residual variance,
# deflates every t, and it is why purple needs |r| ~ 0.80 to clear BH at n = 18.
#
# The symptom is visible without fitting anything: purple's marginal p-value
# histogram RISES in the last decile (11.8% against 10% expected) instead of being
# flat. A well-behaved mixture of null and signal is flat or decreasing. A rising
# right tail means the null is not uniform, which is what an unmodelled variance
# component of that size does.
#
# WHAT THIS IS NOT. It is not a looser threshold. The FDR level and the |r| floor
# are the SAME values 07 uses; only the model changes.
#
# THE MODELS
#   sugarcane   expr ~ genotype + segment + N      n=48, resid df 42
#   purple      expr ~ genotype + N                n=18, resid df 15
#
# and for sugarcane one control, because its 48 libraries are NOT 48 replicates:
#
#   PSEUDOREPLICATION. The 48 sugarcane libraries are 12 plants x 4 leaf segments
#   (base0/base/mid/tip), i.e. repeated measures on the same plant -- the trailing
#   _1/_2/_3 of the sample name is the plant within a genotype x N cell. The
#   marginal test's df = 46 therefore overstates the independent replication for
#   the nitrogen contrast roughly 4-fold. Putting segment in as a fixed block
#   removes the segment MEANS but does not account for within-plant correlation,
#   so a PLANT-LEVEL run (segments averaged, n = 12, resid df 9) is computed
#   alongside and written to the same table. If the two disagree, the blocked
#   p-values are not to be trusted, and saying so is the point of carrying it.
#
# WHY SPEARMAN FOR PURPLE. Its trait is ORDINAL -- a 0/2/6 mM dose, and 2 mM is
# the control of a stress-control-stress design. Pearson reads those spacings
# literally, i.e. it asks whether the 2->6 response is exactly twice the 0->2 one.
# Nothing in the design justifies that. This is the same argument that already
# moved the MODULE level to Spearman (19_module_trait_spearman.r:8-27, carried into
# 53_module_trait_blocked.r). Both statistics are computed and both are written, so
# the choice is visible in the table rather than buried in the code.
#
# ---------------------------------------------------------------------------
# THE GENE UNIVERSE IS THE PEARSON-ONLY NETWORK'S NODE SET, and that is a change
# from the first version of this script. It used to filter on
# conserved_genes_<study>_FULL.txt -- the genes carrying a conserved edge in the
# MERGED (Pearson + MI) graph, 39,226 / 44,118 of them. The analysis no longer uses
# that graph. The universe is now network_<study>_node_metrics.tsv, i.e. every node
# of the unpruned Pearson-only network (101,990 / 170,135), which is the same gene
# universe the modules are built on. It raises the BH denominator 2.6x / 3.9x, so
# the counts this writes are NOT comparable with the ones the old table held.
# ---------------------------------------------------------------------------
#
# THE MARGINAL FIT IS COMPUTED HERE, on these genes, NOT read from 07's table. 07
# corrects over a different gene set, so joining the two would compare two
# analyses and call the difference an effect of blocking. 53_module_trait_blocked.r
# learned this the hard way at module level (docs/results.md, "One error caught").
# Same solver, no blocks, identical rows.
#
# RUN: through run.sh  ->  ./run.sh traitblocked purple
# =============================================================================
suppressPackageStartupMessages({library(data.table)})
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE),
                                                 value = TRUE)[1])), "lib", "common.R"))

STUDY       <- env_req("CLEAN_STUDY")
VST_PREFIX  <- env_req("CLEAN_VST_PREFIX")
META_FILE   <- env_req("CLEAN_META")
TRAIT_SPEC  <- env_req("CLEAN_TRAITS")
OUT_FILE    <- env_req("CLEAN_OUT_FILE")
GENE_FILTER <- env_opt("CLEAN_GENE_FILTER")
TRAIT       <- env_opt("CLEAN_SELECT_TRAIT", "treatment")
BLOCK       <- strsplit(trimws(env_opt("CLEAN_BLOCK", "genotype")), "[ ,]+")[[1L]]
STAT        <- env_opt("CLEAN_BLOCKED_STAT", "pearson")
PLANT_FROM  <- env_opt("CLEAN_PLANT_FROM", "")   # non-empty -> plant-level control
R_THR       <- env_num("CLEAN_TRAIT_R_THR", 0.6)
PADJ_THR    <- env_num("CLEAN_TRAIT_PADJ_THR", 0.05)
setDTthreads(as.integer(env_num("CLEAN_CORES", 100)))

if (!STAT %in% c("pearson", "spearman"))
  stop("CLEAN_BLOCKED_STAT must be pearson or spearman (got '", STAT, "')", call. = FALSE)

banner(sprintf("blocked gene-trait test: %s  [expr ~ %s + %s, %s]",
               STUDY, paste(BLOCK, collapse = " + "), TRAIT, STAT))

# --- trait encoding, the same parser 07, 19, 31 and 53 use -------------------
parse_traits <- function(spec) {
  out <- list()
  for (blk in strsplit(spec, ";", fixed = TRUE)[[1L]]) {
    blk <- trimws(blk)
    if (!nzchar(blk)) next
    nm   <- sub(":.*$", "", blk)
    body <- sub("^[^:]*:", "", blk)
    kv   <- strsplit(strsplit(body, ",", fixed = TRUE)[[1L]], "=", fixed = TRUE)
    if (any(lengths(kv) != 2L))
      stop("cannot parse trait spec '", blk, "'", call. = FALSE)
    v <- as.numeric(vapply(kv, `[`, "", 2L))
    names(v) <- trimws(vapply(kv, `[`, "", 1L))
    out[[trimws(nm)]] <- v
  }
  out
}
TRAITS <- parse_traits(TRAIT_SPEC)
if (!TRAIT %in% names(TRAITS))
  stop("CLEAN_SELECT_TRAIT '", TRAIT, "' is not encoded", call. = FALSE)
enc <- TRAITS[[TRAIT]]
say("trait: ", TRAIT, "  ", paste(sprintf("%s=%g", names(enc), enc), collapse = ", "))

# THE SOLVER IS lib/common.R's. `fit_blocked`, `verify_against_lm` and
# `row_midranks` used to be private copies in this file; they now live in
# common.R:198 onward, shared with 53_module_trait_blocked.r. One implementation
# and one verification, so the gene level and the module level cannot drift apart.

# --- inputs -------------------------------------------------------------------
say("reading VST from ", basename(VST_PREFIX), ".f32")
vst <- read_vst(VST_PREFIX)
say("  ", fmt_n(nrow(vst)), " genes x ", ncol(vst), " samples")

if (nzchar(GENE_FILTER)) {
  vst <- restrict_to_universe(vst, GENE_FILTER)
} else {
  say("NO gene filter -- testing every VST gene. The network node set is the ",
      "universe the modules use; correcting over a different one makes the counts ",
      "incomparable with them.")
}
say("genes to test: ", fmt_n(nrow(vst)))

meta <- fread(META_FILE, header = TRUE)
setnames(meta, tolower(names(meta)))
m <- meta[match(colnames(vst), sample)]
if (anyNA(m$sample)) stop("VST samples missing from ", basename(META_FILE), call. = FALSE)
for (cn in c(TRAIT, BLOCK))
  if (!cn %in% names(m)) stop("column '", cn, "' not in ", basename(META_FILE),
                              "\n  available: ", paste(names(m), collapse = ", "),
                              call. = FALSE)

y <- unname(enc[as.character(m[[TRAIT]])])
if (anyNA(y))
  stop("samples with an unencoded ", TRAIT, ": ",
       paste(unique(as.character(m[[TRAIT]])[is.na(y)]), collapse = ", "), call. = FALSE)
blocks <- setNames(lapply(BLOCK, function(b) factor(m[[b]])), BLOCK)

say("design: ", ncol(vst), " samples")
for (b in BLOCK) say("  ", b, ": ", paste(sprintf("%s=%d", names(table(blocks[[b]])),
                                                  table(blocks[[b]])), collapse = "  "))
say("  ", TRAIT, ": ", paste(sprintf("%g=%d", as.numeric(names(table(y))), table(y)),
                             collapse = "  "))

row_var <- rowSums((vst - rowMeans(vst))^2)
ok <- which(is.finite(row_var) & row_var > 0)
if (length(ok) < nrow(vst))
  say("dropping ", fmt_n(nrow(vst) - length(ok)), " zero-variance genes")
V <- vst[ok, , drop = FALSE]
genes <- rownames(V)

# --- both statistics, always ---------------------------------------------------
# Spearman here is the blocked model fitted to MIDRANKS on both sides, which is
# what "Spearman with a block" means: rank, then residualise. Ties are handled by
# rank()'s default "average", as 19_module_trait_spearman.r:138 requires.
say("")
say("fitting expr ~ ", paste(BLOCK, collapse = " + "), " + ", TRAIT)
fit_p <- fit_blocked(V, y, blocks)
say("  pearson : residual df = ", fit_p$df)

Rk <- row_midranks(V)                            # genes x samples, midranks
fit_s <- fit_blocked(Rk, rank(y), blocks)
say("  spearman: residual df = ", fit_s$df)

prim <- if (STAT == "spearman") fit_s else fit_p

# THE MARGINAL FIT, ON THESE GENES. Same solver, no blocks, identical rows -- so
# blocked-vs-marginal is two MODELS of one gene set. It is NOT read from 07's
# table: that table corrects over a different gene universe, so a join would
# attribute a change of denominator to the change of model.
fit_m <- if (STAT == "spearman") {
  fit_blocked(Rk, rank(y), list())
} else {
  fit_blocked(V, y, list())
}
say("  marginal (no blocks): residual df = ", fit_m$df)

out <- data.table(
  gene = genes, trait = TRAIT,
  model = sprintf("expr ~ %s + %s [%s]", paste(BLOCK, collapse = " + "), TRAIT, STAT),
  r_partial = prim$r, t = prim$t, pval = prim$p,
  n = ncol(V), df = prim$df,
  r_pearson = fit_p$r, p_pearson = fit_p$p,
  r_spearman = fit_s$r, p_spearman = fit_s$p,
  r_marginal = fit_m$r, p_marginal = fit_m$p, df_marginal = fit_m$df)
out[, padj := p.adjust(pval, method = "BH")]
out[, padj_marginal := p.adjust(p_marginal, method = "BH")]
out[, responsive          := padj <= PADJ_THR & abs(r_partial) >= R_THR]
out[, responsive_marginal := padj_marginal <= PADJ_THR & abs(r_marginal) >= R_THR]
out[!is.finite(r_partial),  responsive := FALSE]
out[!is.finite(r_marginal), responsive_marginal := FALSE]

# --- sugarcane's plant-level control ------------------------------------------
# Averaging the 4 segments of a plant collapses the repeated measures into 12
# genuinely independent units. It costs df (9 instead of 42) and it cannot see a
# segment-specific response, but it cannot be accused of pseudoreplication, which
# is exactly what it is here to answer.
out[, `:=`(r_plant = NA_real_, p_plant = NA_real_, padj_plant = NA_real_)]
PLANT_DF <- NA_integer_; PLANT_R <- NA_real_
if (nzchar(PLANT_FROM)) {
  # plant identity = genotype x trait level x the replicate suffix of the sample
  # name (e.g. 10_R_270N_B_3 and 9_R_270N_B0_3 are segments of the same plant).
  rep_id <- sub("^.*_", "", m$sample)
  if (!all(grepl("^[0-9]+$", rep_id)))
    stop("cannot read a replicate number off the sample names; ",
         "CLEAN_PLANT_FROM assumes a trailing _<n>", call. = FALSE)
  plant <- paste(as.character(m[[PLANT_FROM]]), as.character(m[[TRAIT]]), rep_id, sep = "|")
  idx <- split(seq_len(ncol(V)), plant)
  say("")
  say("plant-level control: ", length(idx), " plants from ", ncol(V), " libraries",
      "  (", paste(sort(unique(lengths(idx))), collapse = "/"), " libraries each)")
  if (length(idx) >= ncol(V))
    stop("plant grouping did not collapse anything -- check CLEAN_PLANT_FROM",
         call. = FALSE)
  P  <- vapply(idx, function(i) rowMeans(V[, i, drop = FALSE]), numeric(nrow(V)))
  yp <- vapply(idx, function(i) y[i[1]], numeric(1))
  # A block that VARIES within a plant (segment) is what the averaging removed,
  # so it must not be carried over -- keeping it would fit a meaningless factor
  # and silently eat degrees of freedom. Only blocks constant within a plant
  # (genotype) survive, and they must still vary between plants to be estimable.
  const_within <- vapply(BLOCK, function(b)
    all(vapply(idx, function(i) length(unique(as.character(m[[b]])[i])) == 1L, logical(1))),
    logical(1))
  bp <- setNames(lapply(BLOCK, function(b)
          factor(vapply(idx, function(i) as.character(m[[b]])[i[1]], character(1)))), BLOCK)
  keepb <- const_within & vapply(bp, nlevels, integer(1)) > 1L
  if (any(!keepb))
    say("  blocks removed by the averaging: ", paste(BLOCK[!keepb], collapse = ", "),
        " (they vary within a plant)")
  say("  blocks kept: ", if (any(keepb)) paste(BLOCK[keepb], collapse = " + ") else "(none)")
  fit_pl <- if (STAT == "spearman") {
    fit_blocked(row_midranks(P), rank(yp), bp[keepb])
  } else {
    fit_blocked(P, yp, bp[keepb])
  }
  PLANT_DF <- fit_pl$df
  say("  residual df = ", fit_pl$df)
  out[, `:=`(r_plant = fit_pl$r, p_plant = fit_pl$p,
             padj_plant = p.adjust(fit_pl$p, method = "BH"))]
  PLANT_R <- cor(out$r_partial, out$r_plant)
  say(sprintf("  agreement with the blocked model: r = %+.4f over %s genes",
              PLANT_R, fmt_n(nrow(out))))
  if (PLANT_R < 0.8)
    say("  WARNING: the two disagree. The blocked p-values are pseudoreplicated ",
        "and should not be read without this column.")
}

# --- the correctness gate, fatal ----------------------------------------------
# The vectorised solver is the whole basis for not looping over 170,135 genes, so
# it is required to reproduce lm() exactly before any number below it is reported.
say("")
verify_against_lm(V,  y,        blocks, fit_p, "blocked pearson")
verify_against_lm(Rk, rank(y),  blocks, fit_s, "blocked spearman")
if (STAT == "spearman") {
  verify_against_lm(Rk, rank(y), list(), fit_m, "marginal spearman")
} else {
  verify_against_lm(V, y, list(), fit_m, "marginal pearson")
}

# And the Spearman path must be the rank correlation it claims to be: with no
# block it is cor.test(method="spearman", exact=FALSE) exactly.
if (length(blocks)) {
  ct <- suppressWarnings(cor.test(V[1, ], y, method = "spearman", exact = FALSE))
  f0 <- fit_blocked(Rk[1, , drop = FALSE], rank(y), list())
  if (abs(f0$p - ct$p.value) > 1e-8)
    stop("the unblocked Spearman path does not reproduce cor.test()", call. = FALSE)
  say("unblocked Spearman reproduces cor.test(exact=FALSE)")
}

# --- out ----------------------------------------------------------------------
setorder(out, pval)
write_tsv(out, OUT_FILE)

sel    <- out[responsive == TRUE, .N]
sel_m  <- out[responsive_marginal == TRUE, .N]
say("")
say(sprintf("genes tested:                    %s", fmt_n(nrow(out))))
say(sprintf("  padj <= %.2f:                    %s", PADJ_THR, fmt_n(out[padj <= PADJ_THR, .N])))
say(sprintf("  |r_partial| >= %.2f:              %s", R_THR, fmt_n(out[abs(r_partial) >= R_THR, .N])))
say(sprintf("  SELECTED, blocked (both):       %s", fmt_n(sel)))
say(sprintf("  SELECTED, marginal (both):      %s", fmt_n(sel_m)))
say(sprintf("  |r| at the BH boundary:         %s",
            { o <- out[padj <= PADJ_THR]; if (!nrow(o)) "n/a" else sprintf("%.3f", min(abs(o$r_partial))) }))
say(sprintf("  raw p <= 0.05 (uncorrected):    %s", fmt_n(out[pval <= 0.05, .N])))

# Blocking is the same rule on a better-specified model, so a gene the MARGINAL
# test called responsive should still be responsive once a variance component is
# removed. If one is lost, the block is absorbing signal rather than noise, and
# that has to be seen rather than averaged away.
lost <- out[responsive_marginal == TRUE & responsive == FALSE, .N]
if (lost) {
  say(sprintf("  NOTE: %s marginal call(s) are NOT selected by the blocked model",
              fmt_n(lost)), " -- the block may be absorbing signal, not only noise")
} else {
  say("  blocked is a strict superset of the marginal rule on these genes")
}

# The histogram is the diagnostic this whole stage was built around: if blocking
# has done its job the rising right tail of the marginal test is gone.
hist_of <- function(p, label) {
  br <- table(cut(p, breaks = seq(0, 1, 0.1)))
  say("  ", label, ": ", paste(sprintf("%.1f%%", 100 * as.numeric(br) / length(p)),
                               collapse = " "))
  100 * as.numeric(br)[10] / length(p)
}
say("  p-value histogram (flat or decreasing = a well-behaved null):")
last_m <- hist_of(out$p_marginal, "marginal")
last_b <- hist_of(out$pval,       "blocked ")
if (last_b > 11)
  say("  WARNING: the blocked last decile is ", sprintf("%.1f%%", last_b),
      " -- the null is still not uniform, so something remains unmodelled.")

summ <- data.table(
  metric = c("study", "model", "statistic", "gene_universe", "genes_tested",
             "n_samples", "residual_df",
             "padj_le_thr", "r_ge_thr", "selected", "bh_boundary_abs_r",
             "raw_p_le_0.05", "padj_threshold", "r_threshold",
             "selected_marginal", "marginal_lost_under_blocking",
             "last_decile_pct_marginal", "last_decile_pct_blocked",
             "plant_level_df", "plant_level_selected", "plant_vs_blocked_r"),
  value = c(STUDY, out$model[1], STAT,
            if (nzchar(GENE_FILTER)) basename(GENE_FILTER) else "all VST genes",
            nrow(out), ncol(V), prim$df,
            out[padj <= PADJ_THR, .N], out[abs(r_partial) >= R_THR, .N], sel,
            { o <- out[padj <= PADJ_THR]; if (!nrow(o)) NA else round(min(abs(o$r_partial)), 4) },
            out[pval <= 0.05, .N], PADJ_THR, R_THR,
            sel_m, lost, round(last_m, 2), round(last_b, 2),
            PLANT_DF,
            if (nzchar(PLANT_FROM)) out[padj_plant <= PADJ_THR & abs(r_plant) >= R_THR, .N] else NA,
            if (nzchar(PLANT_FROM)) round(PLANT_R, 4) else NA))
write_tsv(summ, sub("\\.tsv$", "_summary.tsv", OUT_FILE))
say("done: ", STUDY)
