#!/usr/bin/env Rscript
# =============================================================================
# 53_module_trait_blocked.r — module eigengene vs trait, with the design in the
# model.
#
# WHAT CHANGES FROM 19. 19_module_trait_spearman.r is a MARGINAL Spearman: rho
# between a module's eigengene and the encoded trait, with nothing else in the
# model. This project's own QC records what that leaves in the residual:
#
#     purple    genotype R2 on PC1 = 0.999      sugarcane genotype R2 = 0.998
#     sugarcane segment  R2 on PC2 = 0.802
#
# so the single largest source of variance in either matrix sits in the residual
# of every module-level nitrogen test, inflating the residual variance and
# deflating every t. At GENE level the same correction moved purple from 30
# responsive genes to 2,331 and flattened a p-value histogram that had been
# RISING in its last decile -- which a mixture of a uniform null and real signal
# cannot do (31_gene_trait_blocked.r, branch blocked-gene-trait).
#
#   sugarcane   eigengene ~ genotype + segment + N      n = 48, resid df 42
#   purple      eigengene ~ genotype + N                n = 18, resid df 15
#
# WHY SPEARMAN IS STILL PRIMARY. Unchanged from 19:8-27 and already accepted by
# this project. Purple's trait is an ORDINAL 0/2/6 mM dose and Pearson reads
# those spacings literally -- it asks whether the 2->6 response is exactly twice
# the 0->2 one, which nothing in the design justifies. Sugarcane's is two-level,
# where Spearman on eigengene ranks is the rank-biserial correlation. "Blocked
# Spearman" is the blocked model fitted to MIDRANKS on both sides: rank, then
# residualise. Blocked Pearson is computed and written beside it so the choice is
# visible in the table rather than buried in the code.
#
# THE MARGINAL RHO IS CARRIED, not discarded. Blocking is not a looser threshold
# -- the FDR level and the |rho| floor are the same values 19 uses, and only the
# model changes -- so the table reports both and the effect of blocking on the
# responsive set can be read off directly.
#
# SUGARCANE'S PSEUDOREPLICATION CONTROL IS NOT OPTIONAL. Its 48 libraries are 12
# plants x 4 leaf segments: repeated measures, not 48 replicates. Putting segment
# in as a fixed block removes the segment MEANS but does not account for
# within-plant correlation, so a PLANT-LEVEL run (segments averaged, n = 12,
# resid df 9) is computed alongside and correlated against the blocked fit. If
# the two disagree the blocked p-values are not to be trusted, and carrying it is
# the point.
#
# THE NULL PERMUTES WITHIN BLOCK. 19 shuffles trait labels freely, which under a
# blocked model would break the design the model conditions on and produce an
# optimistic null. Labels are permuted within each block stratum instead, so the
# null respects the same structure the test does.
#
# Writes the same contract 19 writes, so every downstream consumer keeps working:
#   module_trait_<study>.tsv          module, trait, n, rho, pval, padj,
#                                     responsive, direction  (+ blocked/marginal
#                                     columns appended)
#   module_trait_<study>.null.tsv     one row per permutation
#   module_trait_<study>.summary.json
#
# RUN: through run.sh  ->  ./run.sh moduletrait sugarcane
# =============================================================================
suppressMessages({ library(data.table) })
source(file.path(dirname(sub("--file=", "",
       grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "lib", "common.R"))

STUDY        <- env_req("CLEAN_STUDY")
EIG_PREFIX   <- env_req("CLEAN_EIGENGENE_PREFIX")
META_FILE    <- env_req("CLEAN_META")
TRAIT_SPEC   <- env_req("CLEAN_TRAITS")
OUT_FILE     <- env_req("CLEAN_OUT_FILE")
SELECT_TRAIT <- env_opt("CLEAN_SELECT_TRAIT", "treatment")
BLOCK        <- strsplit(trimws(env_opt("CLEAN_BLOCK", "genotype")), "[ ,]+")[[1L]]
STAT         <- env_opt("CLEAN_BLOCKED_STAT", "spearman")
PLANT_FROM   <- env_opt("CLEAN_PLANT_FROM", "")
R_THR        <- env_num("CLEAN_MODULE_R_THR", 0.6)
PADJ_THR     <- env_num("CLEAN_MODULE_PADJ_THR", 0.05)
N_PERM       <- as.integer(env_num("CLEAN_MODULE_PERM", 1000))
SEED         <- as.integer(env_num("CLEAN_SEED", 1))
setDTthreads(as.integer(env_num("CLEAN_CORES", 32)))

if (!STAT %in% c("pearson", "spearman"))
  stop("CLEAN_BLOCKED_STAT must be pearson or spearman (got '", STAT, "')", call. = FALSE)
set.seed(SEED)

banner(sprintf("blocked module-trait: %s  [eigengene ~ %s + %s, %s primary]",
               STUDY, paste(BLOCK, collapse = " + "), SELECT_TRAIT, STAT))

# --- trait encoding, the same parser 07, 19 and 31 use -----------------------
parse_traits <- function(spec) {
  out <- list()
  for (blk in strsplit(spec, ";", fixed = TRUE)[[1L]]) {
    blk <- trimws(blk); if (!nzchar(blk)) next
    nm <- sub(":.*$", "", blk); body <- sub("^[^:]*:", "", blk)
    kv <- strsplit(strsplit(body, ",", fixed = TRUE)[[1L]], "=", fixed = TRUE)
    if (any(lengths(kv) != 2L)) stop("cannot parse trait spec '", blk, "'", call. = FALSE)
    v <- as.numeric(vapply(kv, `[`, "", 2L))
    names(v) <- trimws(vapply(kv, `[`, "", 1L))
    out[[trimws(nm)]] <- v
  }
  out
}
TRAITS <- parse_traits(TRAIT_SPEC)
if (!SELECT_TRAIT %in% names(TRAITS))
  stop("trait '", SELECT_TRAIT, "' not in CLEAN_TRAITS", call. = FALSE)
enc <- TRAITS[[SELECT_TRAIT]]

# --- inputs -------------------------------------------------------------------
# read_vst reads the eigengene container too: 14 writes <prefix>.f32 +
# <prefix>.genes.txt, exactly the shape the VST uses, which is why the blocked
# fit needs no reader of its own.
E <- read_vst(EIG_PREFIX)
modules <- rownames(E)
say(fmt_n(nrow(E)), " eigengenes x ", ncol(E), " samples")

meta <- fread(META_FILE, header = TRUE)
setnames(meta, tolower(names(meta)))
m <- meta[match(colnames(E), sample)]
if (anyNA(m$sample)) stop("eigengene samples missing from ", basename(META_FILE), call. = FALSE)
for (cn in c(SELECT_TRAIT, BLOCK))
  if (!cn %in% names(m))
    stop("column '", cn, "' not in ", basename(META_FILE),
         "\n  available: ", paste(names(m), collapse = ", "), call. = FALSE)

y <- unname(enc[as.character(m[[SELECT_TRAIT]])])
if (anyNA(y))
  stop("samples with an unencoded ", SELECT_TRAIT, ": ",
       paste(unique(as.character(m[[SELECT_TRAIT]])[is.na(y)]), collapse = ", "), call. = FALSE)
blocks <- setNames(lapply(BLOCK, function(b) factor(m[[b]])), BLOCK)

say("design: ", ncol(E), " samples")
for (b in BLOCK)
  say("  ", b, ": ", paste(sprintf("%s=%d", names(table(blocks[[b]])), table(blocks[[b]])),
                           collapse = "  "))
say("  ", SELECT_TRAIT, ": ",
    paste(sprintf("%g=%d", as.numeric(names(table(y))), table(y)), collapse = "  "))

# --- both statistics, always --------------------------------------------------
say("")
say("fitting eigengene ~ ", paste(BLOCK, collapse = " + "), " + ", SELECT_TRAIT)
fit_p <- fit_blocked(E, y, blocks)
say("  blocked pearson : residual df = ", fit_p$df)
Rk <- row_midranks(E)
fit_s <- fit_blocked(Rk, rank(y), blocks)
say("  blocked spearman: residual df = ", fit_s$df)

verify_against_lm(E,  y,       blocks, fit_p, "blocked pearson")
verify_against_lm(Rk, rank(y), blocks, fit_s, "blocked spearman")

# The Spearman path must be the rank correlation it claims to be: with no block
# it is cor.test(method="spearman", exact=FALSE) exactly.
if (length(blocks)) {
  ct <- suppressWarnings(cor.test(E[1, ], y, method = "spearman", exact = FALSE))
  f0 <- fit_blocked(Rk[1, , drop = FALSE], rank(y), list())
  if (abs(f0$p - ct$p.value) > 1e-8)
    stop("the unblocked Spearman path does not reproduce cor.test()", call. = FALSE)
  say("unblocked Spearman reproduces cor.test(exact = FALSE)")
}

prim <- if (STAT == "spearman") fit_s else fit_p

# THE MARGINAL FIT, ON THESE EIGENGENES. It is the same solver with NO blocks, so
# blocked-vs-marginal compares two MODELS of one module set.
#
# It is computed here rather than read from the previous run's table, and that is
# not a convenience: module names are positional (Module_%03d, largest first), so
# "Module_001" under a different clustering is a different set of genes. Joining
# the two tables on that name silently compares unrelated modules -- it runs
# cleanly and answers nothing. Measured on this data: old Module_001 held 19,604
# genes, the new one holds 23,439.
fit_m <- if (STAT == "spearman") {
  fit_blocked(Rk, rank(y), list())
} else fit_blocked(E, y, list())
say("  marginal (no blocks): residual df = ", fit_m$df)
padj_m <- p.adjust(fit_m$p, method = "BH")

# --- the plant-level control (sugarcane) --------------------------------------
plant_rho <- rep(NA_real_, nrow(E)); plant_df <- NA_integer_
if (nzchar(PLANT_FROM)) {
  # Plant id = every design column EXCEPT the one that varies within plant.
  # The trailing _1/_2/_3 of the sample name is the plant within a genotype x N
  # cell, so grouping on (all blocks but the within-plant one) + trait + that
  # suffix recovers the 12 plants.
  within <- setdiff(BLOCK, PLANT_FROM)
  rep_id <- sub("^.*_([0-9]+)$", "\\1", colnames(E))
  key <- paste(m[[PLANT_FROM]], y, rep_id, sep = "|")
  say("")
  say("plant-level control: ", length(unique(key)), " plants from ", ncol(E), " libraries",
      if (length(within)) paste0(" (averaging over ", paste(within, collapse = ", "), ")") else "")
  if (length(unique(key)) < ncol(E)) {
    idx <- split(seq_len(ncol(E)), key)
    Ep  <- vapply(idx, function(j) rowMeans(E[, j, drop = FALSE]), numeric(nrow(E)))
    yp  <- vapply(idx, function(j) y[j[1]], numeric(1))
    bp  <- setNames(lapply(PLANT_FROM, function(b)
             factor(vapply(idx, function(j) as.character(m[[b]])[j[1]], character(1)))),
             PLANT_FROM)
    fp <- if (STAT == "spearman") {
      fit_blocked(row_midranks(Ep), rank(yp), bp)
    } else fit_blocked(Ep, yp, bp)
    plant_rho <- fp$r; plant_df <- fp$df
    ok <- is.finite(plant_rho) & is.finite(prim$r)
    say(sprintf("  plant-level resid df = %d; agreement with the blocked fit: r = %+.4f",
                fp$df, cor(prim$r[ok], plant_rho[ok])))
  } else say("  sample names do not collapse into plants; control skipped")
}

# --- the null: permute WITHIN block -------------------------------------------
strata <- if (length(blocks)) {
  interaction(as.data.frame(blocks), drop = TRUE)
} else factor(rep(1L, ncol(E)))
say("")
say("null: ", fmt_n(N_PERM), " permutations of ", SELECT_TRAIT,
    " WITHIN ", nlevels(strata), " block strata")
grp <- split(seq_along(y), strata)
null_counts <- integer(N_PERM)
for (i in seq_len(N_PERM)) {
  yp <- y
  for (g in grp) yp[g] <- y[sample(g)]
  f <- if (STAT == "spearman") {
    fit_blocked(Rk, rank(yp), blocks)
  } else fit_blocked(E, yp, blocks)
  padj <- p.adjust(f$p, method = "BH")
  null_counts[i] <- sum(padj <= PADJ_THR & abs(f$r) >= R_THR, na.rm = TRUE)
}

# --- out ----------------------------------------------------------------------
padj_b <- p.adjust(prim$p, method = "BH")
out <- data.table(
  module = modules, trait = SELECT_TRAIT, n = ncol(E),
  rho = prim$r, pval = prim$p, padj = padj_b,
  responsive = padj_b <= PADJ_THR & abs(prim$r) >= R_THR,
  direction = ifelse(prim$r > 0, "up", "down"),
  # everything below is reported beside the primary call, not used to make it
  blocked_stat = STAT, blocked_df = prim$df,
  rho_blocked_spearman = fit_s$r, padj_blocked_spearman = p.adjust(fit_s$p, "BH"),
  rho_blocked_pearson  = fit_p$r, padj_blocked_pearson  = p.adjust(fit_p$p, "BH"),
  rho_plant = plant_rho, plant_df = plant_df)
out[!is.finite(rho), `:=`(responsive = FALSE, direction = NA_character_)]

# The marginal columns come from fit_m above -- same modules, same samples.
out[, `:=`(rho_marginal = fit_m$r, padj_marginal = padj_m,
           responsive_marginal = padj_m <= PADJ_THR & abs(fit_m$r) >= R_THR)]
out[!is.finite(rho_marginal), responsive_marginal := FALSE]
n_marg <- sum(out$responsive_marginal, na.rm = TRUE)

out <- out[order(padj, -abs(rho))]
write_tsv(out, OUT_FILE)
write_tsv(data.table(perm = seq_len(N_PERM), n_responsive = null_counts),
          sub("\\.tsv$", ".null.tsv", OUT_FILE))

n_resp <- sum(out$responsive, na.rm = TRUE)
jf <- sub("\\.tsv$", ".summary.json", OUT_FILE)
writeLines(sprintf(paste0(
  '{\n  "study": "%s",\n  "trait": "%s",\n  "model": "%s",\n',
  '  "primary_stat": "%s",\n  "residual_df": %d,\n  "n_samples": %d,\n',
  '  "n_modules": %d,\n  "n_responsive": %d,\n  "n_responsive_marginal": %d,\n',
  '  "r_threshold": %g,\n  "padj_threshold": %g,\n  "n_perm": %d,\n',
  '  "null_mean": %.3f,\n  "null_max": %d,\n  "null_ge_observed": %d,\n',
  '  "empirical_p": %.4g,\n  "plant_df": %s\n}'),
  STUDY, SELECT_TRAIT, paste("eigengene ~", paste(c(BLOCK, SELECT_TRAIT), collapse = " + ")),
  STAT, prim$df, ncol(E), nrow(out), n_resp, n_marg,
  R_THR, PADJ_THR, N_PERM, mean(null_counts), max(null_counts),
  sum(null_counts >= n_resp), (sum(null_counts >= n_resp) + 1) / (N_PERM + 1),
  ifelse(is.na(plant_df), "null", as.character(plant_df))), jf)

say("")
say(STUDY, ": ", fmt_n(n_resp), " responsive modules of ", fmt_n(nrow(out)),
    sprintf(" (%.2f%%)", 100 * n_resp / nrow(out)))
say("  null over ", fmt_n(N_PERM), " within-block permutations: mean ",
    sprintf("%.2f", mean(null_counts)), ", max ", max(null_counts))
say(sprintf("  permutations reaching %s: %d of %s  (empirical p <= %.4g)",
            fmt_n(n_resp), sum(null_counts >= n_resp), fmt_n(N_PERM),
            (sum(null_counts >= n_resp) + 1) / (N_PERM + 1)))
if (n_marg > 0) {
  both <- sum(out$responsive & out$responsive_marginal, na.rm = TRUE)
  say("  marginal responsive: ", fmt_n(n_marg), "; blocked: ", fmt_n(n_resp),
      "; in both: ", fmt_n(both))
  say("  blocked is a superset of marginal: ",
      ifelse(both == n_marg, "YES", paste0("NO -- ", fmt_n(n_marg - both),
             " marginal calls do not survive the design")))
}
say("done: ", STUDY)
