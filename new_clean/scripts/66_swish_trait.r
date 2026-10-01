#!/usr/bin/env Rscript
# =============================================================================
# 66_swish_trait.r -- the nitrogen test again, with quantification uncertainty
# propagated, run BESIDE the blocked fit rather than instead of it.
#
# WHY A SECOND TEST. 31_gene_trait_blocked.r is the primary nitrogen test and
# stays so: it is design-aware, verified against lm() to 1e-8, and every
# responsive count in docs/ comes from it. What it cannot do is know how sure
# salmon was about a gene's counts. In an 8-12x polyploid that is not a detail --
# 45.5% of purple copies and 19.6% of sugarcane copies have NO unique 31-mer
# (docs/dnds.md), so a "responsive" call can in principle be made of reads that
# could equally have been assigned to a sibling.
#
# Swish (Zhu, Srivastava, Ibrahim, Patro & Love, NAR 2019) answers the same
# question over the 30 Gibbs replicates, so a gene passing BOTH tests is
# responsive in a way that does not depend on how the reads were split. THE
# USEFUL OUTPUT IS THE AGREEMENT TABLE, not a new gene list.
#
# IT MATTERS MOST FOR PURPLE. Its blocked p-value histogram is still not flat --
# the last decile sits at 10.2% and the shape dips to 5.4% and rises again, so
# genotype plus nitrogen does not account for everything at n = 18. Quantification
# uncertainty is one candidate for part of that residual, and it is the one the
# blocked fit cannot model.
#
# THE DESIGNS MIRROR THE BLOCKED FIT AS CLOSELY AS SWISH ALLOWS
#
#   purple      x = treatment (0/2/6 mM, numeric), cov = genotype, cor = spearman
#               The blocked test is Spearman here for the same reason: the dose is
#               ORDINAL in a stress-control-stress design where 2 mM is the
#               control, so Pearson's literal spacing is unjustified.
#   sugarcane   x = treatment (high/low), cov = genotype:segment interaction
#               Swish takes ONE stratifying covariate, so the two blocks of the
#               blocked model (genotype + segment) are crossed into 8 strata of 6
#               libraries. That is a stratified Wilcoxon, the natural analogue.
#
# THE UNIVERSES DIFFER AND THAT IS REPORTED, NOT HIDDEN. labelKeep() drops genes
# with fewer than minCount in minN samples, so Swish tests fewer genes than the
# blocked fit's network-node universe. The agreement table is computed on the
# INTERSECTION and the gap is printed, because a gene absent from one test is not
# evidence against the other.
#
# RUN: through run.sh -> ./run.sh swishtrait <study>
# =============================================================================

suppressMessages({
  library(tximport)
  library(fishpond)
  library(SummarizedExperiment)
  library(data.table)
})

source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE),
                                                 value = TRUE)[1])), "lib", "common.R"))

STUDY      <- env_req("CLEAN_STUDY")
QUANT_DIR  <- env_req("CLEAN_QUANT_DIR")
META_FILE  <- env_req("CLEAN_META")
TRAIT_FILE <- env_req("CLEAN_GENETRAIT")
TRAITS     <- env_req("CLEAN_TRAITS")
OUT_FILE   <- env_req("CLEAN_OUT_FILE")
TRAIT      <- env_opt("CLEAN_SELECT_TRAIT", "treatment")
BLOCK      <- env_list("CLEAN_BLOCK", "genotype")
STRIP_VER  <- env_flag("CLEAN_STRIP_VERSION")
CORMODE    <- env_opt("CLEAN_SWISH_COR", "none")
NPERMS     <- as.integer(env_num("CLEAN_SWISH_NPERMS", 100))
QTHR       <- env_num("CLEAN_SWISH_QTHR", 0.05)
SEED       <- as.integer(env_num("CLEAN_SEED", 1))

if (!CORMODE %in% c("none", "spearman", "pearson"))
  stop("CLEAN_SWISH_COR must be none|spearman|pearson (got '", CORMODE, "')",
       call. = FALSE)

banner(paste("swish nitrogen test:", STUDY))
set.seed(SEED)

# --- the trait encoding, read from the same spec the blocked fit uses -------
# Same mini-DSL as 31_gene_trait_blocked.r:92-107 -- "trait:LEVEL=value,...;..."
parse_traits <- function(spec) {
  out <- list()
  for (part in strsplit(spec, ";", fixed = TRUE)[[1]]) {
    if (!nzchar(trimws(part))) next
    kv <- strsplit(part, ":", fixed = TRUE)[[1]]
    lv <- strsplit(kv[2], ",", fixed = TRUE)[[1]]
    pairs <- do.call(rbind, lapply(lv, function(x)
      strsplit(x, "=", fixed = TRUE)[[1]]))
    out[[trimws(kv[1])]] <- setNames(as.numeric(pairs[, 2]), trimws(pairs[, 1]))
  }
  out
}
enc_all <- parse_traits(TRAITS)
if (!TRAIT %in% names(enc_all))
  stop("trait '", TRAIT, "' is not in CLEAN_TRAITS", call. = FALSE)
enc <- enc_all[[TRAIT]]

meta <- fread(META_FILE, header = TRUE)
setnames(meta, tolower(names(meta)))
samples <- as.character(meta$sample)
y_num <- unname(enc[as.character(meta[[TRAIT]])])
if (anyNA(y_num))
  stop("samples with an unencoded ", TRAIT, ": ",
       paste(unique(meta[[TRAIT]][is.na(y_num)]), collapse = ", "), call. = FALSE)

# One stratifying covariate: the blocks crossed, since swish takes only one.
strat <- {
  if (length(BLOCK) == 1L) {
    factor(meta[[BLOCK]])
  } else {
    droplevels(interaction(lapply(BLOCK, function(b) factor(meta[[b]])),
                           drop = TRUE, sep = ":"))
  }
}
say(sprintf("%d samples | trait %s -> %s | %d strata of %s",
            length(samples), TRAIT,
            paste(sort(unique(y_num)), collapse = "/"),
            nlevels(strat), paste(table(strat), collapse = "/")))
if (any(table(strat, y_num) == 0))
  stop("a stratum is missing a trait level -- swish cannot stratify on it",
       call. = FALSE)

# --- load the inferential replicates ----------------------------------------
qf <- file.path(QUANT_DIR, samples, "quant.sf"); names(qf) <- samples
if (!all(file.exists(qf)))
  stop("missing quant.sf for: ",
       paste(samples[!file.exists(qf)], collapse = ", "), call. = FALSE)
t2g <- fread(file.path(QUANT_DIR, "salmon.merged.tx2gene.tsv"), header = TRUE)

say("tximport with inferential replicates (the slow step)")
txi <- tximport(qf, type = "salmon", tx2gene = as.data.frame(t2g[, 1:2]),
                countsFromAbundance = "lengthScaledTPM", dropInfReps = FALSE)
if (is.null(txi$infReps))
  stop("no inferential replicates in the quantification", call. = FALSE)
N_REP <- ncol(txi$infReps[[1]])

genes <- rownames(txi$counts)
if (STRIP_VER) {
  genes <- strip_version(genes)
  if (anyDuplicated(genes))
    stop("stripping versions collides gene ids", call. = FALSE)
}

# computeInfRV/swish want infRep1..K as genes x samples; tximport gives
# genes x replicates per SAMPLE, so transpose the nesting once.
mk <- function(mats, ids) {
  setNames(lapply(seq_len(N_REP), function(k) {
    A <- vapply(mats, function(M) M[, k], numeric(nrow(mats[[1]])))
    dimnames(A) <- list(ids, samples); A
  }), paste0("infRep", seq_len(N_REP)))
}
cm <- txi$counts; dimnames(cm) <- list(genes, samples)
lm_ <- txi$length; dimnames(lm_) <- list(genes, samples)
se <- SummarizedExperiment(
  assays = c(list(counts = cm, length = lm_), mk(txi$infReps, genes)),
  colData = DataFrame(row.names = samples, sample = samples,
                      trait = y_num, stratum = strat))
rm(txi); invisible(gc())
say(sprintf("SummarizedExperiment: %s genes x %d samples x %d replicates",
            fmt_n(nrow(se)), ncol(se), N_REP))

# --- the swish flow, in its published order ---------------------------------
se <- scaleInfReps(se, quiet = TRUE)
se <- labelKeep(se)
n_all <- nrow(se)
se <- se[mcols(se)$keep, ]
say(sprintf("labelKeep: %s of %s genes tested (%s dropped as low-count)",
            fmt_n(nrow(se)), fmt_n(n_all), fmt_n(n_all - nrow(se))))

# HOW THE DESIGN ENTERS, AND THE ONE PLACE SWISH CANNOT MATCH THE BLOCKED FIT.
# swish's own guards (fishpond::swish) are:
#     if (correlation & !paired) stopifnot(!cov_given)
#     if (correlation & paired)  stopifnot(cov_given)
# so a CONTINUOUS trait admits a stratifying covariate only in a PAIRED design.
# Purple is 2 genotypes x 3 doses x 3 replicates -- not paired -- so the dose
# correlation cannot be blocked on genotype the way 31_gene_trait_blocked.r blocks
# it. Rather than quietly drop the genotype term, purple is run three ways and all
# three are written: POOLED over both genotypes, and once WITHIN each genotype.
# "Responsive in both genotypes" is then a visible, stronger claim than the pooled
# call, and it is the same shape as the per-genotype estimates
# 30_gene_trait_ushape.r already carries. Sugarcane's trait is two-level, so
# correlation is off and cov = genotype:segment stratifies as intended.
pull <- function(obj, tag) {
  mc <- mcols(obj)
  keep <- intersect(c("stat", "log2FC", "pvalue", "locfdr", "qvalue"), colnames(mc))
  d <- as.data.table(as.data.frame(mc[, keep, drop = FALSE]))
  if (nzchar(tag)) setnames(d, keep, paste0(keep, "_", tag))
  d[, gene := rownames(obj)]
  d
}

if (CORMODE == "none") {
  say(sprintf("swish: x=trait, cov=stratum, two-group, %d permutations", NPERMS))
  se <- swish(se, x = "trait", cov = "stratum", nperms = NPERMS, quiet = TRUE)
  sw <- pull(se, "")
} else {
  say(sprintf("swish: x=trait, cor=%s, POOLED (swish forbids cov with cor), %d permutations",
              CORMODE, NPERMS))
  se <- swish(se, x = "trait", cor = CORMODE, nperms = NPERMS, quiet = TRUE)
  sw <- pull(se, "")
  for (lv in levels(strat)) {
    sel <- which(strat == lv)
    say(sprintf("  within %s (n=%d)", lv, length(sel)))
    s1 <- swish(se[, sel], x = "trait", cor = CORMODE, nperms = NPERMS, quiet = TRUE)
    sw <- merge(sw, pull(s1, make.names(lv)), by = "gene", all.x = TRUE)
  }
  qcols <- grep("^qvalue_", names(sw), value = TRUE)
  if (length(qcols) >= 2L) {
    sw[, n_genotypes_responsive :=
         rowSums(as.data.table(lapply(.SD, function(q) is.finite(q) & q <= QTHR))),
       .SDcols = qcols]
    say(sprintf("  responsive in BOTH genotypes: %s genes",
                fmt_n(sw[n_genotypes_responsive == length(qcols), .N])))
  }
}
if (!"gene" %in% names(sw)) sw[, gene := rownames(se)]
sw[, meanInfRV := as.numeric(mcols(se)$meanInfRV)[match(gene, rownames(se))]]
setcolorder(sw, "gene")
sw[, swish_responsive := is.finite(qvalue) & qvalue <= QTHR]
say(sprintf("swish calls %s genes at qvalue <= %.2f", fmt_n(sw[, sum(swish_responsive)]), QTHR))

# --- the agreement table, on the intersection -------------------------------
bl <- fread(TRAIT_FILE, select = c("gene", "r_partial", "padj", "responsive"))
bl[, gene := strip_version(gene)]
both <- merge(sw, bl, by = "gene")
say(sprintf("\nuniverses: swish %s | blocked %s | intersection %s",
            fmt_n(nrow(sw)), fmt_n(nrow(bl)), fmt_n(nrow(both))))
say(sprintf("  blocked-responsive genes absent from the swish universe: %s",
            fmt_n(bl[responsive %in% c(TRUE, "TRUE") & !gene %chin% sw$gene, .N])))

agree <- both[, .N, by = .(blocked = responsive %in% c(TRUE, "TRUE"),
                           swish = swish_responsive)][order(-blocked, -swish)]
banner("agreement")
for (i in seq_len(nrow(agree)))
  say(sprintf("  blocked %-5s swish %-5s %9s", agree$blocked[i], agree$swish[i],
              fmt_n(agree$N[i])))
nb <- both[(responsive %in% c(TRUE, "TRUE")), .N]
ns <- both[(swish_responsive), .N]
nboth <- both[(responsive %in% c(TRUE, "TRUE")) & (swish_responsive), .N]
say(sprintf("\nCONFIDENT SET -- called by both: %s", fmt_n(nboth)))
if (nb > 0)
  say(sprintf("  that is %.1f%% of the %s blocked calls and %.1f%% of the %s swish calls",
              100 * nboth / nb, fmt_n(nb), 100 * nboth / max(ns, 1), fmt_n(ns)))
if (both[, sum(is.finite(r_partial) & is.finite(stat))] > 2)
  say(sprintf("  rank correlation of the two statistics: %.4f",
              cor(both$r_partial, both$stat, method = "spearman",
                  use = "complete.obs")))

both[, confident := (responsive %in% c(TRUE, "TRUE")) & swish_responsive]
ensure_dir(dirname(OUT_FILE))
write_tsv(both[order(qvalue, -abs(r_partial))], OUT_FILE)
write_tsv(agree, sub("\\.tsv$", "_agreement.tsv", OUT_FILE))
write_tsv(data.table(
  metric = c("study", "n_samples", "n_replicates", "trait", "cor_mode",
             "n_strata", "nperms", "qvalue_threshold",
             "n_genes_quantified", "n_genes_swish_tested", "n_genes_blocked",
             "n_intersection", "n_swish_responsive", "n_blocked_responsive",
             "n_confident_both", "pct_of_blocked_confirmed",
             "blocked_responsive_outside_swish_universe",
             "spearman_rpartial_vs_swishstat"),
  value = c(STUDY, ncol(se), N_REP, TRAIT, CORMODE, nlevels(strat), NPERMS,
            sprintf("%.2f", QTHR), fmt_n(n_all), fmt_n(nrow(sw)), fmt_n(nrow(bl)),
            fmt_n(nrow(both)), fmt_n(ns), fmt_n(nb), fmt_n(nboth),
            if (nb > 0) sprintf("%.2f", 100 * nboth / nb) else "NA",
            fmt_n(bl[responsive %in% c(TRUE, "TRUE") & !gene %chin% sw$gene, .N]),
            sprintf("%.4f", cor(both$r_partial, both$stat, method = "spearman",
                                use = "complete.obs")))),
  sub("\\.tsv$", "_summary.tsv", OUT_FILE))

say("done")
