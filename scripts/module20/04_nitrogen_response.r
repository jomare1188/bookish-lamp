# ============================================================================
# 04_nitrogen_response.r — are the Module-20 orthologs nitrogen-responsive in
# OUR re-quantification, in each network's own data?
#
# The two designs differ, so the test differs:
#
#   sugarcane (Muñoz, 48 libs): 2 genotypes x 2 N levels (High/Low) x leaf.
#     The three leaf segments (Leaf Apex / Leaf Base / Leaf) are treated as ONE
#     leaf tissue per the project's decision; segment is kept only as a
#     sensitivity covariate. With TWO N levels there is no U-shape problem, so a
#     per-genotype Welch test is the right instrument — but it must be run
#     WITHIN genotype, because Muñoz's whole point is a responsive vs
#     non-responsive genotype contrast that pooling would cancel.
#
#   purple (Kiet, 18 leaf libs): 2 genotypes x (0N, 2N, 6N).
#     Three levels, and the companion study describes a non-monotonic response,
#     so this reuses the scripts/myb61/08 machinery: per-genotype ANOVA plus a
#     U-shape contrast c(+1,-2,+1) and a linear contrast.
#
# Independence: the mapped genes are polyploid haplotype copies of only 3-4
# distinct loci, so per-gene p-values are NOT independent. Results are therefore
# also summarised per LOCUS, and that is the level to interpret.
#
# Figures live in 05_module20_heatmaps.r (ComplexHeatmap); this step is tests only.
#
# RUN: /home/genomics/miniconda3/envs/r_env/bin/Rscript 04_nitrogen_response.r
# ============================================================================

suppressMessages(library(data.table))

base   <- "/dados04/jorge/comparative_saccharum"
outdir <- file.path(base, "GET_TFS/new/results/module20")
COL    <- c(sugarcane = "#1b7837", purple = "#762a83")

tab <- fread(file.path(outdir, "module20_network_table.tsv"))

# ---------------------------------------------------------------------------
# expression matrices
# ---------------------------------------------------------------------------
load_expr <- function(tpmf, metaf, sp) {
  tpm <- fread(tpmf); setnames(tpm, 1, "gene")
  tpm[, gene := sub("\\.v[0-9.]+$", "", gene)]
  meta <- fread(metaf)
  if (sp == "sugarcane") {
    m <- meta[, .(sample, genotype, treatment, segment = tissue)]
    m[, nlev := fifelse(grepl("Low", treatment), "LowN", "HighN")]
    m[, nlev := factor(nlev, levels = c("LowN", "HighN"))]
    m[, genotype := fifelse(genotype == "RB975375", "RB975375 (responsive)",
                                                    "RB937570 (non-responsive)")]
  } else {
    m <- meta[tissue == "leaf", .(sample, genotype, treatment)]
    m[, segment := "leaf"]
    m[, nlev := factor(treatment, levels = c("0N", "2N", "6N"))]
  }
  samp <- intersect(names(tpm), m$sample)
  stopifnot(length(samp) == nrow(m))
  E <- melt(tpm[gene %in% tab[species == sp, gene], c("gene", samp), with = FALSE],
            id.vars = "gene", variable.name = "sample", value.name = "tpm")
  E[, sample := as.character(sample)]
  E <- merge(E, m, by = "sample")
  E[, l2 := log2(tpm + 1)][, species := sp]
  E[]
}

E_sc <- load_expr(file.path(base, "run1/salmon/salmon.merged.gene_tpm.tsv"),
                  file.path(base, "samplesheet.csv"), "sugarcane")
E_pu <- load_expr(file.path(base, "china/run2_onlyL/salmon/salmon.merged.gene_tpm.tsv"),
                  file.path(base, "china/samplesheet_china.csv"), "purple")

cat(sprintf("sugarcane: %d genes x %d libs | purple: %d genes x %d libs\n",
            uniqueN(E_sc$gene), uniqueN(E_sc$sample),
            uniqueN(E_pu$gene), uniqueN(E_pu$sample)))

# ---------------------------------------------------------------------------
# tests
# ---------------------------------------------------------------------------
contrast_test <- function(y, lev, cvec) {
  m <- tapply(y, lev, mean); n <- tapply(y, lev, length)
  s2 <- sum(tapply(y, lev, function(v) sum((v - mean(v))^2))) / (length(y) - nlevels(lev))
  if (!is.finite(s2) || s2 <= 0) return(list(est = NA_real_, p = NA_real_))
  est <- sum(cvec * m); se <- sqrt(s2 * sum(cvec^2 / n))
  list(est = est, p = 2 * pt(-abs(est / se), df = length(y) - nlevels(lev)))
}

# --- sugarcane: per-genotype High vs Low, plus genotype x N interaction -----
res_sc <- rbindlist(lapply(unique(E_sc$gene), function(g) {
  rbindlist(lapply(unique(E_sc$genotype), function(gt) {
    d <- E_sc[gene == g & genotype == gt]
    ok <- uniqueN(d$nlev) == 2 && sd(d$l2) > 0
    tt <- if (ok) tryCatch(t.test(l2 ~ nlev, data = d), error = function(e) NULL) else NULL
    data.table(species = "sugarcane", gene = g, genotype = gt,
               mean_tpm = mean(d$tpm),
               log2FC_HighminusLow = if (ok) diff(tapply(d$l2, d$nlev, mean)) else NA_real_,
               p = if (!is.null(tt)) tt$p.value else NA_real_)
  }))
}))
inter_sc <- rbindlist(lapply(unique(E_sc$gene), function(g) {
  d <- E_sc[gene == g]
  p <- if (sd(d$l2) > 0) tryCatch(anova(lm(l2 ~ genotype * nlev, data = d))$`Pr(>F)`[3],
                                  error = function(e) NA_real_) else NA_real_
  # sensitivity: same model with leaf segment as a blocking covariate
  ps <- if (sd(d$l2) > 0) tryCatch(anova(lm(l2 ~ segment + genotype * nlev, data = d))$`Pr(>F)`[4],
                                   error = function(e) NA_real_) else NA_real_
  data.table(gene = g, interaction_p = p, interaction_p_segment_adj = ps)
}))
res_sc <- merge(res_sc, inter_sc, by = "gene")

# --- purple: per-genotype ANOVA + U-shape + linear contrasts ----------------
res_pu <- rbindlist(lapply(unique(E_pu$gene), function(g) {
  rbindlist(lapply(unique(E_pu$genotype), function(gt) {
    d <- E_pu[gene == g & genotype == gt]; d[, nlev := droplevels(nlev)]
    ok <- uniqueN(d$nlev) == 3 && sd(d$l2) > 0
    a <- if (ok) tryCatch(anova(lm(l2 ~ nlev, data = d))$`Pr(>F)`[1], error = function(e) NA_real_) else NA_real_
    u <- if (ok) contrast_test(d$l2, d$nlev, c(1, -2, 1)) else list(est = NA_real_, p = NA_real_)
    l <- if (ok) contrast_test(d$l2, d$nlev, c(-1, 0, 1)) else list(est = NA_real_, p = NA_real_)
    data.table(species = "purple", gene = g, genotype = gt, mean_tpm = mean(d$tpm),
               anova_p = a, U_est = u$est, U_p = u$p, lin_est = l$est, lin_p = l$p)
  }))
}))

# BH within species x genotype
res_sc[, padj := p.adjust(p, "BH"), by = genotype]
res_pu[, anova_padj := p.adjust(anova_p, "BH"), by = genotype]
res_pu[, U_padj := p.adjust(U_p, "BH"), by = genotype]

key <- unique(tab[, .(species, gene, locus, family, our_family, in_network, deg_pct)])
res_sc <- merge(res_sc, key[species == "sugarcane"], by = c("species", "gene"), all.x = TRUE)
res_pu <- merge(res_pu, key[species == "purple"],    by = c("species", "gene"), all.x = TRUE)

fwrite(res_sc, file.path(outdir, "module20_nitrogen_sugarcane.tsv"), sep = "\t")
fwrite(res_pu, file.path(outdir, "module20_nitrogen_purple.tsv"),    sep = "\t")

# ---------------------------------------------------------------------------
# summaries
# ---------------------------------------------------------------------------
cat("\n=== sugarcane: High vs Low N, within genotype ===\n")
print(res_sc[, .(genes = uniqueN(gene), expressed_TPM_gt1 = sum(mean_tpm > 1),
                 p_lt_0.05 = sum(p < 0.05, na.rm = TRUE),
                 padj_lt_0.05 = sum(padj < 0.05, na.rm = TRUE),
                 median_abs_log2FC = round(median(abs(log2FC_HighminusLow), na.rm = TRUE), 2)),
             by = genotype], row.names = FALSE)
cat("\nper locus (the independent unit) — genes with padj<0.05:\n")
print(res_sc[padj < 0.05, .N, by = .(locus, genotype)][order(locus)], row.names = FALSE)
cat(sprintf("\ngenotype x N interaction padj<0.05: %d of %d genes (with leaf-segment covariate: %d)\n",
            sum(p.adjust(unique(res_sc[, .(gene, interaction_p)])$interaction_p, "BH") < 0.05, na.rm = TRUE),
            uniqueN(res_sc$gene),
            sum(p.adjust(unique(res_sc[, .(gene, interaction_p_segment_adj)])$interaction_p_segment_adj, "BH") < 0.05, na.rm = TRUE)))

cat("\n=== purple: 0N / 2N / 6N, within genotype ===\n")
print(res_pu[, .(genes = uniqueN(gene), expressed_TPM_gt1 = sum(mean_tpm > 1),
                 anova_p_lt_0.05 = sum(anova_p < 0.05, na.rm = TRUE),
                 anova_padj_lt_0.05 = sum(anova_padj < 0.05, na.rm = TRUE),
                 U_padj_lt_0.05 = sum(U_padj < 0.05, na.rm = TRUE),
                 lin_p_lt_0.05 = sum(lin_p < 0.05, na.rm = TRUE)),
             by = genotype], row.names = FALSE)

cat(sprintf("\nOutputs in %s\n", outdir))
