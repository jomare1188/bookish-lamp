# ============================================================================
# 08_myb61_expression_test.r — does Kiet et al.'s MYB61 signal exist in the data?
#
# WHY THIS EXISTS
# Our network-level test (07) found 0 MYB61 copies "N-correlated". That test is
# a PEARSON CORRELATION against treatment coded 0/2/6, pooled over both
# genotypes (scripts/gene_trait_cor.r:55). Kiet et al. describe ScMYB61.1 as
#   (a) NON-MONOTONIC  — highest at BOTH low and high N, lowest at normal, and
#   (b) OPPOSITE between the two genotypes.
# A linear correlation has ~zero power against a U-shape, and pooling genotypes
# that move in opposite directions cancels the signal. So the earlier null
# result is partly a blind spot of the test, not evidence against their claim.
#
# This script runs the test their claim actually predicts, on the same leaf
# samples our purple network was built from (18 = 2 genotypes x 3 N x 3 reps):
#   * per-genotype one-way ANOVA across 0N / 2N / 6N
#   * U-SHAPE contrast  c(+1, -2, +1) on (0N, 2N, 6N)  <- Kiet's pattern
#   * LINEAR contrast   c(-1,  0, +1)                  <- what 07 could see
#   * genotype x N interaction (their "opposite patterns" claim)
#   * the pooled Pearson r, to show explicitly why 07 saw nothing
# and draws a heatmap of every copy across all 18 leaf samples.
#
# Power is limited: n=3 per genotype x level, so ANOVA is F(2,6). Interpret
# accordingly — this can support or fail to support, not settle.
#
# RUN: /home/genomics/miniconda3/envs/r_env/bin/Rscript 08_myb61_expression_test.r
# ============================================================================

suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })

base   <- "/dados04/jorge/comparative_saccharum"
outdir <- "/dados04/jorge/comparative_saccharum/new_clean/results/readouts/myb61"
# work/ is cached sequence output from steps 01-06 and lives in the original
# tree, not under results/ -- steps 01-06 are network-independent and are not
# re-run here.
cache  <- file.path(base, "GET_TFS/new/results/myb61")
tpmf   <- file.path(base, "china/run2_onlyL/salmon/salmon.merged.gene_tpm.tsv")
tabf   <- file.path(outdir, "MYB61_network_table.tsv")
purple_pep <- file.path(cache, "work/purple.clean.faa")

KIET_LOCUS <- "09G0002230"   # the published locus, kept as a negative control

# ---------------------------------------------------------------------------
# expression matrix -> long, with genotype / N level parsed from sample names
# ---------------------------------------------------------------------------
tpm <- fread(tpmf)
setnames(tpm, 1, "gene")
samp <- setdiff(names(tpm), c("gene", "gene_name"))

meta <- data.table(sample = samp)
meta[, genotype := fifelse(grepl("51N03", sample), "51NG3", "TAGZ")]
meta[, nlev := fifelse(grepl("_ON-", sample), "0N",
               fifelse(grepl("_2N-", sample), "2N", "6N"))]
meta[, nnum := c("0N" = 0, "2N" = 2, "6N" = 6)[nlev]]
meta[, nlev := factor(nlev, levels = c("0N", "2N", "6N"))]
stopifnot(nrow(meta) == 18, all(table(meta$genotype, meta$nlev) == 3))

# ---------------------------------------------------------------------------
# gene sets: confirmed MYB61 copies + the published locus as negative control
# ---------------------------------------------------------------------------
tab   <- fread(tabf)
myb61 <- tab[species == "purple", .(gene, clade, in_network, degree)]

kiet <- grep(KIET_LOCUS, sub("^>", "", grep("^>", readLines(purple_pep), value = TRUE)),
             value = TRUE)
genes <- rbind(myb61[, .(gene, clade, in_network, degree, set = "MYB61 (this work)")],
               data.table(gene = kiet, clade = "published locus", in_network = NA,
                          degree = NA, set = "Kiet id 09G0002230"))

E <- melt(tpm[gene %in% genes$gene, c("gene", samp), with = FALSE],
          id.vars = "gene", variable.name = "sample", value.name = "tpm")
E[, sample := as.character(sample)]
E <- merge(E, meta, by = "sample")
E <- merge(E, genes, by = "gene")
E[, l2 := log2(tpm + 1)]

cat(sprintf("Genes: %d MYB61 copies + %d copies at the published locus\n",
            uniqueN(genes[set != "Kiet id 09G0002230", gene]), length(kiet)))

# ---------------------------------------------------------------------------
# TESTS
# ---------------------------------------------------------------------------
# contrasts on the ordered levels (0N, 2N, 6N)
C_U   <- c(1, -2, 1)    # U-shape: low + high vs normal   (Kiet's description)
C_LIN <- c(-1, 0, 1)    # monotonic low -> high

contrast_test <- function(y, lev, cvec) {
  m  <- tapply(y, lev, mean); n <- tapply(y, lev, length)
  s2 <- sum(tapply(y, lev, function(v) sum((v - mean(v))^2))) / (length(y) - nlevels(lev))
  if (!is.finite(s2) || s2 <= 0) return(list(est = NA_real_, p = NA_real_))
  est <- sum(cvec * m)
  se  <- sqrt(s2 * sum(cvec^2 / n))
  tt  <- est / se
  list(est = est, p = 2 * pt(-abs(tt), df = length(y) - nlevels(lev)))
}

res <- rbindlist(lapply(unique(E$gene), function(g) {
  rbindlist(lapply(c("51NG3", "TAGZ"), function(gt) {
    d <- E[gene == g & genotype == gt]
    d[, nlev := droplevels(nlev)]
    if (uniqueN(d$nlev) < 3 || all(d$tpm == 0)) {
      return(data.table(gene = g, genotype = gt, mean_tpm = mean(d$tpm),
                        anova_p = NA_real_, U_est = NA_real_, U_p = NA_real_,
                        lin_est = NA_real_, lin_p = NA_real_))
    }
    a  <- tryCatch(anova(lm(l2 ~ nlev, data = d))$`Pr(>F)`[1], error = function(e) NA_real_)
    u  <- contrast_test(d$l2, d$nlev, C_U)
    l  <- contrast_test(d$l2, d$nlev, C_LIN)
    data.table(gene = g, genotype = gt, mean_tpm = mean(d$tpm),
               anova_p = a, U_est = u$est, U_p = u$p, lin_est = l$est, lin_p = l$p)
  }))
}))

# genotype x N interaction ("opposite patterns between genotypes")
inter <- rbindlist(lapply(unique(E$gene), function(g) {
  d <- E[gene == g]
  if (all(d$tpm == 0)) return(data.table(gene = g, interaction_p = NA_real_))
  p <- tryCatch(anova(lm(l2 ~ genotype * nlev, data = d))$`Pr(>F)`[3],
                error = function(e) NA_real_)
  data.table(gene = g, interaction_p = p)
}))

# the pooled linear Pearson — i.e. what step 07 / gene_trait_cor.r measured
pooled <- E[, .(pearson_pooled = suppressWarnings(cor(l2, nnum)),
                pearson_p = tryCatch(cor.test(l2, nnum)$p.value, error = function(e) NA_real_)),
            by = gene]

res <- merge(res, inter, by = "gene")
res <- merge(res, pooled, by = "gene")
res <- merge(res, genes, by = "gene")

# BH within each test family, over the MYB61 copies only (the a-priori set)
for (col in c("anova_p", "U_p", "lin_p")) {
  res[set == "MYB61 (this work)", paste0(sub("_p$", "", col), "_padj") :=
        p.adjust(get(col), method = "BH"), by = genotype]
}
res[set == "MYB61 (this work)", interaction_padj := p.adjust(interaction_p, "BH"), by = genotype]

setcolorder(res, c("set", "gene", "clade", "in_network", "degree", "genotype",
                   "mean_tpm", "anova_p", "anova_padj", "U_est", "U_p", "U_padj",
                   "lin_est", "lin_p", "lin_padj",
                   "interaction_p", "interaction_padj", "pearson_pooled", "pearson_p"))
setorder(res, set, clade, gene, genotype)
fwrite(res, file.path(outdir, "MYB61_expression_anova.tsv"), sep = "\t")

# ---------------------------------------------------------------------------
# summary to console
# ---------------------------------------------------------------------------
m <- res[set == "MYB61 (this work)"]
cat("\n=== per-genotype ANOVA across 0N / 2N / 6N (MYB61 copies) ===\n")
print(m[, .(n_tested = sum(!is.na(anova_p)),
            expressed_mean_tpm_gt1 = sum(mean_tpm > 1, na.rm = TRUE),
            anova_p_lt_0.05 = sum(anova_p < 0.05, na.rm = TRUE),
            anova_padj_lt_0.05 = sum(anova_padj < 0.05, na.rm = TRUE),
            U_p_lt_0.05 = sum(U_p < 0.05, na.rm = TRUE),
            U_padj_lt_0.05 = sum(U_padj < 0.05, na.rm = TRUE),
            lin_p_lt_0.05 = sum(lin_p < 0.05, na.rm = TRUE)),
         by = genotype], row.names = FALSE)

cat("\n=== copies with ANY nominal evidence (p<0.05, either genotype) ===\n")
hit <- m[anova_p < 0.05 | U_p < 0.05, unique(gene)]
if (length(hit)) {
  print(m[gene %in% hit, .(gene, clade, genotype, mean_tpm = round(mean_tpm, 1),
                           anova_p = signif(anova_p, 3), U_est = round(U_est, 2),
                           U_p = signif(U_p, 3), U_padj = signif(U_padj, 3),
                           pooled_r = round(pearson_pooled, 3))], row.names = FALSE)
} else cat("  none\n")

cat("\n=== genotype x N interaction (Kiet's 'opposite patterns') ===\n")
print(unique(m[, .(gene, clade, interaction_p = signif(interaction_p, 3),
                   interaction_padj = signif(interaction_padj, 3))])[
        order(interaction_p)][1:5], row.names = FALSE)

cat("\n=== the published locus, as negative control ===\n")
print(res[set == "Kiet id 09G0002230",
          .(gene, genotype, mean_tpm = round(mean_tpm, 2))], row.names = FALSE)

cat(sprintf("\npooled Pearson |r| over MYB61 copies: max=%.3f (this is what step 07 measured)\n",
            max(abs(m$pearson_pooled), na.rm = TRUE)))

# ---------------------------------------------------------------------------
# HEATMAP
# ---------------------------------------------------------------------------
E[, sample_lab := paste0(genotype, "\n", nlev, "-", sub(".*-", "", sample))]
ord <- meta[order(genotype, nlev, sample), sample]
E[, sample := factor(sample, levels = ord)]

# gene order: clade block, then by mean expression
gord <- E[, .(mu = mean(l2)), by = .(set, clade, gene)][order(set, clade, -mu), gene]
E[, gene := factor(gene, levels = rev(gord))]

# row z-score for the pattern panel
E[, z := (l2 - mean(l2)) / fifelse(sd(l2) > 0, sd(l2), NA_real_), by = gene]

lab_strip <- function(p) p +
  facet_grid(set + clade ~ genotype, scales = "free", space = "free", switch = "y") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7),
        panel.grid = element_blank(),
        strip.text.y.left = element_text(angle = 0, size = 7.5),
        strip.text.x = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        legend.position = "bottom", legend.key.width = unit(1.4, "cm"))

## A — absolute expression: sequential, one hue
pA <- lab_strip(
  ggplot(E, aes(sample, gene, fill = l2)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    scale_fill_gradient(low = "#f7f4f9", high = "#4d004b", name = "log2(TPM + 1)") +
    labs(title = "A  MYB61 copies — absolute expression, Kiet leaf samples",
         subtitle = "18 leaf libraries: 2 genotypes x (0N, 2N, 6N) x 3 replicates",
         x = NULL, y = NULL))

## B — per-gene z-score: diverging, two hues + neutral midpoint
pB <- lab_strip(
  ggplot(E[!is.na(z)], aes(sample, gene, fill = z)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    scale_fill_gradient2(low = "#2166ac", mid = "#f0f0f0", high = "#b2182b",
                         midpoint = 0, name = "row z-score") +
    labs(title = "B  Same data, scaled per copy — shape of the N response",
         subtitle = "Kiet's claim predicts high at 0N and 6N, low at 2N, and opposite between genotypes",
         x = NULL, y = NULL))

fig <- pA / pB + plot_layout(heights = c(1, 1))
ggsave(file.path(outdir, "MYB61_expression_heatmap.png"), fig,
       width = 30, height = 34, units = "cm", dpi = 300)
ggsave(file.path(outdir, "MYB61_expression_heatmap.pdf"), fig,
       width = 30, height = 34, units = "cm")

cat(sprintf("\nOutputs:\n  %s\n  %s\n",
            file.path(outdir, "MYB61_expression_anova.tsv"),
            file.path(outdir, "MYB61_expression_heatmap.png")))
