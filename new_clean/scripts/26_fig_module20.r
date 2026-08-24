# =============================================================================
# 26_fig_module20.r — Muñoz's Module 20 in purple, against sugarcane.
#
# Muñoz-Perez et al. build their nitrogen claim on Module 20: 12 transcripts,
# ~75% MYB or MYB-related, central in their own network. Figure 2 showed the
# module's nitrogen response reproduces in their own data. This figure asks the
# comparative question -- what is left of it in the other species -- and the
# answer is not "less of the same". It is a different response shape, in one
# copy, that every monotone test in this project is blind to.
#
#   A  THE MODULE AS MUNOZ REPORT IT, reproduced -- all 33 mapped sugarcane genes
#      across all 48 libraries, split by nitrogen. Small, because it is the
#      premise rather than the result: it is here so the rest of the figure has
#      something to be compared against.
#
#   B  WHAT THE 12 MEMBERS ACTUALLY ARE. They map to only three Arabidopsis
#      anchors, and the mapped genes are haplotype COPIES of those few loci
#      rather than 12 independent genes -- 27 purple genes from 6 members, 33
#      sugarcane genes from 10. The panel also shows the MYB attrition: the
#      AtMYB59 anchor carries 9 copies in sugarcane and 3 in purple, and the
#      other two anchors are not MYB at all under our own domain rules, however
#      the source study labelled them.
#
#   C  THE SAME MODULE IN PURPLE, as a per-gene Z-SCORE and split by GENOTYPE as
#      well as nitrogen. Both changes are deliberate. The z-score is what makes
#      the one responding copy visible at all -- on an absolute scale it is a
#      bright row among dim ones and the SHAPE of its response is invisible.
#      Splitting by genotype is what makes the dominant structure legible: in
#      purple these genes vary far more between the two species than across
#      nitrogen, and an absolute heatmap split by nitrogen alone hides that.
#      Only the responding copy is labelled; the rest would be 26 unreadable ids.
#
#   D  THAT COPY ON ITS OWN. Soffic.09G0001580-9H, called MYB by the source study
#      AND by our own pipeline, is the only AtMYB59-anchor copy in purple that is
#      expressed at all (60.6 TPM against 0.01 and 0.00 for the other two). Its
#      response is a U centred on the CONTROL: highest under nitrogen starvation,
#      lowest at 2 N, high again under excess. Since 2 N is the control and 0 N
#      and 6 N are stresses in opposite directions, this is a gene induced by
#      nitrogen stress in EITHER direction -- which is why the Spearman test used
#      everywhere else in this project cannot see it, and why the U-shape
#      contrast exists as a separate test.
#
# NOTHING IS WRITTEN ON THE FIGURE THAT BELONGS IN THE LEGEND.
#
# RUN: through run.sh  ->  ./run.sh figmodule20   (then ./run.sh legends)
# =============================================================================

suppressMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(scico); library(grid); library(svglite); library(scales)
})

source(file.path(dirname(sub("--file=", "",
       grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "lib", "common.R"))

M20_DIR    <- env_req("CLEAN_M20_DIR")
TPM_SUG    <- env_req("CLEAN_TPM_SUGARCANE")
TPM_PUR    <- env_req("CLEAN_TPM_PURPLE")
META_SUG   <- env_req("CLEAN_META_SUGARCANE")
META_PUR   <- env_req("CLEAN_META_PURPLE")
OUT_PREFIX <- env_req("CLEAN_OUT_PREFIX")
FIG        <- env_opt("CLEAN_FIG_NUM", "7")
FOCUS      <- env_opt("CLEAN_M20_FOCUS", "Soffic.09G0001580-9H")
FOCUS_ANCH <- env_opt("CLEAN_M20_ANCHOR", "AT5G59780")
setDTthreads(as.integer(env_num("CLEAN_CORES", 8)))

banner(sprintf("Figure %s — Munoz Module 20 in purple", FIG))

ANCHOR_LAB <- c(AT5G59780 = "AT5G59780\nAtMYB59",
                AT3G46130 = "AT3G46130\nAtMYB48",
                AT3G55960 = "AT3G55960\nuncharacterised",
                AT4G10770 = "AT4G10770\nOPT7 transporter")
NSTATUS_LEVELS <- c("low N", "control", "high N")
NSTATUS <- c("Low Nitrogen" = "low N", "High Nitrogen" = "high N",
             "0N" = "low N", "2N" = "control", "6N" = "high N")
PAL_NSTATUS <- setNames(rev(scico(3, palette = "managua", begin = 0.08, end = 0.92)),
                        NSTATUS_LEVELS)
PAL_STUDY <- setNames(scico(3, palette = "batlow")[1:2], c("sugarcane", "purple"))

theme_f <- theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(face = "bold", size = 7.2),
        legend.key.size = unit(3.5, "mm"),
        legend.title = element_text(face = "bold", size = 7.5),
        plot.margin = margin(2, 4, 2, 2))

TAB  <- fread(file.path(M20_DIR, "module20_network_table.tsv"))
NPUR <- fread(file.path(M20_DIR, "module20_nitrogen_purple.tsv"))
NSUG <- fread(file.path(M20_DIR, "module20_nitrogen_sugarcane.tsv"))

load_tpm <- function(tpmf, metaf, sp) {
  tpm <- fread(tpmf); setnames(tpm, 1, "gene")
  tpm[, gene := strip_version(gene)]
  m <- fread(metaf); setnames(m, tolower(names(m)))
  samp <- intersect(names(tpm), m$sample)
  m <- m[match(samp, sample)]
  m[, status := factor(NSTATUS[as.character(treatment)], levels = NSTATUS_LEVELS)]
  genes <- TAB[species == sp, unique(gene)]
  X <- as.matrix(tpm[gene %chin% genes, ..samp])
  rownames(X) <- tpm[gene %chin% genes, gene]
  list(X = X, meta = m)
}
PU <- load_tpm(TPM_PUR, META_PUR, "purple")
SU <- load_tpm(TPM_SUG, META_SUG, "sugarcane")

# =============================================================================
# A — the module as Munoz report it, reproduced in sugarcane
# =============================================================================
zrow <- function(M) {
  Z <- t(scale(t(M))); Z[!is.finite(Z)] <- 0; Z
}
XS <- SU$X; ms <- SU$meta
ordS <- order(ms$status, ms$genotype,
              if ("segment" %in% names(ms)) ms$segment else rep("", nrow(ms)), ms$sample)
XS <- XS[, ordS, drop = FALSE]; ms <- ms[ordS]
ZS <- zrow(XS)
ZS <- ZS[order(-rowMeans(ZS[, ms$status == "low N", drop = FALSE])), , drop = FALSE]

dA <- data.table(gene = factor(rep(rownames(ZS), times = ncol(ZS)),
                               levels = rev(rownames(ZS))),
                 sample = factor(rep(colnames(ZS), each = nrow(ZS)),
                                 levels = colnames(ZS)),
                 z = as.vector(ZS))
dA[, status := ms$status[match(as.character(sample), ms$sample)]]
ZLIM <- as.numeric(quantile(abs(dA$z), 0.98, na.rm = TRUE))
say(sprintf("panel A: %d sugarcane genes x %d libraries", nrow(ZS), ncol(ZS)))

PAL_Z <- rev(scico(256, palette = "roma"))
pA <- ggplot(dA, aes(sample, gene, fill = z)) +
  geom_raster() +
  facet_grid(. ~ status, scales = "free_x", space = "free_x") +
  scale_fill_gradientn(colours = PAL_Z, limits = c(-ZLIM, ZLIM), oob = squish,
                       name = "z", guide = guide_colourbar(
                         barheight = unit(1.4, "cm"), barwidth = unit(2.4, "mm"))) +
  labs(x = sprintf("%d sugarcane libraries", ncol(ZS)),
       y = sprintf("%d Module-20 genes", nrow(ZS))) +
  theme_f +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(), panel.spacing = unit(0.8, "mm"),
        axis.title = element_text(size = 6.5))

# =============================================================================
# B — 12 members are copies of three anchors
# =============================================================================
anch <- unique(TAB[, .(species, gene, at_anchor, our_family)])
anch[, is_myb := grepl("MYB", our_family)]
cnt <- anch[, .(genes = .N, myb = sum(is_myb)), by = .(species, at_anchor)]
cnt <- melt(cnt, id.vars = c("species", "at_anchor"),
            measure.vars = c("myb", "genes"), variable.name = "kind", value.name = "n")
# "genes" is the total, so the non-MYB part is total - myb; stack them.
cnt <- dcast(cnt, species + at_anchor ~ kind, value.var = "n")
cnt[, `not MYB` := genes - myb]
cnt <- melt(cnt[, .(species, at_anchor, `called MYB` = myb, `not MYB`)],
            id.vars = c("species", "at_anchor"), variable.name = "call", value.name = "n")
cnt[, species := factor(species, levels = c("sugarcane", "purple"))]
cnt[, anchor := factor(ANCHOR_LAB[at_anchor], levels = unname(ANCHOR_LAB))]
cnt[, call := factor(call, levels = c("called MYB", "not MYB"))]
say("panel B: mapped genes per Arabidopsis anchor")
print(dcast(cnt, anchor ~ species + call, value.var = "n", fill = 0), row.names = FALSE)

# Horizontal: the anchor labels are two lines each and collide on an x axis at
# any legible size.
pB <- ggplot(cnt, aes(n, anchor, fill = call)) +
  geom_col(width = 0.62, colour = "white", linewidth = 0.3) +
  facet_wrap(~ species, nrow = 1) +
  scale_fill_manual(values = c(`called MYB` = unname(PAL_STUDY[1]), `not MYB` = "grey78"),
                    name = NULL) +
  # 0/5/10 only: at the facet boundary the "15" of one panel and the "0" of the
  # next print on top of each other.
  scale_x_continuous(breaks = c(0, 5, 10),
                     expand = expansion(mult = c(0, 0.14))) +
  labs(y = NULL, x = "mapped genes (haplotype copies)") +
  theme_f + theme(axis.text.y = element_text(size = 5.8, lineheight = 0.95))

# =============================================================================
# C — every purple copy, z-scored, split by genotype as well as nitrogen
# =============================================================================
rdB <- unique(TAB[species == "purple", .(gene, at_anchor, our_family)])
rdB[, anchor := factor(ANCHOR_LAB[at_anchor], levels = unname(ANCHOR_LAB))]
X <- PU$X; mp <- PU$meta
# Genotype first, then nitrogen: in purple these genes vary far more between the
# two species than across nitrogen, and splitting by nitrogen alone hides it.
ordc <- order(mp$genotype, mp$status, mp$sample)
X <- X[, ordc, drop = FALSE]; mp <- mp[ordc]
rdB <- rdB[match(rownames(X), gene)]
rdB[, mean_tpm := rowMeans(X)]
setorder(rdB, anchor, -mean_tpm)
X <- X[rdB$gene, , drop = FALSE]
# Per-gene z. An absolute scale makes the one responding copy a bright row among
# dim ones and hides the SHAPE of its response, which is the thing to see.
ZP <- zrow(X)

dC <- data.table(gene = factor(rep(rownames(ZP), times = ncol(ZP)),
                               levels = rev(rdB$gene)),
                 sample = factor(rep(colnames(ZP), each = nrow(ZP)),
                                 levels = colnames(ZP)),
                 z = as.vector(ZP))
dC[, anchor := rdB$anchor[match(as.character(gene), rdB$gene)]]
dC[, status := mp$status[match(as.character(sample), mp$sample)]]
dC[, genotype := mp$genotype[match(as.character(sample), mp$sample)]]
ZLIM_P <- as.numeric(quantile(abs(dC$z), 0.98, na.rm = TRUE))
say(sprintf("panel C: %d purple copies x %d libraries; %d copies below 1 TPM",
            nrow(ZP), ncol(ZP), sum(rdB$mean_tpm < 1)))

pC <- ggplot(dC, aes(sample, gene, fill = z)) +
  geom_raster() +
  facet_grid(anchor ~ genotype + status, scales = "free", space = "free",
             switch = "y") +
  scale_fill_gradientn(colours = PAL_Z, limits = c(-ZLIM_P, ZLIM_P), oob = squish,
                       name = "z", guide = guide_colourbar(
                         barheight = unit(1.6, "cm"), barwidth = unit(2.4, "mm"))) +
  # Only the responding copy is named; the other 26 ids would be unreadable.
  scale_y_discrete(breaks = FOCUS, labels = FOCUS) +
  labs(x = "18 purple libraries", y = NULL) +
  theme_f +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(),
        axis.text.y = element_text(size = 5.6, face = "bold"),
        panel.grid = element_blank(), panel.spacing = unit(0.7, "mm"),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, size = 5.4, lineheight = 0.95),
        strip.text.x = element_text(size = 5.8))

# =============================================================================
# D — that copy on its own
# =============================================================================
vC <- data.table(sample = colnames(PU$X), tpm = as.numeric(PU$X[FOCUS, ]))
vC <- merge(vC, PU$meta[, .(sample, genotype, status)], by = "sample")
uu <- NPUR[gene == FOCUS, .(genotype, U_est, U_p, anova_p, mean_tpm)]
say("panel D: the focus copy"); print(uu, row.names = FALSE)
print(vC[, .(mean_tpm = round(mean(tpm), 1)), by = .(genotype, status)][order(genotype, status)],
      row.names = FALSE)

pD <- ggplot(vC, aes(status, tpm, colour = genotype, group = genotype)) +
  stat_summary(fun = mean, geom = "line", linewidth = 0.6) +
  geom_point(size = 1.5, alpha = 0.85) +
  scale_colour_manual(values = c(`51NG3` = unname(PAL_STUDY[1]),
                                 TAGZ = unname(PAL_STUDY[2])), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.12))) +
  labs(x = NULL, y = sprintf("%s\nTPM", FOCUS)) +
  theme_f + theme(axis.title.y = element_text(size = 6.5, lineheight = 1.1))

# =============================================================================
# compose
# =============================================================================
# A is the premise rather than the result and gets the narrow column; C carries
# 27 rows of heatmap and gets the wide one.
fig <- (pA | pB) / (pC | pD) +
  plot_layout(widths = c(1, 1.05), heights = c(0.75, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12))

W <- 22; H <- 15
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

nA <- function(sp, a) sum(cnt[species == sp & at_anchor == a, n])
mA <- function(sp, a) cnt[species == sp & at_anchor == a & call == "called MYB", n]
fm <- function(gt, col) uu[genotype == gt][[col]]
mn <- function(gt, st) round(vC[genotype == gt & status == st, mean(tpm)], 1)
sil <- rdB[at_anchor == FOCUS_ANCH & gene != FOCUS, mean_tpm]

legend <- paste0(
"Figure ", FIG, ". Muñoz-Perez et al.'s Module 20 in purple, against sugarcane. Figure 2 showed ",
"that the module's nitrogen response reproduces in the study's own data. This figure asks what ",
"survives in the other species, and the answer is not less of the same thing -- it is a different ",
"response shape, in one copy, that every monotone test in this work is blind to.\n",
"\n",
"(A) The module as reported, reproduced: all ", fmt_n(nrow(ZS)),
" mapped sugarcane genes across all ", fmt_n(ncol(ZS)), " libraries as per-gene z-scores, split ",
"by nitrogen. It is small because it is the premise rather than the result -- the rest of the ",
"figure is what happens to this when it is carried into the other species.\n",
"\n",
"(B) What the 12 published members actually are. They are TransDecoder transcripts from a de novo ",
"assembly; mapped into each reference proteome by reciprocal DIAMOND blastp through Arabidopsis ",
"they resolve to only THREE anchors, and the mapped genes are haplotype COPIES of those few loci ",
"rather than 12 independent genes -- ", fmt_n(sum(cnt[species == "sugarcane", n])),
" sugarcane genes from 10 members and ", fmt_n(sum(cnt[species == "purple", n])),
" purple genes from 6. Fill is our own Pfam/PlnTFDB domain call, which agrees with the source ",
"study on their proteins but not on the orthologs we map to: only the AtMYB59 anchor carries a ",
"MYB call, at ", fmt_n(mA("sugarcane", FOCUS_ANCH)), " of ", fmt_n(nA("sugarcane", FOCUS_ANCH)),
" copies in sugarcane and ", fmt_n(mA("purple", FOCUS_ANCH)), " of ",
fmt_n(nA("purple", FOCUS_ANCH)), " in purple, while the AT3G55960 and OPT7 anchors -- which ",
"supply most of the mapped genes in both species -- are not transcription factors at all. The ",
"AtMYB48 anchor maps to nothing that is in either network. So the module's MYB identity is thin ",
"where it survives and absent where it does not.\n",
"\n",
"(C) The same module in purple, as a per-gene Z-SCORE and split by GENOTYPE as well as nitrogen; ",
"rows are grouped by Arabidopsis anchor and ordered by mean expression, and only the responding ",
"copy is labelled because the other 26 ids would be unreadable at this size. Both choices are ",
"deliberate. The Z-SCORE is what makes that copy legible at all: on an absolute scale it is a ",
"bright row among dim ones and the SHAPE of its response cannot be seen. Splitting by GENOTYPE is ",
"what makes the dominant structure visible -- in purple these genes vary far more between the two ",
"species than across nitrogen, and a panel split by nitrogen alone hides that and invites the ",
"reader to attribute the variation to the treatment. Most copies are flat or near-silent: ",
fmt_n(sum(rdB$mean_tpm < 1)), " of the ", fmt_n(nrow(ZP)), " run below 1 TPM, and the AtMYB59 ",
"block is nearly empty -- of its ", fmt_n(nA("purple", FOCUS_ANCH)), " copies two sit at ",
sprintf("%.2f", max(sil)), " and ", sprintf("%.2f", min(sil)), " TPM, i.e. they are not expressed ",
"in leaf at all.\n",
"\n",
"(D) That copy on its own; points are libraries, lines join the means. ", FOCUS,
" is the only AtMYB59-anchor copy expressed in purple (", sprintf("%.1f", fm("51NG3", "mean_tpm")),
" TPM against ", sprintf("%.2f", max(sil)), " and ", sprintf("%.2f", min(sil)),
" for the other two) and is called MYB both by the source study and by our own domain rules. Its ",
"response is a U CENTRED ON THE CONTROL: in 51NG3 it runs at ", mn("51NG3", "low N"),
" TPM under nitrogen starvation, ", mn("51NG3", "control"), " at the 2 N control and ",
mn("51NG3", "high N"), " under excess -- an eleven-fold drop from starvation to control and a ",
"four-fold rise again beyond it (U-shape contrast c(+1,-2,+1) p = ",
sprintf("%.4f", fm("51NG3", "U_p")), ", one-way p = ", sprintf("%.4f", fm("51NG3", "anova_p")),
"). TAGZ shows the same shape more weakly (", mn("TAGZ", "low N"), ", ", mn("TAGZ", "control"),
", ", mn("TAGZ", "high N"), " TPM).\n",
"\n",
"Because 2 N is the CONTROL and 0 N and 6 N are stresses in opposite directions, this is a gene ",
"induced by nitrogen stress in EITHER direction. That is exactly the shape a monotone rank ",
"correlation cannot see -- it is why the module-level analysis reports this gene\'s module as ",
"unresponsive, and why the U-shape contrast exists as a separate test. For contrast, the nine ",
"expressed copies at the same anchor in SUGARCANE are every one monotonically REPRESSED by ",
"nitrogen (log2 fold change -1.8 to -3.4, BH-adjusted p ~ 1e-08 in the responsive genotype ",
"RB975375) -- but sugarcane\'s design has no control level, only Low and High, so what it can ",
"observe is one arm of a curve. The two studies are therefore consistent rather than ",
"contradictory: sugarcane sees the low-nitrogen arm of the same response, and only purple\'s ",
"three-level design shows that the high-nitrogen side turns back up.\n",
"\n",
"THE LIMIT. This rests on ONE expressed copy in one species, at a p-value that does not clear BH ",
"over the ", fmt_n(nrow(NPUR[genotype == "51NG3"])), " purple Module-20 genes tested. It is a ",
"candidate worth following, not a finding.")

legend_file <- paste0(OUT_PREFIX, "_legend.txt")
writeLines(wrap_at(legend), legend_file)
say("wrote ", basename(legend_file))

stats <- rbindlist(list(
  cnt[, .(panel = "B", study = as.character(species),
          quantity = sprintf("%s: %s", at_anchor, call), value = fmt_n(n))],
  vC[, .(panel = "D", study = "purple",
         quantity = sprintf("%s, %s, %s", FOCUS, genotype, status),
         value = sprintf("%.1f TPM", mean(tpm))), by = .(genotype, status)][
           , .(panel, study, quantity, value)],
  uu[, .(panel = "D", study = "purple",
         quantity = sprintf("%s U-contrast, %s", FOCUS, genotype),
         value = sprintf("est %+.2f, p %.4f", U_est, U_p))]))
write_tsv(stats, paste0(OUT_PREFIX, "_stats.tsv"))
say("done")
