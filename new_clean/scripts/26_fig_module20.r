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
#   A  WHAT THE 12 MEMBERS ACTUALLY ARE. They map to only three Arabidopsis
#      anchors, and the mapped genes are haplotype COPIES of those few loci
#      rather than 12 independent genes -- 27 purple genes from 6 members, 33
#      sugarcane genes from 10. The panel also shows the MYB attrition: the
#      AtMYB59 anchor carries 9 copies in sugarcane and 3 in purple, and the
#      other two anchors are not MYB at all under our own domain rules, however
#      the source study labelled them.
#
#   B  ALL 27 PURPLE COPIES ACROSS ALL 18 LIBRARIES. Most are silent or flat.
#      Drawn as absolute log2(TPM+1) rather than a z-score precisely because a
#      z-score would rescale the silent copies to look as structured as the one
#      that matters.
#
#   C  THE ONE COPY THAT RESPONDS. Soffic.09G0001580-9H, called MYB by the source
#      study AND by our own pipeline, is the only AtMYB59-anchor copy in purple
#      that is expressed at all (60.6 TPM against 0.01 and 0.00 for the other
#      two). Its response is a U centred on the CONTROL: highest under nitrogen
#      starvation, lowest at 2 N, high again under excess. Since 2 N is the
#      control and 0 N and 6 N are stresses in opposite directions, this is a
#      gene induced by nitrogen stress in EITHER direction -- which is why the
#      Spearman test used everywhere else in this project cannot see it, and why
#      the U-shape contrast exists.
#
#   D  THE SAME ANCHOR IN SUGARCANE, for contrast: nine expressed copies, every
#      one of them monotonically repressed by nitrogen. Sugarcane's design has no
#      control level, only Low and High, so what it can see is one arm of a curve.
#      Read C and D together and the two studies are consistent rather than
#      contradictory: sugarcane observes the low-nitrogen arm of the same shape.
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
# A — 12 members are copies of three anchors
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
say("panel A: mapped genes per Arabidopsis anchor")
print(dcast(cnt, anchor ~ species + call, value.var = "n", fill = 0), row.names = FALSE)

# Horizontal: the anchor labels are two lines each and collide on an x axis at
# any legible size.
pA <- ggplot(cnt, aes(n, anchor, fill = call)) +
  geom_col(width = 0.62, colour = "white", linewidth = 0.3) +
  facet_wrap(~ species, nrow = 1) +
  scale_fill_manual(values = c(`called MYB` = unname(PAL_STUDY[1]), `not MYB` = "grey78"),
                    name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.10))) +
  labs(y = NULL, x = "mapped genes (haplotype copies)") +
  theme_f + theme(axis.text.y = element_text(size = 5.8, lineheight = 0.95))

# =============================================================================
# B — every purple copy across every purple library
# =============================================================================
rdB <- unique(TAB[species == "purple", .(gene, at_anchor, our_family)])
rdB[, anchor := factor(ANCHOR_LAB[at_anchor], levels = unname(ANCHOR_LAB))]
X <- PU$X; mp <- PU$meta
ordc <- order(mp$status, mp$genotype, mp$sample)
X <- X[, ordc, drop = FALSE]; mp <- mp[ordc]
rdB <- rdB[match(rownames(X), gene)]
rdB[, mean_tpm := rowMeans(X)]
setorder(rdB, anchor, -mean_tpm)
X <- X[rdB$gene, , drop = FALSE]

dB <- data.table(gene = factor(rep(rownames(X), times = ncol(X)),
                               levels = rev(rdB$gene)),
                 sample = factor(rep(colnames(X), each = nrow(X)),
                                 levels = colnames(X)),
                 l2 = as.vector(log2(X + 1)))
dB[, anchor := rdB$anchor[match(as.character(gene), rdB$gene)]]
dB[, status := mp$status[match(as.character(sample), mp$sample)]]
say(sprintf("panel B: %d purple copies x %d libraries; %d copies below 1 TPM",
            nrow(X), ncol(X), sum(rdB$mean_tpm < 1)))

pB <- ggplot(dB, aes(sample, gene, fill = l2)) +
  geom_raster() +
  facet_grid(anchor ~ status, scales = "free", space = "free", switch = "y") +
  scale_fill_scico(palette = "davos", direction = -1,
                   name = expression(log[2]*"(TPM+1)")) +
  labs(x = "18 purple libraries", y = NULL) +
  theme_f +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(), panel.spacing = unit(0.8, "mm"),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, size = 5.6, lineheight = 0.95),
        legend.key.width = unit(2.6, "mm"), legend.key.height = unit(1.5, "cm"))

# =============================================================================
# C — the one copy that responds, in purple
# =============================================================================
vC <- data.table(sample = colnames(PU$X), tpm = as.numeric(PU$X[FOCUS, ]))
vC <- merge(vC, PU$meta[, .(sample, genotype, status)], by = "sample")
uu <- NPUR[gene == FOCUS, .(genotype, U_est, U_p, anova_p, mean_tpm)]
say("panel C: the focus copy"); print(uu, row.names = FALSE)
print(vC[, .(mean_tpm = round(mean(tpm), 1)), by = .(genotype, status)][order(genotype, status)],
      row.names = FALSE)

pC <- ggplot(vC, aes(status, tpm, colour = genotype, group = genotype)) +
  stat_summary(fun = mean, geom = "line", linewidth = 0.6) +
  geom_point(size = 1.5, alpha = 0.85) +
  scale_colour_manual(values = c(`51NG3` = unname(PAL_STUDY[1]),
                                 TAGZ = unname(PAL_STUDY[2])), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.12))) +
  labs(x = NULL, y = sprintf("%s\nTPM", FOCUS)) +
  theme_f + theme(axis.title.y = element_text(size = 6.5, lineheight = 1.1))

# =============================================================================
# D — the same anchor in sugarcane
# =============================================================================
gD <- unique(TAB[species == "sugarcane" & at_anchor == FOCUS_ANCH, gene])
gD <- intersect(gD, rownames(SU$X))
dD <- rbindlist(lapply(gD, function(g)
  data.table(gene = g, sample = colnames(SU$X), tpm = as.numeric(SU$X[g, ]))))
dD <- merge(dD, SU$meta[, .(sample, status)], by = "sample")
dD <- dD[, .(tpm = mean(tpm)), by = .(gene, status)]
say(sprintf("panel D: %d sugarcane copies at the %s anchor", length(gD), FOCUS_ANCH))
print(dcast(dD, gene ~ status, value.var = "tpm"), row.names = FALSE)

pD <- ggplot(dD, aes(status, tpm, group = gene)) +
  geom_line(linewidth = 0.5, colour = "grey55") +
  geom_point(aes(colour = status), size = 1.6) +
  scale_colour_manual(values = PAL_NSTATUS, guide = "none") +
  scale_y_log10(breaks = 10^(0:3),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "l", linewidth = 0.2,
                      short = unit(0.4, "mm"), mid = unit(0.7, "mm"),
                      long = unit(1.1, "mm")) +
  labs(x = NULL, y = sprintf("sugarcane %s copies\nmean TPM", "AtMYB59")) +
  theme_f + theme(axis.title.y = element_text(size = 6.5, lineheight = 1.1))

# =============================================================================
# compose
# =============================================================================
fig <- (pA | pB) / (pC | pD) +
  plot_layout(heights = c(1, 0.8)) +
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
"(A) What the 12 published members actually are. They are TransDecoder transcripts from a de novo ",
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
"(B) Every purple copy across every purple library, as absolute log2(TPM + 1) rather than a ",
"z-score: a z-score would rescale a silent copy to look exactly as structured as an expressed ",
"one, which is the specific error this panel exists to prevent. Rows are grouped by anchor and ",
"ordered by mean expression, columns by nitrogen status. Most copies are low or flat, and the ",
"AtMYB59 block is nearly empty -- of its ", fmt_n(nA("purple", FOCUS_ANCH)), " copies, two run at ",
sprintf("%.2f", max(sil)), " and ", sprintf("%.2f", min(sil)), " TPM on average, i.e. they are ",
"not expressed in leaf at all.\n",
"\n",
"(C) The one that is. ", FOCUS, " is the only AtMYB59-anchor copy expressed in purple (",
sprintf("%.1f", fm("51NG3", "mean_tpm")), " TPM) and is called MYB both by the source study and ",
"by our own domain rules. Points are libraries, lines join the means. Its response is a U CENTRED ",
"ON THE CONTROL: in 51NG3 it runs at ", mn("51NG3", "low N"), " TPM under nitrogen starvation, ",
mn("51NG3", "control"), " at the 2 N control and ", mn("51NG3", "high N"),
" under excess -- an eleven-fold drop from starvation to control and a four-fold rise again ",
"beyond it (U-shape contrast c(+1,-2,+1) p = ", sprintf("%.4f", fm("51NG3", "U_p")),
", one-way p = ", sprintf("%.4f", fm("51NG3", "anova_p")), "; it does not clear BH over the ",
fmt_n(nrow(NPUR[genotype == "51NG3"])), " purple Module-20 genes tested, so it is a candidate ",
"rather than a finding). TAGZ shows the same shape more weakly (", mn("TAGZ", "low N"), ", ",
mn("TAGZ", "control"), ", ", mn("TAGZ", "high N"), " TPM).\n",
"\n",
"Because 2 N is the CONTROL and 0 N and 6 N are stresses in opposite directions, this is a gene ",
"induced by nitrogen stress in EITHER direction. That is exactly the shape a monotone rank ",
"correlation cannot see -- it is the reason the module-level analysis reports this gene's module ",
"as unresponsive, and the reason the U-shape contrast exists as a separate test.\n",
"\n",
"(D) The same anchor in sugarcane, on a log axis: ", fmt_n(length(gD)), " expressed copies, every ",
"one monotonically REPRESSED by nitrogen (log2 fold change -1.8 to -3.4, BH-adjusted p ~ 1e-08 in ",
"the responsive genotype RB975375). Sugarcane's design has no control level -- only Low and High ",
"-- so what it can observe is one arm of a curve. Read C and D together and the two studies are ",
"consistent rather than contradictory: sugarcane sees the low-nitrogen arm of the same response, ",
"and only purple's three-level design reveals that the high-nitrogen side turns back up. The ",
"caveat that limits how far this can be pushed is that it rests on ONE expressed copy in one ",
"species, at a p-value that does not survive correction.")

legend_file <- paste0(OUT_PREFIX, "_legend.txt")
writeLines(wrap_at(legend), legend_file)
say("wrote ", basename(legend_file))

stats <- rbindlist(list(
  cnt[, .(panel = "A", study = as.character(species),
          quantity = sprintf("%s: %s", at_anchor, call), value = fmt_n(n))],
  vC[, .(panel = "C", study = "purple",
         quantity = sprintf("%s, %s, %s", FOCUS, genotype, status),
         value = sprintf("%.1f TPM", mean(tpm))), by = .(genotype, status)][
           , .(panel, study, quantity, value)],
  uu[, .(panel = "C", study = "purple",
         quantity = sprintf("%s U-contrast, %s", FOCUS, genotype),
         value = sprintf("est %+.2f, p %.4f", U_est, U_p))],
  dD[, .(panel = "D", study = "sugarcane",
         quantity = sprintf("%s, %s", gene, status),
         value = sprintf("%.1f TPM", tpm))]))
write_tsv(stats, paste0(OUT_PREFIX, "_stats.tsv"))
say("done")
