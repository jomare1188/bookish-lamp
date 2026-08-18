# =============================================================================
# 15_module_profile.r — one row per module: response, coherence, TF content
#
# Joins the module-trait result (linear + non-linear, from 12_gene_trait_mi.py
# run on the eigengene matrix) with module size, eigengene coherence, and a
# transcription-factor enrichment test, into a single table per study.
#
# THE QUESTION THIS EXISTS TO ANSWER: are the modules that respond to nitrogen --
# and in particular the ones that respond NON-LINEARLY, which a correlation
# cannot see -- the TF-rich ones? A regulatory module should be.
#
# TF ENRICHMENT is a hypergeometric test per module against the network's own
# node universe:
#
#     phyper(q - 1, m, n - m, k, lower.tail = FALSE)
#       q = TF genes in the module      m = TF genes in the network
#       k = genes in the module         n = genes in the network
#
# BH across the tested modules. The universe is the NETWORK, not the genome, for
# the same reason the GO background is: it isolates "this module is TF-rich" from
# the generic "co-expressed genes differ from the genome" effect.
#
# One detail borrowed from 11_readouts/get_tfs/05_tf_characterization.r: a gene
# can appear several times in TF_in_network.tsv when different isoforms hit
# different families. Collapse to one row per gene FIRST, or a gene with three
# hits is counted three times and every enrichment is inflated.
#
# RUN: through run.sh  ->  ./run.sh moduleprofile sugarcane
# =============================================================================

suppressMessages(library(data.table))

source(file.path(dirname(sub("--file=", "",
       grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "lib", "common.R"))

STUDY      <- env_req("CLEAN_STUDY")
TRAIT_FILE <- env_req("CLEAN_MODULE_TRAIT")
MEMBERSHIP <- env_req("CLEAN_MEMBERSHIP")
PC1_FILE   <- env_req("CLEAN_PC1_VARIANCE")
TF_FILE    <- env_req("CLEAN_TF_FILE")
NODES      <- env_req("CLEAN_NODE_METRICS")
OUT_FILE   <- env_req("CLEAN_OUT_FILE")
ALPHA      <- env_num("CLEAN_PADJ_THR", 0.05)
# Effect-size floor on the LINEAR side, matching the gene-level TRAIT_R_THR.
# Without it a module with |r| = 0.36 -- 13% of the eigengene's variance -- counts
# as "responsive" at n = 48, and the module count is not comparable to the
# gene-level one that did apply this floor.
R_THR      <- env_num("CLEAN_MODULE_R_THR", 0.6)
setDTthreads(as.integer(env_num("CLEAN_CORES", 100)))

banner(paste("module profile:", STUDY))

for (f in c(TRAIT_FILE, MEMBERSHIP, PC1_FILE, NODES))
  if (!file.exists(f)) stop("missing ", basename(f), call. = FALSE)

# --- module response ---------------------------------------------------------
# 12_gene_trait_mi.py names its first column `gene`; here those are module names.
tr <- fread(TRAIT_FILE)
setnames(tr, "gene", "module")
say("modules tested: ", fmt_n(nrow(tr)))
say("  response: ", paste(sprintf("%s %s", tr[, .N, by = finding]$finding,
                                  fmt_n(tr[, .N, by = finding]$N)), collapse = " | "))

pc1 <- fread(PC1_FILE)
prof <- merge(tr, pc1, by = "module", all.x = TRUE)

# --- TF content --------------------------------------------------------------
mem <- fread(MEMBERSHIP, select = c("gene", "module_name"))
mem[, gene := strip_version(gene)]
setnames(mem, "module_name", "module")
mem <- mem[module %chin% prof$module]

universe <- unique(strip_version(fread(NODES, select = "gene")$gene))
n_univ <- length(universe)

if (file.exists(TF_FILE)) {
  tf <- fread(TF_FILE)
  tf[, gene := strip_version(gene)]
  # collapse isoform-driven multi-family calls to ONE row per gene
  tf <- tf[gene %chin% universe,
           .(families = paste(sort(unique(Family)), collapse = ";")), by = gene]
  m_tf <- nrow(tf)
  say("TF genes in the network: ", fmt_n(m_tf), " of ", fmt_n(n_univ),
      sprintf(" (%.2f%%)", 100 * m_tf / n_univ))

  mem[, is_tf := gene %chin% tf$gene]
  agg <- mem[, .(n_genes_in_module = .N, n_tf = sum(is_tf)), by = module]
  agg[, tf_frac := n_tf / n_genes_in_module]
  # hypergeometric: P(X >= q)
  agg[, tf_p := phyper(n_tf - 1L, m_tf, n_univ - m_tf, n_genes_in_module,
                       lower.tail = FALSE)]
  agg[, tf_padj := p.adjust(tf_p, method = "BH")]

  fam <- merge(mem[is_tf == TRUE, .(gene, module)], tf, by = "gene")
  top_fam <- fam[, .(top_tf_families = paste(head(names(sort(table(
                        unlist(strsplit(families, ";", fixed = TRUE))),
                        decreasing = TRUE)), 3), collapse = ";")), by = module]
  agg <- merge(agg, top_fam, by = "module", all.x = TRUE)
  prof <- merge(prof, agg, by = "module", all.x = TRUE)
} else {
  say("NOTE: ", basename(TF_FILE), " not found — skipping TF enrichment")
  prof[, `:=`(n_genes_in_module = NA_integer_, n_tf = NA_integer_,
              tf_frac = NA_real_, tf_p = NA_real_, tf_padj = NA_real_,
              top_tf_families = NA_character_)]
}

# --- effect-size floor, on BOTH sides ---------------------------------------
# Flooring only the linear side is not a neutral choice: it moves every module
# with a modest linear response out of `both`/`pearson_only` and into `mi_only`,
# which then reads as "non-linear" when it is nothing of the kind. Measured here:
# flooring Pearson alone took mi_only from 253 to 786 in sugarcane.
#
# So MI needs an equivalent floor. The Gaussian identity the network layer uses
# to convert nats to |r| does NOT apply -- it assumes two continuous variables,
# and here the trait is discrete with MI capped at H(trait). Instead the floor is
# CALIBRATED EMPIRICALLY on this data: among modules whose linear response sits
# at |r| ~ R_THR, what MI do they carry? That median is the MI equivalent of the
# linear cut, at this n, this class structure, and this estimator.
#
# It is a calibration, not a theoretical equivalence, and it is reported so the
# reader can see the number rather than trust the label.
null_file <- sub("\\.tsv$", ".null.tsv", TRAIT_FILE)
H <- NA_real_
if (file.exists(null_file)) {
  nl <- fread(null_file, header = TRUE, fill = TRUE)
  hit <- nl[[1]] == "trait_entropy"
  if (any(hit, na.rm = TRUE)) H <- as.numeric(nl[[2]][which(hit)[1]])
}
prof[, mi_norm := if (is.finite(H) && H > 0) mi / H else NA_real_]

near <- prof[abs(abs(pearson) - R_THR) < 0.03 & is.finite(mi)]
if (nrow(near) >= 30) {
  MI_THR <- median(near$mi)
  say(sprintf("MI floor calibrated on %s modules at |r| ~ %.2f: MI >= %.3f (%.2f of H)",
              fmt_n(nrow(near)), R_THR, MI_THR, MI_THR / H))
} else {
  MI_THR <- 0
  say(sprintf("WARNING: only %d modules near |r| = %.2f -- cannot calibrate an MI",
              nrow(near), R_THR))
  say("         floor, so the MI side keeps no effect-size cut and `mi_only` will")
  say("         include modules whose linear response merely fell below the floor.")
}

prof[, finding_stat := finding]                 # the unfloored, padj-only call
pass_pearson <- prof$pearson_padj <= ALPHA & abs(prof$pearson) >= R_THR
pass_mi      <- prof$padj <= ALPHA & prof$mi >= MI_THR
prof[, finding := fifelse(pass_pearson & pass_mi, "both",
                  fifelse(pass_mi, "mi_only",
                  fifelse(pass_pearson, "pearson_only", "neither")))]

say(sprintf("effect-size floors: |r| >= %.2f (linear), MI >= %.3f (non-linear)",
            R_THR, MI_THR))
say("  padj only:    ", paste(sprintf("%s %s", prof[, .N, by = finding_stat]$finding_stat,
                                fmt_n(prof[, .N, by = finding_stat]$N)), collapse = " | "))
say("  with floors:  ", paste(sprintf("%s %s", prof[, .N, by = finding]$finding,
                                fmt_n(prof[, .N, by = finding]$N)), collapse = " | "))

prof[, responsive := finding != "neither"]
prof[, tf_enriched := !is.na(tf_padj) & tf_padj <= ALPHA]
setorder(prof, padj, pearson_padj)
setcolorder(prof, c("module", "n_genes", "pc1_var_pct", "finding", "finding_stat",
                    "pearson", "pearson_padj", "mi", "mi_norm", "padj",
                    "n_tf", "tf_frac", "tf_p", "tf_padj", "top_tf_families"))
write_tsv(prof, OUT_FILE)

# --- the headline cross-tab --------------------------------------------------
banner("response x TF enrichment")
ct <- prof[, .(modules = .N, tf_enriched = sum(tf_enriched, na.rm = TRUE)), by = finding]
ct[, pct_tf_enriched := round(100 * tf_enriched / modules, 2)]
setorder(ct, -modules)
print(ct, row.names = FALSE)

base_rate <- prof[finding == "neither", mean(tf_enriched, na.rm = TRUE)]
say("")
for (f in c("mi_only", "pearson_only", "both")) {
  r <- prof[finding == f]
  if (!nrow(r)) next
  obs <- mean(r$tf_enriched, na.rm = TRUE)
  # Fisher against the non-responsive modules -- the honest comparison, since
  # every module here passed the same size and coherence filters
  tab <- matrix(c(sum(r$tf_enriched, na.rm = TRUE), nrow(r) - sum(r$tf_enriched, na.rm = TRUE),
                  prof[finding == "neither", sum(tf_enriched, na.rm = TRUE)],
                  prof[finding == "neither", .N - sum(tf_enriched, na.rm = TRUE)]), nrow = 2)
  ft <- fisher.test(tab)
  say(sprintf("  %-13s TF-enriched %5.2f%% vs %5.2f%% in non-responsive   OR %.2f  p %.3g",
              f, 100 * obs, 100 * base_rate, ft$estimate, ft$p.value))
}
say("")
say("done: ", STUDY)
