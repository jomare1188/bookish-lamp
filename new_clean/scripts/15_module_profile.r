# =============================================================================
# 15_module_profile.r — one row per module: response, coherence, TF content
#
# Joins the module-trait result (Spearman, from 19_module_trait_spearman.r) with
# module size, eigengene coherence, and a transcription-factor enrichment test,
# into a single table per study.
#
# THE QUESTION THIS EXISTS TO ANSWER: are the modules that respond to nitrogen
# the TF-rich ones? A regulatory module should be.
#
# THE RESPONSE CALL IS NOT MADE HERE. It used to be: this script owned the
# effect-size floors and the pearson_only/mi_only/both classification. With one
# statistic there is one rule -- padj <= MODULE_PADJ_THR and |rho| >= MODULE_R_THR
# -- and it lives in 19 with the statistic it applies to, so there is a single
# source of truth for what "responsive" means. This script reads `responsive`
# and `direction` and does not recompute them.
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
# Alpha for the TF hypergeometric only; the response call arrives pre-made.
ALPHA      <- env_num("CLEAN_PADJ_THR", 0.05)
setDTthreads(as.integer(env_num("CLEAN_CORES", 100)))

banner(paste("module profile:", STUDY))

for (f in c(TRAIT_FILE, MEMBERSHIP, PC1_FILE, NODES))
  if (!file.exists(f)) stop("missing ", basename(f), call. = FALSE)

# --- module response ---------------------------------------------------------
tr <- fread(TRAIT_FILE)
need <- c("module", "rho", "pval", "padj", "responsive", "direction")
if (length(setdiff(need, names(tr))))
  stop(basename(TRAIT_FILE), " is missing column(s): ",
       paste(setdiff(need, names(tr)), collapse = ", "),
       "\n  found: ", paste(names(tr), collapse = ", "),
       "\n  this stage expects the 19_module_trait_spearman.r schema; re-run",
       "\n  ./run.sh moduletrait <study>", call. = FALSE)
say("modules tested: ", fmt_n(nrow(tr)))
say(sprintf("  responsive: %s  (%s positive, %s negative)",
            fmt_n(tr[responsive == TRUE, .N]),
            fmt_n(tr[direction == "positive", .N]),
            fmt_n(tr[direction == "negative", .N])))

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

prof[, tf_enriched := !is.na(tf_padj) & tf_padj <= ALPHA]
prof <- prof[order(padj, -abs(rho))]
setcolorder(prof, c("module", "n_genes", "pc1_var_pct", "responsive", "direction",
                    "rho", "pval", "padj",
                    "n_tf", "tf_frac", "tf_p", "tf_padj", "top_tf_families"))
write_tsv(prof, OUT_FILE)

# --- the headline cross-tab --------------------------------------------------
# One statistic means one comparison: responsive vs not. The direction split is
# reported underneath because "modules that rise with nitrogen" and "modules that
# fall" are different biology, but it is a description, not a second test family.
banner("response x TF enrichment")
ct <- prof[, .(modules = .N, tf_enriched = sum(tf_enriched, na.rm = TRUE)),
           by = .(group = fifelse(responsive, "responsive", "not responsive"))]
ct[, pct_tf_enriched := round(100 * tf_enriched / modules, 2)]
setorder(ct, -modules)
print(ct, row.names = FALSE)

# Fisher against the non-responsive modules -- the honest comparison, since every
# module here passed the same size and coherence filters.
base <- prof[responsive == FALSE]
base_rate <- mean(base$tf_enriched, na.rm = TRUE)
say("")
for (g in list(list("responsive", prof[responsive == TRUE]),
               list("  positive", prof[direction == "positive"]),
               list("  negative", prof[direction == "negative"]))) {
  r <- g[[2]]
  if (!nrow(r)) next
  tab <- matrix(c(sum(r$tf_enriched, na.rm = TRUE), nrow(r) - sum(r$tf_enriched, na.rm = TRUE),
                  sum(base$tf_enriched, na.rm = TRUE), nrow(base) - sum(base$tf_enriched, na.rm = TRUE)),
                nrow = 2)
  ft <- fisher.test(tab)
  say(sprintf("  %-13s %5s modules, TF-enriched %5.2f%% vs %5.2f%% in non-responsive   OR %.2f  p %.3g",
              g[[1]], fmt_n(nrow(r)), 100 * mean(r$tf_enriched, na.rm = TRUE),
              100 * base_rate, ft$estimate, ft$p.value))
}
say("")
say("done: ", STUDY)
