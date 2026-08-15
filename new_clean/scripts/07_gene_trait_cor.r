# =============================================================================
# 07_gene_trait_cor.r — per-gene expression vs trait, for one study
#
# Writes:
#   gene_trait_correlations_<study>.tsv   gene, trait, pearson, pval, padj
#   selected_genes_<trait>_<study>.tsv    the responsive subset
#
# WHY PER-GENE AND NOT PER-MODULE: a module eigengene (PC1 over all its genes)
# is a deliberate summary, but it dilutes a real nitrogen signal carried by a
# subset of genes inside an otherwise heterogeneous module. This correlates each
# gene's own expression against the trait instead, restricted to genes sitting
# on at least one conserved cross-species edge. The output is therefore, by
# construction, the intersection of "nitrogen-responsive" and "conserved".
#
# WHERE THE EXPRESSION COMES FROM — this is the fix that motivated the rebuild.
# The old script loaded a DESeq2 object directly, and it loaded a DIFFERENT one
# from the one the purple network was built from (china/run2_onlyL vs
# china/run1 + Group1=="L"). Nothing tied the two together, so the mismatch
# survived for months. This reads results/<study>/vst/<study>.f32 — the exact
# matrix both network layers were computed from. The two cannot diverge because
# there is only one.
#
# RUN: through run.sh  ->  ./run.sh trait sugarcane
# =============================================================================

suppressMessages(library(data.table))

source(file.path(dirname(sub("--file=", "",
       grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "lib", "common.R"))

STUDY       <- env_req("CLEAN_STUDY")
VST_PREFIX  <- env_req("CLEAN_VST_PREFIX")
META_FILE   <- env_req("CLEAN_META")
TRAIT_SPEC  <- env_req("CLEAN_TRAITS")
OUT_DIR     <- ensure_dir(env_req("CLEAN_OUT_DIR"))
GENE_FILTER <- env_opt("CLEAN_GENE_FILTER")     # "" = correlate every gene
SELECT_TRAIT <- env_opt("CLEAN_SELECT_TRAIT", "treatment")
R_THR       <- env_num("CLEAN_TRAIT_R_THR", 0.6)
PADJ_THR    <- env_num("CLEAN_TRAIT_PADJ_THR", 0.05)
setDTthreads(as.integer(env_num("CLEAN_CORES", 100)))

banner(paste("gene-trait correlation:", STUDY))

# --- trait encoding ----------------------------------------------------------
# "genotype:A=1,B=0;treatment:High Nitrogen=1,Low Nitrogen=0"
parse_traits <- function(spec) {
  out <- list()
  for (block in strsplit(spec, ";", fixed = TRUE)[[1L]]) {
    block <- trimws(block)
    if (!nzchar(block)) next
    nm   <- sub(":.*$", "", block)
    body <- sub("^[^:]*:", "", block)
    kv   <- strsplit(strsplit(body, ",", fixed = TRUE)[[1L]], "=", fixed = TRUE)
    if (any(lengths(kv) != 2L))
      stop("cannot parse trait spec '", block, "' — expected LEVEL=value pairs",
           call. = FALSE)
    v <- as.numeric(vapply(kv, `[`, "", 2L))
    names(v) <- trimws(vapply(kv, `[`, "", 1L))
    out[[trimws(nm)]] <- v
  }
  out
}
TRAITS <- parse_traits(TRAIT_SPEC)
say("traits: ", paste(names(TRAITS), collapse = ", "))
for (tr in names(TRAITS))
  say("  ", tr, ": ", paste(sprintf("%s=%g", names(TRAITS[[tr]]), TRAITS[[tr]]),
                            collapse = ", "))

# --- expression --------------------------------------------------------------
say("reading VST from ", basename(VST_PREFIX), ".f32")
vst <- read_vst(VST_PREFIX)
say("  ", fmt_n(nrow(vst)), " genes x ", ncol(vst), " samples")

if (nzchar(GENE_FILTER)) {
  if (!file.exists(GENE_FILTER))
    stop("gene filter not found: ", GENE_FILTER,
         "\n  run  ./run.sh conserve  first, or unset CLEAN_GENE_FILTER",
         call. = FALSE)
  say("restricting to ", basename(GENE_FILTER))
  want <- strip_version(readLines(GENE_FILTER, warn = FALSE))
  keep <- intersect(want, rownames(vst))
  say("  ", fmt_n(length(keep)), " of ", fmt_n(length(want)),
      " conserved genes are in the VST matrix")
  if (length(keep) == 0L)
    stop("no conserved gene matched the VST matrix — id conventions differ",
         call. = FALSE)
  vst <- vst[keep, , drop = FALSE]
} else {
  say("no gene filter — correlating every gene")
}
genes <- rownames(vst)
say("genes to correlate: ", fmt_n(length(genes)))

# --- metadata ----------------------------------------------------------------
meta <- fread(META_FILE, header = TRUE)
setnames(meta, tolower(names(meta)))
if (!"sample" %in% names(meta))
  stop(basename(META_FILE), " has no `sample` column", call. = FALSE)

meta_matched <- meta[match(colnames(vst), sample)]
if (any(is.na(meta_matched$sample)))
  stop("these VST samples are missing from ", basename(META_FILE), ": ",
       paste(head(colnames(vst)[is.na(meta_matched$sample)], 5), collapse = ", "),
       call. = FALSE)

for (tr in names(TRAITS)) {
  if (!tr %in% names(meta_matched))
    stop("trait column '", tr, "' is not in ", basename(META_FILE),
         "\n  available: ", paste(names(meta_matched), collapse = ", "),
         call. = FALSE)
}

# --- correlate ---------------------------------------------------------------
# Vectorised: one cor() call per trait over the whole samples x genes matrix,
# then the same analytic t-test the network's linear layer uses.
say("")
vst_t <- t(vst)

cor_results <- rbindlist(lapply(names(TRAITS), function(tr) {
  enc       <- TRAITS[[tr]]
  levels_in <- as.character(meta_matched[[tr]])
  trait_vec <- unname(enc[levels_in])
  ok        <- !is.na(trait_vec)
  n_ok      <- sum(ok)

  if (n_ok < length(trait_vec)) {
    dropped <- unique(levels_in[!ok])
    say("  ", tr, ": dropping ", length(trait_vec) - n_ok,
        " samples with unencoded level(s): ", paste(dropped, collapse = ", "))
  }
  if (n_ok < 3L)
    stop("trait '", tr, "' has only ", n_ok, " usable samples", call. = FALSE)

  df <- n_ok - 2L
  r  <- as.numeric(cor(vst_t[ok, , drop = FALSE], trait_vec[ok]))
  t_stat <- r * sqrt(df / (1 - r^2 + 1e-15))
  say(sprintf("  %-10s n = %2d, df = %2d", tr, n_ok, df))

  data.table(gene = genes, trait = tr, pearson = round(r, 4),
             pval = 2 * pt(-abs(t_stat), df = df))
}))

cor_results[, padj := signif(p.adjust(pval, method = "BH"), 4), by = trait]
setorder(cor_results, trait, padj)
setcolorder(cor_results, c("gene", "trait", "pearson", "pval", "padj"))
write_tsv(cor_results, file.path(OUT_DIR,
          sprintf("gene_trait_correlations_%s.tsv", STUDY)))

for (tr in names(TRAITS))
  say(sprintf("  %-10s %s genes with padj <= %.2f", tr,
              fmt_n(cor_results[trait == tr & padj <= PADJ_THR, .N]), PADJ_THR))

# --- selection ---------------------------------------------------------------
if (!SELECT_TRAIT %in% names(TRAITS))
  stop("CLEAN_SELECT_TRAIT '", SELECT_TRAIT, "' is not one of the encoded ",
       "traits (", paste(names(TRAITS), collapse = ", "), ")", call. = FALSE)

selected <- cor_results[trait == SELECT_TRAIT &
                        abs(pearson) >= R_THR &
                        padj <= PADJ_THR][order(-abs(pearson))]
say("")
say(sprintf("%s-responsive genes (|r| >= %.2f, padj <= %.2f): %s",
            SELECT_TRAIT, R_THR, PADJ_THR, fmt_n(nrow(selected))))
write_tsv(selected, file.path(OUT_DIR,
          sprintf("selected_genes_%s_%s.tsv", SELECT_TRAIT, STUDY)))

say("")
say("These genes are already restricted to the conserved-edge set, so the")
say("selection is the intersection of '", SELECT_TRAIT, "-responsive' and")
say("'cross-species conserved'. For EDGES rather than nodes, filter")
say("conserved_edges_*_FULL.tsv to rows where both genes are in this list.")
say("done: ", STUDY)
