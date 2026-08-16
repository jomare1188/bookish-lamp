# =============================================================================
# 14_module_eigengene.r — one eigengene per MCL module
#
# Summarises each module by the first principal component of its member genes'
# expression, so a module's nitrogen response can be tested once instead of once
# per gene. That is the point: the gene-level analysis died on a testing burden
# of 39,226 / 44,118 genes, and there are only ~6,500 modules.
#
# WHAT THE EIGENGENE IS
#
#   prcomp(t(vst_sub), center = TRUE, scale. = TRUE)   -> PCA on the gene-gene
#                                                         CORRELATION matrix
#   pc1 <- pca$x[, 1]                                  -> sample scores
#   orient so cor(pc1, mean module expression) > 0
#   z-score across samples
#
# Scaling genes to unit variance before the PCA is the WGCNA convention and it
# matters here: without it PC1 chases the few highest-variance members instead of
# summarising the module.
#
# THE OUTPUT IS DELIBERATELY THE VST EXPORT'S FORMAT. `<prefix>.f32` +
# `.genes.txt` + `.meta.json`, with module names where gene names normally go.
# That is what lets 12_gene_trait_mi.py run on modules completely unchanged, so
# the module-level linear/non-linear call comes from the SAME validated code path
# as the gene-level one and the two answers cannot diverge for methodological
# reasons.
#
# THREE FIXES relative to the dead scripts/eigengene.r this is lifted from:
#   1. expression comes from read_vst(), i.e. the exact matrix both network
#      layers were built from -- not a separately loaded dds, which is how the
#      old tree let the network and the trait analysis drift apart
#   2. the "too few genes" guard runs AFTER the zero-variance filter. The
#      original checked before (:99 vs :51), so a 2-gene module with one flat
#      gene reached prcomp as a single row
#   3. eigengenes are z-scored. The original kept raw PC1 scores, whose magnitude
#      scales with module size -- Module_001 ran +-120 against +-2 for a 5-gene
#      module, which makes modules incomparable on a plot. Pearson and rank-based
#      MI are both invariant to this, so it costs nothing statistically.
#
# RUN: through run.sh  ->  ./run.sh eigengene sugarcane
# =============================================================================

suppressMessages(library(data.table))

source(file.path(dirname(sub("--file=", "",
       grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "lib", "common.R"))

STUDY      <- env_req("CLEAN_STUDY")
VST_PREFIX <- env_req("CLEAN_VST_PREFIX")
MEMBERSHIP <- env_req("CLEAN_MEMBERSHIP")
OUT_PREFIX <- env_req("CLEAN_OUT_PREFIX")
MIN_SIZE   <- as.integer(env_num("CLEAN_MIN_MODULE_SIZE_EIGEN", 3))
setDTthreads(as.integer(env_num("CLEAN_CORES", 100)))

banner(paste0("module eigengenes: ", STUDY, "  (modules with >= ", MIN_SIZE, " genes)"))

# --- inputs ------------------------------------------------------------------
if (!file.exists(MEMBERSHIP))
  stop("missing ", basename(MEMBERSHIP), "\n  run  ./run.sh mcl ", STUDY,
       "  first", call. = FALSE)

say("reading VST")
vst <- read_vst(VST_PREFIX)
say("  ", fmt_n(nrow(vst)), " genes x ", ncol(vst), " samples")

mem <- fread(MEMBERSHIP, select = c("gene", "module_name"))
mem[, gene := strip_version(gene)]
mem <- mem[module_name != "Unassigned" & gene %chin% rownames(vst)]
sizes <- mem[, .N, by = module_name]
keep <- sizes[N >= MIN_SIZE, module_name]
say("modules: ", fmt_n(nrow(sizes)), " named | ", fmt_n(length(keep)),
    " with >= ", MIN_SIZE, " genes present in the VST (",
    fmt_n(mem[module_name %chin% keep, .N]), " genes)")
if (length(keep) == 0L) stop("no module meets the size cutoff", call. = FALSE)

mem <- mem[module_name %chin% keep]
by_module <- split(mem$gene, mem$module_name)

# --- the eigengene -----------------------------------------------------------
compute_module_eigengene <- function(vst_sub) {
  # zero-variance genes carry no information and break scale. = TRUE
  vst_sub <- vst_sub[apply(vst_sub, 1L, var, na.rm = TRUE) > 0, , drop = FALSE]
  # guard AFTER the filter, not before -- see the header
  if (nrow(vst_sub) < 2L) return(NULL)

  pca <- prcomp(t(vst_sub), center = TRUE, scale. = TRUE)
  pc1 <- pca$x[, 1L]
  var_pct <- summary(pca)$importance[2L, 1L] * 100

  # PC1's sign is arbitrary; orient it so "up" means "more expressed"
  if (cor(pc1, colMeans(vst_sub)) < 0) pc1 <- -pc1

  s <- stats::sd(pc1)
  if (!is.finite(s) || s == 0) return(NULL)
  list(pc1 = (pc1 - mean(pc1)) / s, var_pct = var_pct, n_used = nrow(vst_sub))
}

say("computing ", fmt_n(length(by_module)), " eigengenes")
modules <- sort(names(by_module))
eig <- vector("list", length(modules))
info <- vector("list", length(modules))
skipped <- character(0)
t0 <- Sys.time()

for (i in seq_along(modules)) {
  m <- modules[i]
  res <- compute_module_eigengene(vst[by_module[[m]], , drop = FALSE])
  if (is.null(res)) { skipped <- c(skipped, m); next }
  eig[[i]]  <- res$pc1
  info[[i]] <- data.table(module = m, n_genes = length(by_module[[m]]),
                          n_genes_used = res$n_used,
                          pc1_var_pct = round(res$var_pct, 3))
  if (i %% 1000 == 0)
    say(sprintf("  %s/%s  (%.1f min)", fmt_n(i), fmt_n(length(modules)),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

ok <- !vapply(eig, is.null, logical(1))
if (length(skipped))
  say("skipped ", fmt_n(length(skipped)),
      " modules with < 2 informative genes after the zero-variance filter")
eig <- eig[ok]; modules <- modules[ok]
info <- rbindlist(info[ok])
say("computed ", fmt_n(length(eig)), " eigengenes")

E <- do.call(rbind, eig)                    # modules x samples
rownames(E) <- modules
colnames(E) <- colnames(vst)

say(sprintf("PC1 variance explained: median %.1f%%  mean %.1f%%  (min %.1f, max %.1f)",
            median(info$pc1_var_pct), mean(info$pc1_var_pct),
            min(info$pc1_var_pct), max(info$pc1_var_pct)))

# --- write, in the VST export's own format -----------------------------------
invisible(ensure_dir(dirname(OUT_PREFIX)))

con <- file(paste0(OUT_PREFIX, ".f32"), "wb")
writeBin(as.vector(t(E)), con, size = 4L)   # row-major, as read_vst expects
close(con)
writeLines(rownames(E), paste0(OUT_PREFIX, ".genes.txt"))

json_esc <- function(x) paste0('"', gsub('"', '\\\\"', x), '"')
writeLines(paste0(
  '{\n',
  '  "label": ',      json_esc(paste0(STUDY, "_eigengenes")),                ',\n',
  '  "n_genes": ',    nrow(E),                                              ',\n',
  '  "n_samples": ',  ncol(E),                                              ',\n',
  '  "unit": ',       json_esc("module eigengene (PC1, z-scored)"),         ',\n',
  '  "min_module_size": ', MIN_SIZE,                                        ',\n',
  '  "source_vst": ', json_esc(VST_PREFIX),                                 ',\n',
  '  "created": ',    json_esc(format(Sys.time(), "%Y-%m-%d %H:%M:%S")),     ',\n',
  '  "samples": [',   paste(json_esc(colnames(E)), collapse = ", "),        ']\n',
  '}'), paste0(OUT_PREFIX, ".meta.json"))

write_tsv(info[order(-n_genes)], paste0(OUT_PREFIX, "_pc1_variance.tsv"))
say("wrote ", basename(OUT_PREFIX), ".{f32,genes.txt,meta.json}  (",
    fmt_n(nrow(E)), " modules x ", ncol(E), " samples)")
say("done: ", STUDY)
