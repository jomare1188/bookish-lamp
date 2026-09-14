# ============================================================================
# GO ENRICHMENT — CONSERVED SUB-NETWORK SANITY CHECK (H0)
#
# Purpose: sanity-check the cross-species conservation pipeline. If the
# OrthoFinder-based conservation (R570 <-> LA purple reference proteomes) is
# capturing real biology and not noise, then the genes sitting on conserved
# co-expression edges should be enriched for universally-conserved,
# co-expression-heavy processes (photosynthesis, translation/ribosome,
# oxidative phosphorylation, etc.).
#
# This is a stripped, generalised version of top_GO_dev2.r:
#   - takes a PLAIN gene list (one gene ID per line), no up/down direction
#   - background = GO-annotated NODES of THAT species' network (NOT the whole
#     genome), so enrichment reflects what is special about the conserved genes
#     relative to the network itself, removing the generic "co-expressed genes
#     differ from the genome" effect
#   - runs BOTH networks (sugarcane, purple) in one pass
#   - NO tx2gene / GTF translation: conserved-gene IDs already match the
#     eggNOG query IDs in each species' annotation, verified:
#       sugarcane: SoffiXsponR570.NNXgNNNNNN  ==  emapper col1
#       purple   : Soffic.*                    ==  emapper col1
#
# RUN INSIDE conda env:  topGO_env
#   conda activate topGO_env
#   Rscript top_GO_conserved.r
#
# Required packages (already in topGO_env): topGO, ggplot2, dplyr, readr
# ============================================================================

library(topGO)
library(ggplot2)
library(dplyr)
library(readr)

# ============================================================================
# PARAMETERS — edit this section
# ============================================================================

# Everything here comes from config.sh via run.sh -- the ontology in particular,
# which used to be a hand-edited variable that had to be set three times.
source(file.path(dirname(sub("--file=", "",
       grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "lib", "common.R"))

RESULTS  <- env_req("CLEAN_RESULTS")
ontology <- env_req("CLEAN_ONTOLOGY")
p_threshold   <- env_num("CLEAN_GO_P", 0.05)   # on the RAW weight01 p
ntop          <- as.integer(env_num("CLEAN_GO_NTOP", 20))

plot_width  <- 40          # cm
plot_height <- 30          # cm
plot_dpi    <- 300

cons_dir <- file.path(RESULTS, "conservation")
# WHICH SIDE OF THE PARTITION. `conserved` is the genes on at least one conserved
# edge; `nonconserved` is the exact complement -- network nodes with NO conserved
# edge. Together they partition the node set (sugarcane 37,867 + 64,123 = 101,990),
# so the two runs are a real contrast against one shared background rather than two
# overlapping sets. The literal alternative, "genes on at least one NON-conserved
# edge", is 99.3% of the network -- with a mean degree near 1,477 almost every gene
# has some non-conserved edge -- so it would be tested against itself.
GENE_SET <- env_opt("CLEAN_GENE_SET", "conserved")
if (!GENE_SET %in% c("conserved", "nonconserved"))
  stop("CLEAN_GENE_SET must be conserved or nonconserved (got '", GENE_SET, "')",
       call. = FALSE)
# The tag rides on every output path AND every filename, so the two runs can never
# overwrite each other or be mistaken for one another downstream.
TAG <- GENE_SET
enrich_dir <- file.path(cons_dir, sprintf("enrichment_%s", TAG))

# One entry per network.
#   gene_list = genes on >=1 conserved edge (06_conservation_join.r) [interest set]
#   nodes     = node-metrics col 1 = every gene in the network. Its GO-annotated
#               subset is the UNIVERSE, so the test asks what is special about the
#               conserved genes RELATIVE TO THE NETWORK, not to the genome. This is
#               why node_metrics must be regenerated whenever the network changes.
#
# node_id_strip is gone: ids are normalised once, in 01_export_vst.r, so every
# file under results/ already carries the bare form.
# WHICH CONSERVED-GENE SET, and it is not a cosmetic switch. `_FULL` is the genes on
# a conserved edge of the MERGED (Pearson + MI) graph, written by 06; `_pearson` is
# the same thing on the Pearson-only graphs, written by 62. They are different gene
# sets from different networks -- 39,226/44,118 against 37,867/42,194 -- so mixing
# them with the current annotation would describe two analyses at once. Default is
# the current one; set CLEAN_CONS_SET=FULL to reproduce the merged-network result.
CONS_SET <- env_opt("CLEAN_CONS_SET", "pearson")
if (!CONS_SET %in% c("pearson", "FULL"))
  stop("CLEAN_CONS_SET must be pearson or FULL (got '", CONS_SET, "')", call. = FALSE)
CONS_GENES_FILE <- function(sp) {
  if (GENE_SET == "conserved") {
    sprintf("conserved_genes_%s_%s.txt", sp, CONS_SET)
  } else {
    sprintf("nonconserved_genes_%s_%s.txt", sp, CONS_SET)
  }
}
cat(sprintf("gene set: %s  (%s edge set)\n", GENE_SET, CONS_SET))

networks <- list(
  list(
    label      = "sugarcane",
    annotation = env_req("CLEAN_EMAPPER_SUGARCANE"),
    gene2go    = env_opt("CLEAN_GENE2GO_SUGARCANE", ""),
    gene_list  = file.path(cons_dir, CONS_GENES_FILE("sugarcane")),
    nodes      = file.path(RESULTS, "sugarcane", "network_sugarcane_node_metrics.tsv"),
    node_id_strip = NULL,
    out_dir    = file.path(enrich_dir, "sugarcane")
  ),
  list(
    label      = "purple",
    annotation = env_req("CLEAN_EMAPPER_PURPLE"),
    gene2go    = env_opt("CLEAN_GENE2GO_PURPLE", ""),
    gene_list  = file.path(cons_dir, CONS_GENES_FILE("purple")),
    nodes      = file.path(RESULTS, "purple", "network_purple_node_metrics.tsv"),
    node_id_strip = NULL,
    out_dir    = file.path(enrich_dir, "purple")
  )
)

# ============================================================================
# HELPER FUNCTIONS  (carried over from top_GO_dev2.r, unchanged logic)
# ============================================================================

#' Parse an eggNOG-mapper annotation file into a gene2GO list.
#' Column 1 = query ID, column 10 = comma-separated GO IDs ("-" if none).
parse_eggnog <- function(annotation_file) {
  raw <- read.table(
    annotation_file, sep = "\t", header = FALSE,
    comment.char = "", quote = "", fill = TRUE, stringsAsFactors = FALSE
  )
  raw <- raw[!grepl("^##",     raw[[1]]), ]
  raw <- raw[!grepl("^#query", raw[[1]]), ]

  gene_col <- raw[[1]]
  go_col   <- raw[[10]]

  go_list <- strsplit(go_col, ",", fixed = TRUE)
  go_list <- lapply(go_list, function(x) {
    x <- trimws(x)
    x[x != "-" & x != "" & grepl("^GO:", x)]
  })

  has_go         <- sapply(go_list, length) > 0
  gene2GO        <- go_list[has_go]
  names(gene2GO) <- gene_col[has_go]

  gene2GO_merged <- tapply(
    seq_along(gene2GO), names(gene2GO),
    function(idx) unique(unlist(gene2GO[idx]))
  )
  as.list(gene2GO_merged)
}

#' topGO caps very small p-values as e.g. "< 1e-30"; keep them numeric.
parse_pval <- function(x) {
  num   <- suppressWarnings(as.numeric(x))
  is_na <- is.na(num)
  if (any(is_na)) {
    extracted <- suppressWarnings(as.numeric(gsub("^<\\s*", "", x[is_na])))
    num[is_na] <- ifelse(is.na(extracted), NA, extracted)
  }
  num
}

#' Run topGO (weight01 Fisher + BH) and return a significant-terms data.frame.
#'
#' ALGORITHM: weight01, topGO's default. It walks the GO DAG and down-weights a
#' term's genes when a more specific child term already explains the signal, so
#' the reported terms are the specific ones rather than every ancestor they
#' inherit significance from. `classic` -- which this pipeline used previously --
#' scores each term independently and therefore reports whole ancestor chains as
#' separate "findings", which inflates the count and makes the term list read as
#' far more informative than it is.
#'
#' THRESHOLD: the RAW weight01 p-value, not an FDR-adjusted one. This is topGO's
#' own convention and the reason is structural, not a shortcut. weight01's
#' p-values are deliberately NOT independent -- the algorithm's whole mechanism
#' is to condition each term on its neighbours in the DAG -- so they are not a
#' family of exchangeable tests and BH's assumptions do not hold on them. The
#' conditioning has already absorbed most of the redundancy that a correction
#' would otherwise be compensating for.
#'
#' A BH-adjusted column is still computed and written as `p.adj`, over ALL tested
#' terms, so the stricter convention remains available in the output without a
#' re-run. It is NOT what selects the terms.
run_topgo <- function(interesting_genes, gene2GO, geneUniverse, ontology) {

  if (length(interesting_genes) == 0) {
    message("    No genes in this set — skipping.")
    return(NULL)
  }

  geneList <- factor(as.integer(geneUniverse %in% interesting_genes))
  names(geneList) <- geneUniverse

  GOdata <- suppressMessages(
    new("topGOdata",
        ontology = ontology,
        allGenes = geneList,
        annot    = annFUN.gene2GO,
        gene2GO  = gene2GO)
  )

  allGO <- usedGO(GOdata)
  if (length(allGO) == 0) {
    message("    No GO terms found for these genes — skipping.")
    return(NULL)
  }

  result <- suppressMessages(
    runTest(GOdata, algorithm = "weight01", statistic = "fisher")
  )

  table_all <- GenTable(
    GOdata, pvalue = result,
    topNodes = length(allGO), orderBy = "pvalue"
  )

  table_all$pvalue <- parse_pval(table_all$pvalue)
  table_all <- table_all[!is.na(table_all$pvalue), ]
  if (nrow(table_all) == 0) return(NULL)

  # Reported for reference only -- see the note above on why it does not select.
  # BH runs over EVERY tested term, not a pre-filtered subset: an earlier version
  # filtered to p < 0.05 first and then corrected, which shrinks the denominator
  # to terms already known to be small and makes the adjusted values
  # anti-conservative. BH's m must be the number of tests performed.
  # signif() not round(): round(1.4e-07,4)=0 but signif(1.4e-07,4)=1.4e-07
  table_all$p.adj <- signif(p.adjust(table_all$pvalue, method = "BH"), digits = 4)
  table_all <- table_all[order(table_all$pvalue), ]

  table_all[table_all$pvalue <= p_threshold, ]
}

#' Save a dot plot of the top GO terms.
save_go_plot <- function(results_df, ntop, label, out_dir) {

  ggdata <- results_df[seq_len(min(ntop, nrow(results_df))), ]
  ggdata <- ggdata[complete.cases(ggdata), ]
  if (nrow(ggdata) == 0) return(invisible(NULL))

  ggdata$pvalue      <- as.numeric(ggdata$pvalue)
  ggdata$Significant <- as.integer(ggdata$Significant)
  ggdata <- ggdata[order(ggdata$pvalue), ]

  n_zeroes  <- sum(ggdata$pvalue == 0, na.rm = TRUE)
  floor_val <- NA
  if (n_zeroes > 0) {
    nonzero_vals <- ggdata$pvalue[ggdata$pvalue > 0]
    floor_val    <- if (length(nonzero_vals) > 0) min(nonzero_vals) / 2 else p_threshold / 1000
    message(sprintf("    NOTE: %d term(s) with p = 0; floored to %.2e for plotting only",
                    n_zeroes, floor_val))
    ggdata$pvalue[ggdata$pvalue == 0] <- floor_val
  }

  ggdata <- ggdata[!duplicated(ggdata$Term), ]
  ggdata$Term <- factor(ggdata$Term, levels = rev(unique(ggdata$Term)))

  plot_title <- sprintf("GO %s - %s sub-network\n%s", ontology, TAG, label)
  xlab_note  <- if (!is.na(floor_val) && n_zeroes > 0)
    sprintf("GO Term  [*%d term(s) with p=0 floored at %.2e for display]",
            n_zeroes, floor_val)
  else "GO Term"

  p <- ggplot(ggdata, aes(x = Term, y = -log10(pvalue), size = Significant)) +
    geom_point(colour = "black") +
    scale_size(range = c(2.5, 12.5)) +
    xlab(xlab_note) +
    ylab(expression(-log[10](p))) +
    labs(title = plot_title, size = "Significant\ngenes") +
    theme_bw(base_size = 14) +
    theme(
      plot.title  = element_text(size = 13, face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 10)
    ) +
    coord_flip()

  base_name <- file.path(out_dir, paste0("GO_", ontology, "_", TAG, "_", label))
  ggsave(paste0(base_name, ".png"), plot = p, device = "png",
         width = plot_width, height = plot_height, dpi = plot_dpi, units = "cm")
  ggsave(paste0(base_name, ".pdf"), plot = p, device = "pdf",
         width = plot_width, height = plot_height, units = "cm")
  invisible(p)
}

# ============================================================================
# MAIN LOOP — one enrichment per network
# ============================================================================

cat("=================================================================\n")
cat(sprintf("GO ENRICHMENT — %s sub-network (background = network nodes)\n", TAG))
cat("=================================================================\n\n")

summary_rows   <- list()
results_by_net <- list()   # keep each network's significant-term table for the cross-comparison

for (net in networks) {

  cat(sprintf("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"))
  cat(sprintf("Network: %s\n", net$label))

  # --- annotation (genome-wide gene -> GO) ---------------------------------
  # CLEAN_GENE2GO_* takes precedence: the adopted full-InterProScan table, which
  # reaches 62.9%/64.4% of network genes against the eggNOG GO column's 8.0%/7.2%.
  # Same switch 18_module_go.r:111-119 makes, and the same parse_gene2go() reader,
  # so the module level and the conserved-set level cannot end up on two
  # annotations. See annotation/README.md for why eggNOG's column was rejected.
  if (nzchar(net$gene2go)) {
    cat(sprintf("  [1] Parsing derived GO table: %s\n", basename(net$gene2go)))
    gene2GO_all <- parse_gene2go(net$gene2go)
  } else {
    cat("  [1] Parsing eggNOG annotation...\n")
    gene2GO_all <- parse_eggnog(net$annotation)
  }
  cat(sprintf("      Genome GO-annotated genes: %d\n", length(gene2GO_all)))

  # --- background: GO-annotated NODES of THIS network ----------------------
  # Universe = genes present in the filtered network (node-metrics col 1) that
  # also carry a GO annotation. Isolates the conservation signal from the
  # generic "co-expressed genes differ from the genome" effect.
  cat("  [2] Building network-node background...\n")
  node_ids <- read.table(net$nodes, sep = "\t", header = TRUE,
                         stringsAsFactors = FALSE, quote = "")[[1]]
  if (!is.null(net$node_id_strip)) {
    node_ids <- sub(net$node_id_strip, "", node_ids)
  }
  node_ids     <- unique(trimws(node_ids))
  geneUniverse <- intersect(node_ids, names(gene2GO_all))
  gene2GO      <- gene2GO_all[geneUniverse]      # restrict annotation to universe
  cat(sprintf("      Network nodes: %d  |  GO-annotated (= background): %d\n",
              length(node_ids), length(geneUniverse)))

  # --- interesting set: genes on conserved edges ---------------------------
  cat("  [3] Loading conserved gene list...\n")
  conserved <- readLines(net$gene_list)
  conserved <- unique(trimws(conserved))
  conserved <- conserved[conserved != ""]
  n_in_universe <- sum(conserved %in% geneUniverse)
  cat(sprintf("      Conserved genes: %d  |  in node background (with GO): %d\n",
              length(conserved), n_in_universe))

  if (n_in_universe == 0) {
    cat("      No annotated conserved genes in background — skipping.\n")
    next
  }

  # --- enrichment ----------------------------------------------------------
  cat("  [4] Running topGO (weight01 Fisher + BH over all tested terms)...\n")
  res <- tryCatch(
    run_topgo(conserved, gene2GO, geneUniverse, ontology),
    error = function(e) { cat("      topGO ERROR:", e$message, "\n"); NULL }
  )

  dir.create(net$out_dir, showWarnings = FALSE, recursive = TRUE)

  n_sig <- if (is.null(res)) 0 else nrow(res)
  if (n_sig == 0) {
    cat("      No significant GO terms after FDR correction.\n")
  } else {
    cat(sprintf("      Enriched GO terms (raw weight01 p <= %.2f): %d\n", p_threshold, n_sig))
    out_csv <- file.path(net$out_dir, paste0("GO_", ontology, "_", TAG, "_", net$label, ".csv"))
    write.csv(res, out_csv, row.names = FALSE)   # quote: GO Term labels contain commas
    cat(sprintf("      Table -> %s\n", out_csv))
    save_go_plot(res, ntop, net$label, net$out_dir)
    cat(sprintf("      Plot  -> %s/GO_%s_%s_%s.{png,pdf}\n",
                net$out_dir, ontology, TAG, net$label))
  }

  results_by_net[[net$label]] <- res   # NULL if no significant terms — handled below

  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    Network              = net$label,
    Conserved_genes      = length(conserved),
    Conserved_in_bg      = n_in_universe,
    Node_background_w_GO = length(geneUniverse),
    Sig_GO_terms         = n_sig
  )
}

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n[4] Saving summary...\n")
if (length(summary_rows) > 0) {
  summary_df <- do.call(rbind, summary_rows)
  rownames(summary_df) <- NULL
  summary_file <- file.path(enrich_dir,
    paste0("GO_", ontology, "_", TAG, "_summary.csv"))
  dir.create(dirname(summary_file), showWarnings = FALSE, recursive = TRUE)
  write.csv(summary_df, summary_file, row.names = FALSE, quote = FALSE)
  cat(sprintf("    Summary -> %s\n\n", summary_file))
  print(summary_df, row.names = FALSE)
}

# ============================================================================
# CROSS-SPECIES TERM COMPARISON — how many enriched terms are shared vs unique
# Matching is on GO.ID (accession), NOT the Term label. Produces 4 files at the
# top of enrichment_conserved/, following the existing naming pattern:
#   GO_<ont>_<tag>_shared_terms.csv        (terms enriched in BOTH species)
#   GO_<ont>_<tag>_unique_<A>.csv          (terms enriched only in A)
#   GO_<ont>_<tag>_unique_<B>.csv          (terms enriched only in B)
#   GO_<ont>_<tag>_comparison_summary.csv  (counts + Jaccard)
#   where <tag> is `conserved` or `nonconserved` -- see CLEAN_GENE_SET above.
# Comparison CSVs are quoted (default) because GO Term labels can contain commas.
# ============================================================================

cat("\n[5] Comparing enriched terms across networks...\n")

cmp_dir <- enrich_dir

labA <- networks[[1]]$label            # sugarcane
labB <- networks[[2]]$label            # purple
resA <- results_by_net[[labA]]
resB <- results_by_net[[labB]]

empty_terms <- function(x) is.null(x) || nrow(x) == 0

if (empty_terms(resA) || empty_terms(resB)) {
  cat(sprintf("    One or both networks have no significant %s terms — skipping comparison.\n",
              ontology))
} else {
  idA <- resA$GO.ID
  idB <- resB$GO.ID

  shared_ids <- intersect(idA, idB)
  uniqA_ids  <- setdiff(idA, idB)
  uniqB_ids  <- setdiff(idB, idA)

  # shared: side-by-side stats from both networks, strongest (worst-case p) first
  a <- resA[match(shared_ids, resA$GO.ID), c("GO.ID", "Term", "Significant", "pvalue")]
  b <- resB[match(shared_ids, resB$GO.ID), c("Significant", "pvalue")]
  colnames(a) <- c("GO.ID", "Term", paste0("Significant_", labA), paste0("pvalue_", labA))
  colnames(b) <- c(paste0("Significant_", labB), paste0("pvalue_", labB))
  shared_df <- cbind(a, b)
  shared_df <- shared_df[order(pmax(shared_df[[paste0("pvalue_", labA)]],
                                    shared_df[[paste0("pvalue_", labB)]])), ]

  uniqA_df <- resA[resA$GO.ID %in% uniqA_ids, ]
  uniqB_df <- resB[resB$GO.ID %in% uniqB_ids, ]

  f_shared <- file.path(cmp_dir, sprintf("GO_%s_%s_shared_terms.csv", ontology, TAG))
  f_uniqA  <- file.path(cmp_dir, sprintf("GO_%s_%s_unique_%s.csv", ontology, TAG, labA))
  f_uniqB  <- file.path(cmp_dir, sprintf("GO_%s_%s_unique_%s.csv", ontology, TAG, labB))
  f_cmpsum <- file.path(cmp_dir, sprintf("GO_%s_%s_comparison_summary.csv", ontology, TAG))

  write.csv(shared_df, f_shared, row.names = FALSE)
  write.csv(uniqA_df,  f_uniqA,  row.names = FALSE)
  write.csv(uniqB_df,  f_uniqB,  row.names = FALSE)

  n_union <- length(union(idA, idB))
  jaccard <- length(shared_ids) / n_union

  cmp_summary <- data.frame(Ontology = ontology, check.names = FALSE,
                            stringsAsFactors = FALSE)
  cmp_summary[[paste0("Terms_",  labA)]] <- length(idA)
  cmp_summary[[paste0("Terms_",  labB)]] <- length(idB)
  cmp_summary[["Shared"]]                <- length(shared_ids)
  cmp_summary[[paste0("Unique_", labA)]] <- length(uniqA_ids)
  cmp_summary[[paste0("Unique_", labB)]] <- length(uniqB_ids)
  cmp_summary[["Union"]]                 <- n_union
  cmp_summary[["Jaccard"]]               <- signif(jaccard, 4)
  write.csv(cmp_summary, f_cmpsum, row.names = FALSE)

  cat(sprintf("    Shared %s terms: %d  |  unique %s: %d  |  unique %s: %d  |  Jaccard: %.3f\n",
              ontology, length(shared_ids), labA, length(uniqA_ids),
              labB, length(uniqB_ids), jaccard))
  cat(sprintf("    Shared  -> %s\n", f_shared))
  cat(sprintf("    Unique %s -> %s\n", labA, f_uniqA))
  cat(sprintf("    Unique %s -> %s\n", labB, f_uniqB))
  cat(sprintf("    Summary -> %s\n", f_cmpsum))
  cat("\n"); print(cmp_summary, row.names = FALSE)
}

cat("\n=================================================================\n")
cat("DONE.\n")
cat("=================================================================\n")
