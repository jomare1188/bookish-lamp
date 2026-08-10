# ============================================================================
# GO ENRICHMENT ANALYSIS — ALL CONTRASTS x DIRECTION
# Uses topGO (classic Fisher + BH correction) with annotations from
# eggNOG-mapper v2.1.3.
#
# ID resolution chain (applied once at annotation parse time):
#   eggNOG query IDs  (XP_018757021.1)
#     → GTF CDS attributes  →  FVEG_* locus tags   [STEP 1c]
#     → geneUniverse now in FVEG_* space
#   DEG gene_id column (XM_018905938.1 or FVEG_09964)
#     → tx2gene lookup       →  FVEG_* locus tags   [STEP 1b, applied in loop]
#   Both sides now match → topGO works correctly.
#
# Fixes vs original:
#   1. run_topgo: signif(4) instead of round(digits=4) for p.adj so that
#      very small values like 1.4e-07 are not zeroed out
#   2. run_topgo: parse_pval() handles "< 1e-30" strings from topGO
#   3. save_go_plot: floor p.adj = 0 with a safe fallback (handles the case
#      where ALL values are 0 after rounding)
#   4. save_go_plot: deduplicate Term labels before building factor levels
#   5. save_go_plot: ASCII hyphen instead of em dash in plot title
#
# Required packages:
#   BiocManager::install("topGO")
#   install.packages(c("ggplot2", "dplyr", "readr"))
# ============================================================================

library(topGO)
library(ggplot2)
library(dplyr)
library(readr)

# ============================================================================
# PARAMETERS — edit this section
# ============================================================================

results_dir     <- "/dados04/jorge/comparative_saccharum/files/purple/new/enrichment"
annotation_file <- "/dados04/jorge/comparative_saccharum/annotation/sugarca/emapper.annotations"
genes_in_conserved_edges <- "/dados04/jorge/comparative_saccharum/files/network_conservation/conserved_genes_purple_FULL.txt"

# tx2gene TSV (columns: transcript_id, gene_id, gene_name).
# Used to translate DEG file IDs (XM_* or FVEG_*) → FVEG_* locus tags.
# Set to NULL if DEG files already contain FVEG_* IDs.
#tx2gene_file <- "/dados02/jorge/israel_rnaseq/rnaseq/run01/star_salmon/salmon.merged.tx2gene.tsv"

# GTF file for F. verticillioides.
# Used to build XP_* protein_id → FVEG_* gene_id lookup so that
# eggNOG annotations (queried with XP_* IDs) end up in FVEG_* space.
# Set to NULL if eggNOG was run directly against FVEG_* sequences.
#gtf_file <- "/dados01/jorge/rnaseq_diatraea/reference_genomes/fusarium_verticillioides/GCF_000149555.1_ASM14955v1_genomic.gtf.gz"

# GO ontology to test: "BP", "MF", or "CC"
ontology <- "BP"

# FDR threshold (BH-corrected p-value)
fdr_threshold <- 0.05

# Number of top terms to show in the dot plot
ntop <- 20

# Plot dimensions (cm)
plot_width  <- 40
plot_height <- 30
plot_dpi    <- 300

# Contrasts to process — must match subdirectory names exactly
contrasts <- c(
  "conserverd_sugarcane"
)

# ============================================================================
# STEP 1b — BUILD tx2gene LOOKUP (optional)
# Translates DEG file IDs → FVEG_* locus tags at analysis time (in the loop).
# ============================================================================

cat("=================================================================\n")
cat("GO ENRICHMENT ANALYSIS PIPELINE\n")
cat("=================================================================\n\n")

if (!is.null(tx2gene_file)) {
  cat("[1b] Loading tx2gene mapping...\n")

  tx2gene_raw <- read.table(
    tx2gene_file, sep = "\t", header = TRUE,
    stringsAsFactors = FALSE, quote = ""
  )

  # Build a lookup that handles both input types transparently:
  #   XM_018905938.1  → FVEG_16698   (transcript → gene)
  #   FVEG_16698      → FVEG_16698   (identity, already correct)
  tx2gene_map <- setNames(
    c(tx2gene_raw$gene_id, tx2gene_raw$gene_id),
    c(tx2gene_raw$transcript_id, tx2gene_raw$gene_id)
  )

  cat(sprintf("  Transcripts: %d  |  Unique gene IDs: %d\n",
              nrow(tx2gene_raw), length(unique(tx2gene_raw$gene_id))))
} else {
  tx2gene_map <- NULL
  cat("[1b] No tx2gene file — DEG IDs used as-is.\n")
}

# ============================================================================
# STEP 1c — BUILD PROTEIN → LOCUS TAG MAP FROM GTF (optional)
# Translates eggNOG query IDs (XP_*) → FVEG_* so geneUniverse matches DEGs.
# ============================================================================

if (!is.null(gtf_file)) {
  cat("[1c] Building XP_* -> FVEG_* map from GTF CDS attributes...\n")

  gtf_raw <- read.table(
    gtf_file, sep = "\t", header = FALSE,
    comment.char = "#", quote = "",
    stringsAsFactors = FALSE, fill = TRUE
  )

  # Keep only CDS lines — these are the only ones with protein_id
  gtf_cds <- gtf_raw[gtf_raw[[3]] == "CDS", , drop = FALSE]

  if (nrow(gtf_cds) == 0) {
    stop("No CDS lines found in GTF — check that column 3 contains 'CDS'.")
  }

  attr_col <- gtf_cds[[9]]   # the attributes string column

  # Vectorised extractor: returns NA for lines where the key is absent
  extract_attr <- function(attrs, key) {
    pattern <- paste0(key, ' "([^"]+)"')
    ifelse(
      grepl(pattern, attrs, perl = TRUE),
      sub(paste0('.*', key, ' "([^"]+)".*'), "\\1", attrs, perl = TRUE),
      NA_character_
    )
  }

  protein_ids <- extract_attr(attr_col, "protein_id")
  gene_ids    <- extract_attr(attr_col, "gene_id")

  # Keep only rows where both attributes were found
  keep <- !is.na(protein_ids) & !is.na(gene_ids) &
          protein_ids != "" & gene_ids != ""

  protein_ids <- protein_ids[keep]
  gene_ids    <- gene_ids[keep]

  # Deduplicate (same XP_* can appear in every exon CDS line)
  dup_mask    <- duplicated(protein_ids)
  protein2gene <- setNames(gene_ids[!dup_mask], protein_ids[!dup_mask])

  cat(sprintf("  Unique XP_* protein IDs: %d  |  Unique FVEG_* gene IDs: %d\n",
              length(protein2gene), length(unique(protein2gene))))
} else {
  protein2gene <- NULL
  cat("[1c] No GTF file — eggNOG IDs used as-is in geneUniverse.\n")
}

# ============================================================================
# STEP 1 — PARSE eggNOG-MAPPER ANNOTATION FILE
# ============================================================================

cat("[1]  Loading and parsing eggNOG-mapper annotations...\n")

raw <- read.table(
  annotation_file,
  sep           = "\t",
  header        = FALSE,
  comment.char  = "",
  quote         = "",
  fill          = TRUE,
  stringsAsFactors = FALSE
)

raw <- raw[!grepl("^##",     raw[[1]]), ]
raw <- raw[!grepl("^#query", raw[[1]]), ]

gene_col <- raw[[1]]
go_col   <- raw[[10]]

# Strip isoform suffixes that eggNOG sometimes appends (e.g. XP_018757021.1)
# Note: the version number (.1) is part of the accession and must be kept.
# Only strip biological suffix patterns like .t1 / -T1 appended by some pipelines.
gene_col_clean <- sub("[-.]t[0-9]+$", "", gene_col, ignore.case = TRUE)

# Translate XP_* protein accessions → FVEG_* locus tags using GTF map.
# This makes geneUniverse use the same ID space as the DEG files.
if (!is.null(protein2gene)) {
  translated <- protein2gene[gene_col_clean]
  n_hit      <- sum(!is.na(translated))
  n_miss     <- sum(is.na(translated))
  gene_col_clean <- ifelse(!is.na(translated), translated, gene_col_clean)
  cat(sprintf("  eggNOG IDs translated XP_* -> FVEG_*: %d hit  |  %d kept as-is\n",
              n_hit, n_miss))
  # Sanity check: show a few examples of each
  if (n_hit > 0) {
    ex_idx <- which(!is.na(translated))[1:min(3, n_hit)]
    cat("  Examples (XP_* -> FVEG_*):\n")
    for (i in ex_idx) {
      cat(sprintf("    %s  ->  %s\n", raw[[1]][i], gene_col_clean[i]))
    }
  }
}

go_list <- strsplit(go_col, ",", fixed = TRUE)
go_list <- lapply(go_list, function(x) {
  x <- trimws(x)
  x[x != "-" & x != "" & grepl("^GO:", x)]
})

has_go         <- sapply(go_list, length) > 0
gene2GO        <- go_list[has_go]
names(gene2GO) <- gene_col_clean[has_go]

gene2GO_merged <- tapply(
  seq_along(gene2GO),
  names(gene2GO),
  function(idx) unique(unlist(gene2GO[idx]))
)
gene2GO      <- as.list(gene2GO_merged)
geneUniverse <- names(gene2GO)

cat(sprintf("  Genes with GO annotations: %d\n", length(geneUniverse)))
cat(sprintf("  Sample geneUniverse IDs: %s\n",
            paste(head(geneUniverse, 4), collapse = "  ")))

# ============================================================================
# STEP 2 — HELPER FUNCTIONS
# ============================================================================

#' Parse topGO p-value strings.
#' topGO caps very small p-values and returns them as e.g. "< 1e-30".
#' as.numeric() converts those to NA and they get silently dropped.
#' This function strips the "< " prefix so 1e-30 is kept as a numeric value.
parse_pval <- function(x) {
  num <- suppressWarnings(as.numeric(x))
  is_na <- is.na(num)
  if (any(is_na)) {
    extracted <- suppressWarnings(as.numeric(gsub("^<\\s*", "", x[is_na])))
    num[is_na] <- ifelse(is.na(extracted), NA, extracted)
  }
  num
}

#' Run topGO and return a result data.frame or NULL if nothing significant.
run_topgo <- function(interesting_genes, gene2GO, geneUniverse, ontology) {

  if (length(interesting_genes) == 0) {
    message("    No genes in this set — skipping.")
    return(NULL)
  }

  geneList <- factor(as.integer(geneUniverse %in% interesting_genes))
  names(geneList) <- geneUniverse

  GOdata <- suppressMessages(
    new("topGOdata",
        ontology  = ontology,
        allGenes  = geneList,
        annot     = annFUN.gene2GO,
        gene2GO   = gene2GO)
  )

  allGO <- usedGO(GOdata)
  if (length(allGO) == 0) {
    message("    No GO terms found for these genes — skipping.")
    return(NULL)
  }

  result_classic <- suppressMessages(
    runTest(GOdata, algorithm = "classic", statistic = "fisher")
  )

  table_all <- GenTable(
    GOdata,
    Classic  = result_classic,
    topNodes = length(allGO),
    orderBy  = "Classic"
  )

  # parse_pval() handles "< 1e-30" strings that as.numeric() would turn to NA
  table_all$Classic <- parse_pval(table_all$Classic)
  table_all <- table_all[!is.na(table_all$Classic), ]

  table_sig <- filter(table_all, Classic < 0.05)
  if (nrow(table_sig) == 0) return(NULL)

  # signif() instead of round(): round(1.4e-07, 4) = 0; signif(1.4e-07, 4) = 1.4e-07
  table_sig$p.adj <- signif(p.adjust(table_sig$Classic, method = "BH"), digits = 4)
  table_sig <- table_sig[order(table_sig$p.adj), ]

  table_sig[table_sig$p.adj <= fdr_threshold, ]
}

#' Save dot plot of top GO terms.
save_go_plot <- function(results_df, ntop, direction, contrast_name, out_dir) {

  ggdata <- results_df[seq_len(min(ntop, nrow(results_df))), ]
  ggdata <- ggdata[complete.cases(ggdata), ]
  if (nrow(ggdata) == 0) return(invisible(NULL))

  ggdata$p.adj       <- as.numeric(ggdata$p.adj)
  ggdata$Significant <- as.integer(ggdata$Significant)
  ggdata <- ggdata[order(ggdata$p.adj), ]

  # Floor p.adj == 0 so -log10() stays finite.
  # After signif() this rarely triggers, kept as safety net.
  n_zeroes <- sum(ggdata$p.adj == 0, na.rm = TRUE)
  floor_val <- NA
  if (n_zeroes > 0) {
    nonzero_vals <- ggdata$p.adj[ggdata$p.adj > 0]
    floor_val    <- if (length(nonzero_vals) > 0) min(nonzero_vals) / 2 else fdr_threshold / 1000
    message(sprintf("    NOTE: %d term(s) with p.adj = 0; floored to %.2e for plotting only",
                    n_zeroes, floor_val))
    ggdata$p.adj[ggdata$p.adj == 0] <- floor_val
  }

  # Deduplicate Term labels (two GO IDs can share the same label)
  ggdata <- ggdata[!duplicated(ggdata$Term), ]
  ggdata$Term <- factor(ggdata$Term, levels = rev(unique(ggdata$Term)))

  direction_label <- if (direction == "upregulated") "Up-regulated" else "Down-regulated"
  dot_colour      <- if (direction == "upregulated") "black" else "black"

  # ASCII hyphen — em dash triggers mbcsToSbcs warnings on some systems
  plot_title <- sprintf("GO %s - %s\n%s genes", ontology, contrast_name, direction_label)

  xlab_note <- if (!is.na(floor_val) && n_zeroes > 0)
    sprintf("GO Term  [*%d term(s) with p.adj=0 floored at %.2e for display]",
            n_zeroes, floor_val)
  else
    "GO Term"

  p <- ggplot(ggdata, aes(x = Term, y = -log10(p.adj), size = Significant)) +
    geom_point(colour = dot_colour) +
    scale_size(range = c(2.5, 12.5)) +
    xlab(xlab_note) +
    ylab(expression(-log[10](p[adj]))) +
    labs(title = plot_title, size = "Significant\ngenes") +
    theme_bw(base_size = 14) +
    theme(
      plot.title  = element_text(size = 13, face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 10)
    ) +
    coord_flip()

  base_name <- file.path(out_dir, paste0("GO_", ontology, "_", direction))
  ggsave(paste0(base_name, ".png"), plot = p, device = "png",
         width = plot_width, height = plot_height, dpi = plot_dpi, units = "cm")
  ggsave(paste0(base_name, ".pdf"), plot = p, device = "pdf",
         width = plot_width, height = plot_height, units = "cm")
  ggsave(paste0(base_name, ".svg"), plot = p, device = "svg",
         width = plot_width, height = plot_height, units = "cm")

  invisible(p)
}

# ============================================================================
# STEP 3 — LOOP OVER CONTRASTS AND DIRECTIONS
# ============================================================================

cat("\n[2/3] Running GO enrichment for each contrast...\n")

summary_rows <- list()

for (contrast in contrasts) {

  contrast_dir <- file.path(results_dir, contrast)
  go_out_dir   <- file.path(contrast_dir, "GO_enrichment")
  dir.create(go_out_dir, showWarnings = FALSE, recursive = TRUE)

  cat(sprintf("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"))
  cat(sprintf("Contrast: %s\n", contrast))

  for (direction in c("upregulated", "downregulated")) {

    cat(sprintf("  Direction: %s\n", direction))

    deg_file <- file.path(
      contrast_dir,
      paste0("DEGs_", direction, "_", contrast, ".csv")
    )

    if (!file.exists(deg_file)) {
      cat(sprintf("    WARNING: file not found — %s\n", deg_file))
      next
    }

    deg_df <- tryCatch(
      read_csv(deg_file, show_col_types = FALSE),
      error = function(e) { cat("    ERROR reading file:", e$message, "\n"); NULL }
    )
    if (is.null(deg_df) || !"gene_id" %in% colnames(deg_df)) {
      cat("    WARNING: gene_id column missing — skipping.\n")
      next
    }

    # Translate DEG IDs → FVEG_* via tx2gene (XM_* or already FVEG_* both work)
    raw_ids <- as.character(deg_df$gene_id)

    if (!is.null(tx2gene_map)) {
      translated        <- tx2gene_map[raw_ids]
      n_hit             <- sum(!is.na(translated))
      n_miss            <- sum(is.na(translated))
      interesting_genes <- ifelse(!is.na(translated), translated, raw_ids)
      cat(sprintf("    tx2gene: %d mapped  |  %d kept as-is\n", n_hit, n_miss))
    } else {
      interesting_genes <- raw_ids
    }

    n_in_universe <- sum(interesting_genes %in% geneUniverse)
    cat(sprintf("    DEGs: %d  |  with GO annotation: %d\n",
                length(interesting_genes), n_in_universe))

    if (n_in_universe == 0) {
      cat("    No annotated DEGs — skipping.\n")
      # If this still fires after the fix, uncomment for quick debugging:
      # cat("    Sample DEG IDs:      ", paste(head(interesting_genes, 4), collapse="  "), "\n")
      # cat("    Sample universe IDs: ", paste(head(geneUniverse,     4), collapse="  "), "\n")
      next
    }

    res <- tryCatch(
      run_topgo(interesting_genes, gene2GO, geneUniverse, ontology),
      error = function(e) { cat("    topGO ERROR:", e$message, "\n"); NULL }
    )

    if (is.null(res) || nrow(res) == 0) {
      cat("    No significant GO terms after FDR correction.\n")
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        Contrast     = contrast,
        Direction    = direction,
        DEGs         = length(interesting_genes),
        DEGs_w_GO    = n_in_universe,
        Sig_GO_terms = 0
      )
      next
    }

    cat(sprintf("    Significant GO terms (FDR <= %.2f): %d\n",
                fdr_threshold, nrow(res)))

    out_csv <- file.path(go_out_dir,
                         paste0("GO_", ontology, "_", direction, ".csv"))
    write.csv(res, out_csv, row.names = FALSE, quote = FALSE)

    save_go_plot(res, ntop, direction, contrast, go_out_dir)

    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      Contrast     = contrast,
      Direction    = direction,
      DEGs         = length(interesting_genes),
      DEGs_w_GO    = n_in_universe,
      Sig_GO_terms = nrow(res)
    )
  }
}

# ============================================================================
# STEP 4 — SAVE SUMMARY TABLE
# ============================================================================

cat("\n[3/3] Saving summary...\n")

if (length(summary_rows) > 0) {
  summary_df   <- do.call(rbind, summary_rows)
  rownames(summary_df) <- NULL
  summary_file <- file.path(results_dir, paste0("GO_", ontology, "_enrichment_summary.csv"))
  write.csv(summary_df, summary_file, row.names = FALSE, quote = FALSE)
  cat(sprintf("  Summary saved to: %s\n", summary_file))
  cat("\nResults overview:\n")
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
  print(summary_df, row.names = FALSE)
}

cat("\n=================================================================\n")
cat("DONE! Results saved inside each contrast's GO_enrichment/ subfolder\n")
cat("=================================================================\n")
