# ============================================================================
# SEMANTIC CONVERGENCE OF THE CONSERVED-SUBNETWORK GO ENRICHMENT
#
# Claim to support: in the topologically-conserved gene set, functional
# enrichment in BOTH species independently points to the SAME basic cellular
# processes. So we cannot merge the two species — we must show each species'
# own enrichment landing on the same regions of GO semantic space.
#
# Design:
#   - load each species' FULL significant term set (per-species tables from
#     top_GO_conserved.r), not just the shared intersection.
#   - build ONE common semantic space from the UNION of terms (GOSemSim Wang,
#     graph-based / IC-free -> no OrgDb needed: godata(ont, computeIC=FALSE)),
#     cluster it once into macro-themes, embed with PCoA.
#   - plot BOTH species in that SAME space, faceted by species: if both cover
#     the same clusters, functional enrichment is convergent.
#   - quantify convergence per theme: sugarcane-only / both / purple-only.
#
# Input : results/conservation/enrichment_conserved/
#             sugarcane/GO_<ONT>_conserved_sugarcane.csv
#             purple/GO_<ONT>_conserved_purple.csv
# Output: results/conservation/enrichment_conserved/semantic/
#
# THE INPUTS CARRY THE ANNOTATION WITH THEM. 09_go_enrichment.r moved onto the
# adopted full-InterProScan table and onto the Pearson-only conserved-gene sets on
# 2026-09-13; this stage reads its output, so re-running it is what brings the
# semantic view onto the same annotation. It has no annotation of its own and must
# never grow one -- two GO sources in one comparison is the failure this project
# spent a week removing. Re-run `./run.sh go <ONT>` before this whenever the
# annotation or the conserved-gene set changes.
#
# WHY Wang AND NOT AN IC-BASED MEASURE: `r_clusterprofiler` has no OrgDb for either
# species, so information content cannot be estimated from a corpus. Wang similarity
# is graph-based and IC-free -- godata(ont, computeIC = FALSE) -- which is what makes
# this runnable at all here. Clustering is base hclust and the embedding is PCoA,
# for the same reason: rrvgo/treemap/pheatmap are not installed in that env.
#
# RUN: through run.sh  ->  ./run.sh gosem        (env: r_clusterprofiler)
# ============================================================================

suppressMessages({
  library(GOSemSim)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(readr)
})

# ============================================================================
# PARAMETERS
# ============================================================================

source(file.path(dirname(sub("--file=", "",
       grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "lib", "common.R"))

in_dir  <- file.path(env_req("CLEAN_RESULTS"), "conservation", "enrichment_conserved")
out_dir <- file.path(in_dir, "semantic")

ontologies       <- env_list("CLEAN_ONTOLOGIES", c("BP", "MF", "CC"))
sim_measure      <- "Wang"                # graph-based, IC-free (no OrgDb needed)
n_macro_clusters <- 12                    # macro-themes on the UNION (visualisation choice)
label_min_size   <- 4                     # label a theme only if it has >= this many union terms

sp_a <- list(label = "sugarcane", file_tag = "sugarcane", colour = "#1b7837")
sp_b <- list(label = "purple",    file_tag = "purple",    colour = "#762a83")

plot_width  <- 32                         # cm (faceted scatter)
plot_height <- 18
plot_dpi    <- 300

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# HELPERS
# ============================================================================

# read a per-species enrichment table -> GO.ID, Term, score(-log10 raw p)
# Raw weight01 p, matching how 09_go_enrichment.r selects terms.
read_species <- function(ont, tag) {
  f <- file.path(in_dir, tag, sprintf("GO_%s_conserved_%s.csv", ont, tag))
  if (!file.exists(f)) return(NULL)
  d <- read_csv(f, show_col_types = FALSE)
  d <- d[!is.na(d$GO.ID), c("GO.ID", "Term", "pvalue")]
  d
}

score_from_padj <- function(p) {
  floor_p <- suppressWarnings(min(p[p > 0], na.rm = TRUE)) / 10
  if (!is.finite(floor_p)) floor_p <- 1e-6
  p[is.na(p) | p == 0] <- floor_p
  -log10(p)
}

analyse_ontology <- function(ont) {

  cat(sprintf("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\nOntology: %s\n", ont))

  a <- read_species(ont, sp_a$file_tag)
  b <- read_species(ont, sp_b$file_tag)
  if (is.null(a) || is.null(b)) { cat("  per-species table(s) missing — skipping.\n"); return(NULL) }
  a$score_a <- score_from_padj(a$pvalue)
  b$score_b <- score_from_padj(b$pvalue)
  cat(sprintf("  Significant %s terms: %s=%d  %s=%d\n",
              ont, sp_a$label, nrow(a), sp_b$label, nrow(b)))

  # --- UNION of terms + per-species membership -----------------------------
  uni <- full_join(a[c("GO.ID","Term","score_a")], b[c("GO.ID","Term","score_b")],
                   by = "GO.ID")
  uni$Term   <- ifelse(is.na(uni$Term.x), uni$Term.y, uni$Term.x)
  uni$in_a   <- !is.na(uni$score_a)
  uni$in_b   <- !is.na(uni$score_b)
  uni$shared <- uni$in_a & uni$in_b
  uni <- uni[, c("GO.ID","Term","score_a","score_b","in_a","in_b","shared")]
  if (nrow(uni) < 3) { cat("  <3 union terms — skipping.\n"); return(NULL) }

  # --- common semantic space (Wang) on the UNION ---------------------------
  semData <- godata(ont = ont, computeIC = FALSE)
  sim <- mgoSim(uni$GO.ID, uni$GO.ID, semData = semData,
                measure = sim_measure, combine = NULL)
  keep <- !is.na(diag(sim))
  if (sum(!keep) > 0) cat(sprintf("  Dropped %d term(s) absent from current GO.db.\n", sum(!keep)))
  sim <- sim[keep, keep, drop = FALSE]
  uni <- uni[match(rownames(sim), uni$GO.ID), ]

  # cluster the union ONCE (shared assignment for both species)
  d  <- as.dist(1 - sim)
  hc <- hclust(d, method = "ward.D2")
  k  <- min(n_macro_clusters, nrow(uni) - 1)
  uni$cluster <- cutree(hc, k = k)
  mds <- cmdscale(d, k = 2); uni$Dim1 <- mds[, 1]; uni$Dim2 <- mds[, 2]

  # representative term per cluster = strongest combined evidence
  uni$score_comb <- pmax(uni$score_a, uni$score_b, na.rm = TRUE)
  reps <- uni %>% group_by(cluster) %>% slice_max(score_comb, n = 1, with_ties = FALSE) %>%
    ungroup() %>% select(cluster, rep_GO = GO.ID, rep_Term = Term)
  uni <- left_join(uni, reps, by = "cluster")

  csize <- table(uni$cluster)
  big   <- as.integer(names(csize[csize >= label_min_size]))

  # ------------------------------------------------------------------------
  # FIGURE 1 — each species in the SAME semantic space (faceted)
  # ------------------------------------------------------------------------
  long <- bind_rows(
    uni[uni$in_a, ] %>% mutate(Species = sprintf("%s  (n=%d)", sp_a$label, nrow(a)),
                               score = score_a),
    uni[uni$in_b, ] %>% mutate(Species = sprintf("%s  (n=%d)", sp_b$label, nrow(b)),
                               score = score_b)
  )
  # labels replicated in both facets so the themes line up visually
  lab <- long
  lab$label <- ifelse(lab$GO.ID %in% reps$rep_GO[reps$cluster %in% big], lab$rep_Term, NA)

  p_sc <- ggplot(long, aes(Dim1, Dim2, colour = factor(cluster), size = score)) +
    geom_point(alpha = 0.8) +
    geom_point(data = long[long$shared, ], shape = 1, colour = "grey15",
               stroke = 0.4, show.legend = FALSE) +
    geom_text_repel(data = lab, aes(label = label), size = 3.1, colour = "grey15",
                    max.overlaps = 18, box.padding = 0.5, min.segment.length = 0,
                    na.rm = TRUE, seed = 1) +
    facet_wrap(~ Species) +
    scale_size(range = c(1.8, 8), name = expression(-log[10](p[adj]))) +
    guides(colour = "none") +
    labs(title = sprintf("Conserved %s enrichment converges on the same themes in both species", ont),
         subtitle = sprintf("union of %d enriched terms in one semantic space; %d ringed terms enriched in BOTH; %d macro-themes",
                            nrow(uni), sum(uni$shared), k),
         x = "PCoA 1", y = "PCoA 2") +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold", size = 13),
          strip.text = element_text(face = "bold"))

  f_sc <- file.path(out_dir, sprintf("GO_%s_conserved_semantic_space_by_species", ont))
  ggsave(paste0(f_sc, ".png"), p_sc, width = plot_width, height = plot_height,
         dpi = plot_dpi, units = "cm")
  ggsave(paste0(f_sc, ".pdf"), p_sc, width = plot_width, height = plot_height, units = "cm")

  # ------------------------------------------------------------------------
  # FIGURE 2 — per-theme convergence: sugarcane-only / both / purple-only
  # ------------------------------------------------------------------------
  conv <- uni %>%
    mutate(status = ifelse(shared, "both",
                    ifelse(in_a, sprintf("%s only", sp_a$label),
                                 sprintf("%s only", sp_b$label)))) %>%
    group_by(rep_Term, status) %>% summarise(n = n(), .groups = "drop")
  ord <- uni %>% group_by(rep_Term) %>% summarise(tot = n(), .groups = "drop") %>%
    arrange(desc(tot)) %>% slice_head(n = 15)
  conv <- conv[conv$rep_Term %in% ord$rep_Term, ]
  conv$rep_Term <- factor(conv$rep_Term, levels = rev(ord$rep_Term))
  conv$status <- factor(conv$status,
                        levels = c(sprintf("%s only", sp_a$label), "both",
                                   sprintf("%s only", sp_b$label)))

  p_bar <- ggplot(conv, aes(rep_Term, n, fill = status)) +
    geom_col() + coord_flip() +
    scale_fill_manual(values = setNames(
      c(sp_a$colour, "#f0a202", sp_b$colour),
      c(sprintf("%s only", sp_a$label), "both", sprintf("%s only", sp_b$label))),
      name = NULL) +
    labs(title = sprintf("Per-theme convergence of conserved %s enrichment", ont),
         subtitle = "terms per macro-theme, split by which species enriched them",
         x = NULL, y = "number of GO terms") +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold", size = 13),
          legend.position = "top")

  f_bar <- file.path(out_dir, sprintf("GO_%s_conserved_semantic_convergence", ont))
  ggsave(paste0(f_bar, ".png"), p_bar, width = 26, height = 18, dpi = plot_dpi, units = "cm")
  ggsave(paste0(f_bar, ".pdf"), p_bar, width = 26, height = 18, units = "cm")

  # --- table ---------------------------------------------------------------
  out_tab <- uni %>%
    mutate(pvalue_sugarcane = 10^(-score_a), pvalue_purple = 10^(-score_b)) %>%
    select(GO.ID, Term, cluster, rep_Term, in_sugarcane = in_a, in_purple = in_b,
           shared, pvalue_sugarcane, pvalue_purple, Dim1, Dim2) %>%
    arrange(cluster, desc(shared))
  write.csv(out_tab,
            file.path(out_dir, sprintf("GO_%s_conserved_semantic_clusters.csv", ont)),
            row.names = FALSE)

  cat(sprintf("  Union: %d terms  |  shared (both species): %d (%.0f%%)  |  themes: %d\n",
              nrow(uni), sum(uni$shared), 100*mean(uni$shared), k))
  cat(sprintf("  -> figures + table in %s\n", out_dir))

  data.frame(Ontology = ont, Sugarcane = nrow(a), Purple = nrow(b),
             Union = nrow(uni), Shared = sum(uni$shared),
             Pct_shared = round(100*mean(uni$shared), 1), Themes = k)
}

# ============================================================================
# MAIN
# ============================================================================

cat("=================================================================\n")
cat("SEMANTIC CONVERGENCE — conserved GO enrichment, both species\n")
cat("=================================================================\n")

summ <- lapply(ontologies, analyse_ontology)
summ <- do.call(rbind, summ[!sapply(summ, is.null)])
if (!is.null(summ)) {
  write.csv(summ, file.path(out_dir, "GO_conserved_semantic_summary.csv"), row.names = FALSE)
  cat("\n"); print(summ, row.names = FALSE)
}

cat("\n=================================================================\n")
cat("DONE.\n")
cat("=================================================================\n")
