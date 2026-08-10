library(data.table)
library(igraph)
library(leidenAlg)
library(ggplot2)

# ==============================================================================
# CONFIGURATION
# ==============================================================================
RESOLUTION           <- 1.0   # Leiden resolution — tune if modules are too big/small
                               #   too many tiny modules → lower (0.5, 0.3)
                               #   too few huge modules  → raise (1.5, 2.0)
MIN_MODULE_SIZE      <- 2     # modules smaller than this → labelled "Unassigned"
                               # set to 2 so virtually all modules are retained
MIN_MODULE_SIZE_PLOT <- 10    # minimum module size shown in the plot
N_ITERATIONS         <- 10    # leiden n.iterations per run
N_STARTS             <- 5     # independent starts — best modularity is kept
RANDOM_SEED          <- 42

OUT_DIR_SUGARCANE <- "/home/genomics/jorge/files/sugarcane/"
OUT_DIR_PURPLE    <- "/home/genomics/jorge/files/purple/new/"

# ==============================================================================
# INPUT FILES — already filtered and weight-normalised
# ==============================================================================
networks <- list(
  sugarcane = list(
    edge_file = "/home/genomics/jorge/files/sugarcane/network_sugarcane_filtered_edges.tsv",
    prefix    = file.path(OUT_DIR_SUGARCANE, "leiden_sugarcane")
  ),
  purple = list(
    edge_file = "/home/genomics/jorge/files/purple/new/network_purple_filtered_edges.tsv",
    prefix    = file.path(OUT_DIR_PURPLE, "leiden_purple")
  )
)

# ==============================================================================
# HELPER: run Leiden N times, keep partition with highest modularity
# Correct signature: leiden.community(graph, resolution, n.iterations)
# ==============================================================================
run_leiden_stable <- function(g, resolution, n_iter, n_starts, seed) {
  best_mod <- -Inf
  best_mem <- NULL   # store membership vector, not the community object

  for (s in seq_len(n_starts)) {
    set.seed(seed + s)
    part <- leiden.community(
      g,
      resolution   = resolution,
      n.iterations = n_iter
    )

    # leiden.community returns a "fakeCommunities" object —
    # extract the integer membership vector first, then compute modularity
    mem <- membership(part)             # named integer vector
    q   <- igraph::modularity(g, mem)  # explicit igraph dispatch

    if (q > best_mod) {
      best_mod <- q
      best_mem <- mem
    }
  }

  cat(sprintf("  Best modularity across %d starts: %.4f\n", n_starts, best_mod))

  # Return a plain named list so the rest of the script stays unchanged
  list(membership = best_mem, modularity = best_mod)
}

# ==============================================================================
# MAIN: cluster one network
# ==============================================================================
cluster_network <- function(group_name, edge_file, prefix,
                            resolution, min_module_size, min_module_size_plot,
                            n_iter, n_starts, seed) {

  cat("\n", strrep("=", 60), "\n", sep = "")
  cat("Clustering: ", group_name, "\n", sep = "")
  cat(strrep("=", 60), "\n", sep = "")

  # ── 1. Read edge list ─────────────────────────────────────────────
  cat("Reading edge list...\n")
  edges <- fread(edge_file, header = TRUE,
                 select = c("gene1", "gene2", "weight"))

  cat(sprintf("  Edges: %s\n",
              format(nrow(edges), big.mark = ",")))
  cat(sprintf("  Weight range: [%.6f, %.6f]\n",
              min(edges$weight), max(edges$weight)))

  # ── 2. Build weighted igraph ──────────────────────────────────────
  cat("Building igraph object...\n")
  g <- graph_from_data_frame(edges[, .(gene1, gene2, weight)],
                              directed = FALSE)
  g <- simplify(g,
                remove.multiple = TRUE,
                remove.loops    = TRUE,
                edge.attr.comb  = list(weight = "max"))

  cat(sprintf("  Nodes: %s | Edges: %s\n",
              format(vcount(g), big.mark = ","),
              format(ecount(g), big.mark = ",")))
  rm(edges); gc()

  # ── 3. Run Leiden ─────────────────────────────────────────────────
  cat(sprintf("Running Leiden (resolution=%.2f, n.iterations=%d, %d starts)...\n",
              resolution, n_iter, n_starts))

  partition     <- run_leiden_stable(g, resolution, n_iter, n_starts, seed)
  n_raw_modules <- length(unique(partition$membership))
  q             <- partition$modularity
  cat(sprintf("  Raw modules: %d | Modularity Q: %.4f\n", n_raw_modules, q))

  # ── 4. Build membership table ─────────────────────────────────────
  membership_dt <- data.table(
    gene       = V(g)$name,
    module_raw = as.integer(partition$membership)
  )

  # Size-rank modules: largest = Module_001, etc.
  # Only modules with N >= min_module_size get a name; rest → "Unassigned"
  mod_sizes  <- membership_dt[, .N, by = module_raw][order(-N)]
  large_mods <- mod_sizes[N >= min_module_size, module_raw]

  rank_map <- data.table(
    module_raw  = large_mods,
    module_name = sprintf("Module_%03d", seq_along(large_mods))
  )

  membership_dt <- merge(membership_dt, rank_map, by = "module_raw", all.x = TRUE)
  membership_dt[is.na(module_name), module_name := "Unassigned"]
  membership_dt[, module_raw := NULL]

  # ── 5. Add node-level metrics ─────────────────────────────────────
  node_str <- strength(g, vids = V(g), weights = E(g)$weight)
  node_deg <- degree(g)

  membership_dt[, strength := node_str[gene]]
  membership_dt[, degree   := node_deg[gene]]
  setorder(membership_dt, module_name, -strength)

  # ── 6. Module summary (all modules, including small ones) ─────────
  summary_dt <- membership_dt[, .N, by = module_name][order(-N)]
  setnames(summary_dt, c("module", "n_genes"))
  summary_dt[, modularity_Q       := round(q, 4)]
  summary_dt[, resolution         := resolution]
  summary_dt[, min_module_size    := min_module_size]
  summary_dt[, shown_in_plot      := n_genes >= min_module_size_plot & module != "Unassigned"]

  n_modules_all    <- sum(summary_dt$module != "Unassigned")
  n_modules_plot   <- sum(summary_dt$shown_in_plot)
  n_unassigned     <- membership_dt[module_name == "Unassigned", .N]
  assigned_modules <- summary_dt[module != "Unassigned"]

  cat(sprintf("  Total modules (>= %d genes): %d\n", min_module_size, n_modules_all))
  cat(sprintf("  Modules shown in plot (>= %d genes): %d\n", min_module_size_plot, n_modules_plot))
  cat(sprintf("  Unassigned genes:    %d\n", n_unassigned))
  cat(sprintf("  Largest module:      %d genes\n", max(assigned_modules$n_genes)))
  cat(sprintf("  Median module size:  %.0f genes\n", median(assigned_modules$n_genes)))

  # ── 7. Hub genes — top 10 per module by weighted strength ─────────
  hubs_dt <- membership_dt[module_name != "Unassigned",
                            head(.SD, 10),
                            by      = module_name,
                            .SDcols = c("gene", "strength", "degree")]
  hubs_dt[, hub_rank := seq_len(.N), by = module_name]

  # ── 8. Plot data table — modules that will appear in the plot ─────
  # This is the table used to generate the plot, saved as a separate file
  # so results are fully reproducible without re-running the analysis.
  plot_dt <- summary_dt[shown_in_plot == TRUE][order(-n_genes)]
  plot_dt[, module := factor(module, levels = module)]

  # ── 9. Save all output files ──────────────────────────────────────
  out_membership <- paste0(prefix, "_membership.tsv")
  out_summary    <- paste0(prefix, "_module_summary.tsv")
  out_hubs       <- paste0(prefix, "_hub_genes.tsv")
  out_plot_data  <- paste0(prefix, "_plot_data.tsv")

  # gene → module mapping (all genes, including Unassigned)
  fwrite(membership_dt[, .(gene, module_name, strength, degree)],
         file = out_membership, sep = "\t", quote = FALSE)

  # module sizes — all modules (includes small ones and Unassigned)
  fwrite(summary_dt[, .(module, n_genes, modularity_Q,
                         resolution, min_module_size, shown_in_plot)],
         file = out_summary, sep = "\t", quote = FALSE)

  # hub genes
  fwrite(hubs_dt,
         file = out_hubs, sep = "\t", quote = FALSE)

  # table used to generate the plot (modules >= min_module_size_plot only)
  fwrite(plot_dt[, .(module = as.character(module), n_genes,
                      modularity_Q, resolution)],
         file = out_plot_data, sep = "\t", quote = FALSE)

  cat("  Saved:", basename(out_membership), "\n")
  cat("  Saved:", basename(out_summary),    "\n")
  cat("  Saved:", basename(out_hubs),       "\n")
  cat("  Saved:", basename(out_plot_data),  "\n")

  # ── 10. Module size distribution plot ────────────────────────────
  # Only modules with >= min_module_size_plot genes are plotted;
  # a dashed line marks the min_module_size_plot threshold.
  p <- ggplot(plot_dt, aes(x = module, y = n_genes)) +
    geom_col(fill = "#2c7bb6", alpha = 0.85) +
    geom_hline(yintercept = min_module_size_plot, linetype = "dashed",
               color = "firebrick", linewidth = 0.6) +
    labs(
      title    = paste("Module size distribution —", group_name),
      subtitle = sprintf(
        "Leiden resolution=%.2f | Q=%.4f | %d total modules | %d shown (≥%d genes)",
        resolution, q, n_modules_all, n_modules_plot, min_module_size_plot
      ),
      x = "Module (sorted by size)",
      y = "Number of genes"
    ) +
    theme_classic(base_size = 13) +
    theme(
      axis.text.x   = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
      plot.title    = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )

  ggsave(paste0(prefix, "_module_sizes.pdf"), plot = p,
         width = 12, height = 5, device = "pdf")
  ggsave(paste0(prefix, "_module_sizes.png"), plot = p,
         width = 12, height = 5, dpi = 300)

  cat("  Saved:", basename(paste0(prefix, "_module_sizes.pdf")), "\n")
  cat("  Saved:", basename(paste0(prefix, "_module_sizes.png")), "\n")

  cat("Clustering complete for:", group_name, "\n")

  invisible(list(
    group      = group_name,
    graph      = g,
    membership = membership_dt,
    summary    = summary_dt,
    plot_data  = plot_dt,
    hubs       = hubs_dt,
    modularity = q
  ))
}

# ==============================================================================
# ENTRY POINT
# ==============================================================================
results <- list()

for (grp in names(networks)) {
  n <- networks[[grp]]
  results[[grp]] <- cluster_network(
    group_name           = grp,
    edge_file            = n$edge_file,
    prefix               = n$prefix,
    resolution           = RESOLUTION,
    min_module_size      = MIN_MODULE_SIZE,
    min_module_size_plot = MIN_MODULE_SIZE_PLOT,
    n_iter               = N_ITERATIONS,
    n_starts             = N_STARTS,
    seed                 = RANDOM_SEED
  )
}

rds_out <- file.path(OUT_DIR_SUGARCANE, "leiden_results_both_networks.rds")
saveRDS(results, file = rds_out)
cat("\nSaved combined RDS:", basename(rds_out), "\n")
cat("All networks clustered.\n")
