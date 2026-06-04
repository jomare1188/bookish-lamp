library(data.table)
library(igraph)
library(ggplot2)

# ==============================================================================
# CONFIGURATION
# ==============================================================================
INFLATION            <- 2.0   # MCL Inflation parameter — tune if modules are too big/small
                              #   too many tiny modules -> lower inflation (1.4, 1.8)
                              #   too few huge modules  -> raise inflation (2.5, 4.0)
MIN_MODULE_SIZE      <- 2     # modules smaller than this -> labelled "Unassigned"
MIN_MODULE_SIZE_PLOT <- 10    # minimum module size shown in the plot
NUM_CORES            <- 250     # Number of CPU threads dedicated to the C MCL execution (-te flag)

OUT_DIR_SUGARCANE <- "/home/genomics/jorge/files/sugarcane/"
OUT_DIR_PURPLE    <- "/home/genomics/jorge/files/purple/new/"

# ==============================================================================
# INPUT FILES — already filtered and weight-normalised
# ==============================================================================
networks <- list(
  sugarcane = list(
    edge_file = "/home/genomics/jorge/files/sugarcane/network_sugarcane_filtered_edges.tsv",
    prefix    = file.path(OUT_DIR_SUGARCANE, "mcl_sugarcane") # changed prefix to match algorithm
  ),
  purple = list(
    edge_file = "/home/genomics/jorge/files/purple/new/network_purple_filtered_edges.tsv",
    prefix    = file.path(OUT_DIR_PURPLE, "mcl_purple")
  )
)

# ==============================================================================
# HELPER: Run native C MCL implementation with multi-threading
# ==============================================================================
run_mcl_multicore <- function(g, inflation, cores) {
  cat(sprintf("  Executing original C MCL with %d threads (Inflation = %.2f)...\n", cores, inflation))
  
  # 1. Extract edge data frame to write out in ABC format (NodeA NodeB Weight)
  edges_df <- as_data_frame(g, what = "edges")
  
  # 2. Configure temporary files for interprocess communication
  tmp_in  <- tempfile(fileext = ".abc")
  tmp_out <- tempfile(fileext = ".mcl")
  
  # 3. Write native edge list format for MCL
  fwrite(data.table(edges_df$from, edges_df$to, edges_df$weight), 
         file = tmp_in, sep = "\t", col.names = FALSE)
  
  # 4. Invoke the C-binary via system call using the multithreading flag (-te)
  status <- system2("mcl", 
                    args = c(tmp_in, "--abc", "-I", inflation, "-te", cores, "-o", tmp_out),
                    stdout = FALSE, stderr = FALSE)
  
  if (status != 0) {
    stop("MCL execution failed. Ensure that the 'mcl' binary is installed and accessible in your system PATH.")
  }
  
  if (!file.exists(tmp_out)) {
    stop("MCL completed but output cluster file was not found.")
  }
  
  # 5. Read resulting clusters (each line contains tab-separated node names)
  lines <- readLines(tmp_out)
  
  # 6. Reconstruct an igraph compatible named membership vector
  node_names <- V(g)$name
  mem <- integer(length(node_names))
  names(mem) <- node_names
  
  for (i in seq_along(lines)) {
    if (trimws(lines[i]) == "") next
    nodes <- strsplit(lines[i], "\t")[[1]]
    mem[nodes] <- i
  }
  
  # Handle fallback edge cases: if any nodes were completely skipped by MCL,
  # assign them to single-element orphan modules
  unassigned_nodes <- names(mem)[mem == 0]
  if (length(unassigned_nodes) > 0) {
    max_cluster <- length(lines)
    for (j in seq_along(unassigned_nodes)) {
      mem[unassigned_nodes[j]] <- max_cluster + j
    }
  }
  
  # Clean up system temporary files
  unlink(c(tmp_in, tmp_out))
  
  # Calculate modularity of the MCL configuration to ensure downstream compatibility
  q <- igraph::modularity(g, mem)
  
  list(membership = mem, modularity = q)
}

# ==============================================================================
# MAIN: cluster one network
# ==============================================================================
cluster_network <- function(group_name, edge_file, prefix,
                            inflation, min_module_size, min_module_size_plot,
                            n_cores) {

  cat("\n", strrep("=", 60), "\n", sep = "")
  cat("Clustering (MCL): ", group_name, "\n", sep = "")
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

  # ── 3. Run MCL ────────────────────────────────────────────────────
  partition     <- run_mcl_multicore(g, inflation, n_cores)
  n_raw_modules <- length(unique(partition$membership))
  q             <- partition$modularity
  cat(sprintf("  Raw modules: %d | Modularity Q: %.4f\n", n_raw_modules, q))

  # ── 4. Build membership table ─────────────────────────────────────
  membership_dt <- data.table(
    gene       = V(g)$name,
    module_raw = as.integer(partition$membership)
  )

  # Size-rank modules: largest = Module_001, etc.
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
  summary_dt[, modularity_Q        := round(q, 4)]
  
  # Kept column name as 'resolution' for file schema compatibility, but populated with inflation value
  summary_dt[, resolution          := inflation] 
  summary_dt[, min_module_size     := min_module_size]
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

  # ── 8. Plot data table ────────────────────────────────────────────
  plot_dt <- summary_dt[shown_in_plot == TRUE][order(-n_genes)]
  plot_dt[, module := factor(module, levels = module)]

  # ── 9. Save all output files ──────────────────────────────────────
  out_membership <- paste0(prefix, "_membership.tsv")
  out_summary    <- paste0(prefix, "_module_summary.tsv")
  out_hubs       <- paste0(prefix, "_hub_genes.tsv")
  out_plot_data  <- paste0(prefix, "_plot_data.tsv")

  fwrite(membership_dt[, .(gene, module_name, strength, degree)],
         file = out_membership, sep = "\t", quote = FALSE)

  fwrite(summary_dt[, .(module, n_genes, modularity_Q,
                        resolution, min_module_size, shown_in_plot)],
         file = out_summary, sep = "\t", quote = FALSE)

  fwrite(hubs_dt,
         file = out_hubs, sep = "\t", quote = FALSE)

  fwrite(plot_dt[, .(module = as.character(module), n_genes,
                     modularity_Q, resolution)],
         file = out_plot_data, sep = "\t", quote = FALSE)

  cat("  Saved:", basename(out_membership), "\n")
  cat("  Saved:", basename(out_summary),    "\n")
  cat("  Saved:", basename(out_hubs),       "\n")
  cat("  Saved:", basename(out_plot_data),  "\n")

  # ── 10. Module size distribution plot ────────────────────────────
  p <- ggplot(plot_dt, aes(x = module, y = n_genes)) +
    geom_col(fill = "#1b9e77", alpha = 0.85) + # Color tweaked slightly to visually identify MCL run
    geom_hline(yintercept = min_module_size_plot, linetype = "dashed",
               color = "firebrick", linewidth = 0.6) +
    labs(
      title    = paste("Module size distribution (MCL) —", group_name),
      subtitle = sprintf(
        "MCL Inflation=%.2f | Q=%.4f | %d total modules | %d shown (≥%d genes)",
        inflation, q, n_modules_all, n_modules_plot, min_module_size_plot
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
    inflation            = INFLATION,
    min_module_size      = MIN_MODULE_SIZE,
    min_module_size_plot = MIN_MODULE_SIZE_PLOT,
    n_cores              = NUM_CORES
  )
}

rds_out <- file.path(OUT_DIR_SUGARCANE, "mcl_results_both_networks.rds")
saveRDS(results, file = rds_out)
cat("\nSaved combined RDS:", basename(rds_out), "\n")
cat("All networks clustered via MCL.\n")
