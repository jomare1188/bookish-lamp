# =============================================================================
# 05_mcl_clustering.r — MCL modules for one network
#
# Writes, per study:
#   mcl_<study>_membership.tsv      gene, module_name, strength, degree
#   mcl_<study>_module_summary.tsv  module, n_genes, modularity_Q, ...
#   mcl_<study>_hub_genes.tsv       top 10 by weighted strength per module
#   mcl_<study>_plot_data.tsv
#   mcl_<study>_module_sizes.{pdf,png}
#
# Clustering is delegated to the native C `mcl` binary via a temporary .abc
# file. That file is roughly the size of the edge list -- tens of GB for purple
# -- so run.sh points TMPDIR at /dados04 rather than letting it land on a small
# system /tmp.
#
# ONE STUDY PER INVOCATION. The old script looped over both and accumulated the
# full igraph objects for each in a list, so peak RAM held sugarcane's graph and
# purple's 710 M-edge graph at the same time, and a failure on the second lost
# the first. Per-study runs also make the "smoke-test on sugarcane first" advice
# actually possible.
#
# RUN: through run.sh  ->  ./run.sh mcl sugarcane
# =============================================================================

suppressMessages({
  library(data.table)
  library(igraph)
  library(ggplot2)
})

source(file.path(dirname(sub("--file=", "",
       grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "lib", "common.R"))

STUDY      <- env_req("CLEAN_STUDY")
EDGE_FILE  <- env_req("CLEAN_EDGES")
PREFIX     <- env_req("CLEAN_PREFIX")
INFLATION  <- env_num("CLEAN_INFLATION", 2)
MIN_SIZE   <- as.integer(env_num("CLEAN_MIN_MODULE_SIZE", 5))
MIN_PLOT   <- as.integer(env_num("CLEAN_MIN_MODULE_SIZE_PLOT", 10))
CORES      <- as.integer(env_num("CLEAN_CORES", 100))
SAVE_RDS   <- env_flag("CLEAN_SAVE_GRAPH_RDS", FALSE)
setDTthreads(CORES)

banner(paste("MCL clustering:", STUDY))
assert_network_schema(EDGE_FILE)

# =============================================================================
# native C MCL
# =============================================================================
run_mcl <- function(g, inflation, cores) {
  say(sprintf("running C mcl: inflation %.2f, %d threads", inflation, cores))
  edges_df <- as_data_frame(g, what = "edges")

  tmp_in  <- tempfile(fileext = ".abc")
  tmp_out <- tempfile(fileext = ".mcl")
  on.exit(unlink(c(tmp_in, tmp_out)), add = TRUE)

  fwrite(data.table(edges_df$from, edges_df$to, edges_df$weight),
         file = tmp_in, sep = "\t", col.names = FALSE)
  rm(edges_df); invisible(gc())
  say("  wrote ", round(file.size(tmp_in) / 1e9, 2), " GB of .abc to ",
      dirname(tmp_in))

  status <- system2("mcl",
                    args = c(tmp_in, "--abc", "-I", inflation,
                             "-te", cores, "-o", tmp_out),
                    stdout = FALSE, stderr = FALSE)
  if (status != 0)
    stop("mcl exited with status ", status,
         ". Is the 'mcl' binary on PATH, and is TMPDIR (", dirname(tmp_in),
         ") big enough?", call. = FALSE)
  if (!file.exists(tmp_out))
    stop("mcl reported success but wrote no cluster file", call. = FALSE)

  lines <- readLines(tmp_out)
  node_names <- V(g)$name
  mem <- integer(length(node_names))
  names(mem) <- node_names
  for (i in seq_along(lines)) {
    if (trimws(lines[i]) == "") next
    mem[strsplit(lines[i], "\t")[[1L]]] <- i
  }

  # Any node mcl skipped becomes its own singleton, so membership stays total.
  orphans <- names(mem)[mem == 0L]
  if (length(orphans))
    mem[orphans] <- length(lines) + seq_along(orphans)

  list(membership = mem, modularity = igraph::modularity(g, mem))
}

# =============================================================================
say("reading ", basename(EDGE_FILE))
edges <- fread(EDGE_FILE, header = TRUE, select = c("gene1", "gene2", "weight"))
say("  ", fmt_n(nrow(edges)), " edges | weight range [",
    sprintf("%.6f", min(edges$weight)), ", ",
    sprintf("%.6f", max(edges$weight)), "]")

g <- graph_from_data_frame(edges, directed = FALSE)
g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE,
              edge.attr.comb = list(weight = "max"))
rm(edges); invisible(gc())
say("nodes ", fmt_n(vcount(g)), " | edges ", fmt_n(ecount(g)))

partition <- run_mcl(g, INFLATION, CORES)
q <- partition$modularity
say(sprintf("raw modules: %d | modularity Q: %.4f",
            length(unique(partition$membership)), q))

# --- membership, ranked so the largest module is Module_001 ------------------
membership_dt <- data.table(gene       = V(g)$name,
                            module_raw = as.integer(partition$membership))
mod_sizes  <- membership_dt[, .N, by = module_raw][order(-N)]
large_mods <- mod_sizes[N >= MIN_SIZE, module_raw]
rank_map   <- data.table(module_raw  = large_mods,
                         module_name = sprintf("Module_%03d", seq_along(large_mods)))
membership_dt <- merge(membership_dt, rank_map, by = "module_raw", all.x = TRUE)
membership_dt[is.na(module_name), module_name := "Unassigned"]
membership_dt[, module_raw := NULL]

membership_dt[, strength := strength(g, vids = V(g), weights = E(g)$weight)[gene]]
membership_dt[, degree   := degree(g)[gene]]
setorder(membership_dt, module_name, -strength)

# --- summary -----------------------------------------------------------------
summary_dt <- membership_dt[, .N, by = module_name][order(-N)]
setnames(summary_dt, c("module", "n_genes"))
summary_dt[, modularity_Q    := round(q, 4)]
summary_dt[, inflation       := INFLATION]
summary_dt[, min_module_size := MIN_SIZE]
summary_dt[, shown_in_plot   := n_genes >= MIN_PLOT & module != "Unassigned"]

assigned <- summary_dt[module != "Unassigned"]
say(sprintf("modules (>= %d genes): %d | shown in plot (>= %d): %d",
            MIN_SIZE, nrow(assigned), MIN_PLOT, sum(summary_dt$shown_in_plot)))
say(sprintf("unassigned genes: %s | largest module: %s | median size: %.0f",
            fmt_n(membership_dt[module_name == "Unassigned", .N]),
            fmt_n(max(assigned$n_genes)), median(assigned$n_genes)))

hubs_dt <- membership_dt[module_name != "Unassigned", head(.SD, 10),
                         by = module_name, .SDcols = c("gene", "strength", "degree")]
hubs_dt[, hub_rank := seq_len(.N), by = module_name]

plot_dt <- summary_dt[shown_in_plot == TRUE][order(-n_genes)]
plot_dt[, module := factor(module, levels = module)]

write_tsv(membership_dt[, .(gene, module_name, strength, degree)],
          paste0(PREFIX, "_membership.tsv"))
write_tsv(summary_dt[, .(module, n_genes, modularity_Q, inflation,
                         min_module_size, shown_in_plot)],
          paste0(PREFIX, "_module_summary.tsv"))
write_tsv(hubs_dt, paste0(PREFIX, "_hub_genes.tsv"))
write_tsv(plot_dt[, .(module = as.character(module), n_genes,
                      modularity_Q, inflation)],
          paste0(PREFIX, "_plot_data.tsv"))

p <- ggplot(plot_dt, aes(x = module, y = n_genes)) +
  geom_col(fill = "#1b9e77", alpha = 0.85) +
  geom_hline(yintercept = MIN_PLOT, linetype = "dashed",
             color = "firebrick", linewidth = 0.6) +
  labs(title    = paste("Module size distribution (MCL) —", STUDY),
       subtitle = sprintf("Inflation=%.2f | Q=%.4f | %d modules | %d shown (>=%d genes)",
                          INFLATION, q, nrow(assigned),
                          sum(summary_dt$shown_in_plot), MIN_PLOT),
       x = "Module (sorted by size)", y = "Number of genes") +
  theme_classic(base_size = 13) +
  theme(axis.text.x   = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
        plot.title    = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
ggsave(paste0(PREFIX, "_module_sizes.pdf"), p, width = 12, height = 5)
ggsave(paste0(PREFIX, "_module_sizes.png"), p, width = 12, height = 5, dpi = 300)

# The old pipeline always saved an .rds holding the full igraph objects for both
# networks -- 8 GB, written into the sugarcane directory regardless of study.
# Nothing in the pipeline reads it (every consumer reads _membership.tsv), so it
# is off unless asked for.
if (SAVE_RDS) {
  rds <- paste0(PREFIX, "_graph.rds")
  saveRDS(list(study = STUDY, graph = g, membership = membership_dt,
               summary = summary_dt, modularity = q), rds)
  say("wrote ", basename(rds), " (", round(file.size(rds) / 1e9, 2), " GB)")
}

say("done: ", STUDY)
