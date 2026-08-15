# =============================================================================
# 06_conservation_join.r — which edges of network A survive in network B?
#
# For every edge (i, j) in network A, project i and j through orthology and ask
# whether ANY resulting ortholog pair is an edge in network B. Writes:
#
#   conserved_edges_<A>_to_<B>_FULL.tsv     gene1, gene2, source, conserved
#   conserved_genes_<A>_FULL.txt            genes touching a conserved edge
#   conservation_summary_<A>_to_<B>_FULL.tsv  totals, overall + BY LAYER
#
# WHY A JOIN AND NOT A MATRIX PRODUCT: the matrix formulation
# (O %*% A_b %*% t(O)) > 0 materialises a sparse product over the entire
# gene-pair space of B. At full-network scale -- purple has ~710 M edges and a
# mean degree in the thousands -- that is not affordable. The chunked join costs
# what the actual ortholog fan-out costs instead, and streams, so memory is
# bounded by CHUNK_SIZE rather than by network size.
#
# THE BY-LAYER BREAKDOWN is the point of the whole augmentation: it answers
# whether NON-LINEAR co-expression (edges only the MI layer found) is conserved
# between the two species at the same rate as linear co-expression. Computing it
# here is free -- the `source` column is already being read.
#
# RUN TWICE, once per direction:
#   ./run.sh conserve sugarcane_to_purple    # cheaper outer loop — run first
#   ./run.sh conserve purple_to_sugarcane    # ~9x more outer edges
# =============================================================================

suppressMessages(library(data.table))

source(file.path(dirname(sub("--file=", "",
       grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "lib", "common.R"))

DIRECTION   <- env_req("CLEAN_DIRECTION")
EDGES_SC    <- env_req("CLEAN_EDGES_SUGARCANE")
EDGES_PU    <- env_req("CLEAN_EDGES_PURPLE")
OUT_DIR     <- ensure_dir(env_req("CLEAN_OUT_DIR"))
ORTHOGROUPS <- env_req("CLEAN_ORTHOGROUPS")
CHUNK_SIZE  <- as.integer(env_num("CLEAN_CHUNK_SIZE", 2e6))
setDTthreads(as.integer(env_num("CLEAN_CORES", 100)))

# OrthoFinder's species names, as they appear in Orthogroups.tsv.
SPECIES_SUGARCANE <- env_opt("CLEAN_OG_SPECIES_SUGARCANE", "sugarcane_one_transcript")
SPECIES_PURPLE    <- env_opt("CLEAN_OG_SPECIES_PURPLE", "one_transcript_purple_proteins")

banner(paste("conservation join:", DIRECTION))
assert_network_schema(EDGES_SC)
assert_network_schema(EDGES_PU)

# =============================================================================
# 1. orthologs
# =============================================================================
say("loading OrthoFinder results")
og <- as.data.table(read_orthogroups_tsv(ORTHOGROUPS))
og[, Gene := sub("\\.p[0-9]+$", "", Gene)]
og[, Gene := strip_version(Gene)]      # network ids are bare; proteome ids are not

if (!SPECIES_SUGARCANE %in% og$Species || !SPECIES_PURPLE %in% og$Species)
  stop("Orthogroups.tsv does not contain the expected species names.\n",
       "  looking for: ", SPECIES_SUGARCANE, " and ", SPECIES_PURPLE, "\n",
       "  found: ", paste(unique(og$Species), collapse = ", "), call. = FALSE)

og[Species == SPECIES_SUGARCANE, Species := "sugarcane"]
og[Species == SPECIES_PURPLE,    Species := "purple"]

ortholog_pairs <- merge(
  og[Species == "sugarcane", .(Orthogroup, sugarcane_gene = Gene)],
  og[Species == "purple",    .(Orthogroup, purple_gene    = Gene)],
  by = "Orthogroup", allow.cartesian = TRUE)
say("  ortholog pairs (many-to-many): ", fmt_n(nrow(ortholog_pairs)))
say(sprintf("  mean orthologs per gene: sugarcane %.2f | purple %.2f",
            ortholog_pairs[, .N, by = sugarcane_gene][, mean(N)],
            ortholog_pairs[, .N, by = purple_gene][, mean(N)]))

# =============================================================================
# 2. direction
# =============================================================================
if (DIRECTION == "sugarcane_to_purple") {
  edge_file_a <- EDGES_SC; edge_file_b <- EDGES_PU
  pairs_ab <- ortholog_pairs[, .(gene_a = sugarcane_gene, gene_b = purple_gene)]
  label_a <- "sugarcane"; label_b <- "purple"
} else if (DIRECTION == "purple_to_sugarcane") {
  edge_file_a <- EDGES_PU; edge_file_b <- EDGES_SC
  pairs_ab <- ortholog_pairs[, .(gene_a = purple_gene, gene_b = sugarcane_gene)]
  label_a <- "purple"; label_b <- "sugarcane"
} else {
  stop("CLEAN_DIRECTION must be sugarcane_to_purple or purple_to_sugarcane",
       call. = FALSE)
}
setkey(pairs_ab, gene_a)

out_edges   <- file.path(OUT_DIR, sprintf("conserved_edges_%s_to_%s_FULL.tsv", label_a, label_b))
out_genes   <- file.path(OUT_DIR, sprintf("conserved_genes_%s_FULL.txt", label_a))
out_summary <- file.path(OUT_DIR, sprintf("conservation_summary_%s_to_%s_FULL.tsv", label_a, label_b))

# =============================================================================
# 3. target adjacency
# =============================================================================
say("loading target network (", label_b, ") into an adjacency lookup")
edges_b <- fread(edge_file_b, header = TRUE, select = c("gene1", "gene2"))
say("  ", fmt_n(nrow(edges_b)), " edges")
adj_b <- rbindlist(list(edges_b[, .(g1 = gene1, g2 = gene2)],
                        edges_b[, .(g1 = gene2, g2 = gene1)]))
rm(edges_b); invisible(gc())
setkey(adj_b, g1, g2)
say("  adjacency rows (both directions): ", fmt_n(nrow(adj_b)))

# An empty intersection between the ortholog table and the network ids means
# the id conventions drifted. The old script printed samples and asked the
# operator to eyeball them; this fails instead, because the silent result is a
# perfectly well-formed file saying nothing is conserved.
probe <- fread(edge_file_a, header = TRUE, nrows = 200000L, select = "gene1")
overlap <- sum(unique(probe$gene1) %in% pairs_ab$gene_a)
if (overlap == 0L)
  stop("none of the first 200k genes in ", label_a, " appear in the ortholog ",
       "table — the id conventions do not match.\n",
       "  network id: ", probe$gene1[1], "\n",
       "  ortholog id: ", pairs_ab$gene_a[1], call. = FALSE)
say(sprintf("  id sanity: %s of %s sampled %s genes are in the ortholog table",
            fmt_n(overlap), fmt_n(uniqueN(probe$gene1)), label_a))
rm(probe)

# =============================================================================
# 4. stream network A
# =============================================================================
n_total <- tryCatch(
  as.integer(system(sprintf("wc -l < %s", shQuote(edge_file_a)), intern = TRUE)) - 1L,
  error = function(e) NA_integer_)
say("streaming network ", label_a, " (",
    if (is.na(n_total)) "unknown total" else fmt_n(n_total), " edges)")

con <- file(edge_file_a, open = "r")
on.exit(close(con), add = TRUE)
hdr <- strsplit(readLines(con, n = 1L, warn = FALSE), "\t", fixed = TRUE)[[1L]]
col_source <- match("source", hdr)
stopifnot(identical(hdr[1:2], c("gene1", "gene2")), !is.na(col_source))

writeLines("gene1\tgene2\tsource\tconserved", out_edges)

gene_chunks     <- list()      # collected, uniqued once — union() per chunk is
by_source       <- list()      # quadratic in the accumulated vector
total_conserved <- 0L
total_processed <- 0L
chunk_i         <- 0L

repeat {
  lines <- readLines(con, n = CHUNK_SIZE, warn = FALSE)
  if (length(lines) == 0L) break
  chunk_i <- chunk_i + 1L

  chunk <- fread(text = lines, header = FALSE,
                 select = c(1L, 2L, col_source),
                 col.names = c("gene1", "gene2", "source"))
  rm(lines)
  chunk[, edge_id := .I]

  map1 <- merge(chunk[, .(edge_id, gene1)], pairs_ab,
                by.x = "gene1", by.y = "gene_a", allow.cartesian = TRUE)
  setnames(map1, "gene_b", "ortho1")
  map2 <- merge(chunk[, .(edge_id, gene2)], pairs_ab,
                by.x = "gene2", by.y = "gene_a", allow.cartesian = TRUE)
  setnames(map2, "gene_b", "ortho2")

  cand <- merge(map1[, .(edge_id, ortho1)], map2[, .(edge_id, ortho2)],
                by = "edge_id", allow.cartesian = TRUE)
  rm(map1, map2)

  matched <- merge(cand, adj_b, by.x = c("ortho1", "ortho2"),
                   by.y = c("g1", "g2"))
  conserved_ids <- unique(matched$edge_id)
  rm(cand, matched)

  chunk[, conserved := edge_id %in% conserved_ids]
  total_conserved <- total_conserved + length(conserved_ids)
  total_processed <- total_processed + nrow(chunk)

  by_source[[chunk_i]] <- chunk[, .(n = .N, n_conserved = sum(conserved)),
                                by = source]

  fwrite(chunk[, .(gene1, gene2, source, conserved)], file = out_edges,
         sep = "\t", quote = FALSE, append = TRUE, col.names = FALSE)

  if (any(chunk$conserved))
    gene_chunks[[length(gene_chunks) + 1L]] <-
      unique(c(chunk[conserved == TRUE, gene1], chunk[conserved == TRUE, gene2]))

  pct <- if (is.na(n_total)) NA_real_ else 100 * total_processed / n_total
  say(sprintf("  chunk %d | %s edges%s | conserved %s (%.2f%%)",
              chunk_i, fmt_n(total_processed),
              if (is.na(pct)) "" else sprintf(" (%.1f%%)", pct),
              fmt_n(total_conserved),
              100 * total_conserved / max(total_processed, 1L)))

  rm(chunk); invisible(gc())
}
close(con); on.exit()

conserved_gene_set <- unique(unlist(gene_chunks, use.names = FALSE))
writeLines(conserved_gene_set, out_genes)
say("wrote ", basename(out_genes), "  (", fmt_n(length(conserved_gene_set)), " genes)")

# =============================================================================
# 5. summary, overall and by layer
# =============================================================================
src <- rbindlist(by_source)[, .(total_edges = sum(n),
                                conserved_edges = sum(n_conserved)), by = source]
src[, conserved_fraction := round(conserved_edges / total_edges, 6)]
setorder(src, -total_edges)

summary_dt <- rbindlist(list(
  data.table(direction = sprintf("%s_to_%s", label_a, label_b),
             layer = "ALL",
             total_edges = total_processed,
             conserved_edges = total_conserved,
             conserved_fraction = round(total_conserved / max(total_processed, 1L), 6),
             conserved_genes = length(conserved_gene_set)),
  data.table(direction = sprintf("%s_to_%s", label_a, label_b),
             layer = src$source,
             total_edges = src$total_edges,
             conserved_edges = src$conserved_edges,
             conserved_fraction = src$conserved_fraction,
             conserved_genes = NA_integer_)
))
write_tsv(summary_dt, out_summary)

banner("conservation by layer")
print(summary_dt, row.names = FALSE)
say("")
say("Read the `mi` row against the `pearson` row: that is whether non-linear")
say("co-expression is conserved between the species at the same rate as linear.")
say("done: ", DIRECTION)
