# =============================================================================
# 13_conservation_null.r — is the observed edge-conservation rate above chance?
#
# 06_conservation_join.r reports that 10.76% of sugarcane edges have a purple
# counterpart and 1.52% of purple edges have a sugarcane one. Those rates are
# NOT interpretable on their own, because both are driven by things that have
# nothing to do with conservation:
#
#   * the target network's DENSITY -- purple's is 0.0484 against sugarcane's
#     0.0143, so a random ortholog pair is 3.4x more likely to "be an edge" in
#     purple no matter what the biology says. This is also the whole reason the
#     two directions are not comparable to each other.
#   * ORTHOLOGY FAN-OUT -- a gene with 5 orthologs gets 25 chances to match, one
#     with 1 ortholog gets 1.
#   * orthology COVERAGE -- ~29% of genes have no ortholog at all and can never
#     be conserved, which drags every rate down by a constant factor.
#
# THE NULL: permute the target-gene column of the ortholog table. That destroys
# the true orthology assignment while preserving exactly:
#   - each source gene's number of orthologs (fan-out)
#   - which source genes have any ortholog at all (coverage)
#   - the multiset of target genes used, hence aggregate exposure to the target
#     network's degree distribution
# so the only thing removed is *which* target gene each source gene maps to.
#
# If the observed rate sits inside this null, the conservation signal is an
# artifact of density and fan-out. If it sits far above, it is real.
#
# The per-LAYER breakdown is computed under the same null, so the mi-vs-pearson
# contrast that the whole augmentation rests on gets its own significance test
# rather than resting on a z-test that assumes edges are independent (they are
# emphatically not -- they share genes).
#
# COST AND SPACE: the target adjacency is built ONCE and reused for every
# replicate, and the source network is Bernoulli-sampled rather than streamed in
# full, so a replicate costs seconds instead of the ~11-48 min a full direction
# takes. Nothing but a small summary table is written -- no per-edge output.
#
# RUN:  ./run.sh conservenull sugarcane_to_purple
#       ./run.sh conservenull purple_to_sugarcane
# =============================================================================

suppressMessages(library(data.table))

source(file.path(dirname(sub("--file=", "",
       grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "lib", "common.R"))

DIRECTION   <- env_req("CLEAN_DIRECTION")
EDGES_SC    <- env_req("CLEAN_EDGES_SUGARCANE")
EDGES_PU    <- env_req("CLEAN_EDGES_PURPLE")
OUT_DIR     <- ensure_dir(env_req("CLEAN_OUT_DIR"))
ORTHOGROUPS <- env_req("CLEAN_ORTHOGROUPS")
N_REPS      <- as.integer(env_num("CLEAN_NULL_REPS", 20))
N_SAMPLE    <- as.numeric(env_num("CLEAN_NULL_SAMPLE", 5e6))
SEED        <- as.integer(env_num("CLEAN_SEED", 1188))
setDTthreads(as.integer(env_num("CLEAN_CORES", 100)))

SPECIES_SUGARCANE <- env_opt("CLEAN_OG_SPECIES_SUGARCANE", "sugarcane_one_transcript")
SPECIES_PURPLE    <- env_opt("CLEAN_OG_SPECIES_PURPLE", "one_transcript_purple_proteins")

set.seed(SEED)
banner(paste("conservation permutation null:", DIRECTION))

# --- orthologs ---------------------------------------------------------------
say("loading orthologs")
og <- as.data.table(read_orthogroups_tsv(ORTHOGROUPS))
og[, Gene := strip_version(sub("\\.p[0-9]+$", "", Gene))]
og[Species == SPECIES_SUGARCANE, Species := "sugarcane"]
og[Species == SPECIES_PURPLE,    Species := "purple"]
op <- merge(og[Species == "sugarcane", .(Orthogroup, sugarcane_gene = Gene)],
            og[Species == "purple",    .(Orthogroup, purple_gene    = Gene)],
            by = "Orthogroup", allow.cartesian = TRUE)

if (DIRECTION == "sugarcane_to_purple") {
  edge_a <- EDGES_SC; edge_b <- EDGES_PU
  pairs  <- op[, .(gene_a = sugarcane_gene, gene_b = purple_gene)]
  label_a <- "sugarcane"; label_b <- "purple"
} else if (DIRECTION == "purple_to_sugarcane") {
  edge_a <- EDGES_PU; edge_b <- EDGES_SC
  pairs  <- op[, .(gene_a = purple_gene, gene_b = sugarcane_gene)]
  label_a <- "purple"; label_b <- "sugarcane"
} else stop("CLEAN_DIRECTION must be sugarcane_to_purple or purple_to_sugarcane",
            call. = FALSE)
say("  ortholog pairs: ", fmt_n(nrow(pairs)))

# --- target adjacency, built once and reused for every replicate -------------
say("building ", label_b, " adjacency (once, reused for all replicates)")
eb <- fread(edge_b, header = TRUE, select = c("gene1", "gene2"))
adj <- rbindlist(list(eb[, .(g1 = gene1, g2 = gene2)],
                      eb[, .(g1 = gene2, g2 = gene1)]))
rm(eb); invisible(gc())
setkey(adj, g1, g2)
say("  adjacency rows: ", fmt_n(nrow(adj)))

# --- Bernoulli sample of the source network ----------------------------------
# Sampling rather than streaming all of A: at 5e6 edges the standard error on a
# 10% rate is 0.013%, far finer than the effect being tested, and it turns a
# ~45 min pass into ~30 s so replicates are affordable. Sampled with awk so the
# full edge file is never held in memory.
n_total <- as.numeric(system(sprintf("wc -l < %s", shQuote(edge_a)), intern = TRUE)) - 1
frac <- min(1, N_SAMPLE / n_total)
say(sprintf("sampling %s of %s %s edges (%.3f)",
            fmt_n(round(frac * n_total)), fmt_n(n_total), label_a, frac))

tmp <- tempfile(fileext = ".tsv")
on.exit(unlink(tmp), add = TRUE)
hdr <- strsplit(readLines(edge_a, n = 1L), "\t", fixed = TRUE)[[1L]]
col_source <- match("source", hdr)
stopifnot(identical(hdr[1:2], c("gene1", "gene2")), !is.na(col_source))
system(sprintf(
  "awk -F'\\t' -v OFS='\\t' -v p=%.10f -v s=%d 'BEGIN{srand(%d)} NR>1 && rand()<p {print $1,$2,$%d}' %s > %s",
  frac, col_source, SEED, col_source, shQuote(edge_a), shQuote(tmp)))
A <- fread(tmp, header = FALSE, col.names = c("gene1", "gene2", "source"))
unlink(tmp)
say("  sampled ", fmt_n(nrow(A)), " edges  (",
    paste(sprintf("%s %s", A[, .N, by = source]$source,
                  fmt_n(A[, .N, by = source]$N)), collapse = " | "), ")")
A[, edge_id := .I]

# --- one pass: how many sampled edges are conserved, overall and by layer ----
conserved_counts <- function(pr) {
  setkey(pr, gene_a)
  m1 <- merge(A[, .(edge_id, gene1)], pr, by.x = "gene1", by.y = "gene_a",
              allow.cartesian = TRUE)[, .(edge_id, o1 = gene_b)]
  m2 <- merge(A[, .(edge_id, gene2)], pr, by.x = "gene2", by.y = "gene_a",
              allow.cartesian = TRUE)[, .(edge_id, o2 = gene_b)]
  cand <- merge(m1, m2, by = "edge_id", allow.cartesian = TRUE)
  rm(m1, m2)
  hit <- unique(merge(cand, adj, by.x = c("o1", "o2"), by.y = c("g1", "g2"))$edge_id)
  rm(cand)
  tot <- A[, .(n = .N), by = source]
  con <- A[edge_id %in% hit, .(k = .N), by = source]
  res <- merge(tot, con, by = "source", all.x = TRUE)
  res[is.na(k), k := 0L]
  rbind(data.table(source = "ALL", n = nrow(A), k = length(hit)), res)
}

say("computing observed")
obs <- conserved_counts(copy(pairs))
print(obs[, .(layer = source, edges = n, conserved = k,
              rate = sprintf("%.4f%%", 100 * k / n))], row.names = FALSE)

# --- the null ----------------------------------------------------------------
say("running ", N_REPS, " permutation replicates")
null_list <- vector("list", N_REPS)
t0 <- Sys.time()
for (r in seq_len(N_REPS)) {
  pr <- copy(pairs)
  pr[, gene_b := sample(gene_b)]        # break the assignment, keep the margins
  nc <- conserved_counts(pr)
  nc[, rep := r]
  null_list[[r]] <- nc
  if (r %% 5 == 0 || r == N_REPS)
    say(sprintf("  rep %d/%d  (%.1f min elapsed)", r, N_REPS,
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
nulls <- rbindlist(null_list)

# --- summarise ---------------------------------------------------------------
summ <- nulls[, .(null_mean_rate = mean(k / n),
                  null_sd_rate   = sd(k / n),
                  null_min_rate  = min(k / n),
                  null_max_rate  = max(k / n)), by = source]
summ <- merge(obs[, .(source, edges = n, conserved = k, obs_rate = k / n)],
              summ, by = "source")
summ[, fold_over_null := obs_rate / null_mean_rate]
summ[, z := (obs_rate - null_mean_rate) / null_sd_rate]
# empirical p: how often did a replicate reach the observed rate?
# Empirical p: how often did a replicate reach the observed rate? With N_REPS
# replicates the smallest attainable value is 1/(N_REPS+1), so "p = 0.048" here
# means "no replicate got close", not a calibrated 4.8%.
summ[, p_emp := {
  o <- obs[match(summ$source, obs$source), k / n]
  vapply(seq_len(.N), function(i) {
    nr <- nulls[source == summ$source[i], k / n]
    (sum(nr >= o[i]) + 1) / (length(nr) + 1)
  }, numeric(1))
}]
setorder(summ, -edges)
summ[, direction := DIRECTION]
setcolorder(summ, c("direction", "source", "edges", "conserved", "obs_rate",
                    "null_mean_rate", "null_sd_rate", "fold_over_null", "z", "p_emp"))

write_tsv(summ, file.path(OUT_DIR, sprintf("conservation_null_%s.tsv", DIRECTION)))
banner("observed vs permuted orthology")
print(summ[, .(layer = source, obs = sprintf("%.4f%%", 100 * obs_rate),
               null = sprintf("%.4f%%", 100 * null_mean_rate),
               sd = sprintf("%.4f%%", 100 * null_sd_rate),
               fold = round(fold_over_null, 2), z = round(z, 1),
               p = signif(p_emp, 3))], row.names = FALSE)
say("")
say("fold > 1 means the true orthology assignment finds more conserved edges")
say("than a random one with the same fan-out and the same target network.")
say("done: ", DIRECTION)
