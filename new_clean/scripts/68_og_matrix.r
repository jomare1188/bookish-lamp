#!/usr/bin/env Rscript
# =============================================================================
# 68_og_matrix.r -- the expression matrix with the ORTHOGROUP as the row, for
# both species at once, on ONE shared set of rows.
#
# WHY THE ORTHOGROUP. 67 measured it: log degree has an ICC of only 0.242
# (sugarcane) and 0.226 (purple) over orthogroups, against 0.940 / 0.939 for
# per-copy omega. Sequence constraint is a property of the orthogroup; network
# position is not. So gene-level cross-species wiring conservation was asking the
# data for something it does not carry, and the ortholog mapping that 61 and 62
# perform is many-to-many -- which is what inflated the conserved-pair count to
# 1.25x by multiplicity.
#
# WHY BOTH SPECIES IN ONE STAGE, breaking the one-study-per-invocation habit of
# every other stage here. The entire point is that the two networks end up on
# IDENTICAL VERTICES, so that edge conservation becomes a direct two-graph
# comparison with no projection at all. A shared row set is a joint property of
# the two species and cannot be produced one study at a time: the CV filter drops
# different groups in each, so filtering per species and hoping would give two
# different row sets. Here the intersection is taken explicitly and the cost is
# printed.
#
# THE AGGREGATION RULE IS SUM OF RAW COUNTS, THEN TRANSFORM. VST is log-like, so
# averaging VST values is not averaging expression and summing them is meaningless.
# Counts add, so they are what is added -- and the sum is also the quantity salmon
# can actually resolve: 65 showed that when the EM cannot separate two copies it
# moves reads between them while THEIR SUM STAYS PUT, which is precisely why
# collapsing removes uncertainty (infrv_gain median 1.43 / 2.28, above 1 for
# 68.7% / 77.2% of multi-copy groups).
#
# THE TRANSFORM IS varianceStabilizingTransformation(), NOT vst(). Measured
# 2026-10-01: varianceStabilizingTransformation reproduces nf-core's stored
# assay(dds, "vst") at max|diff| = 0 and cor = 1.0000000000, while the faster
# vst() approximation differs by a MEDIAN of 1.706 units. 01_export_vst.r reads
# the stored assay and never recomputes; this stage must recompute, because the
# rows are new, so it has to recompute the same way.
#
# WHAT THIS IS NOT. Not a replacement for the gene-level matrix: results/<study>/
# vst/ is untouched and stays the primary analysis. This writes a PARALLEL tree.
#
# RUN: through run.sh -> ./run.sh ogmatrix
# =============================================================================

suppressMessages({
  library(DESeq2)
  library(data.table)
  library(jsonlite)
})

source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE),
                                                 value = TRUE)[1])), "lib", "common.R"))

STUDIES    <- env_list("CLEAN_STUDIES", c("sugarcane", "purple"))
ORTHO_FILE <- env_req("CLEAN_ORTHOGROUPS")
MIN_CV     <- env_num("CLEAN_MIN_CV", 15)
BOTH_ONLY  <- env_flag("CLEAN_OG_BOTH_SPECIES_ONLY", TRUE)

cfg_for <- function(study, key) {
  v <- Sys.getenv(sprintf("CLEAN_%s_%s", key, toupper(study)))
  if (!nzchar(v))
    stop("CLEAN_", key, "_", toupper(study), " is not set -- launch through run.sh",
         call. = FALSE)
  v
}

banner("orthogroup expression matrix, both species, one shared row set")

# --- the orthogroup partition, per species ----------------------------------
og <- as.data.table(read_orthogroups_tsv(ORTHO_FILE))
og[, Gene := strip_version(sub("\\.p[0-9]+$", "", Gene))]

raw <- list(); genes <- list(); samples <- list(); map <- list()
for (s in STUDIES) {
  dds_file <- cfg_for(s, "DDS")
  col_f    <- Sys.getenv(sprintf("CLEAN_COLS_%s", toupper(s)))
  strip    <- Sys.getenv(sprintf("CLEAN_STRIP_VERSION_%s", toupper(s))) == "1"
  e <- new.env()
  load(dds_file, envir = e)
  nm <- ls(e)[vapply(ls(e), function(x) is(get(x, envir = e), "DESeqDataSet"), TRUE)]
  if (!length(nm)) stop("no DESeqDataSet in ", dds_file, call. = FALSE)
  d <- get(nm[1], envir = e)
  m <- counts(d, normalized = FALSE)
  if (nzchar(col_f)) {
    keep <- grep(col_f, colnames(m))
    if (!length(keep)) stop("column filter '", col_f, "' matched nothing", call. = FALSE)
    m <- m[, keep, drop = FALSE]
  }
  g <- rownames(m)
  if (strip) {
    g <- strip_version(g)
    if (anyDuplicated(g))
      stop(s, ": stripping versions collides ", fmt_n(sum(duplicated(g))), " ids",
           call. = FALSE)
    rownames(m) <- g
  }
  raw[[s]] <- m; genes[[s]] <- g; samples[[s]] <- colnames(m)
  col <- cfg_for(s, "OG_SPECIES")
  mp <- og[Species == col & Gene %chin% g, .(Orthogroup, Gene)]
  map[[s]] <- unique(mp)
  say(sprintf("%-10s %s genes x %d libraries | %s genes in an orthogroup (%.1f%%) | %s groups",
              s, fmt_n(nrow(m)), ncol(m), fmt_n(uniqueN(mp$Gene)),
              100 * uniqueN(mp$Gene) / nrow(m), fmt_n(uniqueN(mp$Orthogroup))))
  unassigned <- setdiff(g, mp$Gene)
  say(sprintf("           %s genes have NO orthogroup and are dropped; they carry %.1f%% of the counts",
              fmt_n(length(unassigned)),
              100 * sum(m[unassigned, , drop = FALSE]) / sum(m)))
}

# --- the shared row set ------------------------------------------------------
per_species <- lapply(map, function(x) unique(x$Orthogroup))
shared <- if (BOTH_ONLY) Reduce(intersect, per_species) else Reduce(union, per_species)
say(sprintf("\northogroups with an expressed gene: %s",
            paste(sprintf("%s %s", names(per_species),
                          vapply(per_species, function(x) fmt_n(length(x)), "")),
                  collapse = " | ")))
say(sprintf("%s: %s groups",
            if (BOTH_ONLY) "present in BOTH species" else "union", fmt_n(length(shared))))

# --- aggregate, filter, transform -------------------------------------------
agg <- list()
for (s in STUDIES) {
  mp <- map[[s]][Orthogroup %chin% shared]
  idx <- match(mp$Gene, rownames(raw[[s]]))
  A <- rowsum(raw[[s]][idx, , drop = FALSE], group = mp$Orthogroup, reorder = TRUE)
  # Every group in the shared set must get a row, even in a species where it has
  # no expressed gene (possible only when BOTH_ONLY is off).
  missing <- setdiff(shared, rownames(A))
  if (length(missing)) {
    Z <- matrix(0, nrow = length(missing), ncol = ncol(A),
                dimnames = list(missing, colnames(A)))
    A <- rbind(A, Z)
  }
  A <- A[sort(shared), , drop = FALSE]
  agg[[s]] <- A
  say(sprintf("%-10s aggregated to %s x %d", s, fmt_n(nrow(A)), ncol(A)))
}

# THE CV FILTER, ON THE AGGREGATED COUNTS, then intersected. MIN_CV=15 is the
# pipeline's only gene filter (01_export_vst.r:76-88) and it is applied to RAW
# counts there, so it is applied to raw aggregated counts here -- but a group kept
# in one species and dropped in the other would break the shared row set, so the
# surviving sets are intersected and both costs are printed.
cv_keep <- list()
for (s in STUDIES) {
  A <- agg[[s]]
  mu <- rowMeans(A)
  sdv <- sqrt(rowSums((A - mu)^2) / (ncol(A) - 1))
  cv <- (sdv / (mu + 1e-6)) * 100
  cv_keep[[s]] <- rownames(A)[is.finite(cv) & cv >= MIN_CV]
  say(sprintf("%-10s CV >= %g keeps %s of %s groups (%.1f%%)",
              s, MIN_CV, fmt_n(length(cv_keep[[s]])), fmt_n(nrow(A)),
              100 * length(cv_keep[[s]]) / nrow(A)))
}
final <- Reduce(intersect, cv_keep)
say(sprintf("\nSHARED ROW SET after the CV filter: %s groups", fmt_n(length(final))))
for (s in STUDIES)
  say(sprintf("  %-10s loses %s groups it would have kept alone",
              s, fmt_n(length(setdiff(cv_keep[[s]], final)))))
if (!length(final)) stop("the shared row set is empty", call. = FALSE)
final <- sort(final)

for (s in STUDIES) {
  A <- agg[[s]][final, , drop = FALSE]
  storage.mode(A) <- "integer"
  d <- DESeqDataSetFromMatrix(A, colData = data.frame(row.names = colnames(A),
                                                      one = rep(1L, ncol(A))),
                              design = ~ 1)
  V <- assay(varianceStabilizingTransformation(d, blind = TRUE))
  if (!identical(rownames(V), final))
    stop("the transform reordered rows", call. = FALSE)
  if (anyNA(V)) stop("the transform produced NA", call. = FALSE)

  prefix <- cfg_for(s, "OG_PREFIX")
  ensure_dir(dirname(prefix))
  con <- file(paste0(prefix, ".f32"), "wb")       # row-major: the engine wants gene-major
  writeBin(as.vector(t(V)), con, size = 4L)
  close(con)
  writeLines(rownames(V), paste0(prefix, ".genes.txt"))
  write(toJSON(list(n_genes = nrow(V), n_samples = ncol(V), min_cv = MIN_CV,
                    unit = "orthogroup_2sp", aggregation = "sum_of_raw_counts",
                    transform = "varianceStabilizingTransformation(blind=TRUE)",
                    shared_row_set = TRUE, source_dds = cfg_for(s, "DDS"),
                    samples = colnames(V)), auto_unbox = TRUE, pretty = TRUE),
        paste0(prefix, ".meta.json"))
  say(sprintf("%-10s wrote %s.{f32,genes.txt,meta.json}  %s x %d",
              s, basename(prefix), fmt_n(nrow(V)), ncol(V)))
}

# --- the assertion the whole design exists for -------------------------------
a <- readLines(paste0(cfg_for(STUDIES[1], "OG_PREFIX"), ".genes.txt"))
for (s in STUDIES[-1]) {
  b <- readLines(paste0(cfg_for(s, "OG_PREFIX"), ".genes.txt"))
  if (!identical(a, b))
    stop("the two row sets differ -- edge conservation would need a projection ",
         "again, which is the one thing this stage exists to avoid", call. = FALSE)
}
say(sprintf("\nverified: both species carry the SAME %s rows, in the same order",
            fmt_n(length(a))))

write_tsv(data.table(
  metric = c("unit", "aggregation", "transform", "min_cv",
             "orthogroups_both_species", "orthogroups_after_cv",
             paste0("cv_kept_", STUDIES), paste0("lost_to_intersection_", STUDIES),
             paste0("libraries_", STUDIES)),
  value = c("orthogroup_2sp", "sum_of_raw_counts",
            "varianceStabilizingTransformation(blind=TRUE)", MIN_CV,
            fmt_n(length(shared)), fmt_n(length(final)),
            vapply(STUDIES, function(s) fmt_n(length(cv_keep[[s]])), ""),
            vapply(STUDIES, function(s) fmt_n(length(setdiff(cv_keep[[s]], final))), ""),
            vapply(STUDIES, function(s) as.character(ncol(agg[[s]])), ""))),
  file.path(dirname(cfg_for(STUDIES[1], "OG_PREFIX")), "og_matrix_summary.tsv"))

say("done")
