# =============================================================================
# lib/common.R — helpers shared by the R stages
#
# Every stage takes its configuration from the environment, which run.sh fills
# in from config.sh. Nothing is hardcoded and nothing is edited between runs.
# =============================================================================

# --- configuration -----------------------------------------------------------

env_req <- function(name) {
  v <- Sys.getenv(name)
  if (!nzchar(v))
    stop("required environment variable ", name, " is unset. ",
         "This stage is meant to be launched through new_clean/run.sh.",
         call. = FALSE)
  v
}

env_opt <- function(name, default = "") {
  v <- Sys.getenv(name)
  if (nzchar(v)) v else default
}

env_num  <- function(name, default = NA_real_) {
  v <- Sys.getenv(name)
  if (nzchar(v)) as.numeric(v) else default
}

env_flag <- function(name, default = FALSE) {
  v <- Sys.getenv(name)
  if (!nzchar(v)) return(default)
  v %in% c("1", "TRUE", "true", "yes", "Y")
}

env_list <- function(name, default = character()) {
  v <- Sys.getenv(name)
  if (!nzchar(v)) return(default)
  strsplit(trimws(v), "[[:space:]]+")[[1L]]
}

# --- logging -----------------------------------------------------------------

say <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

banner <- function(x) {
  cat("\n", strrep("=", 70), "\n", x, "\n", strrep("=", 70), "\n", sep = "")
}

fmt_n <- function(x) format(x, big.mark = ",", scientific = FALSE)

# --- gene ids ----------------------------------------------------------------

# Ids are normalised once, in 01_export_vst.r, so everything under
# new_clean/results/ already carries the bare form. This exists for joining
# against files OUTSIDE that tree -- proteomes, orthogroups, TF tables, the
# Munoz supplementary sheets -- which still carry the annotation version.
# Stripping an id that has no suffix is a no-op, so it is always safe to apply.
strip_version <- function(x) sub("\\.v[0-9]+(\\.[0-9]+)*$", "", x)

# --- edge files --------------------------------------------------------------

# The network schema, as written by 03_merge_layers.py. Asserting it here means
# a downstream stage fails loudly on a schema change rather than reading the
# wrong column -- the failure mode that produced the old tree's mixed layers.
NETWORK_COLS <- c("gene1", "gene2", "stat", "pval", "padj",
                  "weight", "pearson_r", "ksg", "source")

assert_network_schema <- function(path, need = NETWORK_COLS) {
  if (!file.exists(path))
    stop("network edge file not found: ", path,
         "\n  build it first with:  ./run.sh build <study>", call. = FALSE)
  hdr <- names(data.table::fread(path, nrows = 0L))
  missing <- setdiff(need, hdr)
  if (length(missing))
    stop("edge file ", basename(path), " is missing column(s): ",
         paste(missing, collapse = ", "),
         "\n  found: ", paste(hdr, collapse = ", "),
         "\n  expected the 03_merge_layers.py schema.", call. = FALSE)
  # Positional readers downstream rely on gene1/gene2 being first.
  if (!identical(hdr[1:2], c("gene1", "gene2")))
    stop("edge file ", basename(path), " does not start with gene1, gene2 ",
         "(found ", paste(hdr[1:2], collapse = ", "), "). Some stages read ",
         "those positionally.", call. = FALSE)
  invisible(hdr)
}

# --- orthogroups -------------------------------------------------------------

# Parse an OrthoFinder Orthogroups.tsv into long form: Orthogroup, Species, Gene.
#
# This replaces cogeqc::read_orthogroups(). cogeqc is installed in three conda
# envs on this machine and none of them is the one the network stages run in, so
# depending on it meant either a fourth interpreter in config.sh or a stage that
# dies on `there is no package called 'cogeqc'` -- which is exactly what
# happened. The file is a header row of species names and one comma-separated
# gene list per cell; that is not worth a cross-env dependency.
#
# Verified against cogeqc::read_orthogroups on the project's own Orthogroups.tsv:
# 360,794 rows, same columns in the same order, and Orthogroup and Gene identical
# element for element. The one difference is that cogeqc returns Species as a
# factor and this returns character. That is irrelevant to both callers -- they
# do `Species == x` then `Species := y`, and data.table adds the new factor level
# transparently, so the two behave the same -- and character is the safer of the
# two for a column that gets relabelled.
read_orthogroups_tsv <- function(path) {
  if (!file.exists(path))
    stop("Orthogroups.tsv not found: ", path, call. = FALSE)
  wide <- data.table::fread(path, header = TRUE, sep = "\t",
                            colClasses = "character", na.strings = NULL)
  species <- setdiff(names(wide), "Orthogroup")
  long <- data.table::rbindlist(lapply(species, function(sp) {
    v <- wide[[sp]]
    keep <- !is.na(v) & nzchar(v)
    if (!any(keep)) return(NULL)
    genes <- strsplit(v[keep], ", ", fixed = TRUE)
    data.table::data.table(
      Orthogroup = rep(wide$Orthogroup[keep], lengths(genes)),
      Species    = sp,
      Gene       = unlist(genes, use.names = FALSE))
  }))
  as.data.frame(long)
}

# --- the VST export ----------------------------------------------------------

# Read the <prefix>.{f32,genes.txt,meta.json} triple written by 01_export_vst.r.
#
# Every stage that needs expression reads THIS, not the dds. That is what makes
# the provenance mismatch structurally impossible: in the old tree the network
# came from one quantification and the gene-trait correlation from another, and
# nothing connected the two well enough to notice. Here there is one matrix per
# study and both the network and the trait correlation are computed from it.
read_vst <- function(prefix) {
  meta_f  <- paste0(prefix, ".meta.json")
  genes_f <- paste0(prefix, ".genes.txt")
  bin_f   <- paste0(prefix, ".f32")
  for (f in c(meta_f, genes_f, bin_f))
    if (!file.exists(f))
      stop("missing ", basename(f), " — run  ./run.sh export <study>  first",
           call. = FALSE)

  txt <- paste(readLines(meta_f, warn = FALSE), collapse = " ")
  grab_num <- function(k)
    as.integer(sub(paste0('.*"', k, '"[[:space:]]*:[[:space:]]*([0-9]+).*'), "\\1", txt))
  n_genes   <- grab_num("n_genes")
  n_samples <- grab_num("n_samples")
  samples <- regmatches(sub('.*"samples"[[:space:]]*:[[:space:]]*\\[([^]]*)\\].*', "\\1", txt),
                        gregexpr('"[^"]*"',
                                 sub('.*"samples"[[:space:]]*:[[:space:]]*\\[([^]]*)\\].*', "\\1", txt)))[[1L]]
  samples <- gsub('"', "", samples)

  genes <- readLines(genes_f, warn = FALSE)
  if (length(genes) != n_genes)
    stop("genes.txt has ", length(genes), " lines, meta says ", n_genes,
         call. = FALSE)
  if (length(samples) != n_samples)
    stop("meta.json lists ", length(samples), " sample names for ", n_samples,
         " samples", call. = FALSE)

  # written row-major (gene-major); R fills column-major, so read transposed
  v <- readBin(bin_f, "numeric", n = n_genes * n_samples, size = 4L)
  m <- matrix(v, nrow = n_samples, ncol = n_genes)
  m <- t(m)
  dimnames(m) <- list(genes, samples)
  m
}

# --- output ------------------------------------------------------------------

ensure_dir <- function(path) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  path
}

write_tsv <- function(x, path) {
  ensure_dir(dirname(path))
  data.table::fwrite(x, path, sep = "\t", quote = FALSE)
  say("wrote ", basename(path), "  (", fmt_n(nrow(x)), " rows)")
  invisible(path)
}

# =============================================================================
# BLOCKED (DESIGN-AWARE) CORRELATION
#
# Shared by the gene level (31_gene_trait_blocked.r, branch blocked-gene-trait)
# and the module level (53_module_trait_blocked.r). One implementation, one
# verification, so the two levels cannot drift apart.
#
# The statistic is the t on the trait coefficient of  y ~ block + trait, which is
# exactly the partial correlation's t-test, so  r = t / sqrt(t^2 + df)  is the
# partial correlation of y and the trait GIVEN the blocks -- the same quantity
# the |r| floor was always cutting on, now conditional on the design.
#
# Vectorised over ROWS via one QR: a per-row loop does not finish over 44,118
# genes, and the module level has 7,493 rows.
# =============================================================================

fit_blocked <- function(Y, y, blocks) {
  # Y: rows (genes or modules) x samples. y: numeric trait. blocks: list of factors.
  n <- ncol(Y)
  X <- if (length(blocks))
         model.matrix(~ ., data = cbind(data.frame(blocks), .trait = y))
       else model.matrix(~ y)
  j   <- ncol(X)                                   # trait is the last column
  qrX <- qr(X)
  k   <- qrX$rank
  df  <- n - k
  if (df <= 0) stop("no residual degrees of freedom", call. = FALSE)
  # A rank-deficient model matrix is what a block CONFOUNDED WITH THE TRAIT looks
  # like -- e.g. a genotype that only ever appears at one nitrogen level. Fitting
  # it anyway would silently attribute the trait effect to the block.
  if (k < ncol(X))
    stop("model matrix is rank-deficient (", k, " < ", ncol(X),
         ") -- a block is confounded with the trait", call. = FALSE)

  B    <- qr.coef(qrX, t(Y))                       # coefficients x rows
  RES  <- t(Y) - X %*% B
  s2   <- colSums(RES^2) / df
  XtXi <- chol2inv(qr.R(qrX))
  est  <- B[j, ]
  se   <- sqrt(s2 * XtXi[j, j])
  tt   <- as.numeric(est / se)
  list(est = as.numeric(est), se = as.numeric(se), t = tt,
       r = tt / sqrt(tt^2 + df), p = 2 * pt(-abs(tt), df = df), df = df)
}

# Refit a spread of rows with lm() and ABORT on disagreement. The vectorised
# solver is fast because it skips every check lm() makes; this is what licenses
# trusting it.
verify_against_lm <- function(Y, y, blocks, fit, label, n_check = 25L) {
  idx  <- unique(round(seq(1, nrow(Y), length.out = min(n_check, nrow(Y)))))
  dd   <- if (length(blocks)) cbind(as.data.frame(blocks), .trait = y)
          else data.frame(.trait = y)
  mx_t <- 0; mx_p <- 0
  for (i in idx) {
    co <- summary(stats::lm(as.numeric(Y[i, ]) ~ ., data = dd))$coefficients
    r  <- co[".trait", ]
    mx_t <- max(mx_t, abs(r["t value"]  - fit$t[i]))
    mx_p <- max(mx_p, abs(r["Pr(>|t|)"] - fit$p[i]))
  }
  say(sprintf("lm() agreement (%s) over %d rows: max |dt| = %.2e, max |dp| = %.2e",
              label, length(idx), mx_t, mx_p))
  if (mx_t > 1e-8 || mx_p > 1e-8)
    stop("the vectorised solver does not reproduce lm() for ", label, call. = FALSE)
  invisible(TRUE)
}

# Midranks, row-wise. rank()'s default ties.method = "average" is what
# 19_module_trait_spearman.r requires: both traits are heavily tied by design.
row_midranks <- function(Y) t(apply(Y, 1L, rank))
