# =============================================================================
# 19_module_trait_spearman.r — module eigengene vs trait, by Spearman only
#
# Replaces the module-level use of 12_gene_trait_mi.py. The module response is
# now called by ONE statistic: Spearman's rho between a module's eigengene and
# the encoded trait, at padj <= CLEAN_MODULE_PADJ_THR and |rho| >= CLEAN_MODULE_R_THR.
#
# WHY SPEARMAN AND WHY ONLY SPEARMAN.
#   * The trait is ORDINAL, not interval. Purple's nitrogen is a dose (0/2/6 mM)
#     and Pearson reads those spacings literally -- it asks whether the response
#     from 2 to 6 mM is exactly twice the response from 0 to 2. Nothing in the
#     design justifies that. Spearman asks the question the experiment actually
#     poses: does the module move monotonically with nitrogen?
#   * Sugarcane's trait is two-level, where Spearman on the eigengene's ranks is
#     the rank-biserial correlation -- a Wilcoxon in correlation clothing, and
#     robust to the eigengene outliers a PC1 can carry.
#   * MI is dropped. It is an OMNIBUS test: it fires on any dependence at all,
#     including a dispersion change or a quirk in one library, and at the module
#     level that is what it was doing. A third of the previous `mi_only` modules
#     had two of 48 samples carrying over a quarter of the eigengene's variance.
#     It also needed an MI effect-size floor that could only be CALIBRATED
#     against the linear one, never derived, because the nats-to-|r| identity
#     assumes two continuous variables and this trait is discrete.
#   * Reporting three statistics also meant reporting three response classes
#     (`pearson_only` / `mi_only` / `both`), and every downstream figure, GO run
#     and TF test had to be read class by class on classes too small to separate.
#     One statistic, one responsive set.
#
# TIES. Both traits are heavily tied by design (sugarcane 24/24, purple 6/6/6),
# so the exact/AS-89 p-value is invalid and is not used. rho is Pearson on
# MIDRANKS and the p-value is the asymptotic t approximation on it -- exactly
# what `cor.test(method = "spearman", exact = FALSE)` computes. Verified against
# cor.test on this data at the end of the run; the check is fatal if it fails.
#
# The BH correction is over every module tested in the study, and `pval` is kept
# alongside `padj` so an uncorrected reading is available without a re-run.
#
# A LABEL-PERMUTATION NULL runs afterwards: the same eigengenes against shuffled
# trait labels, CLEAN_MODULE_PERM times, counting how many modules clear the same
# two thresholds. It is what lets the responsive count be reported as more than
# an artifact of testing 6,576 correlated summaries.
#
# Writes, replacing what the MI stage used to write to the same paths:
#   module_trait_<study>.tsv           module, trait, n, rho, pval, padj,
#                                      responsive, direction
#   module_trait_<study>.null.tsv      one row per permutation
#   module_trait_<study>.summary.json  counts and thresholds
#
# RUN: through run.sh  ->  ./run.sh moduletrait sugarcane
# =============================================================================

suppressMessages(library(data.table))

source(file.path(dirname(sub("--file=", "",
       grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "lib", "common.R"))

STUDY        <- env_req("CLEAN_STUDY")
EIG_PREFIX   <- env_req("CLEAN_EIGENGENE_PREFIX")
META_FILE    <- env_req("CLEAN_META")
TRAIT_SPEC   <- env_req("CLEAN_TRAITS")
SELECT_TRAIT <- env_opt("CLEAN_SELECT_TRAIT", "treatment")
OUT_FILE     <- env_req("CLEAN_OUT_FILE")
R_THR        <- env_num("CLEAN_MODULE_R_THR", 0.6)
PADJ_THR     <- env_num("CLEAN_MODULE_PADJ_THR", 0.05)
N_PERM       <- as.integer(env_num("CLEAN_MODULE_PERM", 1000))
SEED         <- as.integer(env_num("CLEAN_SEED", 1))
setDTthreads(as.integer(env_num("CLEAN_CORES", 100)))

banner(paste("module-trait Spearman:", STUDY))

# --- trait encoding ----------------------------------------------------------
# Same spec and same parser as 07_gene_trait_cor.r:
#   "genotype:A=1,B=0;treatment:High Nitrogen=1,Low Nitrogen=0"
# Only the ordering of the encoded values matters to Spearman, not the spacing,
# which is the point -- purple's 0/2/6 and a 0/1/2 coding give identical rho.
parse_traits <- function(spec) {
  out <- list()
  for (block in strsplit(spec, ";", fixed = TRUE)[[1L]]) {
    block <- trimws(block)
    if (!nzchar(block)) next
    nm   <- sub(":.*$", "", block)
    body <- sub("^[^:]*:", "", block)
    kv   <- strsplit(strsplit(body, ",", fixed = TRUE)[[1L]], "=", fixed = TRUE)
    if (any(lengths(kv) != 2L))
      stop("cannot parse trait spec '", block, "' — expected LEVEL=value pairs",
           call. = FALSE)
    v <- as.numeric(vapply(kv, `[`, "", 2L))
    names(v) <- trimws(vapply(kv, `[`, "", 1L))
    out[[trimws(nm)]] <- v
  }
  out
}
TRAITS <- parse_traits(TRAIT_SPEC)
if (!SELECT_TRAIT %in% names(TRAITS))
  stop("CLEAN_SELECT_TRAIT '", SELECT_TRAIT, "' is not one of the encoded ",
       "traits (", paste(names(TRAITS), collapse = ", "), ")", call. = FALSE)
enc <- TRAITS[[SELECT_TRAIT]]
say("trait: ", SELECT_TRAIT, "  ",
    paste(sprintf("%s=%g", names(enc), enc), collapse = ", "))

# --- eigengenes --------------------------------------------------------------
# 14_module_eigengene.r writes the eigengene matrix in the VST export's own
# format, so read_vst() reads it unchanged; the rows are modules, not genes.
E <- read_vst(EIG_PREFIX)
say("modules with an eigengene: ", fmt_n(nrow(E)), " x ", ncol(E), " samples")

meta <- fread(META_FILE, header = TRUE)
setnames(meta, tolower(names(meta)))
if (!"sample" %in% names(meta))
  stop(basename(META_FILE), " has no `sample` column", call. = FALSE)
m <- meta[match(colnames(E), sample)]
if (anyNA(m$sample))
  stop("these eigengene samples are missing from ", basename(META_FILE), ": ",
       paste(head(colnames(E)[is.na(m$sample)], 5), collapse = ", "), call. = FALSE)
if (!SELECT_TRAIT %in% names(m))
  stop("trait column '", SELECT_TRAIT, "' is not in ", basename(META_FILE),
       "\n  available: ", paste(names(m), collapse = ", "), call. = FALSE)

trait_vec <- unname(enc[as.character(m[[SELECT_TRAIT]])])
ok <- !is.na(trait_vec)
if (sum(!ok))
  say("dropping ", sum(!ok), " samples with unencoded level(s): ",
      paste(unique(as.character(m[[SELECT_TRAIT]])[!ok]), collapse = ", "))
n <- sum(ok)
if (n < 4L) stop("trait '", SELECT_TRAIT, "' has only ", n, " usable samples",
                 call. = FALSE)
df <- n - 2L
say(sprintf("usable samples: n = %d, df = %d,  level counts: %s", n, df,
            paste(sprintf("%g:%d", sort(unique(trait_vec[ok])),
                          tabulate(match(trait_vec[ok], sort(unique(trait_vec[ok])))) ),
                  collapse = "  ")))

# --- Spearman ----------------------------------------------------------------
# Rank once, then one vectorised cor() over the whole modules x samples matrix.
# Midranks on both sides; `rank()`'s default ties.method = "average" is what
# Spearman requires.
Eok  <- E[, ok, drop = FALSE]
Rk   <- t(apply(Eok, 1L, rank))            # modules x samples, midranks
ty   <- rank(trait_vec[ok])

spear <- function(Rk, ty) {
  r <- suppressWarnings(as.numeric(cor(t(Rk), ty)))
  r[!is.finite(r)] <- 0                    # a constant eigengene has no rank order
  r
}
rho <- spear(Rk, ty)

t_stat <- rho * sqrt(df / (1 - rho^2 + 1e-15))
pval   <- 2 * pt(-abs(t_stat), df = df)
padj   <- p.adjust(pval, method = "BH")

res <- data.table(
  module     = rownames(E),
  trait      = SELECT_TRAIT,
  n          = n,
  rho        = round(rho, 4),
  pval       = signif(pval, 4),
  padj       = signif(padj, 4))
res[, responsive := padj <= PADJ_THR & abs(rho) >= R_THR]
res[, direction  := fifelse(!responsive, "none",
                    fifelse(rho > 0, "positive", "negative"))]

# --- the check that this really is cor.test ----------------------------------
# Cheap, and it is the whole basis for not writing a p-value by hand: take a
# spread of modules and confirm rho and pval match R's own implementation.
idx <- unique(round(seq(1, nrow(E), length.out = min(25L, nrow(E)))))
mx_r <- 0; mx_p <- 0
for (i in idx) {
  ct <- suppressWarnings(cor.test(Eok[i, ], trait_vec[ok],
                                  method = "spearman", exact = FALSE))
  mx_r <- max(mx_r, abs(unname(ct$estimate) - rho[i]))
  mx_p <- max(mx_p, abs(ct$p.value - pval[i]))
}
say(sprintf("cor.test agreement over %d modules: max |drho| = %.2e, max |dp| = %.2e",
            length(idx), mx_r, mx_p))
if (mx_r > 1e-8 || mx_p > 1e-8)
  stop("this stage does not reproduce cor.test(method='spearman', exact=FALSE)",
       call. = FALSE)

# --- label-permutation null --------------------------------------------------
# The same eigengenes against shuffled trait labels. The eigengenes stay exactly
# as correlated with each other as they really are, so this is a null for the
# TRAIT association and not for the module structure.
set.seed(SEED)
null_rows <- rbindlist(lapply(seq_len(N_PERM), function(p) {
  r  <- spear(Rk, sample(ty))
  ts <- r * sqrt(df / (1 - r^2 + 1e-15))
  pa <- p.adjust(2 * pt(-abs(ts), df = df), method = "BH")
  data.table(perm = p, max_abs_rho = max(abs(r)),
             n_padj = sum(pa <= PADJ_THR),
             n_responsive = sum(pa <= PADJ_THR & abs(r) >= R_THR))
}))
write_tsv(null_rows, sub("\\.tsv$", ".null.tsv", OUT_FILE))
say(sprintf("null over %s permutations: responsive modules mean %.2f, max %d; max |rho| %.3f",
            fmt_n(N_PERM), mean(null_rows$n_responsive), max(null_rows$n_responsive),
            max(null_rows$max_abs_rho)))

# --- out ---------------------------------------------------------------------
res <- res[order(padj, -abs(rho))]
write_tsv(res, OUT_FILE)

n_resp <- res[responsive == TRUE, .N]
say("")
say(sprintf("modules tested:            %s", fmt_n(nrow(res))))
say(sprintf("  padj <= %.2f:              %s", PADJ_THR, fmt_n(res[padj <= PADJ_THR, .N])))
say(sprintf("  |rho| >= %.2f:             %s", R_THR, fmt_n(res[abs(rho) >= R_THR, .N])))
say(sprintf("  RESPONSIVE (both):       %s  (%s positive, %s negative)",
            fmt_n(n_resp), fmt_n(res[direction == "positive", .N]),
            fmt_n(res[direction == "negative", .N])))
say(sprintf("  for reference, raw p <= %.2f and |rho| >= %.2f: %s (uncorrected)",
            PADJ_THR, R_THR, fmt_n(res[pval <= PADJ_THR & abs(rho) >= R_THR, .N])))

jf <- sub("\\.tsv$", ".summary.json", OUT_FILE)
writeLines(sprintf(paste0(
  '{\n  "study": "%s",\n  "statistic": "spearman",\n  "trait": "%s",\n',
  '  "n_modules": %d,\n  "n_samples": %d,\n  "trait_levels": [%s],\n',
  '  "rho_threshold": %g,\n  "padj_threshold": %g,\n',
  '  "n_padj": %d,\n  "n_rho": %d,\n  "n_responsive": %d,\n',
  '  "n_positive": %d,\n  "n_negative": %d,\n',
  '  "n_perm": %d,\n  "null_responsive_mean": %.4f,\n  "null_responsive_max": %d,\n',
  '  "null_max_abs_rho": %.4f,\n  "out_file": "%s"\n}'),
  STUDY, SELECT_TRAIT, nrow(res), n, paste(sort(unique(trait_vec[ok])), collapse = ", "),
  R_THR, PADJ_THR, res[padj <= PADJ_THR, .N], res[abs(rho) >= R_THR, .N], n_resp,
  res[direction == "positive", .N], res[direction == "negative", .N],
  N_PERM, mean(null_rows$n_responsive), max(null_rows$n_responsive),
  max(null_rows$max_abs_rho), OUT_FILE), jf)
say("wrote ", basename(jf))
say("done: ", STUDY)
