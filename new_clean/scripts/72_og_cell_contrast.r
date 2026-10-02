#!/usr/bin/env Rscript
# =============================================================================
# 72_og_cell_contrast.r -- what kind of gene family lets its copies drift apart,
# and what kind holds them together.
#
# WHAT THE TWO SETS ARE. 67 sorted every multi-copy orthogroup into a 2x2: do its
# copies share an expression profile (uniform / divergent), and can the quantifier
# tell them apart (separable / inseparable). Two cells carry the biology:
#
#   divergent/separable  13,902 sugarcane / 16,271 purple -- copies drifted apart,
#                        and 65's read-stealing diagnostic says that is real
#   uniform/separable    11,740 / 8,742 -- copies still move in lockstep
#
# This stage describes them. It does not test function -- that is 73 -- it prepares
# the sets, measures every covariate already on disk, and settles the two things
# that would otherwise make 73's enrichment uninterpretable.
#
# CONFOUND ONE: COPY NUMBER, AND IT RUNS THE WRONG WAY. Uniform groups have MORE
# copies (mean 3.55 vs 2.82 in sugarcane, Mann-Whitney p = 6e-236). That is almost
# certainly a power effect rather than biology: with more members, concordance is
# easier to prove above the size-matched null, so "uniform" partly means "enough
# copies to show it". Every test here is therefore run TWICE -- once on the raw
# cells and once on a 1:1 EXACT-COPY-NUMBER-MATCHED pair of sets -- and both are
# written, so a reader sees which conclusions survive the matching.
#
# CONFOUND TWO: ANNOTATION COVERAGE, AND MATCHING DOES NOT FIX IT. GO coverage is
# 57.5% for divergent against 67.5% for uniform (sugarcane), and the gap persists
# INSIDE every copy-number band -- ~10 points at n = 2, 3, 4-5 and 6-20 alike.
# InterPro shows the same gap at a higher baseline (79.1% vs 88.2%). So it is a
# property of the sets, not of the annotation source, and it cannot be matched
# away. It is reported as a finding, and 73's matching happens WITHIN the annotated
# orthogroups so the enrichment at least compares like with like.
#
# WHAT THIS STAGE WRITES FOR 73. There is no orthogroup-level GO table anywhere in
# this tree -- gene2go is gene-keyed. One is built here, with the column named
# `gene` even though it holds an orthogroup id, because that is what makes it drop
# into parse_gene2go(), 63's reader and annFUN.gene2GO with no code change, and it
# is what results_og/<study>/network_<study>_node_metrics.tsv already does.
#
# AN ORTHOGROUP INHERITS THE UNION of its members' terms. The alternative --
# requiring a fraction of copies to carry a term -- would discard exactly the
# signal this analysis is about, since a term held by one copy and not its sibling
# IS copy divergence. The union's hazard is size: OG0000000 has 1,233 purple
# copies and would accumulate the union of 1,233 genes' annotations. The classified
# set is capped at 20 members by 67, so the two cells are safe; the giant groups
# only reach 73 through the all-orthogroup background, which caps them there.
#
# RUN: through run.sh -> ./run.sh ogcellcontrast <study>
# =============================================================================

suppressMessages(library(data.table))
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE),
                                                 value = TRUE)[1])), "lib", "common.R"))

STUDY      <- env_req("CLEAN_STUDY")
UNIF_FILE  <- env_req("CLEAN_UNIFORMITY")
OTHER_UNIF <- env_req("CLEAN_UNIFORMITY_OTHER")
OTHER_NAME <- env_req("CLEAN_OTHER_STUDY")
ORTHO_FILE <- env_req("CLEAN_ORTHOGROUPS")
OG_COL     <- env_req("CLEAN_OG_SPECIES")
GENE2GO    <- env_req("CLEAN_GENE2GO")
IPS_FILE   <- env_req("CLEAN_IPS_TSV")
CROSSWALK  <- env_req("CLEAN_CROSSWALK")
TF_FILE    <- env_req("CLEAN_TF_FILE")
OG2GO_OUT  <- env_req("CLEAN_OG2GO_OUT")
OUT_DIR    <- env_req("CLEAN_OUT_DIR")
OG_NODES   <- env_opt("CLEAN_OG_NODE_METRICS")
MAX_COPIES <- as.integer(env_num("CLEAN_MAX_OG_COPIES", 20))
NULL_REPS  <- as.integer(env_num("CLEAN_NULL_REPS", 200))
SEED       <- as.integer(env_num("CLEAN_SEED", 1))

CELL_DIV <- "divergent/separable"
CELL_UNI <- "uniform/separable"

banner(paste("orthogroup cell contrast:", STUDY))
ensure_dir(OUT_DIR)
set.seed(SEED)

# --- the orthogroup partition ------------------------------------------------
og <- as.data.table(read_orthogroups_tsv(ORTHO_FILE))
og[, Gene := strip_version(sub("\\.p[0-9]+$", "", Gene))]
og <- unique(og[Species == OG_COL, .(og_2sp = Orthogroup, gene = Gene)])
say(sprintf("%s orthogroups carry %s %s genes",
            fmt_n(uniqueN(og$og_2sp)), fmt_n(nrow(og)), STUDY))

# --- the orthogroup-level GO table, which does not exist until now -----------
g2g <- fread(GENE2GO, sep = "\t", header = TRUE,
             select = c("gene", "go_id", "source"), colClasses = "character")
g2g <- g2g[grepl("^GO:", go_id)]
j <- merge(og, g2g, by = "gene", allow.cartesian = TRUE)
og2go <- j[, .(source = "og_union", n_members_with_term = uniqueN(gene)),
           by = .(gene = og_2sp, go_id)]
setcolorder(og2go, c("gene", "go_id", "source", "n_members_with_term"))
if (anyDuplicated(og2go, by = c("gene", "go_id")))
  stop("og2go has duplicate (orthogroup, term) rows", call. = FALSE)
ensure_dir(dirname(OG2GO_OUT))
write_tsv(og2go[order(gene, go_id)], OG2GO_OUT)
say(sprintf("og2go: %s orthogroups, %s distinct terms, %s rows",
            fmt_n(uniqueN(og2go$gene)), fmt_n(uniqueN(og2go$go_id)), fmt_n(nrow(og2go))))

# --- per-orthogroup annotation presence --------------------------------------
# InterPro streamed, not read whole: the table is 338 MB (sugarcane) / 421 MB.
say("streaming the InterProScan table for IPR accessions")
ipr_cmd <- sprintf("awk -F'\\t' 'NF>=13 && $12!=\"-\" && $12!=\"\" {print $1\"\\t\"$12}' %s | sort -u",
                   shQuote(IPS_FILE))
ipr <- fread(cmd = ipr_cmd, header = FALSE, col.names = c("gene", "ipr"),
             colClasses = "character")
og_ipr <- unique(merge(og, ipr, by = "gene", allow.cartesian = TRUE)[, .(og_2sp, ipr)])
say(sprintf("IPR: %s genes, %s orthogroups, %s distinct accessions",
            fmt_n(uniqueN(ipr$gene)), fmt_n(uniqueN(og_ipr$og_2sp)),
            fmt_n(uniqueN(og_ipr$ipr))))

tf <- fread(TF_FILE, sep = "\t", header = TRUE, colClasses = "character")
tf <- unique(tf[, .(gene = strip_version(gene), Family)])     # isoforms hit several families
og_tf <- unique(merge(og, tf, by = "gene")[, .(og_2sp, Family)])

ann <- og[, .(n_genes = .N), by = og_2sp]
ann[, has_go  := og_2sp %chin% unique(og2go$gene)]
ann[, has_ipr := og_2sp %chin% unique(og_ipr$og_2sp)]
ann[, has_tf  := og_2sp %chin% unique(og_tf$og_2sp)]
ann <- merge(ann, og2go[, .(n_go_terms = .N), by = .(og_2sp = gene)], by = "og_2sp", all.x = TRUE)
ann <- merge(ann, og_ipr[, .(n_ipr = .N), by = og_2sp], by = "og_2sp", all.x = TRUE)
ann <- merge(ann, og_tf[, .(tf_families = paste(sort(unique(Family)), collapse = ";")),
                        by = og_2sp], by = "og_2sp", all.x = TRUE)
for (cl in c("n_go_terms", "n_ipr")) ann[is.na(get(cl)), (cl) := 0L]
ann[is.na(tf_families), tf_families := ""]
write_tsv(ann[order(og_2sp)], file.path(OUT_DIR, sprintf("og_annotation_%s.tsv", STUDY)))

# --- assemble the two cells with every covariate on disk ---------------------
u <- fread(UNIF_FILE)
# `cell` reads back as "" for unclassified rows, NOT NA. Filtering on !is.na()
# silently admits ~49,000 single-copy groups; this is the filter that does not.
d <- u[cell %in% c(CELL_DIV, CELL_UNI)]
d[, divergent := cell == CELL_DIV]
say(sprintf("\n%s: %s divergent/separable | %s uniform/separable",
            STUDY, fmt_n(sum(d$divergent)), fmt_n(sum(!d$divergent))))

cw <- fread(CROSSWALK)
d <- merge(d, cw[, .(og_2sp,
                     omega = as.numeric(get(sprintf("omega_median_%s", STUDY))),
                     frac_unique = as.numeric(get(sprintf("frac_unique_min_%s", STUDY))),
                     dup_class = get(sprintf("class_%s", STUDY)))],
           by = "og_2sp", all.x = TRUE)
d <- merge(d, ann[, .(og_2sp, has_go, has_ipr, has_tf, n_go_terms, n_ipr)],
           by = "og_2sp", all.x = TRUE)
if (nzchar(OG_NODES) && file.exists(OG_NODES)) {
  # First column is named `gene` but holds orthogroup ids, as every results_og
  # table does. This is the ORTHOGROUP network's degree -- og_uniformity's
  # degree_median is the GENE-level one and the two must not be conflated.
  ogn <- fread(OG_NODES, select = c("gene", "degree"))
  setnames(ogn, c("og_2sp", "og_degree"))
  d <- merge(d, ogn, by = "og_2sp", all.x = TRUE)
}

# --- 1:1 exact-copy-number matching ------------------------------------------
# Exact n, not bands: banding leaves a residual gradient inside each band, and the
# whole point is that the matched sets be indistinguishable on copy number.
match_pairs <- function(dt) {
  out <- vector("list", 0L)
  for (k in sort(unique(dt$n_members))) {
    a <- dt[divergent == TRUE  & n_members == k, og_2sp]
    b <- dt[divergent == FALSE & n_members == k, og_2sp]
    m <- min(length(a), length(b))
    if (m == 0L) next
    out[[length(out) + 1L]] <- data.table(
      og_2sp = c(sample(a, m), sample(b, m)),
      divergent = rep(c(TRUE, FALSE), each = m), n_members = k)
  }
  rbindlist(out)
}
mp <- match_pairs(d)
d[, matched := og_2sp %chin% mp$og_2sp]
say(sprintf("matched 1:1 on exact copy number: %s pairs (%s of %s divergent, %s of %s uniform)",
            fmt_n(sum(mp$divergent)), fmt_n(d[divergent & matched, .N]),
            fmt_n(d[(divergent), .N]), fmt_n(d[!divergent & matched, .N]),
            fmt_n(d[!(divergent), .N])))

# And the same, among GO-annotated orthogroups only -- the set 73 must use, since
# the annotation gap survives copy-number matching and topGO conditions on being
# annotated anyway.
mg <- match_pairs(d[has_go == TRUE])
d[, matched_go := og_2sp %chin% mg$og_2sp]
mi <- match_pairs(d[has_ipr == TRUE])
d[, matched_ipr := og_2sp %chin% mi$og_2sp]
say(sprintf("  within GO-annotated: %s pairs | within IPR-annotated: %s pairs",
            fmt_n(sum(mg$divergent)), fmt_n(sum(mi$divergent))))

# AND THE SAME AGAIN, HOMEOLOG CLASS ONLY -- the control that decides whether a
# domain is a POLYPLOIDY result. `homeolog` means one chromosome, several
# haplotypes, i.e. copies the recent hybrid polyploidy retained; `dispersed` and
# `unplaced` include transposon families, whose many near-identical scattered
# members OrthoFinder groups into one orthogroup. Those are 2.6% of purple's cells
# and 0.8% of sugarcane's, but 78% and 73% of them land in the divergent cell,
# while only 11.5% and 17.3% are homeologs against 47.1% and 31.3% of the rest.
# Measured 2026-10-02: restricting to homeologs removes EVERY transposase,
# reverse-transcriptase and RNase-H domain from the divergent enrichment, and
# removes F-box as well -- so neither was a polyploidy finding. PPR, TPR and the E
# motif survive in both species.
mgh <- match_pairs(d[has_go == TRUE & dup_class == "homeolog"])
mih <- match_pairs(d[has_ipr == TRUE & dup_class == "homeolog"])
d[, matched_go_homeolog  := og_2sp %chin% mgh$og_2sp]
d[, matched_ipr_homeolog := og_2sp %chin% mih$og_2sp]
say(sprintf("  homeolog only: %s GO pairs | %s IPR pairs",
            fmt_n(sum(mgh$divergent)), fmt_n(sum(mih$divergent))))

# --- the test accumulator ----------------------------------------------------
RES <- list()
add <- function(test, statistic, value, n, note = "") {
  RES[[length(RES) + 1L]] <<- data.table(
    study = STUDY, test = test, statistic = statistic,
    value = if (is.numeric(value)) value else NA_real_, n = n, note = note)
}

num_test <- function(dt, col, label, scope) {
  a <- suppressWarnings(as.numeric(dt[divergent == TRUE][[col]]))
  b <- suppressWarnings(as.numeric(dt[divergent == FALSE][[col]]))
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (length(a) < 30L || length(b) < 30L) {
    add(sprintf("%s: %s", scope, label), "skipped", NA, length(a) + length(b),
        sprintf("too few finite values (%d divergent / %d uniform)", length(a), length(b)))
    return(invisible())
  }
  w <- wilcox.test(a, b)
  add(sprintf("%s: %s", scope, label), "wilcox_p", w$p.value, length(a) + length(b),
      sprintf("median %.4f divergent (n=%d) vs %.4f uniform (n=%d)",
              median(a), length(a), median(b), length(b)))
}

cat_test <- function(dt, col, label, scope) {
  x <- dt[!is.na(get(col))]
  if (!nrow(x)) return(invisible())
  tb <- table(x$divergent, as.logical(x[[col]]))
  if (nrow(tb) != 2L || ncol(tb) != 2L) return(invisible())
  ft <- fisher.test(tb)
  add(sprintf("%s: %s", scope, label), "fisher_p", ft$p.value, nrow(x),
      sprintf("%.2f%% divergent vs %.2f%% uniform | OR %.3f",
              100 * mean(as.logical(x[divergent == TRUE][[col]])),
              100 * mean(as.logical(x[divergent == FALSE][[col]])), ft$estimate))
}

NUMS <- list(n_members = "copies per group", omega = "omega vs sorghum",
             expr_concordance = "expression concordance",
             sign_concordance = "N-response sign concordance",
             degree_median = "gene-level degree (median)",
             infrv_member_mean = "member InfRV", infrv_gain = "infrv_gain",
             frac_unique = "frac_unique_min", n_go_terms = "GO terms per group",
             n_ipr = "IPR accessions per group")
if ("og_degree" %in% names(d)) NUMS$og_degree <- "orthogroup-network degree"

for (scope in c("1 unmatched", "2 copy-number matched")) {
  dt <- if (scope == "1 unmatched") d else d[matched == TRUE]
  for (cl in names(NUMS)) num_test(dt, cl, NUMS[[cl]], scope)
  dt2 <- copy(dt)
  dt2[, any_responsive := n_responsive > 0]
  cat_test(dt2, "any_responsive", "has an N-responsive member", scope)
  cat_test(dt2, "has_tf", "contains a TF", scope)
  cat_test(dt2, "has_go", "has any GO term", scope)
  cat_test(dt2, "has_ipr", "has any IPR accession", scope)
  hom <- dt[dup_class %in% c("homeolog", "dispersed")]
  if (nrow(hom) > 60L) {
    hom[, is_hom := dup_class == "homeolog"]
    cat_test(hom, "is_hom", "duplication class is homeolog", scope)
  }
}

# --- the matching has to have worked, or everything above it is decoration ---
md <- d[matched == TRUE]
w_cop <- wilcox.test(md[divergent == TRUE, n_members], md[divergent == FALSE, n_members])
add("3 matching check: copies after matching", "wilcox_p", w_cop$p.value, nrow(md),
    sprintf("mean %.3f vs %.3f", mean(md[divergent == TRUE, n_members]),
            mean(md[divergent == FALSE, n_members])))
if (w_cop$p.value < 0.5)
  stop("after matching, copy number still differs (p = ", signif(w_cop$p.value, 3),
       ") -- the matching did not work and every matched test below it is ",
       "uninterpretable", call. = FALSE)
say(sprintf("matching check: copies indistinguishable after matching (p = %.3f)",
            w_cop$p.value))
for (a in c("go", "ipr")) {
  cl <- paste0("matched_", a); hs <- paste0("has_", a)
  mm <- d[get(cl) == TRUE]
  g1 <- 100 * mean(mm[divergent == TRUE][[hs]]); g0 <- 100 * mean(mm[divergent == FALSE][[hs]])
  add(sprintf("3 matching check: %s coverage in the matched-%s sets", toupper(a), a),
      "pct_point_gap", abs(g1 - g0), nrow(mm),
      sprintf("%.1f%% divergent vs %.1f%% uniform (both 100%% by construction)", g1, g0))
}

# --- a model, to see which axis is doing the work ---------------------------
mv <- d[matched == TRUE & is.finite(omega) & is.finite(degree_median) &
          degree_median > 0 & is.finite(infrv_member_mean)]
if (nrow(mv) > 300L) {
  fit <- glm(divergent ~ log(n_members) + log(degree_median) + omega +
               infrv_member_mean + has_go, data = mv, family = binomial())
  ci <- suppressMessages(confint.default(fit))
  co <- summary(fit)$coefficients
  for (v in rownames(co)[-1]) {
    add(sprintf("4 multivariable: %s", v), "log_odds", co[v, "Estimate"], nrow(mv),
        sprintf("CI [%.4f, %.4f], p = %.3g", ci[v, 1], ci[v, 2], co[v, "Pr(>|z|)"]))
  }
  add("4 multivariable model", "n", nrow(mv), nrow(mv),
      sprintf("omega is the thinnest axis: it covers %.1f%% of the matched set",
              100 * mean(is.finite(d[matched == TRUE, omega]))))
}

# --- is divergence shared between the species? ------------------------------
o <- fread(OTHER_UNIF)[cell %in% c(CELL_DIV, CELL_UNI), .(og_2sp, other = cell)]
both <- merge(d[, .(og_2sp, this = cell)], o, by = "og_2sp")
if (nrow(both) > 100L) {
  tb <- table(this_divergent = both$this == CELL_DIV,
              other_divergent = both$other == CELL_DIV)
  ft <- fisher.test(tb)
  add("5 cross-species: divergence shared", "fisher_or", ft$estimate, nrow(both),
      sprintf("p = %.3g | %d divergent in both, %d in neither",
              ft$p.value, tb["TRUE", "TRUE"], tb["FALSE", "FALSE"]))
  nullor <- replicate(NULL_REPS, {
    t2 <- table(both$this == CELL_DIV, sample(both$other == CELL_DIV))
    tryCatch(fisher.test(t2)$estimate, error = function(e) NA_real_)
  })
  add("5 cross-species: label-permutation null", "null_mean_or",
      mean(nullor, na.rm = TRUE), NULL_REPS,
      sprintf("sd %.4f over %d permutations; the observed OR must sit far above 1",
              sd(nullor, na.rm = TRUE), NULL_REPS))
  say(sprintf("\ncross-species: OR %.3f (p %.3g) on %s orthogroups classified in both | null OR %.3f",
              ft$estimate, ft$p.value, fmt_n(nrow(both)), mean(nullor, na.rm = TRUE)))
  write_tsv(as.data.table(as.data.frame(table(this = both$this, other = both$other))),
            file.path(OUT_DIR, sprintf("og_cell_cross_species_%s_vs_%s.tsv",
                                       STUDY, OTHER_NAME)))
}

# --- write -------------------------------------------------------------------
tests <- rbindlist(RES)
write_tsv(tests, file.path(OUT_DIR, sprintf("og_cell_tests_%s.tsv", STUDY)))
keep <- c("og_2sp", "cell", "divergent", "n_members", "matched", "matched_go",
          "matched_ipr", "matched_go_homeolog", "matched_ipr_homeolog",
          "has_go", "has_ipr", "has_tf", "n_go_terms", "n_ipr",
          "omega", "frac_unique", "dup_class", "expr_concordance",
          "sign_concordance", "n_responsive", "degree_median",
          "infrv_member_mean", "infrv_gain",
          if ("og_degree" %in% names(d)) "og_degree")
write_tsv(d[, ..keep][order(-divergent, og_2sp)],
          file.path(OUT_DIR, sprintf("og_cell_sets_%s.tsv", STUDY)))

banner("the contrast")
for (i in seq_len(nrow(tests)))
  with(tests[i], say(sprintf("  %-52s %-12s %12.4g  %s", test, statistic, value, note)))
say("done")
