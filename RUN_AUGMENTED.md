# Runbook — downstream pipeline on the augmented networks

Prepared 2026-08-09. Stage 0 and all script edits are **done**; what follows is
the compute, none of which has been run. Plan: `~/.claude/plans/quiet-growing-tower.md`.
Method and threshold rationale: `MI_for_review.md`.

## Already done (no compute needed)

- **Augmented edge files built** — `files/{sugarcane,purple/new}/network_*_augmented_edges.tsv`
  (76.2 M and 710.1 M edges; 866,148 and 28,963,013 MI-only edges).
- **Stage 0 archive** — 91 GB of Pearson-derived results moved to
  `files/pearson_baseline/` (see its `README.md`). Inputs left in place.
- **`scripts/general_stats.r` repaired** — it did not parse (orphan paren at the
  old `:29-30`) and read the edge file positionally on 5 columns. It now detects
  the format from the header, and gained `COMPUTE_TRANSITIVITY` (default `FALSE`).
  Smoke-tested end to end on 200 k edges.
- **Edge paths repointed** — `mcl_clustering.r:25,29` and
  `network_conservation_join.r:32,33`.

## Interpreter

`general_stats.r` needs igraph + hexbin, which `cor_env` does **not** have:

```bash
RSCRIPT=/home/genomics/miniconda3/envs/r_net_env/bin/Rscript
```

`r_net_env`, `R_popstat_jorge` and `r_clusterprofiler` all satisfy the full
dependency set; `comparative_network` and `matrix2-*` lack `hexbin`.

---

## Stage 1 — node + global metrics

```bash
cd /dados04/jorge/comparative_saccharum
$RSCRIPT scripts/general_stats.r
```

Runs both studies in one pass. Sugarcane minutes; purple ~1–2 h. It will **not**
rewrite `network_*_filtered_edges.tsv` — the augmented file is already the
filtered, weighted edge list, so it only produces metrics and plots.

Transitivity is off by default. Turning it on (`COMPUTE_TRANSITIVITY <- TRUE`,
`general_stats.r:21`) costs ~15 h on purple and is only ever *reported*, never
tested — `get_tfs/05:112` and `myb61/07:93`. Leave it off for the first pass.

Verify:
```bash
wc -l files/sugarcane/network_sugarcane_node_metrics.tsv \
      files/pearson_baseline/sugarcane/network_sugarcane_node_metrics.tsv
# expect ~103,364 vs 102,020 — the +1,344 nodes the MI layer adds
```

## Stage 2 — MCL clustering

```bash
export TMPDIR=/dados04/jorge/tmp && mkdir -p $TMPDIR
$RSCRIPT scripts/mcl_clustering.r
```

`TMPDIR` matters: the script writes a temporary `.abc` file for the MCL binary
(`:44-45`) that will be tens of GB for purple. Needs `mcl` on `PATH` (`:57`).
`INFLATION <- 2`, `NUM_CORES <- 100`.

Verify the augmentation actually changed the modules:
```bash
head -3 files/sugarcane/mcl_sugarcane_module_summary.tsv
head -3 files/pearson_baseline/sugarcane/mcl_sugarcane_module_summary.tsv
```

## Stage 3 — conservation join (run twice)

Edit `scripts/network_conservation_join.r:44 DIRECTION`, cheap direction first:

```bash
# 1) DIRECTION <- "sugarcane_to_purple"
$RSCRIPT scripts/network_conservation_join.r
# 2) DIRECTION <- "purple_to_sugarcane"     (~9x more outer edges)
$RSCRIPT scripts/network_conservation_join.r
```

If memory presses, lower `:46 CHUNK_SIZE <- 2e6` before anything else.

Do **not** run `full_matrix_edge_comparation.r` — it carries an OOM warning at
the old size, purple grew 19%, and its only consumer (step 14) is out of scope.

## Stage 4 — gene–trait correlations (run twice)

Swap the mutually exclusive comment blocks at `scripts/gene_trait_cor.r:43-56`
(currently purple), once per species. Fast.

⚠️ Pre-existing provenance mismatch: this reads `china/run2_onlyL` while the
purple network and the MI layer both came from `china/run1` + `Group1=="L"`
(68-gene difference, `scripts/README.md:540-542`). Not caused by the
augmentation, but it now propagates into the augmented results too — worth
resolving or stating explicitly in the methods.

## Stage 5 — conserved correlated genes

```bash
$RSCRIPT scripts/conserved_cor_genes.r      # no edits, trivial runtime
```

## Stage 6 — GO enrichment

```bash
conda activate topGO_env
# edit scripts/top_GO_conserved.r:41 ontology, once each for BP, MF, CC
Rscript scripts/top_GO_conserved.r

conda activate r_clusterprofiler
Rscript scripts/GO_semanthinc_enrichment.r  # loops all three itself
```

The new background universe matters most here: 1,344 extra sugarcane nodes
change the enrichment denominator.

## Stage 7 — H1 readouts

No path edits needed — the archive-and-swap layout means every path is already
correct. Skip the cached sequence-only steps.

```bash
scripts/get_tfs/04_postprocess.sh sugarcane
scripts/get_tfs/04_postprocess.sh purple
$RSCRIPT scripts/get_tfs/05_tf_characterization.r     # skip 01-03: hmmsearch is cached

scripts/myb61/run_all.sh 07                            # runs 07, 08; 01-06 cached
scripts/module20/run_all.sh 03                         # runs 03, 04, 05; 01-02b cached
```

Sanity check before this stage — the merges use `all.x = TRUE`, so an id
mismatch fails **silently as NAs**, not as an error:

```r
nodes <- data.table::fread("files/sugarcane/network_sugarcane_node_metrics.tsv")
tf    <- data.table::fread("GET_TFS/new/results/sugarcane/TF_in_network.tsv")
sum(sub("\\.v[0-9.]+$","",tf$gene) %in% sub("\\.v[0-9.]+$","",nodes$gene))  # must be > 0
```

(The new `node_metrics` ids have no `.v2.1` suffix because they come from the
augmented edge file; the readouts strip defensively, and stripping an absent
suffix is a no-op, so the join is safe — but check it anyway.)

---

## The analysis this was all for

Every augmented result can be split by how the edge was found:

```bash
awk -F'\t' 'NR>1{c[$9]++} END{for(k in c) print k, c[k]}' \
    files/sugarcane/network_sugarcane_augmented_edges.tsv
```

Questions worth asking of the re-run that the Pearson network cannot answer:

- Do the new modules recruit MI-only edges preferentially, or are they spread evenly?
- Do Module 20 or the MYB61 copies gain neighbours that were invisible to Pearson?
- Does the conserved-edge set become more or less conserved when non-monotone
  edges are admitted — i.e. is non-linear co-expression conserved between the
  two species at the same rate as linear co-expression?

## Out of scope

`eigengene.r` (8), `module_trait_cor.r` (9), `comparative_networks2.r` (10),
`full_matrix_edge_comparation.r`, `conserverd_edges_treatment.r` (14).

If you later want the module-trait branch, fix `module_trait_cor.r:183` first —
it is `fread(EDGE_FILE)` with no `select=`, i.e. the whole 69.8 GB purple file,
all 9 columns, into RAM.
