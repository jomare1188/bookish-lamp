# Decisions

Every non-obvious choice in this pipeline, with the evidence behind it. Entries
are dated; the date is when the decision was made, not when it was written down.

---

## 2026-08-13 — Purple uses `china/run2_onlyL`, not `china/run1`

**The problem.** The old pipeline was internally inconsistent about which purple
quantification it used. The network came from `china/run1` (36 libraries, leaf
and root) narrowed to the 18 leaf libraries with `Group1 == "L"`. But
`gene_trait_cor.r:50`, `myb61/08` and `module20/04` all read
`china/run2_onlyL` — a separate nf-core run of the leaf libraries alone.

**Why it matters.** These are not two views of one quantification. Measured:

| | run1 leaf subset | run2_onlyL |
|---|---|---|
| raw counts | — | differ, max abs 6,398 |
| VST | — | differ, mean abs 0.071, max 2.34 |
| genes at CV ≥ 15 | 170,736 | 170,740 (68 in symmetric difference) |

Size factors and the dispersion fit differ because DESeq2 saw 36 samples in one
case and 18 in the other, so every VST value for the leaf libraries differs.

**Decision.** `run2_onlyL` everywhere. It is the dedicated leaf-only run, it is
*already* exactly the 18 libraries, and using it removes the `Group1 == "L"`
filter as a source of drift.

**Consequence.** n = 18 is unchanged, so every threshold in
[thresholds.md](thresholds.md) that depends on n carries over untouched. The
purple network had to be rebuilt from the new VST, which is affordable only
because of the next decision.

---

## 2026-08-13 — No dense correlation matrices; both layers come from one engine

**The problem.** `pearson_cor.r` wrote an n × n correlation matrix and a second
one of p-values — 402 GB per study, 804 GB for the two — and `build_edgelist.r`
streamed them once to produce a 105 GB intermediate edge list, which was then
thresholded down to a few GB. Rebuilding purple through that path would have
cost days and another ~290 GB.

**Decision.** Add `pearson` as an estimator to the existing GPU engine
(`02_network_engine.py`), so the thresholded edge list comes straight out of the
sweep. The intermediates never exist.

**Result.** The sugarcane linear layer — all 1.46 × 10¹⁰ pairs — takes **5.5
minutes** end to end, at 1.6 billion pairs/s. The old path took hours and 220 GB.

**Bonus, and not a small one:** both layers now read the same
`<study>.f32`, so their gene sets are identical *by construction*. The old
`00_export_vst.r` had to verify its gene list against the header of the 200 GB
Pearson matrix and abort on mismatch. That check is gone because the situation
it guarded against can no longer arise.

---

## 2026-08-13 — The old networks contained ~47,000 edges below their own threshold

**Found while validating the new linear layer against the old one.** The new
sugarcane layer has 75,333,769 edges; the old `network_sugarcane_filtered_edges.tsv`
has 75,380,961. The shortfall was 47,192 (0.063%) — and it *grew* with |r|:

| cut | old | new | deficit |
|---|---|---|---|
| \|r\| ≥ 0.80 | 75,380,961 | 75,333,769 | 0.063% |
| \|r\| ≥ 0.85 | 36,443,207 | 36,415,235 | 0.077% |
| \|r\| ≥ 0.90 | 12,795,207 | 12,779,788 | 0.120% |
| \|r\| ≥ 0.95 | 1,856,204 | 1,851,315 | 0.263% |

That pattern is not float precision — precision only moves edges within ~10⁻⁶ of
a cut. So it was checked against ground truth: r computed in **float64** from the
DESeq2 VST for the first 6,000 genes, thresholded at 0.8 ≤ |r| ≤ 0.9999.

```
TRUTH (float64):  73,404 edges
new layer:        73,404 edges   truth-only: 0   new-only: 0     <- exact match
old network:      73,453 edges   truth-only: 0   old-only: 49    <- 49 spurious
```

The new layer reproduces float64 ground truth **exactly**. The old network has
extra edges and misses none. Every one of the 49 has its r reported as exactly
±0.8.

**Cause.** `matrix_sugarcane_pearson.tsv` stored r rounded to 4 decimal places:

```
1  0.5748  0.0368  0.006  0.2628  -0.1126  -0.0528  0.699  -0.2165 ...
```

A true r of 0.79996 was written as `0.8`, then passed `abs_r >= 0.8` in
`general_stats.r`. The error is strictly one-directional — it can only admit
edges, never drop them — and its relative size grows at higher cuts because the
|r| distribution decays steeply, so a fixed 5 × 10⁻⁵ rounding window is a larger
fraction of an ever-smaller tail.

**Consequence.** This is a defect in the *old* results, inherited by everything
downstream of them, and it is fixed here as a side effect of dropping the dense
matrices: the new pipeline thresholds the float value directly and never
round-trips r through text.

---

## 2026-08-13 — Gene ids are normalised once, at export

**The problem.** R570 ids carry a `.v2.1` suffix; the proteomes, orthogroups and
TF tables use the bare form. In the old tree the strip happened in different
places at different times — `network_sugarcane_node_metrics.tsv` (Jun 2) had the
suffix while `network_sugarcane_filtered_edges.tsv` (Jun 25) did not — and every
readout carried a defensive `sub()` to cope. Those merges use `all.x = TRUE`, so
a mismatch produces silent NAs, not an error.

**Decision.** Strip in `01_export_vst.r`, once, with a collision assert. Every
file under `results/` carries the bare id. The defensive strips downstream become
harmless no-ops (stripping an absent suffix does nothing), so they stay as
insurance against files from outside this tree.

---

## 2026-08-13 — `--match-pearson 0.8`, not 0.7

Carried over from the previous method note, restated because it is easy to get
wrong. `build_edgelist.r` used `PEARSON_THRESHOLD <- 0.7`, but that was only a
pre-filter for an intermediate file. The networks actually analysed came from
`general_stats.r`, whose `PEARSON_MIN <- 0.8` is the real cut — verified against
the data: both `network_*_filtered_edges.tsv` bottom out at exactly |r| = 0.8.

In this pipeline 0.7 survives only as `CAND_PEARSON`, the candidate cut, which
exists so BH is exact over the full rejection set rather than only over the edges
finally kept.

---

## 2026-08-13 — Every toggle became an argument

The old tree needed at least 8 manual file edits to complete one pass: species
comment-block swaps in five scripts (`gene_trait_cor.r` set `NETWORK_NAME` twice,
at `:26` and `:48`), plus `DIRECTION`, `ontology`, `TRAIT_ENCODING` and
`SPLIT_BY`. An edit left in the wrong state produces results that look fine and
are wrong. All of it now comes from `config.sh` through `run.sh <stage> <arg>`.

---

## 2026-08-13 — Things deliberately left off

- **Local transitivity** (`COMPUTE_TRANSITIVITY=0`). O(Σ deg²), ~15 h on purple,
  and nothing in the pipeline *tests* it — it appears only as a reported column
  in the TF and MYB61 readouts. Set to 1 to backfill.
- **The 8 GB `mcl_results_both_networks.rds`** (`MCL_SAVE_GRAPH_RDS=0`). It held
  the full igraph objects for both networks and was written into the sugarcane
  directory regardless of study. Nothing reads it; every consumer reads
  `mcl_*_membership.tsv`.
- **`full_matrix_edge_comparation.r`.** Carries its own OOM warning at the old
  network size. Its only consumer is `conserverd_edges_treatment.r`, which is out
  of scope. Note for whoever picks it up: that script reads
  `conserved_edges_full_*` (this script's output) rather than
  `conserved_edges_*_FULL.tsv` (the join's output that actually completes) — the
  names differ by one underscore and they are different files.
- **`eigengene.r`, `module_trait_cor.r`, `comparative_networks2.r`.** The
  module-trait branch. **Superseded, not skipped** — this is now stages 14–18
  (`14_module_eigengene.r` … `18_module_go.r`), rebuilt rather than ported. Three
  bugs in `eigengene.r` did not survive the rewrite and are listed in
  `14_module_eigengene.r`'s header. `module_trait_cor.r` was not ported at all:
  its `module_trait_cor.r:183` reads the whole purple edge file, all 9 columns,
  into RAM with no `select=`, and the module-level test does not need edges —
  `moduletrait` runs on the eigengene matrix instead.
- **Nothing outside `new_clean/` was modified or deleted.** The 804 GB of dense
  matrices and 105 GB of intermediate edge lists are now unreferenced and
  reclaimable, but this pipeline does not touch them.
