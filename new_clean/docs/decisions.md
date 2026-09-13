# Decisions

Every non-obvious choice in this pipeline, with the evidence behind it. Entries
are dated; the date is when the decision was made, not when it was written down.

---

## 2026-09-10 — the module analysis moves to per-species inflation, and the design goes into the module test

**Supersedes the `-I 2` entry below for the module analysis.** That entry asked
which single inflation to use for both networks and answered `-I 2`, on the
grounds that it keeps 98.1% of each species' peak modularity and keeps the two
networks comparable. The analysis now uses **each species' own optimum** —
`-I 1.5` sugarcane, `-I 3.5` purple — on the unpruned Pearson-only graphs.

**Why the change is defensible.** The two questions are not the same. A single
inflation matters when the two partitions are being compared *to each other*, as
the conservation and cross-species stages do. The module analysis is run *within*
a species: eigengenes, module-trait response, module GO. There, taking each
network's own optimum costs nothing in comparability and gains the granularity the
criterion actually selected. `-I 2` remains the right answer for anything that
joins the two species, and the figure-12 argument stands as written.

**What it costs, stated plainly.** Purple's optimum is a finer partition with a
long tail: 29,624 clusters of which 15,600 are singletons (9.17% of genes, now
"Unassigned"), against 276 singletons at `-I 2`. The analysable core is unaffected
— 7,493 modules get an eigengene and 1,063 hold ten genes or more.

**The module-trait test is now blocked, and blocked is primary.**

```
sugarcane   eigengene ~ genotype + segment + N      n = 48, residual df 42
purple      eigengene ~ genotype + N                n = 18, residual df 15
```

Genotype carries R² = 0.999 of PC1 in purple and 0.998 in sugarcane, and leaf
segment 0.802 of PC2, so the largest variance component in either matrix had been
in the residual of every marginal test. Spearman stays primary — the ordinal-dose
argument is unchanged — as the blocked fit on midranks.

The marginal statistic is computed on the *same* eigengenes and reported beside
the blocked one, so the effect of the model is visible in one table. It more than
doubles the responsive set in both species (sugarcane 251 → 588, purple 96 → 182)
and **loses nothing**: blocked is a strict superset. A test that only loosened a
threshold could not have that property, and a block absorbing signal rather than
noise would have broken it.

**What licenses trusting it:** the solver reproduces `lm()` to 4.4e-16; the model
matrix is rank-checked against a block confounded with the trait; the unblocked
Spearman path reproduces `cor.test(exact = FALSE)`; sugarcane's plant-level control
(12 plants × 4 segments are repeated measures, not 48 replicates) agrees at
r = +0.9486; and the permutation null runs **within block**, because free shuffling
would break the structure the model conditions on. No permutation of 1,000 reached
the observed count in either species.

`fit_blocked` and `verify_against_lm` live in `scripts/lib/common.R` so the gene
level (`31_gene_trait_blocked.r`, branch `blocked-gene-trait`) and the module level
share one verified implementation; when that branch is merged, 31 should be
switched to the shared copy.

---

## 2026-09-10 — `-I 2` stays, on evidence; Louvain and Leiden-modularity are ruled out

**The problem.** The pipeline had always clustered at `-I 2`, a value nobody
chose, on a graph nobody had compared against another algorithm. Two open
questions, both fair to ask of any paper.

**What was done.** Both studies' unpruned Pearson-only graphs, every partition
scored by `clm info` against the same matrix so only the partition varies: an
MCL inflation ladder (1.2 to 13 in sugarcane, 1.2 to 6 in purple), Leiden with
the CPM objective over a resolution ladder, Leiden with the modularity objective,
and Louvain. Plus an independent biological check — PFAM Sørensen–Dice
homogeneity above a size-matched null.

**The decision, and why.**

* **Keep `-I 2`.** `mf`, `af` and `eff` are all monotone over MCL's usable range,
  so each picks a grid boundary and none can choose an inflation. Modularity is
  the only criterion with an interior optimum. Resampled densely around the peak,
  it lands at `-I 1.5` (sugarcane, plateau 1.4-1.7) and `-I 3.5` (purple, plateau
  2.7-4.5) — so `-I 2` is not the optimum in either species, but it keeps **98.1%
  of the peak in both**, which is the argument for one inflation across two
  networks. The setting was never justified before; it is now.
* **Never Louvain or Leiden-modularity on these graphs.** They post the highest
  modularity of anything tested and reach it by building giant modules — purple's
  optimum is **50 clusters for 170,135 genes**, the resolution limit at
  2m = 1.35e9. Their PFAM excess is +0.0002 and +0.0004, i.e. no annotation
  signal at all above a random partition of the same shape. This also settles
  what the giant module is: not an MCL artefact, but what modularity wants.
* **Leiden CPM is the option, and only if the giant module must go.** It is the
  only method that dissolves it, and holds the highest `eff` on either graph with
  a genuine interior optimum (gamma 0.1 sugarcane, 0.2 purple). But sorted by
  area fraction the two methods interleave on **one curve**, with MCL the higher
  wherever they overlap and carrying far fewer singletons, and MCL also wins the
  PFAM check at matched module count in both species. Leiden CPM's advantage is
  reach, not quality: MCL's area fraction bottoms out because inflation
  underflows before it can fragment further.

**The caveat that turned out not to be one.** MCL prunes each node's neighbour
list while computing (`-scheme 7` keeps 1,200, against mean degrees of 1,477 and
7,946) and grades that pruning "deplorable" to "abominable". Measured rather than
assumed: at `-S 4000` and `-S 10000` sugarcane's `-I 1.7` is identical and within
0.7% of the default, and purple's `-I 3` at `-S 10000` moves every statistic by
under 0.3%. The grade says how much was discarded, not whether it mattered.

**Three silent tool failures were found doing this, and are now guarded**
(see `docs/results.md`): `clm info`'s `eff` and `mf` depend on which other
clusterings share the call; mcl ignores `-I` above 30 and uses the default
instead; and mcl underflows at high inflation, returning plausible-looking
partitions computed on a graph a third of which has gone to zero.

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
  module-trait branch. **Superseded, not skipped** — this is now stages 14–19
  (`14_module_eigengene.r` … `19_module_trait_spearman.r`), rebuilt rather than
  ported. Three
  bugs in `eigengene.r` did not survive the rewrite and are listed in
  `14_module_eigengene.r`'s header. `module_trait_cor.r` was not ported at all:
  its `module_trait_cor.r:183` reads the whole purple edge file, all 9 columns,
  into RAM with no `select=`, and the module-level test does not need edges —
  `moduletrait` runs on the eigengene matrix instead.
- **Nothing outside `new_clean/` was modified or deleted.** The 804 GB of dense
  matrices and 105 GB of intermediate edge lists are now unreferenced and
  reclaimable, but this pipeline does not touch them.

---

## 2026-08-22 — The module response is Spearman only: not Pearson, not MI

**The problem.** The module-level nitrogen call was made by running
`12_gene_trait_mi.py` on the eigengene matrix, which returned Pearson *and*
mutual information and sorted every module into `pearson_only` / `mi_only` /
`both`. Reusing the gene-level script was a genuinely good structural decision —
one validated code path, no second implementation of the statistics — but it
imported two statistics that do not suit this question, and a three-way
classification that nothing downstream could use.

**Pearson reads an ordinal trait as an interval one.** Purple's nitrogen levels
are 0, 2 and 6 mM. Pearson on `{0, 2, 6}` asks whether a module moves exactly
twice as far from 2 to 6 mM as from 0 to 2 — an arithmetic claim about the
response curve that the design never made and the experiment cannot test.

(Recorded 2026-08-22, after this decision: those three levels are **not** a dose
series at all. 2 N is the control and 0 N / 6 N are opposite stresses, so a
monotonic test of any kind asks about nitrogen *supply*, not nitrogen *stress*.
The measured cost of that is one module — see
[results.md](results.md#what-responds) — but the wording of every purple module
claim should say "tracks nitrogen supply".) A
module that responds and then saturates above 2 mM is a real monotone nitrogen
response, and Pearson discounts it. Spearman asks only the question the gradient
poses: does this module move monotonically with nitrogen? On sugarcane's
two-level trait it becomes the rank-biserial correlation, which is the same test
with a different name and is robust to the outliers a PC1 can carry.

The cost of the wrong choice was measured on the identical eigengenes at the
identical thresholds (padj ≤ 0.05, |r| ≥ 0.6):

| | Pearson | Spearman | Spearman only | Pearson only |
|---|---|---|---|---|
| sugarcane | 408 | 465 | 71 | 14 |
| **purple** | **38** | **79** | **45** | 4 |

**Purple's responsive set more than doubled**, and purple is exactly the study
with the three-level gradient. This is the largest single effect of the change.

**MI at this level was an omnibus test finding artifacts.** MI fires on any
dependence at all, including a dispersion change or a distributional quirk in one
or two libraries. That is a virtue at the edge level, where the question is
whether non-linear co-expression exists. At the module level it produced a
`mi_only` class of 239 sugarcane modules of which 35% had two of 48 samples
carrying over a quarter of the eigengene's variance — worse than the
non-responsive background at 14%. The Spearman-selected set runs the other way:
20% against 50% in the modules it did not select.

**And the MI threshold could not be derived, only calibrated.** `padj ≤ 0.05`
alone admitted modules down to |r| = 0.356, so an effect-size floor was needed.
Flooring the linear side alone was worse than no floor — it pushed every
modestly-linear module into `mi_only`, 253 → 786, where it then read as
"non-linear" when it was nothing of the kind. But the network's nats-to-|r|
identity could not supply the MI equivalent, because it assumes two continuous
variables and this trait is discrete with MI capped at H(trait). The floor was
therefore set to the median MI among modules sitting at |r| ≈ 0.6 — calibrated
against the very statistic it was meant to be independent of. A threshold that
can only be defined by reference to the alternative is not an independent test.

**What changed.** `19_module_trait_spearman.r` replaces `12_gene_trait_mi.py` in
the `moduletrait` stage. One statistic, one rule — `padj <= MODULE_PADJ_THR` and
`|rho| >= MODULE_R_THR` — owned by that script alone; 15, 16, 17 and 18 read
`responsive` and `direction` and do not recompute anything. The three response
classes are gone; the figures and the per-module GO now split by the sign of rho,
which is different biology rather than two ways of detecting the same thing.

**What it bought, beyond purple's modules.** The TF result came back. Under the
three-class call, TF enrichment among responsive modules was a null spread thin:
`both` at OR 1.96, p = 0.061, the "strongest" class `pearson_only` at 3 of 41
modules. One responsive set of 465 gives **OR 2.65, p = 0.0016**, concentrated in
the modules that fall with nitrogen. The signal was there; three classes too
small to separate had been hiding it. Purple's testable module count for GO also
went from 3 to 10.

**What it did not change.** Every headline module survives: Module_026, _100,
_440, _267, _469, _1820, _009, _514, all rising with nitrogen, and the BP↔MF
cross-check (nitrate assimilation ↔ nitrate reductase activity; ammonia
assimilation ↔ glutamate synthase) holds unchanged. Module_117 (urea transport)
drops out at |rho| = 0.43.

**Scope.** This is the module level only. The gene-level selection
(`07_gene_trait_cor.r` Pearson, `12_gene_trait_mi.py` MI) is untouched, because
the conserved-response, edge-level and GO results all depend on it and changing
it would invalidate them. The ordinal argument applies there too and is recorded
as an open item in [results.md](results.md#open-items), not as a settled choice.

**The one thing to read carefully.** Purple's permutation null has a long right
tail — 1,000 label shuffles gave a median of 0 responsive modules but a maximum
of 244, against 79 observed (empirical p = 0.002). Purple's network is a single
dense component, so its eigengenes are strongly correlated and one lucky shuffle
lights up many at once. The 79 are a real signal and are **not** 79 independent
findings.

---

## 2026-08-23 — GO is asked at three grains, and the middle one is new

**The problem.** `18_module_go.r` tested one module at a time and, in its global
figure, all responsive modules pooled. Between those two sits the question a
reader asks first — *what do the modules that go UP with nitrogen do, as a set,
and is it different from the ones that go DOWN?* — and neither grain answers it.
Pooling averages the two directions together; per-module testing answers about
individual modules and does so for very few of them.

**The annotation gate is why "very few" is the right description.** A module needs
`MODULE_GO_MIN_ANNOTATED` GO-annotated members for the test to be defined. Of
sugarcane's 465 responsive modules **49** clear that, and 10 of purple's 79 — so
the per-module grain describes ~11% of the responsive set, biased toward the
large modules, and cannot be read as characterising it.

**The fix, and why it is not just convenience.** The union of every rising module
is tested against the union of every falling one, on the same `topGOdata` object,
the same network-node background, the same weight01 statistic and the same raw-p
threshold, so all three grains stay directly comparable. Pooling by direction
dissolves the gate: modules too small to test alone all contribute, and the two
sets carry 383 and 161 annotated genes drawn from 242 and 223 modules.

**What it found.** Of sugarcane's 103 enriched BP terms, **74 are enriched only
among the modules that rise, 23 only among those that fall, and 6 in both**
(purple: 7 / 23 / 3). The split is not one functional programme cut in half, so
pooling the directions was averaging two distinct programmes into one list. And
the two are interpretable in a way neither other grain recovered: nitrate
assimilation, nitric oxide, proline/asparagine/spermidine biosynthesis and a
defence block go UP; flavonoid biosynthesis, raffinose-family oligosaccharides,
triglyceride biosynthesis and cold response go DOWN.

**What it does not license.** The direction labels mean "tracks nitrogen supply",
not "responds to nitrogen stress". In purple the three levels are
stress-control-stress, so a module moved the same way by both extremes is
invisible to a monotone test and therefore absent from both direction sets. See
`Soffic.09G0001580-9H` in [results.md](results.md#module-20-in-purple--one-copy-responds-and-not-monotonically)
for a gene that does exactly that.

---

## 2026-08-26 — A stochastic block model is fitted alongside MCL, and neither is chosen yet

**Why.** MCL's partition of the sugarcane network has a shape that makes the
module-level analysis awkward: 10,309 modules with a **median of 3 genes** and a
single module holding 19,604 (19% of the network). A median of 3 is below the
GO annotation gate almost by construction, which is why only 49 of the 465
responsive modules could be tested for function at all. A stochastic block model
was fitted to the same network to see whether a different partition removes that.

**How it is wired, and why nothing was rewritten.** Every module-level stage
reads the clustering through exactly two files. `27_sbm_membership.r` writes
those two files from the SBM fit in MCL's exact schema, so **no consumer changed**
— the alternative is carried through the identical eigengene, Spearman,
TF-hypergeometric and topGO code. A `CLUSTERING=mcl|sbm` switch selects which
clustering is read and sends its outputs to a separate directory, so the MCL
results are untouched by construction rather than by care. Level 0 of the nested
fit is used: 1,009 blocks, median 51 genes, no giant block. Levels 1 and up
coarsen fast (294 blocks by level 1, 49 by level 3) and stop being a module set.

**They are not the same structure re-parameterised.** Adjusted Rand index between
the two partitions is **0.0156** over all 103,336 genes. Newman modularity is
0.1005 for MCL and 0.0081 for the SBM — expected, since an SBM minimises a
description length and does not optimise modularity, so this is the two methods
answering different questions rather than one failing.

**What the comparison shows** (`./run.sh figclustering`, sugarcane):

| | MCL | SBM level 0 |
|---|---|---|
| modules / median size / largest | 10,309 / 3 / 19,604 | 995 / 51 / 1,366 |
| genes left unassigned | 401 | 14 |
| modules tested | 6,576 | 979 |
| responsive | **465** | **23** |
| median PC1 variance explained | 81.1% | 77.6% |
| responsive modules GO-testable | **49 (11%)** | **20 (87%)** |
| BP terms returned | 246 | **325** |

The trade is stark. MCL finds twenty times more responsive modules, of which
**nine in ten cannot be tested for function**. The SBM finds 23, of which almost
all can, and they return more GO terms than MCL's 465 do. Eigengene coherence is
close, slightly favouring MCL — but its median is dominated by very small
modules, where one component explains most of the variance almost by
construction, so the two medians are not measured on comparable objects.

**Not decided, and deliberately waiting.** Which partition is preferable depends
on whether the next step is to name modules or to count them — and that question
cannot be settled on sugarcane alone. ~~The purple fit is running on another
machine as of 2026-08-27~~ — **this was wrong, corrected 2026-09-03**: the purple
fit never ran. `sbm/output/purple/` holds only two `.gt` files, nothing since
2026-08-25, and `sbm/code/purple/sbm_3.6_saveall.py` calls `lowmem.load_gt`
(`:915`) and `check_weighted_model_lowmem` (`:924`) without ever importing
`lowmem` — a `NameError` on both lines. Until a purple fit lands the comparison
covers one species, and a change that would regenerate every module-level figure
should not rest on one species.

**The modularity comparison in this entry is also superseded.** MCL's 0.1005 was
unweighted and the SBM's 0.0081 weighted; see the 2026-09-03 entry on that.

The specific thing purple will decide is whether the trade above is a property of
the METHOD or of this network. Sugarcane's MCL partition is pathological in a
particular way — a median module of 3 genes and one holding 19% of the nodes —
and if purple's MCL partition is better behaved, the SBM's advantage may shrink
to nothing there. Purple is also the harder case for the SBM: its network is a
single dense component of 170,736 nodes with mean degree 8,265, where MCL already
produces 9,881 modules with a median of 3, so whether a block model finds
structure in it at all is an open question rather than a formality.

MCL remains the default until then. Switching is one line (`CLUSTERING=sbm`), but
figures 5–7, the Module-20 readout and the module sections of results.md all
describe the MCL partition and would all need regenerating — so the switch is a
decision about the paper, not just about the pipeline.

**What is already in place for when it lands.** `sbmclust` takes a study argument
like every other stage, `SBM_DIR_purple` is already configured, and the comparison
figure takes a study argument too — so purple needs `./run.sh sbmclust purple`,
the module chain under `CLUSTERING=sbm`, and `./run.sh figclustering purple`.
Nothing else has to change.

---

## 2026-09-03 — The giant module is a node-degree artefact, and MCL said so all along

**The problem.** MCL at `-I 2` gives sugarcane one module holding 19,604 genes
(19.0% of the network) and purple one holding 47,887 (28.0%), both with a median
module of 3. The instinct was to tune inflation.

**Inflation is not the handle, and the measurement says so.** Sweeping `-I` from
1.4 to 6 on the unreduced sugarcane graph moves the largest module from 25.3% to
12.6% — it never stops being a giant — while stranding 12,905 genes as singletons
at the top end. mclfaq(7) §7.3 names the real cause: *"Preferably the network
should not have nodes of very high degree... Such nodes tend to obscure cluster
structure and contribute to coarse clusters."* clmprotocols(5) puts a number on
it for co-expression graphs: *"the median node degree should be at most one
hundred neighbours."*

| | nodes | median degree | p90 | max |
|---|---|---|---|---|
| sugarcane | 103,336 | 52 | 6,666 | 16,019 |
| purple | 170,736 | **859** | **33,071** | **40,960** |

**MCL had been grading its own work as unusable, into `/dev/null`.**
`05_mcl_clustering.r:67` passed `stderr = FALSE`, which discarded the jury
pruning synopsis — mcl's own verdict on whether its resource scheme kept enough
of each node's neighbourhood. Recovered, on the unreduced sugarcane graph:

| -I | 1.4 | **2 (production)** | 3 | 4 | 6 |
|---|---|---|---|---|---|
| jury | 34.8 woeful | **39.2 deplorable** | 43.4 poor | 45.5 dodgy | 47.2 shabby |

The clustering every module-level result in this project rests on was computed
under pruning mcl itself calls **deplorable**, and no inflation value fixes it.
The cause is arithmetic: the default `-scheme 7` tracks at most ~1,200 neighbours
per node, and purple's top decile has more than 33,000.

**Decision.** Apply mcl's documented remedy, a k-NN reduction (`-tf '#knn(k)'`),
to the matrix mcl clusters. **The network is not touched** — edges, conservation,
and every published figure are unchanged; only what mcl is handed is reduced. k
is chosen by the protocol's own criterion, measured with `mcx query -vary-knn`:
the largest k that still meets the median-degree target, so the reduction is the
gentlest one that does the job.

Sugarcane, same graph, same inflation, reduction the only difference:

| | largest module | modularity Q | jury |
|---|---|---|---|
| `-I 2`, no reduction | **18.97%** | 0.0819 | 39.2 deplorable |
| `-I 2`, `#knn(180)` | **0.57%** | 0.2692 | **70.0 adequate** |
| `-I 1.4`, `#knn(180)` | 3.07% | **0.5506** | 54.6 tolerable |

Modularity rises **6.7-fold**. The community structure was there; the hub edges
were masking it. The cost is singletons: 401 at the current setting, 16,144 at
`#knn(180) -I 2`.

**Supporting infrastructure.** Each network is now loaded once into mcl's native
binary format (`34_mcl_load.sh`), replacing the 36 GB text `.abc` file
`05_mcl_clustering.r` rewrote on every run. Sugarcane: 4.67 GB and ~4 min becomes
1.2 GB and 2.5 min, peak RSS 2.4 GB; purple's 705,571,723 edges load in 11 GB at
22.5 GB peak. That is what makes a sweep affordable at all.

`--stream-mirror` is required and is not optional: the edge table stores each
undirected edge once, so without it `mcxload` builds a directed matrix in which
roughly half the nodes have out-degree 0, and `#knn` — which needs both
endpoints' neighbour lists — reduces the graph to zero edges. `34` asserts zero
degree-0 nodes for this reason.

**Judged on four criteria, not one.** `clm info` (efficiency, mass fraction, and
**area fraction**, which is the sum of squared cluster sizes over N² and
therefore the giant-module statistic), `clm dist` for where the partition stops
moving, the jury grade, and an annotation-homogeneity score against a
size-matched null.

**cogeqc was considered and not used.** It is a comparative-genomics QC package —
BUSCO, orthogroup inference, synteny — with no co-expression module functions.
Its `calculate_H()` is the right *metric* (mean pairwise Sørensen–Dice over gene
annotations) but its implementation cannot serve here: `max_size = 200` makes it
return nothing for larger groups, so it would silently refuse to score exactly
the giant modules being diagnosed, and it enumerates pairs with `combn()` in an R
loop. `37_cluster_homogeneity.r` computes the same score from a sparse binary
gene × domain matrix and **subsamples** large modules rather than skipping them.
The package is installed only in `comparative_network`, `matrix2-core` and
`matrix2-expression`, none of which this pipeline uses.

**Nothing is adopted yet.** The sweep is diagnostic; `MCL_KNN_*` and
`MCL_INFLATION_*` stay empty, so `./run.sh mcl` still reproduces the shipped
partition exactly. Changing them regenerates `mcl_*_membership.tsv`, which every
module-level stage reads, and that is a decision about the paper rather than
about the pipeline.

---

## 2026-09-03 — Modularity was being compared between two different statistics

`05_mcl_clustering.r:89` computed `igraph::modularity(g, mem)` — **unweighted**.
`27_sbm_membership.r:125` computed `modularity(g, comm, weights = E(g)$weight)` —
**weighted**. igraph 2.1.4 does not pick up the `weight` edge attribute when
`weights` is left NULL (verified directly). So the comparison quoted throughout
this project — MCL 0.1005 against SBM 0.0081, a 12-fold gap — was between two
different quantities, while `decisions.md`, `methods.md`, `config.sh` and the
figure-8 legend all asserted the two were measured the same way.

Both stages now compute and write **both**, as `modularity_Q` (unweighted) and
`modularity_Q_weighted`. `28_fig_clustering_compare.r` reads one column and says
so. The corrected like-for-like comparison replaces the old one wherever it is
quoted.

The adjusted Rand index in `28_fig_clustering_compare.r:104-113` was checked at
the same time and is **correct**: hand-rolled Hubert–Arabie, 0.0155872455 against
`clue::cl_agreement(method = "cRand")`'s 0.0155872455, difference 0. It is kept.

---

## 2026-09-03 — |r| >= 0.9 tested as a parallel network, not adopted

**Why.** The giant module was diagnosed as a node-degree artefact and treated with
a k-NN reduction applied after the network is built. A stricter correlation
threshold attacks the same cause earlier — it changes what counts as an edge —
so it is the more principled version of the same idea and worth measuring against
the incumbent.

**What was known before building anything.** `mcx query --vary-correlation` on
the existing networks already answers most of it:

| at \|r\| ~ 0.90, from the survey | sugarcane | purple |
|---|---|---|
| median node degree | 52 -> **3** | 859 -> **112** |
| genes left isolated | 30.4% | 0.3% |

The degree numbers held up. **The gene-loss numbers did not, and the reason is
the trap recorded below.** `mcx query --vary-correlation` cuts the *weight*
scale, which is derived from `stat`, and `stat` for an MI edge is its
`r_eq` — about 0.9 by construction at the 0.8-matched floor. So every MI edge
survives a weight cut and keeps its nodes attached. A real rebuild re-matches the
MI floor to the new threshold, the MI layer collapses, and the nodes it was
holding in fall out. Measured after the rebuild:

| purple, genes with at least one edge | |
|---|---|
| survey predicted at \|r\| ~ 0.90 | 170,212 (0.3% lost) |
| **actual 0.9 rebuild** | **147,446 (13.6% lost)** |

A 45-fold error in the predicted cost, in the optimistic direction. The survey is
still the right tool for choosing a **k-NN** cut, where no re-matching happens;
it is the wrong tool for predicting a threshold change, and this entry exists
partly to say so.

**And one consequence that runs against the intuition.** A common threshold makes
the two studies *less* comparable, not more. The p implied by the cut moves
9.0e-12 -> 3.3e-18 in sugarcane but only 6.7e-05 -> 3.7e-07 in purple, so the
specificity gap widens from 7.5e6-fold to 1.1e11-fold. The threshold fixes
purple's degree distribution; it does not fix the asymmetry between the studies,
and raising it is not a route to doing so.

**Decision: build it as a separate tree, adopt nothing.** `RESULTS=…/results_r09
STAT_MIN=0.9 ./run.sh build <study>`. The main tree is untouched.

**Why a rebuild and not a filter of the existing edge table.** For sugarcane the
two are equivalent, and provably so: `n_perm` is capped at `MAX_PERM = 2e8` under
both thresholds and the candidate cut is unchanged, so the KSG null is
bit-identical and filtering `pval <= p(|r| = 0.9)` reproduces the rebuild. It was
used as an independent prediction and the rebuild matched it to **0.013%** —
12,778,116 edges against 12,779,788 predicted, the residual being the `%.5f`
printing width at the boundary.

**For purple it is not equivalent, and this is the reason the shortcut was
rejected.** `n_perm` auto-sizes from the target p (`02_network_engine.py:907-911`),
so purple's permutation null goes from **1,000,000 draws at 0.8 to 81,986,326 at
0.9** — a different null sample, a different GPD tail fit, a different upper
endpoint (2.4376 -> 2.2157 nats), and therefore different p-values for the same
MI value. A filtered purple table would have carried 0.8-calibrated p-values
under a 0.9 label. Measured after the rebuild, all as predicted.

**Two traps recorded so they are not rediscovered.**

- **Never filter on `stat`.** MI edges enter the 0.8 table at `r_eq ~ 0.907`, so
  `stat >= 0.9` keeps the *entire* MI layer regardless of threshold: 797,576 of
  sugarcane's 866,575 MI-only edges pass a `stat` filter that a real 0.9 rebuild
  rejects. Only `pval` is on the matched scale, because `MATCH_PEARSON` follows
  `STAT_MIN` and both layers are cut at the same target p.
- **`03_merge_layers.py:117-124` stops a half-done job.** Feeding a 0.8-matched
  MI layer into a `--stat-min 0.9` merge exits with an error, which is correct:
  both layers must be rebuilt together or the union is not licensed.

**A hazard fixed on the way, which is the real reason this entry exists.**
`RESULTS`, `STAT_MIN`, `STAT_MAX`, `MATCH_PEARSON`, `CAND_PEARSON` and
`MCL_WORK_DIR` were all unconditional assignments, and `layer_out()` /
`network_tsv()` carry no tag or suffix. `STAT_MIN=0.9 ./run.sh build purple`
would have been silently reset to 0.8 — and had it taken effect, it would have
overwritten the 0.8 network **in place**: 77 GB of edge tables and ~5.5 h of GPU,
with nine `run.sh` call sites picking up the new file without comment. All six
now use the `${VAR:-default}` form, and `network`/`merge` read the threshold out
of the summary json beside their target and refuse to overwrite an output built
at a different one. `FORCE=1` overrides.

**What the layers came out at.**

| layer | 0.8 | 0.9 |
|---|---:|---:|
| sugarcane pearson | 75,333,769 | **12,778,116** |
| purple pearson | 675,955,918 | **212,252,625** |
| purple ksg (floor, nats) | 79,795,793 (0.98822) | **3,054,956** (1.29011) |

Purple's MI layer contributes **26x fewer** edges at the stricter cut. That is the
power argument in `thresholds.md` §3.3 showing up directly: at n = 18 the KSG
estimator is noisy, and demanding p <= 3.66e-07 rather than 6.7e-05 removes most
of what it had to offer.

---

## 2026-09-08 — k-NN is chosen by a graph model, not by hand; and the network is Pearson-only

**Why.** `#knn(k)` is what breaks the giant module — purple went 28.0% → 4.2% at
`#knn(160)`, stranding fewer genes than the current setting does. But k was set
by a heuristic, and a hardcoded k is not defensible. k is now chosen per network
as the value whose graph is closest to a Barabási–Albert model under
`statGraph::graph.model.selection`, swept over k = 50…600.

**Pearson only.** The MI layer is 1.14% of sugarcane's edges and 4.20% of
purple's. Dropping it costs sugarcane 1,346 genes (103,336 → 101,990 nodes in
the network) and simplifies the object being reduced. No rebuild was needed:
`02_network_engine.py` already writes the Pearson layer separately, and every
column `03_merge_layers.py` adds for a pearson-source edge is a per-row function
of what is already there, so `40_pearson_only_network.sh` is one streaming awk.
Verified: 75,333,769 edges, exactly the layer summary's count.

### Four things about statGraph that had to be found by measuring

**Its default path cannot run here.** `graph.model.selection` defaults to
`method="diag"`, a full dense eigendecomposition — and not one: 502 grid points ×
50 simulated graphs = **25,100** of them. `eigen()` is cleanly cubic on this
machine (fitted exponent 3.10). At n = 170,736 a single call is **34.6 h and
233 GB** (466 GB peak, i.e. OOM), so the default is **~99 years**.

**`method="fast"` scales, and by a lot.** It derives the spectral density from
the degree distribution instead. Measured, one density on Barabási–Albert graphs
of mean degree 40:

| n | 1,000 | 10,000 | 50,000 | 170,736 |
|---|---|---|---|---|
| unique degrees | 117 | 262 | 452 | 694 |
| seconds | 2.32 | 1.55 | 2.43 | **5.11** |

Cost tracks unique degrees, not node count.

**Its `numCores` path deadlocks.** It builds and tears down a PSOCK cluster per
spectral density; across dozens of GIC calls that leaks connections and hangs —
observed at `numCores=32` with the main process pinned at 0% CPU, fine at 4–16.
statGraph is therefore called **single-threaded**, and all parallelism sits
outside, across k. With 12 k values × 2 studies on 256 cores that is the better
split regardless.

**And a trap that silently corrupts the answer.** statGraph's default ER grid is
`seq(0, 1, 0.01)`. At p = 0 the model graph is empty and the whole call dies with
"missing value where TRUE/FALSE needed"; at p = 1e-06 it returns a spuriously
**negative** GIC that beats every real model — so a graph generated as
Barabási–Albert is reported as Erdős–Rényi. Every parameter range in
`42_knn_model_selection.r` is passed explicitly and excludes the degenerate end,
with ER centred on the observed mean degree rather than spanning [0, 1].

### What the validation gate found, and how to read the result

`39_statgraph_validate.r` runs before the criterion is allowed to choose
anything. On synthetic graphs of known type:

| truth | n = 1,000 | n = 10,000 |
|---|---|---|
| BA | **BA** ✓ | **BA** ✓ |
| ER | ER ✓ | WS ✗ |
| WS | ER ✗ | ER ✗ |

**BA is recognised every time it is the truth; ER and WS are confused in both
directions.** That split is the whole verdict. Choosing k by argmin GIC(BA) needs
only the first property, so the k choice is licensed — but `selected_model`
should be read as *BA vs not-BA*, and which of ER/WS placed second must not be
quoted.

Against the exact method on real co-expression subgraphs (snowball-sampled, so
the sample keeps the density and clustering that make the check meaningful),
`fast` and `diag` **agree**. An earlier apparent disagreement came from sampling
600 nodes uniformly out of 98,297, which induces a 133-node, 78-edge fragment
with none of the structure under test — a reminder that a subgraph of a network
is not a small version of it.

Both criteria are computed for every k regardless: `gic_ba`, and a
Clauset–Shalizi–Newman power-law fit (`pl_alpha`, `pl_ks_stat`) as the agreed
fallback. Computing both costs nothing and means a failed gate needs no re-run.

### Cost, measured

One k on the real Pearson-only sugarcane graph at k = 200 (96,335 nodes,
1,724,999 edges): **939 s**, selecting BA with a clear margin — GIC 0.990 against
2.089 (ER) and 2.347 (WS), BA exponent 1.016. Twelve k values run concurrently,
so a study is one batch, not twelve.

**Nothing is adopted.** `MCL_KNN_*` stay empty and the main tree is untouched;
this lives in `results_pearson/` on branch `pearson-knn-ba`.


## 2026-09-12 — The gene level moves to the network's node universe, blocked only, nodes only

Four decisions, taken together so the gene level describes the same graph the module
level does.

**The gene universe is the Pearson-only network's node set** — 101,990 sugarcane and
170,135 purple — read from `network_<study>_node_metrics.tsv` via `node_list()` and
`restrict_to_universe()`. The alternatives were every VST gene (170,790 / 170,740) and
the old `conserved_genes_<study>_FULL.txt` (39,226 / 44,118). The old one was rejected
because it is the genes carrying a conserved edge in the **merged** graph, which this
analysis no longer uses; the whole VST matrix was rejected because it tests genes the
network analysis never sees. The consequence is explicit: the BH denominator rose 2.6x
and 3.9x, so the responsive counts are **not** comparable with the ones the previous
blocked table held, and the write-up says so rather than presenting 8,737 as an
improvement on 4,137.

The universes are asymmetric — purple's network holds 99.6% of its VST genes,
sugarcane's 59.7% — and that asymmetry is a property of the data at |r| >= 0.8, not a
choice. It is why a fourth status (`ortholog_not_a_node`) had to be added: a gene can
have a perfectly good ortholog that is not in the other species' graph.

**Purple's primary gene statistic is blocked Spearman on midranks.** Its trait is an
ordinal 0/2/6 mM dose and Pearson reads that spacing literally. This is the argument
already accepted at module level (`19`, `53`) and recorded in `Open items` as never
having been applied at gene level; it is applied now. Sugarcane keeps Pearson because
a two-level trait makes the two the same test up to a monotone relabelling. Both are
computed and written in either case, so the choice is visible in the table.

**The non-monotone tier is included for purple, as a second test family.** Each family
is BH-corrected within itself over the *same* genes and the selections are unioned;
`61` aborts if the two gene sets differ, because a union across two denominators is
not a union. It contributes **one** gene, whose ortholog is not sugarcane-responsive,
so it adds no conserved pair. Including it anyway was the right call: the question
"how much does the monotone-only rule cost in a three-level design" now has a measured
answer at gene level instead of an untested assumption, and the answer is ~nothing that
survives FDR against a 2.06x aggregate excess that does not localise.

Sugarcane's two-level design cannot express curvature at all, so the tier is purple-only
and the asymmetry is declared: an admitted pair says "purple responds non-monotonically
and its sugarcane ortholog responds monotonically", never "both species share a
non-monotone response". That question is unanswerable with these two designs.

**Node level only.** The edge level would need `conserved_edges_*_FULL.tsv`, which
describes the merged graph. `08_conserved_cor_genes.r` is left untouched as the record
of that analysis rather than being edited into something half-rebuilt, and `61` is a
new script beside it.

### Two reporting rules that came out of the numbers

**The pair count is never quoted without the gene count.** 392 conserved pairs beat
their ortholog-shuffle null 1.25x (p 0.000999), but the 326 *genes* behind them beat it
only 1.08x (p 0.071), and under the directed design 0.94x — below the null. Orthology is
many-to-many (`OG0001395`: 6 pairs from 3 x 2 genes), so most of the pair excess is
multiplicity. Reporting the significant row alone would be a real overstatement of a
result that does not hold at gene level.

**Sign concordance is tested at orthogroup level, not only over pairs.** The pair-level
binomial treats 392 non-independent trials as independent — the same pseudoreplication
sugarcane's plant-level control exists to answer at gene level. One vote per orthogroup,
majority sign, ties abstaining: 59.5% of 269 at p 0.0022 genome-wide, 57.7% of 418 at
p 0.0020 directed, against 62.2% / p 1.4e-06 over pairs. The claim survives the
correction, weakened as it should be, and both are reported.

### One thing blocking did not fix

Purple's blocked p-value histogram still is not flat: the last decile falls from 12.5%
to 10.2% but the shape dips to 5.4% and rises again. Genotype and nitrogen do not
account for everything at n = 18. The script warns above 11% and this sits below it, so
the warning does not fire — which is exactly why it is written down here instead.


## 2026-09-13 — One GO source, no merge, no tier; and PANNZER is out

The plan this closes was four annotation sources, tiered by evidence and merged, each
adopted only if it improved module coherence. Three were built and scored. The
decision is to ship **one**: full InterProScan 5.78, all 17 member databases, through
the pinned interpro2go/pfam2go and normalised to most-specific terms.

**The merge was rejected by its own judge, which is why this is a result and not a
simplification.** On a fixed gene set -- identical genes, identical modules, so only
the annotation varies -- the union of InterPro and eggNOG scores +0.0671 in sugarcane
against +0.0748 for InterPro alone, and +0.0264 against +0.0361 in purple. Merging
lowers coherence in both species. A tiered table would therefore have cost a `source`
and `tier` column on every pair, a provenance caveat on every downstream result, and
bought negative signal. The user's instinct that merging was "unnecessary complexity"
is what the measurement says too.

**Two things the numbers say that the headline does not.** eggNOG's raw homogeneity is
2.6x InterPro's while its excess over the null is LOWER -- without a size-matched null
it would have been adopted on sight, and what its raw H measures is term density (16.7
terms per gene after normalisation, against 2.1), not shared function. And the 5-DB
and 17-DB InterPro tables are TIED on coherence, 0.9% apart and inside the noise: the
17-DB scan was adopted for REACH, not because it is the more coherent annotation. It
annotates 7,204 more sugarcane and 15,094 more purple network genes, which is what
turns into 37 more testable sugarcane modules and 21 more terms clearing BH. The
write-up says this rather than "InterPro full is the best annotation", which the data
does not support.

**MF and CC are now run on the adopted table**, and MF turns out to be the stronger
ontology: 272 terms clear cross-module BH against BP's 164 in sugarcane, 17 against 9
in purple. That is the expected shape for a domain-derived annotation -- an InterPro
signature names what a protein does far more sharply than what process it is in. BP
stays primary anyway, because the question is about a nitrogen response and that is a
process; the choice is now made knowing it costs power, instead of by default. CC is
weak in both species and near-empty in purple (one term) and nothing should rest on it.

**PANNZER2 is out**, by decision, and its 6.1 GB of output was deleted. It was also
never finished: purple completed (242/242 chunks, 5.2M predictions over 170,151
proteins) but sugarcane stopped at 174 of 195, with 21 chunks killed by ConnectTimeout
to the public SANS service at Helsinki. The runner and the commit that hardened it
(b283178 -- quarantine partial chunks, bounded retry passes) stay in the repo, so
re-running is hours rather than reconstruction, but it depends on a service that is
not ours. It was never scored against the coherence judge and no claim is made about
whether it would have helped.

The other three tables are kept and suffixed `.notused`, with `annotation/README.md`
recording which table is live and one line per source saying why it is not. They are
the measurement that justifies the single-source decision, and deleting them would
leave the decision unsupported.

**What this pass does NOT touch.** `09_go_enrichment.r` and `10_go_semantic.r` still
run on the original 7.7% eggNOG GO column and its background. They are conserved-set
and semantic-clustering stages, not module stages, and re-running them is a separate
piece of work; they stay marked stale.
