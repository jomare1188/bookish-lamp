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
cannot be settled on sugarcane alone. **The purple fit is running on another
machine as of 2026-08-27**; until it lands the comparison covers one species, and
a change that would regenerate every module-level figure should not rest on one
species.

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
