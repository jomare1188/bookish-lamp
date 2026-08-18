# Results

Numbers from the current build. Update this file as stages complete — it is the
one place where counts live, so that [thresholds.md](thresholds.md) and
[methods.md](methods.md) can stay stable.

Build started 2026-08-13.

---

## Inputs

| study | dds | libraries | genes before filter | genes at CV ≥ 15 |
|---|---|---|---|---|
| sugarcane | `run1/salmon/deseq2_qc` | 48 | 190,973 | **170,790** |
| purple | `china/run2_onlyL/salmon/deseq2_qc` | 18 | 215,183 | **170,740** |

Purple is now the leaf-only quantification — see [decisions.md](decisions.md).
Its previous gene count, from `china/run1` + `Group1=="L"`, was 170,736.

The sugarcane VST export is **byte-identical** (same md5) to the one the previous
MI sweep used, so that side of the rebuild is a controlled comparison: only the
method changed, not the data.

---

## Layer construction

| study | layer | edges | candidate edges | floor | null | runtime |
|---|---|---|---|---|---|---|
| sugarcane | pearson | **75,333,769** | 213,790,679 | \|r\| = 0.8 | analytic | 5.5 min |
| sugarcane | ksg | **2,871,152** | 24,190,102 | 0.98911 nats | 2e8 perm | 169.7 min |
| purple | pearson | **675,955,918** | 1,364,134,888 | \|r\| = 0.8 | analytic | 58.7 min |
| purple | ksg | **79,795,793** | | 0.98822 nats | 1e6 perm | 45.5 min |

The two KSG floors — 0.98911 at n = 48 and 0.98822 at n = 18 — sit five decimal
places from the previous build's 0.98910 / 0.98827, despite purple having a
different VST, a different gene set and an independently drawn null. The floor is
fixed by the target p-value and n, both unchanged, so the input change moves it
only in the fifth decimal. That stability is the matching argument holding.

Purple's null needed only 1e6 permutations against sugarcane's 2e8: at n = 18 the
target p is 6.72e-5, which 1e6 draws resolve directly, so the auto-sizer does not
reach for the ceiling. It also reports **0 edges beyond the fitted null's upper
endpoint**, i.e. every purple MI p-value is calibrated rather than floored;
sugarcane has 69,390 such edges (2.4%), the price of extrapolating to 9e-12.

Measured KSG null means, read by the merge from each layer's own `null.tsv`
rather than a lookup table: **+0.12428** (n = 48) and **+0.17968** (n = 18),
against the +0.124 / +0.180 that were previously hardcoded.

The KSG floor is the MI value with the same per-edge false-positive rate as
\|r\| = 0.8 at n = 48 (p = 9.02e-12) — see [thresholds.md](thresholds.md).

At n = 18 the engine auto-raises the candidate cut from the requested
p = 1.219e-3 (\|r\| = 0.7) to 1.345e-3 (\|r\| = 0.696), so the 0.8 floor clears the
candidate boundary by the required 20×. BH then stays exact up to 392 M
rejections, far above what purple needs.

All 1.46 × 10¹⁰ sugarcane pairs were swept in **12 seconds** (1.62 × 10⁹
pairs/s); the rest of the 5.5 min is the BH pass and writing 6.1 GB of text.

---

## Validation of the linear layer

Sugarcane's input is unchanged from the previous build, so the new linear layer
is directly comparable to the old `network_sugarcane_filtered_edges.tsv`.

Ground truth: r computed in **float64** from the DESeq2 VST for the first 6,000
genes, thresholded at 0.8 ≤ |r| ≤ 0.9999.

| | edges | vs truth |
|---|---|---|
| float64 truth | 73,404 | — |
| **new layer** | **73,404** | **0 missing, 0 extra** |
| old network | 73,453 | 0 missing, **49 extra** |

The new layer reproduces float64 ground truth exactly. Across the full network it
has 75,333,769 edges against the old 75,380,961 — **the old network contained
47,192 edges whose true \|r\| was below its own 0.8 threshold**, admitted because
the dense correlation matrix stored r rounded to 4 decimal places. See
[decisions.md](decisions.md) for the full diagnosis.

Independent checks that passed:

- `./run.sh validate` — all checks, including GPU r vs `numpy.corrcoef`
  (max 2.4e-07), analytic p vs `scipy.stats.t` (exact), sign preservation for
  negative correlations, and `--match-pearson 0.8` round-tripping to exactly 0.8.
- The KSG path is **bit-identical** to the original engine on a 2,000-gene
  subset — edge list and null file both compare equal — confirming the three
  estimator hooks are true no-ops for the MI estimators.
- `lib/common.R:read_vst()` round-trips the DESeq2 VST to 9.5e-07 (pure float32).
- Measured KSG null at n = 48: mean **+0.12428**, sd 0.10004 — matching the
  +0.124 used for the `r_eq` bias correction.

---

## The network

| study | total edges | pearson only | both | mi only | size |
|---|---|---|---|---|---|
| sugarcane | **76,200,344** | 73,329,192 (96.23%) | 2,004,577 (2.63%) | **866,575 (1.14%)** | 8.27 GB |
| purple | **705,571,723** | 625,775,930 (88.69%) | 50,179,988 (7.11%) | **29,615,805 (4.20%)** | ~65 GB |

**30.2% of the sugarcane MI edges and 37.1% of the purple MI edges are invisible
to the linear layer** — 866,575 and 29,615,805 edges that a correlation cannot
represent. That is the return on the method.

Do **not** read purple's larger MI share as more non-linear biology. At n = 18 the
p-value implied by \|r\| = 0.8 is 6.72e-5 against 9.02e-12 at n = 48, so purple's
matched floor is far softer in power terms. The asymmetry is sample size, not
species. See [thresholds.md](thresholds.md) §3.3.

### Sugarcane vs the previous build

Same input (the VST export is md5-identical), so this is a controlled comparison
of the *method* alone:

| | previous | current | change |
|---|---|---|---|
| pearson only | 73,375,808 | 73,329,192 | **−46,616** |
| both | 2,005,153 | 2,004,577 | −576 |
| mi only | 866,148 | 866,575 | **+427** |
| total | 76,247,109 | 76,200,344 | −46,765 |

The story is internally consistent: the linear layer loses ~46.6 k edges — the
ones admitted by the old 4-decimal rounding — and the MI-only count goes *up* by
427, because some edges previously counted as `both` had a spurious Pearson
partner and are now correctly attributed to the MI layer alone.

The MI layer itself changed by only **149 edges** (2,871,301 → 2,871,152,
0.005%). Same seed and the same 2 × 10⁸ permutations, so the null is identical;
the two differ only in how the floor is reached — the old file was an `awk`
filter on `pval`, the new one inverts the GPD fit at the same p. Agreement to
0.005% is an independent check that baking the threshold into the sweep
reproduces the post-hoc filter.

---

## Topology

| study | nodes | edges | components | giant component | density |
|---|---|---|---|---|---|
| sugarcane | **103,336** | 76,200,344 | 939 | 101,253 (98.0%) | 1.43e-2 |
| purple | **170,736** | 705,571,723 | 1 | 170,736 (100.0%) | 0.04840871 |

The node count is the load-bearing number: against the Pearson-only baselines of
102,020 and 170,103, the MI layer contributes ~1,300 (sugarcane) and ~630
(purple) genes that have **no linear edge at all** at |r| >= 0.8. Purple's
network now covers 170,736 of the 170,740 genes in its VST — only four genes are
isolated. Those change the GO enrichment denominator and the
degree-matched permutation nulls in the MYB61 and Module-20 readouts, which is
why `node_metrics` is regenerated rather than reused.

`igraph::simplify()` in the clustering stage reported the same 103,336 / 76,200,344
after removing multi-edges and self-loops — i.e. it removed none. That is a check
on the merge's edge key (`min(i,j) * n_genes + max(i,j)`): a collision there
would have surfaced as a collapsed edge count.

**Purple is one connected component.** All 170,736 nodes, mean degree ~8,264,
density 0.0484 against sugarcane's 0.0143. The Pearson-only baseline already had
44 components with a 170,012-node giant, so this is not created by the MI layer —
625.8 M of the 705.6 M edges are Pearson-only — but it is the single most
important caveat on every purple result downstream:

- It produces a top-heavy partition: one module with 28% of the genes and 14,594
  genes in modules too small to name. **It does not produce worse modularity** —
  purple's Q of 0.1631 is higher than sugarcane's 0.1005. Density and modularity
  are not the same axis, and an earlier note here predicting otherwise was
  wrong: Q is measured against a degree-preserving null, so a denser graph is not
  penalised for being dense.
- The threshold is soft here. At n = 18. At that sample size \|r\| >= 0.8 is p = 6.72e-5 — six
  orders of magnitude softer than the same \|r\| at n = 48 — so the threshold
  admits far more of the correlation distribution. It is FDR-honest (BH over
  1.46e10 tests) but not selective.
- The design compounds it: 2 genotypes x 3 nitrogen levels x 3 replicates means
  most genes share a coarse genotype-or-treatment response, which is real
  covariation but not specific co-regulation.

If purple's modules turn out to be uninformative, the fix is a stricter threshold
for purple specifically — not a different clustering algorithm. That would break
the equal-specificity matching with sugarcane, so it is a deliberate trade to
make and document, not a tuning knob.

Local transitivity is NA by design (`COMPUTE_TRANSITIVITY=0`).

---

## Modules

MCL, inflation 2, `min_module_size = 2` (the pipeline's historical value).

| study | modules | largest | median size | modularity Q | unassigned |
|---|---|---|---|---|---|
| sugarcane | **10,309** | 19,604 | 3 | 0.1005 | 401 |
| purple | **9,881** | 47,887 | 3 | 0.1631 | 14,594 |

Against the Pearson-only baselines:

| | sugarcane old → new | purple old → new |
|---|---|---|
| modules | 8,691 → 10,309 | 7,464 → 9,881 |
| largest | — → 19,604 | 44,050 → 47,887 |
| Q | 0.1036 → 0.1005 | 0.1726 → 0.1631 |
| unassigned | 105 → 401 | 273 → 14,594 |

Both networks yield **more, smaller modules** after augmentation, with modularity
essentially unchanged (down ~0.003 and ~0.010). The extra edges fragment the
partition slightly rather than reorganising it.

Purple's first module holds 47,887 genes — 28% of the network — and 14,594 genes
land in modules too small to name. Sugarcane's largest is 19,604 (19%) with only
401 unassigned. Purple's partition is the more top-heavy of the two, which is
consistent with its density, though it is *not* the less modular one by Q.

`min_module_size` was briefly run at 5 during this build. It does **not** affect
the clustering — the raw partition (10,710 modules) and Q (0.1005) came out
identical either way; it only decides how small a module may be and still get a
`Module_NNN` name rather than "Unassigned" (3,307 named / 18,935 unassigned at 5,
versus 10,309 / 401 at 2). Reverted to 2 to match the original pipeline.

---

## Conservation

Each direction streams one network and looks its edges up in the other. **The two
directions therefore test different MI layers** — `sugarcane_to_purple` tests
sugarcane's 866,575 MI-only edges, `purple_to_sugarcane` tests purple's
29,615,805. They are not two measurements of one quantity.

### sugarcane → purple

| layer | edges | conserved | fraction | vs pearson |
|---|---|---|---|---|
| ALL | 76,200,344 | 8,197,303 | 10.758% | — |
| pearson only | 73,329,192 | 7,836,068 | 10.686% | 1.00 |
| both | 2,004,577 | 257,564 | 12.849% | **1.20** |
| mi only | 866,575 | 103,671 | 11.963% | **1.12** |

39,226 sugarcane genes sit on at least one conserved edge.

### purple → sugarcane

| layer | edges | conserved | fraction | vs pearson |
|---|---|---|---|---|
| ALL | 705,571,723 | 10,722,788 | 1.520% | — |
| pearson only | 625,775,930 | 9,393,280 | 1.501% | 1.00 |
| both | 50,179,988 | 909,550 | 1.813% | **1.21** |
| mi only | 29,615,805 | 419,958 | 1.418% | **0.95** |

The overall rates (10.8% vs 1.5%) are not comparable across directions: purple has
9x more edges to test against a sparser target, so a purple edge is far less
likely to find an ortholog partner. Only the **within-direction** contrast between
layers is meaningful.

### Is any of it above chance? — the permutation null

`13_conservation_null.r` permutes the target column of the ortholog table,
destroying the true assignment while preserving fan-out, coverage and aggregate
exposure to the target network's degree distribution. 20 replicates, 5 M sampled
source edges.

| direction | layer | observed | null | **fold** | z |
|---|---|---|---|---|---|
| sc → pu | pearson | 10.6902% | 4.1668% | **2.57** | 78 |
| | both | 12.7968% | 4.6337% | **2.76** | 75 |
| | mi | 12.0698% | 4.6796% | **2.58** | 63 |
| pu → sc | pearson | 1.4996% | 0.6086% | **2.46** | 68 |
| | both | 1.8398% | 0.6232% | **2.95** | 55 |
| | mi | 1.4135% | 0.6230% | **2.27** | 41 |

**Edge conservation is real.** Both directions run ~2.5x above a null with the
same orthology structure, and no replicate came close (p = 1/21, the floor at 20
replicates). The absolute rates mean something after all.

The null also explains the asymmetry: the null rate is 4.17% one way and 0.61%
the other — a 6.8x difference that tracks the density gap almost exactly. The
*fold over null* is 2.57 and 2.46, i.e. **the two directions agree once the
denominator is removed.** That is the symmetric quantity; the raw rates are not.

### The MI advantage does not survive the null

This corrects the earlier reading of the raw rates.

| | raw rate vs pearson | **fold-over-null vs pearson** |
|---|---|---|
| sc → pu, `both` | 1.197 | **1.076** |
| sc → pu, **`mi`** | **1.129** | **1.005** |
| pu → sc, `both` | 1.227 | **1.198** |
| pu → sc, **`mi`** | **0.943** | **0.921** |

MI-only edges have a *higher null rate* than Pearson-only edges (4.68% vs 4.17%;
0.623% vs 0.609%) — they connect genes with more orthologs and better-connected
partners, so they had more chances to match by accident. Once that opportunity is
removed, **MI-only edges are conserved at essentially the Pearson rate** (1.005)
in sugarcane and below it (0.921) in purple.

This resolves the puzzle in the raw numbers rather than deepening it. The raw
1.129 and 0.943 looked like a contradiction between directions; corrected to
1.005 and 0.921 they are consistent, and both say the same thing: **the MI layer
finds edges that are no better conserved than the linear layer's, and the earlier
"12% better" was opportunity bias, not biology.**

What *does* survive is the `both` layer: 1.076 and 1.198 over null, elevated in
both directions. Edges that two independent estimators agree on are genuinely
better conserved. That was the most robust result before the null and it remains
the most robust one after.

---

## Gene–trait correlation and node-level conservation

Per-gene expression vs trait, restricted to genes on a conserved edge.

| study | genes tested | padj ≤ 0.05 (treatment) | selected (also \|r\| ≥ 0.6) |
|---|---|---|---|
| sugarcane | 39,226 | 6,302 | **1,361** |
| purple | 44,118 | **62** | **62** |

**The two selections are not defined by the same constraint**, which has to be
said before any comparison of them:

| | sugarcane | purple |
|---|---|---|
| n / df | 48 / 46 | 18 / 16 |
| \|r\| at the BH boundary | 0.378 | **0.799** |
| \|r\| needed by the single best gene | 0.635 | 0.884 |
| what actually binds | the \|r\| ≥ 0.6 effect-size cut | the FDR cut |

At n = 18 with BH over 44,118 genes, a gene needs \|r\| ≈ 0.80 just to clear the
FDR, so `TRAIT_R_THR = 0.6` never binds on purple and every one of its 62 genes
is far above it. On sugarcane the FDR admits \|r\| ≥ 0.378 and the 0.6 cut does
the work. Same nominal thresholds, different operative ones. Verified there is
no encoding bug: both traits used all 18 / 48 samples, no sample was dropped as
unencoded, and the samplesheet's `treatment` values (`0N`/`2N`/`6N`) match the
config exactly.

### Conservation of the nitrogen response — node and edge level

`08_conserved_cor_genes.r` takes a `selection` argument (`pearson` | `mi` |
`union`, default `union`) and reads both statistics from
`gene_trait_mi_<study>.tsv`, so the two are computed on the same samples and the
same VST and the rules are directly comparable. Genes are first restricted to the
conserved-edge gene set, which is built from conserved edges of **all** layers.

| | Pearson | MI | union |
|---|---|---|---|
| responsive sugarcane genes | 1,361 | 3,221 | 3,265 |
| responsive purple genes | **30** | **5** | **32** |
| **conserved correlated ortholog pairs** | **1** | **0** | **2** |
| expected by chance | ~2.1 | ~0.8 | ~5.5 |

**Node level: at or below chance under every rule.** 1,360 of sugarcane's 1,361
fall in `ortholog_not_correlated`, not `no_ortholog`, so it is not an orthology-
coverage artifact — the orthologs exist and are simply not responsive on the
other side.

### Edge level — the question the MI layer was built for

A conserved edge counts only if **both** its genes are responsive in **both**
species. On the sugarcane side alone there are plenty, and the MI layer supplies a
real share of them:

| selection | conserved edges, both endpoints responsive in sugarcane | `pearson` | `both` | **`mi`** |
|---|---|---|---|---|
| Pearson | 5,894 | 4,761 | 1,065 | **68** |
| MI | 34,898 | 25,507 | 6,801 | **2,590** |
| union | 35,761 | 26,170 | 6,993 | **2,598** |

**Both species: 0, under every selection rule.**

The reason is arithmetic rather than a weak signal. An edge needs *two* genes
responsive on both sides, and there are 1 (Pearson) or 2 (union) in the entire
analysis — `SoffiXsponR570.10os1g012500`, plus `SoffiXsponR570.05Cg220600` when MI
is included — and they are not connected to each other.

So the funnel closes at **purple's 30–32 responsive genes**, not at the edges and
not at the MI layer. Adding MI moves the sugarcane side from 5,894 to 35,761
qualifying edges (2,598 of them invisible to Pearson) and moves the purple side
*down*, 30 → 5. The n = 18 selection is the binding constraint everywhere.

This is worth stating as a result rather than an absence: **the MI layer did
supply non-linear conserved edges between nitrogen-responsive genes — 2,598 of
them — in the study that has the samples to find them.** What it could not do is
supply the matching purple side.

---

### The directed test, and what it shows

Correcting purple genome-wide (m = 44,118) is the wrong burden for a comparative
question. `08_conserved_cor_genes.r` therefore takes a **directed** mode
(`./run.sh conscor 1`): purple's p-values are corrected over only the orthologs of
sugarcane's responsive genes, which is the hypothesis actually being tested.

| | |
|---|---|
| candidate purple orthologs | 4,745 of 44,118 |
| BH denominator | 44,118 → 4,745 |
| \|r\| required of the best gene | ≈0.884 → ≈0.844 |
| **passing, same candidates, genome-wide BH** | **2** |
| **passing, same candidates, directed BH** | **1** |

**It found fewer, not more.** BH is adaptive — a gene benefits from sitting in a
set containing many discoveries — and the candidate set is depleted of purple's
strongest signal, so the smaller denominator does not compensate.

That depletion is the result:

| | |
|---|---|
| purple's genome-wide responsive genes | 32 |
| …that are orthologs of a sugarcane-responsive gene | **2** |
| purple's top 30 genes by p, that are candidates | **2** (≈3.2 expected if independent) |

**The two species' nitrogen-responsive gene sets are independent in ortholog
space, marginally below chance.** This is a substantive negative result, not a
missing-power caveat.

Nor is the signal merely below a threshold. Among all 4,745 candidates:

| best \|r\| | required | \|r\| ≥ 0.8 | \|r\| ≥ 0.6 |
|---|---|---|---|
| 0.854 | 0.844 | **3** | 128 |

Three genes clear 0.8. An edge needs two *connected* genes responsive on both
sides, so no denominator choice produces one. Reported uncorrected (raw p ≤ 0.05,
\|r\| ≥ 0.6) the candidate count is 128 — the most generous available reading, and
it should carry that label.

**The zero survives three selection rules, two correction burdens, and both
conservation directions.**

---

## Gene–trait mutual information (`12_gene_trait_mi.py`)

Added to attack the 62-gene purple bottleneck: the trait is discrete, so this
uses the **Ross (2014)** estimator rather than KSG, with a label-permutation null
that — because every gene is rank-transformed — is shared across all genes and
therefore cheap to make deep (10⁷ permutations, resolving p to 1e-7 empirically).

| study | n | trait | MI sig. | Pearson sig. | both | **mi only** | pearson only |
|---|---|---|---|---|---|---|---|
| sugarcane | 48 | binary High/Low | 9,400 | 20,331 | 7,103 | **2,297** | 13,228 |
| purple | 18 | 3-level dose | 21 | 85 | 14 | **7** | 71 |

*(over all genes, BH across the whole transcriptome; `07_gene_trait_cor.r`
restricts to conserved genes first, hence its different counts)*

### It is correctly calibrated

- **Negative control:** re-running sugarcane with the trait labels shuffled gives
  **0 significant genes by MI and 0 by Pearson**. The FDR is controlled.
- Validation suite: the batched kernel matches a naive loop reference to 5.6e-07;
  perfect class separation approaches H(trait) from below (0.657 vs ln 2; 0.908
  vs ln 3); p-values under H₀ are uniform (0.0495 / 0.0096 / 0.00098 against
  nominal 0.05 / 0.01 / 0.001).
- 656 of 170,790 sugarcane genes exceed the empirical null's maximum and get
  GPD-extrapolated p-values (97 of those are floored). Those p-values **rank**
  correctly but their magnitude is not calibrated — a reported padj of 1e-20
  means "p < 1e-7", nothing more. The output table flags each row's
  `p_source` as `empirical` or `extrapolated`.

### But it does not find what it was added to find

**On sugarcane it adds real signal — of a different kind than advertised.** The
2,297 MI-only genes are *not* non-linear dose responses:

| property of the 2,297 MI-only genes | |
|---|---|
| significant by Wilcoxon (a location shift) | **2.4%** |
| variance ratio between N groups > 2× | 52.4% |
| median \|AUC − 0.5\| (group separation) | 0.146 |
| median share of variance explained by N | 5.4% (vs 20.4% for Pearson-only) |

By contrast **100%** of Pearson-only genes are Wilcoxon-significant: those are
clean location shifts. So MI is picking up genes whose expression *distribution*
differs by nitrogen — dispersion and shape — rather than genes whose mean shifts.
That is real (the negative control rules out noise) and arguably interesting, but
it is a broader phenomenon than the saturating dose response that motivated the
work, and it should not be described as "non-linear nitrogen response" without
qualification.

Checked and excluded: this is **not** a genotype × nitrogen interaction. A
two-way ANOVA gives the interaction 0.7% of the variance in MI-only, Pearson-only
and non-significant genes alike.

**On purple it does not help at all.** 21 significant genes against Pearson's 85,
adding 7. At n = 18 with three classes of six, the MI null has sd 0.130 and a
maximum of 0.908 against a ceiling of ln 3 = 1.099, so significance requires
near-perfect three-way separation — a *higher* bar than the \|r\| ≈ 0.8 Pearson
needs. An omnibus test spends its power detecting any of the ways three classes
might differ; when the effect is monotone, the directed test wins.

### The bottom line: it does not change the biological conclusion

Redoing the node-level conservation with each selection rule, restricted to the
conserved-edge genes as `08_conserved_cor_genes.r` does:

| selection | sugarcane genes | purple genes | conserved correlated pairs | expected by chance |
|---|---|---|---|---|
| Pearson (as run) | 1,361 | 30 | 1 | 2.1 |
| MI | 3,221 | 5 | 0 | 0.8 |
| **union** | 3,265 | 32 | **2** | **5.5** |

Two observed against ~5.5 expected. **Still no evidence of a shared node-level
nitrogen response**, and the purple side is still the bottleneck — MI made it
*smaller*, not larger.

This was proposed as the highest-value next step. It was worth doing — the
sugarcane signal is real and the machinery is validated and reusable — but as a
fix for the purple bottleneck it failed, and the honest reading is that **n = 18
is the binding constraint, not the choice of statistic.** No estimator recovers
information that 18 libraries do not contain.

---

## Module-level nitrogen response

One eigengene per MCL module (PC1 of member genes' VST, z-scored, oriented to
mean expression), for modules with ≥ 3 genes: 6,576 sugarcane / 6,318 purple,
covering 92% and 87% of each network's genes. Median PC1 variance explained
**81.1%** (sugarcane) and **65.0%** (purple).

Each eigengene is tested by **both** statistics, by running `12_gene_trait_mi.py`
unchanged on the eigengene matrix — the same validated code path as the gene
level, so the two answers cannot diverge for methodological reasons.

### Effect-size floors, and why both sides need one

`padj <= 0.05` alone admits sugarcane modules down to **|r| = 0.356** (13% of the
eigengene's variance) and is not comparable to the gene-level selection, which
required `|r| >= 0.6`. Applying that floor to the linear side **only** is worse
than not applying it: it moves every modestly-linear module out of `both` and
into `mi_only`, which then reads as "non-linear" when it is nothing of the kind.
Measured: `mi_only` went 253 → **786**.

So MI needs an equivalent floor. The nats-to-|r| identity the network layer uses
does **not** apply — it assumes two continuous variables, and here the trait is
discrete with MI capped at H(trait). The floor is instead **calibrated on this
data**: among the 192 sugarcane modules sitting at |r| ≈ 0.60, the median MI is
**0.218 nats (0.31 of H)**. That is the MI equivalent of the linear cut at this
n and class structure. A calibration, not a theoretical equivalence, and the
number is printed on every run.

| sugarcane | padj only | **both floors** |
|---|---|---|
| neither | 4,614 | 5,929 |
| both | 939 | **367** |
| mi_only | 253 | **239** |
| pearson_only | 770 | **41** |

**647 of 6,576 modules respond; 239 only non-linearly.** Permuting the trait
labels against the same eigengenes gives **0** significant modules.

### TF enrichment is weak once the floors are applied

| class | modules | TF-enriched | OR | Fisher p |
|---|---|---|---|---|
| pearson_only | 41 | 7.32% | 6.16 | 0.016 |
| both | 367 | 2.45% | 1.96 | 0.061 |
| mi_only | 239 | 1.67% | 1.33 | 0.55 |

An earlier unfloored run had `both` at OR 3.23, p = 4.7e-06, and that no longer
holds. The strongest enrichment is now in `pearson_only`, but that is **3 of 41
modules**, `both` is marginal, and `mi_only` is null. The honest reading: TF
enrichment among responsive modules is weak and the classes are too small to
separate confidently. Do not report "TFs concentrate in the modules both
statistics agree on" — that was an artifact of the missing effect-size floor.

### A third of the `mi_only` modules are sample-driven

Visible in `module_summary_sugarcane.png` as vertical streaks in the `mi_only`
block, and quantified:

| class | median top-2 sample share of variance | R² of the High/Low split | 2 samples carry >25% |
|---|---|---|---|
| both | 0.182 | **0.484** | 14% |
| pearson_only | 0.185 | 0.384 | — |
| **mi_only** | **0.218** | **0.229** | **35%** |

(under a flat null one of 48 samples carries 0.021)

`mi_only` modules explain **less than half** as much of their variance by the
nitrogen split as `both` modules do, and **35% have two samples carrying over a
quarter of the eigengene's variance**. An omnibus dependence test will fire on a
distributional quirk in one or two libraries, and some of these are that.

It is a minority — 65% are not obviously sample-driven — but the class should be
filtered on this diagnostic before any individual `mi_only` module is followed
up, not taken as a list of non-linear responders.

### Purple (n = 18)

**38 responsive: 37 `pearson_only`, 1 `both`, 0 `mi_only`**, unchanged by the
floors (all 38 already cleared |r| = 0.6 — FDR binds that hard at n = 18). No TF
enrichment, 0 of 37.

### Figures

- `module_summary_<study>.{png,pdf}` — every responsive module's eigengene in one
  panel, split by response class, with PC1 variance explained and module size as
  row annotations. This is the figure that made the sample-driven `mi_only`
  problem visible. Purple's shows textbook monotone dose responses across
  0N/2N/6N, splitting cleanly into modules that fall and rise with nitrogen.
  Sample labels are dropped (the design is carried by the annotation bars), which
  lets the cells shrink to 2 mm.

  **Columns group replicates.** Sorting by sample name alone interleaves
  sugarcane's three leaf positions -- the names run `B0_1, B_1, M_1, P_1, B0_2,
  ...`, so replicates of one tissue sit four columns apart and the tissue effect
  reads as vertical striping right across the figure. Ordering by treatment,
  genotype, tissue, then sample (`HEATMAP_GROUP_BY`) puts replicates adjacent and
  the treatment blocks resolve.
- `heatmaps/module_<id>_<study>.{png,pdf}` — per-module gene-level z-scores, top
  20 per response class. Height scales at 30.5 px per gene (r = 0.9998), so a
  3-gene and a 100-gene module are both legible.

### What the responsive modules are for — GO BP, one test per module

One topGO BP enrichment per responsive module, response classes pooled, against
the same network-node background the gene-level GO and the TF hypergeometric use
(`./run.sh modulego <study>`).

**The annotation gate is the binding constraint, and it is severe.**

| | sugarcane | purple |
|---|---|---|
| responsive modules | 647 | 38 |
| GO-annotated network nodes (background) | 8,251 of 103,336 | 12,255 of 170,736 |
| median GO-annotated members per responsive module | **0** | **0** |
| **modules testable (≥ 3 annotated members)** | **71 (11%)** | **3 (8%)** |
| modules with ≥ 1 enriched term | 65 of 71 | 3 of 3 |
| terms written | 359 | 20 |
| terms clearing cross-module BH ≤ 0.05 | 91 | 5 |

**89% of the responsive modules cannot be tested at all.** Two facts multiply:
only 8% of sugarcane's network nodes carry any eggNOG GO annotation, and the
median responsive module holds 5 genes. **401 of the 647 (62%) have *zero*
annotated members**; purple's figure is 27 of 38. So this section describes 71 modules, not 647, and the 71 are biased
toward the large ones — nothing here should be read as characterising the
responsive set as a whole.

Within that 11%, the signal is strong and it is about nitrogen:

| module | class | genes | ann. | top BP term | p | global BH |
|---|---|---|---|---|---|---|
| Module_267 | both | 26 | 12 | response to chitin | 1.2e-26 | 1.0e-21 |
| Module_514 | both | 16 | 10 | monoterpene biosynthetic process | 5.0e-27 | 8.6e-22 |
| Module_469 | both | 16 | 5 | spermine / spermidine biosynthesis | 1.0e-17 | 2.0e-13 |
| Module_1820 | both | 7 | 6 | proline biosynthetic process | 1.6e-17 | 2.8e-13 |
| Module_440 | both | 17 | 4 | ammonia assimilation cycle | 4.1e-11 | 3.5e-07 |
| Module_100 | both | 78 | 14 | **nitrate assimilation** | 9.6e-08 | 3.5e-04 |
| Module_026 | both | 324 | 34 | nitric oxide biosynthetic process | 1.1e-06 | 3.0e-03 |
| Module_009 | both | 521 | 96 | photosynthesis, light harvesting in PSI | 3.7e-10 | 2.3e-06 |

Nitrogen assimilation recurs across independent modules — nitrate assimilation
(Module_026, Module_100), the ammonia assimilation cycle and glutamate
biosynthesis (Module_440, Module_106), ammonium ion metabolism and the polyamines
(Module_469), proline and asparagine biosynthesis (Module_1820, Module_009,
Module_146, Module_191), urea transport (Module_117), nitric oxide, and cellular
response to nitrogen starvation (Module_417). These are different modules finding
different parts of the same pathway, which is the shape a real result has.

Module_009 (521 genes, the largest responsive module) is photosynthesis and
translation; Module_026 (324 genes, |r| = 0.96) carries both nitrate assimilation
and nitric oxide biosynthesis. That the two largest, strongest-responding modules
return interpretable primary metabolism is the check that the gene→GO join is
sound.

**By class:** 45 of 367 `both`, 22 of 239 `mi_only` and 4 of 41 `pearson_only`
modules were testable. `mi_only` is not enriched for anything the other classes
lack; with 22 testable modules it could not be shown either way. No class-level
claim is supported here.

**Purple is 3 modules.** 35 of its 38 responsive modules have fewer than 3
annotated members. The three that survive (Module_334, Module_4168, Module_097)
give phenylpropanoid biosynthesis, anion and phosphate transport, and acetyl-CoA
biosynthesis. That is a report of three gene sets, not a characterisation of
purple's nitrogen response, and it should not be compared with sugarcane's.

---

## Open items

- The 47,192 spurious sugarcane edges were present in every downstream result of
  the previous build. The affected artifacts are archived under
  `files/pearson_baseline/`; nothing there has been corrected, only superseded.
- ~~Gene–trait association is Pearson-only~~ — **done** (`12_gene_trait_mi.py`),
  and it did not fix the purple bottleneck. See above. What it did produce is a
  set of 2,297 sugarcane genes whose expression *dispersion* tracks nitrogen;
  whether that is biologically interesting is an open question, not a settled one.
- **Purple needs more libraries, not a better statistic.** Every underpowered
  result in this build traces back to n = 18: the soft MI floor for the network
  layer, the 62-gene trait selection, and the failure of the omnibus test to
  improve on it. That is worth stating in the methods rather than working around.
- **Purple's MI layer is matched on p, not on power.** At n = 18 the matched floor
  admits a 34x larger MI-only set than sugarcane's and those edges conserve
  slightly *worse* than linear ones. Worth testing whether an absolute-stringency
  purple MI layer recovers a conservation ratio above 1.
- Local transitivity was never computed (`COMPUTE_TRANSITIVITY=0`). It is only a
  reported column, but the TF and MYB61 tables carry NA for it.
