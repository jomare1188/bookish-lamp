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

Each eigengene is tested against nitrogen by **Spearman's rho alone**
(`19_module_trait_spearman.r`), at the same two thresholds the gene level uses:

> **padj ≤ 0.05 (BH over every module in the study) and |rho| ≥ 0.6.**

### Why one rank correlation, and not the three statistics this used to report

The module response was previously called by running `12_gene_trait_mi.py` on the
eigengene matrix, which returned Pearson *and* mutual information and sorted every
module into `pearson_only` / `mi_only` / `both`. That is gone. Three reasons, in
order of how much they mattered:

**1. The trait is ordinal, and Pearson reads it as interval.** Purple's nitrogen
is a dose — 0, 2 and 6 mM — and Pearson asks whether a module moves exactly twice
as far from 2 to 6 mM as it does from 0 to 2. Nothing in the design justifies
that arithmetic. A module whose response saturates above 2 mM is a perfectly good
monotone nitrogen response and Pearson penalises it. Measured on the identical
eigengenes at the identical thresholds:

| | Pearson | **Spearman** | in both | Spearman only | Pearson only |
|---|---|---|---|---|---|
| sugarcane (n = 48) | 408 | **465** | 394 | 71 | 14 |
| **purple (n = 18)** | **38** | **79** | 34 | **45** | 4 |

**Purple's responsive set more than doubles.** That is the single largest effect
of this change, and it is exactly where the ordinal argument predicts it: purple
is the study with the three-level gradient. Sugarcane's trait is two-level, where
Spearman is the rank-biserial correlation, so the gain there is smaller (and is
mostly robustness to the outliers a PC1 can carry).

**2. MI at this level was an omnibus test finding artifacts.** MI fires on *any*
dependence, including a dispersion change or a quirk in one or two libraries.
Under the old call, 35% of the 239 `mi_only` sugarcane modules had two of 48
samples carrying over a quarter of the eigengene's variance, against 14% for
`both`. Under the Spearman call the responsive set is *cleaner than the
background it is drawn from*:

| sugarcane | median top-2 sample share of eigengene variance | 2 samples carry > 25% | median R² of the nitrogen split |
|---|---|---|---|
| responsive (465) | **0.190** | **20%** | **0.447** |
| not responsive (6,111) | 0.249 | 50% | 0.040 |
| (flat null, 2 of 48 samples) | 0.042 | — | — |

The old `mi_only` class sat at 0.218 and 35%, i.e. *worse* than the non-responsive
background. Selecting on a monotone rank association removes that class of
artifact instead of concentrating it.

**3. The MI floor could only be calibrated, never derived.** `padj ≤ 0.05` alone
admitted modules down to |r| = 0.356, so an effect-size floor was needed; but
flooring only the linear side pushed every modestly-linear module into `mi_only`
(253 → 786), and the network's nats-to-|r| identity could not supply the MI
equivalent because it assumes two continuous variables and this trait is
discrete. The floor was therefore set to the median MI among modules at
|r| ≈ 0.6 — a calibration against the statistic it was supposed to be
independent of. One statistic needs one threshold and no calibration.

Everything downstream is simpler for it: no response classes to read the TF test,
the figures or the per-module GO by, and one rule for what "responsive" means,
owned by `19_module_trait_spearman.r` alone.

### What responds

| | sugarcane (n = 48) | purple (n = 18) |
|---|---|---|
| modules tested | 6,576 | 6,318 |
| padj ≤ 0.05 | 1,792 | **79** |
| \|rho\| ≥ 0.6 | **465** | 274 |
| **responsive (both)** | **465** | **79** |
| — rises with nitrogen | 242 | 32 |
| — falls with nitrogen | 223 | 47 |
| median \|rho\| (range) | 0.695 (0.602–0.866) | 0.787 (0.734–0.944) |
| median PC1 var | 80.6% | 79.9% |
| median module size | 5 genes | 5 genes |
| genes covered | 4,824 | 662 |

**The two studies are limited by different things, and it shows in which
threshold binds.** At n = 48 an |rho| of 0.6 already implies p ≈ 6e-06, so BH
never binds for sugarcane — the effect-size floor is the *only* active
constraint, and the corrected and uncorrected readings are the same 465 modules.
At n = 18 it is the reverse: 274 modules clear |rho| ≥ 0.6 and only 79 survive
BH. Purple's uncorrected count is therefore 274, and should be labelled that way
if it is ever used.

### The permutation null, and one caveat it raises

The same eigengenes against 1,000 shuffled trait label sets, counting how many
clear both thresholds:

| | observed | null median | null p99 | null max | permutations ≥ observed |
|---|---|---|---|---|---|
| sugarcane | 465 | 0 | 1 | 16 | **0 / 1,000** |
| purple | 79 | 0 | 10 | **244** | **2 / 1,000** (p = 0.002) |

Sugarcane is unambiguous. **Purple is significant but its null has a heavy tail**
— one shuffle in a thousand produced 244 "responsive" modules, three times the
observed count. That is what a single dense giant component does to a
permutation test: purple's eigengenes are strongly correlated with each other, so
one lucky label assignment lights up many modules at once and BH, which is
adaptive, then loosens for all of them. Purple's 79 modules are a real signal
(p = 0.002) but they are **not 79 independent findings**, and no per-module claim
from purple should be made without looking at the module itself.

### TF enrichment: real in sugarcane, and it is the falling modules

Hypergeometric per module against the network node universe, BH across modules,
then Fisher against the non-responsive modules (`15_module_profile.r`):

| sugarcane | modules | TF-enriched | OR | Fisher p |
|---|---|---|---|---|
| **responsive** | 465 | **3.23%** | **2.65** | **0.0016** |
| — falls with nitrogen | 223 | 4.04% | 3.34 | 0.0028 |
| — rises with nitrogen | 242 | 2.48% | 2.02 | 0.13 |
| not responsive | 6,111 | 1.24% | — | — |

This is a **restored** result, and the restoration is the point. Under the old
three-class call the same data gave `both` at OR 1.96, p = 0.061 and the
strongest class was `pearson_only` at 3 of 41 modules — a null spread across
three classes too small to separate. One responsive set of 465 recovers it:
**nitrogen-responsive modules are TF-enriched at OR 2.65**, and the effect sits
in the modules that go *down* with nitrogen. Do not over-read the direction split
— 9 enriched modules against 6 — but the pooled result is solid.

**Purple: 0 of 79, against 8 of 6,239 non-responsive (0.13%).** Nothing, as
before, and its TF-enrichment rate is an order of magnitude below sugarcane's
everywhere in the network.

### Figures

- `module_summary_<study>.{png,pdf}` — every responsive module's eigengene in one
  panel, **split by the sign of rho**: modules that rise with nitrogen above,
  modules that fall below. That replaces the old split by response class, and it
  is the more useful cut — the two halves are different biology rather than two
  ways of detecting the same thing. PC1 variance explained and module size ride
  as row annotations. Sugarcane draws the top 125 per direction of its 465;
  purple draws all 79. The vertical streaking that exposed the sample-driven
  `mi_only` block is gone from both, which is the visual form of the diagnostic
  table above.

  **Columns group replicates.** Sorting by sample name alone interleaves
  sugarcane's four leaf segments — the names run `B0_1, B_1, M_1, P_1, B0_2,
  ...`, so replicates of one segment sit four columns apart and the segment
  effect reads as vertical striping right across the figure. Ordering by
  treatment, genotype, **`segment`**, then sample (`HEATMAP_GROUP_BY`) puts
  replicates adjacent and the treatment blocks resolve. This grouped on the
  sheet's `tissue` column until 2026-08-22, which collapses `base0` and `base`
  into one label and so left those two interleaved inside a block of six — the
  striping the setting exists to remove. Leaf segment explains 0.802 of
  sugarcane's PC2 on the four real levels against 0.450 on the collapsed three. Sample labels are dropped (the design is carried
  by the annotation bars), which lets the cells shrink to 2 mm.
- `heatmaps/module_<id>_<study>.{png,pdf}` — per-module gene-level z-scores, top
  20 per direction. Height scales at 30.5 px per gene (r = 0.9998), so a 3-gene
  and a 100-gene module are both legible.

### What the responsive modules are for — GO, one test per module

One topGO enrichment per responsive module against the same network-node
background the gene-level GO and the TF hypergeometric use
(`./run.sh modulego <study> [BP|MF|CC]`). With one responsive set there are no
classes to pool; `direction` rides through as a column but does not partition
the run.

**The annotation gate is the binding constraint, and it is severe.**

| | sugarcane | purple |
|---|---|---|
| responsive modules | 465 | 79 |
| GO-annotated network nodes (background) | 8,251 of 103,336 | 12,255 of 170,736 |
| median GO-annotated members per responsive module | **0** | **0** |
| modules with *zero* annotated members | 293 (63%) | 51 (64%) |
| **modules testable (≥ 3 annotated members)** | **49 (11%)** | **10 (13%)** |

**89% of sugarcane's responsive modules cannot be tested at all.** Two facts
multiply: only 8% of sugarcane's network nodes carry any eggNOG GO annotation,
and the median responsive module holds 5 genes. So this section describes 49
modules, not 465, and the 49 are biased toward the large ones — nothing here
should be read as characterising the responsive set as a whole.

**Purple gains here.** It had 3 testable modules under the Pearson call and has
**10** under Spearman, because the responsive set it draws from more than
doubled. Ten gene sets is still a report of ten gene sets, not a
characterisation of purple's nitrogen response, but it is no longer a footnote.

| | BP | MF | CC |
|---|---|---|---|
| sugarcane: modules with ≥ 1 term | 44 of 49 | 40 of 49 | 34 of 49 |
| sugarcane: terms written | 246 | 137 | 66 |
| sugarcane: clearing cross-module BH | 73 | 36 | 14 |
| purple: modules with ≥ 1 term | 10 of 10 | 10 of 10 | 7 of 10 |
| purple: terms written | 50 | 30 | 13 |
| purple: clearing cross-module BH | 7 | 4 | 1 |

Within sugarcane's 11%, the signal is strong and it is about nitrogen:

| module | dir | genes | ann. | top BP term | p | global BH |
|---|---|---|---|---|---|---|
| Module_514 | ↑ | 16 | 10 | monoterpene biosynthetic process | 5.0e-27 | 5.9e-22 |
| Module_267 | ↑ | 26 | 12 | response to chitin | 1.2e-26 | 7.0e-22 |
| Module_469 | ↑ | 16 | 5 | spermine / spermidine biosynthesis | 1.0e-17 | 1.4e-13 |
| Module_1820 | ↑ | 7 | 6 | proline biosynthetic process | 1.6e-17 | 1.9e-13 |
| Module_613 | ↑ | 14 | 10 | protein folding in the ER | 4.8e-17 | 5.2e-13 |
| Module_469 | ↑ | 16 | 5 | ammonium ion metabolic process | 1.3e-14 | 1.2e-10 |
| Module_162 | ↓ | 45 | 7 | cellular response to cold | 4.3e-13 | 3.7e-09 |
| Module_215 | ↓ | 31 | 6 | flavonoid biosynthetic process | 7.4e-12 | 5.1e-08 |
| Module_440 | ↑ | 17 | 4 | ammonia assimilation cycle | 2.3e-11 | 2.3e-07 |
| Module_100 | ↑ | 78 | 14 | **nitrate assimilation** | 9.6e-08 | 2.8e-04 |
| Module_026 | ↑ | 324 | 34 | nitrate assimilation / nitric oxide | 1.1e-06 | 2.5e-03 |

Nitrogen assimilation recurs across independent modules — nitrate assimilation
(Module_026, Module_100), the ammonia assimilation cycle and glutamate
biosynthesis (Module_440), ammonium ion metabolism and the polyamines
(Module_469), proline and asparagine biosynthesis (Module_1820, Module_009,
Module_146, Module_191), and nitric oxide. These are different modules finding
different parts of the same pathway, which is the shape a real result has. **All
of them rise with nitrogen**, which is the direction they should.

**Every one of the previous build's headline modules survives the change of
statistic**, all in the positive direction: Module_026 (rho = 0.87, the strongest
responsive module in the study), Module_100 (0.81), Module_440 (0.78), Module_267
(0.75), Module_469 (0.73), Module_1820 (0.70), Module_009 (0.64), Module_514
(0.63). One drops out: **Module_117** (urea transmembrane transport) reaches only
|rho| = 0.43 and is now below the effect-size floor.

#### MF and CC as an internal check

The universe and the testable set are identical across the three ontologies (49
sugarcane modules, 10 purple), so they are comparable module by module — and
they are separate term sets scored in separate runs, which makes MF an
independent check on BP:

| module | BP said | MF said | global BH (MF) |
|---|---|---|---|
| Module_100 | nitrate assimilation | **nitrate reductase (NADH) activity** | 2.0e-04 |
| Module_026 | nitrate assimilation | **nitrate reductase (NADH) activity** | 4.0e-03 |
| Module_440 | ammonia assimilation cycle | **glutamate synthase activity** | 3.4e-07 |
| Module_1820 | proline biosynthetic process | **glutamate-5-semialdehyde dehydrogenase** | 4.2e-13 |
| Module_514 | monoterpene biosynthesis | **S-linalool synthase activity** | 5.3e-22 |
| Module_215 | flavonoid biosynthesis | **naringenin 3-dioxygenase activity** | 1.9e-10 |

Module_440's pairing is GS/GOGAT — the primary ammonium assimilation route — and
Module_100/Module_026's is the nitrate reductase step upstream of it. Different
ontologies, different tests, same modules.

CC behaves as CC usually does: its most recurrent term is *nucleus* in 8 modules,
which says little. Where it is sharp it corroborates — Module_009 is *chloroplast
thylakoid membrane* at p = 9.6e-47 against BP's photosynthesis and light-harvesting
call. MF's most recurrent terms are water channel activity and magnesium ion
binding (4 modules each), the former matching BP's *water transport*, which is
also BP's most recurrent term (4 modules).

**Figures.** `module_GO_<ONT>_<study>_global.{png,pdf}` ranks terms by how many
modules share them. Per-module panels live in `module_go/modules/<Module>/`, one
directory per tested module, with terms that survive the cross-module BH marked
in red so a panel cannot be over-read.

**Purple's 10 modules** give phenylpropanoid biosynthesis and protein
arginylation (Module_334, falls with nitrogen; MF: **phenylalanine ammonia-lyase
activity**, BH 2.8e-06), cellulose biosynthesis and primary cell wall biogenesis
(Module_455, rises; MF: **cellulose synthase activity**, BH 7.2e-06), sulfate
assimilation (Module_636; MF: phosphoadenylyl-sulfate reductase), and iron and
manganese homeostasis (Module_053, falls). Sulfate assimilation moving with
nitrogen is the one that connects to the nitrogen literature; the rest is a
report of gene sets. It should not be compared with sugarcane's — different n,
different annotation depth, and purple's permutation null has the heavy tail
described above.

---

## The reproduction figure — both source studies reproduce (with one asymmetry)

Before any cross-species claim, the re-quantification has to recover what each
study found on its own data. `./run.sh figrepro` draws that check
(`results/figures/figure2_reproduction.{png,pdf,svg}`); `./run.sh legends`
assembles its generated legend into `figures_legends.txt` at the repo root, and
the numbers are also written to `figure2_reproduction_stats.tsv`. The paper's Figure 1 is the
dataset, QC and quantification figure (`./run.sh figdataset`), whose panel D
carries a result worth reading before anything else: **PC1 is genotype in both
studies, at R² = 0.998 and 0.999**, and nitrogen loads on neither first
component. The nitrogen response has to be found against genotype and tissue,
not read off a leading component.

| | Muñoz-Perez 2025 — Module 20 | Ta Quang Kiet 2025 — MYB61 |
|---|---|---|
| what they claim | Module 20 is nitrogen-responsive | MYB61 responds non-monotonically, oppositely between genotypes |
| what we test | 33 mapped genes x 48 sugarcane libraries | 16 confirmed copies x 18 purple libraries |
| result | **24 of 32 testable respond in RB975375**, and **all 24 go down at high N** | **U-shape contrast padj < 0.05 in 10 of 13 copies in 51NG3; 0 in TAGZ** |
| variance | nitrogen 20.0% vs genotype 12.1% (median); N > genotype in 18 of 33 | — |
| **verdict** | **reproduces outright** | **reproduces in shape, inverted in sign** |

**Muñoz's claim reproduces cleanly.** Every significant gene moves the same way —
down at high nitrogen — across all four independent Module-20 loci, and nitrogen
beats genotype as a variance component in most of the set. This is the strongest
single validation the project has that the re-quantification is sound.

**Kiet's claim reproduces in shape but not in sign.** A significant non-monotonic
response exists and is restricted to one genotype, which is the structure they
describe. But all 10 significant copies **peak at 2N and fall at both 0N and 6N**,
where the paper reports maxima at the extremes. Two caveats belong with it: these
copies are barely expressed in leaf (most under 1 TPM), and the paper emphasises
roots, which this dataset does not contain.

The figure also carries the copies at the published id `Soff.09G0002230` as a
negative control. That id does not resolve to a MYB at all, and the genes at that
locus run one to two orders of magnitude hotter than any true MYB61 copy — the
likely source of the signal the paper attributes to MYB61. The full argument,
including the cloning-primer evidence that settles which gene they actually
worked on, is in `scripts/11_readouts/myb61/README.md`.

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
- **The gene level is still Pearson + MI; only the module level moved to
  Spearman.** The ordinal-trait argument that motivated the module change applies
  just as much to purple's 0/2/6 mM gradient at the gene level, and the module
  result suggests it would matter there too — Spearman found 45 purple modules
  Pearson missed. Re-running the gene selection on Spearman would, however,
  invalidate the conserved-response, edge-level and GO results that depend on it,
  so it is a deliberate open item rather than an oversight.
- **Purple's 79 responsive modules are not 79 independent findings.** Its
  permutation null reaches 244 in the worst of 1,000 shuffles because the network
  is one dense component. Any per-module claim from purple needs the module
  looked at, not just the padj.
