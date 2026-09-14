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

~~If purple's modules turn out to be uninformative, the fix is a stricter
threshold for purple specifically — not a different clustering algorithm.~~
**Superseded 2026-09-03.** The fix is neither: it is a k-NN degree reduction
applied to the matrix MCL clusters, leaving the network untouched. See
[the giant module is a node-degree artefact](#the-giant-module-is-a-node-degree-artefact-and-mcl-had-been-saying-so)
below. A stricter purple threshold remains a separate, unmade decision — it would
break the equal-specificity matching with sugarcane, so it is a deliberate trade
to make and document, not a tuning knob, and the evidence for it is now recorded
in `mcl_threshold_survey_purple.tsv`.

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


### The giant module is a node-degree artefact, and MCL had been saying so

The partition above was produced by `mcl <file>.abc --abc -I 2 -te 100 -o <out>`
and nothing else — no transform, and `stderr` discarded. Both of those turn out
to matter.

**MCL grades its own pruning, and the grade was going to `/dev/null`.** The
default `-scheme 7` tracks at most ~1,200 neighbours per node during computation;
purple's top decile of nodes has more than 33,000. Recovering the jury synopsis
from stderr, on the unreduced sugarcane graph:

| -I | 1.4 | **2 — the shipped partition** | 3 | 4 | 6 |
|---|---|---|---|---|---|
| jury | 34.8 woeful | **39.2 deplorable** | 43.4 poor | 45.5 dodgy | 47.2 shabby |
| largest module | 25.3% | 19.0% | 16.0% | 14.4% | 12.6% |

**No inflation value fixes either problem.** Across the whole range the largest
module never drops below 12.6% of the network, and mcl never rates its own
computation better than "shabby". This is what `mclfaq(7)` §7.3 predicts —
*"nodes of very high degree... tend to obscure cluster structure and contribute
to coarse clusters"* — and the degrees are extreme: purple's median is 859 with
p90 = 33,071 and a worst hub of 40,960, against the ≤ 100 median that
`clmprotocols(5)` recommends for co-expression graphs.

**The documented remedy works, on the same graph, at the same inflation.**
`-tf '#knn(k)'` requires an edge to be in the top k for *both* endpoints; k is
chosen by `mcx query -vary-knn` (`./run.sh mclsurvey`) as the gentlest reduction
meeting the degree target. Sugarcane, `#knn(180)`:

| | largest module | modules | median | singletons | Q | jury |
|---|---|---|---|---|---|---|
| `-I 2`, no reduction | **18.97%** | 10,710 | 3 | 401 | 0.0819 | 39.2 deplorable |
| `-I 2`, `#knn(180)` | **0.57%** | 32,240 | 1 | 16,144 | 0.2692 | **70.0 adequate** |
| `-I 1.4`, `#knn(180)` | 3.07% | 12,139 | 2 | 5,235 | **0.5506** | 54.6 tolerable |
| `-I 1.4`, `#knn(240)` | 3.53% | 9,704 | 3 | 3,048 | **0.5624** | 48.9 off colour |

Modularity rises **6.7-fold**. The community structure was always there; the hub
edges were masking it. `clm info`'s area fraction — the sum of squared cluster
sizes over N², i.e. the giant-module statistic — falls from 0.0402 to 0.0004.

**Purple is worse to start with and improves more.** Its shipped partition draws
the lowest jury grade in either study, and inflation does even less for it than
for sugarcane:

| purple | largest module | modules | median | singletons | Q | jury |
|---|---|---|---|---|---|---|
| `-I 2`, no reduction — **shipped** | **28.05%** | 24,475 | 1 | 14,594 | 0.1357 | **19.2 pathetic** |
| `-I 6`, no reduction | 20.15% | 101,101 | 1 | 96,605 | 0.1207 | 30.1 lousy |
| `-I 1.4`, `#knn(160)` | **4.24%** | 4,874 | **9** | **1,296** | **0.5249** | 54.2 tolerable |
| `-I 2`, `#knn(160)` | 0.36% | 100,291 | 1 | 82,960 | 0.2801 | 71.0 adequate |
| `-I 1.4`, `#knn(440)` | 9.31% | 2,522 | **13** | **372** | **0.5319** | 32.9 rotten |

`#knn(160) -I 1.4` improves purple on **every** axis at once, including the one
that is usually a trade: it strands **1,296** genes where the shipped partition
strands 14,594. Median module size goes from 1 to 9, and modularity nearly
quadruples. Note that k = 440 — the gentlest reduction meeting the median-degree
target — gives a *worse* jury grade than k = 160, because more surviving edges
means more pruning during the computation.

**The cost is real and is not hidden: singletons.** In sugarcane, 401 genes are
unassigned at the current setting; `#knn(180) -I 2` strands 16,144, and `-I 6`
strands 36,687. `#knn(180) -I 1.4` is the gentler trade at 5,235. Purple's curve
is kinder throughout, which is why k is chosen per study from its own survey
rather than shared.

**Are the smaller modules meaningful, or just rubble?** Homogeneity alone rises
trivially as modules shrink, so each clustering is scored against a size-matched
random partition (`37_cluster_homogeneity.r`; the metric is cogeqc's Sørensen–Dice
over eggNOG PFAM domains, reimplemented sparse because cogeqc's own version skips
groups above 200 genes and would refuse to score the giant module at all).
Excess over that null, sugarcane:

| clustering | modules scored | H | null | excess |
|---|---|---|---|---|
| shipped (`-I 2`, no knn) | 8,598 | 0.0835 | 0.0042 | +0.0793 |
| `#knn(180) -I 1.4` | 6,106 | 0.0457 | 0.0044 | +0.0413 |
| `#knn(180) -I 2` | 13,317 | 0.0910 | 0.0040 | **+0.0870** |
| `#knn(180) -I 3` | 15,970 | 0.1235 | 0.0034 | **+0.1200** |

The excess **rises** with inflation on the reduced graph, and the reduced
clustering beats the unreduced one at **every** inflation from 2 upward — 2:
+0.087 vs +0.079; 3: +0.120 vs +0.097; 4: +0.130 vs +0.106; 6: +0.139 vs +0.115.
The one cell that does worse is `-I 1.4`, the coarsest. So the extra granularity
is not rubble — but note that the reduced clusterings score
*fewer genes* overall, because a singleton cannot be scored, and that asymmetry
should be read alongside the excess rather than after it.

**Nothing is adopted.** `MCL_KNN_*` and `MCL_INFLATION_*` are empty, so
`./run.sh mcl` still reproduces the partition described above exactly. Changing
them regenerates `mcl_*_membership.tsv`, which every module-level stage reads.
Full grid: `results/<study>/mcl_sweep_<study>.tsv`; figure
`results/figures/figure10_clustering_choice.{png,pdf,svg}`.

**A correction to the sentence above.** This section previously said that if
purple's modules turned out to be uninformative "the fix is a stricter threshold
for purple specifically — not a different clustering algorithm". The first half
is half right and the second is not the issue: the fix is a degree reduction,
which is neither a threshold change nor a different algorithm. The
threshold evidence is recorded anyway — `mcl_threshold_survey_<study>.tsv`, from
`mcx query --vary-correlation` — because purple's `|r| >= 0.8` is p = 6.7e-5 at
n = 18 where sugarcane's is p = 9e-12, so the two networks were never built at
equal stringency. That remains a separate, unmade decision.


### Testing a stricter network: |r| >= 0.9

The threshold above is 0.8. A stricter cut was tested as a **separate network**
(`results_r09/`, `RESULTS=… STAT_MIN=0.9 ./run.sh build <study>`); nothing in this
document above or below describes it, and the 0.8 tree is untouched.

**What it does to the networks.** Both layers re-match to the new threshold —
`MATCH_PEARSON` follows `STAT_MIN`, which is what keeps the union licensed — so
the MI floor moves with the Pearson one:

| | sugarcane 0.8 | sugarcane 0.9 | purple 0.8 | purple 0.9 |
|---|---:|---:|---:|---:|
| pearson layer | 75,333,769 | **12,778,116** | 675,955,918 | **212,252,625** |
| ksg layer | 2,871,152 | **308,632** | 79,795,793 | **3,054,956** |
| ksg floor (nats) | 0.98911 | **1.13049** | 0.98822 | **1.29011** |
| **merged edges** | 76,200,344 | **12,971,011** | 705,571,723 | **213,892,607** |
| of which MI-only | 866,575 | **192,895** | 29,615,805 | **1,639,982** |
| nodes with an edge | 103,336 | **55,103** | 170,736 | **147,446** |

**The MI layer is what pays.** In purple it goes from 11.3% of the network
(both + MI-only) to 1.4%, an 18-fold drop in MI-only edges. That is
[the power argument in thresholds.md §3.3](thresholds.md) arriving in the data:
at n = 18 the KSG estimator is noisy, and demanding p <= 3.7e-07 rather than
6.7e-05 removes most of what it had to contribute. **If the MI layer is part of
what the network is for, a stricter threshold is expensive in purple.**

**A prediction that was wrong, and why it is worth recording.** The
`mcx query --vary-correlation` survey predicted purple would lose 0.3% of its
genes at |r| ~ 0.9. It lost **13.6%** (170,736 -> 147,446 nodes with an edge).
The survey cuts the *weight* scale, which derives from `stat`, and `stat` for an
MI edge is its r-equivalent — about 0.9 by construction at the 0.8-matched floor.
So every MI edge survives a weight cut and holds its nodes in; a real rebuild
re-matches the MI floor and they fall out. The survey remains the right tool for
choosing a **k-NN** cut, where nothing re-matches. It is the wrong tool for
predicting a threshold change.

**A common threshold does not make the two studies comparable.** It makes them
less so. The p implied by the cut moves 9.0e-12 -> 3.3e-18 in sugarcane and
6.7e-05 -> 3.7e-07 in purple, widening the specificity gap from 7.5e6-fold to
1.1e11-fold. Raising the threshold is not a route to equal stringency between
n = 48 and n = 18.

**And what it does to the modules — the question the test was for.** Purple:

| purple | \|r\| >= 0.8 | \|r\| >= 0.9 |
|---|---:|---:|
| nodes in the network | 170,736 | 147,446 |
| genes with no edge at all | 0 | **23,290** |
| named modules | 9,881 | 16,069 |
| largest module | 47,887 | **26,801** |
| — as % of its own network | 28.0% | **18.2%** |
| median module size | 3 | 3 |
| genes in no named module | 14,594 (8.5%) | **25,304 (14.8%)** |
| modularity Q | **0.1631** | **0.1299** |

**It helps, and it is not enough.** The giant module shrinks from 28.0% to 18.2%
of the network — real, but still a giant — while modularity goes *down*, the
median module stays at 3, and the share of genes in no named module nearly
doubles. The stricter threshold buys a smaller giant by discarding the genes that
were holding it together, rather than by revealing structure.

**Against the k-NN reduction on the unchanged 0.8 network, it is not close:**

| purple, largest module | genes in no named module | Q |
|---|---|---|
| 0.8 baseline — 28.0% | 8.5% | 0.163 |
| **0.9 rebuild — 18.2%** | **14.8%** | **0.130** |
| **0.8 + `#knn(160)` `-I 1.4` — 4.24%** | **0.8%** | **0.525** |

k-NN gives a giant module four times smaller, a modularity four times higher, and
strands an order of magnitude fewer genes, on a network that was never rebuilt.
Sugarcane tells the same story more bluntly: 46.7% of its nodes gone, modularity
halved, giant module still 14.1%.
**The answer to "was 0.8 too loose" is: not in the way that matters for
clustering.** Too many edges per *node* is the problem, and a global threshold is
a blunt instrument against it — it removes weak edges everywhere, including the
ones that were a sparse gene's only connections, while the hubs keep enough of
their thousands to go on dominating.

**Sugarcane is the harder case, and it fails the test outright.**

| sugarcane | \|r\| >= 0.8 | \|r\| >= 0.9 |
|---|---:|---:|
| nodes in the network | 103,336 | **55,103** |
| median node degree | 52 | 8 |
| p90 node degree | 6,666 | 1,672 |
| named modules | 10,309 | 9,334 |
| largest module | 19,604 | **7,782** |
| — as % of its own network | 19.0% | **14.1%** |
| median module size | 3 | 3 |
| modularity Q | **0.1005** | **0.0524** |

**It loses 46.7% of the network's nodes and halves modularity to buy a giant
module that is still 14.1%.** Of the 170,790 genes in the VST, 60.5% were in the
0.8 network and only 32.3% survive in the 0.9 one. The median module is still 3.
Nothing about the clustering is better; a great deal of the data is gone.

**The prediction check passed exactly, which is worth recording separately from
the verdict.** Filtering the 0.8 table forecast **12,970,859** edges; the rebuild
produced **12,971,011** — 152 edges apart, 0.001%. The forecast even got the
layer reassignment right: 115,754 edges predicted to stay `both`, 115,737 actual;
192,911 predicted to become MI-only, 192,895 actual. Sugarcane's permutation null
is capped at `MAX_PERM` under both thresholds, so the filter and the rebuild are
the same computation, and they agree to four decimal places.

`min_module_size` was briefly run at 5 during this build. It does **not** affect
the clustering — the raw partition (10,710 modules) and Q (0.1005) came out
identical either way; it only decides how small a module may be and still get a
`Module_NNN` name rather than "Unassigned" (3,307 named / 18,935 unassigned at 5,
versus 10,309 / 401 at 2). Reverted to 2 to match the original pipeline.

---

### Choosing k-NN by a graph model: the criterion runs, and it cannot choose

Branch `pearson-knn-ba`. Separate results tree (`results_pearson/`, gitignored),
separate MCL work dir. Nothing in `results/` was touched. **This analysis is
recorded and closed** — the bulk data it produced has been deleted; the summary
tables it rests on are in `docs/data/knn_ba/`. The work going forward uses
**unpruned** graphs.

The question: `#knn(k)` demonstrably breaks the giant module, but k was picked by
hand, and a hardcoded k is not defensible. The proposal was to choose k per
network as the value whose reduced graph is closest to a Barabási–Albert model
under `statGraph::graph.model.selection` — a scale-free network being the shape a
co-expression network is supposed to have. **The answer is that this criterion
cannot pick k, for a reason that is a property of the statistic and not of our
networks.**

#### The network the sweep ran on: Pearson only

MI edges are a small minority, so the hairball is essentially all linear and
dropping them costs little:

| | sugarcane | purple |
|---|---|---|
| nodes | 101,990 | 170,135 |
| edges | 75,333,769 | 675,955,918 |
| genes lost vs the merged network | 1,346 | 601 |
| MI-only share of the merged network | 1.14% | 4.20% |

Built by one streaming `awk` over the existing Pearson layer
(`scripts/40_pearson_only_network.sh`), reproducing the merge's weight formula
exactly; no re-run of the correlation engine.

#### `method="fast"` is the only affordable path, and it was validated first

`graph.model.selection` defaults to `method="diag"`, a full dense
eigendecomposition per grid point. Measured scaling on this machine is cubic
(exponent 3.10): at n = 170,736 one `eigen()` is **34.6 h** on a **233 GB**
matrix, and a default selection wants 25,100 of them — about **99 years**.
`method="fast"` derives the spectral density analytically from the degree
distribution and touches no dense matrix: **2.32 s** at n = 1,000 rising only to
**5.11 s** at n = 170,736, even as unique degrees go 117 → 694.

The validation gate (`scripts/39_statgraph_validate.r`,
`docs/data/knn_ba/statgraph_validation.tsv`) returned **PARTIAL**:

- **BA recovery 2/2.** On synthetic BA graphs at n = 1,000 and n = 10,000, fast
  named BA, with GIC(BA) ~ 0.0025 against ~0.34 for ER and WS. The one model we
  actually score on is the one it recovers.
- **All-model recovery 3/6.** It confuses ER with WS in both directions. That is
  a real limitation, but ER-vs-WS is not the discrimination this analysis needs.
- **fast vs exact on real subgraphs.** They agree on a sparse snowball-sampled
  subgraph (mean degree 4.4) and disagree on a dense one (mean degree 30.7: fast
  says BA, exact says WS). This is the finding that explains everything below.

Two implementation traps found and worked around, both worth keeping:
`numCores > 1` **deadlocks** (a PSOCK cluster is created per spectral density and
leaks connections; the main process sits at 0% CPU) — so statGraph runs
single-threaded and parallelism goes across k instead. And statGraph's default ER
grid `seq(0, 1, 0.01)` is unusable: p = 0 crashes, and p = 1e-06 returns a
**negative** GIC that beats every real model, which labelled a synthetic BA graph
as ER. All grids are given explicitly with the degenerate ends excluded.

#### The grids

k = 50…600 in steps of 50, plus `none` as control (`scripts/41_knn_ba_sweep.sh`,
exact counts from `mcx query` on each reduced matrix). Selected rows:

**Sugarcane** (unreduced: 101,990 nodes / 75,333,769 edges)

| k | nodes | edges | median deg | singletons |
|---|---|---|---|---|
| 50 | 85,050 | 448,140 | 6 | 16,940 |
| 200 | 96,335 | 1,724,999 | 13 | 5,655 |
| 400 | 100,335 | 2,991,867 | 20 | 1,655 |
| 600 | 100,771 | 3,896,333 | 24 | 1,219 |

**Purple** (unreduced: 170,135 nodes / 675,955,918 edges)

| k | nodes | edges | median deg | singletons |
|---|---|---|---|---|
| 50 | 156,796 | 1,360,219 | 14 | 13,339 |
| 200 | 157,590 | 4,614,459 | 47 | 12,545 |
| 400 | 159,222 | 7,922,149 | 78 | 10,913 |
| 600 | 162,776 | 11,358,388 | 106 | 7,359 |

Full grids with GIC and power-law columns:
`docs/data/knn_ba/knnba_selection_{sugarcane,purple}.tsv`.

#### Why the criterion cannot choose k

**BA was selected at every single k in both species — 16/16 sugarcane, 12/12
purple.** Against a criterion of "closest to BA" that sounds like success. It is
not, because GIC(BA) is **monotone in sparsity**: its minimum always sits on the
sparsest edge of whatever grid you offer it.

Sugarcane GIC(BA) climbs 0.394 (k=50) → 0.677 (k=100) → 0.990 (k=200) → 2.181
(k=600). Purple climbs 0.0714 (k=50) → 0.153 (k=200) → 0.207 (k=400) → 0.232
(k=600). The argmin is k=50 in both — the boundary of the grid. Extending
sugarcane *below* the requested range settles it: k=40 gives 0.321, k=30 gives
0.239, k=20 gives 0.154, k=10 gives **0.0714**. It keeps falling. There is no
interior optimum to find; the criterion answers "as sparse as you will let me",
whatever grid it is handed.

The mechanism is check B above. `method="fast"` assumes a locally tree-like
graph. As k rises the reduced graph gets denser and more clustered, the
approximation degrades, and GIC rises — so GIC is tracking **how well the
approximation holds**, not how BA-like the graph is. Its ordering across k is not
a comparison of graphs.

**The power-law fallback did not rescue it.** Run at every k regardless
(`igraph::fit_power_law`, `implementation="plfit"`), the best KS statistic lands
at k=300 in sugarcane and k=50 in purple — inconsistent across two species built
by the same pipeline — and α wanders implausibly (purple k=150 returns
α = 31.7, a fit failure, not an exponent).

#### What k-NN does achieve, and the mistake in reading it

The reduction does exactly what it was brought in to do. Giant module and MCL's
own jury grade of the resource scheme:

| | sugarcane | purple |
|---|---|---|
| unreduced | 19,840 (19.45%), jury 39.2 *deplorable* | 43,816 (25.75%), jury 20.1 *awful* |
| k = 50 | 920 (0.90%), jury 95.5 *fabulous* | 1,110 (0.65%), jury 91.9 *cracking* |
| k = 400 | 6,046 (5.93%), jury 78.9 *groovy* | 7,020 (4.13%), jury 63.1 *fairish* |

**I initially reported the modularity rise (0.081 → 0.645 in sugarcane) as
evidence that reduction improves structure. That was wrong, and the user's
question — "how do you measure the structure to say that it peaks at k=400" —
is what exposed it.** Each k's Q had been computed on its *own* reduced graph, so
those numbers are not comparable to one another: a sparser graph makes any
partition look more modular. The apparent purple "peak" at k=200 also beat k=150
by 0.00004, which is noise, not a peak.

Scored properly — every partition evaluated against **one fixed reference**, the
unreduced Pearson network (`scripts/45_score_vs_reference.sh`, `clm info`) — the
finding **inverts**:

| | sugarcane Q_own | sugarcane Q_ref | sugarcane mass frac | purple Q_ref | purple mass frac |
|---|---|---|---|---|---|
| k = 50 | 0.645 | 0.014 | 0.340 | 0.012 | 0.187 |
| k = 400 | 0.732 | 0.057 | 0.586 | 0.068 | 0.311 |
| k = 600 | 0.720 | 0.059 | 0.606 | 0.086 | 0.362 |
| **unreduced** | **0.081** | **0.081** | **0.687** | **0.150** | **0.634** |

> **Correction (2026-09-09).** The mass-fraction column above is unreliable and
> should not be quoted. `clm info` was called with all ~25 partitions in one
> invocation, and its `eff` and `mf` depend on which other clusterings share the
> call — measured later on the same tool: one cluster file scored alone gives
> `eff=0.47281 mf=0.51576`, and in a batch of eight `eff=0.37821 mf=0.46878`.
> `mod` and `af` are unaffected, so **the modularity column, which is what the
> conclusion rests on, stands**; `Q_own` and `Q_ref` are both modularity. The
> partitions were deleted with the k-NN bulk data, so the mass fractions cannot
> be recomputed. `45_score_vs_reference.sh` now scores one partition per call.

Against the real network the unreduced partition has the highest modularity
**and** — on the modularity evidence — the partition the reduced ones fail to
improve on, in both species. k-NN buys a smaller giant
module by discarding the edges that made it, and what it discards is signal the
partition then fails to explain. Q_ref rises monotonically with k precisely
because larger k throws less away.

The one metric that still favours small k is PFAM annotation homogeneity above a
size-matched null (sugarcane k50.I2 +0.1049 vs unreduced I2 +0.0779; purple
k50.I2 +0.0254). That is expected and not decisive: a partition of many tiny
modules scores well on homogeneity almost mechanically, and the null correction
does not fully absorb it.

#### Verdict

**k-NN reduction is not adopted.** The objective criterion asked for does not
exist in usable form — GIC(BA) under `method="fast"` measures its own
approximation error, not BA-likeness, and the power-law fallback is inconsistent
across species. And when the reduced clusterings are scored like-for-like against
the real network, reduction is a net loss. The giant module is real structure in
a dense graph, not an artefact k-NN should be tuned to erase. Subsequent work
uses the unpruned graphs.

Machinery kept and reusable: `scripts/39_statgraph_validate.r`,
`40_pearson_only_network.sh`, `41_knn_ba_sweep.sh`, `42_knn_model_selection.r`,
`43_knn_ba_collect.r`, `44_run_knn_ba.sh`, `45_score_vs_reference.sh`
(`run.sh` stages `sgvalidate pearsononly knnsweep knnselect knncollect scoreref`).
`45_score_vs_reference.sh` is the generally useful one: it scores any set of
clusterings against one fixed graph, which is the only way to compare partitions
of differently-pruned networks.

---

## Choosing the clustering: inflation, and whether MCL is the right tool

Branch `clustering-methods`, tree `results_cluster/`. The graph is the
**unpruned Pearson-only network** — `|r| >= 0.8`, weights rescaled to
`[0.01, 1]`, no k-NN reduction (sugarcane 101,990 nodes / 75,333,769 edges,
mean degree 1,477; purple 170,135 / 675,955,918, mean degree 7,946). Built by
streaming the Pearson layer straight into `mcxload` (`46_pearson_mci.sh`), which
skips the 70 GB nine-column table whose only consumer was the loader; edge counts
match the layer summaries exactly and no node has degree 0.

Two questions, neither previously answered. `-I 2` was inherited, never chosen.
And MCL was never compared against anything.

### Everything is scored on one graph, by one tool, and that had to be verified

Every partition — MCL's, Leiden's, Louvain's — is scored by `clm info` against
the same unpruned matrix, so only the partition varies. Leiden's memberships are
written in mcl's own cluster format to make this possible.

That introduces one place where a silent error would invalidate everything, so it
has its own gate (`verify_cluster_roundtrip.py`). igraph's `Read_Ncol` assigns
vertex ids by **order of first appearance**, not by the integers in the file —
measured: a file beginning `5 3` makes node 5 into vertex 0. Writing such a
membership against the `.mci` node domain would put every gene in the wrong
cluster, and `clm info` would score it happily. The gate pushes a real MCL
partition through the Leiden writing path and requires identical statistics, and
carries a **negative control** that drops the permutation deliberately and
requires the numbers to move. Both species pass on all seven statistics
(sugarcane `cls.knone.I12`, purple `cls.knone.I2`).

**A defect in `clm info` itself was found on the way, and it matters beyond this
section.** `clm info <graph> <cls1> <cls2> ...` is documented to take many
clusterings at once, but `eff` and `mf` depend on which others share the call:

| | eff | mf |
|---|---|---|
| `cls.knone.I6` scored alone | 0.47281 | 0.51576 |
| the same file in a batch of 8 | 0.37821 | 0.46878 |

`mod` and `af` are unaffected. The batched columns are not merely offset, they
are **not monotone**: batched, `mf` rises from 0.469 at `-I 6` to 0.479 at
`-I 8` and then falls to 0.410 at `-I 10`, while cluster count, largest module,
`af` and `mod` move smoothly through the same cells — a discontinuity that reads
as a feature of the data and is an artefact of the call. All three scorers now
invoke `clm info` once per partition. **This also affects the k-NN section
above**, whose mass-fraction column came from a single call holding ~25
partitions; that section carries a correction, and its conclusion rests on
modularity, which is unaffected.

### The MCL inflation ladder

**Sugarcane** (`-I 20` is excluded from every conclusion: it underflowed,
zeroing 35,044 of 101,990 vectors, and returned a plausible-looking 34.8% giant):

| -I | clusters | largest | singletons | eff | mf | af | mod | jury |
|---|---|---|---|---|---|---|---|---|
| 1.2 | 2,350 | 31.43% | 0 | 0.1146 | 0.8851 | 0.1233 | 0.0802 | 33.0 rotten |
| 1.7 | 7,064 | 21.02% | 37 | 0.2768 | 0.7629 | 0.0510 | **0.0818** | 37.1 deplorable |
| 2 | 8,795 | 19.45% | 105 | 0.3258 | 0.7200 | 0.0428 | 0.0806 | 39.2 deplorable |
| 6 | 25,538 | 15.56% | 7,315 | 0.4728 | 0.5158 | 0.0247 | 0.0600 | 47.0 shabby |
| 13 | 35,573 | 14.27% | 15,344 | **0.4960** | 0.4455 | 0.0206 | 0.0548 | 47.9 shabby |

**Purple**:

| -I | clusters | largest | singletons | eff | mf | af | mod | jury |
|---|---|---|---|---|---|---|---|---|
| 1.2 | 716 | 38.34% | 0 | 0.1484 | 0.7879 | 0.1773 | 0.1377 | 16.1 terrible |
| 2 | 7,722 | 25.75% | 276 | 0.3859 | 0.6527 | 0.0817 | 0.1503 | 20.1 awful |
| 3 | 22,545 | 23.60% | 9,383 | 0.4699 | 0.5805 | 0.0665 | **0.1529** | 24.1 miserable |
| 4 | 35,561 | 22.67% | 21,218 | 0.4884 | 0.5428 | 0.0610 | **0.1529** | 26.8 abominable |
| 6 | 50,417 | 22.37% | 36,132 | **0.4902** | 0.4968 | 0.0580 | 0.1457 | 29.8 abominable |

**Inflation cannot break the giant module in either species.** Sugarcane goes
31.4% -> 14.3% and purple 38.3% -> 22.4%, and the cost is 15,344 and 36,132
singletons respectively. Purple's giant is the harder one and never drops below
22%.

**Most of the criteria cannot choose a setting, because they are monotone.** `mf`
falls monotonically (it picks the coarsest cell; the one-cluster baseline scores
`mf = 1.0`, which is why it is in the table). `af` falls monotonically (it picks
the finest). `eff` rises monotonically across MCL's entire usable range in both
species and is still climbing at the last valid cell — the argmax sits on the
grid boundary, exactly the failure that sank the Barabási–Albert criterion in the
section above. **Only modularity has an interior optimum**, and the ladder was resampled densely
around it to make sure the peak is real rather than an artefact of a coarse grid:
sugarcane peaks at **`-I 1.5`** (0.08212, on a broad plateau from 1.4 to 1.7,
where the whole range varies by under 0.5%), purple at **`-I 3.5`** (0.15324, on a
plateau from 2.7 to 4.5). So `-I 2` is **not** the optimum in either species — it
is the value that sits between the two optima, and it keeps **98.1% of the peak in
both**. That, rather than a peak at 2 that does not exist, is the argument for
using one inflation across two networks. **Figure 12**, whose legend
(`figure12_cluster_methods_legend.txt`) is generated from the same table.

The ladder cannot simply be extended to find `eff`'s peak, because mcl breaks
first, in two ways that are silent unless stderr is kept:

* `-I` above 30 is out of range. mcl warns to stderr and **uses the default 2.0**
  instead. Measured: `-I 40` returned `-I 2`'s partition byte-for-byte.
* around `-I 15` and above it stops converging. `-I 15` ran 2,027 iterations with
  mcl auto-escalating the inflation to 109 before it was killed; `-I 20`
  completed but underflowed 34% of the graph.

### Modularity maximisation is the wrong objective on these graphs

Leiden with the modularity objective, and Louvain, both score far higher
modularity than anything else tested — and reach it by building giant modules:

| | clusters | largest | eff | mod | PFAM excess |
|---|---|---|---|---|---|
| sugarcane leiden-mod | 1,107 | 42.91% | 0.0719 | **0.1399** | +0.0471 |
| sugarcane louvain | 1,084 | 40.79% | 0.0698 | 0.1386 | +0.0448 |
| sugarcane best MCL | 7,064 | 21.02% | 0.2768 | 0.0818 | +0.0695 |
| purple leiden-mod | **55** | 27.73% | 0.1384 | **0.2122** | **+0.0004** |
| purple louvain | **50** | 27.36% | 0.1374 | 0.2121 | **+0.0002** |
| purple best MCL | 22,545 | 23.60% | 0.4699 | 0.1529 | +0.0211 |

Purple's modularity optimum is **fifty clusters for 170,135 genes**. This is the
resolution limit doing exactly what theory says: modularity cannot resolve
communities below ~sqrt(2m) edges, and 2m here is 1.35e9. The independent check
settles it — those partitions carry **no PFAM signal above a size-matched random
partition of the same shape** (+0.0002, +0.0004). They are the highest-modularity
and the least biologically meaningful partitions in the table, simultaneously.

**So the giant module was never an MCL artefact. It is what modularity
maximisation actively wants.** Louvain and Leiden-modularity should not be used
on these networks.

### MCL and Leiden CPM trace the same curve; only the reach differs

CPM has no resolution limit, and Leiden-CPM does dissolve the giant module —
sugarcane 17.2% at gamma 0.005 down to 0.31% at gamma 0.9, purple 27.4% down to
2.3%. It also attains the highest `eff` of any method, and — unlike MCL — its
`eff` has a genuine **interior** optimum: gamma = 0.1 in sugarcane (0.5504),
gamma = 0.2 in purple (0.5611).

But that is not because it is the better algorithm at a given granularity. Sorted
by area fraction the two methods **interleave on one curve**, and in sugarcane
MCL is consistently the higher of the two where they overlap:

| af | method | eff | singletons |
|---|---|---|---|
| 0.0205 | Leiden CPM gamma 0.025 | 0.4760 | 14,556 |
| 0.0206 | MCL -I 13 | **0.4960** | 15,344 |
| 0.0271 | Leiden CPM gamma 0.01 | 0.4015 | 11,495 |
| 0.0282 | MCL -I 4 | **0.4488** | **3,239** |
| 0.0331 | Leiden CPM gamma 0.005 | 0.3435 | 10,621 |
| 0.0358 | MCL -I 2.5 | **0.3808** | **507** |

The same holds for the biological check. At matched module count MCL wins in
**both** species — sugarcane ~19.5k modules +0.1051 (`-I 4`) against +0.0680
(gamma 0.01), ~31.5k +0.1187 (`-I 9`) against +0.1009 (gamma 0.05); purple ~22k
+0.0211 (`-I 3`) against +0.0166 (gamma 0.1), ~36k +0.0251 (`-I 4`) against
+0.0223 (gamma 0.2). Leiden CPM's higher raw homogeneity comes from fragmenting
past where MCL can go, and its apparent peak there is hollow: at sugarcane gamma
0.7 only 3,982 of 90,070 modules have enough annotation to score at all, falling
to 1,174 at gamma 0.9, which is why the number then drops.

**Leiden CPM's advantage is reach, not quality.** MCL's area fraction bottoms out
at 0.021 (sugarcane) and 0.058 (purple) because inflation underflows before it
can fragment further; the peak of the shared curve lies below that, where only
CPM can operate.

On modularity specifically the ranking is **species-dependent inside the overlap**
and is not claimed either way: at ~19,300 clusters sugarcane's Leiden CPM is above
MCL (0.0799 against 0.0664), while purple's MCL is above Leiden CPM throughout.
What holds in both is that **MCL's maximum exceeds Leiden CPM's maximum** (0.0821
against 0.0803; 0.1532 against 0.1426), reached at a coarser granularity than the
overlap. The per-method numbers are in
`docs/data/clustering/cluster_methods_{sugarcane,purple}.tsv`; the only axis on
which the methods can be plotted together is cluster count, since inflation and
gamma are not commensurable.

### MCL's internal pruning is not distorting the answer

mcl prunes each node's neighbour list *while computing* — `-scheme 7`, the
default and the highest preset, keeps 1,200 neighbours against mean degrees of
1,477 and 7,946 — and grades its own pruning "deplorable" to "abominable"
throughout. That confound had to be measured rather than assumed:

| sugarcane `-I 1.7` | clusters | largest | eff | mod | jury |
|---|---|---|---|---|---|
| `-S 1200` (default) | 7,064 | 21.02% | 0.27685 | 0.08184 | 37.1 deplorable |
| `-S 4000` | 7,106 | 20.39% | 0.27964 | 0.08239 | 56.2 acceptable |
| `-S 10000` | 7,106 | 20.39% | 0.27964 | 0.08239 | 56.2 acceptable |

`-S 4000` and `-S 10000` are **identical**: the pruning saturates, and the
converged answer is within 0.7% of the default.

Purple is the case that matters, since its grades are the worst:

| purple `-I 3` | clusters | largest | eff | mod | jury |
|---|---|---|---|---|---|
| `-S 1200` (default) | 22,545 | 23.60% | 0.46992 | 0.15294 | 24.1 miserable |
| `-S 10000` | 22,508 | 23.52% | 0.47108 | 0.15324 | 30.8 lousy |

Eight times the headroom moves **every statistic by under 0.3%**. The grade rises
but stays poor, which it must: at mean degree 7,946 even `S = 10000` is only
1.26x a typical neighbour list, whereas sugarcane's 1,477 is comfortably covered.
The partition has nonetheless converged.

So **the jury grade is not a guide to whether a result is affected** — it grades
how much was discarded, not whether discarding it mattered. The confound flagged
at the start of this section is real in the grade and absent in the answer, in
both species.

### The hierarchical family

Classical hierarchical clustering was not attempted, and the reason is not
squeamishness about cost: it needs a dense `170,135 x 170,135` dissimilarity
matrix — 232 GB, ~1.4e10 pairs, with O(n^2 log n) linkage on top — built from
precisely the correlations this pipeline thresholded away at `|r| >= 0.8`. It
would be answering a different question from every other row in the comparison.

The graph-native member of the family that still yields a dendrogram is
fast-greedy (Clauset–Newman–Moore), run on sugarcane only, under a 4 h cap.
**It is not adopted, for a reason that has nothing to do with runtime.** The
sugarcane graph has **958 connected components**, so igraph builds only 101,032
merges for 101,990 vertices and no cut below 958 clusters exists; CNM's own
`optimal_count` is **338**. Its preferred partition is unreachable on this graph.
And CNM optimises modularity, which the section above shows is the wrong
objective here regardless.

On cost the evidence is mixed and is reported as such: one run built the
dendrogram in under 2 h (and then raised on the cut, which is how the component
floor was found), while a later run on the same graph and code exceeded the 4 h
cap with the worker's CPU share falling from 34% to 24% under contention. So it
is affordable on an idle machine and marginal on a busy one — but the
component floor and the objective settle the question either way, and purple, at
nine times the edges, was not attempted.

### What to use

* **`-I 2` is defensible.** The modularity optima are `-I 1.5` (sugarcane) and
  `-I 3.5` (purple), and `-I 2` retains **98.1% of the peak in both** — it is the
  one setting that is near-optimal for two networks at once. Modularity is also
  the only criterion with an interior optimum for MCL, so it is the only one that
  can make this argument at all. The inherited setting survives the test.
* **Do not use Louvain or Leiden-modularity here.** They win the metric they
  optimise and lose everything else, catastrophically in purple.
* **Leiden CPM at gamma 0.1–0.2 is the option if the giant module must go.** It
  is the only method that dissolves it, and it holds the highest `eff` on either
  graph. The price is 26,539 singletons in sugarcane (26% of genes) and 19,476 in
  purple, and slightly worse annotation coherence per module than MCL at matched
  granularity.
* **Raising inflation is not a way to break the giant module.** It costs genes far
  faster than it costs the giant, and above `-I 13` mcl stops being numerically
  trustworthy.

Tables: `docs/data/clustering/`. Figure: `figure12_cluster_methods` (one panel:
modularity vs inflation, peak marked), with its legend in
`figure12_cluster_methods_legend.txt`.
Machinery: `scripts/46`–`50`, `verify_cluster_roundtrip.py`, `run.sh` stages
`pearsonmci mclladder leidensweep clustercompare`.

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

### And it survives a non-monotone test, which is the one blind spot that mattered

Every selection rule above is **monotone** — Pearson, mutual information, their
union. Purple's design is stress–control–stress (0 mM and 6 mM are deficiency and
excess, 2 mM the control), so a gene moved the *same way* by both stresses is
structurally invisible to all of them. Until now the U-shape contrast
`c(+1,−2,+1)` had only ever been run on a handful of named genes and once at
module level, so the negative above was, strictly, a claim about monotone
responses only.

It has now been run genome-wide at gene level (`./run.sh ushape`, then
`./run.sh conscor 0 ushape` and `conscor 1 ushape`). Blocked model
`expr ~ genotype + treatment`, genotype a block rather than a pool because 51NG3
is *S. robustum* and TAGZ *S. officinarum*; 170,740 genes, residual df 14.

| | |
|---|---|
| genes tested | 170,740 |
| **survive BH at 0.05** | **1** (a peak at the control, not a trough) |
| raw p ≤ 0.05, uncorrected | 17,528 |
| expected by chance | 8,537 |
| **fold over chance** | **2.05×** |
| **conserved responsive ortholog pairs, genome-wide burden** | **0** |
| **conserved responsive ortholog pairs, directed burden (BH over 4,745)** | **0** |

**Both readings matter, and they say different things.** The 2.05× enrichment at
raw p — roughly 9,000 genes' worth of excess — says purple really does contain
non-monotone nitrogen responses that the monotone tests cannot see. The single BH
survivor says that at n = 18 over 170,740 tests, essentially none of them can be
individually resolved. The blind spot was real; closing it changes nothing,
because the wall is power, not test shape.

Under the directed burden the U-contrast selects **0 of 4,745** candidate
orthologs, against the monotone rules' 1–2. So the non-monotone test does not
merely fail to add conserved pairs, it adds fewer candidates than the monotone
one at the same burden.

**`Soffic.09G0001580-9H`, calibrated at last.** (Full per-gene dossier:
[gene_Soffic_09G0001580_9H.md](gene_Soffic_09G0001580_9H.md).) The Module-20 copy that motivated
this — the one concrete case of a gene the monotone tests miss — sits at **rank
890 of 170,740** (percentile **0.52**), u_est **+3.67** (a trough: up at both 0 N
and 6 N), raw p **0.0037**, **padj 0.485**. It is genuinely in the top half-percent
genome-wide and it is nowhere near surviving correction. It remains an anecdote,
but now a numbered one.

The p-value histogram is skewed toward zero (19% in the first decile against 10%
expected), consistent with that diffuse real signal rather than with a flat null.
That is a reason to believe the 2.05×, and equally a reason not to read anything
into the top of the list gene by gene.

**What this bounds.** The sentence "the two species' responsive sets are
independent in ortholog space" no longer carries an unstated *"among monotone
responders"*. It survives a non-monotone test under both correction burdens. What
it does not survive is a power argument — and that was already the finding.

One asymmetry is forced and should be stated wherever this is quoted: sugarcane's
design has **two** nitrogen levels and cannot express a U-shape, so under
`SELECTION=ushape` purple is selected by the contrast while sugarcane keeps its
monotone rule. The question answered is *"does purple have non-monotone responders
whose orthologs are sugarcane-responsive?"* — not *"do the two species share a
non-monotone response?"*, which these two designs cannot ask.

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

> **Superseded on 2026-09-10 by "The module analysis, rebuilt" below.** Everything
> in this section was computed on the **merged (Pearson + MI) network at `-I 2`**
> with a **marginal** Spearman. The analysis now runs on the unpruned Pearson-only
> graphs at each species' modularity optimum, with the experimental design in the
> model. The numbers here are kept because the reasoning that produced them —
> especially why one rank correlation and not three statistics — still stands and
> is not repeated below. **Do not quote these counts against the new modules**, and
> note that module names are positional, so `Module_001` is not the same gene set
> in the two analyses. See `results/STALE.md`.


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
levels are 0, 2 and 6 mM, and Pearson asks whether a module moves exactly twice
as far from 2 to 6 mM as it does from 0 to 2. Nothing in the design justifies
that arithmetic. A module whose response saturates above 2 mM is a perfectly good
monotone response and Pearson penalises it. Measured on the identical eigengenes
at the identical thresholds:

| | Pearson | **Spearman** | in both | Spearman only | Pearson only |
|---|---|---|---|---|---|
| sugarcane (n = 48) | 408 | **465** | 394 | 71 | 14 |
| **purple (n = 18)** | **38** | **79** | 34 | **45** | 4 |

**Purple's responsive set more than doubles.** That is the single largest effect
of this change, and it is exactly where the ordinal argument predicts it: purple
is the study with three levels. Sugarcane's trait is two-level, where Spearman is
the rank-biserial correlation, so the gain there is smaller (and is mostly
robustness to the outliers a PC1 can carry).

> **What a monotonic test can and cannot ask of purple.** Purple's three levels
> are not a dose-response gradient: **2 N is the control, and 0 N and 6 N are
> stresses in opposite directions** (deficiency and excess). So Spearman on this
> coding asks *does the module track nitrogen supply?* — a real question, but not
> *does it respond to nitrogen stress?* A module moved the same way by both
> stresses is invisible to it, and that is the shape the design predicts.
>
> The cost is measured rather than assumed. A U-shape contrast `c(+1, -2, +1)`
> over the three levels finds **one** significant purple module at padj ≤ 0.05,
> and Spearman misses that one; Spearman finds 79. So the monotonic test is not
> leaving a large non-monotonic set behind — at n = 18 the U-shape test has almost
> no power, the same wall everything else in purple hits. **The 79 should be
> described as tracking nitrogen supply, not as stress responders.**

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

### GO by response direction — the grain that was missing

The per-module test above describes 49 of 465 modules. Pooling every responsive
module into one set — the other grain the pipeline had — answers a question
nobody asked, because the two response directions are opposite in expression and
there is no reason to expect them to be the same in function.

So the modules that **rise** with nitrogen are now tested as one gene set against
the modules that **fall**, on the same `topGOdata` object, the same background and
the same raw-p threshold (`./run.sh modulego <study> [ONT]`, written to
`module_GO_<ONT>_<study>_bydirection.tsv`).

**Pooling dissolves the annotation gate.** Modules too small to test alone all
contribute:

| | rises with N | falls with N |
|---|---|---|
| sugarcane: modules / genes / GO-annotated | 242 / 3,211 / **383** | 223 / 1,613 / **161** |
| sugarcane: BP terms at raw p ≤ 0.05 | **80** (34 clear BH) | **29** (4 clear BH) |
| purple: modules / genes / GO-annotated | 32 / 202 / 23 | 47 / 460 / 47 |
| purple: BP terms | 10 | 26 |

**The two directions are almost disjoint in function.** Of sugarcane's 103
enriched BP terms, **74 are found only among the modules that rise, 23 only among
those that fall, and 6 in both**; purple gives 7 / 23 / 3. The split is not one
programme cut in half, and pooling the directions — which is what the previous
analysis did — averages two distinct programmes into one list.

**And what they are is interpretable.** Sugarcane:

| rises with nitrogen | falls with nitrogen |
|---|---|
| **nitrate assimilation** (p = 3.2e-09) | flavonoid biosynthesis (2.1e-05) |
| **nitric oxide biosynthesis** (3.2e-09) | raffinose-family oligosaccharide biosynthesis (1.6e-04) |
| **proline biosynthesis** (1.7e-12) | triglyceride biosynthesis (9.1e-04) |
| spermidine biosynthesis (8.4e-08) | cellular response to cold (5.7e-06) |
| reactive oxygen species biosynthesis (3.2e-09) | skotomorphogenesis (7.3e-05) |
| response to fungus / bacterium / chitin (1.7e-14 …) | microtubule bundle formation (2.0e-05) |
| monoterpene biosynthesis (8.5e-13) | flowering-time regulation (2.5e-04) |

Nitrogen assimilation and the metabolism that consumes it go **up**; the
low-nitrogen carbon-storage and anthocyanin-type programme goes **down**. That is
the textbook shape of a nitrogen response, read at module resolution, and
**neither of the other two grains recovers it** — the per-module grain sees 11%
of the modules and the pooled grain averages the two directions away.

Purple's direction sets are much smaller (23 and 47 annotated genes) and give a
different picture: cell-wall construction up (plant-type primary cell wall
biogenesis and cellulose biosynthesis, both p = 2.3e-08), phenylpropanoid
biosynthesis, phosphate and anion transport and shoot morphogenesis down. Read
that as a report of what those sets contain rather than as purple's nitrogen
response.

> **The direction labels mean "tracks nitrogen supply", not "responds to nitrogen
> stress".** In purple the three levels are stress-control-stress, so a monotone
> test cannot see a module moved the same way by both extremes. See
> [Module 20 in purple](#module-20-in-purple--one-copy-responds-and-not-monotonically)
> for a gene that does exactly that.

---

---

## The module analysis, rebuilt: Pearson-only graphs, per-species inflation, blocked trait test

Two changes land together: the clustering moves to each species' own modularity
optimum on the unpruned Pearson-only graphs, and the module-trait test gains the
experimental design.

### The clustering, adopted rather than recomputed

The inflation ladder had already written the partition at every setting, so the
two chosen cells were adopted directly (`51_mcl_membership_from_cls.r`, following
`27_sbm_membership.r`'s contract) instead of re-running mcl — which would in any
case need the 9-column edge table these networks no longer have.

| | sugarcane `-I 1.5` | purple `-I 3.5` |
|---|---|---|
| clusters | 5,653 | 29,624 |
| named modules (≥ 2 genes) | 5,650 | 14,024 |
| largest | 23,439 (22.98%) | 39,230 (23.06%) |
| Unassigned | 3 | 15,600 (9.17%) |
| modularity Q | 0.08212 | 0.15324 |
| eigengenes (≥ 3 genes) | 3,627 | 7,493 |
| PC1 variance, median | 80.8% | 88.4% |

Cluster count and largest module were cross-checked against `mcl_sweep_<study>.tsv`
at that inflation, and every gene is accounted for (101,990 and 170,135 exactly).
Because an mcl cell tag strips the decimal point — `-I 1.5` and `-I 15` would share
a file — the adopter refuses any cell whose `.inflation` sidecar disagrees with the
value requested. That guard was tested firing.

### The design goes into the module test

`53_module_trait_blocked.r`. The project's own QC puts genotype at R² = 0.999 of
PC1 in purple and 0.998 in sugarcane, and leaf segment at 0.802 of PC2 — so the
largest variance component in either matrix had been sitting in the residual of
every module-level nitrogen test.

```
sugarcane   eigengene ~ genotype + segment + N      n = 48, residual df 42
purple      eigengene ~ genotype + N                n = 18, residual df 15
```

Spearman stays primary (both traits are ordinal or two-level; the argument is in
the superseded section above and is unchanged), implemented as the blocked fit on
midranks. Blocked Pearson is written beside it.

| | sugarcane | purple |
|---|---|---|
| modules tested | 3,627 | 7,493 |
| responsive, **marginal** | 251 | 96 |
| responsive, **blocked** | **588** (16.2%) | **182** (2.4%) |
| in both | 251 | 96 |
| up / down | 276 / 312 | 65 / 117 |
| null mean (1,000 within-block permutations) | 0.14 | 0.39 |
| permutations reaching the observed count | **0 / 1,000** | **0 / 1,000** |

**Blocking is a strict superset in both species — nothing the marginal test found
is lost.** That is the result that says the blocks absorbed noise rather than
signal, and it is the same behaviour the change produced at gene level. Neither
count is reachable by chance.

Three checks the numbers rest on, all fatal if they fail:

* the vectorised solver reproduces `lm()` to **4.4e-16** (t) and **2.2e-16** (p);
* the model matrix is rank-checked, which is what a block confounded with the
  trait would look like, and the unblocked Spearman path reproduces
  `cor.test(method = "spearman", exact = FALSE)`;
* **sugarcane's 48 libraries are 12 plants × 4 leaf segments — repeated measures,
  not 48 replicates.** Segment as a fixed block removes the segment means but not
  the within-plant correlation, so a plant-level run (n = 12, residual df 9) is
  computed alongside. It agrees with the blocked fit at **r = +0.9486**.

The null permutes **within block**. Shuffling labels freely would break the
structure the model conditions on and give an optimistic null.

TF enrichment among responsive modules reproduces the earlier asymmetry:
sugarcane 3.57% against 1.51% in non-responsive (OR 2.41, p = 0.00207); purple
0.55% against 0.16% (OR 3.36, p = 0.274, one module — not significant).

### One error caught, and the trap it came from

The first version of 53 carried the marginal rho by joining the **previous run's**
`module_trait_<study>.tsv` on `module`. Module names are positional — `Module_%03d`,
largest first — so that joined different gene sets across two clusterings: old
`Module_001` held 19,604 sugarcane genes, the new one holds 23,439. It ran cleanly
and reported that 244 marginal calls had been "lost under the design". The marginal
fit is now computed on the same eigengenes, and the true answer is the opposite —
blocked is a strict superset. `results/STALE.md` records the trap.

### What is stale

`results/` now holds Pearson-only modules beside the **merged** network's edge
table, because the Pearson-only edge tables were deleted (70 GB whose only consumer
was `mcxload`; the live graph is the `.mci`). So the topology figure, both
conservation stages, `08_conserved_cor_genes.r` and the gene-level trait results
still describe the merged network at `-I 2` and must not be quoted against the new
modules. Listed in `results/STALE.md`; module GO is pending its topGO run.


## The gene analysis, rebuilt: blocked nitrogen correlation on the Pearson-only nodes

The module level moved to the Pearson-only graphs with the design in the model; the
gene level had not moved with it, and `results/STALE.md` said so. It has now, with
three changes that land together.

### One gene universe, and it is the network's

Every gene-level test now corrects over the **node set of the unpruned Pearson-only
network** — 101,990 sugarcane and 170,135 purple genes, the same sets the modules are
built on (`node_list()` in `config.sh`, `restrict_to_universe()` in `lib/common.R`).
It replaces `conserved_genes_<study>_FULL.txt`, the genes carrying a conserved edge in
the **merged** graph, which is what the first blocked run corrected over: 39,226 and
44,118 genes. The BH denominator therefore rose 2.6× and 3.9×, so **the counts below
are not comparable with the ones that table used to hold** (4,137 / 2,331).

The two universes are asymmetric and that is a property of the data, not a choice:
purple's network holds 99.6% of its VST genes, sugarcane's only 59.7%.

### The design goes into the gene test

`31_gene_trait_blocked.r` (brought forward from branch `blocked-gene-trait`, commit
`052fcd6`, and refactored onto `lib/common.R`'s shared solver so the gene and module
levels cannot drift apart):

    sugarcane   expr ~ genotype + segment + N            n = 48, resid df 42, Pearson
    purple      expr ~ genotype + N       on midranks    n = 18, resid df 15, Spearman

Purple is Spearman because its trait is an **ordinal 0/2/6 mM dose** and Pearson reads
that spacing literally — the argument already accepted at module level, and recorded
in `Open items` as never having been applied at gene level. It is now. Sugarcane's
trait is two-level, where Spearman on midranks is the same test up to a monotone
relabelling, so it keeps Pearson; both statistics are computed and written either way.

| | sugarcane | purple |
|---|---|---|
| genes tested | 101,990 | 170,135 |
| responsive, **blocked** | **8,737** | **3,854** |
| responsive, marginal, *same genes* | 3,238 | 882 |
| \|r\| at the BH boundary | 0.366 | 0.720 |
| plant-level control | resid df 9, agrees at **r = +0.9761** | n/a |

Blocking is a **strict superset** of the marginal rule in both species — no marginal
call is lost, so the block is removing noise rather than absorbing signal. The
marginal column is computed **in the same run on the same rows**, never read from
`07`'s table, which corrects over a different universe; joining the two would
attribute a change of denominator to the change of model. That is the mistake
`53_module_trait_blocked.r` made at module level and it is not repeated here.

The diagnostic that motivated the whole stage, p-value histogram by decile:

    purple  marginal   14.7  8.9 13.3  8.2 11.2  7.2  8.4  8.4  7.2 12.5
            blocked    24.7 10.8 12.8  9.0  8.6  6.3  6.3  5.4  5.8 10.2

The **rising last decile is the symptom** — a mixture of a uniform null and real
signal cannot rise. Blocking removes most of it (12.5% → 10.2%) but **not all**: the
blocked histogram is still not monotone, so something remains unmodelled in purple at
n = 18. Worth stating rather than declaring the problem solved. Sugarcane's falls from
6.6% to 4.1%.

### The non-monotone tier: included, and it contributes one gene

Purple has three nitrogen levels, so a gene moved the same way by deficiency and
excess is invisible to every monotone test. `30_gene_trait_ushape.r` (the blocked
quadratic contrast, `expr ~ genotype + nlev`, resid df 14) now runs on the **same**
gene universe, so its BH and the monotone test's are over identical genes and the two
selections can legitimately be unioned — `61` asserts that.

Over 170,135 genes it selects **1** gene at padj ≤ 0.05 (a peak at the control).
Raw p ≤ 0.05 gives 17,507 against 8,507 expected, a 2.06× aggregate excess, so
non-monotone signal is *present* and the design cannot resolve which genes carry it.
That one gene's ortholog is not sugarcane-responsive, so the tier contributes **0**
conserved pairs. This is the same shape as the module-level answer, and it closes the
open item that said the U-shape had never been measured genome-wide at gene level.

`Soffic.09G0001580-9H`, the anecdote, lands at raw p = 0.0037, padj = 0.484, rank 890
of 170,135 (0.52nd percentile) — a strong trough, nowhere near surviving correction.

### The conserved nitrogen response, at node level

`61_conserved_blocked_nodes.r`. `08_conserved_cor_genes.r` is left alone as the record
of the merged network under the marginal `pearson|mi|union` rules; it has an edge level
this rebuild has no inputs for, since conserved-edge tables describe the merged graph.

The universe is ortholog pairs whose **both** sides are nodes of their own network:
275,047 unrestricted pairs → **114,681** over 43,390 orthogroups.

| | genome-wide | directed (discover in sugarcane) |
|---|---|---|
| responsive sugarcane / purple | 8,737 / 3,855 | 8,737 / 492 of 8,725 candidates |
| conserved correlated **pairs** | **392** over 270 orthogroups | **616** over 419 |
| null (ortholog shuffle) | 314.3 ± 16.5, **fold 1.25**, p 0.000999 | 573.8 ± 13.8, fold 1.07, p 0.002 |
| sugarcane **genes** in both | **326** | 502 |
| null | 301.9 ± 15.8, **fold 1.08**, p **0.071** | 536.3 ± 13.3, fold **0.94**, p 0.995 |
| sign concordance, pairs | 244/392 = 62.2%, p 1.4e-06 | 373/616 = 60.6%, p 1.8e-07 |
| sign concordance, **orthogroups** | 160/269 = 59.5%, p **0.0022** | 241/418 = 57.7%, p **0.0020** |

**What this establishes, and what it does not.** The *gene* count does not beat its own
null — 1.08× at p = 0.071 genome-wide, and 0.94× (below the null) under the directed
design. The excess in the *pair* count is therefore mostly ortholog multiplicity, not
more conserved genes: `OG0001395` alone contributes 6 pairs from 3 sugarcane × 2 purple
genes. **Counting pairs and calling them findings would overstate the result** and the
pair row must never be quoted without the gene row beside it.

What does hold up is **direction**. Of the concordant pairs 94 are up in both species
and 150 down in both, and the agreement survives the pseudoreplication correction: one
vote per orthogroup, majority sign, ties abstaining — 59.5% of 269 at p = 0.0022, and
57.7% of 418 under the directed design. Weaker than the pair-level p of 1e-06 as it
should be, and significant in both designs. The pair-level binomial is reported beside
it, not instead of it.

### Four causes of "not conserved", kept apart

Collapsing them would make a coverage gap look like biological absence. Of sugarcane's
8,737 responsive genes:

| status | n |
|---|---|
| `ortholog_not_correlated` | 5,558 |
| `no_ortholog` | **2,775** (31.8% — could never have been called conserved) |
| `conserved_correlated` | 326 |
| `ortholog_not_a_node` | 78 |

`ortholog_not_a_node` is new, and it exists because the universe is now the network
rather than the transcriptome: a gene can have a perfectly good ortholog that is simply
not in the other species' graph. Purple: 1,736 / 1,217 / 311 / 591.

## What hubs are for, and what the periphery is for

`63_degree_go.r` ranks every GO-annotated network node by degree and tests each GO term
for concentration at one end, by a Kolmogorov–Smirnov statistic under topGO's `weight01`
algorithm. The ranking is used whole: degree spans four orders of magnitude (sugarcane's
10th percentile is 2, median 53, 90th 6,677), so any decile cut is arbitrary and
discards the middle. `weight01` is kept rather than moving to `fgsea` because it
decorrelates the GO DAG — without it a parent and its children score on the same genes
and the table fills with near-duplicates.

**The two species agree**, across networks built from unrelated experiments:

| | hubs | median degree / background |
|---|---|---|
| both | photosynthesis, light harvesting | 5.7× / 9.8× |
| both | photosynthesis | 5.1× / 10.7× |
| both | translation | 2.6× / 4.2× |
| both | mRNA splicing, intracellular protein transport, protein deubiquitination | 2.6–7.4× |

| | periphery | median degree / background |
|---|---|---|
| sugarcane | regulation of DNA-templated transcription | 0.55× |
| both | RNA modification | 0.52× / 0.41× |
| both | protein phosphorylation | 0.69× / 0.76× |
| both | hydrogen peroxide catabolism | 0.77× / 0.49× |
| sugarcane | trehalose biosynthesis | **0.15×** |

The densely connected core of a co-expression network is the housekeeping machinery;
signalling, regulation and specialised metabolism sit at its edge.

**An effect size is reported beside every p, and it is not decoration.** Over 1,711
tested terms a KS test reaches significance on small shifts: sugarcane's *carbohydrate
metabolic process* clears p = 2.6e-08 with its genes' median degree at **1.04×** the
background — no shift at all — while *trehalose biosynthesis* sits at 0.15×. The figure
encodes the ratio as point size, so a term that is significant but flat is visibly small.

### This is the same contrast as figure 4 panel B, reached another way

Genes on a conserved edge **are** hubs:

| | median degree, on a conserved edge | no conserved edge |
|---|---|---|
| sugarcane | 228 | 18 |
| purple | 6,916 | 376 |

63.3% of sugarcane's top-degree decile sits on a conserved edge against 7.1% of its
bottom. **Part of that is arithmetic**: at a 10.4% per-edge conservation rate, one
conserved edge out of 6,677 is near-certain and out of 2 is unlikely.

So the conserved/non-conserved GO split and the hub/periphery GO split are **entangled,
and neither is independent evidence for the other**. Both legends now say so, with the
numbers. What the pair does establish is that the housekeeping-core/regulatory-periphery
organisation is visible from two directions in two species; what it does not establish is
that cross-species edge conservation selects housekeeping genes *over and above* their
being hubs. Separating those would need a degree-matched comparison, which is not done
here.

One caveat measured rather than assumed: GO coverage is mildly degree-dependent — 59.7%
in sugarcane's bottom degree decile against 66.8% in its top (62.6% / 65.4% in purple) —
so the periphery is slightly the less well annotated end.

## Edge-level conservation, on the Pearson-only networks

The node level moved to the Pearson-only graphs with `61`; the **edge** level had not,
so every edge-conservation number in this project still described the merged
(Pearson + MI) graph. `62_conserved_edges_pearson.py` closes that, in both directions
and for the nitrogen-correlated subsets. `06_conservation_join.r` and
`13_conservation_null.r` are untouched and remain the record of the merged analysis.

### The input, since the edge tables were deleted

`06` streams `network_<study>_edges.tsv`, and the Pearson-only versions were deleted —
70 GB whose only consumer was `mcxload`. Nothing had to be regenerated: `mcxdump` had
already written `<study>.pairs` for the Leiden sweep, one row per undirected edge as
`idx1 idx2 weight`, and both files verify **exactly** against the matrices
(75,333,769 and 675,955,918 edges). The script asserts that before doing anything, so
a stale dump cannot quietly shrink a network and depress every rate.

Everything runs in integer index space — gene names are materialised only when
conserved edges are written — which is what makes 676 M edges against a 676 M-edge
adjacency affordable: one sorted `int64` key array and a binary search.

### Both directions

| | sugarcane → purple | purple → sugarcane |
|---|---|---|
| edges | 75,333,769 | 675,955,918 |
| conserved | **7,806,379 (10.36%)** | **10,196,172 (1.51%)** |
| fold over null | **1.69** | **1.71** |
| empirical p | 0.0099 | 0.0099 |
| genes on a conserved edge | 37,867 | 42,194 |

These land where they should: the merged network's *Pearson-layer* rows were 10.69%
and 1.50%, and the graphs differ by only 1.1% and 4.2% of edges.

**The two directions are not comparable to each other** and never should be. Purple's
density is 3.2× sugarcane's (0.0467 against 0.0145), so a random ortholog pair is far
likelier to land on an edge there. Each is read only against its own null.

**The fold is lower than the merged network's 2.57, and that is a stricter null, not a
weaker signal.** This null draws only from ortholog pairs whose **both** sides are
network nodes — the same 114,681-pair universe `61` uses, asserted equal at run time.
The old null drew from all 275,047 pairs including genes that could never have been
matched, which inflates the fold by giving the permutation impossible assignments to
fail at.

### Conservation rises with edge strength, in both species

This replaces the old by-layer breakdown ("are MI edges conserved as often as Pearson
edges?"), which is vacuous on a single-layer graph. Rank-based deciles, so each bin
holds a tenth of the edges:

| decile | sugarcane → purple | fold | purple → sugarcane | fold |
|---|---|---|---|---|
| D1 (weakest) | 9.72% | 1.57 | 1.17% | 1.39 |
| D5 | 10.20% | 1.66 | 1.39% | 1.59 |
| D10 (strongest) | **11.30%** | **1.85** | **2.16%** | **2.30** |

**Monotone across all ten deciles in both directions, in rate and in fold.** The fold
rising is the part that matters: if strong edges were merely joining better-annotated
genes, the rate would rise and the fold would not. Purple's effect is the larger one —
its strongest decile is conserved 1.8× as often as its weakest.

Weight ranges are reported per bin because the [0.01, 1] rescaling is per-study: a
decile in sugarcane and a decile in purple do not stand for the same |r|.

### The nitrogen-correlated edges — and what the null does to them

Who is nitrogen-correlated is read from `61`, never recomputed, so the node and edge
levels cannot drift.

| | edges | conserved | rate | fold | p |
|---|---|---|---|---|---|
| **sugarcane → purple** | | | | | |
| both ends responsive in sugarcane | 407,700 | 44,114 | 10.82% | 1.44 | 0.0099 |
| both ends responsive in **both species** | 945 | 295 | **31.22%** | **1.58** | **0.0198** |
| **purple → sugarcane** | | | | | |
| both ends responsive in purple | 849,681 | 12,615 | 1.49% | 1.27 | 0.0198 |
| both ends responsive in **both species** | 5,692 | 337 | **5.92%** | **1.07** | **0.37** |

Read the raw rates alone and this looks like a 3–4× enrichment in both species: edges
joining two genes that respond to nitrogen on both sides are conserved at 31.2% and
5.92% against backgrounds of 10.4% and 1.5%.

**Against the proper null, most of that is orthology, and in purple all of it is.**
Genes responsive in both species are by construction genes with good orthologs, and
the ortholog shuffle prices that in. Sugarcane keeps a modest real excess (1.58×,
p = 0.02); purple's vanishes entirely (1.07×, p = 0.37). The honest statement is that
**a shared nitrogen response predicts edge conservation in sugarcane and not in
purple**, and that the eye-catching raw rates are mostly a selection effect.

### A defect in the first version of this null, and why the fix changed the answer

The first run computed every stratum's null on a 5 M-edge Bernoulli sample. That is
right for a 676 M-edge network and **wrong for a rare stratum**: the sample caught 46
of purple's 5,692 `resp_both` edges, none of which happened to be conserved, so the
null divided by nothing and reported **fold 0.0, p 1.0** against a true rate of 5.92%.
The null now runs on the **complete** stratum for both responsive sets (5,692 and
849,681 edges — trivially affordable) and samples only for `all`; every output row
carries a `null_basis` column naming the set its fold came from, so a sampled fold and
a complete-stratum fold cannot be compared by accident.

It was not only purple that moved: sugarcane's `resp_both` fold was 1.92 on the
sampled null and is 1.58 on the complete stratum. The sampled null was optimistic in
both species.

### What is written

Only **conserved** edges, with their weights — 7.8 M and 10.2 M rows, ~1 GB. The old
schema wrote every edge with a `conserved` TRUE/FALSE column (~36 GB here) and every
consumer immediately filtered to TRUE; the totals those FALSE rows carried are in the
summary.

Outputs carry a `_pearson` suffix. Note that `08_conserved_cor_genes.r` also writes
`*_pearson.tsv` files, where the word means its **pearson selection rule** on the
merged graph; here it means the **Pearson-only network**. The base names differ
(`conserved_correlated_*` for `08`, `conserved_edges_*_to_*` and
`conservation_summary_*_to_*` here) so no file is overwritten, but the two must not be
read as one series. Nothing here touches a `_FULL` path — `conserved_genes_<study>_FULL.txt`
is still read by `run.sh trait`, `08` and `09_go_enrichment.r`, and overwriting it
would have made three superseded stages describe a mixture of two graphs.

### What transfers is housekeeping; what does not is regulation

`09_go_enrichment.r` now runs on either side of a partition: the genes on at least one
conserved edge, and the **exact complement** — network nodes with no conserved edge.
They partition the node set (sugarcane 37,867 + 64,123 = 101,990; purple 42,194 +
127,941 = 170,135) and share a background, so this is one gene set split two ways
rather than two unrelated tests.

| conserved edges | no conserved edge |
|---|---|
| translation | regulation of DNA-templated transcription |
| intracellular protein transport | response to auxin |
| protein folding | protein ubiquitination |
| vesicle-mediated transport | positive regulation of transcription by RNA pol II |
| photosynthesis, light harvesting | mitotic cell cycle phase transition |
| mRNA splicing, via spliceosome | ethylene-activated signalling pathway |

*(top 6 BP terms per set, ranked by the worse of the two species' p-values, so each is
a term both species agree on. 111 shared terms in the conserved set, 61 in the
complement.)*

The co-expression that survives between the two species is **core cellular machinery**;
what does not is the **regulatory layer** — transcriptional control, hormone signalling,
ubiquitination, the cell cycle.

**One caveat that is structural, not incidental.** A gene with no ortholog at all cannot
have a conserved edge, so it lands in the complement by construction — 2,775 sugarcane
and 1,217 purple genes among the nitrogen-responsive alone. The regulatory signal is
therefore partly a statement about which gene families have clean orthology between
these two genomes, not purely about which are conserved in co-expression. Transcription
factors and hormone-signalling families are exactly the ones that expand and diverge, so
the two explanations are not separable here.

**The literal alternative was degenerate and was not used.** "Genes on at least one
non-conserved edge" is 101,256 of 101,990 sugarcane nodes — 99.3% — because at a mean
degree near 1,477 almost every gene has some non-conserved edge. It overlaps the
conserved set by 37,133 genes and would have been tested against a background that is
essentially itself.

### The join was verified independently

A fast vectorised join that is subtly wrong returns plausible numbers — that is how
the k-NN result went wrong once. So the join was checked by a different route
entirely: resolve gene **names** through a fresh parse of `Orthogroups.tsv`, build
candidate index pairs as plain tuples, and confirm membership with one `awk` pass over
the raw 14.5 GB `purple.pairs`.

    POSITIVE  4,000 of 4,000 conserved edges confirmed
    NEGATIVE  0 of 4,000 non-conserved edges wrongly match

The negative case is the one with teeth: 4,000 edges the script called *not*
conserved, all with both endpoints mappable through orthology, and none of their
36,011 candidate pairs is a real purple edge. A join that over-matches passes the
positive test and fails that one. The check also reproduced 61,846 mappable sugarcane
genes independently.

## The module enrichment analysis, on one annotation

GO is now derived from **one source**: a full local InterProScan 5.78 over all 17
member databases, mapped to GO through the GO Consortium's pinned `interpro2go` and
`pfam2go`, then normalised to most-specific terms (`60_normalise_gene2go.r`). Four
sources were built and scored; three are not used. That is a result, not a shortcut,
and the reason is measured.

### Why one source and not four

The judge is `55_go_coherence.r`: Sørensen–Dice annotation homogeneity above a
size-matched null, asking whether co-expressed genes share function more than chance.
Scored on a **fixed gene set** — identical genes, identical modules — so the only
thing varying is the annotation:

| annotation | genes scored | terms | H | H null | **excess** |
|---|---|---|---|---|---|
| nfcore5db (InterPro, 5 DBs) | 34,056 | 2,005 | 0.0978 | 0.0223 | **+0.0755** |
| **ipsfull (InterPro, 17 DBs)** | 34,056 | 2,114 | 0.1009 | 0.0261 | **+0.0748** |
| union (everything merged) | 34,056 | 7,936 | 0.2560 | 0.1889 | +0.0671 |
| eggnog_auto | 34,056 | 7,719 | 0.2617 | 0.1967 | +0.0650 |

Purple orders identically: ipsfull +0.0361, eggNOG +0.0316, union +0.0264.

**Merging sources measurably lowers coherence, in both species.** The union scores
below either InterPro table on its own, so a tiered or merged annotation would buy
negative signal at the cost of a `source` and `tier` column on every pair and a
provenance caveat on every result. There is no merged table for that reason.

Two things in that table are worth not misreading. First, **eggNOG's raw H is 2.6×
ipsfull's and its excess is lower** — without the size-matched null it would have
looked like the best annotation by a wide margin. What the raw H measures is term
density (16.7 terms/gene after normalisation, against 2.1), not shared function.
Second, **the two InterPro tables are tied**: 0.9% apart, inside the noise. ipsfull
was adopted for **reach**, not coherence — on the same genes it is no more coherent,
but it annotates 7,204 more sugarcane and 15,094 more purple network genes, and that
is what turns into testable modules.

| | GO on network genes | | testable responsive modules | | terms clearing cross-module BH | |
|---|---|---|---|---|---|---|
| | 5 DBs | **17 DBs** | 5 DBs | **17 DBs** | 5 DBs | **17 DBs** |
| sugarcane | 56,974 (55.9%) | **64,178 (62.9%)** | 442 | **479** | 143 | **164** |
| purple | 94,497 (55.5%) | **109,591 (64.4%)** | 106 | **118** | 8 | **9** |

Coverage is quoted per **network node set**, never per proteome — the two differ
(194,593 sugarcane proteins against 101,990 network nodes) and mixing them is how
these numbers got muddled earlier in this work.

The gain from 8.0%/7.2% (the original eggNOG GO column) to 62.9%/64.4% is the whole
arc of this stage. It is worth stating what it did **not** buy: the first large jump
raised testable modules 4.6× and yet *fewer* terms cleared correction, because the
test family grew from 229,140 to 996,268 module × term pairs and the old 8%
background was biased toward well-studied genes, which had inflated apparent
enrichment. Coverage and significance do not move together, and this section reports
both every time.

### Three ontologies, on the adopted annotation

`./run.sh modulego <study> [BP|MF|CC]`, one topGO run per ontology per species
against the same network-node background, all on the same 588 and 182 responsive
modules from the blocked test:

| | ontology | modules with ≥ 1 term | terms written | **clearing cross-module BH** |
|---|---|---|---|---|
| **sugarcane** | BP | 425 of 479 | 1,233 | **164** |
| | MF | 457 of 479 | 1,766 | **272** |
| | CC | 175 of 479 | 312 | **22** |
| **purple** | BP | 109 of 118 | 291 | **9** |
| | MF | 110 of 118 | 375 | **17** |
| | CC | 44 of 118 | 70 | **1** |

The annotation gate is identical across ontologies by construction — 479 and 118
testable modules — because it counts genes with *any* GO, not genes with a term in
that ontology.

**MF outperforms BP**, and by a lot: 272 terms clearing FDR against 164 in sugarcane,
17 against 9 in purple. That is the expected shape for a domain-derived annotation —
InterPro signatures name what a protein *does* far more sharply than what process it
participates in — and it is worth knowing that the ontology this project reports as
primary is not the one with the most power. BP stays primary because the question is
about a nitrogen *response*, and figure 6 draws BP for that reason; MF and CC are
reported beside it, not promoted.

**CC is weak in both species and near-empty in purple** (one term). Nothing should be
built on it.

### The honest negative, again

Purple clears cross-module BH on **9 BP terms across 118 testable modules**. More
annotation did not fix that: it went 14 → 8 → 9 across three successive annotations
while purple's coverage went 7.2% → 55.5% → 64.4%. Sugarcane's 164 is not the whole
picture and quoting it alone would misrepresent the comparison. The constraint in
purple is n = 18 and a network that is one dense component, not the annotation.

## The paper figures

Seven figures, each generated by one script and each carrying a **generated
legend** — written by the script that draws it, from the same variables, so a
number cannot disagree between a figure and its caption. `./run.sh legends`
assembles them into `figures_legends.txt` at the repo root. Never hand-edit that
file; re-run the figure, then `legends`.

| | stage | what it shows |
|---|---|---|
| **1** | `figdataset` | the two designs, library QC, the gene funnel, PCA per study |
| **2** | `figrepro` | each source study's own finding, reproduced here |
| **3** | `figtopology` | degree and module-size CCDFs, and what hubs vs the periphery are enriched for |
| **4** | `figconservation` | conservation vs null (both directions), conserved-set GO, conservation vs edge strength, the funnel |
| **5** | `figmodules` | module selection, the Spearman gain, the responsive eigengenes |
| **6** | `figmodulego` | the annotation gate, then GO by response direction per species |
| **7** | `figmodule20` | Module 20 in sugarcane with its AtMYB59 copies marked, in purple with a TF column, and the one copy with a U-shape response |

Figures carry a panel letter and the labels the data needs, and nothing else — no
titles, no subtitles, no statistics printed on the panel. Everything else is in
the legend. Figure NUMBERS live in `config.sh` (`FIG_*`), so renumbering the
paper is one edit.

Two results in Figure 1 panel D are worth carrying into any reading of the rest:
**PC1 is genotype in both studies** (R² = 0.998 and 0.999) and **nitrogen loads
on neither first component**. The nitrogen response is real but it is nowhere
near the dominant structure in either dataset, and it has to be found against
genotype and tissue rather than read off a leading component.

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

## Module 20 in purple — one copy responds, and not monotonically

Muñoz-Perez et al. build their nitrogen claim on Module 20. Its response
reproduces in their own data (see the reproduction figure above). This is the
comparative half: what is left of it in purple.

**The 12 published members are not 12 genes.** They are TransDecoder transcripts
from a de novo assembly; mapped into each reference proteome by reciprocal
DIAMOND blastp through Arabidopsis they resolve to **three anchors**, and the
mapped genes are haplotype copies of those few loci — 33 sugarcane genes from 10
members, 27 purple genes from 6.

| Arabidopsis anchor | sugarcane copies (MYB by our rules) | purple copies (MYB) |
|---|---|---|
| `AT5G59780` **AtMYB59** | 9 (**9**) | 3 (**2**) |
| `AT3G55960` uncharacterised | 14 (0) | 12 (0) |
| `AT4G10770` **OPT7** transporter | 10 (0) | 12 (0) |
| `AT3G46130` **AtMYB48** | 0 | 0 |

Only the AtMYB59 anchor carries a MYB call under our own Pfam/PlnTFDB rules. The
two anchors that supply most of the mapped genes are not transcription factors at
all, whatever the source study labelled the transcripts — and its own
classification reproduces 12/12 on its own proteins, so the disagreement is about
the orthologs we map to, not about how TFs are called.

**One purple copy is expressed, and it responds in a U.**
`Soffic.09G0001580-9H` is the only AtMYB59-anchor copy expressed in purple leaf —
60.6 TPM against 0.01 and 0.00 for the other two — and is called MYB both by the
source study and by our pipeline.

| 51NG3 | 0 N (starvation) | 2 N (control) | 6 N (excess) |
|---|---|---|---|
| mean TPM | **124.1** | **10.8** | **47.4** |
| TAGZ | 65.5 | 40.0 | 75.6 |

An **11-fold drop** from starvation to the control and a **4-fold rise** again
beyond it (U-shape contrast `c(+1, −2, +1)` p = 0.0078, one-way p = 0.0126).
Because 2 N is the control and 0 N and 6 N are stresses in opposite directions,
this is a gene **induced by nitrogen stress in either direction**.

**That shape is invisible to every monotone test in this project.** It is why the
module carrying this gene is not in purple's 79 responsive modules, and it is the
concrete case for the U-shape contrast existing as a separate test.

**Sugarcane sees one arm of the same curve.** The nine expressed copies at the
same anchor are every one monotonically repressed by nitrogen (log2FC −1.8 to
−3.4, padj ≈ 1e-08 in the responsive genotype RB975375). Sugarcane's design has
no control level, only Low and High — so what it *can* observe is the
low-nitrogen arm. Read together, the two studies are **consistent rather than
contradictory**: only purple's three-level design shows that the high-nitrogen
side turns back up.

> **The limit.** This rests on one expressed copy in one species, at a p-value
> that does **not** clear BH over the 27 purple Module-20 genes tested. It is a
> candidate worth following, not a finding. It also sits against the module's
> overall verdict, which is unchanged: Module 20's *nitrogen response* transfers,
> its *MYB identity* barely reaches our references, and its *network position*
> does not transfer at all.

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
- ~~**A U-shape test is run nowhere except the two readouts.**~~ **CLOSED.** It now
  runs genome-wide at gene level on the network universe, as a second test family
  unioned with the monotone one. The cost is measured in both places and is almost
  nothing: **one module** and **one gene**, and that gene's ortholog is not
  sugarcane-responsive, so it adds no conserved pair. `Soffic.09G0001580-9H` sits at
  padj 0.484, the 0.52nd percentile — a real trough the design cannot resolve. The
  aggregate excess (17,507 raw hits against 8,507 expected, 2.06x) says the signal
  exists and n = 18 cannot localise it.
- ~~**The gene level is still Pearson + MI; only the module level moved to
  Spearman.**~~ **CLOSED.** The gene level is now blocked, Spearman for purple, on the
  Pearson-only network's node universe, and it did matter: 3,854 responsive purple
  genes against 882 under the marginal rule, a strict superset. The conserved-response
  answer was rebuilt with it at **node** level (`61_conserved_blocked_nodes.r`) and at
  **edge** level (`62_conserved_edges_pearson.py`), the latter reading the mcxdump
  edge stream since the Pearson-only edge tables were deleted. `conscor`'s numbers
  remain the merged-network record and are marked stale.
- **Why does conservation rise with edge strength?** Monotone across all ten weight
  deciles in both directions, in rate AND in fold over null (sugarcane 1.57 → 1.85,
  purple 1.39 → 2.30). The fold rising rules out the simplest explanation — that
  strong edges merely join better-annotated genes. Whether it reflects stronger
  selective constraint on tight co-expression, or the weakest decile sitting nearer
  the |r| ≥ 0.8 threshold where edge membership is noisiest in both species, is not
  answered here and is a real question.
- **Purple's blocked p-value histogram is still not flat.** Blocking took the last
  decile from 11.8% (old universe) / 12.5% (new) down to 10.2%, but the shape dips to
  5.4% and then rises again. Genotype plus nitrogen does not account for everything at
  n = 18, and any purple gene-level claim should be read with that in mind.
- **The conserved GENE count does not beat chance.** 1.08x at p = 0.071 genome-wide,
  0.94x under the directed design. Only the *pair* count does (1.25x), and that is
  ortholog multiplicity. The conservation claim rests on **direction** (59.5% of
  orthogroups, p 0.0022), not on how many genes respond in both species.
- **Purple's 79 responsive modules are not 79 independent findings.** Its
  permutation null reaches 244 in the worst of 1,000 shuffles because the network
  is one dense component. Any per-module claim from purple needs the module
  looked at, not just the padj.
