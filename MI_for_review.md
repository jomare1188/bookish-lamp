# Mutual information as a non-linear layer on the co-expression networks

*Method note for review — comparative Saccharum project (Muñoz-Perez 2025 "sugarcane",
n = 48 leaf libraries, R570 reference; Ta Quang Kiet 2025 "purple", n = 18 leaf
libraries, LA-purple reference).*

Pipeline: `scripts/mutual_information/` — see its `README.md` for how to run it.

---

## 0. The threshold, settled

Two correlation cut-offs appear in the code and only one of them built the
networks:

| where | value | what it is |
|---|---|---|
| `scripts/build_edgelist.r:22` | `PEARSON_THRESHOLD <- 0.7` | **pre-filter only**, produces the intermediate `edgelist_*_pearson.tsv` |
| `scripts/general_stats.r:11` | `PEARSON_MIN <- 0.8` | **the real cut**, produces `network_*_filtered_edges.tsv` |

Confirmed against the data rather than the code. Both network files are sorted
by `|r|` descending and their last rows are:

```
sugarcane   |r| = 0.8   p = 9.063e-12   weight = 0.01
purple      |r| = 0.8   p = 6.737e-05   weight = 0.01
```

`weight = 0.01` is the floor of the min–max normalisation, i.e. the weakest edge
retained. **The analysed networks are `|r| >= 0.8`.** Everything in this note is
matched to 0.8, and `config.sh:MATCH_PEARSON` has been set to 0.8 accordingly.

---

## 1. Why

Both networks are currently built from Pearson correlation. Pearson is fast,
robust and well understood, but it measures only **linear** co-variation. A gene
pair related by a saturating, threshold-like or non-monotone function is
invisible to it. In a nitrogen-response context — where saturation and dose
thresholds are the expected shape of a regulatory response — that is a
systematic blind spot, not a rare edge case.

Mutual information has no such restriction: it is zero if and only if two
variables are independent, whatever the shape of the dependence.

The purpose of this note is to explain (a) what the estimator actually computes,
and (b) exactly what the phrase "the MI threshold is matched to Pearson 0.8"
does and does not claim.

---

## 2. What the KSG estimator is

`--estimator ksg` is the **Kraskov–Stögbauer–Grassberger** estimator
(Kraskov, Stögbauer & Grassberger, *Phys. Rev. E* **69**, 066138, 2004),
"algorithm 1". Implementation: `custom_mutual_information.py:166`.

### 2.1 The problem it solves

Mutual information is defined through densities:

```
I(X;Y) = ∫∫ p(x,y) · log [ p(x,y) / (p(x)·p(y)) ]
```

With 48 samples those densities cannot be estimated directly. A 2-D histogram of
48 points is mostly empty bins and the resulting estimate depends strongly on
the arbitrary bin width. KSG never estimates a density at all. It uses
k-nearest-neighbour distances as an adaptive local ruler, chosen so that the
biases of the joint and marginal terms cancel against one another.

### 2.2 The algorithm

For one gene pair the data are `n` points — one per library — at coordinates
(rank of gene A in that library, rank of gene B in that library). For each
point *i*:

1. Find its **k-th nearest neighbour in the joint 2-D space** under the
   Chebyshev (maximum) metric, `d(i,j) = max(|xi−xj|, |yi−yj|)`. Call that
   distance `eps_i`. Geometrically this is an axis-aligned **square box** around
   point *i* containing exactly k neighbours.
2. Count `nx(i)` = points falling within `eps_i` **on the x-axis alone**,
   ignoring y — i.e. the vertical strip of that same width.
3. Count `ny(i)` = the same for the horizontal strip.

Average over all n points:

```
I = psi(k) + psi(N) − < psi(nx+1) + psi(ny+1) >          psi = digamma
```

Two details are load-bearing rather than cosmetic:

- **The Chebyshev metric is required.** It makes the joint neighbourhood a box
  whose projection onto each axis is exactly the strip counted in steps 2–3.
  A Euclidean ball would not have this property and the bias cancellation would
  fail.
- **Digamma, not log.** `psi` is the exact finite-count correction;
  `psi(m) ≈ ln m − 1/2m`. At the counts involved here (n = 18–48) the difference
  is not negligible.

### 2.3 Why it works

```
   INDEPENDENT                        DEPENDENT (y ≈ f(x))
   y                                  y
   │  ·   ·  ·   ·                    │           ··
   │ ·  ┌───┐ ·                       │        ·┌┐
   │  · │ · │  ·  ·                   │      ··│·│
   │ ·  └───┘   ·                     │    ··  └┘
   │   ·  ·  ·   ·                    │  ··
   └──────────────── x                └──────────────── x

     the box is small, but the         the box is small AND the x-strip
     x-strip still catches ~√(kN)      catches only ~k points — everything
     points from all over the range    else is far away along the curve

     psi terms large   →  I ≈ 0        psi terms small  →  I large
```

`eps` is set by how tightly the points cluster **jointly**; `nx` and `ny`
measure how many points that same width captures **marginally**. Under
independence the two agree and the formula collapses to zero. Under dependence
the points hug a curve, the strips are nearly empty relative to expectation, and
that shortfall *is* the mutual information.

The curve does not have to be a line. `y = x²`, `y = |x|` and `y = sin 3x`
concentrate points into thin strips just as effectively as `y = x`. This is the
entire motivation: Pearson scores `sin 3x` at |r| = 0.03, KSG recovers it
(test group 3 in `99_validate.py`).

### 2.4 How it is applied here

Every gene is first replaced by its **ranks across libraries**, ties broken at
random. Mutual information is invariant under monotone transformations of the
marginals, so this discards no information, and it buys three things:

- each gene becomes an exact permutation of 0…n−1, so under independence a gene
  *pair* is a random permutation pair. The null therefore depends only on (n, k)
  and **not on the genes** — it can be built once and reused for all ~1.46e10
  pairs. Per-pair permutation testing at this scale is otherwise impossible.
- rank differences are small integers, so fp16 arithmetic is bit-exact (verified
  against fp32 in `99_validate.py`), which halves both memory and runtime.
- no tied values, which KSG handles badly.

`k = 3` is the standard bias/variance compromise. Smaller k is less biased but
noisier; larger k smooths the neighbourhood until it is no longer local. At
n = 48, k = 3 means the ruler spans about 6% of the data.

Realised cost for the full sweeps (RTX A4500, fp16): **31 min** for sugarcane
(170,790 genes, 1.46e10 pairs) and **65 min** for purple (170,736 genes).

Output is in **nats**. For a bivariate Gaussian the identity
`I = −½·ln(1 − r²)` holds exactly; this is the anchor used later to place MI
edges on a correlation-equivalent scale.

### 2.5 Two properties that shape everything downstream

**The estimator can return negative values, and its null is not centred on
zero.** Measured on shuffled rank data:

| n | null mean | null sd |
|---|---|---|
| 48 | **+0.124** nats | 0.151 |
| 18 | **+0.180** nats | 0.146 |

This upward small-sample bias is why the pipeline never assumes `I = 0` under
independence and instead builds an explicit permutation null. (It is also why
`ksg − gcmi` is not a clean "non-linearity score" until the ~0.113 nat offset
between the two estimators' nulls at n = 48 is subtracted.)

**The estimator is noisy.** At n = 48 its standard deviation is ~0.15 nats
against a true signal of 0.51 nats for r = 0.8. That ratio drives everything in
section 3.

---

## 3. What `--match-pearson 0.8` does

The flag does exactly one thing: **it finds the MI value that admits false edges
at the same rate that `|r| >= 0.8` does.**

```
|r| = 0.8 at n=48  ──t-test──▶  p = 9.02e-12
                                    │
                                    │  "under independence, 9 in a trillion
                                    │   random gene pairs reach |r| >= 0.8"
                                    ▼
          find the MI value with that same tail probability
          under MI's own permutation null
                                    ▼
                            MI = 0.9891 nats
```

The MI null has no closed form, so it is built empirically: shuffle two rank
vectors, compute KSG, repeat (auto-sized, capped at 2e8). A target of p = 9e-12
lies far beyond what any affordable number of permutations resolves directly, so
a generalised-Pareto distribution is fitted to the extreme tail of the empirical
null (probability-weighted moments, Hosking & Wallis) and extrapolated to the
required quantile. Fitted parameters, the fit's upper endpoint, and the number of
edges falling beyond it are all recorded in `<out>.summary.json`.

### 3.1 The production thresholds

| study | n | p at \|r\|=0.8 | MI floor (nats) | r-equivalent | MI edges | Pearson edges | MI / Pearson |
|---|---|---|---|---|---|---|---|
| sugarcane | 48 | 9.02e-12 | **0.9891** | 0.907 | 2,871,301 | 75,380,961 | 3.8% |
| purple | 18 | 6.72e-05 | **0.9883** | 0.895 | 79,058,982 | 681,090,855 | 11.6% |

The r-equivalent column converts the floor through the Gaussian identity after
subtracting the null bias: `r_eq = sqrt(1 − exp(−2·(MI − bias)))`. Both land near
0.90, comfortably inside the pipeline's existing `[0.8, 0.9999]` weighting
window, so merging MI edges into `network_*_filtered_edges.tsv` does not shift
the weights of the existing Pearson edges.

### 3.2 The sense in which the two networks are comparable

**Equal specificity.** Both cuts admit noise at the same per-edge rate. Across
1.46e10 tested pairs, a Pearson network at `|r| >= 0.8` and an MI network at its
matched floor contain approximately the same (very small) number of pure-noise
edges.

This is the property that makes it legitimate to **take the union** of the two
edge sets: adding the MI layer does not loosen the error rate of the existing
Pearson network. Any other matching rule would.

An incidental but useful consequence: the two studies' floors land at **0.9891
and 0.9883 nats — 0.08% apart** — even though the p-value targets they were
derived from differ by seven orders of magnitude. **The MI layer is therefore
far more symmetric between the two studies than the Pearson layer is**, where
the same nominal `|r| = 0.8` means p = 9e-12 in sugarcane and p = 6.7e-5 in
purple. For a conservation comparison that is a real point in its favour, and it
is worth stating explicitly in the methods.

### 3.3 What is NOT matched

Everything else — and this must be stated plainly in any methods section.

Viewed from the *signal* side rather than the noise side, at n = 48:

| threshold matched at p = 9e-12 | Pearson \|r\| >= 0.8 | KSG >= 0.9891 nats |
|---|---|---|
| distance from the null | 7.37 sd | ~5.7 sd |
| distance from the centre of a **true r = 0.8** pair's distribution | **0.00 sd** | **~+1.9 sd** |
| power at true r = 0.8 | **50%** | **~3%** |

The middle row is the crux. Pearson's cut at 0.8 sits exactly *on top of* the
distribution of genuine r = 0.8 pairs — by construction, half of them land above
it. The matched MI floor sits nearly two standard deviations *past* the centre of
the corresponding distribution, so almost none of them clear it.

The cause is estimator noise, not a flaw in the matching. At r = 0.8 the KSG
estimate separates from its null by ~3.5 sd, where Fisher-z separates by 7.4 sd —
roughly half the discriminating power per sample. Buying equal specificity from
a noisier estimator is paid for in sensitivity, and at n = 48 the bill is ~94% of
the power.

At n = 18 power at r = 0.8 is ~9%, better only because the Pearson p-value there
is a far softer 6.7e-05.

The practical consequence is that **the MI network is much smaller than the
Pearson network and does not contain most of it.** Measured on the complete
networks by the merge step (section 5.1):

| | sugarcane (n=48) | | purple (n=18) | |
|---|---|---|---|---|
| Pearson only | 73,375,808 | 96.23% | 630,994,886 | 88.87% |
| **both** | **2,005,153** | 2.63% | **50,095,969** | 7.06% |
| **MI only** | **866,148** | 1.14% | **28,963,013** | 4.08% |
| union | 76,247,109 | +1.15% | 710,053,868 | +4.25% |

Read two ways:

- **30.2% of sugarcane MI edges and 36.6% of purple MI edges are invisible to
  Pearson.** That is the return on the method — those edges are exactly what a
  correlation cannot represent. Both are better than the 14.2% measured earlier
  at the looser 0.7-matched floor: raising the threshold discards linear edges
  faster than non-linear ones, so the surviving MI set is proportionally more
  novel.
- **97.3% (sugarcane) and 92.6% (purple) of Pearson edges do not clear the MI
  floor.** That is the cost, and it is fully explained by the power argument
  above. It is not evidence that those edges are false.

**An independent check on the power model.** The fraction of Pearson edges that
the MI layer recovers should be predicted by the simulated power at the
threshold. Predicted power at r = 0.8 was 3.1% (n = 48) and 8.5% (n = 18);
observed recovery was **2.66%** and **7.36%**. Right magnitude, right ordering
between the two studies, both slightly below prediction — which is what one
expects if real co-expression pairs are somewhat less clean than the bivariate
Gaussians the power curve simulates (outliers and heteroscedasticity both pull
the KSG estimate down). The power argument in this section is therefore not
merely theoretical; it predicts the observed behaviour of the real data to
within about half a percentage point.

### 3.4 There is no setting that removes this trade-off

Lowering the MI floor to match Pearson's *effect size* instead of its
false-positive rate — i.e. to `I = 0.511` nats, the true MI of an r = 0.8
Gaussian pair — puts the threshold **below the null's 99.9th percentile** at
these sample sizes. The resulting network would be predominantly noise. The
trade-off is a property of n = 18–48, not of the implementation.

---

## 4. How the result should be interpreted

**Correct framing.** The Pearson network, plus a discovery layer of non-monotone
edges detected at equal specificity, tagged by provenance. Each edge in the
merged network carries a `source` field (`pearson`, `mi`, `both`) so every
downstream result can be split by how the edge was found.

**Incorrect framing.** "An MI network, which is better than the Pearson network."
It is not a replacement: at these sample sizes it recovers a minority of what
Pearson already finds. Reporting it as a standalone network invites the reader to
interpret its smaller size as less co-expression, which would be wrong.

Two specific cautions:

- **Do not compare MI edge counts between the two studies.** n = 48 versus n = 18
  changes both the null and the power, so the fact that purple retains 79 M MI
  edges and sugarcane 2.9 M is a statistical artefact of sample size, not
  biology.
- **Do not read a smaller MI network as evidence of weaker co-expression.** It is
  evidence of a stricter threshold on a noisier statistic.

If a *density-matched* network is wanted instead — same edge count as the Pearson
network, threshold chosen so that network structure is comparable rather than
error rate — that is `--min-value` set by hand, and the methods must say the
threshold was chosen for equal density, not equal significance.

---

## 5. Applying the 0.8 threshold to the completed sweeps

Both sweeps were run before this threshold question was settled and therefore
used the 0.7-matched floor (0.8376 nats for sugarcane, 0.7744 for purple). **They
do not need re-running.** The floor is applied only at assembly time, after the
per-edge p-values are computed, and both edge lists were written well below the
0.8-matched floor. Tightening the cut is a pure filter on the `pval` column that
the edge lists already carry:

```bash
# sugarcane: |r| = 0.8 at n = 48  ->  p <= 9.02e-12
awk -F'\t' 'NR==1 || $4+0 <= 9.02e-12' \
    files/sugarcane/mi/sugarcane_ksg.edgelist.tsv \
  > files/sugarcane/mi/sugarcane_ksg.p08.edgelist.tsv

# purple: |r| = 0.8 at n = 18  ->  p <= 6.72e-05
awk -F'\t' 'NR==1 || $4+0 <= 6.72e-05' \
    files/purple/new/mi/purple_ksg.edgelist.tsv \
  > files/purple/new/mi/purple_ksg.p08.edgelist.tsv
```

The `padj` values are unaffected: BH was computed over all 1.46e10 tests from
the candidate set, independently of where the effect-size floor was placed, so
raising the floor only removes rows.

Verification that this reproduces the intended thresholds — the minimum retained
KSG value after filtering was **0.98910** (sugarcane) and **0.98827** (purple),
matching the generalised-Pareto quantiles to five and four decimal places
respectively. Both filters have been run; the results are
`*_ksg.p08.edgelist.tsv` (2,871,301 and 79,058,982 edges).

## 5.1 Merging the two networks

`scripts/mutual_information/10_augment_network.py` takes the union and writes it
in exactly the format `general_stats.r` produces, so **every downstream script
consumes it unchanged** — `mcl_clustering.r`, `eigengene.r`, `module_trait_cor.r`,
`comparative_networks2.r`, `full_matrix_edge_comparation.r` and
`network_conservation_join.r` all read only `gene1`/`gene2`/`weight`, or just the
gene pair.

```
gene1  gene2  stat  pval  padj  weight  pearson_r  ksg  source
```

- `stat` — effective strength on the **correlation** scale: `|r|` for Pearson
  edges, `r_eq` for MI-only edges, the larger of the two for edges both found.
  `general_stats.r` reads column 3 positionally, so re-running it on this file
  works.
- `weight` — `0.01 + (stat−0.8)/(0.9999−0.8)·0.99`, i.e. `general_stats.r`'s
  minmax with the bounds fixed at `PEARSON_MIN`/`PEARSON_MAX`. The script
  **verifies against the input's own weight column** and reports the maximum
  discrepancy; on sugarcane it was **6.66e-16**, so existing Pearson edge weights
  are preserved bit-for-bit and no re-weighting artefact is introduced.
- `source` — `pearson` | `mi` | `both`. This is the point of the whole exercise:
  every downstream result can be split by how the edge was found.
- `pearson_r` is `NA` for MI-only edges and `ksg` is `NA` for Pearson-only edges,
  so neither column can be silently misread as zero.

**`r_eq` is not a correlation.** For a non-monotone edge it means "as strong as a
linear r of this size", not "its Pearson r is this". Never interpret it without
`source`.

Both have been run:

```bash
cd scripts/mutual_information
python 10_augment_network.py --study sugarcane   #  6.5 min,   8.3 GB   [DONE]
python 10_augment_network.py --study purple      # 73.9 min,  69.8 GB   [DONE]
```

| output | edges | size | weight self-check |
|---|---|---|---|
| `files/sugarcane/network_sugarcane_augmented_edges.tsv` | 76,247,109 | 8.3 GB | 6.66e-16 |
| `files/purple/new/network_purple_augmented_edges.tsv` | 710,053,868 | 69.8 GB | 5.55e-16 |

One naming detail the script handles and asserts on: the sugarcane Pearson
network carries bare locus ids (`SoffiXsponR570.01Ag000100`) while the MI gene
list carries the annotation version suffix (`…v2.1`). The suffix is stripped and
the result checked for collisions — all 170,790 ids stay unique. Purple ids are
unaffected. If the two gene sets ever genuinely diverge the script aborts rather
than silently dropping edges.

---

## 6. Reproducing the numbers in this note

```bash
cd scripts/mutual_information

# correctness suite: GPU kernel vs naive CPU KSG, analytic Gaussian MI,
# non-monotone detection, rectangular tiles, fp16 vs fp32, null + GPD fit
./run_all.sh validate

# nulls, matched floors, power curves, both n
/home/genomics/miniconda3/envs/docling/bin/python 98_power_curve.py \
    --match-pearson 0.8 --n 48 18
```

Fresh sweeps at the corrected threshold, should they ever be wanted, are
`./run_all.sh sugarcane` and `./run_all.sh purple` — `config.sh` now defaults to
`MATCH_PEARSON=0.8`, and both always run with `--resume`.

## 7. References

- Kraskov A., Stögbauer H., Grassberger P. (2004) Estimating mutual information.
  *Phys. Rev. E* **69**, 066138.
- Ince R.A.A. *et al.* (2017) A statistical framework for neuroimaging data
  analysis based on mutual information estimated via a Gaussian copula.
  *Hum. Brain Mapp.* **38**, 1541–1573. *(the `gcmi` estimator)*
- Chatterjee S. (2021) A new coefficient of correlation. *JASA* **116**,
  2009–2022. *(the `xi` estimator)*
- Hosking J.R.M., Wallis J.R. (1987) Parameter and quantile estimation for the
  generalized Pareto distribution. *Technometrics* **29**, 339–349.
  *(the tail extrapolation)*
- Benjamini Y., Hochberg Y. (1995) Controlling the false discovery rate.
  *JRSS B* **57**, 289–300.
