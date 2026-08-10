# Mutual-information co-expression networks (GPU)

An alternative to `scripts/pearson_cor.r` + `scripts/build_edgelist.r` that measures
**non-linear** dependence. Same gene sets, same output format, so the existing
downstream scripts consume it unchanged.

```
00_export_vst.r              dds -> flat float32 matrix (+ gene list, meta)
custom_mutual_information.py the GPU engine: all-pairs sweep -> edgelist
98_power_curve.py            what dependence strength is detectable at this n
99_validate.py               correctness suite (run this after any edit)
run_all.sh                   driver;  config.sh  holds every path and parameter
```

```bash
./run_all.sh export      # both VST matrices (verifies gene sets vs Pearson)
./run_all.sh validate    # ~2 min, must print ALL CHECKS PASSED
./run_all.sh sugarcane   # ~2.5 h on the A4500
./run_all.sh purple      # ~45 min
```

Sweeps always run with `--resume`; re-running after a crash picks up at the first
unfinished tile.

---

## Read this before using the output

**At these sample sizes, MI is not a drop-in "better Pearson". It is a stricter
network that trades most of Pearson's edges for the ability to see non-monotone
ones.** Two measurements, both from this pipeline:

On the first 20,000 sugarcane genes, with the MI floor matched to Pearson `|r| = 0.7`:

| | |
|---|---|
| KSG edges also visible to Pearson (`\|r\| >= 0.7`) | 265,210 (85.8%) |
| **KSG edges invisible to Pearson** (`\|r\| < 0.7`) | **43,819 (14.2%)** |
| …of those, strongly non-monotone (`\|r\| < 0.3`) | 1,062 |
| **Pearson `\|r\| >= 0.7` pairs absent from the MI network** | **2,510,478 (90.4%)** |

The 14.2% is the prize: those edges are exactly what Pearson cannot represent.
The 90.4% is the price, and `98_power_curve.py` explains it:

```
n = 48   null mean +0.124 (KSG is biased upward at small n)
         |r| = 0.7  ->  p = 3.08e-08  ->  MI floor 0.832 nats

  |r|   true MI   KSG mean   KSG sd    power
 0.70    0.3367     0.4889   0.1506     1.6%
 0.80    0.5108     0.6699   0.1609    15.6%
 0.90    0.8304     0.9913   0.1696    82.6%
```

`--match-pearson` equalises the **false-positive rate**, not the effect size.
Because KSG at n = 48 has a standard deviation of ~0.15 nats against a signal of
~0.34, holding the per-edge false-positive rate at Pearson's level forces a floor
that behaves like a Pearson cut at **|r| ≈ 0.9**. Power at the nominal |r| = 0.7 is
1.6%. At n = 18 the null is wider still (sd 0.146) and power at |r| = 0.7 is 11.5%.

**The production setting is `--match-pearson 0.8`, not 0.7.** `build_edgelist.r`'s
`PEARSON_THRESHOLD <- 0.7` is only a pre-filter; the analysed networks come from
`general_stats.r`, whose `PEARSON_MIN <- 0.8` is the real cut — verified against
the data, both `network_*_filtered_edges.tsv` bottom out at |r| = 0.8. At 0.8 the
floors are 0.9891 nats (n = 48, 2,871,301 edges) and 0.9883 nats (n = 18,
79,058,982 edges). See `MI_for_review.md` at the repository root.

There is no setting that fixes this. Lowering the floor to match Pearson's *effect
size* (MI = 0.337 nats for r = 0.7) puts it **below the null's 99.9th percentile**,
so the network would be mostly noise. The honest summary is:

- **Use MI to find what Pearson misses**, i.e. treat the ~14% of MI edges with low
  `|r|` as the result, and interpret them directly.
- **Do not** read a smaller MI network as evidence of less co-expression, and do
  not compare MI edge counts between the two studies — n = 48 vs n = 18 changes
  both the null and the power.
- If you want a density-matched network for a structural comparison, set
  `--min-value` explicitly rather than `--match-pearson`, and say in the methods
  that the threshold was chosen for equal density, not equal significance.

---

## Estimators (`--estimator`)

| | what it sees | cost | use for |
|---|---|---|---|
| `ksg` | any dependence, monotone or not | ~1.6 M pairs/s at n=48 | the real run |
| `gcmi` | monotone only (rank correlation in disguise) | one matmul | fast baseline; `ksg >> gcmi` isolates the non-linear pairs |
| `xi` | any dependence, incl. non-monotone | cheap | Chatterjee's ξ; noisy at small n — on the synthetic test it produced 8 false positives where KSG produced 1 |

`gcmi` is a proven lower bound on the true MI, so **`ksg` minus `gcmi` is a
usable "how non-linear is this edge" score**. Running both and taking the
difference is the cheapest way to rank candidate non-linear relationships.

## How it works

Every gene is replaced by its **ranks** across libraries, ties broken at random.
This costs no information (MI is invariant under monotone marginal transforms),
removes the VST's residual mean–variance trend, and eliminates the KSG estimator's
degenerate tie case. Most importantly it makes the **null data-independent**: after
ranking, an independent pair *is* a random permutation pair, so the null depends
only on (n, k). The permutation null is therefore built once and reused for all
1.5e10 pairs — per-pair permutation testing at this scale is otherwise impossible.

`--n-perm 0` (the default) sizes that null automatically to resolve the strictest
threshold empirically, capped by `--max-perm`. Below the empirical range a
generalised-Pareto tail (PWM fit) extrapolates; edges past the fitted endpoint get
`pval` at the floor, which means "beyond the fit", not a calibrated number. The
count is reported and stored in `summary.json`.

**GPU**: the KSG kernel evaluates an `(a, b, n, n)` tensor of joint Chebyshev
distances in fp16 — exact, since every value is a small integer rank difference,
and verified bit-identical to fp32 by the test suite. `--gpu-mem-gb` sets the tile
size; measured throughput saturates past ~3 GB, so the default 5 GB is ample on a
20 GB card.

**RAM**: the expression matrix is 33 MB and stays on the GPU. The host never holds
anything of order `n_genes²`. Candidate edges stream to 12-byte binary shards that
double as checkpoints. Peak host RAM is 8 bytes per candidate edge during the FDR
step (purple is the demanding one, ~15 GB).

## Thresholds

Two independent cuts, both applied:

- `--cand-p` (default 1e-3) — what gets written to a shard, and the range over
  which BH stays exact (up to `cand_p * n_tests / alpha` rejections). Raised
  automatically if it would not clear the effect-size floor by 20×, which is what
  happens at n = 18 where `|r| = 0.7` is only p = 1.2e-3.
- `--match-pearson` (default 0.8, matching `general_stats.r:PEARSON_MIN`) — the
  effect-size floor, as described above. `--min-value` overrides it with a raw
  statistic value. Because the floor is applied only at assembly time, a sweep
  run at a looser setting can be tightened afterwards by filtering its `pval`
  column; it does not have to be repeated.

BH runs over **all** `n_genes*(n_genes-1)/2` tests, not just the candidates; this is
exact because every unstored pair has a larger p than every stored one. If the
final cut ever reaches the weakest stored candidate the run says `TRUNCATED` —
that is the one condition that invalidates the edge list.

## Provenance note: which `dds` feeds which network

The purple Pearson network was built from **`china/run1`** (36 libraries, leaf and
root) narrowed to the 18 leaf libraries by `Group1 == "L"` — *not* from
`china/run2_onlyL`, which the `module20` and `myb61` expression scripts use. The
two quantifications disagree about which genes are flat, and their CV≥15 gene sets
differ by 68 genes. `00_export_vst.r` reproduces the Pearson selection exactly and
**aborts** if the gene set ever stops matching the Pearson matrix header, so this
cannot drift silently.

## Environment

`config.sh:PYTORCH` points at the `docling` conda env — the only one on this
machine with a CUDA build (torch 2.11 + cu130, Python 3.14). It is borrowed, not
purpose-built; a dedicated env is:

```bash
conda create -n mi_gpu python=3.12 numpy scipy
conda run -n mi_gpu pip install torch --index-url https://download.pytorch.org/whl/cu124
```

then point `PYTORCH` at it. Nothing else in the pipeline depends on the env.

## Validation

`99_validate.py` checks, in order: the GPU kernel against a deliberately naive
per-pair CPU KSG; recovery of the analytic MI of a bivariate Gaussian; detection of
`y = x²`, `y = |x|`, `y = sin 3x` where Pearson collapses; rectangular tiles
(the last tile of every sweep — a real bug lived here); the paired kernel used for
the null against the diagonal of the block kernel; null calibration and the GPD
tail; and fp16 vs fp32. Run it after any change to the kernels.
