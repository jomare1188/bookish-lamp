# Comparative gene co-expression networks in sugarcane, from reused public RNA-seq

Two independently published sugarcane nitrogen experiments, re-quantified through
one pipeline and turned into two comparable gene co-expression networks. The
question the project asks is whether the **nitrogen response is conserved between
them** — at the level of individual genes, of network edges, and of biological
processes.

The secondary goal, and the reason the data are public ones, is to show what can
be recovered from reusing existing sugarcane RNA-seq rather than generating new
libraries.

Networks are built from **two layers**: Pearson correlation (linear) and
Kraskov–Stögbauer–Grassberger mutual information (non-linear), thresholded at the
same per-edge false-positive rate so their union is legitimate, with every edge
labelled by which layer found it. That label is what makes the central
methodological question answerable — *does non-linear co-expression behave
differently from linear co-expression?*

---

## Start here

| you want | go to |
|---|---|
| to run the pipeline | [`new_clean/README.md`](new_clean/README.md) |
| the results, with caveats | [`new_clean/docs/results.md`](new_clean/docs/results.md) |
| why it is built this way | [`new_clean/docs/decisions.md`](new_clean/docs/decisions.md) |
| what each stage does | [`new_clean/docs/methods.md`](new_clean/docs/methods.md) |
| the MI threshold argument | [`new_clean/docs/thresholds.md`](new_clean/docs/thresholds.md) |

**`new_clean/` is the pipeline.** Everything at the top level of `scripts/` is the
earlier version, kept for history — see [Repository layout](#repository-layout).

---

## The data

| Network | Study | Genotype | Trait | Samples | Organism | Reference |
|---|---|---|---|---|---|---|
| sugarcane | Muñoz-Perez et al. 2025 | RB937570 | non-responsive | 24 | *Saccharum* hybrid | R570 |
| sugarcane | Muñoz-Perez et al. 2025 | RB975375 | responsive | 24 | *Saccharum* hybrid | R570 |
| purple | Ta Quang Kiet et al. 2025 | 51NG3 | responsive | 9 | *S. robustum* | LA purple |
| purple | Ta Quang Kiet et al. 2025 | TAGZ | non-responsive | 9 | *S. officinarum* | LA purple |

- Muñoz-Perez et al. 2025 — <https://onlinelibrary.wiley.com/doi/10.1111/ppl.70612>
  — two genotypes of contrasting nitrogen-use efficiency, high vs low nitrogen,
  sampled across leaf segments. Counts from Sweet Recycler precomputed matrices.
- Ta Quang Kiet et al. 2025 — <https://www.sciencedirect.com/science/article/pii/S0926669025017017>
  — *S. officinarum* and *S. robustum* across a nitrogen **gradient** (0/2/6 mM).
  Downloaded from the SugarCane multi-Omics Database and quantified with the same
  nf-core/rnaseq + salmon workflow.

**The sample-size difference is the single most important fact about this
project.** n = 48 against n = 18 sets what each study can support, and nearly
every underpowered result below traces back to it.

## The pipeline

```
dds ──► VST export ──┬──► GPU Pearson  (analytic t null)  ──┐
                     └──► GPU KSG      (permutation null)  ──┴──► merge ──► network
                                                                              │
     topology ──► MCL modules ──► cross-species conservation ──► trait ──► GO ──► H1 readouts
```

All-pairs dependence runs on the GPU straight to a thresholded, FDR-corrected
edge list — the n × n matrix is never materialised. Both studies, VST to merged
edge tables, take about **5.5 hours**; everything through clustering about **9**.

## Results

### Networks

| | sugarcane (n=48) | purple (n=18) |
|---|---|---|
| genes after CV ≥ 15 filter | 170,790 | 170,740 |
| **edges** | **76,200,344** | **705,571,723** |
| nodes | 103,336 | 170,736 |
| components / giant | 939 / 101,253 | 1 / 170,736 |
| pearson only | 73,329,192 (96.2%) | 625,775,930 (88.7%) |
| both layers | 2,004,577 (2.6%) | 50,179,988 (7.1%) |
| **MI only** | **866,575 (1.1%)** | **29,615,805 (4.2%)** |
| MCL modules | 10,309 | 9,881 |
| modularity Q | 0.1005 | 0.1631 |

Purple's larger MI share is **not** more non-linear biology: at n = 18 the
p-value implied by |r| = 0.8 is 6.7e-05 against 9.0e-12 at n = 48, so its matched
floor is far softer. Purple is also a single connected component with mean degree
8,265 — dense enough that its module structure should be read cautiously.

### Cross-species edge conservation

Both directions run **~2.5× above** a permutation null that preserves orthology
fan-out and coverage, so conservation is real:

| direction | observed | null | fold | 
|---|---|---|---|
| sugarcane → purple | 10.76% | 4.18% | **2.57** |
| purple → sugarcane | 1.52% | 0.61% | **2.46** |

The raw rates differ 7× only because purple has 9.3× more edges. **The fold over
null is the comparable quantity** and the two directions agree on it.

**Edges found by both estimators are better conserved** than Pearson-only edges
(1.08× and 1.20× over null). **MI-only edges are not** — 1.005 and 0.921 once
opportunity is accounted for. An apparent MI advantage in the raw rates turned
out to be opportunity bias: MI-only edges connect genes with more orthologs, so
they had more chances to match by accident.

### Nitrogen response

| | sugarcane | purple |
|---|---|---|
| trait-correlated genes (of conserved set) | 1,361 | **62** |
| \|r\| needed to clear FDR | 0.378 (the 0.6 cut binds) | **0.799** (FDR binds) |
| **conserved correlated ortholog pairs** | **1** (≈4.4 expected by chance) | |

**No evidence of a shared node-level nitrogen response.** The bottleneck is
purple's 62 genes, and it is a power problem: at n = 18 over 44,118 genes a gene
needs |r| ≈ 0.80 just to clear the FDR. A mutual-information gene–trait test was
built to attack exactly this and **did not fix it** — 21 significant genes against
Pearson's 85. n = 18 is the binding constraint, not the choice of statistic.

### Function and H1 readouts

- **GO** (topGO weight01, raw p ≤ 0.05, background = each network's own nodes):
  **69 shared BP terms** (106 sugarcane, 155 purple), 39 MF, 24 CC. Shared
  processes despite ~no shared responsive genes — an ordinary evolutionary
  pattern, but note the power gap between the two tests before reading it as
  agreement.
- **Transcription factors**: 7,183 (sugarcane) and 12,197 (purple) TF genes in
  the networks, 68 and 67 families.
- **MYB61**: purple's 15 copies sit on cross-species conserved edges far more than
  expected — 60% vs 25.8% background, **degree-matched permutation p = 0.0010**,
  and *not* a hub artifact (their median degree is below average). Sugarcane's 8
  copies show nothing (p = 0.67).
- **Muñoz Module 20**: nothing significant in either species. The module's
  nitrogen response transfers between studies; its network position does not.

## Repository layout

```
new_clean/            THE PIPELINE — scripts, config, docs, and results
  run.sh              one command per stage; nothing is edited between runs
  config.sh           every path and parameter
  scripts/            01 export ... 13 conservation null, + H1 readouts
  docs/               decisions, methods, thresholds, results
  results/            all output (gitignored — regenerable)

scripts/              the EARLIER pipeline, kept for history. Superseded by
                      new_clean/; see new_clean/docs/decisions.md for what
                      changed and why.
  old/                superseded even within that generation

files/                correlation matrices, edge lists, networks (gitignored,
                      ~1.2 TB). files/pearson_baseline/ holds the pre-augmentation
                      results for comparison.
GET_TFS/              TF identification: HMM databases and cached hmmsearch output
annotation/           eggNOG-mapper annotations for both references
china/, run1/         nf-core/rnaseq output for the two studies (gitignored)
```

## Reproducing

```bash
cd new_clean
./run.sh validate                # correctness suite — must pass first
./run.sh build sugarcane         # VST -> both layers -> network
./run.sh build purple
# then: stats, mcl, conserve, conservenull, trait, traitmi, conscor, go, gosem
```

`new_clean/README.md` has the full stage list with measured runtimes. Every stage
logs to `new_clean/logs/`, and the GPU sweeps checkpoint per tile so a killed run
resumes.

## Things a reader should know before using these results

1. **n = 18 limits purple everywhere.** Its MI network threshold is soft, its
   trait selection is FDR-bound rather than effect-size-bound, and its module
   structure is one dense component.
2. **The previous networks contained ~47,000 sugarcane edges below their own
   0.8 threshold**, admitted because the correlation matrix was stored as text
   rounded to 4 decimals. Fixed here; the affected results are archived under
   `files/pearson_baseline/` and were never corrected, only superseded.
3. **`r_eq` is not a correlation.** For a non-monotone edge it means "as strong
   as a linear r of this size". Never read it without the `source` column.
4. **Raw conservation rates are not comparable between directions.** Use the fold
   over null.
5. **The MI layer's value is not settled.** It adds ~0.9 M sugarcane edges the
   linear layer cannot see, but those edges are not better conserved. What is
   robust is that edges *both* estimators find are better conserved than either
   alone.
