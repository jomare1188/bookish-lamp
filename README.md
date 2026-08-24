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
| the paper figures and their legends | [`figures_legends.txt`](figures_legends.txt) |

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
     topology ──► MCL modules ──┬─ cross-species conservation ─► trait ─► GO ─► H1 readouts
                                │
                                └─ eigengenes ─► module response ─► TF ─► per-module GO
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

"Trait-responsive" here means: a gene sitting on **at least one conserved edge**
(any layer) that is also associated with nitrogen. Association is tested by
Pearson **and** by mutual information, and the selection rule is configurable
(`TRAIT_SELECTION`, default `union`) — because a Pearson-only search returning
nothing cannot distinguish "no shared response" from "no *linear* shared
response".

| | sugarcane | purple |
|---|---|---|
| genes on conserved edges | 39,226 | 44,118 |
| responsive — Pearson | 1,361 | **30** |
| responsive — MI | 3,221 | **5** |
| responsive — union | 3,265 | **32** |
| \|r\| needed to clear FDR | 0.378 (the 0.6 cut binds) | **0.799** (FDR binds) |

**Node level — conserved correlated ortholog pairs: 1 (Pearson), 0 (MI), 2 (union)**,
against ≈4.4 expected by chance.

**Edge level — 0, under every selection rule.** Conserved edges joining two genes
responsive in *sugarcane* are plentiful, and the MI layer contributes a real share
of them:

| selection | conserved edges, both endpoints responsive in sugarcane | of which `pearson` / `both` / **`mi`** |
|---|---|---|
| Pearson | 5,894 | 4,761 / 1,065 / **68** |
| union | 35,761 | 26,170 / 6,993 / **2,598** |

But **not one of them has a purple counterpart whose endpoints are also
responsive**, because an edge needs *two* genes responsive on both sides and there
are only 1 (Pearson) or 2 (union) such genes in the whole analysis — and they are
not connected.

So the funnel does not close at the edges, or at the MI layer. It closes at
**purple's 30–32 responsive genes**, which is a power result: at n = 18 over
44,118 genes a gene needs |r| ≈ 0.80 merely to clear the FDR. A mutual-information
gene–trait test was built specifically to attack this and **did not fix it** — it
*shrinks* purple's set (30 → 5 alone).

### The two species' responsive genes are independent in ortholog space

Testing purple genome-wide is the wrong burden for a comparative question, so the
test was repeated **directed**: correct purple's p-values over only the 4,745
orthologs of sugarcane-responsive genes rather than over all 44,118
(`./run.sh conscor 1`). That drops the |r| a gene must reach from ≈0.88 to ≈0.84.

It found **fewer**, not more — 1 gene against 2 under the genome-wide burden on
the identical candidate set. The reason is the finding:

> Of purple's **32** genome-wide nitrogen-responsive genes, only **2** are
> orthologs of a sugarcane-responsive gene. Of purple's top 30 genes by p-value,
> 2 are candidates against ~3.2 expected if the two sets were independent.

**The nitrogen-responsive gene sets of the two species are, in ortholog space,
independent — marginally below chance overlap.** The candidate set is therefore
*depleted* of purple's strongest signal, and BH is adaptive: shrinking the
denominator does not compensate for losing that company.

And the signal is genuinely absent rather than hidden by a threshold. Among all
4,745 candidates the best |r| is **0.854** against a required 0.844, and only
**3** reach |r| ≥ 0.8 at all. No choice of denominator turns three genes into a
conserved responsive edge, which needs two connected genes responsive on both
sides. (At raw p ≤ 0.05 and |r| ≥ 0.6 — no multiple-testing correction — 128
candidates qualify; that is the most generous reading available and should be
labelled uncorrected.)

**The zero survives three selection rules, two correction burdens and both
conservation directions.** n = 18 with 0/2/6 mM in triplicate is a design limit,
not a method limit.

### Module-level nitrogen response

The gene level dies on a testing burden of 39,226 / 44,118 genes. One eigengene
per MCL module cuts that to ~6,500 tests, and a module signal can survive where a
single gene cannot. Each eigengene is tested against nitrogen by **Spearman's rho
alone** — `padj <= 0.05` and `|rho| >= 0.6`, the same two thresholds the gene
level uses.

**Not Pearson, and not MI.** The trait is ordinal — purple's nitrogen levels are
0/2/6 mM and Pearson reads that spacing as arithmetic the design never claimed.
MI at this resolution was an omnibus test firing on two-library quirks, and its
effect-size floor could only be calibrated against the linear one, never derived.
On the identical eigengenes at identical thresholds:

| | Pearson | **Spearman** | Spearman only | Pearson only |
|---|---|---|---|---|
| sugarcane (n=48) | 408 | **465** | 71 | 14 |
| **purple (n=18)** | **38** | **79** | **45** | 4 |

**Purple's responsive set more than doubles** — and purple is the study with the
gradient, which is exactly where the ordinal argument predicts the gain.

| | sugarcane | purple |
|---|---|---|
| modules tested | 6,576 | 6,318 |
| **responsive** | **465** (242 up, 223 down) | **79** (32 up, 47 down) |
| which threshold binds | \|rho\| ≥ 0.6 (BH never binds at n=48) | padj (274 clear \|rho\| ≥ 0.6) |
| median \|rho\| | 0.695 | 0.787 |
| label-permutation null | 0 of 1,000 shuffles reach 465 | 2 of 1,000 reach 79 (p = 0.002) |
| TF-enriched | **3.23% vs 1.24%, OR 2.65, p = 0.0016** | 0 of 79 |
| GO-testable (≥3 annotated genes) | 49 (11%) | 10 (13%) |

Two things came back that the three-statistic call had buried. **TF enrichment**
was a null spread across classes too small to separate (`both` at OR 1.96,
p = 0.061); pooled into one responsive set it is OR 2.65, p = 0.0016, and it sits
in the modules that go *down* with nitrogen. And the responsive set is now
*cleaner* than the background it is drawn from — 20% of responsive modules have
two of 48 samples carrying over a quarter of the eigengene variance, against 50%
of non-responsive ones. The old `mi_only` class ran the other way, at 35%.

Sugarcane's testable modules find nitrogen assimilation repeatedly and
independently — nitrate assimilation (Module_026, Module_100), the ammonia
assimilation cycle and glutamate biosynthesis (Module_440), ammonium metabolism
and the polyamines (Module_469), proline and asparagine biosynthesis — **all
rising with nitrogen**, and MF independently names the enzymes BP inferred from
the processes (nitrate reductase, glutamate synthase). Purple's 10 testable
modules give phenylpropanoid biosynthesis, cellulose biosynthesis and sulfate
assimilation; that is a report of ten gene sets, not a characterisation of its
nitrogen response.

**Purple's 79 modules are not 79 independent findings.** Its network is one dense
component, so its eigengenes are strongly correlated and one lucky label shuffle
in a thousand produced 244 "responsive" modules against 79 observed. The signal
is real at p = 0.002; no single purple module should be believed on its padj alone.

### What the responsive modules do — and the two directions do different things

One topGO enrichment per responsive module names individual modules, but the
annotation gate is brutal: **49 of 465** sugarcane modules carry enough
GO-annotated members to be testable at all (10 of 79 in purple), biased toward
the large ones.

So the modules that **rise** with nitrogen are also tested as one set against the
modules that **fall**. Pooling dissolves the gate — every module contributes —
and the two directions turn out to be **almost disjoint in function**:

| sugarcane BP terms | rises with N | falls with N | both |
|---|---|---|---|
| enriched at raw p ≤ 0.05 | **74** | **23** | 6 |

| rises with nitrogen | falls with nitrogen |
|---|---|
| **nitrate assimilation**, nitric oxide biosynthesis | flavonoid biosynthesis |
| **proline**, asparagine, spermidine biosynthesis | raffinose-family oligosaccharides |
| reactive oxygen species biosynthesis | triglyceride biosynthesis, cold response |
| defence: fungus, bacterium, chitin, monoterpenes | skotomorphogenesis, flowering time |

Nitrogen assimilation and the metabolism that consumes it go up; the low-nitrogen
carbon programme goes down. **Neither of the other two grains recovers this** —
the per-module grain sees 11% of the modules, and pooling all responsive modules
averages the two directions away.

### Module 20 in purple — one copy responds, and not monotonically

Muñoz's 12 published members are not 12 genes: they resolve to **three
Arabidopsis anchors** and the mapped genes are haplotype copies — 33 sugarcane
genes from 10 members, 27 purple from 6. Only the AtMYB59 anchor carries a MYB
call under our own domain rules (9 of 9 copies in sugarcane, 2 of 3 in purple);
the two anchors supplying most of the mapped genes are not TFs at all.

**One purple copy is expressed, and it responds in a U.**
`Soffic.09G0001580-9H` — MYB by both classifications, 60.6 TPM against 0.01 and
0.00 for the other two AtMYB59 copies:

| 51NG3 | 0 N (starvation) | 2 N (control) | 6 N (excess) |
|---|---|---|---|
| mean TPM | **124.1** | **10.8** | **47.4** |

An 11-fold drop to the control and a 4-fold rise beyond it (U-contrast
p = 0.0078). Since 2 N is the control and both extremes are stresses, it is
**induced by nitrogen stress in either direction** — a shape every monotone test
here is blind to. Sugarcane's nine copies at the same anchor are all
monotonically repressed by nitrogen, but sugarcane has no control level, so it
can only observe one arm of the same curve. The two studies are consistent, not
contradictory.

It rests on one expressed copy at a p that does **not** clear BH over the 27
purple Module-20 genes: a candidate, not a finding.

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
- **Muñoz Module 20**: no network-position signal in either species. The module's
  nitrogen response transfers between studies; its network position does not.
  See above for the one purple copy that does respond.

## Repository layout

```
new_clean/            THE PIPELINE — scripts, config, docs, and results
  run.sh              one command per stage; nothing is edited between runs
  config.sh           every path and parameter
  scripts/            01 export ... 13 conservation null, 14-18 module level,
                      + H1 readouts
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
# module level, per study (~5 min):
#   eigengene, moduletrait, moduleprofile, modulego, moduleheatmap, modulesummary
```

`new_clean/README.md` has the full stage list with measured runtimes. Every stage
logs to `new_clean/logs/`, and the GPU sweeps checkpoint per tile so a killed run
resumes.

## Things a reader should know before using these results

1. **Purple's three nitrogen levels are not a dose series.** 2 N is the
   **control**; 0 N and 6 N are stresses in opposite directions. Every monotone
   test in this project therefore asks whether something tracks nitrogen
   *supply*, not whether it responds to nitrogen *stress* — a gene moved the same
   way by both extremes is invisible to it. At module level that costs one
   module; `Soffic.09G0001580-9H` is a concrete example of what it misses.
2. **n = 18 limits purple everywhere.** Its MI network threshold is soft, its
   trait selection is FDR-bound rather than effect-size-bound, and its module
   structure is one dense component.
3. **The previous networks contained ~47,000 sugarcane edges below their own
   0.8 threshold**, admitted because the correlation matrix was stored as text
   rounded to 4 decimals. Fixed here; the affected results are archived under
   `files/pearson_baseline/` and were never corrected, only superseded.
4. **`r_eq` is not a correlation.** For a non-monotone edge it means "as strong
   as a linear r of this size". Never read it without the `source` column.
5. **Raw conservation rates are not comparable between directions.** Use the fold
   over null.
6. **MI is a network layer, not a universal choice.** It earns its place at the
   edge level, where the question is whether non-linear co-expression exists at
   all. At the module level, where the question is whether one eigengene tracks
   one ordinal trait, it was dropped for Spearman — see
   `new_clean/docs/decisions.md`. The two are consistent; the resolution differs.
7. **The MI layer's value is not settled.** It adds ~0.9 M sugarcane edges the
   linear layer cannot see, but those edges are not better conserved. What is
   robust is that edges *both* estimators find are better conserved than either
   alone.
8. **GO annotation coverage, not statistics, limits the module-level function
   results.** Only 8,251 of sugarcane's 103,336 network nodes and 12,255 of
   purple's 170,736 carry any eggNOG GO term. Combined with a median responsive
   module of 5 genes, that leaves 11% of responsive modules testable (49 of 465;
   63% have *zero* annotated members). The per-module GO section describes those
   49, not the responsive set.
9. **Purple's module count is not a count of independent findings.** Its
   permutation null reaches 244 in the worst of 1,000 shuffles against 79
   observed, because the network is a single dense component. Look at the module,
   not only the padj.
