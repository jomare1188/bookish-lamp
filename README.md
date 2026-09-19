# Comparative gene co-expression networks in sugarcane, from reused public RNA-seq

Two independently published sugarcane nitrogen experiments, re-quantified through
one pipeline and turned into two comparable gene co-expression networks. The
question the project asks is whether the **nitrogen response is conserved between
them** — at the level of individual genes, of network edges, and of biological
processes.

The secondary goal, and the reason the data are public ones, is to show what can
be recovered from reusing existing sugarcane RNA-seq rather than generating new
libraries.

Networks are **Pearson-only**, at |r| ≥ 0.8 with no degree reduction, built from
the variance-stabilised matrix of each study.

They were not always. The project began with **two layers** — Pearson and
Kraskov–Stögbauer–Grassberger mutual information, thresholded at the same per-edge
false-positive rate so their union was legitimate, with every edge labelled by the
layer that found it. That label existed to answer the central methodological
question: *does non-linear co-expression behave differently from linear
co-expression?*

**It was answered, and the answer was no.** MI-only edges were 1.1% of sugarcane's
network and 4.2% of purple's, and they were conserved across species no better than
linear ones. The layer was dropped rather than quietly kept. The threshold argument
that licensed the union is still in
[`docs/thresholds.md`](new_clean/docs/thresholds.md), because a negative result is
only readable if the thing that produced it is still on record.

---

## Start here

| you want | go to |
|---|---|
| to run the pipeline | [`new_clean/README.md`](new_clean/README.md) |
| the results, with caveats | [`new_clean/docs/results.md`](new_clean/docs/results.md) |
| why it is built this way | [`new_clean/docs/decisions.md`](new_clean/docs/decisions.md) |
| what each stage does | [`new_clean/docs/methods.md`](new_clean/docs/methods.md) |
| the MI threshold argument (historical) | [`new_clean/docs/thresholds.md`](new_clean/docs/thresholds.md) |
| which GO annotation is live, and why | [`annotation/README.md`](annotation/README.md) |
| what in `results/` is stale | [`new_clean/docs/results_tree_STALE.md`](new_clean/docs/results_tree_STALE.md) |
| whether sequence evolution follows network conservation | [`new_clean/docs/dnds.md`](new_clean/docs/dnds.md) |
| the U-shaped MYB, gene by gene | [`new_clean/docs/gene_Soffic_09G0001580_9H.md`](new_clean/docs/gene_Soffic_09G0001580_9H.md) |
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
| **nodes** | **101,990** | **170,135** |
| **edges** | **75,333,769** | **675,955,918** |
| edge density | 0.0145 | 0.0467 |
| components / giant | 958 / 99,851 | **43** / 170,046 (99.95%) |
| mean normalised weight | 0.288 | 0.380 |
| adopted MCL inflation | **-I 1.5** | **-I 3.5** |
| modules (≥ 1 gene / ≥ 3 genes) | 5,650 / 3,627 | 14,024 / 7,493 |
| largest module | 23,439 (23.0%) | 39,230 (23.1%) |

Only 59.7% of sugarcane's CV-filtered genes reach the network, against 99.6% of
purple's. That asymmetry runs through everything downstream and is the reason
per-study percentages are always quoted against the **node** set, never the
proteome.

**Purple is not one connected component.** It has 43 on Pearson edges alone, with
the giant holding 99.95% of nodes; the single-component claim in earlier versions
of this file was true of the merged graph, and the MI layer was what joined those
pieces.

### The clustering, and why the giant module is not a bug

Earlier versions of this file called the MCL partition defective, on the strength
of mcl's own jury synopsis grading its pruning "deplorable" and of a `#knn(180)`
degree reduction that cut the largest module from 18.97% to 0.57%. **That is no
longer the reading.**

Inflation was swept independently at every setting on the unpruned Pearson-only
graphs and scored by Newman modularity — the only one of `clm info`'s four criteria
with an interior optimum, since mass fraction, area fraction and efficiency are each
maximised at a grid boundary. Modularity peaks at **-I 1.5** in sugarcane and
**-I 3.5** in purple, and those are what the pipeline adopts. `-I 2`, used
throughout the earlier work, keeps 98.1% of peak modularity in both.

Leiden (CPM and modularity objectives) and Louvain were run on the same graphs and
scored by the same `clm info` calls. They do not remove the giant: Louvain reaches
Q = 0.1386 and Leiden-modularity 0.1399 in sugarcane, both with giants of 41–43%.
**A large module is what modularity maximisation wants on a graph this dense**, not
an artefact of MCL. The k-NN reduction that appeared to fix it was scoring each `k`
against its own reduced graph; scored against one fixed graph the finding inverts.

Two mcl behaviours are worth knowing if you re-run this: it **silently ignores
`-I` above 30** (`-I 40` returns `-I 2`'s partition byte for byte), and it
**underflows above ~`-I 13`** — 35,044 of 101,990 vectors zeroed at `-I 20` — while
still returning a plausible-looking partition.

### Cross-species edge conservation

An edge is **conserved** when both its genes have orthologs in the other species and
those orthologs are joined by an edge there. Scored against a permutation null that
destroys only the orthology *assignment*, preserving every gene's fan-out, which genes
have any ortholog, and the multiset of targets.

| direction | edges | conserved | fold over null | p |
|---|---|---|---|---|
| sugarcane → purple | 75,333,769 | 7,806,379 (**10.36%**) | **1.69** | 0.0099 |
| purple → sugarcane | 675,955,918 | 10,196,172 (**1.51%**) | **1.71** | 0.0099 |

The raw rates differ ~7× because purple's edge density is 3.2× sugarcane's. **Each
direction is read only against its own null**; they are not comparable to each other.

The folds are lower than the 2.5× earlier versions of this file reported, and that is a
**stricter null, not a weaker signal**: it now draws only from ortholog pairs with both
sides in-network (114,681 pairs over 43,390 orthogroups), where the old one included
genes that could never have matched.

**Conservation rises with edge strength, monotonically, in both species** — and so does
the fold over null, which is what rules out the alternative that strong edges merely join
better-orthologued genes:

| weight decile | sugarcane → purple | purple → sugarcane |
|---|---|---|
| D1 (weakest) | 9.72% (fold 1.57) | 1.17% (fold 1.39) |
| D10 (strongest) | **11.30%** (fold 1.85) | **2.16%** (fold 2.30) |

On a relative reading purple's effect is the larger: +84% across the deciles against
sugarcane's +16%.

**What transfers is housekeeping; what does not is regulation.** GO on the genes carrying
a conserved edge against the exact complement — the two partition the node set and share
a background — gives translation, protein folding and transport, photosynthesis and
splicing on one side, and transcriptional regulation, auxin and ethylene signalling,
ubiquitination and the cell cycle on the other.

> Read that with the caveat: genes on a conserved edge **are hubs** (median degree 228
> against 18 in sugarcane), partly by arithmetic, since one conserved edge out of 6,677
> is near-certain at a 10.4% per-edge rate and unlikely out of 2. Ranking the same genes
> by degree instead reproduces nearly the same split, so the two are one contrast seen
> twice and neither is independent evidence for the other.

### Nitrogen response

The gene-level test puts the **experimental design in the model** rather than in the
residual. The project's own QC records genotype at R² = 0.999 of PC1 in purple and 0.998
in sugarcane, and leaf segment at 0.802 of PC2 — so the largest variance component in
either matrix used to sit in the residual of every nitrogen test.

    sugarcane   expr ~ genotype + segment + N            n = 48, resid df 42, Pearson
    purple      expr ~ genotype + N        on midranks   n = 18, resid df 15, Spearman

Purple is Spearman because its trait is an **ordinal 0/2/6 mM dose** and Pearson reads
that spacing literally. The universe is the network's node set, so the BH denominator is
101,990 and 170,135.

| | sugarcane | purple |
|---|---|---|
| genes tested | 101,990 | 170,135 |
| responsive, **blocked** | **8,737** | **3,854** |
| responsive, marginal, same genes | 3,238 | 882 |
| \|r\| at the BH boundary | 0.366 | 0.720 |

Blocking is a **strict superset** in both species — no marginal call is lost — so the
block is removing noise rather than absorbing signal. The diagnostic that motivated it:
purple's marginal p-value histogram *rises* in its last decile (12.5%), which a mixture
of a uniform null and real signal cannot do. Blocking takes it to 10.2% — better, and
still not flat.

Purple also gets a **non-monotone tier**, because three nitrogen levels let a gene
respond to deficiency and excess alike, which every monotone test is blind to. Over
170,135 genes it selects **one**. The aggregate excess is real (17,507 raw hits against
8,507 expected) and n = 18 cannot say which genes carry it.

### Do the two species respond in the same genes?

Over the 114,681 ortholog pairs with both sides in-network:

| | genome-wide | directed (discover in sugarcane) |
|---|---|---|
| conserved correlated **pairs** | 392 (270 orthogroups) | 616 (419) |
| null | 314.3, **fold 1.25**, p 0.001 | 573.8, fold 1.07 |
| sugarcane **genes** responsive in both | **326** | 502 |
| null | 301.9, **fold 1.08**, p **0.071** | 536.3, fold **0.94** |

**The pair count beats chance; the gene count does not.** Orthology is many-to-many —
one orthogroup contributes 6 pairs from 3 × 2 genes — so most of the pair excess is
multiplicity. The pair row must never be quoted without the gene row beside it.

What does hold up is **direction**. Of the concordant pairs 94 are up in both species and
150 down in both, and the agreement survives a per-orthogroup test that treats each
orthogroup as one vote rather than counting non-independent pairs: **59.5% of 269
orthogroups, p = 0.0022** (57.7% of 418 under the directed design).

### Module-level nitrogen response

One eigengene per MCL module of ≥ 3 genes cuts the testing burden from 101,990 and
170,135 genes to **3,627 and 7,493**, and a module signal can survive where a single
gene cannot.

The statistic is a **blocked Spearman partial correlation** — eigengene and trait both
reduced to midranks, then residualised on the design — at `padj ≤ 0.05` and
`|rho| ≥ 0.6`, the same thresholds the gene level uses. Only the model changes.

| | sugarcane | purple |
|---|---|---|
| modules tested | 3,627 | 7,493 |
| responsive, **blocked** | **588** | **182** |
| responsive, marginal, same eigengenes | 527 | 110 |

Blocking is again a strict superset, and purple gains most — 110 → 182 — which is where
the ordinal-trait argument predicts it. Sugarcane's 48 libraries are 12 plants × 4 leaf
segments, so a **plant-level control** (segments averaged, n = 12) is computed alongside
and agrees with the blocked fit at r = +0.9486.

### What the responsive modules do

GO comes from **one source**: a full local InterProScan over all 17 member databases,
mapped through the GO Consortium's pinned `interpro2go`/`pfam2go`. Coverage of network
nodes went from 8.0% / 7.2% under the original eggNOG GO column to **62.9% / 64.4%**.

**Merging sources was rejected by its own judge.** On a fixed gene set — identical genes,
identical modules — the union of InterPro and eggNOG scores +0.0671 in sugarcane against
+0.0748 for InterPro alone, and +0.0264 against +0.0361 in purple. Merging *lowers*
module coherence in both species, so there is no tiered table. eggNOG's raw homogeneity
is 2.6× InterPro's while its excess over the null is lower — what the raw figure measures
is term density (16.7 terms/gene against 2.1), not shared function.

| | responsive | testable | BP | MF | CC |
|---|---|---|---|---|---|
| sugarcane | 588 | **479** | 164 | **272** | 22 |
| purple | 182 | **118** | 9 | **17** | 1 |

*(terms clearing cross-module BH)*

The annotation gate that once left only ~11% of responsive modules testable now leaves
81% and 65%. **MF is the stronger ontology** — expected for a domain-derived annotation,
which names what a protein *does* more sharply than what process it is in. BP stays
primary because the question is about a response.

**Purple's 9 BP terms are the honest negative of this section.** More annotation did not
fix it: the count went 14 → 8 → 9 across three successive annotations while coverage went
7.2% → 55.5% → 64.4%. The constraint there is n = 18 and a network that is one dense
component, not the annotation.

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
purple Module-20 genes, and the U-shape reaches significance in only **one of the two
purple genotypes** (51NG3 p = 0.008; TAGZ p = 0.142, same direction). Run genome-wide
over every purple network gene, the same contrast returns one gene past FDR and this is
not it — it sits at padj 0.48, in the top 0.5% by raw p.

**As a network object it does not corroborate the expression evidence either.** It is not
a hub — 635 neighbours, the **48.2nd** degree percentile, below purple's median — while
Muñoz define Module 20 as high-betweenness on their own network.

The nine sugarcane copies are one locus only by Arabidopsis anchor: OrthoFinder splits
them across at least four orthogroups, and only `06Ag012300` (OG0088735) pairs 1:1 with
the purple gene. Projecting each copy's neighbourhood through orthology and intersecting
with the purple gene's 635 neighbours, against a degree-matched null, **nothing clears
BH** — and the 1:1 ortholog is *specifically* indistinguishable from chance (11 shared
genes against 8.4 expected, z = 0.53, p = 0.29). No copy recovers more than 1.9%.

The measure is not at fault: within sugarcane those same copies share up to **83%** of
the smaller neighbourhood. And the neighbourhoods do different things — the purple gene's
is stress and specialised metabolism (tryptophan and terpene synthase, innate immune
response, wounding), its sugarcane ortholog's is core housekeeping (vesicle fusion,
translation initiation, TCA cycle).

So: **a candidate worth following, not a finding** — and the reasons for caution are now
independent, the expression evidence being thin *and* the network evidence not
corroborating it.

### Function and H1 readouts

- **GO of the conserved gene sets** (topGO weight01, raw p ≤ 0.05, background = each
  network's own GO-annotated nodes, on the adopted annotation): **111 shared BP terms**
  (165 sugarcane, 207 purple), 119 MF, 31 CC. Semantic clustering puts every one of the
  12 BP macro-themes in both species — protein folding, photosynthesis, organelle
  organisation, vesicle transport, splicing. The convergence is about **themes, not
  terms**: the shared fraction is under half in every ontology, so the two species rarely
  enrich the identical term and land in the same regions of GO space.
- **Transcription factors**: 7,183 (sugarcane) and 12,197 (purple) TF genes in
  the networks, 68 and 67 families.
- **MYB61**: purple's 15 copies sit on cross-species conserved edges far more than
  expected — 60% vs 25.8% background, **degree-matched permutation p = 0.0010**,
  and *not* a hub artifact (their median degree is below average). Sugarcane's 8
  copies show nothing (p = 0.67).
- **Muñoz Module 20**: no network-position signal in either species. The module's
  nitrogen response transfers between studies; its network position does not.
  See above for the one purple copy that does respond.

### An alternative clustering, tested on the previous network

> **This section describes the PRE-PEARSON clustering** — MCL at `-I 2` on the merged
> network, 10,309 modules. It is kept because the comparison was real and the
> annotation-gate argument is what motivated the GO work that followed. The clustering
> question was settled differently in the end: see
> [The clustering, and why the giant module is not a bug](#the-clustering-and-why-the-giant-module-is-not-a-bug)
> above. The SBM was never re-fitted on the Pearson-only graphs, and the numbers below
> are **not** comparable with anything in the sections above.

Every module-level result at the time rested on MCL's partition, and that partition had
an awkward shape: **10,309 modules with a median of 3 genes**, one holding 19% of
the network. A median of 3 sits below the GO annotation gate almost by
construction, which is why only 49 of the 465 responsive modules could be tested
for function at all.

A stochastic block model has been fitted to the same sugarcane network and
carried through the *identical* downstream analysis — same eigengenes, same
Spearman test at the same thresholds, same TF and GO code:

| | MCL | SBM |
|---|---|---|
| modules / median size / largest | 10,309 / 3 / 19,604 | 995 / **51** / 1,366 |
| genes left unassigned | 401 | 14 |
| responsive modules | **465** | 23 |
| of those, GO-testable | 49 (**11%**) | 20 (**87%**) |
| BP terms returned | 246 | **325** |

The two partitions barely overlap — adjusted Rand index **0.0156**, verified
against `clue::cl_agreement` — so this is a real alternative, not a
re-parameterisation. (The modularity figures once quoted beside this table
compared an unweighted MCL Q with a weighted SBM one; both stages now write both,
and the purple SBM fit referred to elsewhere never actually ran.) The trade is sharp: MCL finds twenty
times more responsive modules, of which nine in ten cannot be tested for
function; the SBM finds 23, of which almost all can, and they return more GO
terms than MCL's 465 do.

**Nothing was switched.** MCL remained the default then and still is — though for a
different reason than "the SBM never finished". Sweeping inflation on the unpruned
Pearson-only graphs and scoring Leiden and Louvain on the same footing showed the giant
module is what modularity maximisation *wants* on a graph this dense, not a defect MCL
introduced. The annotation gate that motivated this whole comparison was then largely
closed from the other end: GO coverage went from 8% to 63%, and 81% of sugarcane's
responsive modules are now testable rather than 11%.

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
sbm/                  stochastic block model fits (gitignored, 24 GB — the
                      pipeline reads one 8 MB node-block table from it)
GET_TFS/              TF identification: HMM databases and cached hmmsearch output
annotation/           eggNOG-mapper annotations for both references
china/, run1/         nf-core/rnaseq output for the two studies (gitignored)
```

## Reproducing

```bash
cd new_clean
./run.sh validate                     # correctness suite — must pass first
./run.sh build sugarcane              # VST -> Pearson layer -> network
./run.sh build purple

# the Pearson-only graph, clustering, and node metrics, per study
./run.sh pearsonmci  <study>          # Pearson layer -> MCL matrix, direct
./run.sh mclladder   <study>          # inflation ladder, scored one cell per call
./run.sh membership  <study>          # adopt the chosen cell (-I 1.5 / -I 3.5)
./run.sh nodemetrics <study>

# gene level (blocked), per study
./run.sh traitblocked <study>
./run.sh ushape purple                # the non-monotone tier

# module level, per study
./run.sh eigengene moduletrait moduleprofile modulego moduleheatmap modulesummary

# conservation, on the Pearson-only graphs
./run.sh consblocked 0                # node level; 1 = the directed variant
./run.sh consedges sugarcane_to_purple
./run.sh consedges purple_to_sugarcane

# function
./run.sh go BP [conserved|nonconserved]
./run.sh gosem
./run.sh degreego <study>             # what hubs vs the periphery are for
```

`conserve`, `conservenull`, `trait` and `conscor` are the **merged-network** stages,
kept as the record of that analysis. They cannot run on the current graphs — the
gene-name edge tables were deleted — and are superseded by `consedges`, `consblocked`
and `traitblocked`.

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
2. **n = 18 limits purple everywhere.** Its trait selection is FDR-bound rather than
   effect-size-bound (|rho| 0.720 at the boundary against sugarcane's 0.366), its
   non-monotone tier resolves one gene out of 170,135, its module GO clears 9 BP terms,
   and its network is effectively one dense component. Almost every thin purple result in
   this repository traces back to that number.
3. **The previous networks contained ~47,000 sugarcane edges below their own
   0.8 threshold**, admitted because the correlation matrix was stored as text
   rounded to 4 decimals. Fixed here; the affected results are archived under
   `files/pearson_baseline/` and were never corrected, only superseded.
4. **`r_eq` is not a correlation.** For a non-monotone edge it means "as strong
   as a linear r of this size". Never read it without the `source` column.
5. **Raw conservation rates are not comparable between directions.** Use the fold
   over null.
6. ~~**The MI layer's value is not settled.**~~ **Settled, and it lost.** MI-only edges
   were 1.1% and 4.2% of the two networks and conserved no better than linear ones. The
   layer is not used anywhere in the current analysis; `docs/thresholds.md` keeps the
   argument that licensed it.
7. ~~**GO annotation coverage limits the module-level function results.**~~ **Largely
   fixed.** Coverage went from 8.0% / 7.2% of network nodes to **62.9% / 64.4%** on a
   full 17-database InterProScan, and testable responsive modules from 11% to **81% and
   65%**. What did *not* follow is significance in purple: 9 BP terms clear cross-module
   BH, and that count went 14 → 8 → 9 across three successive annotations. The limit
   there is n = 18, not the annotation.
8. **Edge conservation and degree are not separated.** Genes on a conserved edge are hubs
   (median degree 228 against 18 in sugarcane), partly by arithmetic. So the
   conserved/non-conserved functional split and the hub/periphery split are one contrast
   seen twice. Whether conservation selects housekeeping genes *over and above* their
   being hubs would need a degree-matched comparison, which is not done.
9. **Purple's module count is not a count of independent findings.** Its
   permutation null reaches 244 in the worst of 1,000 shuffles against 79
   observed, because the network is a single dense component. Look at the module,
   not only the padj.


10. **The conserved *gene* count does not beat chance.** 1.08× at p = 0.071 genome-wide,
    and 0.94× — below the null — under the directed design. Only the *pair* count beats
    it (1.25×), and that is ortholog multiplicity. The conservation claim rests on
    **direction** (59.5% of orthogroups agree in sign, p = 0.0022), not on how many genes
    respond in both species.
11. **Purple's blocked p-value histogram is still not flat.** Blocking took its last
    decile from 12.5% to 10.2%, but the shape dips to 5.4% and rises again. Genotype plus
    nitrogen does not account for everything at n = 18, and any purple gene-level claim
    should be read with that in mind.
12. **`Soffic.09G0001580-9H` is a candidate, not a result.** Its U-shaped response holds
    in one of two purple genotypes, does not survive genome-wide correction, and its
    network neighbourhood does not corroborate it — see
    [Module 20 in purple](#module-20-in-purple--one-copy-responds-and-not-monotonically).
