# `Soffic.09G0001580-9H` — a MYB with a U-shaped nitrogen response

A per-gene dossier, assembled from the pipeline's own outputs. It exists because
this gene is the clearest single illustration in the project of the difference
between *"no signal"* and *"the test cannot see this shape"* — and because it is
the concrete case that motivated running the U-shape contrast genome-wide
([results.md](results.md#and-it-survives-a-non-monotone-test-which-is-the-one-blind-spot-that-mattered)).

**Read the caveats at the bottom before quoting any of it.**

---

## Identity

| | |
|---|---|
| annotation (eggNOG) | Myb-like DNA-binding domain, PFAM `Myb_DNA-binding` |
| COG category | **K** — transcription |
| seed ortholog | `4577.GRMZM2G013581_P01` (maize), e = 1.27e-112 |
| this project's TF call | **MYB**, family-confirmed (`get_tfs/purple`) |
| GO terms | none assigned |
| sorghum anchor | `LOC110430376` |
| aligned length | 292 codons |
| Module-20 locus id | `SCA2_2__c76820f1p01331` |

## The expression pattern

Mean TPM per genotype and nitrogen level, with fold-change against the 2 mM
control (`china/run2_onlyL` salmon quantification, 3 replicates per cell):

| genotype | 0 N (deficiency) | 2 N (control) | 6 N (excess) |
|---|---|---|---|
| **51NG3** (*S. robustum*) | **124.1  (11.5x)** | 10.8 | **47.4  (4.4x)** |
| TAGZ (*S. officinarum*) | 65.5  (1.6x) | 40.0 | 75.6  (1.9x) |

In 51NG3 it is induced eleven-fold by nitrogen starvation and four-fold by
nitrogen excess, with its minimum at the control. Expression is comfortably above
the noise floor throughout, so this is a shape, not a detection artefact.

## Why every monotone test missed it

| test | result | verdict |
|---|---|---|
| Pearson (the pipeline's main rule) | r = **-0.0895**, p = 0.72, padj = **0.9916** | invisible: the two arms cancel |
| mutual information | 0.322, p = 0.0165, padj = 0.6337 | a hint that dies under correction |
| the pipeline's own `finding` field | **"neither"** | |
| **U-shape contrast** `c(+1,-2,+1)` | u_est **+3.669**, p = **0.003676**, padj = 0.4849 | **rank 890 / 170,740 — top 0.52%** |
| per genotype | 51NG3 +5.292 (p = 0.0115); TAGZ +2.046 (p = 0.0929) | the response is genotype-restricted |
| among purple TFs only | | **rank 73 / 12,197 — top 0.60%** |

Pearson scores it at r = -0.09, indistinguishable from a flat line. That single
number is the argument for why the U-shape contrast had to be run at all.

## Network position

**Purple** — degree **681** (47.6th percentile; the network median is 859),
strength 133.8, `Module_003`, **not a hub**. 12 of its 681 edges are conserved
(1.76%). `Module_003` holds 5,558 genes and is **not** nitrogen-responsive
(rho 0.17, padj 0.85), so the module level does not rescue it either.

**Its 1:1 sugarcane partner**, `SoffiXsponR570.06Ag012300` — also a **MYB-related
TF**, and structurally a very different gene: degree **2,349** (85th percentile,
against a median of 52), strength 464.5, local clustering 0.608, coreness 1,888,
and **flagged as a hub**. Its nitrogen response is monotone repression:
**r = -0.491, padj = 0.0069**.

## Why the conserved pair still fails — twice, for different reasons

Both partners sit on conserved edges, so the data are present. It fails on
thresholds, on both sides at once:

- **purple**: U-shape padj **0.485** — a top-0.5% candidate, not correctable at n = 18
- **sugarcane**: r = **-0.491** — clears padj, but sits under the **|r| >= 0.6** effect floor

There is a third layer. Under the three-species orthology this locus is
**1 copy in purple, 5 in sugarcane**, and three of the other sugarcane copies are
*strongly* nitrogen-repressed — r = **-0.871** (padj 1.96e-12), **-0.842**
(6.3e-11), **-0.716** (1.7e-06) — but **none of them lies on a conserved edge**,
so the conservation funnel excludes them before any test is applied.

## Sequence evolution

| | |
|---|---|
| sugarcane <-> purple | dN = **0**, dS = **0** (identical over all 292 codons) |
| omega vs sorghum | **0.2247** — percentile 61.7, unremarkable |
| GC3 | **0.853** (the GC-rich grass class) |
| **constraint_score** | **+1.322 — percentile 98.7** |

Two things here. The sugarcane and purple proteins are **100% identical**, so the
divergence between the two study species tells us nothing about this locus.

And it is a good demonstration of why the rate-corrected readout exists
([dnds.md](dnds.md#what-actually-drives-omega-here-mostly-the-denominator)). Its raw
omega of 0.22 looks ordinary. But GC3 = 0.85 inflates dS and *should* have pushed
omega much lower; once that is corrected the gene sits in the **top 1.3% fastest
evolving** in the analysis — accumulating far more amino-acid change than its
synonymous rate and composition predict.

## What this supports, and what it does not

**Supports.** A MYB transcription factor, single-copy in purple and five-copy in
sugarcane, induced eleven-fold by nitrogen starvation *and* four-fold by nitrogen
excess in *S. robustum* — a stress-response shape that every monotone test in
this pipeline scores as "neither". Its sugarcane counterpart is a network hub and
a nitrogen-*repressed* MYB. The proteins are sequence-identical between the two
study species, while the locus is among the fastest-evolving once composition is
accounted for.

**Does not support, and must not be written as.**

1. **It is not a significant finding.** padj = 0.485 under the U-contrast. It is
   the best-characterised member of a class this design cannot resolve, and that
   is the sentence to use. The genome-wide run found exactly **1** gene surviving
   BH out of 170,740, and it was a different gene and the opposite shape.
2. **It is not evidence of a conserved response.** The purple copy is U-shaped
   and the sugarcane copies are monotonically repressed. That is the same TF
   family doing *different* things in the two species — arguably the more
   interesting claim, and the one the data actually carry.
3. **The genotype restriction is real but underpowered.** The response is clear in
   51NG3 (p = 0.0115) and weak in TAGZ (p = 0.0929), on n = 9 each.

## How to regenerate every number here

```bash
./run.sh ushape                    # genome-wide U-shape contrast, purple
./run.sh conscor 0 ushape          # conserved pairs, genome-wide burden
./run.sh conscor 1 ushape          # conserved pairs, directed burden
```

Sources: `results/purple/gene_trait_ushape_purple.tsv`,
`gene_trait_mi_purple.tsv`, `network_purple_node_metrics.tsv`,
`mcl_purple_membership.tsv`, `module_trait_purple.tsv`,
`results/dnds/dnds_gene_table.tsv`, `results/dnds/centrality_sugarcane.tsv`,
`results/readouts/get_tfs/purple/TF_genes_families.tsv`,
`annotation/purple/emapper.annotations`.
