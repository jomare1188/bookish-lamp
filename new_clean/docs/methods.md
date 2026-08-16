# Methods — what each stage does

One row per stage: what it computes, what it reads, what it writes, roughly what
it costs, and which interpreter it needs. `run.sh` picks the interpreter; the
column is for reference and for when something has to be debugged by hand.

Costs are **measured** on this machine (RTX A4500, ~100 usable cores) from the
2026-08-13 build, not estimated. Sugarcane is ~9x smaller than purple at every
stage — run it first, always.

The whole network build (both studies, VST to merged edge table) is about
**5.5 hours**, and everything through clustering about **9 hours**. Earlier drafts
of this file guessed "days" for MCL on purple; that was extrapolated from stale
directory timestamps and was wrong by an order of magnitude.

---

## 01 · export — `./run.sh export <study>`

DESeq2 object → the flat matrix everything else is computed from.

| | |
|---|---|
| reads | `deseq2.dds.RData` (the `vst` assay and raw `counts`) |
| writes | `results/<study>/vst/<study>.{f32,genes.txt,meta.json}` |
| cost | minutes |
| env | `cor_env` (DESeq2, matrixStats) |

Gene filter: coefficient of variation on **raw** counts, `sd/(mean+1e-6)·100 ≥ 15`,
applied after any library subset. Character-for-character what the original
`pearson_cor.r` did.

Two things happen here and nowhere else:

- **The gene set is decided.** Both network layers read this one file, so they
  cannot disagree about which genes exist.
- **Gene ids are normalised.** The `.v2.1` suffix is stripped, with a collision
  assert. Everything under `results/` carries the bare id.

The `.f32` is row-major float32, `n_genes × n_samples` — 33 MB for sugarcane.
`lib/common.R:read_vst()` reads it back; it round-trips the DESeq2 VST to 9.5e-07.

---

## 02 · network — `./run.sh network <study> pearson|ksg`

All-pairs dependence on the GPU, straight to a thresholded, FDR-corrected edge
list. The n × n matrix is never materialised.

| | |
|---|---|
| reads | `results/<study>/vst/<study>.f32` |
| writes | `results/<study>/layers/<study>_<estimator>.{edgelist.tsv,null.tsv,summary.json}` |
| cost | pearson: **5.5 min** sugarcane, **58.7 min** purple · ksg: **169.7 min** sugarcane, **45.5 min** purple |
| env | `docling` (torch 2.11 + cu130) |

Edge list schema, identical for both estimators — which is what lets the merge
join them without knowing which is which:

```
gene1  gene2  <estimator>  pval  padj
```

**`--estimator pearson`** — the linear layer. Rows centred and scaled to unit
norm, so a tile of correlations is one matrix multiply. Runs at ~1.6 × 10⁹
pairs/s. Signed r is stored; thresholds and p-values apply to |r|. Its null is
the exact two-sided t-test (df = n−2) — no permutations, no tail extrapolation.

**`--estimator ksg`** — the non-linear layer. Kraskov algorithm 1 on rank-
transformed data, Chebyshev metric, fp16 distances (exact: rank differences are
small integers). Its null is built once from up to 2 × 10⁸ permutations and
reused for all 1.5 × 10¹⁰ pairs — legitimate because after ranking, an
independent pair *is* a random permutation pair, so the null depends only on
(n, k). A generalised-Pareto tail is fitted above the 99th percentile to reach
p ≈ 9e-12, which no affordable permutation count resolves empirically.

Two cuts, both p-values:

- `CAND_PEARSON` (0.7) — what gets written to a shard. Sets the disk cost and
  the range over which BH stays exact.
- `STAT_MIN` (0.8) — the network threshold, via `--match-pearson`.

Shards double as checkpoints: `--resume` skips any tile with an `.ok` marker, so
a killed multi-hour run resumes.

---

## 03 · merge — `./run.sh merge <study>`

The two layers → the network. See [thresholds.md](thresholds.md) §5.1 for the
schema and the guards.

| | |
|---|---|
| reads | both layer edge lists, `<study>.genes.txt`, the KSG `null.tsv` |
| writes | `results/<study>/network_<study>_edges.tsv` + `.summary.json` |
| cost | **5.7 min** sugarcane, **51.5 min** purple |
| env | `docling` |

The MI side is held in RAM as sorted int64 gene-index pairs (~30 B/edge) and the
Pearson file is streamed past it once in 4 M-row chunks.

---

## 04 · stats — `./run.sh stats <study>`

| | |
|---|---|
| reads | the network (`gene1`, `gene2`, `weight`, `source`) |
| writes | `network_<study>_{node,global}_metrics.tsv`, degree and strength plots |
| cost | **3 min** sugarcane, **65 min** purple |
| env | `r_net_env` |

`node_metrics` is the most load-bearing file downstream: it is the **background
universe** for the GO enrichment and for the degree-matched permutation nulls in
the MYB61 and Module-20 readouts. A stale copy biases those tests silently, which
is why it is always regenerated rather than carried over.

Local transitivity is off (`COMPUTE_TRANSITIVITY=0`) — O(Σ deg²), ~15 h on purple,
and nothing tests it. Also prints the per-layer edge counts, since it is reading
`source` anyway.

---

## 05 · mcl — `./run.sh mcl <study>`

| | |
|---|---|
| reads | the network (`gene1`, `gene2`, `weight`) |
| writes | `mcl_<study>_{membership,module_summary,hub_genes,plot_data}.tsv`, size plot |
| cost | **13 min** sugarcane, **2 h** purple |
| env | `r_net_env`, plus the native `mcl` binary on `PATH` |

Clustering is delegated to the C `mcl` via a temporary `.abc` file that is
roughly the size of the edge list — tens of GB for purple. `run.sh` points
`TMPDIR` at `MCL_TMPDIR` (`/dados04/jorge/tmp`) so it does not land on a small
system `/tmp`. `INFLATION=2`, `-te 100`.

One study per invocation, so a purple failure cannot lose the sugarcane result.

---

## 06 · conserve — `./run.sh conserve <direction>`

For each edge of network A, project both genes through orthology and ask whether
any resulting ortholog pair is an edge in network B.

| | |
|---|---|
| reads | both networks, `Orthogroups.tsv` |
| writes | `results/conservation/conserved_{edges,genes}_*_FULL.*`, `conservation_summary_*` |
| cost | **~11 min** (`sugarcane_to_purple`); the reverse streams 9x more edges |
| env | `r_net_env` |

Run **twice**, `sugarcane_to_purple` first — it has ~9× fewer outer edges.

The summary is broken down **by layer**, which is the payoff of the whole
augmentation: comparing the `mi` row against the `pearson` row answers whether
non-linear co-expression is conserved between the species at the same rate as
linear co-expression.

Two guards: the ortholog table must overlap the network ids (an empty
intersection produces a perfectly well-formed file saying nothing is conserved),
and the edge file's first two columns must be `gene1`, `gene2` (the streamed read
is positional). Measured overlap on this build: 2,225 of 3,119 sampled sugarcane
genes, i.e. ~71% — a real orthology rate, not an id mismatch.

Orthogroups are parsed by `lib/common.R:read_orthogroups_tsv()`, not
`cogeqc::read_orthogroups()`. cogeqc is installed in three conda envs here and
none is the one the network stages use, so the dependency meant either a fourth
interpreter or a stage that dies on `there is no package called 'cogeqc'`. The
replacement was verified against cogeqc on this project's own file: 360,794 rows,
`Orthogroup` and `Gene` identical element for element.

Lower `CHUNK_SIZE` first if memory presses.

---

## 13 · conservenull — `./run.sh conservenull <direction>`

Permutation null for the conservation rate. Answers whether the rates reported by
stage 06 beat chance, and makes the two directions comparable.

| | |
|---|---|
| reads | both networks, `Orthogroups.tsv` |
| writes | `results/conservation/conservation_null_<direction>.tsv` |
| cost | **13 min** sc→pu (9 of them building the adjacency), **7 min** pu→sc |
| env | `r_net_env` |

The null permutes the ortholog table's target column, preserving fan-out,
coverage and aggregate exposure to the target's degree distribution, so only
*which* gene maps to *which* is destroyed. The target adjacency is built once and
reused across replicates and the source network is Bernoulli-sampled, so a
replicate costs ~7 s rather than the 11–48 min a full pass takes. Only a summary
table is written — no per-edge output, so disk stays flat.

**The fold over null, not the raw rate, is the comparable quantity.** Raw rates
differ 7x between directions purely because purple has 9.3x more edges; the folds
are 2.57 and 2.46. See [results.md](results.md), including what the null does to
the MI-vs-Pearson contrast.

---

## 07 · trait — `./run.sh trait <study>`

Per-gene expression vs trait, restricted to genes on at least one conserved edge.

| | |
|---|---|
| reads | `results/<study>/vst/<study>.f32`, the samplesheet, `conserved_genes_<study>_FULL.txt` |
| writes | `gene_trait_correlations_<study>.tsv`, `selected_genes_treatment_<study>.tsv` |
| cost | fast |
| env | `r_net_env` |

**This stage reads the same VST the network was built from.** In the old tree it
loaded a different DESeq2 object than the purple network came from, and nothing
connected the two well enough for the mismatch to surface — see
[decisions.md](decisions.md).

Trait encoding comes from `TRAITS_<study>` in `config.sh`. Both studies are
encoded so higher = more nitrogen, so the signs are comparable across species.

---

## 12 · traitmi — `./run.sh traitmi <study>`

Gene-vs-trait mutual information. The trait is **discrete**, so KSG does not
apply (it assumes both variables continuous); this uses **Ross (2014)**, KSG
adapted to a discrete class variable, with a tie correction that is mandatory
here because ranking makes every distance a small integer.

| | |
|---|---|
| reads | `results/<study>/vst/<study>.f32`, the samplesheet |
| writes | `gene_trait_mi_<study>.{tsv,null.tsv,summary.json}` |
| cost | **0.5 min** sugarcane, **0.1 min** purple |
| env | `docling` |

Each gene is rank-transformed, so the statistic depends only on which trait label
sits at which rank position. The permutation null is therefore **shared by every
gene**, exactly as in the network layer — which is what makes 10⁷ permutations
affordable, deep enough to resolve the p ≈ 3e-7 that BH needs without
extrapolating.

Output carries both MI and Pearson on the same samples, a `finding` column
(`both` / `mi_only` / `pearson_only` / `neither`), and a `p_source` column
marking whether each p-value came from the empirical null or the GPD tail.
**Read `p_source`.** Genes past the empirical maximum get extrapolated p-values
that rank correctly but whose magnitude is meaningless — a padj of 1e-20 means
"p < 1e-7".

What it found is in [results.md](results.md): real signal on sugarcane, but
distributional (dispersion) rather than the non-linear dose response it was
built for, and no help at all on purple.

---

## 08 · conscor — `./run.sh conscor`

Node-level conservation of the nitrogen response: which ortholog pairs are
trait-correlated on both sides, and do they move in the same direction?

| | |
|---|---|
| reads | `selected_genes_*` from both studies, `Orthogroups.tsv` |
| writes | `conserved_correlated_ortholog_pairs.tsv`, `*_correlated_conservation_status.tsv`, summary |
| cost | fast |
| env | `r_net_env` |

The per-gene `status` column separates `no_ortholog` from
`ortholog_not_correlated` from `conserved_correlated`, so a small result cannot
be misread as absent shared response when it is really absent orthology coverage.

---

## 14-16 · module-level analysis

```
./run.sh eigengene    <study>          # PC1 per module -> a VST-format matrix
./run.sh moduletrait  <study>          # 12_gene_trait_mi.py, unchanged, on it
./run.sh moduleprofile <study>         # + TF hypergeometric per module
./run.sh moduleheatmap <study> [mods]  # heatmaps for the responsive ones
```

| | |
|---|---|
| cost | eigengene ~1 min · moduletrait 0.5 min · profile seconds · heatmaps ~2 min |
| env | `r_net_env`, except heatmaps (`r_env`: ComplexHeatmap, circlize, scico) and moduletrait (`docling`) |

**Why the eigengene matrix is written in the VST export's format.** `.f32` +
`.genes.txt` + `.meta.json`, module names where gene names go. That lets
`12_gene_trait_mi.py` run on modules with no modification, so the module-level
linear/non-linear classification is the *same code path* as the gene-level one.
Its log says "genes" throughout; read "modules".

The eigengene is `prcomp(t(vst_sub), center=TRUE, scale.=TRUE)`, PC1 scores,
oriented so it correlates positively with mean module expression, then z-scored.
Scaling genes before the PCA is the WGCNA convention and matters — without it PC1
chases the highest-variance members instead of summarising the module. Three
fixes relative to the dead `scripts/eigengene.r` it was lifted from are listed in
that script's header.

`MIN_MODULE_SIZE_EIGEN=3` keeps 92% / 87% of each network's genes. Raising it
does not change the eigengenes of larger modules, so it is a cheap re-run.

TF enrichment is a hypergeometric per module against the **network** node
universe, BH across modules; a gene's isoform-driven multi-family calls are
collapsed to one row first or every enrichment is inflated.

Heatmaps draw paired absolute-VST and per-gene z-score panels — a z-score panel
alone rescales a gene varying by 0.01 VST units to look as structured as one
varying by 5. Ramps are percentile-clipped (99th absolute, 98th of |z|). Modules
above `HEATMAP_MAX_GENES` are subset by intramodular strength, stated in the
subtitle; the largest are 19,604 and 47,887 genes.

---

## 09 · go — `./run.sh go BP|MF|CC`

topGO **weight01** Fisher, once per ontology, thresholded on the **raw**
weight01 p-value (`GO_P`, default 0.05) — not on an FDR.

| | |
|---|---|
| reads | eggNOG annotations, `conserved_genes_*_FULL.txt`, `network_*_node_metrics.tsv` |
| writes | `results/conservation/enrichment_conserved/<study>/GO_<ont>_*` + comparison files |
| cost | minutes |
| env | `topGO_env` |

The universe is the **GO-annotated nodes of that network**, not the genome. That
isolates the conservation signal from the generic "co-expressed genes differ from
the genome" effect — and is why `node_metrics` must be current.

Two deliberate choices, both departures from what this pipeline did before:

- **weight01, not `classic`.** `classic` scores each term independently, so a
  specific term's signal propagates up the GO DAG and every ancestor is reported
  as a separate finding. weight01 conditions on the DAG and reports the specific
  terms.
- **Raw p, not BH.** weight01's p-values are deliberately not independent — the
  conditioning is the mechanism — so they are not an exchangeable family and BH
  does not apply to them. This is topGO's own convention. A `p.adj` column is
  still written for reference but does not select terms.

An earlier version also ran BH over only the terms already at `p < 0.05`, which
makes the adjusted values anti-conservative (BH's *m* must be the number of tests
performed). Fixed; see [results.md](results.md) for what the counts did.

---

## 10 · gosem — `./run.sh gosem`

GO semantic-similarity clustering of the enriched terms (Wang, no OrgDb); loops
all three ontologies itself.

| | |
|---|---|
| reads | `enrichment_conserved/<study>/GO_<ont>_conserved_<study>.csv` |
| writes | `enrichment_conserved/semantic/` |
| cost | minutes |
| env | `r_clusterprofiler` |

---

## 11 · readouts

The H1 readouts (TF identification, MYB61 copies, Muñoz Module 20). None of them
read an edge table — only `node_metrics`, `mcl_*_membership`,
`conserved_genes_*_FULL`, `gene_trait_correlations_*` and TPM matrices.

Skip the cached sequence-only steps: `get_tfs/01–03` (hmmsearch, 1–2 h/species),
`myb61/01–06`, `module20/01–02b`. They are network-independent.

`myb61/07` and `module20/03` run 20,000-draw degree-matched permutation tests
against `node_metrics` — the tests most sensitive to a network rebuild.

Their merges use `all.x = TRUE`, so an id mismatch fails as silent NAs rather
than an error. Check the join before trusting the output:

```r
nodes <- data.table::fread("results/sugarcane/network_sugarcane_node_metrics.tsv")
tf    <- data.table::fread("GET_TFS/new/results/sugarcane/TF_in_network.tsv")
sum(sub("\\.v[0-9.]+$", "", tf$gene) %in% nodes$gene)   # must be > 0
```
