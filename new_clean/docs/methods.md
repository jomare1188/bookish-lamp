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

## 05b · an alternative clustering — `./run.sh sbmclust <study>`

A stochastic block model was fitted to the sugarcane network (`sbm/` at the repo
root, outside this pipeline) as an alternative to MCL. `27_sbm_membership.r`
converts that fit into **the two files every module-level stage reads**, in MCL's
exact schema, so the alternative can be carried through the identical downstream
analysis and the two compared on their biology.

| | |
|---|---|
| cost | ~6 min, almost all of it the modularity pass over the 76 M-row edge table |
| env | `r_net_env` |
| inputs | `<SBM_DIR>/<SBM_TAG>_node_blocks.tsv`, the node-metric table, the edge table |

**The contract.** `<prefix>_membership.tsv` (`gene, module_name, strength,
degree`) and `<prefix>_module_summary.tsv` (`module, n_genes, modularity_Q, …`),
with `module_name` ranked largest-first as `Module_%03d` or the literal
`Unassigned` — every convention copied from `05_mcl_clustering.r`, because
`22_fig_topology.r` separates real modules from the leftover pool with
`grepl("^Module_", module)` and `14_module_eigengene.r` drops `Unassigned` by name.

**No edge table is read for `strength`/`degree`.** MCL's own columns are
*whole-graph* igraph values, not intramodular ones, which makes them identical to
the columns already in `network_<study>_node_metrics.tsv`. Only one consumer uses
`strength` at all — `16_module_heatmaps.r`, to take the top N genes of a module
too large to draw.

**Modularity is computed anyway**, on the same graph, because it is the number
`05_mcl_clustering.r` reports for MCL and a comparison needs both sides measured
the same way.

> **Correction, 2026-09-03: they were not measured the same way.** This stage
> passed `weights = E(g)$weight` and `05_mcl_clustering.r` did not, and igraph
> does not pick up the `weight` attribute when `weights` is NULL — so the quoted
> "SBM 0.0081 against MCL 0.1005" compared a weighted modularity with an
> unweighted one. Both stages now write **both**, as `modularity_Q` (unweighted)
> and `modularity_Q_weighted`; read one column, never a mixture. The
> qualitative point stands and is unaffected: an SBM minimises a description
> length and does not optimise modularity, so a lower Q is expected rather than
> a defect.

`SBM_COMPUTE_Q=0` skips the slow pass and writes NA.

The script **fails** rather than warns if the fit's gene set and the network's
disagree, because a partition of a different graph would silently produce modules
that are not the ones the fit found.

### The `CLUSTERING` switch

`CLUSTERING=mcl|sbm` in `config.sh`, honoured from the environment
(`CLUSTERING="${CLUSTERING:-mcl}"`) so `CLUSTERING=sbm ./run.sh eigengene
sugarcane` works. Two helpers derive every path:

```sh
clus_prefix()   # which membership/summary the module stages READ
module_dir()    # where their outputs are WRITTEN
```

With `mcl` both are exactly what they have always been, so the existing results
cannot be disturbed; with `sbm` the outputs land in `results/<study>/sbm/`.

> **run.sh announces the active clustering** on every module-level stage. That
> line exists because the first attempt at this went wrong silently: `config.sh`
> assigned `CLUSTERING=mcl` unconditionally, so an exported `CLUSTERING=sbm` was
> clobbered on sourcing and three stages re-ran MCL into the MCL location while
> appearing to do SBM work. Nothing was lost — the outputs are deterministic and
> reproduced byte-identically — but only a checksum proved it. A stage using the
> wrong clustering writes plausible output to a plausible place, and being told
> which one is active before the work starts is the only cheap defence.

> **The Module-20 readout stays on MCL.** `module20/03_network_readout.r` asks
> which of *our* modules hold the Module-20 genes — a question about one specific
> clustering — and `run.sh` invokes that stage with no environment at all. Its
> membership paths are now `Sys.getenv` with the MCL files as defaults, so it
> cannot drift with the switch; point `CLEAN_M20_MODS_SUGARCANE` at another
> membership file to re-run it against one.

## 14-19 · module-level analysis

```
./run.sh eigengene    <study>          # PC1 per module -> a VST-format matrix
./run.sh moduletrait  <study>          # Spearman rho vs the trait, on it
./run.sh moduleprofile <study>         # + TF hypergeometric per module
./run.sh moduleheatmap <study> [mods]  # per-module gene heatmaps
./run.sh modulesummary <study>         # one figure: all responsive modules
./run.sh modulego     <study> [ONT]    # topGO per responsive module (default BP)
```

| | |
|---|---|
| cost | eigengene ~1 min · moduletrait ~10 s (incl. 1,000 permutations) · profile seconds · heatmaps ~1 min · modulego ~45 s |
| env | `r_net_env`, except heatmaps and modulesummary (`r_env`: ComplexHeatmap, circlize, scico) and modulego (`topGO_env`) |

**Why the eigengene matrix is written in the VST export's format.** `.f32` +
`.genes.txt` + `.meta.json`, module names where gene names go. Anything that can
read the VST export can read the eigengenes — `read_vst()` in `lib/common.R`
does, unchanged, in three of these five stages. (It also let the module response
be computed by `12_gene_trait_mi.py` when that was the statistic; it no longer
is, but the format is the reason swapping the statistic cost one script and no
plumbing.)

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

### The response call — `19_module_trait_spearman.r`

**One statistic: Spearman's rho.** `padj <= MODULE_PADJ_THR` (0.05, BH over every
module in the study) **and** `|rho| >= MODULE_R_THR` (0.6), mirroring the gene
level's `TRAIT_PADJ_THR` / `TRAIT_R_THR` so the two sets of counts are comparable.
The rule lives in this one script, with the statistic it applies to; nothing
downstream recomputes it, they read `responsive` and `direction`.

Pearson and mutual information are gone from this level. The trait is ORDINAL —
purple's nitrogen is a dose (0/2/6 mM) and Pearson reads that spacing literally —
and MI is an omnibus test that at the module level was firing on two-library
quirks. The measured consequences of both, including the 38 → 79 jump in purple's
responsive set, are in [results.md](results.md#why-one-rank-correlation-and-not-the-three-statistics-this-used-to-report).

**Ties.** Both traits are heavily tied by design (sugarcane 24/24, purple 6/6/6),
so the exact/AS-89 p-value is invalid and is not used. `rho` is Pearson on
MIDRANKS and the p-value is the asymptotic t approximation on it — exactly what
`cor.test(method = "spearman", exact = FALSE)` computes. The script verifies that
against `cor.test` on 25 modules every run and dies if it disagrees beyond 1e-8.

**The null** is 1,000 trait-label permutations against the same eigengenes,
counting how many clear both thresholds. It is a null for the TRAIT association
only: the eigengenes stay exactly as correlated with each other as they really
are, which is why purple's null has a long right tail worth reading before any
single purple module is believed.

`pval` is written alongside `padj` so an uncorrected reading needs no re-run.

### Figures — 16 and 17

Heatmaps show **per-gene z-scores only**. Ramps are clipped at the 98th
percentile of |z|. Device size is derived from the actual body size
(`CELL_W x CELL_H` per cell plus fixed overhead), so height scales with gene
count instead of every figure being forced into one frame. Modules above
`HEATMAP_MAX_GENES` are subset by intramodular strength, stated in the subtitle.
`HEATMAP_TOP_N` modules are drawn per DIRECTION, so the smaller of "rises with
nitrogen" and "falls with nitrogen" is not crowded out by the larger.

**Column order groups replicates, on both figures.** Sorting by sample name alone
interleaves sugarcane's three leaf positions — the names run `B0_1, B_1, M_1,
P_1, B0_2, ...`, so the three replicates of one tissue land four columns apart
and the tissue effect reads as vertical striping right across the panel,
obscuring the treatment pattern the figure exists to show. Columns are ordered by
the split trait, then the other trait, then `HEATMAP_GROUP_BY` (tissue), then the
sample name.

`modulesummary` draws every responsive module's eigengene in one panel, **split by
the sign of rho** — rises with nitrogen above, falls below — with PC1 variance
explained, module size and TF enrichment as row annotations. Sample labels are
dropped (the design is carried by the annotation bars), which lets the cells
shrink to 2 mm. Above `SUMMARY_MAX_MODULES` rows the top N per direction by
significance are shown and the title says so.

### 18 · per-module GO — `./run.sh modulego <study> [BP|MF|CC]`

`09_go_enrichment.r` asks what the *conserved gene set* is for, one test per
species. It cannot name an individual module. `18_module_go.r` does: **one topGO
run per responsive module**, so `Module_100` can be reported as nitrate
assimilation rather than only as `positive, rho = 0.81`.

**Every responsive module goes in as one set.** There are no response classes to
pool any more — the Spearman call produces one responsive set. `direction` (rises
or falls with nitrogen) rides through as a column so it can be split afterwards,
but it does not partition the run: the question is what responsive modules *do*,
not what separates the two halves.

Background, algorithm and threshold are **identical to 09** and deliberately so:
GO-annotated nodes of that species' network as the universe (the same denominator
as 15's TF hypergeometric), `weight01` Fisher, selection on the **raw** p at
`GO_P = 0.05`. Keeping all three aligned is what lets "this module is TF-rich"
and "this module is enriched for nitrate assimilation" be statements about one
population. `parse_eggnog()` is carried over from 09 unchanged rather than
reimplemented.

**One `topGOdata` object, reused via `updateGenes()`.** 09 rebuilds the object
per gene set, which re-runs the DAG mapping every time. Reusing it is not only
the difference between ~45 seconds and hours — it means **every module is scored
against an identical term universe** (2,412 BP terms sugarcane, 2,561 purple), so
the per-module results are comparable to each other. Rebuilding would let the
tested term set drift with the gene set. Verified bit-identical to a rebuild on
three modules before adopting it.

The loop reads `score(runTest(...))` and `termStat()` directly instead of
`GenTable()`. Same numbers, without the round trip through topGO's formatted
strings (`"< 1e-30"`, 4 significant digits) that 09 has to parse back.

**Three p-value columns, one of which selects.**

| column | what it is |
|---|---|
| `pvalue` | raw weight01. **This selects**, at `GO_P`. |
| `p.adj` | BH within the module, over every tested term — 09's convention |
| `p.adj_global` | BH across **all** module × term tests in the study |

`p.adj_global` is the burden per-module testing introduces and 09 never faced:
171,252 tests in sugarcane, 7,683 in purple. It is reported, not used to select —
weight01 p-values are not an exchangeable family and pooling them across modules
compounds that — but a term surviving it is on far firmer ground. Computing it
requires the *complete* p-vector of every module, not just the retained terms,
because BH takes a cumulative minimum from the largest p downwards.

`MODULE_GO_MIN_ANNOTATED = 3` gates which modules are testable. It is
deliberately low: a 3-gene module *can* reach p < 0.05 against an 8,251-gene
background, so the gate only skips modules where the test is undefined. **In
practice it is the binding constraint on the whole stage** — see
[results.md](results.md). `n_annotated` is written for every module, gated ones
included, so filtering harder needs no re-run.

**Two grains of output.** The joined table answers "does the responsive set have
a coherent function"; the per-module directories answer "is *this* module worth
following up". Different questions, so each gets its own artefact and its own
figure, in `results/<study>/module_go/`:

```
module_GO_<ONT>_<study>.tsv           joined — one row per term per module
module_GO_<ONT>_<study>_summary.tsv   one row per responsive module
module_GO_<ONT>_<study>_global.{png,pdf}      GLOBAL figure
modules/<Module_NNN>/GO_<ONT>_<Module_NNN>_<study>.tsv
modules/<Module_NNN>/GO_<ONT>_<Module_NNN>_<study>.{png,pdf}   GRAIN figure
```

The summary carries **one row per responsive module, including the ones that
returned nothing and the ones the gate excluded** — the denominator is the result
here, so it cannot be dropped. A directory exists for every *tested* module, so
its presence means "tested" and a missing figure inside it means "nothing cleared
the threshold", without consulting the summary.

Both figures are dot plots in 09's idiom, with its p = 0 flooring carried over
(`-log10(0)` is infinite and silently drops the best term off the panel).

- **Grain**: one module, terms by `-log10(p)`, sized by how many of the module's
  genes carry the term, and **coloured by whether the term also survives the
  cross-module BH** — the honest part of the panel, because most do not.
**All three ontologies run on the same universe.** `parse_eggnog` keeps every GO
id regardless of ontology, so `geneUniverse` — and therefore which modules clear
`MODULE_GO_MIN_ANNOTATED` — is identical for BP, MF and CC (71 sugarcane, 3
purple). topGO filters to the requested ontology when it builds the graph, exactly
as in 09. Because the outputs are ontology-tagged, `modules/<Module_NNN>/` ends up
holding all three side by side, which is the useful unit: the BP, MF and CC panels
for one module answer "what process", "what enzyme" and "where" about the same
gene set.

- **Global**: terms **ranked by how many distinct modules** they are enriched in —
  a term found once is a lead, a term found in six independent modules is a
  pathway. Recurrence gets the size and colour channels; `-log10(best p)` is the
  x axis rather than recurrence itself, because when nothing recurs (purple, where
  3 modules share no term) a recurrence axis collapses onto one value and the
  figure says nothing. In that case the degenerate legend is dropped and the
  subtitle states it.

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

---

## 20-26 · Paper figures

### The rule every paper figure follows

**Nothing goes on the figure that belongs in the legend.** A panel carries a bold
letter and the labels the data needs to be read — split titles, annotation names
that are not already named by a colour key, axis labels, the keys themselves —
and nothing else. No titles, no subtitles, no statistics printed onto the panel.

**Scripts are named for what they draw, not for their figure number.**
`20_fig_reproduction.r`, `21_fig_dataset_qc.r`. Figure order is editorial and has
already moved once — the reproduction figure opened the paper until the
dataset/QC figure took the front. The numbers live in `config.sh`
(`FIG_DATASET`, `FIG_REPRODUCTION`) and they name the output files and open each
generated legend, so renumbering the paper is one edit in one file.

**Every figure script generates its own legend**, from the same variables that
drew the panels, to `<prefix>_legend.txt`. `./run.sh legends` concatenates those
in figure order into `$FIGURE_LEGENDS` (`figures_legends.txt` at the repo root).
The point of generating rather than writing them is that a number cannot disagree
between a figure and its legend: change a threshold, re-run the figure, re-run
`legends`, and the prose follows. Do not hand-edit `figures_legends.txt`.

It lives outside `results/` deliberately — `results/` is gitignored as
regenerable output, and legends are manuscript text.

**`results/figures/` is the one tracked part of `results/`.** The seven figures'
PNGs and PDFs are committed, along with each figure's generated legend and stats
table, so they render on GitHub and a reader can check any number without a
re-run. SVG is excluded on purpose: it regenerates in seconds and figure 5's
alone is 1.7 MB of XML that would churn on every redraw.

**All three devices are cairo** (`png(type="cairo")`, `cairo_pdf`, `svglite`).
The base `pdf()` device transliterates UTF-8, which turned *Muñoz* into *Munoz*
and em dashes into hyphens.

### `figdataset` — the dataset, its QC and its quantification

```
./run.sh figdataset     # -> results/figures/figure<N>_dataset_qc.{png,pdf,svg}
```

| | |
|---|---|
| cost | ~30 s |
| env | `r_env` (ggplot2, patchwork, scico, svglite) |
| inputs | both sample sheets, both `multiqc_salmon.txt`, both TPM matrices, both VST exports, both node-metric tables |

The paper's opening figure. Every later result rests on two things a reader
cannot check from the text: that the two reused designs can be compared at all,
and that re-quantifying somebody else's libraries actually worked.

**A — the design.** Libraries per study × genotype × nitrogen × leaf segment.
Note `facet_wrap`, not `facet_grid`: a grid draws a cell for every study ×
genotype pair including the six that do not exist. Fill carries the study, not
the count — with every cell at 3 a continuous colour bar is decoration.

> **The sugarcane sheet's `tissue` column is not the design.** The study sampled
> **four** leaf segments — `base0`, `base`, `mid`, `tip`, 12 libraries each,
> carried in the library names as `B0/B/M/P`. `tissue` collapses `base0` and
> `base` into one "Leaf Base" of 24, calls `mid` "Leaf" and `tip` "Leaf Apex".
> Use the **`segment`** column, which carries the four real levels. `tissue` is
> left in place so the original sheet's provenance is intact; nothing should read
> it as the segment factor.

**B — depth against mapping rate**, one point per library from the salmon logs.
Plotted against each other rather than as two distributions because the failure
worth seeing is a library that is both shallow and poorly mapped, which is a
position on this plane and not a value on either axis.

**C — the gene funnel**: annotated and quantified → surviving the CV filter →
nodes in the network, with the fraction retained. A reader who sees only
"103,336 nodes" cannot tell whether that is most of the annotation or a tenth.

**D — PCA per study** on the 2,000 most variable genes (`PCA_NTOP`).
**What each axis tracks is measured, not assumed** — the script regresses each
component on each design factor and the legend quotes the R². That check earned
its place twice: the first draft of the legend asserted that sugarcane's leading
axis was leaf segment (PC1 is genotype at R² = 0.998 in *both* studies, and
nitrogen loads on neither first component), and the second measured segment on
the collapsed `tissue` column. On the four real segments, leaf segment explains
**0.802** of sugarcane's PC2 against 0.450 on the collapsed three — which is also
the quantitative case for `HEATMAP_GROUP_BY="segment"`.

### `figclustering` — MCL against the SBM

```
./run.sh figclustering [study]   # -> figure<N>_clustering_<study>.{png,pdf}
```

Not a paper figure — evidence for a decision that has not been made, numbered
past the seven so it cannot be mistaken for one. Both partitions carried through
the identical downstream analysis: **A** module-size CCDFs with the adjusted Rand
index, **B** eigengene coherence, **C** the module–trait test at each threshold,
**D** what share of each partition's responsive modules can be *named*.

Panel D is a share, not a count, on purpose: counts would put 465 modules beside
325 GO terms on one axis and invite reading the partition with more modules as
the better one.

### `figrepro` — reproducing each source study's own finding

```
./run.sh figrepro       # -> results/figures/figure<N>_reproduction.{png,pdf,svg}
./run.sh legends        # -> figures_legends.txt
```

| | |
|---|---|
| cost | ~40 s |
| env | `r_env` (ComplexHeatmap, circlize, scico, svglite) |
| inputs | the module20 and myb61 readouts under `results/readouts/`, plus both salmon TPM matrices |

**Why this figure exists.** Everything else in this project is a comparison
between two studies that were never designed to be compared. That comparison is
only worth reading if the re-quantification first reproduces what each study
found *on its own data* — so before any cross-species result, the paper shows
Muñoz's Module 20 and Kiet's MYB61 recovered through this pipeline, from these
references.

**Panel A** — Module 20, sugarcane, 33 mapped genes x 48 libraries. Columns split
by **nitrogen only**, genotype as an annotation bar. That is the layout that
separates "tracks nitrogen" from "tracks genotype" by eye: a genotype-driven set
splits into two sub-clusters *inside* each nitrogen block rather than differing
*between* blocks.

**Panel B** — MYB61, purple, 16 confirmed copies x 18 libraries. Columns split by
**genotype**, ordered 0N -> 2N -> 6N inside each, because Kiet's claim is that the
two genotypes behave differently. The copies at the published id `09G0002230`
ride along as a labelled **negative control** block — that id does not resolve to
a MYB (see `11_readouts/myb61/README.md`) and those genes run one to two orders
of magnitude hotter than any true MYB61 copy.

**The panels are not symmetric and the legend does not pretend they are.**
Muñoz's result reproduces outright. Kiet's reproduces in *shape* — a significant
non-monotonic, genotype-restricted response — but with the sign **inverted**, and
at expression levels where most copies sit under 1 TPM in leaf. The generated
legend says so, and the numbers are also written to
`figure2_reproduction_stats.tsv` for checking without re-reading the prose.

**Where OrthoFinder does and does not enter.** It is *not* how panel A was built:
Module-20 members are de novo assembly ORFs, mapped into the R570 proteome by
reciprocal DIAMOND blastp through Arabidopsis. It enters panel B as one of four
independent lines of evidence identifying the purple MYB61 copies (alongside the
reciprocal best-hit search, an independent Myb-domain call and a MAFFT/FastTree
phylogeny), and it is the same orthology — OrthoFinder v3.1.3, `-S diamond -M msa
-A famsa -T fasttree`, 94,273 orthogroups over both proteomes — that underlies
every cross-species comparison in this work. The legend states it in those terms
rather than crediting it with the whole figure.

**Rows are per-gene z-scores** in both panels, because the question is the shape
of the response and these genes span three orders of magnitude in absolute
expression. Mean `log2(TPM+1)` is a row-annotation barplot instead, so a row that
is structured but silent cannot be mistaken for a result — the failure mode a
z-score heatmap has on its own. Panel B clips its ramp at the **95th** percentile
of |z| rather than the 98th the other figures use: most MYB61 copies are near
zero in leaf, so their z-scores are spiky by construction and a 98% clip lets
three cells set the whole scale.

**Expression is TPM**, not VST, from the same salmon quantification the networks
were built from — these are per-study descriptive panels reproducing per-study
claims, and TPM is what both source papers report.

**Titles are drawn by hand** into their own layout rows rather than passed to
`draw(column_title=)`. A ComplexHeatmap column title is centred on the heatmap
*body* and is not wrapped, so a two-sentence caption runs off both edges. The two
heatmaps also share neither rows nor columns — different species, different
libraries — so they go into separate viewports of one grid layout rather than
being combined with `%v%`, which would force them onto one axis.

### `figtopology` — what the two networks are shaped like

```
./run.sh figtopology    # -> results/figures/figure<N>_topology.{png,pdf,svg}
```

| | |
|---|---|
| cost | ~20 s |
| env | `r_env` (ggplot2, patchwork, scico, svglite, scales) |
| inputs | node metrics, global metrics, MCL module summaries, the merge step's `edges.summary.json` |

**No edge table is read.** The two edge files are 76 M and 706 M rows; every
number here comes from the small per-study summaries.

Its job is to stop a reader importing intuitions from sparse biological
networks. These are dense thresholded correlation graphs — mean degree ~1,500
and ~8,300, edge density 1.4% and 4.8%, and purple is a **single** connected
component. Every later claim about hubs, modules or centrality has to be read
against that.

**A and D are complementary CDFs**, not frequency histograms. At this size the
upper tail of a histogram holds one node per bin and reads as noise; a CCDF is
monotone by construction so the tail is legible. Evaluated on a 300-point
log-spaced grid (`TOPO_GRID`) rather than at every unique value, which keeps the
SVG small and draws the same curve.

**B plots mean weight PER EDGE, not node strength.** `strength()` is the *sum* of
a node's edge weights, so it tracks degree trivially (Spearman ρ = 0.98 / 0.99)
and plotting it answers nothing. Dividing by degree asks the real question — is a
high-degree gene well connected, or connected to many weak partners? It uses a
**coarser 60-point grid** than A and D: a CCDF stays smooth at 300 points but a
per-bin *mean* does not, and the narrow mid-range bins spiked on a handful of
nodes each.

> The answer turned out to be neither of the two shapes I first wrote into the
> legend. Per-edge weight is flat across the first nine degree deciles
> (0.174–0.203 in sugarcane) and rises only in the tenth (0.321). A two-decile
> summary reads as a gradient and the raw ρ reads as a trend; both are artifacts
> of the top decile. The script now computes the whole decile profile and the
> legend quotes its flat span.

### `figconservation` — what transfers between the two networks

```
./run.sh figconservation  # -> results/figures/figure<N>_conservation.{png,pdf,svg}
```

| | |
|---|---|
| cost | ~15 s |
| env | `r_env` (ggplot2, patchwork, scico, svglite, scales) |
| inputs | `conservation_summary_*_FULL.tsv`, `conservation_null_*.tsv`, `conserved_correlated_summary_<selection>.tsv`, both node-metric tables |

**No edge table is read.** `conserved_edges_*_FULL.tsv` are 5 GB and 39 GB; every
number comes from the summary and null tables beside them.

The figure has to carry two findings that pull in opposite directions and are
both true: edge conservation is real and ~2.5× above chance in both directions,
and the nitrogen response nevertheless fails to transfer at all.

**A — observed against the permutation null, per direction.** The raw rates
differ ~7× between directions and that difference is *meaningless* on its own —
purple has 9.3× more edges, so a sugarcane edge has far more opportunity. Error
bars are the full null range over `NULL_REPS` permutations. The axis carries the
percentage, the bars carry the **absolute edge counts**, and the axis label
carries the denominator, so all three are readable without arithmetic. The null
bar's count is the number its rate implies over the same full edge set.

**B — what the conserved genes are *for*.** Panel A is a structural statement;
this asks whether the structure that transfers does anything recognisable. One
topGO enrichment per species over the genes on ≥1 conserved edge, against that
network's own GO-annotated nodes, showing the terms enriched in **both**. Ranked
by the *worse* of the two p-values, so the panel shows agreement rather than
terms one species drives. 69 of 192 BP terms are shared (Jaccard 0.36), and two
of them speak directly to the trait — *glutamate biosynthetic process* and the
*ammonia assimilation cycle*.

> `CONS_GO_NTERMS=12`, not 10, because those two nitrogen terms rank 11th and
> 12th by agreement and drop off the panel at 10.

> **topGO truncates term names** at 40 characters with an ellipsis, so 18 of the
> 69 shared BP terms arrive cut in the CSV. `GO.db` has the full names but lives
> in `topGO_env`, not the plotting env, so `run.sh` dumps an id → term cache once
> (`conservation/go_term_names.tsv`) and the figure expands from it, falling back
> to the truncated string if an id is missing.

The fold-over-null-by-layer panel this replaced is not lost: panel A's legend
carries the overall folds (2.57× and 2.49×) and panel C carries the per-layer
comparison. The empirical p is 1/21 for every layer, which is the floor set by 20
permutations, not a precise p — the legend says so and quotes z instead.

**C — edge composition by layer**, stacked as a share per network with the
counts on the bars and the two minority classes under the axis. Panel C used to
plot each non-linear layer relative to Pearson-only under two normalisations;
the composition is the more direct statement of the same case, and the
fold-over-null numbers that carried the argument now live in the legend:
MI-only conserves at 2.58× and 2.27× against Pearson-only's 2.57× and 2.46× —
parity in one direction, **below** Pearson in the other — while `both`-estimator
edges reach 2.76× and 2.95×.

> This duplicates Figure 3's panel C. Deliberate: in Figure 3 it describes the
> networks, here it sits against the conservation result, which is where the
> case on the MI layer actually gets made.

> Small-slice labels do not fit inside a 1–7% band and collide however they are
> justified, so the minority classes are labelled under the axis as percentages;
> counts stay in the legend and `_stats.tsv`.

**D — the funnel**, log scale: nodes → on a conserved edge → nitrogen-responsive
→ responsive on both sides. It closes at purple's 32 responsive genes, not at
orthology or conservation.

> Note the panel-C wording. The first draft said MI-only edges "look better
> conserved in the raw rate" in both directions; in purple → sugarcane they are
> already below parity raw (0.94×). The legend now states each direction
> separately, which makes the conclusion stronger, not weaker: there is no
> direction in which MI-only edges are better conserved once chance is accounted
> for.

### `figmodules` — the nitrogen response at module level

```
./run.sh figmodules       # -> results/figures/figure<N>_modules.{png,pdf,svg}
```

| | |
|---|---|
| cost | ~20 s |
| env | `r_env` (ggplot2, patchwork, scico, svglite, scales) |
| inputs | both `module_profile_*`, both `module_trait_*.null.tsv`, both eigengene matrices, both sample sheets |

**A — selection**, |rho| against −log10(padj) with both thresholds drawn. Its job
is to show *which threshold binds*, and it is a different one in each study: at
n = 48 an |rho| of 0.6 already implies p ≈ 6e-06 so BH never binds for sugarcane
and the effect-size floor is the only active constraint; at n = 18 it is the
reverse (274 clear |rho|, 79 survive BH).

**B — Spearman against Pearson** on the same eigengenes at identical thresholds:
408 → 465 and **38 → 79**.

**C and D — every responsive module's eigengene**, one panel per species. Rows
are modules split by the sign of rho and ordered by |rho| within each block;
columns are libraries split by nitrogen status and ordered by genotype, then leaf
segment, then library. Both splits label themselves, so no annotation strip is
needed. **One colour scale across both panels** (clipped at the 98th percentile
of |z| over both sets together), with the colourbar drawn once — a cell in C then
means what a cell in D means.

Drawn with `geom_raster`, not ComplexHeatmap: these compose with two ggplots in
one patchwork, and a heatmap object would have to be grabbed into a grob to get
there. Rows are ordered by rho rather than clustered — deterministic, and the
ordering is itself interpretable.

A mean profile would show the same two shapes in four points. What the heatmaps
add is whether the *set* is coherent, and it is: the nitrogen split accounts for
a median R² of 0.44/0.45 of each eigengene's variance in sugarcane's two blocks
and 0.56/0.62 in purple's.

**Two things live in the legend rather than costing a panel.** The permutation
null is two numbers (0/1000 and 2/1000), plus purple's heavy right tail — worst
of 1,000 shuffles: 244 responsive modules against 79 observed. And the TF
enrichment result (OR 2.65, p = 0.0016, concentrated in the modules that fall
with nitrogen at OR 3.34) is quoted in full there after its panel was replaced.

### `figmodulego` — module-level GO, at three grains

```
./run.sh modulego <study> [BP|MF|CC]   # produces the tables, incl. by direction
./run.sh figmodulego                   # -> figure<N>_module_go.{png,pdf,svg}
```

**The direction-level test is new** (added 2026-08-23 to `18_module_go.r`). There
were two grains — per module, and all responsive modules pooled — and the
question a reader asks first fell between them: *what do the modules that go UP
with nitrogen do, as a set, and is it different from the ones that go DOWN?*

It also breaks the constraint that dominates the per-module grain. A module needs
`MODULE_GO_MIN_ANNOTATED` annotated members for the test to be defined, and only
**49 of 465** sugarcane modules clear that (10 of 79 in purple), biased toward
the large ones. Pooled by direction every module contributes: the two sets carry
383 and 161 annotated genes drawn from 242 and 223 modules.

Same `topGOdata` object, background, statistic and raw-p threshold as the
per-module run, so the grains are directly comparable. Writes
`module_GO_<ONT>_<study>_bydirection.tsv` and a `_comparison.tsv` classifying
each term as up-only / down-only / shared.

**The result is worth the panel.** In sugarcane, 74 terms are enriched only among
the modules that rise, 23 only among those that fall, and **6 in both** — so the
split is not one programme cut in half. What rises is nitrogen assimilation and
what consumes it (nitrate assimilation, nitric oxide, proline/asparagine/
spermidine biosynthesis) plus a defence block; what falls is the low-nitrogen
carbon programme (flavonoid, raffinose-family oligosaccharide, triglyceride
biosynthesis, cold response). Neither of the other two grains recovers it.

**Three panels.** A is the annotation gate; **B and C are the same plot for the
two species**, side by side so the direction split can be compared. Their x
scales are free — sugarcane reaches −log10(p) ≈ 14 and purple ≈ 8, and a shared
axis would flatten purple into nothing.

The term-overlap counts and the per-module recurrence that earlier drafts drew as
panels are in the legend instead: the first is six numbers, and the second is a
weak signal drawn from the ~11% of modules that are individually testable.

`MODULE_GO_FIG_MAIN` still names which species leads (sugarcane), and the other
follows in C.

### `figmodule20` — Muñoz's Module 20 in purple

```
./run.sh figmodule20      # -> results/figures/figure<N>_module20.{png,pdf,svg}
```

| | |
|---|---|
| cost | ~15 s |
| env | `r_env` |
| inputs | the module20 readout tables, both TPM matrices, both sample sheets |

The comparative question about Module 20, and the answer is not "less of the
same". **A** the module as Muñoz report it, reproduced — 33 mapped sugarcane
genes × 48 libraries, small, because it is the premise rather than the result.
**B** the 12 published members resolve to three Arabidopsis anchors and the
mapped genes are haplotype *copies* of those loci (33 sugarcane genes from 10
members, 27 purple from 6); only the AtMYB59 anchor carries a MYB call.
**C** all 27 purple copies as a **per-gene z-score**, split by **genotype** as
well as nitrogen, with only the responding copy labelled. **D**
`Soffic.09G0001580-9H` on its own.

> Two deliberate choices in C, both reversals of an earlier draft. The **z-score**
> is what makes the responding copy legible: on an absolute scale it is a bright
> row among dim ones and the *shape* of its response cannot be seen. Splitting by
> **genotype** is what makes the dominant structure visible — in purple these
> genes vary far more between the two species than across nitrogen, and a panel
> split by nitrogen alone hides that and invites the reader to credit the
> treatment.

> **The finding.** That one copy responds in a **U centred on the control** —
> 124.1 TPM under starvation, 10.8 at the 2 N control, 47.4 under excess in
> 51NG3. Since 2 N is the control and both extremes are stresses, it is induced
> by nitrogen stress *in either direction*, which is exactly what a monotone rank
> correlation cannot see. Sugarcane's nine copies at the same anchor are all
> monotonically repressed by nitrogen — but sugarcane has no control level, so it
> can only observe one arm. The two studies are consistent, not contradictory.
>
> It rests on one expressed copy at U-contrast p = 0.0078 that does **not** clear
> BH over the 27 purple Module-20 genes. The legend calls it a candidate.

The focus gene and anchor are `M20_FOCUS_GENE` / `M20_FOCUS_ANCHOR` in
`config.sh`, not hardcoded in the script.
