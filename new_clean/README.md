# new_clean — comparative co-expression networks, linear + non-linear

One coherent pipeline for the two sugarcane RNA-seq studies, from DESeq2 objects
to the H1 readouts. It builds **one network per study** whose edges come from two
layers — Pearson correlation (linear) and Kraskov–Stögbauer–Grassberger mutual
information (non-linear) — thresholded at the *same per-edge false-positive rate*
so their union is legitimate, with every edge labelled by which layer found it.

| | study | reference | libraries |
|---|---|---|---|
| `sugarcane` | Muñoz-Perez et al. 2025 | R570 | 48 leaf |
| `purple` | Ta Quang Kiet et al. 2025 | LA purple | 18 leaf |

Nothing outside this directory is read for results or written to. The previous
tree is untouched.

## Why it exists

Pearson sees only linear co-variation, so saturating and non-monotone regulatory
relationships — the expected shape of a nitrogen dose response — are structurally
invisible to it. The MI layer finds those. `source ∈ {pearson, mi, both}` on
every edge is what makes the comparison possible:

- do the modules recruit MI-only edges preferentially, or evenly?
- do Module 20 or the MYB61 copies gain neighbours Pearson could not see?
- **is non-linear co-expression conserved between the two species at the same
  rate as linear co-expression?** (`06_conservation_join.r` reports this directly)

For what changed relative to the old tree and why, read
[docs/decisions.md](docs/decisions.md). It documents, among other things, that
the old networks contained ~47,000 edges below their own stated threshold.

## Running it

Every stage is one command. **No script is ever edited between runs** — species,
direction and ontology are arguments, and everything else lives in `config.sh`.

```bash
cd /dados04/jorge/comparative_saccharum/new_clean

./run.sh export sugarcane          # dds -> VST binary                  minutes
./run.sh export purple
./run.sh validate                  # must print ALL CHECKS PASSED        ~2 min

./run.sh build sugarcane           # both layers + merge                   3.1 h
./run.sh build purple              #                                       2.6 h

./run.sh stats    sugarcane        # node + global metrics                 3 min
./run.sh stats    purple           #                                      65 min
./run.sh mcl      sugarcane        # modules                              13 min
./run.sh mcl      purple           #                                        2 h
./run.sh conserve sugarcane_to_purple    # cheaper direction — run first  11 min
./run.sh conserve purple_to_sugarcane    # 9x more outer edges           ~75 min
./run.sh conservenull sugarcane_to_purple ; ./run.sh conservenull purple_to_sugarcane
./run.sh trait    sugarcane ; ./run.sh trait purple
./run.sh traitmi  sugarcane ; ./run.sh traitmi purple   # gene-trait MI  <1 min
./run.sh conscor
./run.sh go BP ; ./run.sh go MF ; ./run.sh go CC
./run.sh gosem

# module level, per study                                              ~5 min
./run.sh eigengene <study>         # PC1 per MCL module                  1 min
./run.sh moduletrait <study>       # module response: Spearman rho      10 s
./run.sh moduleprofile <study>     # + TF hypergeometric per module
./run.sh modulego <study>          # topGO BP per responsive module      45 s
./run.sh moduleheatmap <study>     # per-module gene heatmaps            1 min
./run.sh modulesummary <study>     # one figure: all responsive modules

# paper figures
./run.sh figure1                   # Fig 1: both source studies reproduced  40 s
./run.sh legends                   # -> figures_legends.txt (repo root)
```

Figures carry a panel letter and the labels the data needs, and nothing else.
Every description lives in the legend, which each figure script GENERATES from
the variables that drew it — so a number cannot disagree between a figure and its
legend. Never hand-edit `figures_legends.txt`; re-run the figure, then `legends`.

`build` is `export` + `network <study> pearson` + `network <study> ksg` +
`merge`; run those individually if you want to watch them.

**The whole thing is an overnight job, not a weekend one** — about 5.5 h to both
networks and ~9 h through clustering. Times above are measured, not estimated;
see [docs/methods.md](docs/methods.md) for the per-stage breakdown. The one long
stage is sugarcane's KSG sweep at 2.8 h, and it is long because of n = 48 (a
2e8-permutation null to resolve p = 9e-12), not because of gene count — purple's
KSG layer, on a 9x denser network, takes 45 min.

Everything is logged to `logs/<stage>_<arg>_<timestamp>.log` as well as stdout.
The GPU sweeps take `--resume` and checkpoint per tile, so a killed run picks up
where it left off.

**Run sugarcane first at every stage.** It is ~9× smaller and catches path and
schema errors in minutes instead of hours.

## The shape of it

```
dds ──► 01_export_vst.r ──► <study>.f32 + .genes.txt + .meta.json
                                  │            one matrix, one gene set,
                    ┌─────────────┴─────────────┐   both layers
                    ▼                           ▼
        02_network_engine.py          02_network_engine.py
          --estimator pearson           --estimator ksg
          (analytic t null)             (permutation null + GPD tail)
                    │                           │
        <study>_pearson.edgelist.tsv   <study>_ksg.edgelist.tsv
                    └─────────────┬─────────────┘
                                  ▼
                        03_merge_layers.py
                                  ▼
                 results/<study>/network_<study>_edges.tsv
                                  │
   04_stats ─► 05_mcl ─► 06_conservation ─► 07_trait ─► 08_conscor ─► 09/10_GO
                     │            │
                     │      13_conservation_null
                     ▼
         14_eigengene ─► 19_spearman ─► 15_profile ─┬─► 16/17 figures
                                                    └─► 18_module_go
                                  │
                            11_readouts/  (TFs, MYB61, Module 20)
```

## The network file

`results/<study>/network_<study>_edges.tsv` is the only edge table. There is
exactly one per study, which is the point — the old tree had a "filtered" and an
"augmented" copy and half the scripts read the stale one.

```
gene1  gene2  stat  pval  padj  weight  pearson_r  ksg  source
```

| column | meaning |
|---|---|
| `stat` | strength on the correlation scale: \|r\| for Pearson edges, `r_eq` for MI-only, the larger of the two for `both` |
| `pval`/`padj` | the more significant of the two tests |
| `weight` | `0.01 + (stat − 0.8)/(0.9999 − 0.8) × 0.99` — what the clustering consumes |
| `pearson_r` | the true **signed** r, `NA` for MI-only edges |
| `ksg` | mutual information in nats, `NA` for Pearson-only edges |
| `source` | `pearson` \| `mi` \| `both` |

`r_eq` converts nats to the correlation scale via the Gaussian identity
`I = −½ln(1−r²)`, after subtracting the KSG small-sample bias. For a
**non-monotone** edge it means "as strong as a linear r of this size", *not*
"its Pearson r is this" — which is why `source` must always be read alongside it.

## Layout

```
config.sh          every path and parameter; the single source of truth
run.sh             the dispatcher
scripts/
  01_export_vst.r        dds -> flat float32 matrix (+ id normalisation)
  02_network_engine.py   GPU all-pairs: pearson | ksg | gcmi | xi
  03_merge_layers.py     the two layers -> the network
  04_network_stats.r     node + global topology
  05_mcl_clustering.r    MCL modules
  06_conservation_join.r cross-species conserved edges, broken down by layer
  07_gene_trait_cor.r    per-gene expression vs trait
  08_conserved_cor_genes.r  node-level conservation of the N response
  09_go_enrichment.r     topGO, network nodes as background
  10_go_semantic.r       GO semantic-similarity clustering
  12_gene_trait_mi.py    gene-vs-trait MI (Ross, discrete trait)
  13_conservation_null.r permutation null for edge conservation
  14_module_eigengene.r  PC1 per MCL module, in the VST export's own format
  15_module_profile.r    module response + coherence + TF hypergeometric
  16_module_heatmaps.r   gene-level heatmap per responsive module
  17_module_summary.r    every responsive module's eigengene in one figure
  18_module_go.r         topGO per responsive module (topGO_env)
  19_module_trait_spearman.r  module eigengene vs trait, Spearman only
  20_figure1_reproduction.r   paper Fig 1: Module 20 + MYB61, png/pdf/svg
  11_readouts/           get_tfs, myb61, module20 (H1); cached sequence work is
                         read from the original GET_TFS tree, output lands here
  validate.py            engine correctness suite
  lib/common.R           config, logging, id handling, VST + orthogroup readers
results/{sugarcane,purple,conservation,figures}/  gitignored
logs/                                       gitignored
docs/
  decisions.md     why the pipeline is the way it is, with evidence
  thresholds.md    the MI threshold argument in full
  methods.md       per stage: inputs, outputs, cost, environment
  results.md       the numbers
```

## Environments

`cor_env` has DESeq2 but **not** igraph or hexbin, so it cannot run the network
scripts. `run.sh` picks the right one per stage; these are for reference.

| use | interpreter |
|---|---|
| VST export (DESeq2) | `/home/genomics/miniconda3/envs/cor_env/bin/Rscript` |
| network R stages | `/home/genomics/miniconda3/envs/r_net_env/bin/Rscript` |
| GPU engine (torch + cu130) | `/home/genomics/miniconda3/envs/docling/bin/python` |
| GO enrichment | `topGO_env`, then `r_clusterprofiler` |
| readout plots | `r_env`; `comparative_network` for the MYB61 tree (`ape`) |

`mcl` does not need to be on `PATH` — it lives in `r_net_env/bin` and `run.sh`
prepends that directory from `MCL_BIN`. Nor is `cogeqc` needed: orthogroups are
parsed by `lib/common.R:read_orthogroups_tsv()`, verified against
`cogeqc::read_orthogroups` on this project's own file. Both were environment
assumptions inherited from the old tree that failed at run time here.
