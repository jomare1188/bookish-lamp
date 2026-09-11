# What in `results/` no longer matches, and why

> Tracked copy of `results/STALE.md`. The results tree is gitignored, so the file
> that records its inconsistencies has to live here too or it disappears with the
> tree it describes.

Written 2026-09-10, when the module analysis moved to the **unpruned Pearson-only**
networks at each species' own modularity optimum (sugarcane `-I 1.5`, purple
`-I 3.5`). Everything displaced is in `_pre_pearson_20260910/`.

## Current — rebuilt on the new networks and clustering

| file | |
|---|---|
| `<study>/network_<study>_node_metrics.tsv` | Pearson-only gene universe (101,990 / 170,135). `transitivity` is NA by design — no consumer reads it |
| `<study>/mcl_<study>_membership.tsv`, `_module_summary.tsv` | the adopted partition at the chosen inflation |
| `<study>/modules/` | eigengenes (3,627 / 7,493) |
| `<study>/module_trait_<study>.tsv` | **blocked** partial correlation, marginal carried beside it |
| `<study>/module_profile_<study>.tsv`, `heatmaps/`, `module_summary_<study>.*` | |
| `figures/figure5_modules.*` | |

## STALE — do not quote against the new modules

These were computed on the **merged (Pearson + MI) network at `-I 2`** and have not
been re-run. They are not wrong; they describe a different graph and a different
clustering, so joining them to anything above compares two analyses.

| file / stage | what it still describes |
|---|---|
| `<study>/network_<study>_edges.tsv` | the MERGED network. The Pearson-only edge tables were deleted (70 GB whose only consumer was `mcxload`); the `.mci` in `/dados04/jorge/tmp/mcl_work_cluster/` is the live Pearson-only graph |
| `<study>/network_<study>_global_metrics.tsv`, degree/strength plots | merged network |
| `figures/figure*_topology*` (`22_fig_topology.r`) | merged network |
| `conservation_*` (`06_conservation_join.r`, `13_conservation_null.r`) | merged networks, both directions |
| `conserved_*` (`08_conserved_cor_genes.r`) | merged networks + the marginal gene-trait rule |
| `<study>/gene_trait_*` (`07_gene_trait_cor.r`) | **marginal** gene-level test. The blocked gene-level version lives on branch `blocked-gene-trait` (`31_gene_trait_blocked.r`) and is not merged here |
| — | **BP is current** (run 2026-09-11 on the new clustering). The CC and MF results in that directory were from 2026-08-23 and describe the OLD clustering; they have been moved to `_pre_pearson_20260910/<study>/module_go_OLD_CLUSTERING/` so one directory does not hold two analyses |

## A loss to record

The backup taken on 2026-09-10 listed `mcl_*`, `modules/`, `network_*_node_metrics.tsv`,
`module_trait_*` and `figures/` — **it did not include `module_go/`**. Running BP on
2026-09-11 therefore overwrote the previous BP results without a copy. They are
reproducible rather than gone: the old membership, node metrics and module-trait
table are all in the backup, so `15_module_profile.r` then `18_module_go.r` pointed
at `_pre_pearson_20260910/<study>/` regenerates them. CC and MF were untouched by
that run and are preserved.

## The trap this file exists to prevent

Module names are **positional** — `Module_%03d`, largest first. `Module_001` under
the old clustering held 19,604 sugarcane genes; under the new one it holds 23,439.
Any join on module name across the two trees runs cleanly and compares unrelated
gene sets. This already happened once during this rebuild and was caught only
because the result looked wrong.
