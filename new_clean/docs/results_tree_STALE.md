# What in `results/` no longer matches, and why

> Tracked copy of `results/STALE.md`. The results tree is gitignored, so the file
> that records its inconsistencies has to live here too or it disappears with the
> tree it describes.

Written 2026-09-10, when the module analysis moved to the **unpruned Pearson-only**
networks at each species' own modularity optimum (sugarcane `-I 1.5`, purple
`-I 3.5`). Updated 2026-09-12, when the **gene** analysis followed it there.
Everything displaced is in `_pre_pearson_20260910/`.

## Current — rebuilt on the new networks and clustering

| file | |
|---|---|
| `<study>/network_<study>_node_metrics.tsv` | Pearson-only gene universe (101,990 / 170,135). `transitivity` is NA by design — no consumer reads it |
| `<study>/mcl_<study>_membership.tsv`, `_module_summary.tsv` | the adopted partition at the chosen inflation |
| `<study>/modules/` | eigengenes (3,627 / 7,493) |
| `<study>/module_trait_<study>.tsv` | **blocked** partial correlation, marginal carried beside it |
| `<study>/gene_trait_blocked_<study>.tsv` | **blocked** gene-level nitrogen test on the Pearson-only node universe (101,990 / 170,135), marginal fit carried beside it on the same genes |
| `purple/gene_trait_ushape_purple.tsv` | the blocked quadratic contrast, on the same universe so the two test families share a denominator |
| `conservation/conserved_correlated_*_blocked_nodes*.tsv` | the conserved nitrogen response at **node** level, blocked rule, with the ortholog-shuffle null |
| `<study>/module_go/` | **BP, MF and CC**, all run 2026-09-13 on the adopted single-source annotation (full InterProScan, 17 DBs) over the blocked responsive modules |
| `figures/figure6_module_go.*` | rebuilt 2026-09-13; panel A's gate numbers read 588/479/425 and 182/118/109 straight from the current tables |
| `<study>/module_profile_<study>.tsv`, `heatmaps/`, `module_summary_<study>.*` | |
| `figures/figure5_modules.*` | |

## STALE — do not quote against the new modules

These were computed on the **merged (Pearson + MI) network at `-I 2`** and have not
been re-run. They are not wrong; they describe a different graph and a different
clustering, so joining them to anything above compares two analyses.

| file / stage | what it still describes |
|---|---|
| `<study>/network_<study>_edges.tsv` | the MERGED network. The Pearson-only edge tables were deleted (70 GB whose only consumer was `mcxload`); the `.mci` in `/dados04/jorge/tmp/mcl_work_cluster/` is the live Pearson-only graph |
| — | `network_<study>_global_metrics.tsv` is **current** (recomputed 2026-09-12 from the Pearson-only matrix: 958 and 43 components). The per-study degree/strength plots are still from the merged build |
| — | **figure 3 is current**: rebuilt 2026-09-12 on the Pearson-only graphs and the new clustering, with panels B (mean weight vs degree) and C (edge composition by layer) removed |
| `conservation_*` (`06_conservation_join.r`, `13_conservation_null.r`) | merged networks, both directions |
| `<study>/gene_trait_correlations_<study>.tsv`, `selected_genes_*` (`07_gene_trait_cor.r`) | the **marginal** gene-level test, over the merged graph's conserved-edge gene set. Superseded at gene level by `gene_trait_blocked_<study>.tsv`, which carries the marginal fit as a column on the current universe |
| `<study>/gene_trait_mi_<study>.tsv` (`12_gene_trait_mi.py`) | the MI statistic. There is no MI layer in the current networks, so nothing selects on it any more |
| `conserved_*` **edge level** (`08_conserved_cor_genes.r`) | merged networks + the marginal rule. `61_conserved_blocked_nodes.r` replaces its NODE level only — the edge level needs `conserved_edges_*_FULL.tsv`, which describes the merged graph, and has not been rebuilt |

## GO: one source, adopted 2026-09-12, closed 2026-09-13

`module_go/` is built on GO derived from a **full local InterProScan 5.78 over all 17
member databases**, via the GO Consortium's pinned interpro2go/pfam2go and normalised
to most-specific terms. `annotation/<study>/gene2go_<study>.tsv` is the only table any
stage reads; every other source was scored and rejected, and their tables now carry a
`.notused` suffix. `annotation/README.md` records which is live and why the others are
not. PANNZER2's output was deleted (6.1 GB, and its sugarcane half never finished).

Network coverage 62.9% (sugarcane) and 64.4% (purple), up from 8.0% / 7.2% under the
original eggNOG GO column — that column used emapper's default `--go_evidence
non-electronic`, which excludes every IEA term, and for grass proteins almost all GO
is IEA.

STILL STALE, and untouched by that work: **`09_go_enrichment.r` and `10_go_semantic.r`
outputs**, which run on the old 7.7% annotation and its background. They are
conserved-set and semantic stages, not module stages. `emapper.annotations` itself is
untouched and is still the source for KEGG and Preferred_name; only its GO column was
rejected. PFAM-based work (`37_cluster_homogeneity.r`) is unaffected — eggNOG's PFAMs
column is 87.6% filled and was never the problem.

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
