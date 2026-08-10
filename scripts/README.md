# Scripts — comparative gene co-expression network pipeline

This folder contains the full analysis pipeline used to build, compare and
mine gene co-expression networks for the two studies described in the main
repository README: **sugarcane** (Muñoz-Perez et al. 2025, "Network 1") and
**purple** (Ta Quang Kiet et al. 2025, *S. officinarum* / *S. robustum*,
"Network 2"). Most scripts process one study at a time and are configured by
editing the `CONFIGURATION` / `INPUT FILES` block at the top (commented-out
blocks are the alternate study — swap which one is active and re-run).

`old/` holds superseded versions of earlier scripts, kept for reference only;
they are not part of the active pipeline.

## Pipeline overview

```
1. run_download.sh / run_rnseq.sh   raw reads -> quantification -> DESeq2 dds (+VST)
2. pca.r                            (QC, optional) sample PCA panel from VST
3. pearson_cor.r                    dds/VST -> all-vs-all gene correlation matrices
4. build_edgelist.r                 correlation matrices -> significant edge list
5. general_stats.r                  edge list -> filtered network + topology stats/plots
6. degree_distribution.r            (optional) combined degree-distribution panel
7. mcl_clustering.r                 filtered network -> MCL modules
8. eigengene.r                      modules -> module eigengenes (PC1 per module)
9. module_trait_cor.r               eigengenes vs traits -> selected modules + subnetwork
10. comparative_networks2.r         subnetworks -> cross-species conserved edges (module scope)
11. network_conservation_join.r     full networks -> cross-species conserved edges (chunked join)
    (full_matrix_edge_comparation.r is the matrix-based alternative to #11)
12. gene_trait_cor.r                per-gene expression vs traits, optionally restricted
                                    to genes on conserved edges -> nitrogen-responsive genes
13. conserved_cor_genes.r           node-level: are nitrogen-responsive genes conserved
                                    (ortholog also nitrogen-responsive)?
14. conserverd_edges_treatment.r    edge-level: STRICT/LOOSE nitrogen-conserved edges
15. top_GO_dev2.r                   (in old/, NOT currently used) DEG-based GO enrichment
16. top_GO_conserved.r              GO enrichment of the conserved-edge genes vs each
                                    network's own nodes -> shared/unique term comparison (H0)
17. GO_semanthinc_enrichment.r      semantic-similarity clustering of the shared conserved
                                    terms (H0 figure) -> "conserved topology = basic processes"
18. get_tfs/run_all.sh              reference proteome -> hmmsearch (TF HMM library) ->
                                    family assignment -> TFs that are network nodes (H1)
    get_tfs/05_tf_characterization.r  in-network TFs x (degree, N-correlation,
                                    conserved-edge) -> big table + MYB comparative figure
19. myb61/run_all.sh                AtMYB61 anchor -> RBH via sorghum/rice -> DIAMOND vs both
                                    proteomes -> RBH + domain-rule filter -> orthogroup bridge
                                    -> tree -> MYB61 copies and their network role (H1)
20. module20/run_all.sh             Munoz Module-20 proteins -> blastp + reciprocal-At filter
                                    -> both networks -> centrality vs MYB background, conserved
                                    edges, module recovery, design-aware N tests (H1)
21. mutual_information/run_all.sh   ALTERNATIVE to steps 3+4: GPU mutual-information networks
                                    (Kraskov kNN / Gaussian-copula / Chatterjee xi) on the same
                                    gene sets -> non-linear edges Pearson cannot see
```

Steps 10 and 11 are two independent, alternative routes to "which co-expression
edges are conserved between species" — 10 operates only on the nitrogen-trait
subnetworks selected in step 9 (cheap), 11 operates on the entire filtered
networks (expensive, needs chunking). Step 12 onward (the nitrogen-response /
conservation analysis) is normally run against the step-11 (full-network)
conserved-gene set.

---

## 1. `run_download.sh` / `run_rnseq.sh`

Nextflow (`nf-core/rnaseq`) launchers that turn raw reads into per-sample
quantifications and a DESeq2 object per study.

- `run_download.sh` — full alignment-based run (STAR + Salmon, `--aligner
  star_salmon`) against a genome + GFF/GTF, used for the sugarcane hybrid
  genome.
- `run_rnseq.sh` — pseudo-alignment run (`--skip_alignment --pseudo_aligner
  salmon`) against the R570 reference, used for the second study.

**Output:** salmon quantifications and a `deseq2.dds.RData` (DESeq2 object
with a VST assay) per study, e.g. `run1/salmon/deseq2_qc/deseq2.dds.RData`
and `china/run2_onlyL/salmon/deseq2_qc/deseq2.dds.RData`. All downstream `.r`
scripts read from these `dds` objects.

## 2. `pca.r` (QC, optional)

Loads both studies directly from their Salmon `quant.sf` files (via
`tximport`), filters genes by coefficient of variation (CV ≥ 15%), runs
`DESeq2::varianceStabilizingTransformation`, and builds a side-by-side PCA
panel colored by genotype. Purely a QC/exploratory step — not consumed by any
later script.

**Output:** `pca_panel/pca_panel_saccharum.png`.

## 3. `pearson_cor.r`

Loads a study's `dds`, filters genes with raw-count CV ≥ 15% (`MIN_CV`),
then computes the full all-vs-all Pearson correlation matrix on VST values
plus the matching p-values (analytic t-test), chunked and parallelized across
`N_CORES`/`CHUNK_ROWS` to stay within memory.

**Output (per study/label):** `matrix_<label>_pearson.tsv`,
`matrix_<label>_pvalues.tsv` (gene × gene matrices).

## 4. `build_edgelist.r`

Streams the pearson/pvalue matrices from step 3 (never loading them fully
into memory: 3-pass strategy — stream candidates to disk, BH-correct
p-values, re-extract surviving rows) and keeps gene pairs with `|r| >=
PEARSON_THRESHOLD` (0.7) and raw `p <= 0.05`, then applies BH FDR correction
and keeps `padj <= 0.05`.

**Output:** `edgelist_<label>_pearson.tsv` (`gene1, gene2, pearson, pval,
padj`).

## 5. `general_stats.r`

Reads the edge list from step 4, keeps edges with `|r|` in `[PEARSON_MIN,
PEARSON_MAX]` = `[0.8, 0.9999]` (the upper bound screens out suspicious
"perfect" correlations) and `padj <= 0.05`, min-max normalises the kept
`|r|` into an edge weight, builds the `igraph` network, and computes
global/per-node topology metrics and diagnostic plots.

**Output:** `network_<label>_filtered_edges.tsv` (the final co-expression
network — the main input to almost everything downstream),
`network_<label>_global_metrics.tsv`, `network_<label>_node_metrics.tsv`,
plus `network_<label>_degree_distribution.{pdf,png}`,
`network_<label>_strength_vs_degree.{pdf,png}`,
`network_<label>_transitivity_vs_degree.{pdf,png}`.

## 6. `degree_distribution.r` (optional)

Combines `network_sugarcane_node_metrics.tsv` and
`network_purple_node_metrics.tsv` (step 5 output for both studies) into one
two-panel degree-distribution figure for direct visual comparison.

**Output:** `degree_distribution_panel/degree_distribution_panel.png`.

## 7. `mcl_clustering.r`

Runs the native MCL binary (`system2("mcl", ...)`, multi-threaded via `-te`)
on each study's `network_<label>_filtered_edges.tsv`, converts the result
into a size-ranked module membership table (`Module_001` = largest), and
computes hub genes (top 10 per module by weighted strength) and module-size
plots.

**Output:** `mcl_<label>_membership.tsv` (gene → module, strength, degree),
`mcl_<label>_module_summary.tsv`, `mcl_<label>_hub_genes.tsv`,
`mcl_<label>_plot_data.tsv`, `mcl_<label>_module_sizes.{pdf,png}`, and a
combined `mcl_results_both_networks.rds`.

## 8. `eigengene.r`

For each study, loads the VST matrix and the module membership from step 7,
and computes each module's eigengene (PC1 across its genes, sign-oriented to
correlate positively with the module's mean expression) — the WGCNA-style
one-value-per-sample summary of a module.

**Output:** `eigengenes_<label>.tsv` (samples × modules),
`eigengenes_<label>_pc1_variance.tsv` (PC1 variance explained per module).

## 9. `module_trait_cor.r`

Correlates each module eigengene (step 8) against numerically-encoded traits
(genotype, nitrogen treatment) from the sample metadata, BH-corrects within
each trait, and selects modules with `|r| >= PEARSON_THR` (0.6) and `padj <=
PADJ_THR` (0.05) for the `treatment` trait. It then pulls out the induced
subnetwork — every edge from step 5 whose both endpoints fall in a selected
module.

**Output:** `module_trait_correlations_<label>.tsv`,
`selected_modules_genotype_<label>.tsv`,
`subnetwork_selected_modules_<label>.tsv` (edges restricted to
nitrogen-associated modules — the input to step 10).

## 10. `comparative_networks2.r`

Cross-species conservation at **module-subnetwork scale**. Maps edges of one
study's selected-module subnetwork (step 9) onto the other study through
OrthoFinder orthogroups (`Orthogroups.tsv`), using a sparse-matrix projection
(`O %*% A_b %*% t(O)`): an edge in network A is "conserved" if its two genes
have orthologs that are themselves connected in network B. Run in both
directions.

**Output:** `network_conservation/conserved_edges_sugarcane_to_purple.tsv`,
`conserved_edges_purple_to_sugarcane.tsv`, `conservation_summary.tsv`
(Jaccard-style conservation fraction per direction).

## 11. `full_matrix_edge_comparation.r` / `network_conservation_join.r`

Same conservation logic as step 10 but over the **entire filtered networks**
(step 5), not just the trait subnetworks — two alternative implementations:

- `full_matrix_edge_comparation.r` — same dense sparse-matrix-product
  approach as step 10. Documented as memory-risky at full scale (especially
  for the high-mean-degree purple network); intended to be tried first /
  monitored, not the default.
- `network_conservation_join.r` — scalable replacement: streams network A in
  chunks, projects each edge's endpoints through orthology, and looks them up
  in a keyed adjacency table of network B, so cost scales with ortholog
  fan-out rather than with B's edge density. This is the version meant to
  actually complete on the full purple/sugarcane networks.

**Output:** `network_conservation/conserved_edges_full_<a>_to_<b>.tsv` (matrix
version) or `conserved_edges_<a>_to_<b>_FULL.tsv` (chunked-join version),
`conserved_genes_<a>_FULL.txt` (genes sitting on ≥1 conserved edge — used as
the gene filter in step 12), and a `conservation_summary*.tsv`.

## 12. `gene_trait_cor.r`

Complements step 9: instead of correlating trait against the module
eigengene (which can dilute a real signal carried by only part of a
heterogeneous module), this correlates each gene's own VST expression
directly against the trait, using the same analytic-t-test method as step 3.
Normally restricted to the gene set from step 11
(`conserved_genes_<label>_FULL.txt`) so the resulting hits are, by
construction, both nitrogen-responsive and cross-species conserved.

**Output:** `gene_trait_correlations_<label>.tsv`,
`selected_genes_treatment_<label>.tsv` (feeds steps 13 and 14).

## 13. `conserved_cor_genes.r`

Node-level question: for a gene that is nitrogen-correlated in one network
(step 12 output), does it have an ortholog in the other network, and is that
ortholog *also* nitrogen-correlated? Classifies every correlated gene as
`no_ortholog`, `ortholog_not_correlated`, or `conserved_correlated`, and
checks sign concordance (same direction of response) for the hits.

**Output:** `network_conservation/conserved_correlated_ortholog_pairs.tsv`,
`sugarcane_correlated_conservation_status.tsv`,
`purple_correlated_conservation_status.tsv`,
`conserved_correlated_summary.tsv`.

## 14. `conserverd_edges_treatment.r`

Edge-level question: combines the conserved-edge table (step 11) with the
per-network nitrogen-responsive gene lists (step 12) to ask whether a
conserved edge's endpoints are trait-correlated on **both** sides — two
readings: **LOOSE** (both endpoints trait-correlated, and each has some
trait-correlated ortholog on the other side) and **STRICT** (LOOSE, plus that
specific ortholog pair must itself be an edge in the other network).

**Output:** `treatment_conserved_pairs_ALL_<a>_to_<b>.tsv`,
`treatment_conserved_edges_LOOSE_<a>.tsv`,
`treatment_conserved_edges_STRICT_<a>.tsv`,
`treatment_conserved_edge_pairs_STRICT_<a>_to_<b>.tsv`,
`treatment_conserved_summary_<a>_to_<b>.tsv`.

## 15. `top_GO_dev2.r` (in `old/`, not currently used)

Moved to `old/` — kept for reference, **not part of the active pipeline**.
The current analysis uses `top_GO_conserved.r` (step 16) instead.

GO term enrichment (`topGO`, classic Fisher test + BH correction) run per
"contrast" (a named gene list of interest, e.g. genes on conserved edges from
step 11) against an eggNOG-mapper functional annotation, with ID translation
between transcript/protein/locus namespaces (tx2gene / GTF CDS attributes) so
DEG IDs and annotation IDs line up. Built for per-contrast up/down DEG CSVs;
still carries the old Fusarium (FVEG) configuration.

**Output:** per-contrast enrichment tables and dot plots under
`<results_dir>/enrichment/<contrast>/`.

## 16. `top_GO_conserved.r`

Conservation sanity check (**H0**): are the genes on cross-species-conserved
edges (step 11) enriched for coherent biology, or is the conservation call
noise? A stripped version of `top_GO_dev2.r` that takes a **plain gene list**
(no up/down direction, no tx2gene/GTF translation — conserved-gene IDs already
match the eggNOG query IDs) and, crucially, uses **each species' own network
nodes as the background** (GO-annotated nodes from `network_<label>_node_metrics.tsv`),
not the whole genome — so enrichment reflects what is special about the
conserved genes *relative to the network*, not the generic "co-expressed genes
differ from the genome" effect. Runs both species in one pass, then compares
the two enriched-term lists (matched on GO.ID) into shared / unique sets.
Per-network `node_id_strip` handles ID-namespace suffixes (sugarcane nodes
carry a `.v2.1` tag the annotation lacks; purple needs none). Set `ontology`
to `BP` / `MF` / `CC`. **Must run inside the `topGO_env` conda env.**

**Output** under `files/network_conservation/enrichment_conserved/`:
per-species tables + dot plots in `sugarcane/` and `purple/`, plus
`GO_<ont>_conserved_{shared_terms,unique_sugarcane,unique_purple,comparison_summary}.csv`.

## 17. `GO_semanthinc_enrichment.r`

Makes the **H0 supporting figure**: does functional enrichment in the
topologically-conserved gene set of **each species independently point to the
same basic cellular processes**? Rather than merging the two species, it takes
each species' FULL significant term set (per-species tables from step 16),
builds **one common GO-semantic space from their union** (GOSemSim **Wang**
method — graph-based / IC-free, so it needs **no OrgDb**, none of which is
installed for R570 / LA purple: `godata(ont, computeIC = FALSE)`), clusters the
union once into **macro-themes** (`n_macro_clusters`, default 12; hclust ward.D2
on `1 - similarity`), embeds with PCoA, and then plots **both species in that
same space, faceted** — if both cover the same clusters, enrichment is
convergent. A second panel quantifies per-theme convergence
(sugarcane-only / both / purple-only). Runs BP/MF/CC. NB the macro-cluster
count is a *visualisation* choice, not a discovered optimum. **Must run inside
the `r_clusterprofiler` conda env** (clusterProfiler 4.14 / GOSemSim 2.32; no
rrvgo/treemap/pheatmap — clustering and plots use base hclust + ggplot2/ggrepel).

**Output** under `files/network_conservation/enrichment_conserved/semantic/`:
`GO_<ont>_conserved_semantic_space_by_species.{png,pdf}` (union of both
species' terms in one PCoA space, faceted by species, coloured by macro-theme,
sized by that species' significance, shared terms ringed),
`GO_<ont>_conserved_semantic_convergence.{png,pdf}` (per-theme
sugarcane-only/both/purple-only stacked bar),
`GO_<ont>_conserved_semantic_clusters.csv` (per union term: theme,
representative, per-species membership + p.adj, PCoA coords), and
`GO_conserved_semantic_summary.csv` (union / shared / %shared per ontology).

---

## 18. `get_tfs/` — transcription-factor identification (H1)

Identifies transcription factors / transcriptional regulators (TAPs) in each
**reference proteome** with the rule-based PlnTFDB method (Riaño-Pachón 2007;
Pérez-Rodríguez 2010), then keeps the TFs that are **nodes of the
co-expression network** — the input set for testing **H1** (is the conserved
core the MYB / TF regulatory layer both papers spotlight?). Reuses the HMM
library and rules from the original `GET_TFS/` run; the driver is a portable,
resumable, locally-parallel rewrite of the old SGE array job (no scheduler / no
`module` / no BioPython — system HMMER 3.4 + perl only).

A cohesive 7-file, config-driven pipeline (one argument = species,
`sugarcane` | `purple`). Run end-to-end:

```
cd scripts/get_tfs
./run_all.sh sugarcane      # ~1–2 h at 32 jobs x 4 cpu
./run_all.sh purple
```

Steps (also runnable individually with the same `<species>` argument):

```
config.sh                per-species paths + id-mapping regexes (edit here only)
00_check_requirements.sh  verify software + input files (executable checklist)
01_split_proteome.sh      proteome FASTA -> N_CHUNKS pieces (awk)
02_run_hmmsearch.sh       hmmsearch --cut_ga --domtblout, MAX_JOBS in parallel, concat
03_assign_families.sh     assign_family_membership.pl + RulesFull -> protein->Family->Type
04_postprocess.sh         drop Orphans, collapse to gene, intersect with network nodes
```

Parallelism is tunable via env vars (`N_CHUNKS`, `MAX_JOBS`, `CPU_PER_JOB`;
defaults 64 / 32 / 4 = 128 cores). `02_...` is resumable (finished chunks
skipped; partial writes atomic). Id reconciliation proteome→gene→node is
encoded in `config.sh` (sugarcane: strip `.<iso>.p` / `.v2.1`; purple ids are
already 1:1 with the network).

**Inputs:** `db/TF.db.hmm` (19,649 Pfam-A + PlnTFDB models, GA cutoffs),
`mytfdb/RulesFull`, `mytfdb/assign_family_membership.pl` (all under
`GET_TFS/`); the OrthoFinder reference proteomes
(`files/fix_orthofinder/{sugarcane,purple}/…`); and each network's
`…node_metrics.tsv`.

**Output** in `GET_TFS/new/results/<species>/`: `family_assignment.tsv`
(protein→Family→Type), `TF_in_network.ids` / `TF_in_network.tsv` (H1 input —
gene · Family · Type for network TFs), `family_counts.tsv` (per-family, all TF
genes vs in-network), `summary.txt`. See `GET_TFS/new/README.md` for the full
method write-up.

### `get_tfs/05_tf_characterization.r` — TF network characteristics (H1, general view)

Joins every in-network TF to its **centrality** (degree, strength,
transitivity — *not* betweenness: exact betweenness is infeasible on the
75M–681M-edge networks), its **nitrogen (`treatment`) correlation**, whether it
sits on a **cross-species conserved edge**, and its ortholog-conservation
status; then summarises per family and draws a MYB / MYB-related comparative
figure. Handles the sugarcane node-id `.v2.1` suffix and collapses the 92
multi-family (multi-isoform) genes to one row. **Run in the `r_env` conda env**
(data.table + ggplot2 + patchwork): `Rscript 05_tf_characterization.r`.

**Output** in `GET_TFS/new/results/characterization/`:
`TF_network_characteristics_{sugarcane,purple,ALL}.tsv` (big per-gene table),
`TF_family_summary.tsv` (per family × species: gene count, median degree /
strength / transitivity, % on-conserved-edge, % N-correlated),
`MYB_TF_network_overview.{png,pdf}` (3 panels: family landscape · degree vs
background · conservation & N-response). NB degree is **not** comparable across
species — the purple network is ~9× denser; compare groups *within* a species.

---

## 19. `myb61/` — locating MYB61 in both references (H1)

Kiet et al. build their nitrogen story on **ScMYB61.1**, published as
`Soff.09G0002230–3D`. That id **cannot be used**: in all 241,263 LA purple gene
ids the haplotype suffix number equals the chromosome number (0 exceptions), so
a chr9 locus with a `-3D` suffix is chimeric; and neither single-field repair
lands on a MYB (the `09G0002230-9*` copies carry `Transposase_21` /
`Transpos_assoc` / `DUF1218` and show **no Myb signal even at `E<10`**;
`03G0002230-3D` is a 230 aa protein with no Arabidopsis hit at all). The paper
provides no sequences. MYB61 is therefore re-derived from sequence.

```bash
cd scripts/myb61 && ./run_all.sh          # ./run_all.sh 03 resumes mid-way
```

Anchor **AtMYB61 = AT1G09540** → reciprocal best hit into sorghum + rice (both
return the expected **pair** of grass-WGD co-orthologues) → DIAMOND blastp into
R570 and LA purple → two independent filters: (1) reciprocal best hit back to
Arabidopsis must return AT1G09540, separating MYB61 from the ~125 other
Arabidopsis R2R3-MYBs, and (2) an independent MYB call from `get_tfs/`.
Then the OrthoFinder bridge and a MAFFT/FastTree clade test with decoys.

**Result**: 12 sugarcane + 16 purple copies in two co-orthologue clades
(chr3-type, chr9-type), joined by **5 shared orthogroups**; all strictly
confirmed copies fall in one 32-tip grass MYB61 clade. In the networks MYB61 is
**not a hub** (median degree percentile 45.7 / 45.1) and **not significantly
N-correlated in either dataset** (0 copies at |r|≥0.6 & padj≤0.05), but in
purple it is enriched on conserved edges (60% vs 24.9%, 2.4×, binomial
p=0.0041; sugarcane 37.5% vs 37.2%, p=0.62).

Two caveats the scripts surface rather than hide: OrthoFinder leaves many
polyploid copies **unassigned** (6/12 sugarcane, 3/16 purple), so absence from
`Orthogroups.tsv` is not evidence of absence (`MYB61_bridge_status.tsv`); and
`SoffiXsponR570.02Ag129700` is a real MYB61 copy whose gene model kept only one
Myb repeat, downgraded to `MYB-related` by the strict rules and recovered by the
tree. Also handled in code: the purple proteome has `.` in 542 sequence lines
(DIAMOND rejects them, HMMER does not) and `Orthogroups.tsv` uses CRLF.

`08_myb61_expression_test.r` then tests Kiet's claim on their own 18 leaf
libraries with the design it predicts — per-genotype ANOVA over 0N/2N/6N, a
**U-shape contrast** `c(+1,-2,+1)`, a linear contrast and a genotype x N
interaction — because the project's standard `gene_trait_cor.r` uses a **linear
Pearson against treatment coded 0/2/6 pooled over genotypes**, which has ~zero
power against the non-monotonic, genotype-opposed pattern they describe. Result:
a significant non-monotonic response does exist in 51NG3 (U-contrast padj<0.05
for 10/13 testable copies) but **inverted** vs their description, with 0 copies
linear; expression is very low (2/16 copies above 1 TPM), so it is suggestive
only. The paper's cloning primer encodes `MGRHSC`, matching all 28 copies found
here and none at the published locus — sequence-level proof the right gene was
identified.

**Output** in `GET_TFS/new/results/myb61/`: `<sp>_MYB61_candidates.tsv` (all
hits + both verdicts), `MYB61_orthologs_<sp>.ids`, `MYB61_orthogroups.tsv`,
`MYB61_bridge_status.tsv`, `phylogeny/MYB61.tree`, `MYB61_network_table.tsv`,
`MYB61_conserved_edge_enrichment.tsv`, `MYB61_network_overview.{png,pdf}`,
`MYB61_expression_anova.tsv`, `MYB61_expression_heatmap.{png,pdf}`.
Full method write-up in `scripts/myb61/README.md`.

---

## 20. `module20/` — Muñoz's Module 20 across both networks (H1)

Muñoz build their nitrogen claim on **Module 20**: 12 genes, ~75% MYB/MYB-related,
high betweenness in *their* de novo network. This tests whether those genes keep
MYB identity, high degree and N-responsiveness in **our** two networks.

```bash
cd scripts/module20 && ./run_all.sh        # ./run_all.sh 03 resumes mid-way
```

Mapping is **protein-level**: Muñoz's proteome
(`raw_sugarcane/transcriptome_munoz/sugar_cane.pep`) is **byte-identical to
`GET_TFS/db/sugar_cane.pep`**, so `blastp` into each reference proteome plus the
**reciprocal-Arabidopsis filter** from §19, done independently per proteome.
Members are TransDecoder ORFs; all ORFs are searched and the best-scoring one is
attributed to the member, because "longest" and "complete" disagree for at least
one member.

**Step 02b settles the classification question**: running our own TF pipeline on
**their** proteins reproduces **12/12** of their calls (9/9 MYB-family, 3/3
`no_tf`). Their "75% MYB" is correct — so every downstream disagreement is about
the orthologs we map to, not about TF calling.

**The independence problem, which dominates interpretation:** 2 of 12 members
have **no predicted ORF** (non-coding fragments); two more are isoforms of one
locus; and the 9 MYB members resolve to only **two** orthologue groups,
`AT5G59780` (**MYB59**) ×5 and `AT3G46130` (**MYB48**) ×3. That is **3–4
independent units**. Deduplicating haplotype copies moved the centrality result
from p=2.8e-05 to p=0.118, so tests are reported at gene / locus / **anchor**
level and the anchor level is the one to read.

**Result — the N response transfers, the MYB topology does not.** What maps in
volume is mostly *not* the MYB part: in-network genes come mainly from
`AT3G55960` (uncharacterised, 11/12) and `AT4G10770` (OPT7, 10/11), while MYB59
gives 9 sugarcane / **2** purple genes and MYB48 reaches neither network — hence
only 30% / 8% of mapped genes carry a MYB call. Centrality is not significant
against the MYB background (anchor-level p = 0.44 sugarcane, 0.72 purple); no
conserved-edge enrichment (degree-matched p = 0.85 / 0.77); the module is partly
recovered in sugarcane (`Module_016`, 474 genes, all 4 loci; 16 vs 24.6 random,
p=0.0003) but not in purple (18 vs 17.0, p=0.74, only the giant module).
**Nitrogen responsiveness is solid**: 24/33 genes padj<0.05 High-vs-Low N in the
responsive genotype (19/33 non-responsive), all four loci, 26/33 above 1 TPM —
Muñoz's core claim reproduces. In purple **nothing survives BH correction**. A
variance partition (step 05) says why: median variance explained is **nitrogen
20.0% vs genotype 12.1%** in sugarcane, but **nitrogen 2.1% vs genotype 40.4%**
in purple — in LA purple these genes vary genotypically, not nutritionally.

Nitrogen tests are **design-aware**: sugarcane has 2 N levels × 2 genotypes with
the three leaf segments pooled as one leaf tissue (segment kept as a sensitivity
covariate), so a per-genotype Welch test is used; purple has 3 levels and reuses
the ANOVA + U-shape machinery from `myb61/08`.

**Output** in `GET_TFS/new/results/module20/`: `module20.faa`,
`module20_orfs.tsv`, `module20_tf_call_comparison.tsv`, `module20_map_<sp>.tsv`,
`module20_mapping_summary.tsv`, `module20_network_table.tsv`,
`module20_vs_background.tsv`, `module20_tests.tsv`,
`module20_module_membership.tsv`, `module20_nitrogen_<sp>.tsv`,
`module20_network_overview.{png,pdf}` and
**`module20_heatmap_{sugarcane,purple}.{png,pdf}`** (step 05: full
ComplexHeatmap figures — absolute + z-score panels on scico `batlow`/`roma`,
columns split by genotype x nitrogen, rows by locus, clustering inside blocks).
Full write-up in `scripts/module20/README.md`.

---

## 21. `mutual_information/` — non-linear network inference (alternative to 3+4)

Replaces `pearson_cor.r` + `build_edgelist.r` with a GPU all-pairs sweep that
measures dependence rather than linear correlation, on the **same** CV≥15 gene
sets (the export step verifies this against the Pearson matrix headers and aborts
on any mismatch). Output columns match `edgelist_*_pearson.tsv`, so steps 5–14
consume it unchanged.

Estimators: `ksg` (Kraskov kNN mutual information — the real one), `gcmi`
(Gaussian-copula MI; monotone only, one matmul — `ksg − gcmi` scores how
non-linear an edge is), `xi` (Chatterjee's coefficient).

Everything is rank-transformed first, which makes the null depend only on
(n, k) rather than on the genes — so one permutation null serves all 1.46e10
pairs. Measured runtime 31 min (sugarcane, n=48) and 65 min (purple, n=18) on
the RTX A4500, resumable at tile granularity.

**The threshold is matched to `|r| = 0.8`**, i.e. `general_stats.r:PEARSON_MIN`,
which is what actually built the analysed networks — `build_edgelist.r`'s
`PEARSON_THRESHOLD <- 0.7` is only a pre-filter. Confirmed against the data:
both `network_*_filtered_edges.tsv` bottom out at |r| = 0.8. Resulting floors:
0.9891 nats (sugarcane, 2,871,301 edges) and 0.9883 nats (purple, 79,058,982).
Those two land 0.08% apart despite p-targets seven orders of magnitude apart, so
the MI layer is *more* symmetric across the two studies than the Pearson layer.

**The key caveat, measured on the first 20,000 sugarcane genes** (at the looser
0.7-matched floor, but the shape holds): ~14% of MI edges are invisible to
Pearson — that is the payoff — but ~90% of Pearson's pairs do not clear the MI
floor. Matching the false-positive rate is not matching effect size: because KSG
at n = 48 has sd ≈ 0.15 nats against a signal of ≈ 0.51 at r = 0.8, the matched
floor behaves like a Pearson cut at `|r| ≈ 0.9`, and power at the nominal
`|r| = 0.8` is ~3% (~9% at n = 18). `98_power_curve.py` reproduces those numbers.
Treat MI as a way to *find non-monotone edges*, not as a better-behaved
replacement for the Pearson networks. Full rationale in `MI_for_review.md` at
the repository root.

Also documents a provenance fact worth knowing: the **purple Pearson network came
from `china/run1` + `Group1 == "L"`**, not from `china/run2_onlyL` that the
module20/myb61 expression scripts use; their CV≥15 gene sets differ by 68 genes.

Full write-up, including the estimator comparison table and the environment
recipe, in `scripts/mutual_information/README.md`.

---

## `old/`

Deprecated/superseded script versions (`build_edgelist.r`, `cluster.r`,
`comparative_networks.r`, `general_stats_net.r`, `parallel_fullTest_xi.r`,
`pearson_cor.r`, `top_GO_dev2.r`) kept for history — not used by the current
pipeline.
