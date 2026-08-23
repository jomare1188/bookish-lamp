# =============================================================================
# config.sh — the single source of truth for the new_clean pipeline
#
# Sourced by run.sh, and exported into the environment for the R and Python
# stages. Nothing downstream hardcodes a path: if a value is not here, it is a
# bug, not a convention.
#
# Per-study values use the `VAR_<study>` naming convention and are read with
# bash indirect expansion, e.g.  eval "dds=\$DDS_$study".
# =============================================================================

BASE="/dados04/jorge/comparative_saccharum"
CLEAN="${BASE}/new_clean"
SCRIPTS="${CLEAN}/scripts"
RESULTS="${CLEAN}/results"
LOGS="${CLEAN}/logs"

STUDIES="sugarcane purple"

# --- interpreters ------------------------------------------------------------
# cor_env has DESeq2 + matrixStats but NOT igraph/hexbin, so it cannot run the
# network scripts. r_net_env has the full set (igraph, data.table, ggplot2,
# hexbin, scales). torch+cu130 exists only in the `docling` env on this machine.
RSCRIPT_DESEQ="/home/genomics/miniconda3/envs/cor_env/bin/Rscript"
RSCRIPT_NET="/home/genomics/miniconda3/envs/r_net_env/bin/Rscript"
RSCRIPT_PLOT="/home/genomics/miniconda3/envs/r_env/bin/Rscript"
PYTORCH="/home/genomics/miniconda3/envs/docling/bin/python"

# GO stages need their own envs and are launched via `conda run`.
CONDA_TOPGO="topGO_env"
CONDA_CLUSTERPROFILER="r_clusterprofiler"

# --- inputs: ONE quantification per study, used by every stage ---------------
#
# sugarcane: Munoz-Perez 2025, R570 reference, all 48 leaf libraries.
DDS_sugarcane="${BASE}/run1/salmon/deseq2_qc/deseq2.dds.RData"
COLS_sugarcane=""
TPM_sugarcane="${BASE}/run1/salmon/salmon.merged.gene_tpm.tsv"
SALMONQC_sugarcane="${BASE}/run1/multiqc/multiqc_report_data/multiqc_salmon.txt"
META_sugarcane="${BASE}/samplesheet.csv"
STRIP_VERSION_sugarcane=1          # R570 ids carry a .v2.1 suffix

# purple: Ta Quang Kiet 2025, LA purple reference, 18 leaf libraries.
#
# run2_onlyL, NOT run1. The old pipeline was internally inconsistent: the
# network came from china/run1 narrowed to the 18 leaf libraries by
# Group1=="L", while gene_trait_cor.r, myb61/08 and module20/04 read
# china/run2_onlyL. Those are different quantifications, not different views of
# one -- measured on 2026-08-13, the raw counts differ by up to 6,398, the VST
# by 0.071 on average (max 2.34), and the CV>=15 gene sets differ by 68 genes.
#
# run2_onlyL is the dedicated leaf-only run and is ALREADY exactly those 18
# libraries, so no column filter is needed and n = 18 is unchanged. Everything
# in docs/thresholds.md that depends on n therefore carries over untouched.
DDS_purple="${BASE}/china/run2_onlyL/salmon/deseq2_qc/deseq2.dds.RData"
COLS_purple=""
TPM_purple="${BASE}/china/run2_onlyL/salmon/salmon.merged.gene_tpm.tsv"
SALMONQC_purple="${BASE}/china/run2_onlyL/multiqc/multiqc_report_data/multiqc_salmon.txt"
META_purple="${BASE}/china/samplesheet_china.csv"
STRIP_VERSION_purple=0             # LA purple ids are already bare

# --- network construction ----------------------------------------------------
MIN_CV=15                 # CV filter on RAW counts; matches the original pipeline
ESTIMATORS="pearson ksg"  # the linear layer and the non-linear layer
KSG_K=3                   # Kraskov k; 3 is the usual small-sample choice
ALPHA=0.05                # BH FDR level

# Edge-strength window, on the correlation scale. STAT_MIN is the real network
# threshold: the old build_edgelist.r used 0.7, but that was only a pre-filter
# for an intermediate file -- the networks actually analysed came from
# general_stats.r, whose PEARSON_MIN was 0.8. STAT_MAX drops suspiciously
# perfect correlations.
STAT_MIN=0.8
STAT_MAX=0.9999

# The MI layer is thresholded at the MI value whose per-edge false-positive rate
# equals that of |r| = STAT_MIN at this n. Equal specificity is what licenses
# taking the union of the two layers. See docs/thresholds.md.
MATCH_PEARSON="${STAT_MIN}"

# Candidate cut, as a Pearson |r|. Everything at least this strong is written to
# a shard and enters the BH correction; the final network keeps only STAT_MIN
# and above. Keeping this looser than STAT_MIN is what makes BH exact over the
# full rejection set rather than only over the edges finally kept.
CAND_PEARSON=0.7

MAX_PERM=200000000        # ceiling on the auto-sized permutation null (KSG only)
GPU_MEM_GB=5              # working set budget; the A4500 has 20 GB
RAM_LIMIT_GB=64           # host RAM for the FDR pass (8 B per candidate edge)

# --- clustering --------------------------------------------------------------
MCL_INFLATION=2

# 2, matching the original mcl_clustering.r. This does NOT affect the clustering
# -- modularity is computed on the raw partition -- it only decides how small a
# module may be and still get a Module_NNN name instead of "Unassigned". At 5,
# sugarcane went from 8,691 named modules / 105 unassigned genes to 3,307 /
# 18,935, which is a big change in how results read for no change in the
# underlying partition. Keep it at the pipeline's historical value; raise it
# deliberately if tiny modules turn out to be noise.
MCL_MIN_MODULE_SIZE=2
MCL_MIN_MODULE_SIZE_PLOT=10
NUM_CORES=100

# The old pipeline always wrote an 8 GB .rds holding the full igraph objects for
# both networks. Nothing reads it -- every consumer reads mcl_*_membership.tsv.
MCL_SAVE_GRAPH_RDS=0

# The mcl binary. It ships inside the r_net_env conda env -- the same one the
# network R stages use -- but that env is never activated, so its bin/ is not on
# PATH and mcl_clustering.r's system2("mcl") would fail. run.sh prepends this
# directory rather than making the caller activate anything.
MCL_BIN="/home/genomics/miniconda3/envs/r_net_env/bin/mcl"

# mcl writes a temporary .abc file that is tens of GB for purple. Do not let it
# land on a small system /tmp.
MCL_TMPDIR="/dados04/jorge/tmp"

# --- topology ----------------------------------------------------------------
# Local transitivity is O(sum deg^2) and took ~15 h on purple. Nothing in the
# pipeline TESTS it -- it appears only as a reported column in the TF and MYB61
# readouts. FALSE writes NA and skips the transitivity-vs-degree plot.
COMPUTE_TRANSITIVITY=0

# --- conservation ------------------------------------------------------------
ORTHOGROUPS="${BASE}/files/fix_orthofinder/proteins/OrthoFinder/Results_Jun04_2/Orthogroups/Orthogroups.tsv"
CHUNK_SIZE=2000000

# Permutation null for the conservation rate (13_conservation_null.r).
# The target adjacency is built once and reused, and the source network is
# Bernoulli-sampled, so replicates are cheap. 5e6 edges gives a standard
# error of ~0.013% on a 10% rate -- far finer than the effect being tested.
NULL_REPS=20
NULL_SAMPLE=5000000

# --- gene-trait correlation --------------------------------------------------
# Trait encoding, per study. Format:
#     trait:LEVEL=value,LEVEL=value;trait:LEVEL=value,...
# These replace the mutually exclusive TRAIT_ENCODING comment blocks that had to
# be swapped by hand in gene_trait_cor.r. Levels are matched against the named
# column of the study's samplesheet; a sample whose level is not listed is
# dropped from that trait's correlation.
#
# `treatment` is the nitrogen axis and is what the gene selection uses. Sugarcane
# is a two-level contrast. Purple is encoded 0/2/6 mM, but READ THE NEXT
# PARAGRAPH BEFORE INTERPRETING ANYTHING THAT USES IT.
#
# Purple's three levels are NOT a dose-response gradient. 2 N is the CONTROL, and
# 0 N and 6 N are stresses in opposite directions -- deficiency and excess. So a
# monotonic test on this coding asks "does expression track nitrogen SUPPLY?",
# which is a real question, but it is NOT "does this gene respond to nitrogen
# stress?" -- a gene moved the same way by both stresses is invisible to it, and
# that is the shape the design predicts.
#
# Measured, so the cost is known rather than feared: at the module level, a
# U-shape contrast c(+1,-2,+1) over the three levels finds ONE significant purple
# module at padj <= 0.05, and Spearman misses that one. Spearman finds 79. So the
# monotonic test is not leaving a large non-monotonic set on the table -- at
# n = 18 the U-shape test has almost no power, which is the same wall everything
# else in purple hits. The 79 should still be described as tracking nitrogen
# supply, not as stress responders.
TRAITS_sugarcane="genotype:RB975375=1,RB937570=0;treatment:High Nitrogen=1,Low Nitrogen=0"
TRAITS_purple="genotype:51NG3=1,TAGZ=0;treatment:0N=0,2N=2,6N=6"

# The trait the gene selection is made on, and its thresholds.
SELECT_TRAIT=treatment

# Which statistic makes a gene "nitrogen-responsive" in 08_conserved_cor_genes.r:
#   pearson | mi | union
# `union` is the default because finding a conserved response that Pearson cannot
# see is the reason the MI layer exists -- a Pearson-only search returning nothing
# does not distinguish "no shared response" from "no LINEAR shared response".
# Overridable per run:  ./run.sh conscor pearson
TRAIT_SELECTION=union

# Discovery species for the DIRECTED test (./run.sh conscor 1). The other
# species' p-values are then corrected over only the orthologs of the
# discovery species' responsive genes, which is the burden the comparative
# question actually implies -- a genome-wide BH over 44,118 genes is testing
# a hypothesis nobody asked. Discovery should be the better-powered study.
TRAIT_DISCOVERY=sugarcane
TRAIT_R_THR=0.6
TRAIT_PADJ_THR=0.05

# Gene-vs-trait mutual information (12_gene_trait_mi.py). The trait is discrete,
# so this uses the Ross (2014) estimator, not KSG. The null is over label
# permutations and is shared by every gene, so it is cheap to make deep: BH over
# ~40,000 genes needs p ~ 1.3e-6 for the strongest gene, and 1e7 permutations
# resolve that empirically rather than by extrapolation.
TRAIT_MI_K=3
TRAIT_MI_PERM=10000000

# R_THR is the H1 readouts' own effect-size cut when they call a gene
# "nitrogen-responsive". PADJ_THR is NOT readout-only despite where it sits: the
# module stages take it as their significance threshold too (`moduletrait
# --alpha`, and `moduleprofile`'s final call). Changing it moves the module
# results as well as the readouts.
R_THR=0.7
PADJ_THR=0.05

# --- module-level analysis ---------------------------------------------------
# Minimum genes for a module to get an eigengene. 3 keeps 6,576 sugarcane /
# 6,318 purple modules, covering 92% and 87% of each network's genes. It only
# modestly eases the testing burden -- purple's strongest module still needs
# |r| ~ 0.85 to clear BH at n = 18, against 0.88 testing every gene -- but the
# eigengenes of larger modules do not change if this is raised later.
MIN_MODULE_SIZE_EIGEN=3

# The module response is called by SPEARMAN ONLY (19_module_trait_spearman.r).
#
# NOT PEARSON: the trait is ordinal, not interval. Purple's nitrogen is a dose
# (0/2/6 mM) and Pearson reads that spacing literally -- it asks whether a module
# moves exactly twice as far from 2 to 6 mM as it does from 0 to 2. Nothing in
# the design justifies that. Spearman asks what the experiment poses: does the
# module move monotonically with nitrogen? Sugarcane's two-level trait makes it
# the rank-biserial correlation, robust to the outliers a PC1 can carry.
#
# NOT MI: it is an omnibus test, firing on any dependence at all including a
# dispersion change in one or two libraries -- which is what a third of the old
# `mi_only` modules turned out to be. It also needed an effect-size floor that
# could only be CALIBRATED against the linear one, never derived, because the
# network's nats-to-|r| identity assumes two continuous variables and this trait
# is discrete. One statistic, one responsive set, no response classes.
#
# The two thresholds mirror the gene-level TRAIT_R_THR / TRAIT_PADJ_THR so the
# module counts stay comparable with the gene ones. padj is BH over every module
# tested in the study; the raw `pval` is kept in the output so an uncorrected
# reading needs no re-run.
MODULE_R_THR=0.6
MODULE_PADJ_THR=0.05

# Label-permutation null for the module response: the same eigengenes against
# shuffled trait labels, counting how many clear both thresholds. This is a null
# for the TRAIT association only -- the eigengenes stay exactly as correlated
# with each other as they really are. 1,000 permutations cost seconds.
MODULE_PERM=1000

# Heatmaps: how many responsive modules to draw, and how many genes of each.
# The largest modules are 19,604 (sugarcane) and 47,887 (purple) genes, which no
# heatmap can render; those are subset to the top genes by intramodular strength.
HEATMAP_TOP_N=20
HEATMAP_MAX_GENES=100

# Summary figure: cap on rows. Above this, the top N per response DIRECTION
# by significance are shown, so the smaller direction is not crowded out.
SUMMARY_MAX_MODULES=250

# Extra sample-metadata columns to group heatmap columns by, after the two
# traits. Sugarcane's four leaf segments otherwise interleave -- its sample names
# run B0_1, B_1, M_1, P_1, B0_2, ... so replicates of one segment sit four columns
# apart and the segment effect reads as striping across the figure.
#
# `segment`, NOT `tissue`. The sheet's own `tissue` column is not the design: it
# collapses base0 and base into one "Leaf Base" of 24 libraries, calls mid "Leaf"
# and tip "Leaf Apex". Grouping on it therefore left base0 and base interleaved
# inside a block of 6 -- the striping this setting exists to remove. `segment`
# carries the four real levels (base0, base, mid, tip; 12 libraries each), and
# the PCA in the dataset figure measures the difference: leaf segment explains
# 0.802 of sugarcane's PC2 on the four true levels against 0.450 on the
# collapsed three. Purple has no segment column and skips this.
HEATMAP_GROUP_BY="segment"

# Per-module GO enrichment (18_module_go.r). One topGO run per responsive module,
# with the response classes pooled -- the question is what responsive modules do,
# not what separates the classes. GO_P below is reused unchanged, same threshold
# and same reasoning.
#
# This is only the DEFAULT ontology; all three are run and kept
# (`./run.sh modulego <study> MF`). The universe is ontology-agnostic, so the same
# modules are testable in each and the three are comparable module by module --
# which is what makes MF's enzyme names a check on BP's process calls.
MODULE_GO_ONTOLOGY=BP

# Minimum GO-annotated members for a module to be tested at all. Deliberately
# low: a 3-gene module CAN reach p < 0.05 against a ~25,000-gene background, so
# this only skips modules where the test is undefined, it does not pre-judge
# small ones. n_annotated is written per module, so filtering harder afterwards
# needs no re-run. The median responsive module holds 5 genes -- most of the
# module set is small, and that fact belongs in the output rather than in a cut.
MODULE_GO_MIN_ANNOTATED=3

# Modules are independent topGO runs over one shared graph, so this forks cleanly
# (copy-on-write: the graph is not duplicated per worker).
MODULE_GO_CORES=16

# --- paper figures -----------------------------------------------------------
# Figures carry a panel letter and the labels the data needs to be read, and
# nothing else -- no titles, no subtitles, no statistics printed on the panel.
# Everything else goes in the legend. Each figure script generates its own
# legend from the variables that drew it, and `./run.sh legends` concatenates
# them here, so a number cannot disagree between a figure and its legend.
#
# This file is OUTSIDE results/ deliberately: results/ is gitignored as
# regenerable output, and the legends are manuscript text.
FIGURE_LEGENDS="${BASE}/figures_legends.txt"

# Figure NUMBERS live here and nowhere else. Scripts are named for what they draw
# (20_fig_reproduction.r, 21_fig_dataset_qc.r) because figure order is editorial
# and has already changed once -- the reproduction figure opened the paper until
# the dataset/QC figure took the front. These numbers name the output files and
# open each generated legend, so renumbering the paper is this one edit.
FIG_DATASET=1
FIG_REPRODUCTION=2
FIG_TOPOLOGY=3
FIG_CONSERVATION=4
FIG_MODULES=5

# Conservation figure, panel B: which ontology the conserved-set GO panel draws,
# and how many shared terms it shows. Terms are ranked by the WORSE of the two
# species' p-values, so the panel shows agreement rather than one species' hits.
# 10 is what fits legibly once the full GO names are wrapped -- topGO truncates
# them at 40 characters and the figure expands them from a GO.db cache.
CONS_GO_ONTOLOGY=BP
# 12 rather than 10 so that `glutamate biosynthetic process` and the `ammonia
# assimilation cycle` -- the two shared terms that speak directly to the trait --
# stay on the panel; they rank 11th and 12th by agreement.
CONS_GO_NTERMS=12

# Genes used for the per-study PCA in the dataset figure: the most variable
# 2,000, the usual DESeq2 plotPCA convention.
PCA_NTOP=2000

# Points per curve in the topology figure's log-log CCDFs. Evaluating at every
# unique degree would put ~41,000 points in the SVG for purple and draw the same
# curve; 300 log-spaced values are indistinguishable and keep the vector small.
TOPO_GRID=300

# --- GO ----------------------------------------------------------------------
ONTOLOGIES="BP MF CC"
# Threshold on the RAW weight01 p-value, topGO's own convention: weight01
# conditions each term on its DAG neighbours, so its p-values are not an
# exchangeable family and BH does not apply to them. A BH column is still
# written to the output for reference; it does not select the terms.
GO_P=0.05
GO_NTOP=20
EMAPPER_sugarcane="${BASE}/annotation/sugarcane/emapper.annotations"
EMAPPER_purple="${BASE}/annotation/purple/emapper.annotations"

# --- derived paths (do not edit) ---------------------------------------------
study_dir()   { echo "${RESULTS}/$1"; }
vst_prefix()  { echo "${RESULTS}/$1/vst/$1"; }
layer_out()   { echo "${RESULTS}/$1/layers/$1_$2"; }
network_tsv() { echo "${RESULTS}/$1/network_$1_edges.tsv"; }
