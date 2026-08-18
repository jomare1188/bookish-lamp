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
# `treatment` is the nitrogen axis and is what the gene selection uses. Purple
# encodes it as a DOSE (0/2/6 mM) because that study is a gradient; sugarcane is
# a two-level contrast.
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

# Used by the H1 readouts when they call a gene "nitrogen-responsive".
R_THR=0.7
PADJ_THR=0.05

# --- module-level analysis ---------------------------------------------------
# Minimum genes for a module to get an eigengene. 3 keeps 6,576 sugarcane /
# 6,318 purple modules, covering 92% and 87% of each network's genes. It only
# modestly eases the testing burden -- purple's strongest module still needs
# |r| ~ 0.85 to clear BH at n = 18, against 0.88 testing every gene -- but the
# eigengenes of larger modules do not change if this is raised later.
MIN_MODULE_SIZE_EIGEN=3

# Effect-size floor on a module's LINEAR response, mirroring TRAIT_R_THR at the
# gene level. Without it padj alone admits modules down to |r| = 0.36 at n = 48,
# and the module counts are not comparable to the gene-level ones. No floor is
# applied to MI: the nats-to-|r| identity the network uses assumes two continuous
# variables and the trait here is discrete, so `mi_norm` (MI / H(trait)) is
# reported instead as a 0-1 effect size.
MODULE_R_THR=0.6

# Heatmaps: how many responsive modules to draw, and how many genes of each.
# The largest modules are 19,604 (sugarcane) and 47,887 (purple) genes, which no
# heatmap can render; those are subset to the top genes by intramodular strength.
HEATMAP_TOP_N=20
HEATMAP_MAX_GENES=100

# Summary figure: cap on rows. Above this, the top N per response class by
# significance are shown, so the smaller non-linear class is not crowded out.
SUMMARY_MAX_MODULES=250

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
