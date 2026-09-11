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
# Overridable so a parameter test can be built into a PARALLEL tree without
# touching the main analysis -- the same pattern as CLUSTERING below. Every path
# helper at the bottom of this file is expressed through RESULTS, so setting it
# redirects the whole pipeline:
#   RESULTS=$PWD/results_r09 STAT_MIN=0.9 ./run.sh build purple
# Without this the assignment was unconditional, an env override was silently
# discarded, and a rebuild overwrote 77 GB of edge tables in place.
RESULTS="${RESULTS:-${CLEAN}/results}"
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
STAT_MIN="${STAT_MIN:-0.8}"
STAT_MAX="${STAT_MAX:-0.9999}"

# The MI layer is thresholded at the MI value whose per-edge false-positive rate
# equals that of |r| = STAT_MIN at this n. Equal specificity is what licenses
# taking the union of the two layers. See docs/thresholds.md.
# Follows STAT_MIN, so raising the network threshold re-matches the MI layer to
# it automatically. That coupling is what licenses the union of the two layers,
# and it is the reason a stricter network must be REBUILT rather than filtered:
# a filtered table would carry MI p-values calibrated to the old threshold.
MATCH_PEARSON="${MATCH_PEARSON:-${STAT_MIN}}"

# Candidate cut, as a Pearson |r|. Everything at least this strong is written to
# a shard and enters the BH correction; the final network keeps only STAT_MIN
# and above. Keeping this looser than STAT_MIN is what makes BH exact over the
# full rejection set rather than only over the edges finally kept.
CAND_PEARSON="${CAND_PEARSON:-0.7}"

MAX_PERM=200000000        # ceiling on the auto-sized permutation null (KSG only)
GPU_MEM_GB=5              # working set budget; the A4500 has 20 GB
RAM_LIMIT_GB=64           # host RAM for the FDR pass (8 B per candidate edge)

# --- clustering --------------------------------------------------------------
# WHICH CLUSTERING THE MODULE-LEVEL STAGES USE: mcl | sbm
#
# `mcl` is the default and every path it produces is byte-identical to what the
# pipeline has always written, so switching to `sbm` and back cannot disturb the
# MCL results. A stochastic block model was fitted to the sugarcane network as an
# alternative (sbm/ at the repo root, not part of this pipeline); this switch is
# what lets the two be carried through the SAME downstream analysis and compared
# on their biology rather than on partition shape alone.
#
# The two clusterings differ far more than a parameter change would suggest: MCL
# gives 10,309 modules with a median of 3 genes and one holding 19% of the
# network, the SBM 1,009 blocks with a median of 50 and no giant block, and their
# adjusted Rand index is 0.0156. Do not read one's results as the other's.
#
# Overridable from the environment -- `CLUSTERING=sbm ./run.sh eigengene
# sugarcane` -- which is the whole point, since the comparison needs the same
# stages run twice. run.sh sourcing this file must NOT clobber an explicit
# choice, so the assignment is conditional.
CLUSTERING="${CLUSTERING:-mcl}"

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

# --- the MCL toolchain: native matrices, degree reduction, granularity -------
# WHY THIS EXISTS. MCL at -I 2 gives purple one module holding 28% of the network
# (47,887 genes) and a median module of 3. That is not an inflation problem.
# mclfaq(7) 7.3 names the cause exactly: "Preferably the network should not have
# nodes of very high degree... Such nodes tend to obscure cluster structure and
# contribute to coarse clusters", and clmprotocols(5) puts a number on it for
# co-expression graphs -- "the median node degree should be at most one hundred
# neighbours". Purple's median degree is 859, its p90 is 33,071 and its worst hub
# has 40,960 neighbours out of 170,736 nodes.
#
# So the handle is a k-NN degree reduction applied to the matrix MCL clusters.
# THE NETWORK ITSELF IS NOT TOUCHED -- edges, conservation and every published
# figure are unchanged; only what MCL is handed is reduced.
#
# Working directory for native binary matrices. Loading each network ONCE with
# mcxload replaces the 36 GB text .abc round-trip that 05_mcl_clustering.r pays
# on every run, which is what made a parameter sweep look unaffordable.
# Overridable too: 34_mcl_load.sh writes <study>.mci here, so two tracks would
# otherwise collide on one filename (purple's matrix is 11.3 GB).
MCL_WORK_DIR="${MCL_WORK_DIR:-${MCL_TMPDIR}/mcl_work}"

# --- choosing k-NN by graph model, not by hand -------------------------------
# k was picked heuristically for the first k-NN pass. This replaces it with a
# criterion: sweep k and keep, per network, the k whose graph is closest to a
# Barabasi-Albert model under statGraph::graph.model.selection.
#
# statGraph is installed in a PROJECT-LOCAL library rather than into any conda
# env the pipeline depends on. rARPACK is on no env here and statGraph needs it.
CLEAN_RLIB="${CLEAN_RLIB:-${CLEAN}/rlibs}"
RSCRIPT_STATGRAPH="${RSCRIPT_STATGRAPH:-/home/genomics/miniconda3/envs/lncadeep2/bin/Rscript}"

# The k grid. 50..600 in steps of 50.
KNN_BA_K_LIST="${KNN_BA_K_LIST:-50 100 150 200 250 300 350 400 450 500 550 600}"

# How many #knn reductions run at once. Each re-reads the source matrix, but it
# stays in page cache, so this is bounded by RAM rather than I/O.
KNN_BA_JOBS="${KNN_BA_JOBS:-6}"

# How many k values are SCORED at once. statGraph runs single-threaded -- its own
# numCores path deadlocks across repeated calls -- so all the parallelism is here.
# 6, not 12. statGraph opens a PSOCK cluster inside every spectral density even
# when called single-threaded, so N concurrent R processes open N clusters; at 12
# the machine ran out of ports and 3 of 16 k values died with
# "Error in serverSocket(port = port)". The jobs are ~15-25 min each, so halving
# the concurrency costs one extra batch and removes the failure mode.
KNN_BA_SELECT_JOBS="${KNN_BA_SELECT_JOBS:-6}"

# Parameter grid resolution per model inside graph.model.selection. 20 points x
# 3 models = 60 GIC evaluations per k. Measured cost of one spectral density at
# n = 170,736 is ~5 s, so a k costs single-digit minutes.
KNN_BA_NPARAM="${KNN_BA_NPARAM:-20}"

# k for the #knn() reduction, per study. EMPTY = no reduction, which is what the
# pipeline has always done. Set from 35_mcl_survey.sh using the author's own
# heuristic: the k that brings median degree toward <=100 without materially
# increasing the number of singletons.
#
# NOTE #knn INTERSECTS neighbour lists -- an edge survives only if it is among the
# top k for BOTH endpoints -- so a small k can empty the graph. #knnj joins
# instead and is the fallback; 36 records which was used.
MCL_KNN_sugarcane=""
MCL_KNN_purple=""

# Per-study inflation, chosen by the sweep. Falls back to MCL_INFLATION when
# empty, so the historical setting stays in force until the sweep has run.
# Set from the modularity optimum measured per species (figure 12): modularity is
# the only criterion with an interior optimum for MCL -- mass fraction, area
# fraction and efficiency are each monotone over the usable range and so are
# maximised at a grid boundary. Densely resampled around the peak, these are it.
# They are NOT equal, and that is the finding: purple needs a much higher
# inflation than sugarcane to reach its own optimum.
MCL_INFLATION_sugarcane=1.5
MCL_INFLATION_purple=3.5

# The work dir holding the Pearson-only matrices and the inflation ladder's
# partitions. `membership` adopts a cell from here rather than re-running mcl.
CLUSTER_WORK_DIR="${CLUSTER_WORK_DIR:-/dados04/jorge/tmp/mcl_work_cluster}"
CLUSTER_SWEEP_TREE="${CLUSTER_SWEEP_TREE:-${CLEAN}/results_cluster}"

# mcl cell tags strip the decimal point, so -I 1.5 lands in "I15". The helper
# reproduces that, and 51 refuses to adopt a cell whose .inflation sidecar
# disagrees with the value asked for.
cls_for() { echo "${CLUSTER_WORK_DIR}/sweep_$1/cls.knone.I$(printf '%s' "$2" | tr -d '.')"; }

# The sweep grid. The inflation values are the FAQ's own starting set (7.2: "A
# good set of values to start with is 1.4, 2 and 6") widened around the
# pipeline's current 2. k values are filled in per study by 35; "none" is the
# no-reduction control and must stay in the grid so the current setting is on the
# plot like any other cell.
MCL_SWEEP_I="${MCL_SWEEP_I:-1.4 2 3 4 6}"

# The unreduced control's ladder, kept separate because mcl on the full matrix is
# ~40x the work of a k-NN-reduced one (sugarcane keeps 2.4% of its arcs at
# k = 180). Measured: a full-matrix sugarcane cell is ~3.5 min, so the whole
# ladder is affordable there and this is the same list. Purple's matrix is 9.3x
# larger; shorten this to "2 4 6" for that study if the full ladder does not fit.
# Three points still answer the control's only question -- does inflation alone
# flatten the giant module?
# Per study, because purple's unreduced matrix is 9.3x sugarcane's: a full-matrix
# sugarcane cell is ~3 min, a purple one was 1 h 12 min in the original run.
MCL_SWEEP_I_NONE_sugarcane="${MCL_SWEEP_I_NONE_sugarcane:-1.4 2 3 4 6}"
MCL_SWEEP_I_NONE_purple="${MCL_SWEEP_I_NONE_purple:-2 6}"
MCL_SWEEP_KNN_RANGE="40/800/40"

# --- comparing clustering methods --------------------------------------------
# WHY. The pipeline has always clustered at -I 2, a value nobody chose. This
# sweeps inflation on the UNPRUNED Pearson-only graph and puts a second,
# graph-native method next to it (Leiden), with every partition scored against
# the SAME graph by clm info so the numbers are comparable. That last point is
# the lesson of the k-NN post-mortem: a modularity computed on each method's own
# graph compares nothing.

# The inflation ladder. mclfaq(7) 7.2 starts at "1.4, 2 and 6"; this widens it
# around the historical 2 so the current setting is one cell among many.
CLUSTER_I_LIST="${CLUSTER_I_LIST:-1.2 1.4 1.7 2 2.5 3 4 6}"

# Extra mcl resource arguments, appended verbatim (e.g. "-S 10000").
#
# THIS IS NOT A TUNING KNOB, IT IS A CONFOUND CONTROL. mcl squares the matrix
# repeatedly and prunes each column every iteration; -scheme 7 (the default AND
# the highest preset) keeps S = 1200 neighbours per node. Unpruned purple has
# mean degree 7,946, so ~85% of each node's list is discarded on the fly, among
# weights that are near-tied by construction (|r| in [0.8, 0.9999] rescaled to
# [0.01, 1]). mcl grades this itself and says "awful" (jury 20.1); sugarcane,
# mean degree 1,477, gets "deplorable" (39.2). Leiden prunes nothing, so a
# comparison at the default would partly be a comparison against mcl's pruner.
# The winning inflation is therefore re-run at higher -S to measure what the
# pruning cost, rather than assuming it cost nothing.
CLUSTER_MCL_RESOURCE="${CLUSTER_MCL_RESOURCE:-}"
CLUSTER_MCL_RESOURCE_PROBE="${CLUSTER_MCL_RESOURCE_PROBE:--S 4000|-S 10000}"

# Leiden CPM resolution ladder. CPM has no resolution limit -- unlike modularity,
# which at mean degree 7,946 would merge anything below sqrt(2m) edges -- and its
# gamma has a direct reading: a community is kept while its internal weighted
# density exceeds gamma. Because gamma is compared against EDGE WEIGHTS, the
# usable range depends on this network's weight distribution and cannot be
# guessed; 47 scouts the four decades below first and this list is set from that.
# Set from the sugarcane scout (17:16, run at 0.01/0.05/0.2/0.5): cluster count
# rose 19,571 -> 76,387 and the largest module fell 15.73% -> 4.29% across that
# range, still moving at the top end, so the ladder extends past 0.5 toward the
# maximum weight of 1.0 and adds a coarser point below 0.01.
CLUSTER_LEIDEN_GAMMA="${CLUSTER_LEIDEN_GAMMA:-0.005 0.01 0.025 0.05 0.1 0.2 0.35 0.5 0.7 0.9}"
CLUSTER_LEIDEN_SCOUT="${CLUSTER_LEIDEN_SCOUT:-0.01 0.05 0.2 0.5}"
CLUSTER_LEIDEN_ITER="${CLUSTER_LEIDEN_ITER:-2}"

# python-igraph 1.0.0 lives in the sbm env; r_net_env has no python igraph and
# leidenalg is installed nowhere. igraph's own C Leiden is used, not leidenalg:
# it reads the edge list at C level, so 676M edges never become Python objects.
CLUSTER_PYTHON="${CLUSTER_PYTHON:-/home/genomics/miniconda3/envs/sbm/bin/python}"

# Wall-clock cap for the fast-greedy (CNM) hierarchical scout, sugarcane only.
# Classical hierarchical clustering is not attempted at all: it needs a dense
# 170,135^2 dissimilarity matrix (232 GB) built from correlations this pipeline
# deliberately thresholded away. fast-greedy is the graph-native substitute that
# still yields a dendrogram; if it does not finish inside the cap, that IS the
# result and purple is not attempted.
CLUSTER_FASTGREEDY_CAP_S="${CLUSTER_FASTGREEDY_CAP_S:-14400}"

# --- clustering quality ------------------------------------------------------
# Sorensen-Dice annotation homogeneity, the metric cogeqc::calculate_H uses for
# orthogroups, applied to modules. cogeqc itself is NOT used: it skips groups
# larger than max_size (200) entirely -- which would silently refuse to score the
# very giant modules being diagnosed -- and it enumerates every pair through
# combn() in an R loop. 37 reimplements the same score sparse and vectorised.
#
# Modules above this size are SUBSAMPLED rather than skipped, with a fixed seed,
# and the subsampling is recorded per module so a score is never mistaken for
# exact.
HOMOGENEITY_MAX_GENES=2000

# Homogeneity rises trivially as modules shrink, and the whole point of raising
# inflation is to make modules smaller. Only the excess over a SIZE-MATCHED
# random partition is evidence, so every clustering is scored against one.
HOMOGENEITY_PERM=100
HOMOGENEITY_SEED=1188

# --- stochastic block model (CLUSTERING=sbm) ---------------------------------
# Output of the graph-tool nested fit, which lives outside this pipeline. The
# per-study directory names are the fit's, not the pipeline's.
SBM_ROOT="${BASE}/sbm/output"
SBM_DIR_sugarcane="${SBM_ROOT}/sugar"
SBM_DIR_purple="${SBM_ROOT}/purple"
SBM_TAG=01_minimize

# Which level of the nested hierarchy is "the modules". 0 is the finest and the
# only one at module granularity -- 1,009 blocks, median 50 genes, max 1,366.
# Level 1 already coarsens to 294 blocks with a median of 276, and by level 3 the
# partition is 49 blocks and no longer a module set.
SBM_LEVEL=0

# Newman modularity of the SBM partition, computed on the same graph so it is
# comparable with the number 05_mcl_clustering.r reports for MCL. It is the only
# slow step (it reads the 76 M-row edge table); set 0 to skip and write NA. An
# SBM minimises description length rather than modularity, so a low Q is the two
# methods optimising different things, not a failure.
SBM_COMPUTE_Q=1

# --- topology ----------------------------------------------------------------
# Local transitivity is O(sum deg^2) and took ~15 h on purple. Nothing in the
# pipeline TESTS it -- it appears only as a reported column in the TF and MYB61
# readouts. FALSE writes NA and skips the transitivity-vs-degree plot.
COMPUTE_TRANSITIVITY=0

# --- conservation ------------------------------------------------------------
ORTHOGROUPS="${BASE}/files/fix_orthofinder/proteins/OrthoFinder/Results_Jun04_2/Orthogroups/Orthogroups.tsv"

# DO NOT REPOINT ORTHOGROUPS. Every conservation number this pipeline reports
# was computed from the two-species run above. The dN/dS stage (scripts/29_dnds)
# adds a Sorghum bicolor outgroup through a SEPARATE three-species OrthoFinder
# run in files/orthofinder_3sp/. The file above stays its analysis unit: sorghum
# is attached to THESE single-copy pairs, so the dN/dS and conservation layers
# describe the same genes. Its own paths live in scripts/29_dnds/config.sh, which
# is where that stage's parameters belong -- this file stays the source of truth
# for the stages run.sh drives directly.
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
MODULE_SEED=1188

# --- the design in the module-trait test -------------------------------------
# The module response is now a BLOCKED partial correlation, not a marginal one.
# The project's own QC puts genotype at R2 = 0.999 of PC1 in purple and 0.998 in
# sugarcane, and leaf segment at 0.802 of PC2, so the largest variance component
# in either matrix was sitting in the residual of every nitrogen test. The same
# correction at gene level moved purple from 30 responsive genes to 2,331.
#
#   sugarcane   eigengene ~ genotype + segment + N     n = 48, resid df 42
#   purple      eigengene ~ genotype + N               n = 18, resid df 15
MODULE_BLOCK_sugarcane="genotype segment"
MODULE_BLOCK_purple="genotype"

# Spearman stays primary at module level -- purple's trait is an ordinal 0/2/6 mM
# dose and sugarcane's is two-level, where Spearman on eigengene ranks is the
# rank-biserial correlation. Blocked Pearson is computed and written beside it.
MODULE_BLOCKED_STAT_sugarcane=spearman
MODULE_BLOCKED_STAT_purple=spearman

# Sugarcane's 48 libraries are 12 plants x 4 leaf segments -- repeated measures,
# not 48 replicates. Segment as a fixed block removes the segment MEANS but not
# the within-plant correlation, so a plant-level run (n = 12, resid df 9) is
# computed alongside and correlated against the blocked fit. Empty = no control.
MODULE_PLANT_FROM_sugarcane="genotype"
MODULE_PLANT_FROM_purple=""

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
FIG_MODULE_GO=6
FIG_MODULE20=7

# The MCL-vs-SBM comparison is not a paper figure yet -- it is the evidence for a
# decision that has not been made. Numbered past the seven so it cannot be
# mistaken for one.
FIG_CLUSTERING=8

# Module-20 figure. The focus copy is the ONLY AtMYB59-anchor copy expressed in
# purple (60.6 TPM against 0.01 and 0.00 for the other two) and the only
# Module-20 gene there with a significant U-shape response to nitrogen.
M20_FOCUS_GENE=Soffic.09G0001580-9H
M20_FOCUS_ANCHOR=AT5G59780

# Module-GO figure: how many terms per direction in panel B, how many recurrent
# terms in panel D, and which study carries those two panels. sugarcane, because
# purple has 10 individually testable modules against sugarcane's 49 and the
# comparison would be between a result and an absence.
MODULE_GO_FIG_NTERMS=8
MODULE_GO_FIG_NRECUR=10
MODULE_GO_FIG_MAIN=sugarcane

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

# Where the module-level stages read their clustering from, and where they write.
# With CLUSTERING=mcl both are exactly what they have always been, so the
# existing results cannot be touched; with sbm the outputs land in a parallel
# subdirectory and the two sets sit side by side.
clus_prefix() { echo "$(study_dir "$1")/${CLUSTERING}_$1"; }
module_dir()  {
  if [ "$CLUSTERING" = "mcl" ]; then echo "$(study_dir "$1")"
  else echo "$(study_dir "$1")/${CLUSTERING}"; fi
}
vst_prefix()  { echo "${RESULTS}/$1/vst/$1"; }
layer_out()   { echo "${RESULTS}/$1/layers/$1_$2"; }
# The Pearson/MI layers are SOURCE data, shared by every track: they are what the
# GPU engine produced and no parallel tree recomputes them. A tree redirected with
# RESULTS= must still read the layer that actually exists, so this helper
# deliberately ignores RESULTS. Use it only for inputs, never for outputs.
main_layer_out() { echo "${CLEAN}/results/$1/layers/$1_$2"; }
network_tsv() { echo "${RESULTS}/$1/network_$1_edges.tsv"; }
