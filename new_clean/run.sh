#!/usr/bin/env bash
# =============================================================================
# run.sh — the one entry point for the new_clean pipeline
#
# Every stage is `./run.sh <stage> [argument]`. No script is ever edited between
# runs: what used to be a comment-block swap or a hand-set DIRECTION/ontology
# variable is now the argument. That is the whole reason this file exists --
# the old tree needed at least 8 manual edits to complete one pass, and an edit
# left in the wrong state produced results that looked fine and were wrong.
#
# Everything is logged to logs/<stage>_<arg>_<timestamp>.log as well as stdout.
#
#   ./run.sh export    sugarcane|purple      dds -> VST binary
#   ./run.sh validate                        engine correctness suite
#   ./run.sh network   <study> pearson|ksg   one layer of one network
#   ./run.sh merge     <study>               layers -> the network
#   ./run.sh stats     <study>               node + global metrics, plots
#   ./run.sh mcl       <study>               MCL modules
#   ./run.sh nodemetrics <study>             gene universe from the .mci
#   ./run.sh membership  <study>             adopt the ladder partition at the chosen -I
#   ./run.sh pearsononly <study>             Pearson layer -> a Pearson-only network
#   ./run.sh sgvalidate                      can statGraph's criterion choose k here?
#   ./run.sh knnsweep  <study>               reduce at each k; nodes/edges/degrees
#   ./run.sh knnselect <study> [k]           score each k against ER/WS/BA
#   ./run.sh knncollect <study>              merge the scores, name the winning k
#   ./run.sh scoreref  <study>               score every k's partition on ONE graph
#   ./run.sh pearsonmci <study>              Pearson layer -> MCL matrix, direct
#   ./run.sh mclload   <study>               network -> native MCL matrix (once)
#   ./run.sh mclsurvey <study>               degree survey; choose the k-NN cut
#   ./run.sh mclsweep  <study> [k-list]      inflation x k-NN grid + clm info/dist
#   ./run.sh mclladder <study> [mcl-args]    inflation ladder, no k-NN reduction
#   ./run.sh leidensweep <study> [scout]     Leiden CPM/modularity + Louvain
#   ./run.sh clustercompare <study>          score every method on ONE graph
#   ./run.sh fastgreedy <study>              CNM hierarchical scout, capped
#   ./run.sh figclustermethods               figure 12: MCL vs Leiden vs Louvain
#   ./run.sh clusterhomog <study>            annotation homogeneity of each cell
#   ./run.sh figclusterchoice                figure 10, the granularity choice
#   ./run.sh sbmclust  <study>               SBM fit -> the same two files
#   ./run.sh conserve  sugarcane_to_purple|purple_to_sugarcane
#   ./run.sh conservenull <direction>        permutation null for the above
#   ./run.sh trait     <study>               gene-trait correlations
#   ./run.sh traitmi   <study>               gene-vs-trait MI (non-linear)
#   ./run.sh eigengene <study>               one eigengene per MCL module
#   ./run.sh moduletrait <study>             module response: Spearman rho vs trait
#   ./run.sh moduleprofile <study>           + TF enrichment per module
#   ./run.sh moduleheatmap <study> [mods]    heatmaps for responsive modules
#   ./run.sh modulesummary <study>           one figure: all responsive modules
#   ./run.sh gene2go   <study>              derive GO from InterPro + Pfam
#   ./run.sh modulego  <study> [BP|MF|CC]    GO enrichment per responsive module
#   ./run.sh conscor   [0|1] [selection]     conserved N response; 1 = directed test
#   ./run.sh ushape    [study]               U-shape contrast, genome-wide (purple)
#   ./run.sh traitblocked <study>            gene-trait with the design in the model
#   ./run.sh consblocked  [0|1]              conserved N response at NODE level, blocked
#   ./run.sh go        BP|MF|CC              GO enrichment
#   ./run.sh gosem                           GO semantic clustering
#   ./run.sh tfs       <study>               TFs in the network (step 04 only)
#   ./run.sh tfchar                          TF characterisation
#   ./run.sh myb61     [from-step]           MYB61 readout (default 07)
#   ./run.sh module20  [from-step]           Module-20 readout (default 03)
#   ./run.sh dnds      [from-step]           dN/dS vs network conservation (default 00)
#   ./run.sh figdataset                      dataset, QC and quantification figure
#   ./run.sh figtopology                     network topology of both networks
#   ./run.sh figconservation                 cross-species edge conservation
#   ./run.sh figmodules                      module-level nitrogen response
#   ./run.sh figmodulego                     module-level GO, incl. by direction
#   ./run.sh figmodule20                     Munoz Module 20 in purple
#   ./run.sh figclustering [study]           MCL vs SBM, same downstream analysis
#   ./run.sh figrepro                        both source studies reproduced
#   ./run.sh legends                         assemble figures_legends.txt
#   ./run.sh build     <study>               export + both layers + merge
#
# Order: export -> validate -> network(x2) -> merge -> stats -> mcl ->
#        conserve(x2) -> trait(x2) -> conscor -> go(x3) -> gosem -> readouts
#
# WHY THE BODY IS INSIDE main(): bash reads a script INCREMENTALLY, from a saved
# byte offset. Edit a script while an invocation of it is still running and that
# shell resumes reading at the old offset in the new file -- i.e. in the middle
# of a different case branch. That really happened here on 2026-08-13: an edit
# during a 2.5 h sweep made the finished run fall into the `gosem` branch and
# exit 1 with the sweep's results perfectly intact. A function body is parsed as
# a single unit before any of it executes, so the whole file is in memory by the
# time anything runs and a mid-run edit cannot reach it.
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/config.sh"

usage() { sed -n '3,39p' "${BASH_SOURCE[0]}" | sed 's/^# \?//'; exit 1; }
die()   { echo "ERROR: $*" >&2; exit 1; }

# Per-study value lookup, e.g. cfg DDS sugarcane -> $DDS_sugarcane
cfg() { local v="$1_$2"; echo "${!v-}"; }

# --- refuse to overwrite a network built at a DIFFERENT threshold -----------
# STAT_MIN and RESULTS are both env-overridable so a parameter test can be built
# into a parallel tree. The failure mode that creates is silent and expensive:
# one forgotten RESULTS= and a 0.9 rebuild lands on top of the 0.8 network --
# 77 GB of edge tables and ~5.5 h of GPU, with nine call sites below silently
# picking up the new file. The threshold is recorded in the summary json that
# sits beside every output, so this checks it rather than trusting the caller.
#
# FORCE=1 overrides, for the case where overwriting really is the intent.
assert_threshold() {   # $1 = summary json, $2 = expected value, $3 = what
  local js="$1" want="$2" what="$3" have
  [ -f "$js" ] || return 0                      # nothing to overwrite
  have=$(sed -n 's/.*"stat_min": *\([0-9.]*\).*/\1/p;s/.*--match-pearson \([0-9.]*\) .*/\1/p' "$js" | head -1)
  [ -n "$have" ] && [ "$have" != "$want" ] || return 0
  cat >&2 <<EOF

REFUSING TO OVERWRITE: $(basename "$js") records ${what} = ${have}, this run has ${want}.
  target : $(dirname "$js")
  You are about to replace a network built at a different threshold with one
  built at this threshold, in place. If that is a parameter test, send it to a
  separate tree instead:
      RESULTS=\$PWD/results_r${want//./} STAT_MIN=${want} ./run.sh ...
  If overwriting really is the intent, re-run with FORCE=1.

EOF
  exit 1
}

check_study() {
  case " $STUDIES " in
    *" $1 "*) ;;
    *) die "unknown study '$1' (expected one of: $STUDIES)" ;;
  esac
}

main() {
  STAGE="${1-}"; [ -n "$STAGE" ] || usage
  ARG="${2-}"
  shift $(( $# > 1 ? 2 : $# ))
  EXTRA=("$@")

  mkdir -p "$LOGS"
  LOG="${LOGS}/${STAGE}${ARG:+_$ARG}_$(date +%Y%m%d_%H%M%S).log"
  exec > >(tee -a "$LOG") 2>&1
  echo "=== ${STAGE} ${ARG} :: $(date '+%F %T') :: log ${LOG}"
  # Announced on every run. A module-level stage silently using the wrong
  # clustering writes plausible-looking output to the wrong place, and the only
  # way to catch it is to be told which one is active before the work starts.
  case "$STAGE" in
    network|merge|build|stats|mcl|mclload|mclsurvey|mclsweep|clusterhomog)
      echo "=== results tree: ${RESULTS}"
      echo "=== threshold   : STAT_MIN=${STAT_MIN}  MATCH_PEARSON=${MATCH_PEARSON}  CAND_PEARSON=${CAND_PEARSON}" ;;
  esac
  case "$STAGE" in
    eigengene|moduletrait|moduleprofile|moduleheatmap|modulesummary|modulego|figtopology|figmodules|figmodulego)
      echo "=== clustering: ${CLUSTERING}  ->  $(module_dir "${ARG:-sugarcane}")" ;;
  esac

  case "$STAGE" in

  # --- 01 VST export ----------------------------------------------------------
  export)
    check_study "$ARG"
    mkdir -p "${RESULTS}/${ARG}/vst"
    CLEAN_DDS="$(cfg DDS "$ARG")" \
    CLEAN_PREFIX="$(vst_prefix "$ARG")" \
    CLEAN_LABEL="$ARG" \
    CLEAN_COL_FILTER="$(cfg COLS "$ARG")" \
    CLEAN_STRIP_VERSION="$(cfg STRIP_VERSION "$ARG")" \
    MIN_CV="$MIN_CV" \
      "$RSCRIPT_DESEQ" "${SCRIPTS}/01_export_vst.r"
    ;;

  validate)
    "$PYTORCH" -u "${SCRIPTS}/validate.py" "${EXTRA[@]}"
    ;;

  # --- 02 one network layer ---------------------------------------------------
  network)
    check_study "$ARG"
    EST="${EXTRA[0]-}"; [ -n "$EST" ] || die "usage: run.sh network <study> pearson|ksg"
    EXTRA=("${EXTRA[@]:1}")
    [ "${FORCE:-0}" = "1" ] || assert_threshold \
      "$(layer_out "$ARG" "$EST").summary.json" "$MATCH_PEARSON" "--match-pearson"
    mkdir -p "${RESULTS}/${ARG}/layers"
    COMMON=(--matrix "$(vst_prefix "$ARG")"
            --out    "$(layer_out "$ARG" "$EST")"
            --estimator "$EST"
            --alpha "$ALPHA"
            --cand-pearson "$CAND_PEARSON"
            --gpu-mem-gb "$GPU_MEM_GB"
            --ram-limit-gb "$RAM_LIMIT_GB"
            --resume)
    case "$EST" in
      pearson)
        # The linear layer's floor IS the network threshold, exactly. Passing it
        # as --match-pearson rather than --min-value is deliberate: it makes the
        # engine invert its own analytic null and land back on STAT_MIN, which is
        # the self-check that the null is wired correctly.
        "$PYTORCH" -u "${SCRIPTS}/02_network_engine.py" "${COMMON[@]}" \
          --match-pearson "$STAT_MIN" --max-value "$STAT_MAX" "${EXTRA[@]}" ;;
      ksg)
        # The MI floor is set to the MI value with the SAME per-edge false-positive
        # rate as |r| = STAT_MIN at this n. Equal specificity is what licenses the
        # union in 03_merge_layers.py -- see docs/thresholds.md.
        "$PYTORCH" -u "${SCRIPTS}/02_network_engine.py" "${COMMON[@]}" \
          --k "$KSG_K" --match-pearson "$MATCH_PEARSON" \
          --max-perm "$MAX_PERM" "${EXTRA[@]}" ;;
      *) die "unknown estimator '$EST' (expected pearson or ksg)" ;;
    esac
    ;;

  # --- 03 merge the layers ----------------------------------------------------
  merge)
    check_study "$ARG"
    [ "${FORCE:-0}" = "1" ] || assert_threshold \
      "$(dirname "$(network_tsv "$ARG")")/network_${ARG}_edges.summary.json" \
      "$STAT_MIN" "stat_min"
    N=$(sed -n 's/.*"n_samples": *\([0-9]*\).*/\1/p' "$(vst_prefix "$ARG").meta.json")
    [ -n "$N" ] || die "could not read n_samples from $(vst_prefix "$ARG").meta.json"
    "$PYTORCH" -u "${SCRIPTS}/03_merge_layers.py" \
      --study "$ARG" --n "$N" \
      --pearson "$(layer_out "$ARG" pearson).edgelist.tsv" \
      --mi      "$(layer_out "$ARG" ksg).edgelist.tsv" \
      --genes   "$(vst_prefix "$ARG").genes.txt" \
      --out     "$(network_tsv "$ARG")" \
      --stat-min "$STAT_MIN" --stat-max "$STAT_MAX" "${EXTRA[@]}"
    ;;

  build)
    check_study "$ARG"
    "$0" export  "$ARG"
    "$0" network "$ARG" pearson
    "$0" network "$ARG" ksg
    "$0" merge   "$ARG"
    ;;

  # --- 04 topology ------------------------------------------------------------
  stats)
    check_study "$ARG"
    CLEAN_STUDY="$ARG" \
    CLEAN_EDGES="$(network_tsv "$ARG")" \
    CLEAN_PREFIX="$(study_dir "$ARG")/network_${ARG}" \
    CLEAN_TRANSITIVITY="$COMPUTE_TRANSITIVITY" \
    CLEAN_CORES="$NUM_CORES" \
      "$RSCRIPT_NET" "${SCRIPTS}/04_network_stats.r"
    ;;

  # --- choosing k by graph model (Pearson-only track) --------------------------
  # See docs/decisions.md. These write ONLY into whatever RESULTS points at, and
  # 40 refuses to run against a tree holding a merged network.
  pearsononly)
    check_study "$ARG"
    CLEAN_STUDY="$ARG" \
    CLEAN_PEARSON_LAYER="${CLEAN_PEARSON_LAYER:-${CLEAN}/results/$ARG/layers/${ARG}_pearson.edgelist.tsv}" \
    CLEAN_OUT="$(network_tsv "$ARG")" \
    CLEAN_STAT_MIN="$STAT_MIN" \
    CLEAN_STAT_MAX="$STAT_MAX" \
    CLEAN_FORCE="${EXTRA[0]:-0}" \
      bash "${SCRIPTS}/40_pearson_only_network.sh"
    ;;

  sgvalidate)
    CLEAN_RLIB="$CLEAN_RLIB" \
    CLEAN_OUT="${RESULTS}/statgraph_validation.tsv" \
    CLEAN_REAL_EDGES="${CLEAN_REAL_EDGES:-}" \
    CLEAN_A_N="${CLEAN_A_N:-1000,10000}" \
    CLEAN_BIG_N="${CLEAN_BIG_N:-50000}" \
      "$RSCRIPT_STATGRAPH" "${SCRIPTS}/39_statgraph_validate.r"
    ;;

  knnsweep)
    check_study "$ARG"
    CLEAN_STUDY="$ARG" \
    CLEAN_WORK_DIR="$MCL_WORK_DIR" \
    CLEAN_MCL_BIN_DIR="$(dirname "$MCL_BIN")" \
    CLEAN_OUT_DIR="$(study_dir "$ARG")" \
    CLEAN_K_LIST="$KNN_BA_K_LIST" \
    CLEAN_JOBS="$KNN_BA_JOBS" \
    CLEAN_CORES="$NUM_CORES" \
      bash "${SCRIPTS}/41_knn_ba_sweep.sh"
    ;;

  scoreref)
    check_study "$ARG"
    CLEAN_STUDY="$ARG" \
    CLEAN_WORK_DIR="$MCL_WORK_DIR" \
    CLEAN_MCL_BIN_DIR="$(dirname "$MCL_BIN")" \
    CLEAN_OUT_DIR="$(study_dir "$ARG")" \
    CLEAN_REF="${CLEAN_REF:-}" \
      bash "${SCRIPTS}/45_score_vs_reference.sh"
    ;;

  knncollect)
    check_study "$ARG"
    CLEAN_STUDY="$ARG" \
    CLEAN_OUT_DIR="$(study_dir "$ARG")" \
    CLEAN_VALIDATION="${RESULTS}/statgraph_validation.tsv" \
      "$RSCRIPT_STATGRAPH" "${SCRIPTS}/43_knn_ba_collect.r"
    ;;

  knnselect)
    check_study "$ARG"
    CLEAN_STUDY="$ARG" \
    CLEAN_RLIB="$CLEAN_RLIB" \
    CLEAN_GRID="$(study_dir "$ARG")/knnba_grid_${ARG}.tsv" \
    CLEAN_OUT_DIR="$(study_dir "$ARG")" \
    CLEAN_K="${EXTRA[0]:-}" \
    CLEAN_NPARAM="$KNN_BA_NPARAM" \
      "$RSCRIPT_STATGRAPH" "${SCRIPTS}/42_knn_model_selection.r"
    ;;

  # --- MCL toolchain: load once, survey degrees, sweep granularity -------------
  # The giant module is a node-degree problem, not an inflation one; see
  # docs/decisions.md. These stages read the network and write only diagnostics --
  # nothing published moves until `mcl` is re-run with the chosen settings.
  pearsonmci)
    check_study "$ARG"
    CLEAN_STUDY="$ARG" \
    CLEAN_LAYER="$(main_layer_out "$ARG" pearson).edgelist.tsv" \
    CLEAN_LAYER_SUMMARY="$(main_layer_out "$ARG" pearson).summary.json" \
    CLEAN_WORK_DIR="$MCL_WORK_DIR" \
    CLEAN_MCL_BIN_DIR="$(dirname "$MCL_BIN")" \
    CLEAN_FORCE="${EXTRA[0]:-0}" \
    CLEAN_CORES="$NUM_CORES" \
      bash "${SCRIPTS}/46_pearson_mci.sh"
    ;;

  mclload)
    check_study "$ARG"
    CLEAN_STUDY="$ARG" \
    CLEAN_EDGES="$(network_tsv "$ARG")" \
    CLEAN_WORK_DIR="$MCL_WORK_DIR" \
    CLEAN_MCL_BIN_DIR="$(dirname "$MCL_BIN")" \
    CLEAN_FORCE="${EXTRA[0]:-0}" \
    CLEAN_CORES="$NUM_CORES" \
      bash "${SCRIPTS}/34_mcl_load.sh"
    ;;

  mclsurvey)
    check_study "$ARG"
    CLEAN_STUDY="$ARG" \
    CLEAN_WORK_DIR="$MCL_WORK_DIR" \
    CLEAN_MCL_BIN_DIR="$(dirname "$MCL_BIN")" \
    CLEAN_OUT_DIR="$(study_dir "$ARG")" \
    CLEAN_KNN_RANGE="$MCL_SWEEP_KNN_RANGE" \
    CLEAN_NODE_METRICS="$(study_dir "$ARG")/network_${ARG}_node_metrics.tsv" \
    CLEAN_FORCE="${EXTRA[0]:-0}" \
    CLEAN_CORES="$NUM_CORES" \
      bash "${SCRIPTS}/35_mcl_survey.sh"
    ;;

  mclsweep)
    check_study "$ARG"
    CLEAN_STUDY="$ARG" \
    CLEAN_WORK_DIR="$MCL_WORK_DIR" \
    CLEAN_MCL_BIN_DIR="$(dirname "$MCL_BIN")" \
    CLEAN_OUT_DIR="$(study_dir "$ARG")" \
    CLEAN_SWEEP_I="$MCL_SWEEP_I" \
    CLEAN_SWEEP_I_NONE="$(cfg MCL_SWEEP_I_NONE "$ARG")" \
    CLEAN_SWEEP_K="${EXTRA[0]:-}" \
    CLEAN_CORES="$NUM_CORES" \
      bash "${SCRIPTS}/36_mcl_sweep.sh"
    ;;

  mclladder)
    # The inflation sweep with the k-NN axis switched OFF: 36 loops "for K in
    # none" when CLEAN_SWEEP_K is empty, so this is a pure inflation ladder on
    # the unpruned matrix. Second argument is an mcl resource string for the
    # pruning probe, e.g.  ./run.sh mclladder sugarcane "-S 10000"
    check_study "$ARG"
    CLEAN_STUDY="$ARG" \
    CLEAN_WORK_DIR="$MCL_WORK_DIR" \
    CLEAN_MCL_BIN_DIR="$(dirname "$MCL_BIN")" \
    CLEAN_OUT_DIR="$(study_dir "$ARG")" \
    CLEAN_SWEEP_I_NONE="${CLEAN_SWEEP_I_NONE:-$CLUSTER_I_LIST}" \
    CLEAN_SWEEP_K="" \
    CLEAN_MCL_RESOURCE="${EXTRA[0]:-$CLUSTER_MCL_RESOURCE}" \
    CLEAN_CORES="$NUM_CORES" \
      bash "${SCRIPTS}/36_mcl_sweep.sh"
    ;;

  leidensweep)
    # Second argument "scout" runs only the calibration gammas.
    check_study "$ARG"
    CLEAN_STUDY="$ARG" \
    CLEAN_WORK_DIR="$MCL_WORK_DIR" \
    CLEAN_MCL_BIN_DIR="$(dirname "$MCL_BIN")" \
    CLEAN_OUT_DIR="$(study_dir "$ARG")" \
    CLEAN_LEIDEN_MODE="${EXTRA[0]:-full}" \
    CLUSTER_LEIDEN_GAMMA="$CLUSTER_LEIDEN_GAMMA" \
    CLUSTER_LEIDEN_SCOUT="$CLUSTER_LEIDEN_SCOUT" \
    CLUSTER_LEIDEN_ITER="$CLUSTER_LEIDEN_ITER" \
      "$CLUSTER_PYTHON" "${SCRIPTS}/47_leiden_sweep.py"
    ;;

  clustercompare)
    check_study "$ARG"
    CLEAN_STUDY="$ARG" \
    CLEAN_WORK_DIR="$MCL_WORK_DIR" \
    CLEAN_MCL_BIN_DIR="$(dirname "$MCL_BIN")" \
    CLEAN_OUT_DIR="$(study_dir "$ARG")" \
    CLEAN_REF="${CLEAN_REF:-}" \
      bash "${SCRIPTS}/48_cluster_compare.sh"
    ;;

  fastgreedy)
    # The hierarchical family, bounded. Classical hierarchical clustering is not
    # offered at all: it needs a dense 170,135^2 dissimilarity (232 GB) built
    # from correlations this pipeline thresholded away.
    check_study "$ARG"
    CLEAN_STUDY="$ARG" \
    CLEAN_WORK_DIR="$MCL_WORK_DIR" \
    CLUSTER_FASTGREEDY_CAP_S="$CLUSTER_FASTGREEDY_CAP_S" \
      "$CLUSTER_PYTHON" "${SCRIPTS}/50_fastgreedy_scout.py"
    ;;

  figclustermethods)
    CLEAN_RESULTS="$RESULTS" \
    CLEAN_OUT_DIR="${RESULTS}/figures" \
    CLEAN_STUDIES="sugarcane purple" \
      "$RSCRIPT_PLOT" "${SCRIPTS}/49_fig_cluster_methods.r"
    ;;

  clusterhomog)
    check_study "$ARG"
    CLEAN_STUDY="$ARG" \
    CLEAN_WORK_DIR="$MCL_WORK_DIR" \
    CLEAN_MCL_BIN_DIR="$(dirname "$MCL_BIN")" \
    CLEAN_OUT_DIR="$(study_dir "$ARG")" \
    CLEAN_EMAPPER="$(cfg EMAPPER "$ARG")" \
    CLEAN_MAX_GENES="$HOMOGENEITY_MAX_GENES" \
    CLEAN_PERM="$HOMOGENEITY_PERM" \
    CLEAN_SEED="$HOMOGENEITY_SEED" \
    CLEAN_CORES="$NUM_CORES" \
      "$RSCRIPT_NET" "${SCRIPTS}/37_cluster_homogeneity.r"
    ;;

  figclusterchoice)
    CLEAN_RESULTS="$RESULTS" \
    CLEAN_OUT_DIR="${RESULTS}/figures" \
    CLEAN_STUDIES="sugarcane purple" \
    CLEAN_WORK_DIR="$MCL_WORK_DIR" \
    CLEAN_TREES="${CLEAN_TREES:-}" \
      "$RSCRIPT_PLOT" "${SCRIPTS}/38_fig_clustering_choice.r"
    ;;

  # --- 05 clustering ----------------------------------------------------------
  nodemetrics)
    # Gene universe for the Pearson-only network, straight from the .mci.
    check_study "$ARG"
    CLEAN_STUDY="$ARG" \
    CLEAN_WORK_DIR="$CLUSTER_WORK_DIR" \
    CLEAN_MCL_BIN_DIR="$(dirname "$MCL_BIN")" \
    CLEAN_OUT_FILE="$(study_dir "$ARG")/network_${ARG}_node_metrics.tsv" \
    CLEAN_CORES="$NUM_CORES" \
      bash "${SCRIPTS}/52_pearson_node_metrics.sh"
    ;;

  membership)
    # Adopt the inflation ladder's partition at this study's chosen inflation.
    # Replaces `./run.sh mcl <study>` for the Pearson-only networks: 05 re-runs
    # mcl from a 9-column edge table these networks no longer have, and the
    # partition is already on disk.
    check_study "$ARG"
    _I="$(cfg MCL_INFLATION "$ARG")"; _I="${_I:-$MCL_INFLATION}"
    CLEAN_STUDY="$ARG" \
    CLEAN_CLS="$(cls_for "$ARG" "$_I")" \
    CLEAN_TAB="${CLUSTER_WORK_DIR}/${ARG}.tab" \
    CLEAN_NODE_METRICS="$(study_dir "$ARG")/network_${ARG}_node_metrics.tsv" \
    CLEAN_SWEEP_TSV="${CLUSTER_SWEEP_TREE}/${ARG}/mcl_sweep_${ARG}.tsv" \
    CLEAN_PREFIX="$(clus_prefix "$ARG")" \
    CLEAN_INFLATION="$_I" \
    CLEAN_MIN_MODULE_SIZE="$MCL_MIN_MODULE_SIZE" \
    CLEAN_MIN_MODULE_SIZE_PLOT="$MCL_MIN_MODULE_SIZE_PLOT" \
    CLEAN_CORES="$NUM_CORES" \
      "$RSCRIPT_NET" "${SCRIPTS}/51_mcl_membership_from_cls.r"
    ;;

  mcl)
    check_study "$ARG"
      # mcl ships inside the same conda env as RSCRIPT_NET but that env is not
    # activated, so its bin/ is not on PATH and the R script's system2("mcl")
    # would fail. Prepend it rather than requiring the caller to activate.
    export PATH="$(dirname "$MCL_BIN"):$PATH"
    command -v mcl >/dev/null || die "mcl not found at $MCL_BIN (set MCL_BIN in config.sh)"
    mkdir -p "$MCL_TMPDIR"
    CLEAN_STUDY="$ARG" \
    CLEAN_EDGES="$(network_tsv "$ARG")" \
    CLEAN_PREFIX="$(study_dir "$ARG")/mcl_${ARG}" \
    CLEAN_INFLATION="$MCL_INFLATION" \
    CLEAN_MIN_MODULE_SIZE="$MCL_MIN_MODULE_SIZE" \
    CLEAN_MIN_MODULE_SIZE_PLOT="$MCL_MIN_MODULE_SIZE_PLOT" \
    CLEAN_SAVE_GRAPH_RDS="$MCL_SAVE_GRAPH_RDS" \
    CLEAN_CORES="$NUM_CORES" \
    TMPDIR="$MCL_TMPDIR" \
      "$RSCRIPT_NET" "${SCRIPTS}/05_mcl_clustering.r"
    ;;

  # Converts a graph-tool nested SBM fit into the two files every module-level
  # stage reads, in MCL's exact schema, so the alternative clustering can be
  # carried through the identical downstream analysis. Writes beside the MCL
  # files; nothing is overwritten.
  sbmclust)
    check_study "$ARG"
    SBM_DIR="$(cfg SBM_DIR "$ARG")"
    [ -d "$SBM_DIR" ] || die "no SBM output for '$ARG' at $SBM_DIR"
    CLEAN_STUDY="$ARG" \
    CLEAN_SBM_NODE_BLOCKS="${SBM_DIR}/${SBM_TAG}_node_blocks.tsv" \
    CLEAN_NODE_METRICS="$(study_dir "$ARG")/network_${ARG}_node_metrics.tsv" \
    CLEAN_EDGES="$(network_tsv "$ARG")" \
    CLEAN_PREFIX="$(study_dir "$ARG")/sbm_${ARG}" \
    CLEAN_SBM_LEVEL="$SBM_LEVEL" \
    CLEAN_SBM_COMPUTE_Q="$SBM_COMPUTE_Q" \
    CLEAN_MIN_MODULE_SIZE="$MCL_MIN_MODULE_SIZE" \
    CLEAN_MIN_MODULE_SIZE_PLOT="$MCL_MIN_MODULE_SIZE_PLOT" \
    CLEAN_CORES="$NUM_CORES" \
      "$RSCRIPT_NET" "${SCRIPTS}/27_sbm_membership.r"
    ;;

  # --- 06 conservation --------------------------------------------------------
  conserve)
    case "$ARG" in
      sugarcane_to_purple|purple_to_sugarcane) ;;
      *) die "usage: run.sh conserve sugarcane_to_purple|purple_to_sugarcane" ;;
    esac
    CLEAN_DIRECTION="$ARG" \
    CLEAN_EDGES_SUGARCANE="$(network_tsv sugarcane)" \
    CLEAN_EDGES_PURPLE="$(network_tsv purple)" \
    CLEAN_OUT_DIR="${RESULTS}/conservation" \
    CLEAN_ORTHOGROUPS="$ORTHOGROUPS" \
    CLEAN_CHUNK_SIZE="$CHUNK_SIZE" \
    CLEAN_CORES="$NUM_CORES" \
      "$RSCRIPT_NET" "${SCRIPTS}/06_conservation_join.r"
    ;;

  conservenull)
    case "$ARG" in
      sugarcane_to_purple|purple_to_sugarcane) ;;
      *) die "usage: run.sh conservenull sugarcane_to_purple|purple_to_sugarcane" ;;
    esac
    CLEAN_DIRECTION="$ARG" \
    CLEAN_EDGES_SUGARCANE="$(network_tsv sugarcane)" \
    CLEAN_EDGES_PURPLE="$(network_tsv purple)" \
    CLEAN_OUT_DIR="${RESULTS}/conservation" \
    CLEAN_ORTHOGROUPS="$ORTHOGROUPS" \
    CLEAN_NULL_REPS="$NULL_REPS" \
    CLEAN_NULL_SAMPLE="$NULL_SAMPLE" \
    CLEAN_CORES="$NUM_CORES" \
      "$RSCRIPT_NET" "${SCRIPTS}/13_conservation_null.r"
    ;;

  # --- 07-08 trait correlation ------------------------------------------------
  trait)
    check_study "$ARG"
    CLEAN_STUDY="$ARG" \
    CLEAN_VST_PREFIX="$(vst_prefix "$ARG")" \
    CLEAN_META="$(cfg META "$ARG")" \
    CLEAN_TRAITS="$(cfg TRAITS "$ARG")" \
    CLEAN_OUT_DIR="$(study_dir "$ARG")" \
    CLEAN_GENE_FILTER="${RESULTS}/conservation/conserved_genes_${ARG}_FULL.txt" \
    CLEAN_SELECT_TRAIT="$SELECT_TRAIT" \
    CLEAN_TRAIT_R_THR="$TRAIT_R_THR" \
    CLEAN_TRAIT_PADJ_THR="$TRAIT_PADJ_THR" \
    CLEAN_CORES="$NUM_CORES" \
      "$RSCRIPT_NET" "${SCRIPTS}/07_gene_trait_cor.r"
    ;;

  traitmi)
    check_study "$ARG"
    "$PYTORCH" -u "${SCRIPTS}/12_gene_trait_mi.py" \
      --matrix "$(vst_prefix "$ARG")" \
      --meta   "$(cfg META "$ARG")" \
      --traits "$(cfg TRAITS "$ARG")" \
      --trait  "$SELECT_TRAIT" \
      --out    "$(study_dir "$ARG")/gene_trait_mi_${ARG}" \
      --k "$TRAIT_MI_K" --n-perm "$TRAIT_MI_PERM" --alpha "$PADJ_THR" \
      "${EXTRA[@]}"
    ;;

  # Purple's design is stress-control-stress, so a gene induced by BOTH nitrogen
  # deficiency and excess is invisible to every monotone test in this pipeline.
  # This is the only genome-wide gene-level test of that shape. Feeds
  # `./run.sh conscor 0 ushape` and `./run.sh conscor 1 ushape`.
  ushape)
    S="${ARG:-purple}"; check_study "$S"
    CLEAN_STUDY="$S" \
    CLEAN_VST_PREFIX="$(vst_prefix "$S")" \
    CLEAN_META="$(cfg META "$S")" \
    CLEAN_OUT_FILE="${RESULTS}/${S}/gene_trait_ushape_${S}.tsv" \
    CLEAN_GENE_FILTER="$(node_list "$S")" \
    CLEAN_SELECT_TRAIT="$SELECT_TRAIT" \
      "$RSCRIPT_NET" "${SCRIPTS}/30_gene_trait_ushape.r"
    ;;

  # The gene-level test with the design in the model -- genotype (and, in
  # sugarcane, leaf segment) as blocks instead of unmodelled residual variance.
  # Same FDR level and same |r| floor as the marginal rule; only the model changes,
  # and 31 computes the marginal fit on the same genes so the two are comparable.
  #
  # The universe is the PEARSON-ONLY NETWORK'S NODE SET, the gene set the modules
  # are built on -- not the merged graph's conserved-edge set the first version of
  # this stage used. Feeds `./run.sh consblocked`.
  traitblocked)
    check_study "$ARG"
    CLEAN_STUDY="$ARG" \
    CLEAN_VST_PREFIX="$(vst_prefix "$ARG")" \
    CLEAN_META="$(cfg META "$ARG")" \
    CLEAN_TRAITS="$(cfg TRAITS "$ARG")" \
    CLEAN_OUT_FILE="${RESULTS}/${ARG}/gene_trait_blocked_${ARG}.tsv" \
    CLEAN_GENE_FILTER="$(node_list "$ARG")" \
    CLEAN_SELECT_TRAIT="$SELECT_TRAIT" \
    CLEAN_BLOCK="$(cfg TRAIT_BLOCK "$ARG")" \
    CLEAN_BLOCKED_STAT="$(cfg TRAIT_BLOCKED_STAT "$ARG")" \
    CLEAN_PLANT_FROM="$(cfg TRAIT_PLANT_FROM "$ARG")" \
    CLEAN_TRAIT_R_THR="$TRAIT_R_THR" \
    CLEAN_TRAIT_PADJ_THR="$TRAIT_PADJ_THR" \
    CLEAN_CORES="$NUM_CORES" \
      "$RSCRIPT_NET" "${SCRIPTS}/31_gene_trait_blocked.r"
    ;;

  # Which ortholog genes are nitrogen-correlated in BOTH species -- NODE level,
  # blocked rule only, on the Pearson-only networks. ARG=1 runs the DIRECTED
  # variant (discover in sugarcane, re-correct purple over the candidate orthologs).
  #
  # It does NOT supersede `conscor`, which is the record of the merged network under
  # the marginal pearson|mi|union rules and still carries the edge level. This has
  # no edge level: that would need the conserved-edge tables, which describe the
  # merged graph.
  consblocked)
    case "${ARG:-0}" in 0|1) ;; *) die "usage: run.sh consblocked [0|1]" ;; esac
    CLEAN_RESULTS="$RESULTS" \
    CLEAN_OUT_DIR="${RESULTS}/conservation" \
    CLEAN_ORTHOGROUPS="$ORTHOGROUPS" \
    CLEAN_NODES_SUGARCANE="$(node_list sugarcane)" \
    CLEAN_NODES_PURPLE="$(node_list purple)" \
    CLEAN_SELECT_TRAIT="$SELECT_TRAIT" \
    CLEAN_DIRECTED="${ARG:-0}" \
    CLEAN_DISCOVERY="$TRAIT_DISCOVERY" \
    CLEAN_USHAPE_TIER="$TRAIT_USHAPE_TIER" \
    CLEAN_TRAIT_R_THR="$TRAIT_R_THR" \
    CLEAN_TRAIT_PADJ_THR="$TRAIT_PADJ_THR" \
    CLEAN_NULL_REPS="$TRAIT_NULL_REPS" \
    CLEAN_SEED="$TRAIT_NULL_SEED" \
    CLEAN_CORES="$NUM_CORES" \
      "$RSCRIPT_NET" "${SCRIPTS}/61_conserved_blocked_nodes.r"
    ;;

  conscor)
    CLEAN_RESULTS="$RESULTS" \
    CLEAN_OUT_DIR="${RESULTS}/conservation" \
    CLEAN_ORTHOGROUPS="$ORTHOGROUPS" \
    CLEAN_SELECT_TRAIT="$SELECT_TRAIT" \
    CLEAN_SELECTION="${EXTRA[0]:-$TRAIT_SELECTION}" \
    CLEAN_DIRECTED="${ARG:-0}" \
    CLEAN_DISCOVERY="$TRAIT_DISCOVERY" \
    CLEAN_TRAIT_R_THR="$TRAIT_R_THR" \
    CLEAN_TRAIT_PADJ_THR="$TRAIT_PADJ_THR" \
    CLEAN_CORES="$NUM_CORES" \
      "$RSCRIPT_NET" "${SCRIPTS}/08_conserved_cor_genes.r"
    ;;

  # --- 14-16 module-level analysis --------------------------------------------
  eigengene)
    check_study "$ARG"
    CLEAN_STUDY="$ARG" \
    CLEAN_VST_PREFIX="$(vst_prefix "$ARG")" \
    CLEAN_MEMBERSHIP="$(clus_prefix "$ARG")_membership.tsv" \
    CLEAN_OUT_PREFIX="$(module_dir "$ARG")/modules/${ARG}_eigengenes" \
    CLEAN_MIN_MODULE_SIZE_EIGEN="$MIN_MODULE_SIZE_EIGEN" \
    CLEAN_CORES="$NUM_CORES" \
      "$RSCRIPT_NET" "${SCRIPTS}/14_module_eigengene.r"
    ;;

  # Spearman only. The module response used to be called by Pearson AND mutual
  # information (12_gene_trait_mi.py run on the eigengene matrix); it is now one
  # rank correlation, because the trait is ordinal and MI at this level was an
  # omnibus test firing on sample-driven quirks. See 19_module_trait_spearman.r.
  moduletrait)
    # BLOCKED (53), not the marginal Spearman (19). 53 also computes the
    # MARGINAL fit on the same eigengenes and carries it in the same table, so
    # the effect of putting the design in the model is a comparison of two models
    # of ONE module set -- never a join against a previous run, whose Module_NNN
    # names are positional and refer to different genes.
    check_study "$ARG"
    CLEAN_STUDY="$ARG" \
    CLEAN_EIGENGENE_PREFIX="$(module_dir "$ARG")/modules/${ARG}_eigengenes" \
    CLEAN_META="$(cfg META "$ARG")" \
    CLEAN_TRAITS="$(cfg TRAITS "$ARG")" \
    CLEAN_SELECT_TRAIT="${SELECT_TRAIT:-treatment}" \
    CLEAN_BLOCK="$(cfg MODULE_BLOCK "$ARG")" \
    CLEAN_BLOCKED_STAT="$(cfg MODULE_BLOCKED_STAT "$ARG")" \
    CLEAN_PLANT_FROM="$(cfg MODULE_PLANT_FROM "$ARG")" \
    CLEAN_OUT_FILE="$(module_dir "$ARG")/module_trait_${ARG}.tsv" \
    CLEAN_MODULE_R_THR="$MODULE_R_THR" \
    CLEAN_MODULE_PADJ_THR="$MODULE_PADJ_THR" \
    CLEAN_MODULE_PERM="$MODULE_PERM" \
    CLEAN_SEED="$MODULE_SEED" \
    CLEAN_CORES="$NUM_CORES" \
      "$RSCRIPT_NET" "${SCRIPTS}/53_module_trait_blocked.r"
    ;;

  moduleprofile)
    check_study "$ARG"
    CLEAN_STUDY="$ARG" \
    CLEAN_MODULE_TRAIT="$(module_dir "$ARG")/module_trait_${ARG}.tsv" \
    CLEAN_MEMBERSHIP="$(clus_prefix "$ARG")_membership.tsv" \
    CLEAN_PC1_VARIANCE="$(module_dir "$ARG")/modules/${ARG}_eigengenes_pc1_variance.tsv" \
    CLEAN_TF_FILE="${RESULTS}/readouts/get_tfs/${ARG}/TF_in_network.tsv" \
    CLEAN_NODE_METRICS="$(study_dir "$ARG")/network_${ARG}_node_metrics.tsv" \
    CLEAN_OUT_FILE="$(module_dir "$ARG")/module_profile_${ARG}.tsv" \
    CLEAN_PADJ_THR="$MODULE_PADJ_THR" \
    CLEAN_CORES="$NUM_CORES" \
      "$RSCRIPT_NET" "${SCRIPTS}/15_module_profile.r"
    ;;

  moduleheatmap)
    check_study "$ARG"
    CLEAN_STUDY="$ARG" \
    CLEAN_VST_PREFIX="$(vst_prefix "$ARG")" \
    CLEAN_MODULE_PROFILE="$(module_dir "$ARG")/module_profile_${ARG}.tsv" \
    CLEAN_MEMBERSHIP="$(clus_prefix "$ARG")_membership.tsv" \
    CLEAN_META="$(cfg META "$ARG")" \
    CLEAN_TRAITS="$(cfg TRAITS "$ARG")" \
    CLEAN_TF_FILE="${RESULTS}/readouts/get_tfs/${ARG}/TF_in_network.tsv" \
    CLEAN_OUT_DIR="$(module_dir "$ARG")/heatmaps" \
    CLEAN_HEATMAP_TOP_N="$HEATMAP_TOP_N" \
    CLEAN_HEATMAP_MAX_GENES="$HEATMAP_MAX_GENES" \
    CLEAN_HEATMAP_MODULES="${EXTRA[0]:-}" \
    CLEAN_HEATMAP_GROUP_BY="$HEATMAP_GROUP_BY" \
    CLEAN_CORES="$NUM_CORES" \
      "$RSCRIPT_PLOT" "${SCRIPTS}/16_module_heatmaps.r"
    ;;

  modulesummary)
    check_study "$ARG"
    CLEAN_STUDY="$ARG" \
    CLEAN_EIGENGENE_PREFIX="$(module_dir "$ARG")/modules/${ARG}_eigengenes" \
    CLEAN_MODULE_PROFILE="$(module_dir "$ARG")/module_profile_${ARG}.tsv" \
    CLEAN_META="$(cfg META "$ARG")" \
    CLEAN_TRAITS="$(cfg TRAITS "$ARG")" \
    CLEAN_OUT_PREFIX="$(module_dir "$ARG")/module_summary_${ARG}" \
    CLEAN_SUMMARY_MAX_MODULES="$SUMMARY_MAX_MODULES" \
    CLEAN_MODULE_R_THR="$MODULE_R_THR" \
    CLEAN_MODULE_PADJ_THR="$MODULE_PADJ_THR" \
    CLEAN_HEATMAP_GROUP_BY="$HEATMAP_GROUP_BY" \
    CLEAN_CORES="$NUM_CORES" \
      "$RSCRIPT_PLOT" "${SCRIPTS}/17_module_summary.r"
    ;;

  # topGO env, like `go` -- the module gene sets and the conserved gene set are
  # tested against the same network-node background, so the two enrichments stay
  # on one denominator.
  gene2go)
    # Derive GO from the protein-annotator output (InterPro + Pfam).
    check_study "$ARG"
    CLEAN_STUDY="$ARG" \
    CLEAN_ANNOT_DIR="$PROTEIN_ANNOT_DIR" \
    CLEAN_GO_MAP_DIR="$GO_MAP_DIR" \
    CLEAN_OUT_FILE="$(gene2go_tsv "$ARG")" \
    CLEAN_NODE_METRICS="$(study_dir "$ARG")/network_${ARG}_node_metrics.tsv" \
    CLEAN_PFAM_EVALUE="$PFAM_EVALUE" \
      bash "${SCRIPTS}/54_build_gene2go.sh"
    ;;

  modulego)
    check_study "$ARG"
    CLEAN_STUDY="$ARG" \
    CLEAN_MODULE_PROFILE="$(module_dir "$ARG")/module_profile_${ARG}.tsv" \
    CLEAN_MEMBERSHIP="$(clus_prefix "$ARG")_membership.tsv" \
    CLEAN_NODE_METRICS="$(study_dir "$ARG")/network_${ARG}_node_metrics.tsv" \
    CLEAN_EMAPPER="$(cfg EMAPPER "$ARG")" \
    CLEAN_GENE2GO="${CLEAN_GENE2GO-$(gene2go_tsv "$ARG")}" \
    CLEAN_OUT_DIR="$(module_dir "$ARG")/module_go" \
    CLEAN_ONTOLOGY="${EXTRA[0]:-$MODULE_GO_ONTOLOGY}" \
    CLEAN_GO_P="$GO_P" \
    CLEAN_GO_NTOP="$GO_NTOP" \
    CLEAN_MODULE_GO_MIN_ANNOTATED="$MODULE_GO_MIN_ANNOTATED" \
    CLEAN_MODULE_GO_CORES="$MODULE_GO_CORES" \
    CLEAN_MODULE_GO_LIMIT="${CLEAN_MODULE_GO_LIMIT:-0}" \
    CLEAN_CORES="$NUM_CORES" \
      conda run --no-capture-output -n "$CONDA_TOPGO" \
        Rscript "${SCRIPTS}/18_module_go.r"
    ;;

  # --- 20-21 paper figures ------------------------------------------------------
  # Scripts are named for what they draw, not for their figure number: figure
  # ORDER is editorial and has already moved once. FIG_* in config.sh carries the
  # number, and it names the output files and opens the generated legend.
  figdataset)
    CLEAN_STUDIES="$STUDIES" \
    CLEAN_FIG_NUM="$FIG_DATASET" \
    CLEAN_LABEL_SUGARCANE="sugarcane" \
    CLEAN_LABEL_PURPLE="purple" \
    CLEAN_META_SUGARCANE="$META_sugarcane" \
    CLEAN_META_PURPLE="$META_purple" \
    CLEAN_SALMONQC_SUGARCANE="$SALMONQC_sugarcane" \
    CLEAN_SALMONQC_PURPLE="$SALMONQC_purple" \
    CLEAN_TPM_SUGARCANE="$TPM_sugarcane" \
    CLEAN_TPM_PURPLE="$TPM_purple" \
    CLEAN_VST_SUGARCANE="$(vst_prefix sugarcane)" \
    CLEAN_VST_PURPLE="$(vst_prefix purple)" \
    CLEAN_NODES_SUGARCANE="$(study_dir sugarcane)/network_sugarcane_node_metrics.tsv" \
    CLEAN_NODES_PURPLE="$(study_dir purple)/network_purple_node_metrics.tsv" \
    CLEAN_PCA_NTOP="$PCA_NTOP" \
    CLEAN_OUT_PREFIX="${RESULTS}/figures/figure${FIG_DATASET}_dataset_qc" \
    CLEAN_CORES="$NUM_CORES" \
      "$RSCRIPT_PLOT" "${SCRIPTS}/21_fig_dataset_qc.r"
    ;;

  figtopology)
    CLEAN_STUDIES="$STUDIES" \
    CLEAN_FIG_NUM="$FIG_TOPOLOGY" \
    CLEAN_NODES_SUGARCANE="$(study_dir sugarcane)/network_sugarcane_node_metrics.tsv" \
    CLEAN_NODES_PURPLE="$(study_dir purple)/network_purple_node_metrics.tsv" \
    CLEAN_GLOBAL_SUGARCANE="$(study_dir sugarcane)/network_sugarcane_global_metrics.tsv" \
    CLEAN_GLOBAL_PURPLE="$(study_dir purple)/network_purple_global_metrics.tsv" \
    CLEAN_MCL_SUGARCANE="$(clus_prefix sugarcane)_module_summary.tsv" \
    CLEAN_MCL_PURPLE="$(clus_prefix purple)_module_summary.tsv" \
    CLEAN_TOPO_GRID="$TOPO_GRID" \
    CLEAN_OUT_PREFIX="${RESULTS}/figures/figure${FIG_TOPOLOGY}_topology" \
    CLEAN_CORES="$NUM_CORES" \
      "$RSCRIPT_PLOT" "${SCRIPTS}/22_fig_topology.r"
    ;;

  figconservation)
    # topGO's GenTable truncates long term names to 40 characters with an
    # ellipsis, and 18 of the 69 shared BP terms come out cut. GO.db has the full
    # names but lives in topGO_env, not the plotting env, so the id -> term map is
    # dumped once here and cached; the figure expands the truncated strings from it.
    GO_NAMES="${RESULTS}/conservation/go_term_names.tsv"
    if [ ! -s "$GO_NAMES" ]; then
      echo "building GO term-name cache -> ${GO_NAMES##*/}"
      conda run --no-capture-output -n "$CONDA_TOPGO" Rscript -e "
        suppressMessages(library(GO.db))
        tt <- AnnotationDbi::Term(GOTERM)
        utils::write.table(data.frame(GO.ID = names(tt), Term_full = unname(tt)),
                           '$GO_NAMES', sep = '\t', quote = FALSE, row.names = FALSE)"
    fi
    CLEAN_GO_NAMES="$GO_NAMES" \
    CLEAN_GO_ONTOLOGY="$CONS_GO_ONTOLOGY" \
    CLEAN_GO_NTERMS="$CONS_GO_NTERMS" \
    CLEAN_FIG_NUM="$FIG_CONSERVATION" \
    CLEAN_CONS_DIR="${RESULTS}/conservation" \
    CLEAN_SELECTION="$TRAIT_SELECTION" \
    CLEAN_NODES_SUGARCANE="$(study_dir sugarcane)/network_sugarcane_node_metrics.tsv" \
    CLEAN_NODES_PURPLE="$(study_dir purple)/network_purple_node_metrics.tsv" \
    CLEAN_OUT_PREFIX="${RESULTS}/figures/figure${FIG_CONSERVATION}_conservation" \
    CLEAN_CORES="$NUM_CORES" \
      "$RSCRIPT_PLOT" "${SCRIPTS}/23_fig_conservation.r"
    ;;

  figmodules)
    CLEAN_STUDIES="$STUDIES" \
    CLEAN_FIG_NUM="$FIG_MODULES" \
    CLEAN_MODULE_R_THR="$MODULE_R_THR" \
    CLEAN_MODULE_PADJ_THR="$MODULE_PADJ_THR" \
    CLEAN_SELECT_TRAIT="$SELECT_TRAIT" \
    CLEAN_PROFILE_SUGARCANE="$(study_dir sugarcane)/module_profile_sugarcane.tsv" \
    CLEAN_PROFILE_PURPLE="$(study_dir purple)/module_profile_purple.tsv" \
    CLEAN_NULL_SUGARCANE="$(study_dir sugarcane)/module_trait_sugarcane.null.tsv" \
    CLEAN_NULL_PURPLE="$(study_dir purple)/module_trait_purple.null.tsv" \
    CLEAN_EIGEN_SUGARCANE="$(study_dir sugarcane)/modules/sugarcane_eigengenes" \
    CLEAN_EIGEN_PURPLE="$(study_dir purple)/modules/purple_eigengenes" \
    CLEAN_META_SUGARCANE="$META_sugarcane" \
    CLEAN_META_PURPLE="$META_purple" \
    CLEAN_TRAITS_SUGARCANE="$TRAITS_sugarcane" \
    CLEAN_TRAITS_PURPLE="$TRAITS_purple" \
    CLEAN_OUT_PREFIX="${RESULTS}/figures/figure${FIG_MODULES}_modules" \
    CLEAN_CORES="$NUM_CORES" \
      "$RSCRIPT_PLOT" "${SCRIPTS}/24_fig_modules.r"
    ;;

  figmodulego)
    CLEAN_STUDIES="$STUDIES" \
    CLEAN_FIG_NUM="$FIG_MODULE_GO" \
    CLEAN_GO_ONTOLOGY="$MODULE_GO_ONTOLOGY" \
    CLEAN_GO_NTERMS="$MODULE_GO_FIG_NTERMS" \
    CLEAN_GO_NRECUR="$MODULE_GO_FIG_NRECUR" \
    CLEAN_GO_MAIN_STUDY="$MODULE_GO_FIG_MAIN" \
    CLEAN_GODIR_SUGARCANE="$(study_dir sugarcane)/module_go" \
    CLEAN_GODIR_PURPLE="$(study_dir purple)/module_go" \
    CLEAN_OUT_PREFIX="${RESULTS}/figures/figure${FIG_MODULE_GO}_module_go" \
    CLEAN_CORES="$NUM_CORES" \
      "$RSCRIPT_PLOT" "${SCRIPTS}/25_fig_module_go.r"
    ;;

  figmodule20)
    CLEAN_FIG_NUM="$FIG_MODULE20" \
    CLEAN_M20_DIR="${RESULTS}/readouts/module20" \
    CLEAN_TPM_SUGARCANE="$TPM_sugarcane" \
    CLEAN_TPM_PURPLE="$TPM_purple" \
    CLEAN_META_SUGARCANE="$META_sugarcane" \
    CLEAN_META_PURPLE="$META_purple" \
    CLEAN_M20_FOCUS="$M20_FOCUS_GENE" \
    CLEAN_M20_ANCHOR="$M20_FOCUS_ANCHOR" \
    CLEAN_OUT_PREFIX="${RESULTS}/figures/figure${FIG_MODULE20}_module20" \
    CLEAN_CORES="$NUM_CORES" \
      "$RSCRIPT_PLOT" "${SCRIPTS}/26_fig_module20.r"
    ;;

  figclustering)
    CLEAN_STUDY="${ARG:-sugarcane}" \
    CLEAN_STUDY_DIR="$(study_dir "${ARG:-sugarcane}")" \
    CLEAN_FIG_NUM="$FIG_CLUSTERING" \
    CLEAN_GO_ONTOLOGY="$MODULE_GO_ONTOLOGY" \
    CLEAN_TOPO_GRID="$TOPO_GRID" \
    CLEAN_OUT_PREFIX="${RESULTS}/figures/figure${FIG_CLUSTERING}_clustering_${ARG:-sugarcane}" \
    CLEAN_CORES="$NUM_CORES" \
      "$RSCRIPT_PLOT" "${SCRIPTS}/28_fig_clustering_compare.r"
    ;;

  figrepro)
    CLEAN_BASE="$BASE" \
    CLEAN_FIG_NUM="$FIG_REPRODUCTION" \
    CLEAN_READOUTS="${RESULTS}/readouts" \
    CLEAN_TPM_SUGARCANE="$TPM_sugarcane" \
    CLEAN_TPM_PURPLE="$TPM_purple" \
    CLEAN_META_SUGARCANE="$META_sugarcane" \
    CLEAN_OUT_PREFIX="${RESULTS}/figures/figure${FIG_REPRODUCTION}_reproduction" \
    CLEAN_CORES="$NUM_CORES" \
      "$RSCRIPT_PLOT" "${SCRIPTS}/20_fig_reproduction.r"
    ;;

  # Assembles every figure's generated legend into the paper-level file. Each
  # figure script writes its own <prefix>_legend.txt from the same variables that
  # drew the panels, so a number cannot disagree between a figure and its legend;
  # this only concatenates them in figure order.
  legends)
    LEG="${FIGURE_LEGENDS}"
    SRC=$(ls "${RESULTS}/figures/"figure*_legend.txt 2>/dev/null | sort -V) || true
    [ -n "$SRC" ] || die "no figure legends found -- run ./run.sh figdataset first"
    {
      echo "Figure legends"
      echo "=============="
      echo
      echo "Generated by new_clean/run.sh legends on $(date '+%F'). Do not edit by"
      echo "hand: each legend is written by the script that draws its figure, from"
      echo "the same variables, so the numbers cannot drift. Re-run the figure, then"
      echo "./run.sh legends."
      echo
      for f in $SRC; do echo; cat "$f"; echo; done
    } > "$LEG"
    echo "wrote $LEG  ($(echo "$SRC" | wc -l) figure legend(s))"
    ;;

  # --- 09-10 GO ---------------------------------------------------------------
  go)
    case "$ARG" in BP|MF|CC) ;; *) die "usage: run.sh go BP|MF|CC" ;; esac
    CLEAN_ONTOLOGY="$ARG" \
    CLEAN_RESULTS="$RESULTS" \
    CLEAN_EMAPPER_SUGARCANE="$EMAPPER_sugarcane" \
    CLEAN_EMAPPER_PURPLE="$EMAPPER_purple" \
    CLEAN_GO_P="$GO_P" CLEAN_GO_NTOP="$GO_NTOP" \
      conda run --no-capture-output -n "$CONDA_TOPGO" \
        Rscript "${SCRIPTS}/09_go_enrichment.r"
    ;;

  gosem)
    CLEAN_RESULTS="$RESULTS" CLEAN_ONTOLOGIES="$ONTOLOGIES" \
      conda run --no-capture-output -n "$CONDA_CLUSTERPROFILER" \
        Rscript "${SCRIPTS}/10_go_semantic.r"
    ;;

  # --- 11 H1 readouts ---------------------------------------------------------
  # None of these read an edge table -- only node_metrics, mcl membership, the
  # conserved gene lists and TPM matrices. The sequence-only steps (get_tfs 01-03,
  # myb61 01-06, module20 01-02b) are network-independent and are NOT re-run;
  # each config.sh reads them from CACHE_DIR in the original GET_TFS tree and
  # writes everything new under results/readouts/.
  tfs)
    check_study "$ARG"
    "${SCRIPTS}/11_readouts/get_tfs/04_postprocess.sh" "$ARG"
    ;;

  tfchar)
    "$RSCRIPT_PLOT" "${SCRIPTS}/11_readouts/get_tfs/05_tf_characterization.r"
    ;;

  myb61)
    "${SCRIPTS}/11_readouts/myb61/run_all.sh" "${ARG:-07}"
    ;;

  module20)
    "${SCRIPTS}/11_readouts/module20/run_all.sh" "${ARG:-03}"
    ;;

  # Sequence-evolution layer. Two of its steps are long and deliberate --
  # 02 (OrthoFinder, 3 species, 2h17m as measured) and 04 -- so it is usually
  # driven step by step from scripts/29_dnds/ rather than through this stage.
  dnds)
    "${SCRIPTS}/29_dnds/run_all.sh" "${ARG:-00}"
    ;;

  *) die "unknown stage '$STAGE'" ;;
  esac

  echo "=== ${STAGE} ${ARG} done :: $(date '+%F %T')"
}

main "$@"
