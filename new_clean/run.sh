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
#   ./run.sh conserve  sugarcane_to_purple|purple_to_sugarcane
#   ./run.sh conservenull <direction>        permutation null for the above
#   ./run.sh trait     <study>               gene-trait correlations
#   ./run.sh traitmi   <study>               gene-vs-trait MI (non-linear)
#   ./run.sh conscor   [pearson|mi|union]    conserved N response (node + edge)
#   ./run.sh go        BP|MF|CC              GO enrichment
#   ./run.sh gosem                           GO semantic clustering
#   ./run.sh tfs       <study>               TFs in the network (step 04 only)
#   ./run.sh tfchar                          TF characterisation
#   ./run.sh myb61     [from-step]           MYB61 readout (default 07)
#   ./run.sh module20  [from-step]           Module-20 readout (default 03)
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

usage() { sed -n '3,34p' "${BASH_SOURCE[0]}" | sed 's/^# \?//'; exit 1; }
die()   { echo "ERROR: $*" >&2; exit 1; }

# Per-study value lookup, e.g. cfg DDS sugarcane -> $DDS_sugarcane
cfg() { local v="$1_$2"; echo "${!v-}"; }

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

  # --- 05 clustering ----------------------------------------------------------
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

  conscor)
    CLEAN_RESULTS="$RESULTS" \
    CLEAN_OUT_DIR="${RESULTS}/conservation" \
    CLEAN_ORTHOGROUPS="$ORTHOGROUPS" \
    CLEAN_SELECT_TRAIT="$SELECT_TRAIT" \
    CLEAN_SELECTION="${ARG:-$TRAIT_SELECTION}" \
    CLEAN_TRAIT_R_THR="$TRAIT_R_THR" \
    CLEAN_TRAIT_PADJ_THR="$TRAIT_PADJ_THR" \
    CLEAN_CORES="$NUM_CORES" \
      "$RSCRIPT_NET" "${SCRIPTS}/08_conserved_cor_genes.r"
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

  *) die "unknown stage '$STAGE'" ;;
  esac

  echo "=== ${STAGE} ${ARG} done :: $(date '+%F %T')"
}

main "$@"
