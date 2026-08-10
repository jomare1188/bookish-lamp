#!/usr/bin/env bash
# ============================================================================
# run_all.sh — build the mutual-information networks for both studies
#
#   ./run_all.sh export            re-export both VST matrices from the dds
#   ./run_all.sh validate          run the correctness suite (fast, ~2 min)
#   ./run_all.sh sugarcane [args]  sweep the Munoz network
#   ./run_all.sh purple    [args]  sweep the Kiet network
#   ./run_all.sh all               export + validate + both sweeps
#
# The sweeps are long (~2.5 h for sugarcane, ~40 min for purple on an A4500),
# so they always pass --resume: re-running after a crash or a reboot picks up
# at the first unfinished tile. Extra arguments are forwarded to the engine,
# e.g.  ./run_all.sh sugarcane --estimator gcmi --out /tmp/quick_check
# ============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

log() { printf '\n=== %s ===\n' "$*"; }

export_one() {
  local label="$1" dds="$2" cols="$3" pearson="$4" prefix="$5"
  mkdir -p "$(dirname "$prefix")"
  MI_DDS="$dds" MI_PREFIX="$prefix" MI_COL_FILTER="$cols" \
  MI_PEARSON="$pearson" MI_LABEL="$label" MIN_CV="$MIN_CV" \
    "$RSCRIPT_DESEQ" 00_export_vst.r
}

sweep_one() {
  local label="$1" prefix="$2" outdir="$3"; shift 3
  mkdir -p "$outdir"
  "$PYTORCH" custom_mutual_information.py \
      --matrix "$prefix" \
      --out    "${outdir}/${label}_${ESTIMATOR}" \
      --estimator "$ESTIMATOR" --k "$KSG_K" \
      --alpha "$ALPHA" --gpu-mem-gb "$GPU_MEM_GB" \
      --match-pearson "$MATCH_PEARSON" --cand-p "$CAND_P" \
      --max-perm "$MAX_PERM" --ram-limit-gb "$RAM_LIMIT_GB" --resume "$@"
}

# No default action: `all` launches multi-hour sweeps, so it has to be asked for.
case "${1:-usage}" in
  export)
    log "exporting sugarcane"
    export_one sugarcane "$DDS_SUGARCANE" "$COLS_SUGARCANE" \
               "$PEARSON_SUGARCANE" "${OUT_SUGARCANE}/sugarcane"
    log "exporting purple"
    export_one purple "$DDS_PURPLE" "$COLS_PURPLE" \
               "$PEARSON_PURPLE" "${OUT_PURPLE}/purple"
    ;;
  validate)
    log "validation suite"
    "$PYTORCH" 99_validate.py
    ;;
  sugarcane)
    shift || true
    log "sugarcane sweep (48 libraries)"
    sweep_one sugarcane "${OUT_SUGARCANE}/sugarcane" "$OUT_SUGARCANE" "$@"
    ;;
  purple)
    shift || true
    log "purple sweep (18 libraries)"
    sweep_one purple "${OUT_PURPLE}/purple" "$OUT_PURPLE" "$@"
    ;;
  all)
    "$0" export
    "$0" validate
    "$0" sugarcane
    "$0" purple
    ;;
  *)
    sed -n '3,13p' "$0" | sed 's/^# \{0,1\}//'; exit 1
    ;;
esac
