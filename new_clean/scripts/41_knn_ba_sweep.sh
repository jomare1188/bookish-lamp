#!/usr/bin/env bash
# =============================================================================
# 41_knn_ba_sweep.sh -- reduce the network at each k, and measure it exactly.
#
# Produces, per k: the reduced matrix, an exact node/edge/degree count, and a
# labelled edge list for the model-selection step. Nothing here decides anything;
# 42 does the choosing.
#
# WHY EXACT COUNTS AND NOT THE SURVEY'S. `mcx query -vary-knn` gives the whole
# grid in one read, but its L/D/R/S columns are node counts and its E column is a
# FRACTION of the ARC count -- so edges have to be back-computed and rounded. The
# reduced matrix gives the real numbers directly, and we need the matrix anyway
# to dump an edge list. The survey is still run first as a cheap cross-check:
# `nodes with >= 1 edge` from both routes must agree at every k.
#
# #knn INTERSECTS neighbour lists -- an edge survives only if it is among the top
# k for BOTH endpoints -- so a small k can empty the graph. #knnj joins instead
# and is the fallback; which one was used is recorded per k.
#
# PARALLELISM. The reductions are independent and the source matrix stays in page
# cache after the first read, so they run concurrently. statGraph itself is left
# SINGLE-THREADED in 42: its numCores path creates and destroys a PSOCK cluster
# per spectral density, which leaks connections and deadlocks at 32 workers
# (measured). Parallelism belongs at this level, across k, not inside it.
#
# RUN: through run.sh  ->  RESULTS=... ./run.sh knnsweep purple
# =============================================================================
set -euo pipefail

STUDY="${CLEAN_STUDY:?}"; WORK="${CLEAN_WORK_DIR:?}"; BIN="${CLEAN_MCL_BIN_DIR:?}"
OUT_DIR="${CLEAN_OUT_DIR:?}"; CORES="${CLEAN_CORES:-16}"
K_LIST="${CLEAN_K_LIST:-50 100 150 200 250 300 350 400 450 500 550 600}"
JOBS="${CLEAN_JOBS:-6}"

MCI="${WORK}/${STUDY}.mci"; TAB="${WORK}/${STUDY}.tab"
SW="${WORK}/knnba_${STUDY}"; mkdir -p "$SW" "$OUT_DIR"
say() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }
[ -s "$MCI" ] || { echo "FATAL: no matrix at $MCI -- run ./run.sh mclload $STUDY" >&2; exit 1; }

# node/edge/degree straight from a matrix: degree column of mcx query
stats_of() {
  "$BIN/mcx" query -imx "$1" 2>/dev/null | awk '
    NR>1 { tot++; if ($2>0) { n++; s+=$2; d[n]=$2 } }
    END { if (!n) { print 0, 0, 0, 0, tot; exit }
          asort(d); printf "%d %d %.2f %d %d\n", n, s/2, s/n, d[int((n+1)/2)], tot-n }'
}

N_TOTAL=$("$BIN/mcx" query -imx "$MCI" 2>/dev/null | awk 'END{print NR-1}')
say "$STUDY: $N_TOTAL nodes in the matrix domain"
say "k values: $K_LIST   (${JOBS} reductions at a time)"

# --- the cheap whole-grid cross-check ----------------------------------------
KMIN=$(echo $K_LIST | tr ' ' '\n' | sort -n | head -1)
KMAX=$(echo $K_LIST | tr ' ' '\n' | sort -n | tail -1)
KSTEP=$(echo $K_LIST | awk '{print ($2>0 && $1>0)? $2-$1 : 50}')
say "cross-check survey: -vary-knn ${KMIN}/${KMAX}/${KSTEP}"
"$BIN/mcx" query -imx "$MCI" -vary-knn "${KMIN}/${KMAX}/${KSTEP}" --output-table -t "$CORES" \
  > "${OUT_DIR}/knnba_survey_${STUDY}.tsv" 2>/dev/null || say "  (survey failed; not fatal)"

# --- one reduction per k, in parallel ----------------------------------------
reduce_one() {
  local K="$1"
  local MAT="${SW}/${STUDY}.knn${K}.mci"
  local MODE=knn
  if [ ! -s "$MAT" ]; then
    "$BIN/mcx" alter -imx "$MCI" -tf "#knn($K)" -o "$MAT" 2>/dev/null
    local A
    A=$("$BIN/mcx" query -imx "$MAT" 2>/dev/null | awk 'NR>1{s+=$2}END{print s+0}')
    if [ "${A:-0}" -lt 1000 ]; then
      "$BIN/mcx" alter -imx "$MCI" -tf "#knnj($K)" -o "$MAT" 2>/dev/null
      MODE=knnj
    fi
  fi
  echo "$MODE" > "${MAT}.mode"
  # labelled edge list for statGraph; --dump-upper gives each edge once
  local EL="${SW}/${STUDY}.knn${K}.edges"
  [ -s "$EL" ] || "$BIN/mcxdump" -imx "$MAT" -tab "$TAB" \
      --dump-pairs --dump-upper --no-loops --no-values -o "$EL" 2>/dev/null
  printf '%s\n' "$K done"
}
export -f reduce_one; export SW STUDY BIN MCI TAB

printf '%s\n' $K_LIST | xargs -P "$JOBS" -I{} bash -c 'reduce_one {}' >/dev/null
say "all reductions complete"

# --- the table ---------------------------------------------------------------
TSV="${OUT_DIR}/knnba_grid_${STUDY}.tsv"
printf 'study\tk\tmode\tnodes\tedges\tmean_degree\tmedian_degree\tsingletons\tedge_list\n' > "$TSV"
for K in $K_LIST; do
  MAT="${SW}/${STUDY}.knn${K}.mci"
  [ -s "$MAT" ] || { say "  k=$K MISSING"; continue; }
  read -r NODES EDGES MEAN MED SGL <<<"$(stats_of "$MAT")"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$STUDY" "$K" "$(cat "${MAT}.mode")" "$NODES" "$EDGES" "$MEAN" "$MED" "$SGL" \
    "${SW}/${STUDY}.knn${K}.edges" >> "$TSV"
done
say "wrote $(basename "$TSV")"
column -t "$TSV"
say "done: $STUDY"
