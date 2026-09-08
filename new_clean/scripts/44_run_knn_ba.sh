#!/usr/bin/env bash
# =============================================================================
# run_knn_ba.sh -- the whole Pearson-only k-NN/BA experiment, start to finish.
#
# Meant to be launched in a screen session and left alone:
#
#     screen -S knnba
#     cd /dados04/jorge/comparative_saccharum/new_clean
#     ./scripts/44_run_knn_ba.sh 2>&1 | tee logs/knn_ba_full.log
#
# IDEMPOTENT. Every stage skips work that is already on disk, so a detach, a
# crash or a Ctrl-C costs only the stage that was in flight. Re-running is safe.
#
# MEASURED COSTS, so the log is readable against an expectation:
#   pearsononly   sugarcane ~3 min (done already), purple ~30 min for its 44 GB layer
#   mclload       sugarcane ~2.5 min (done already), purple ~35 min
#   knnsweep      one #knn reduction is ~2 min on sugarcane, ~15 min on purple;
#                 KNN_BA_JOBS of them run at once
#   knnselect     ~939 s per k on sugarcane at k=200 (96,335 nodes, 1.7 M edges);
#                 purple is larger. All k values run concurrently, so a study is
#                 one batch rather than twelve serial jobs.
# Sugarcane's Pearson-only network and MCL matrix are ALREADY BUILT, so this
# picks up at its k sweep and does purple from scratch.
#
# NOTHING TOUCHES THE MAIN ANALYSIS. RESULTS and MCL_WORK_DIR are pinned to the
# parallel tree below; 40_pearson_only_network.sh additionally refuses to write
# into a tree that holds a merged network, and run.sh's threshold guard refuses
# cross-tree writes.
# =============================================================================
set -uo pipefail
cd /dados04/jorge/comparative_saccharum/new_clean

export RESULTS="$PWD/results_pearson"
export MCL_WORK_DIR=/dados04/jorge/tmp/mcl_work_pearson
STUDIES="${STUDIES:-sugarcane purple}"
SELECT_JOBS="${KNN_BA_SELECT_JOBS:-12}"

say() { printf '\n=== [%s] %s\n' "$(date +%H:%M:%S)" "$*"; }
mkdir -p "$MCL_WORK_DIR"

# --- 0. is the criterion allowed to choose? ----------------------------------
if [ ! -s "$RESULTS/statgraph_validation.tsv" ]; then
  say "validating statGraph's criterion before it chooses anything"
  ./run.sh sgvalidate || say "validation returned non-zero -- read the table before trusting gic_ba"
fi

for S in $STUDIES; do
  say "$S: Pearson-only network"
  ./run.sh pearsononly "$S" || { say "$S FAILED at pearsononly"; continue; }

  say "$S: topology"
  ./run.sh stats "$S" || say "$S stats failed (not fatal for the sweep)"

  say "$S: load into MCL's native format"
  ./run.sh mclload "$S" || { say "$S FAILED at mclload"; continue; }

  say "$S: reduce at every k, and measure"
  ./run.sh knnsweep "$S" || { say "$S FAILED at knnsweep"; continue; }

  say "$S: score each k against ER / WS / BA"
  # statGraph is single-threaded on purpose (its own numCores path deadlocks
  # across repeated calls), so the k values are parallelised here instead.
  KS=$(awk -F'\t' 'NR>1{print $2}' "$RESULTS/$S/knnba_grid_${S}.tsv")
  printf '%s\n' $KS | xargs -P "$SELECT_JOBS" -I{} \
    env RESULTS="$RESULTS" ./run.sh knnselect "$S" {} \
    > "$RESULTS/$S/knnselect_${S}.log" 2>&1
  say "$S: per-k scoring done -- see knnselect_${S}.log"

  say "$S: the winning k"
  ./run.sh knncollect "$S"
done

say "sweep complete. Cluster at the winner with:"
for S in $STUDIES; do
  W=$(awk -F'\t' 'NR>1 && $10!="NA" {if(b==""||$10<b){b=$10;k=$2}} END{print k}' \
       "$RESULTS/$S/knnba_selection_${S}.tsv" 2>/dev/null)
  [ -n "${W:-}" ] && echo "    MCL_KNN_${S}=${W}  RESULTS=\$PWD/results_pearson ./run.sh mcl $S"
done
