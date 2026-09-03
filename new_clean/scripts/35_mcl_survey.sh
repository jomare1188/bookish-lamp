#!/usr/bin/env bash
# =============================================================================
# 35_mcl_survey.sh -- choose the k-NN cut from the data, not by taste.
#
# WHY. MCL's giant module is a node-degree artefact. mclfaq(7) 7.3: "Preferably
# the network should not have nodes of very high degree... Such nodes tend to
# obscure cluster structure and contribute to coarse clusters." The documented
# remedy is a k-NN reduction, and clmprotocols(5) gives the target for
# co-expression graphs: "the median node degree should be at most one hundred
# neighbours", with the number of orphan nodes under ten percent.
#
# It also gives the rule for picking k: "a good heuristic is to choose a value
# that does not significantly change the number of singletons in the input
# graph." Both networks start with ZERO isolated nodes, so every singleton a
# reduction creates is a node it removed from the analysis entirely -- which is
# exactly the cost that has to be traded against the degree target.
#
# #knn INTERSECTS neighbour lists: an edge survives only if it is among the top k
# for BOTH endpoints. That is far stricter than it sounds on a graph with hubs --
# a hub's edges rarely make its partners' top-k -- which is why k values here are
# hundreds rather than the tens the protocol's BLAST example used.
#
# WRITES ONLY DIAGNOSTICS. No clustering is produced and no network is modified.
#
# RUN: through run.sh  ->  ./run.sh mclsurvey purple
# =============================================================================
set -euo pipefail

STUDY="${CLEAN_STUDY:?}"; WORK="${CLEAN_WORK_DIR:?}"; BIN_DIR="${CLEAN_MCL_BIN_DIR:?}"
OUT_DIR="${CLEAN_OUT_DIR:?}"; RANGE="${CLEAN_KNN_RANGE:-100/1000/100}"
NODE_METRICS="${CLEAN_NODE_METRICS:?}"; CORES="${CLEAN_CORES:-16}"
# The two thresholds the choice is made on, both from clmprotocols(5).
MAX_SINGLETON_PCT="${CLEAN_MAX_SINGLETON_PCT:-5}"
TARGET_MEDIAN_DEG="${CLEAN_TARGET_MEDIAN_DEG:-100}"

MCI="${WORK}/${STUDY}.mci"
say() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }
[ -s "$MCI" ] || { echo "FATAL: no matrix at $MCI -- run ./run.sh mclload $STUDY first" >&2; exit 1; }

# --- where the graph starts --------------------------------------------------
read -r N_NODES BASE_MED BASE_MEAN BASE_P90 BASE_MAX <<EOF
$(awk -F'\t' 'NR>1{d[NR]=$2; s+=$2; if($2>mx)mx=$2}
  END{n=NR-1; asort(d); printf "%d %d %.1f %d %d", n, d[int(n/2)], s/n, d[int(n*0.9)], mx}' "$NODE_METRICS")
EOF
say "$STUDY as built: $N_NODES nodes, median degree $BASE_MED, mean $BASE_MEAN, p90 $BASE_P90, max $BASE_MAX"
say "target: median degree <= $TARGET_MEDIAN_DEG with <= ${MAX_SINGLETON_PCT}% singletons (clmprotocols(5))"

# --- the sweep ---------------------------------------------------------------
# One matrix read covers every k in the range.
RAW="${OUT_DIR}/mcl_knn_survey_${STUDY}.raw"
OUT="${OUT_DIR}/mcl_knn_survey_${STUDY}.tsv"
if [ -s "$OUT" ] && [ "${CLEAN_FORCE:-0}" != "1" ]; then
  say "reusing $(basename "$OUT") -- CLEAN_FORCE=1 to re-measure"
else
say "mcx query -vary-knn $RANGE"
"${BIN_DIR}/mcx" query -imx "$MCI" -vary-knn "$RANGE" --output-table -t "$CORES" \
  > "$RAW" 2> "${OUT_DIR}/mcl_knn_survey_${STUDY}.stderr"

# --output-table gives raw node COUNTS in L/D/R/S and a fraction in E.
awk -F'\t' -v n="$N_NODES" -v OFS='\t' '
  NR==1 { print "k","edges_retained_frac","n_singletons","singleton_pct",
                "median_degree","mean_degree","pct_in_largest_component"; next }
  { k=$15; gsub(/[ \t]/,"",k)   # the table is space-padded; an unstripped k breaks every later match
    printf "%d\t%.6f\t%d\t%.3f\t%.1f\t%.2f\t%.2f\n", k, $5, $4, 100*$4/n, $11, $10, 100*$1/n }
' "$RAW" | sort -k1,1n > "$OUT"
say "wrote $(basename "$OUT")"
fi
column -t "$OUT"

# --- what to sweep -----------------------------------------------------------
# NOT one "best" k. The survey cannot know whether a reduction breaks the giant
# module -- only the sweep can -- so its job is to propose a defensible RANGE and
# let 36 test it. Three points are proposed, at roughly 1%, 3% and 5% singleton
# cost, and the sweep always also runs "none" so the current setting stays on the
# plot as a control.
# Median degree RISES with k, so the gentlest reduction that still meets the
# protocol's target is the LARGEST k satisfying it. Taking the smallest instead
# would shred a graph that needed only a trim: on purple, k=440 reaches median 96
# for 0.13% singletons, where the smallest-k rule proposed 40 and a 3.6% cost.
k_meeting_target() { awk -F'\t' -v t="$TARGET_MEDIAN_DEG" 'NR>1 && $5<=t {k=$1} END{if(k) print k}' "$OUT"; }
# ...and the aggressive end, for contrast: smallest k whose singleton cost is
# still acceptable.
k_at_singleton_cost() { awk -F'\t' -v lim="$1" 'NR>1 && $4<=lim {print $1; exit}' "$OUT"; }

K_TARGET=$(k_meeting_target)
K_TIGHT=$(k_at_singleton_cost 5)
if [ "$BASE_MED" -le "$TARGET_MEDIAN_DEG" ]; then
  # The degree criterion does not bind -- this network already meets it -- so it
  # cannot choose a k, and proposing the largest passing value would just be the
  # top of the grid. Fall back to the singleton-cost points.
  K_TARGET=""
  K_MID=$(k_at_singleton_cost 2)
else
  K_MID=""
fi
KSET=$(printf '%s\n' "$K_TARGET" "$K_TIGHT" "$K_MID" | awk 'NF && !seen[$0]++' | sort -n | tr '\n' ' ')
KSET=${KSET% }

say ""
if [ -z "$KSET" ]; then
  say "NO k in this range keeps singletons under 5%."
  say "  Widen MCL_SWEEP_KNN_RANGE upward, or accept a higher cost deliberately."
else
  say "PROPOSED k for the sweep: $KSET   (plus 'none', the current setting)"
  [ -n "$K_TARGET" ] && say "  k=$K_TARGET is the largest k still meeting the median-degree target"
  for k in $KSET; do
    awk -F'\t' -v k="$k" -v b="$BASE_MED" -v t="$TARGET_MEDIAN_DEG" '$1==k{
      printf "  k=%-5s singletons %5.2f%%   median degree %5.1f (was %d, target <= %d)   edges kept %.4f%s\n",
             k, $4, $5, b, t, $2, ($5+0 > t+0 ? "   [above target]" : "")
    }' "$OUT"
  done
fi

# Sugarcane's median degree is already 52, under the target -- its problem is the
# hub TAIL (p90 = 6,666), which the median cannot see. Purple's median is 859 and
# fails the target outright. So a k that looks unnecessary on the median may still
# be doing the work, and that is what the sweep measures.
awk -v b="$BASE_MED" -v t="$TARGET_MEDIAN_DEG" -v p="$BASE_P90" 'BEGIN{
  if (b+0 <= t+0)
    printf "\nNOTE: this network already meets the median-degree target (%d <= %d).\n      The case for reducing it rests on the hub tail instead: p90 = %d.\n", b, t, p
}'
say ""
say "set  MCL_SWEEP_K_${STUDY}=\"$KSET\"  (or pass as an argument), then ./run.sh mclsweep $STUDY"

# --- for the record only -----------------------------------------------------
# What correlation cut-off WOULD have produced a sane graph. This changes nothing
# here -- the networks are not being rebuilt -- but purple's |r| >= 0.8 is only
# p = 6.7e-5 at n = 18 where sugarcane's is p = 9e-12, so the two were never
# built at equal stringency, and this is the evidence a future rebuild needs.
say "recording the correlation-threshold curve (diagnostic only)"
"${BIN_DIR}/mcx" query -imx "$MCI" --vary-correlation --output-table -t "$CORES" \
  > "${OUT_DIR}/mcl_threshold_survey_${STUDY}.tsv" 2>/dev/null || \
  say "  (--vary-correlation failed; not fatal, it is diagnostic)"
say "done: $STUDY"
