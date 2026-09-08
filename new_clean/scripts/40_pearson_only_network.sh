#!/usr/bin/env bash
# =============================================================================
# 40_pearson_only_network.sh -- a Pearson-only network, from the layer we have.
#
# WHY. The MI layer is a minority of the hairball -- 4.20% of purple's edges and
# 1.14% of sugarcane's -- and the k-NN work is about reducing degree, not about
# testing non-linearity. Dropping it makes the graph smaller and the story
# simpler, and costs almost nothing.
#
# NO REBUILD IS NEEDED. 02_network_engine.py already wrote the Pearson layer as
# its own file; 03_merge_layers.py's only job for a pearson-source edge is to add
# four columns that are pure per-row functions of what is already there. So this
# is one streaming awk, not 3 hours of GPU:
#
#   stat      = |r|                                     (the merged `stat` for a
#                                                        pearson-source edge)
#   weight    = 0.01 + (|r| - STAT_MIN)/(STAT_MAX - STAT_MIN) * 0.99
#   pearson_r = r, SIGNED -- purple has negative correlations and the merged
#               table keeps the sign here while `stat` takes the absolute value
#   ksg       = NA
#   source    = pearson
#
# THIS MUST NOT LAND IN THE MAIN TREE. network_tsv() has no estimator tag, so a
# Pearson-only network written to the default RESULTS would silently replace the
# merged network. Always run with RESULTS pointed elsewhere; the caller is
# responsible and the guard below checks it.
#
# RUN: through run.sh  ->  RESULTS=$PWD/results_pearson ./run.sh pearsononly purple
# =============================================================================
set -euo pipefail

STUDY="${CLEAN_STUDY:?}"; LAYER="${CLEAN_PEARSON_LAYER:?}"; OUT="${CLEAN_OUT:?}"
LO="${CLEAN_STAT_MIN:-0.8}"; HI="${CLEAN_STAT_MAX:-0.9999}"
FORCE="${CLEAN_FORCE:-0}"
say() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }

[ -s "$LAYER" ] || { echo "FATAL: no Pearson layer at $LAYER" >&2; exit 1; }

# --- refuse to overwrite a merged network ------------------------------------
# The merged summary records n_mi_only; a Pearson-only one records 0. If the
# target already exists and was NOT Pearson-only, this is the main tree and we
# are about to destroy it.
SUM="$(dirname "$OUT")/network_${STUDY}_edges.summary.json"
if [ -s "$SUM" ] && [ "$FORCE" != "1" ]; then
  MI=$(sed -n 's/.*"n_mi_only": *\([0-9]*\).*/\1/p' "$SUM" | head -1)
  if [ -n "$MI" ] && [ "$MI" != "0" ]; then
    echo "REFUSING: $SUM records n_mi_only = $MI, so $(dirname "$OUT") holds a MERGED" >&2
    echo "  network, not a Pearson-only one. Point RESULTS at a separate tree." >&2
    exit 1
  fi
fi
if [ -s "$OUT" ] && [ "$FORCE" != "1" ]; then
  say "already built: $(basename "$OUT") -- FORCE=1 to rebuild"; exit 0
fi

# The layer's own summary carries the exact edge count, so there is no need to
# read the 44 GB file just to count it.
LSUM="${LAYER%.edgelist.tsv}.summary.json"
N_EXPECT=$(sed -n 's/.*"n_significant_edges": *\([0-9]*\).*/\1/p' "$LSUM" | head -1)
say "$STUDY: ${N_EXPECT:-?} Pearson edges -> $(basename "$OUT")"
say "  weight = 0.01 + (|r| - $LO)/($HI - $LO) * 0.99"

HDR=$(head -1 "$LAYER")
[ "$(printf '%s' "$HDR" | cut -f1-3)" = "$(printf 'gene1\tgene2\tpearson')" ] || {
  echo "FATAL: unexpected layer header: $HDR" >&2; exit 1; }

mkdir -p "$(dirname "$OUT")"
awk -F'\t' -v OFS='\t' -v lo="$LO" -v hi="$HI" '
  NR==1 { print "gene1","gene2","stat","pval","padj","weight","pearson_r","ksg","source"; next }
  { r=$3+0; a=(r<0? -r : r);
    printf "%s\t%s\t%.5f\t%s\t%s\t%.6f\t%.5f\tNA\tpearson\n",
           $1, $2, a, $4, $5, 0.01 + (a-lo)/(hi-lo)*0.99, r }
' "$LAYER" > "${OUT}.part"

N_GOT=$(( $(wc -l < "${OUT}.part") - 1 ))
if [ -n "$N_EXPECT" ] && [ "$N_GOT" != "$N_EXPECT" ]; then
  echo "FATAL: wrote $N_GOT edges, the layer summary says $N_EXPECT" >&2
  rm -f "${OUT}.part"; exit 1
fi
mv "${OUT}.part" "$OUT"
say "  wrote $N_GOT edges, verified against the layer summary"

# A companion summary so 22_fig_topology.r's n_pearson_only/n_both/n_mi_only
# parse still works, and so the guard above can recognise this tree next time.
NS=$(sed -n 's/.*"n_samples": *\([0-9]*\).*/\1/p' "$LSUM" | head -1)
cat > "$SUM" <<EOF
{
  "study": "$STUDY",
  "n_samples": ${NS:-null},
  "stat_min": $LO,
  "stat_max": $HI,
  "n_pearson_only": $N_GOT,
  "n_both": 0,
  "n_mi_only": 0,
  "n_total": $N_GOT,
  "layers": "pearson only -- the MI layer is deliberately excluded",
  "source_layer": "$LAYER",
  "out_file": "$OUT"
}
EOF
say "  wrote $(basename "$SUM")"
say "done: $STUDY"
