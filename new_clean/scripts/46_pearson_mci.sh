#!/usr/bin/env bash
# =============================================================================
# 46_pearson_mci.sh -- the Pearson-only network as a native MCL matrix, direct.
#
# WHY THIS EXISTS RATHER THAN 40 + 34. The previous route wrote a 9-column edge
# table (7.7 GB sugarcane, 62 GB purple) whose only consumer was mcxload, which
# read columns 1, 2 and 6 and threw the rest away. That is 70 GB of disk and two
# full passes over it to move three columns. Here the Pearson layer is streamed
# through one awk straight into mcxload and the table is never written.
#
# The 9-column table is still the right artefact when the network itself is the
# product -- stats, figures, conservation all read it. It is the wrong artefact
# when the only question is how to cluster.
#
# WHAT "UNPRUNED" MEANS HERE: |r| >= 0.8, every surviving Pearson edge, no k-NN
# reduction. The MI layer is excluded -- it is 1.14% of sugarcane's edges and
# 4.20% of purple's, so the hairball is essentially all linear (see docs/results.md,
# "Choosing k-NN by a graph model").
#
# THE WEIGHT is the merge's own formula, reproduced exactly so these numbers stay
# comparable to every earlier run:
#     weight = 0.01 + (|r| - LO) / (HI - LO) * 0.99
# with LO = final_cut and HI = max_value read from the layer's summary json, not
# hardcoded. pearson_r keeps its sign in the network table; |r| is what MCL sees,
# because mcl requires positive entries and an anticorrelation is still a strong
# relationship.
#
# RUN: through run.sh  ->  RESULTS=$PWD/results_cluster \
#                          MCL_WORK_DIR=/dados04/jorge/tmp/mcl_work_cluster \
#                          ./run.sh pearsonmci sugarcane
# =============================================================================
set -euo pipefail

STUDY="${CLEAN_STUDY:?}"
LAYER="${CLEAN_LAYER:?}"          # <study>_pearson.edgelist.tsv
SUMMARY="${CLEAN_LAYER_SUMMARY:?}" # <study>_pearson.summary.json
WORK="${CLEAN_WORK_DIR:?}"
BIN="${CLEAN_MCL_BIN_DIR:?}"
CORES="${CLEAN_CORES:-16}"
FORCE="${CLEAN_FORCE:-0}"

MCI="${WORK}/${STUDY}.mci"
TAB="${WORK}/${STUDY}.tab"
say() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }

# --- guard: never write into the main tracks ---------------------------------
# purple.mci is 11 GB and takes half an hour to build. Two tracks sharing one
# filename would silently hand the wrong graph to whichever ran second.
case "$WORK" in
  */mcl_work|*/mcl_work/) echo "FATAL: $WORK is the MAIN work dir. Set MCL_WORK_DIR." >&2; exit 1;;
esac
mkdir -p "$WORK"

if [ -s "$MCI" ] && [ -s "$TAB" ] && [ "$FORCE" != "1" ]; then
  say "already built: $(basename "$MCI") ($(du -h "$MCI" | cut -f1)) -- CLEAN_FORCE=1 to rebuild"
  exit 0
fi

[ -s "$LAYER" ]   || { echo "FATAL: no Pearson layer at $LAYER" >&2; exit 1; }
[ -s "$SUMMARY" ] || { echo "FATAL: no layer summary at $SUMMARY" >&2; exit 1; }

# --- the header must be what we think it is ----------------------------------
# A reordered column would mean clustering on p-values, and every number
# downstream would be about a different graph without anything complaining.
HDR=$(head -1 "$LAYER")
C1=$(printf '%s' "$HDR" | cut -f1)
C2=$(printf '%s' "$HDR" | cut -f2)
C3=$(printf '%s' "$HDR" | cut -f3)
if [ "$C1" != "gene1" ] || [ "$C2" != "gene2" ] || [ "$C3" != "pearson" ]; then
  echo "FATAL: unexpected layer header." >&2
  echo "  expected columns 1,2,3 = gene1, gene2, pearson" >&2
  echo "  got                    = $C1, $C2, $C3" >&2
  exit 1
fi
say "header verified: 1=$C1  2=$C2  3=$C3"

# --- the weight scale, from the layer's own summary --------------------------
jnum() { sed -n "s/.*\"$1\"[[:space:]]*:[[:space:]]*\([0-9.eE+-]*\).*/\1/p" "$SUMMARY" | head -1; }
LO=$(jnum final_cut)
HI=$(jnum max_value)
N_EXP=$(jnum n_significant_edges)
if [ -z "$LO" ] || [ -z "$HI" ] || [ -z "$N_EXP" ]; then
  echo "FATAL: could not read final_cut / max_value / n_significant_edges from $SUMMARY" >&2
  exit 1
fi
awk -v lo="$LO" -v hi="$HI" 'BEGIN{ if (hi <= lo) { exit 1 } }' \
  || { echo "FATAL: max_value ($HI) <= final_cut ($LO)" >&2; exit 1; }
say "weight scale: 0.01 + (|r| - $LO)/($HI - $LO) * 0.99"
say "expecting $(printf "%'d" "$N_EXP") edges"

# --- stream, never materialise -----------------------------------------------
# --stream-mirror IS REQUIRED. The layer stores each undirected edge once, so
# without it mcxload builds a DIRECTED matrix and every gene appearing only in
# column 2 gets out-degree 0. The degree-0 assertion below is what catches that.
say "loading -> $(basename "$MCI")"
/usr/bin/time -v -o "${WORK}/${STUDY}.pearson.load.time" \
  bash -c "awk -F'\t' -v lo='$LO' -v hi='$HI' '
      NR>1 { r = \$3 < 0 ? -\$3 : \$3
             printf \"%s\t%s\t%.6g\n\", \$1, \$2, 0.01 + (r - lo)/(hi - lo) * 0.99 }
    ' '$LAYER' \
    | '${BIN}/mcxload' -abc - --stream-mirror --write-binary -o '$MCI' -write-tab '$TAB'"

PEAK=$(awk '/Maximum resident set size/{print $NF}' "${WORK}/${STUDY}.pearson.load.time")
WALL=$(awk -F': ' '/Elapsed \(wall clock\)/{print $NF}' "${WORK}/${STUDY}.pearson.load.time")
say "loaded: matrix $(du -h "$MCI" | cut -f1), tab $(wc -l < "$TAB") entries, peak RSS $(( PEAK / 1024 )) MB, wall $WALL"

# --- the matrix must BE the layer --------------------------------------------
DIMS=$("${BIN}/mcx" query -imx "$MCI" --dim 2>&1 | tail -1)
say "mcx query --dim: $DIMS"

# Arc count. mcx query lists per-node degrees; arcs = sum, edges = arcs/2.
ARCS=$("${BIN}/mcx" query -imx "$MCI" -t "$CORES" 2>/dev/null | awk 'NR>1{s+=$2} END{print s+0}')
EDGES=$(( ARCS / 2 ))
say "arcs $(printf "%'d" "$ARCS")  ->  edges $(printf "%'d" "$EDGES")  (expected $(printf "%'d" "$N_EXP"))"
if [ "$EDGES" != "$N_EXP" ]; then
  echo "FATAL: the matrix holds $EDGES edges, the layer reports $N_EXP." >&2
  echo "  A mismatch means mcxload dropped or merged entries -- most likely" >&2
  echo "  duplicate gene pairs in the layer, which mcxload silently collapses." >&2
  exit 1
fi

# Zero-degree nodes are impossible in a thresholded correlation network: every
# node is in the matrix because it had at least one surviving edge. Any here mean
# the load is wrong, and this is the check that catches a missing --stream-mirror.
NZERO=$("${BIN}/mcx" query -imx "$MCI" -t "$CORES" 2>/dev/null | awk 'NR>1 && $2==0 {n++} END{print n+0}')
say "nodes with degree 0: $NZERO"
if [ "$NZERO" != "0" ]; then
  echo "FATAL: $NZERO nodes have degree 0, so the loaded matrix is not the network." >&2
  exit 1
fi

# --- what the clusterers will be handed --------------------------------------
NODES=$(wc -l < "$TAB")
say ""
say "$STUDY unpruned Pearson graph:"
say "  nodes        $(printf "%'d" "$NODES")"
say "  edges        $(printf "%'d" "$EDGES")"
say "  mean degree  $(awk -v a="$ARCS" -v n="$NODES" 'BEGIN{printf "%.0f", a/n}')"
say ""
say "NOTE: mcl -scheme 7 keeps at most 1200 neighbours per node WHILE COMPUTING."
say "      Above that mean degree, the inflation ladder is partly measuring mcl's"
say "      pruner. That is what the -S probe in stage E is for."
say "done: $STUDY"
