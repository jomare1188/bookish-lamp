#!/usr/bin/env bash
# =============================================================================
# 34_mcl_load.sh -- load one network into MCL's native binary format, once.
#
# WHY. 05_mcl_clustering.r hands MCL a text .abc file it writes fresh on every
# run: 4.67 GB for sugarcane, 36.02 GB for purple, and for purple the write alone
# takes ~4 minutes on top of ~11 minutes of igraph work before MCL even starts.
# That cost is per-run, which is what makes a parameter sweep look unaffordable.
#
# Loading once into native binary format removes it. Everything downstream --
# the k-NN survey, the inflation sweep, clm info, clm dist -- reads this one
# matrix. The edge table is streamed straight into mcxload, so the text file is
# never materialised at all.
#
# WHAT IS AND IS NOT TOUCHED. This reads network_<study>_edges.tsv and writes to
# MCL_WORK_DIR. It writes nothing under results/ and modifies no input.
#
# THE WEIGHT COLUMN. Column 6 of the merged edge table is `weight`, the min-max
# normalisation of `stat` onto [0.01, 1.0] that 03_merge_layers.py computes. That
# is the value the clustering has always used, and using anything else here would
# silently change what is being clustered.
#
# RUN: through run.sh  ->  ./run.sh mclload sugarcane
# =============================================================================
set -euo pipefail

STUDY="${CLEAN_STUDY:?CLEAN_STUDY not set}"
EDGES="${CLEAN_EDGES:?CLEAN_EDGES not set}"
WORK="${CLEAN_WORK_DIR:?CLEAN_WORK_DIR not set}"
BIN_DIR="${CLEAN_MCL_BIN_DIR:?CLEAN_MCL_BIN_DIR not set}"
FORCE="${CLEAN_FORCE:-0}"

MCI="${WORK}/${STUDY}.mci"
TAB="${WORK}/${STUDY}.tab"
mkdir -p "$WORK"

say() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }

if [ -s "$MCI" ] && [ -s "$TAB" ] && [ "$FORCE" != "1" ]; then
  say "already loaded: $(basename "$MCI") ($(du -h "$MCI" | cut -f1)) -- FORCE=1 to rebuild"
  exit 0
fi

[ -s "$EDGES" ] || { echo "FATAL: no edge table at $EDGES" >&2; exit 1; }

# The header must be what we think it is. A silently reordered column would mean
# clustering on p-values instead of weights, and nothing downstream would notice.
HDR=$(head -1 "$EDGES")
C1=$(printf '%s' "$HDR" | cut -f1); C2=$(printf '%s' "$HDR" | cut -f2)
C6=$(printf '%s' "$HDR" | cut -f6)
if [ "$C1" != "gene1" ] || [ "$C2" != "gene2" ] || [ "$C6" != "weight" ]; then
  echo "FATAL: unexpected edge-table header." >&2
  echo "  expected columns 1,2,6 = gene1, gene2, weight" >&2
  echo "  got                    = $C1, $C2, $C6" >&2
  exit 1
fi
say "header verified: 1=$C1  2=$C2  6=$C6"

N_EDGES=$(( $(wc -l < "$EDGES") - 1 ))
say "loading $STUDY: $(printf "%'d" "$N_EDGES") edges -> $(basename "$MCI")"

# Stream, never materialise. mcxload reads '-' as stdin (verified).
#
# --stream-mirror IS REQUIRED AND IS NOT OPTIONAL. The edge table stores each
# undirected edge once, as `gene1 gene2 weight`. Without this flag mcxload builds
# a DIRECTED matrix from those arcs, so a gene that only ever appears in column 2
# ends up with out-degree 0. Measured on sugarcane without it: roughly half of
# the 103,336 nodes reported degree 0, and #knn -- which needs both endpoints'
# neighbour lists -- reduced the graph to ZERO edges at every k. The clustering
# would have been computed on a graph that is not this network.
/usr/bin/time -v -o "${WORK}/${STUDY}.load.time" \
  bash -c "awk -F'\t' 'NR>1{print \$1\"\t\"\$2\"\t\"\$6}' '$EDGES' \
    | '${BIN_DIR}/mcxload' -abc - --stream-mirror --write-binary -o '$MCI' -write-tab '$TAB'"

PEAK=$(awk '/Maximum resident set size/{print $NF}' "${WORK}/${STUDY}.load.time")
say "loaded. matrix $(du -h "$MCI" | cut -f1), tab $(wc -l < "$TAB") entries, peak RSS $(( PEAK / 1024 )) MB"

# --- the matrix must BE the network ------------------------------------------
# A mismatch here means every number downstream is about a different graph.
DIMS=$("${BIN_DIR}/mcx" query -imx "$MCI" --dim 2>&1 | tail -1)
say "mcx query --dim: $DIMS"

# Zero-degree nodes are impossible in a thresholded correlation network -- every
# row of node_metrics has degree >= 1 -- so any here mean the load is wrong.
# This is the check that catches a missing --stream-mirror.
NZERO=$("${BIN_DIR}/mcx" query -imx "$MCI" -t "${CLEAN_CORES:-8}" 2>/dev/null \
        | awk 'NR>1 && $2==0 {n++} END{print n+0}')
say "nodes with degree 0 in the loaded matrix: $NZERO"
if [ "$NZERO" != "0" ]; then
  echo "FATAL: $NZERO of the nodes have degree 0, so the loaded matrix is not" >&2
  echo "  this network. The usual cause is a missing --stream-mirror." >&2
  exit 1
fi
say "done: $STUDY"
