#!/usr/bin/env bash
# ============================================================================
# 14_centrality.sh -- local clustering coefficient for the sugarcane network
#
# Reuses the graph-tool object the SBM stage already wrote, so the 8 GB edge
# table is never re-parsed. Runs under the sbm_3.7 env's python (graph-tool),
# not the pipeline's R interpreters.
#
#   COMPUTE_CORENESS=1 ./14_centrality.sh    also emit k-core (~5 min extra)
# ============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh
[ -s "$GRAPH_sugarcane" ] || { echo "FATAL: missing $GRAPH_sugarcane"; exit 1; }
"$GT_PYTHON" ./14_centrality.py
