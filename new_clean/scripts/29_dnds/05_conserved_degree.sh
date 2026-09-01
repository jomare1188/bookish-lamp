#!/usr/bin/env bash
# ============================================================================
# 05_conserved_degree.sh -- per-gene edge-conservation fraction
#
# The predictor the hypothesis is actually about does not exist yet as a
# per-gene number. results/conservation/ has the binary
# conserved_genes_<study>_FULL.txt (>=1 conserved edge, yes/no) and the full
# per-edge table with a conserved flag, but nothing in between -- and "what
# FRACTION of this gene's neighbourhood survives the species jump" is a far
# better graded predictor than a yes/no that 38% of the network satisfies.
#
# One streaming awk pass. Sugarcane is 5 GB / 76 M edges (~5 min); purple is
# 39 GB / 705 M edges (~45 min), so it is opt-in.
#
#   ./05_conserved_degree.sh              sugarcane only (default)
#   ./05_conserved_degree.sh purple       purple only
#   ./05_conserved_degree.sh both
# ============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

WHICH="${1:-sugarcane}"
[ "$WHICH" = "both" ] && STUDIES="sugarcane purple" || STUDIES="$WHICH"

for study in $STUDIES; do
  eval "edges=\$CONS_EDGES_${study}"
  eval "metrics=\$NODE_METRICS_${study}"
  out="${OUTDIR}/conserved_degree_${study}.tsv"
  [ -s "$edges" ] || { echo "FATAL: missing $edges"; exit 1; }

  echo "== ${study}: $(du -h "$edges" | cut -f1) -> $(basename "$out")"
  awk -F'\t' -v OFS='\t' '
    NR == 1 {
      for (i = 1; i <= NF; i++) h[$i] = i
      g1 = h["gene1"]; g2 = h["gene2"]; cv = h["conserved"]
      if (!g1 || !g2 || !cv) { print "FATAL: unexpected header" > "/dev/stderr"; exit 1 }
      next
    }
    {
      n[$g1]++; n[$g2]++
      if ($cv == "TRUE") { c[$g1]++; c[$g2]++ }
    }
    END {
      print "gene", "n_edges", "n_conserved", "frac_conserved"
      for (g in n) print g, n[g], (g in c ? c[g] : 0), (g in c ? c[g] / n[g] : 0)
    }' "$edges" > "$out"

  echo "   genes: $(( $(wc -l < "$out") - 1 ))"

  # cross-check against the degree the network stage recorded: a mismatch means
  # the edge table and the node metrics came from different builds, and every
  # join downstream would be quietly wrong.
  if [ -s "$metrics" ]; then
    bad=$(awk -F'\t' 'NR==FNR{if(FNR>1) d[$1]=$2; next}
                      FNR>1 && ($1 in d) && d[$1] != $2 {n++}
                      END{print n+0}' "$metrics" "$out")
    if [ "$bad" -eq 0 ]; then
      echo "   degree cross-check vs node_metrics: OK"
    else
      echo "   WARNING: degree differs from node_metrics for ${bad} genes"
    fi
  fi
done
