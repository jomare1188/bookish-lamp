#!/usr/bin/env bash
# ============================================================================
# run_all.sh — Muñoz Module-20 cross-network test, end to end
#
#   ./run_all.sh          run every step
#   ./run_all.sh 03       run from step 03 onwards
#
# All local (no internet); reuses the DIAMOND databases built by scripts/myb61.
# Step 01 streams a 2.7 GB FASTA once; the rest is minutes.
# ============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

FROM="${1:-01}"

run() {
  local n="$1" script="$2" interp="${3:-}"
  if [ "$n" \< "$FROM" ]; then echo "== skip $script"; return; fi
  echo; echo "======================================================================"
  echo "== $script"
  echo "======================================================================"
  if [ "$interp" = "conda-topgo" ]; then
    conda run --no-capture-output -n topGO_env Rscript "$script"
  elif [ -n "$interp" ]; then "$interp" "$script"; else "./$script"; fi
}

run 01 01_extract_module20.sh
run 02 02_map_to_references.sh
run 02b 02b_tf_call_on_munoz.sh
run 03 03_network_readout.r    "$RSCRIPT_PLOT"
run 04 04_nitrogen_response.r  "$RSCRIPT_PLOT"
run 05 05_module20_heatmaps.r  "$RSCRIPT_PLOT"

# 06-08: the AtMYB59 copies as NETWORK objects -- degree, cross-species
# neighbourhood overlap, and what the neighbourhoods are enriched for. 06 streams
# both .pairs dumps (~2 min) and needs numpy/pandas; 07 needs topGO_env.
run 06 06_myb59_neighbourhoods.py "$PYTORCH"
run 07 07_myb59_neighbour_go.r    "conda-topgo"
run 08 08_myb59_readout_figure.r  "$RSCRIPT_PLOT"

echo
echo "======================================================================"
echo "Done. Key outputs in ${OUTDIR}:"
echo "  module20_members.tsv                 the 12 published members"
echo "  module20_orfs.tsv / module20.faa      member ORFs from the Munoz proteome"
echo "  module20_tf_call_comparison.tsv      our TF call on THEIR proteins vs theirs"
echo "  module20_map_<sp>.tsv                every hit + reciprocal verdict"
echo "  module20_mapping_summary.tsv         genes mapped per member"
echo "  module20_network_table.tsv           per-gene network readout"
echo "  module20_vs_background.tsv           vs all nodes / TFs / MYBs"
echo "  module20_tests.tsv                   centrality, conservation, modules"
echo "  module20_module_membership.tsv       which of our modules they fall in"
echo "  module20_nitrogen_<sp>.tsv           design-aware N tests"
echo "  module20_network_overview.{png,pdf}  network figure"
echo "  module20_heatmap_<sp>.{png,pdf}      full ComplexHeatmap figure per species"
echo "  myb59_neighbour_overlap.tsv          cross-species neighbourhood overlap + null"
echo "  myb59_within_sugarcane_overlap.tsv   copy-vs-copy overlap (they are not independent)"
echo "  myb59_neighbour_go.tsv               BP/MF/CC enrichment of each neighbourhood"
echo "  myb59_copies_readout.{png,pdf}       the four-panel readout"
echo "======================================================================"
