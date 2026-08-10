#!/usr/bin/env bash
# ============================================================================
# run_all.sh — MYB61 orthologue search, end to end
#
#   ./run_all.sh            run every step
#   ./run_all.sh 03         run from step 03 onwards
#
# Steps 01-02 need internet (Ensembl Plants). Everything else is local.
# Total runtime is a few minutes; the DIAMOND databases are cached in
# results/myb61/work/ and reused.
# ============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

FROM="${1:-01}"

run() {  # run <step-number> <script> [interpreter]
  local n="$1" script="$2" interp="${3:-}"
  if [ "$n" \< "$FROM" ]; then echo "== skip $script"; return; fi
  echo
  echo "======================================================================"
  echo "== $script"
  echo "======================================================================"
  if [ -n "$interp" ]; then "$interp" "$script"; else "./$script"; fi
}

run 01 01_fetch_references.sh
run 02 02_build_query.sh
run 03 03_search_proteomes.sh
run 04 04_rbh_filter.sh
run 05 05_orthogroup_check.sh
run 06 06_phylogeny.sh
run 07 07_network_readout.r "$RSCRIPT_PLOT"
run 08 08_myb61_expression_test.r "$RSCRIPT_PLOT"

echo
echo "======================================================================"
echo "Done. Key outputs in ${OUTDIR}:"
echo "  MYB61_query.faa                       At/Sb/Os anchor set"
echo "  <sp>_MYB61_candidates.tsv             every hit + both verdicts"
echo "  MYB61_orthologs_<sp>.ids              confirmed orthologue genes"
echo "  MYB61_orthogroups.tsv                 OrthoFinder bridge view"
echo "  MYB61_bridge_status.tsv               per-copy bridge status"
echo "  phylogeny/MYB61.tree                  clade confirmation"
echo "  MYB61_network_table.tsv               per-copy network readout"
echo "  MYB61_conserved_edge_enrichment.tsv   enrichment test"
echo "  MYB61_network_overview.{png,pdf}      figure"
echo "  MYB61_expression_anova.tsv            per-genotype ANOVA / U-shape test"
echo "  MYB61_expression_heatmap.{png,pdf}    expression across the 18 leaf samples"
echo "======================================================================"
