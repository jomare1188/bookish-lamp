#!/usr/bin/env bash
# ============================================================================
# run_all.sh -- dN/dS layer, end to end
#
#   ./run_all.sh            every step
#   ./run_all.sh 03         from step 03 onwards
#   DNDS_SUBSET=200 ./run_all.sh 02b    a 200-orthogroup smoke test
#
# TWO STEPS ARE LONG AND ARE MEANT TO BE LAUNCHED DELIBERATELY:
#   02  OrthoFinder, three species, ~3-4 h
#   04  codeml over every alignment
# run_all.sh runs them like any other step; run them on their own if you would
# rather watch them.
#
# 05 (purple) is skipped by default: the purple direction is a 37 GB edge table
# and ~40 min. Run `./05_conserved_degree.sh both` when you want the purple
# replicate of test 1.
# ============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

FROM="${1:-00}"

run() {  # run <step> <script> [interpreter]
  local n="$1" script="$2" interp="${3:-}"
  if [[ "$n" < "$FROM" ]]; then echo "== skip $script"; return; fi
  echo
  echo "======================================================================"
  echo "== $script"
  echo "======================================================================"
  if [ -n "$interp" ]; then "$interp" "$script"; else "./$script"; fi
}

run 00  00_setup.sh
run 01  01_prepare_cds.sh
run 02  02_orthofinder_3sp.sh
run 02b 02b_build_triplets.sh
run 03  03_align.sh
run 04  04_codeml.sh
run 04b 04b_yn00_check.sh
run 05  05_conserved_degree.sh
run 06  06_join_and_test.r "$RSCRIPT_NET"
run 07  07_plots.r         "$RSCRIPT_PLOT"

echo
echo "======================================================================"
echo "Done. Key outputs in ${OUTDIR}:"
echo "  triplets.tsv                     1:1:1 orthologs, gated on the 2-species run"
echo "  dnds_raw.tsv                     codeml pairwise + exact NG86 counts"
echo "  yn00_crosscheck.tsv              independent-estimator agreement"
echo "  conserved_degree_<study>.tsv     per-gene edge-conservation fraction"
echo "  dnds_gene_table.tsv              everything joined, one row per orthogroup"
echo "  dnds_filter_report.tsv           what each filter removed"
echo "  dnds_tests.tsv                   the six tests, in reading order"
echo "  dnds_by_degree_decile.tsv        the stratified check"
echo "  dnds_saccharum_pair_binned.tsv   R570 vs LA purple, counts summed per bin"
echo "  fig_dnds.{png,pdf,svg}           the figure"
echo "======================================================================"
