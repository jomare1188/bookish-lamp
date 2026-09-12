#!/usr/bin/env bash
set -uo pipefail
cd /dados04/jorge/comparative_saccharum/new_clean
P=/dados04/jorge/comparative_saccharum/files/fix_orthofinder/proteins
A=/dados04/jorge/comparative_saccharum/annotation
for pair in "sugarcane:$P/sugarcane_one_transcript.fa:Saccharum officinarum" \
            "purple:$P/one_transcript_purple_proteins.faa:Saccharum officinarum"; do
  s="${pair%%:*}"; rest="${pair#*:}"; f="${rest%%:*}"; sp="${rest#*:}"
  echo "=== $(date '+%F %T') $s ==="
  CLEAN_STUDY="$s" CLEAN_PROTEOME="$f" \
  CLEAN_OUT_DIR="$A/$s/pannzer" \
  CLEAN_PZ_DIR=/dados04/jorge/tmp/pannzer_install/SANSPANZ.3 \
  CLEAN_PZ_PYTHON=/dados04/jorge/tmp/pannzer_install/pz_env/bin/python \
  CLEAN_PZ_SPECIES="$sp" CLEAN_PZ_CHUNK=1000 CLEAN_PZ_PARALLEL=2 \
    bash scripts/59_run_pannzer.sh
  echo "=== $(date '+%F %T') $s exit $? ==="
done
