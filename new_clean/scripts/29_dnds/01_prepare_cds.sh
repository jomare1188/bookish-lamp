#!/usr/bin/env bash
# ============================================================================
# 01_prepare_cds.sh -- one CDS and one protein record per gene, per species,
# headered with the gene id the rest of the project uses.
#
# The gate that matters: translate(CDS) must reproduce the protein. pal2nal
# fails SILENTLY on a frame mismatch, so genes that do not translate cleanly
# are dropped here and counted in work/seqs/cds_qc_<species>.tsv.
# ============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

python3 ./01_prepare_cds.py

echo
echo "01 done. Sequence sets in ${SEQDIR}:"
for f in "$SC_CDS_CLEAN" "$SC_PEP_CLEAN" "$PU_CDS_CLEAN" "$PU_PEP_CLEAN" "$SB_CDS_CLEAN" "$SB_PEP_CLEAN"; do
  printf "  %-24s %8d records\n" "$(basename "$f")" "$(grep -c '^>' "$f")"
done
