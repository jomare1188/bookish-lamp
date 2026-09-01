#!/usr/bin/env bash
# ============================================================================
# 03_align.sh -- protein alignment (MAFFT L-INS-i) -> codon alignment (pal2nal)
#
# NOT reusing OrthoFinder's MultipleSequenceAlignments/: FAMSA output that has
# been through OrthoFinder's column trimming no longer maps 1:1 onto the CDS,
# which is exactly the correspondence pal2nal requires. Three short sequences
# realign in milliseconds, so there is nothing to save by reusing them.
# ============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh
PAL2NAL="$(cat "${WORKDIR}/pal2nal.path")"
export PAL2NAL MAFFT MIN_CODONS

[ -s "$TRIPLETS" ] || { echo "FATAL: no ${TRIPLETS} -- run 02b first"; exit 1; }

python3 ./03_explode.py

LOG="${WORKDIR}/align_status.tsv"
echo "== aligning with ${DNDS_THREADS} workers"
if command -v parallel >/dev/null; then
  parallel --colsep '\t' -j "$DNDS_THREADS" ./_align_one.sh {1} {2} \
    :::: "${WORKDIR}/og_list.txt" > "$LOG"
else
  xargs -a "${WORKDIR}/og_list.txt" -P "$DNDS_THREADS" -L1 ./_align_one.sh > "$LOG"
fi

echo "== codon alignments"
awk -F'\t' '{split($2,a,"_"); print a[1]}' "$LOG" | sort | uniq -c | sort -rn | sed 's/^/   /'
echo "   ok: $(grep -c $'\tok_' "$LOG") -> ${WORKDIR}/og/"
