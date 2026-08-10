#!/bin/bash
# ============================================================================
# 01_split_proteome.sh — split the proteome FASTA into $N_CHUNKS pieces.
# Dependency-free (awk only); replaces the old split_fasta.py / Bio.SeqIO.
# Usage:  ./01_split_proteome.sh <species>
# ============================================================================
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${here}/config.sh" "${1:?usage: 01_split_proteome.sh <species>}"

mkdir -p "$CHUNKDIR"
rm -f "${CHUNKDIR}"/chunk_*.fa

total=$(grep -c '^>' "$PROTEOME")
per=$(( (total + N_CHUNKS - 1) / N_CHUNKS ))     # ceil(total / N_CHUNKS)
echo "[split] ${total} proteins -> ${N_CHUNKS} chunks (~${per} seqs each)"

awk -v per="$per" -v dir="$CHUNKDIR" '
  /^>/ { if (n % per == 0) { c++; f = sprintf("%s/chunk_%d.fa", dir, c) } ; n++ }
  { print > f }
' "$PROTEOME"

echo "$total" > "${OUTDIR}/.nseq"                # remembered for hmmsearch -Z
echo "[split] wrote $(ls "${CHUNKDIR}"/chunk_*.fa | wc -l) chunk(s) in ${CHUNKDIR}"
