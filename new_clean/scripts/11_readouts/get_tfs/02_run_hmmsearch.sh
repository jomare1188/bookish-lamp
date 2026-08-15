#!/bin/bash
# ============================================================================
# 02_run_hmmsearch.sh — hmmsearch every chunk in parallel, then concatenate.
#
# hmmsearch flags mirror the original run_tf.sh:
#   --cut_ga     use each model's Gathering-threshold (bit-score) cutoff
#   --domtblout  per-domain tabular output (input to assign_family_membership.pl)
#   -Z <N>       set search space to the FULL proteome size so E-values are
#                comparable across chunks (advised in mytfdb/README when splitting)
#
# Resumable: a chunk whose .domtbl already exists is skipped; partial writes go
# to .tmp and are only renamed on success.
# Usage:  ./02_run_hmmsearch.sh <species>
# ============================================================================
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${here}/config.sh" "${1:?usage: 02_run_hmmsearch.sh <species>}"

[ -d "$CHUNKDIR" ] || { echo "ERROR: run 01_split_proteome.sh first"; exit 1; }
mkdir -p "$DOMTBLDIR"
Z=$(cat "${OUTDIR}/.nseq")
echo "[hmm] DB=$(basename "$HMM_DB")  Z=${Z}  parallel=${MAX_JOBS}x${CPU_PER_JOB}cpu"

run_one () {
  local chunk="$1"
  local base; base=$(basename "$chunk" .fa)
  local dom="${DOMTBLDIR}/${base}.domtbl"
  if [ -s "$dom" ]; then echo "[hmm] skip ${base} (done)"; return 0; fi
  hmmsearch --cpu "$CPU_PER_JOB" --cut_ga -Z "$Z" \
            --domtblout "${dom}.tmp" -o /dev/null "$HMM_DB" "$chunk" \
    && mv "${dom}.tmp" "$dom" \
    && echo "[hmm] done ${base}"
}

for chunk in "${CHUNKDIR}"/chunk_*.fa; do
  run_one "$chunk" &
  while [ "$(jobs -r | wc -l)" -ge "$MAX_JOBS" ]; do wait -n; done
done
wait

n_done=$(ls "${DOMTBLDIR}"/chunk_*.domtbl 2>/dev/null | wc -l)
n_chunks=$(ls "${CHUNKDIR}"/chunk_*.fa | wc -l)
[ "$n_done" -eq "$n_chunks" ] || { echo "ERROR: only ${n_done}/${n_chunks} chunks finished — rerun to resume"; exit 1; }

cat "${DOMTBLDIR}"/chunk_*.domtbl > "${OUTDIR}/concatenated.domtbl"
echo "[hmm] all ${n_done} chunks -> ${OUTDIR}/concatenated.domtbl"
