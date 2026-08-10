#!/bin/bash
# ============================================================================
# 03_assign_families.sh — classify proteins into TF/TAP families.
# Wraps mytfdb/assign_family_membership.pl (RulesFull + domtblout).
# Output: family_assignment.tsv  (protein_id <tab> Family <tab> Type[TFF/OTR/Orphans])
# Usage:  ./03_assign_families.sh <species>
# ============================================================================
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${here}/config.sh" "${1:?usage: 03_assign_families.sh <species>}"

domtbl="${OUTDIR}/concatenated.domtbl"
out="${OUTDIR}/family_assignment.tsv"
[ -s "$domtbl" ] || { echo "ERROR: ${domtbl} missing — run 02_run_hmmsearch.sh"; exit 1; }

perl "$ASSIGN_PL" --pfam "$domtbl" --rules "$RULES" --species "$SPECIES" --out "$out"
echo "[assign] wrote ${out}  ($(($(wc -l < "$out") - 1)) proteins classified)"
