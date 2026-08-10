#!/bin/bash
# ============================================================================
# 00_check_requirements.sh — verify software + input files before running.
# Usage:  ./00_check_requirements.sh <species>      (sugarcane | purple)
# Exits non-zero if anything essential is missing.
# ============================================================================
set -uo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${here}/config.sh" "${1:?usage: 00_check_requirements.sh <species>}"

ok=0; bad=0
chk_cmd () { if command -v "$1" >/dev/null 2>&1; then printf "  [ok]  %-14s %s\n" "$1" "$(command -v "$1")"; ok=$((ok+1)); else printf "  [--]  %-14s MISSING\n" "$1"; bad=$((bad+1)); fi; }
chk_file() { if [ -s "$2" ]; then printf "  [ok]  %-14s %s\n" "$1" "$2"; ok=$((ok+1)); else printf "  [--]  %-14s MISSING: %s\n" "$1" "$2"; bad=$((bad+1)); fi; }

echo "=== software ==="
chk_cmd hmmsearch
chk_cmd perl
chk_cmd awk
chk_cmd sort
echo "    hmmsearch: $(hmmsearch -h 2>/dev/null | sed -n 2p | sed 's/^# *//')"
echo "    cores available (nproc): $(nproc)   |  plan: ${MAX_JOBS} jobs x ${CPU_PER_JOB} cpu = $((MAX_JOBS*CPU_PER_JOB))"

echo "=== shared pipeline files ==="
chk_file "HMM library" "$HMM_DB"
chk_file "rules"       "$RULES"
chk_file "assign.pl"   "$ASSIGN_PL"

echo "=== ${SPECIES} inputs ==="
chk_file "proteome"      "$PROTEOME"
chk_file "node metrics"  "$NODE_METRICS"
if [ -s "$PROTEOME" ]; then
  echo "    proteins in proteome : $(grep -c '^>' "$PROTEOME")"
  echo "    nodes in network     : $(($(wc -l < "$NODE_METRICS") - 1))"
fi

echo "----------------------------------------------------------------"
echo "  OK: $ok   MISSING: $bad"
[ "$bad" -eq 0 ] && echo "  -> ready to run run_all.sh ${SPECIES}" || echo "  -> resolve the MISSING items above first"
exit "$bad"
