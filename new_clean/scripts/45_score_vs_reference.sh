#!/usr/bin/env bash
# =============================================================================
# 45_score_vs_reference.sh -- score every k's partition against ONE graph.
#
# WHY THIS EXISTS. 36_mcl_sweep.sh calls `clm info` with the REDUCED matrix, so
# each k's modularity is measured on its own graph. Comparing that column across
# k compares "best partition of graph A" with "best partition of graph B" -- not
# two partitions of one network. Newman's Q also drifts upward with sparsity, the
# same bias that made the Barabasi-Albert criterion useless here, so a modularity
# that rises as k falls proves nothing on its own.
#
# Here every clustering is scored against the SAME reference: the unreduced
# Pearson-only network. Only the partition changes, so the numbers are
# comparable and an argmax means something.
#
# WHAT THE COLUMNS MEAN, all from clm info against the reference graph:
#   eff  efficiency -- captured edge mass balanced against cluster footprint
#   mod  Newman modularity
#   mf   mass fraction -- share of the ORIGINAL network's edge weight that falls
#        inside modules. This is the honest "how much of the real network does
#        this partition explain" number, and a reduction that throws away signal
#        will show it here.
#   af   area fraction -- sum of squared cluster sizes / N^2, the giant-module
#        statistic
#
# RUN: through run.sh  ->  RESULTS=... ./run.sh scoreref sugarcane
# =============================================================================
set -euo pipefail

STUDY="${CLEAN_STUDY:?}"; WORK="${CLEAN_WORK_DIR:?}"; BIN="${CLEAN_MCL_BIN_DIR:?}"
OUT_DIR="${CLEAN_OUT_DIR:?}"
REF="${CLEAN_REF:-${WORK}/${STUDY}.mci}"
SW="${WORK}/sweep_${STUDY}"
say() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }

[ -s "$REF" ] || { echo "FATAL: no reference matrix at $REF" >&2; exit 1; }
mapfile -t CLS < <(ls "$SW"/cls.k*.I* 2>/dev/null | grep -v '\.secs$' | sort -V)
[ "${#CLS[@]}" -gt 0 ] || { echo "FATAL: no clusterings in $SW" >&2; exit 1; }

say "$STUDY: scoring ${#CLS[@]} partitions against $(basename "$REF")"
say "  every partition measured on the SAME graph, so the columns are comparable"

INFO="${SW}/info_vs_reference.txt"
"$BIN/clm" info "$REF" "${CLS[@]}" > "$INFO" 2>/dev/null

TSV="${OUT_DIR}/mcl_vs_reference_${STUDY}.tsv"
printf 'study\tk\tinflation\tn_clusters\tlargest\tsingletons\teff_ref\tmod_ref\tmass_frac_ref\tarea_frac_ref\n' > "$TSV"
while read -r line; do
  case "$line" in *src=*) ;; *) continue ;; esac
  f=$(sed -n 's/.*src=\([^ ]*\).*/\1/p' <<<"$line"); b=$(basename "$f")
  k=$(sed -n 's/^cls\.k\([^.]*\)\..*/\1/p' <<<"$b")
  i=$(sed -n 's/^cls\.k[^.]*\.I\(.*\)$/\1/p' <<<"$b")
  # I14 -> 1.4, I2 -> 2
  case "$i" in ??) i="${i:0:1}.${i:1:1}";; esac
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$STUDY" "$k" "$i" \
    "$(sed -n 's/.*ncl=\([0-9]*\).*/\1/p' <<<"$line")" \
    "$(sed -n 's/.*max=\([0-9]*\).*/\1/p' <<<"$line")" \
    "$(sed -n 's/.*sgl=\([0-9]*\).*/\1/p' <<<"$line")" \
    "$(sed -n 's/.*eff=\([0-9.-]*\).*/\1/p' <<<"$line")" \
    "$(sed -n 's/.*mod=\([0-9.-]*\).*/\1/p' <<<"$line")" \
    "$(sed -n 's/.*mf=\([0-9.-]*\).*/\1/p' <<<"$line")" \
    "$(sed -n 's/.*af=\([0-9.-]*\).*/\1/p' <<<"$line")" >> "$TSV"
done < "$INFO"

say "wrote $(basename "$TSV")  ($(( $(wc -l < "$TSV") - 1 )) partitions)"
column -t "$TSV"
say "done: $STUDY"
