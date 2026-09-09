#!/usr/bin/env bash
# =============================================================================
# 48_cluster_compare.sh -- every method, every setting, scored on ONE graph.
#
# WHY THIS IS THE WHOLE POINT. 36 calls `clm info` with whatever matrix the
# clustering came from, which is right when the matrix never changes and wrong
# the moment it does. The k-NN experiment was decided by exactly that mistake:
# each k's modularity was computed on its own reduced graph, "structure" appeared
# to improve as the graph got sparser, and on one fixed graph the finding
# inverted (docs/results.md). Here there is one reference -- the unpruned
# Pearson-only matrix -- and only the partition varies.
#
# It also puts Leiden and MCL in the same table. clm info takes any cluster file
# against any graph, so once 47 has written Leiden's partitions in mcl's format
# (and verify_cluster_roundtrip.py has proved that writer correct), clm info
# scores both on identical definitions of eff/mf/af/mod. No per-method metric
# implementation, nothing to get subtly inconsistent. One call PER PARTITION --
# see the note above the scoring loop, which is not optional.
#
# THE COLUMNS
#   eff  efficiency -- captured edge mass balanced against cluster footprint.
#        THE PRIMARY CRITERION. It is the one clm's author designed to be
#        maximised; mf and af each have a trivial optimum, eff does not.
#   mf   mass fraction -- share of the graph's edge weight falling inside
#        clusters. REPORTED, NEVER USED ALONE: the one-cluster partition scores
#        mf = 1.0. That baseline is emitted into this table on purpose so the
#        trap is visible rather than explained.
#   af   area fraction -- sum of squared cluster sizes / N^2, the giant-module
#        statistic. Also trivially optimised, in the other direction: all
#        singletons scores ~0.
#   mod  Newman modularity, computed by clm on the same footing for every row.
#   jury mcl's own grade of its INTERNAL pruning (mcl rows only). Leiden prunes
#        nothing, so the column is empty for it -- that asymmetry is the point.
#
# RUN: through run.sh  ->  ./run.sh clustercompare sugarcane
# =============================================================================
set -euo pipefail

STUDY="${CLEAN_STUDY:?}"; WORK="${CLEAN_WORK_DIR:?}"; BIN="${CLEAN_MCL_BIN_DIR:?}"
OUT_DIR="${CLEAN_OUT_DIR:?}"
REF="${CLEAN_REF:-${WORK}/${STUDY}.mci}"
SW="${WORK}/sweep_${STUDY}"
MANIFEST="${SW}/partitions.tsv"
say() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }

[ -s "$REF" ]      || { echo "FATAL: no reference matrix at $REF" >&2; exit 1; }
[ -s "$MANIFEST" ] || { echo "FATAL: no manifest at $MANIFEST -- run mclladder / leidensweep first" >&2; exit 1; }
mkdir -p "$OUT_DIR"

N_NODES=$("$BIN/mcx" query -imx "$REF" --dim 2>&1 | grep -oE '[0-9]+ x' | head -1 | grep -oE '[0-9]+')
say "$STUDY: reference $(basename "$REF"), $(printf "%'d" "$N_NODES") nodes"

# --- the baseline that makes the mass-fraction trap visible ------------------
# One cluster holding everything: mf = 1.0, af = 1.0. Any argument that picks a
# setting on mass fraction alone picks this row, and it is not a clustering.
BASE="${SW}/cls.baseline.onecluster"
if [ ! -s "$BASE" ]; then
  say "writing the one-cluster baseline"
  { printf '(mclheader\nmcltype matrix\ndimensions %dx1\n)\n(mclmatrix\nbegin\n0 ' "$N_NODES"
    seq 0 $(( N_NODES - 1 )) | tr '\n' ' '
    printf '$\n)\n'; } > "$BASE"
fi

# --- collect every partition named in the manifest ---------------------------
mapfile -t FILES < <(awk -F'\t' 'NR>1 && $4 != "" {print $4}' "$MANIFEST" | sort -u)
FILES+=("$BASE")
KEEP=(); MISSING=0
for f in "${FILES[@]}"; do
  if [ -s "$f" ]; then KEEP+=("$f"); else say "  MISSING, skipped: $f"; MISSING=$(( MISSING + 1 )); fi
done
[ "${#KEEP[@]}" -gt 0 ] || { echo "FATAL: no partition files exist" >&2; exit 1; }
say "scoring ${#KEEP[@]} partitions against the reference ($MISSING missing)"

# ---------------------------------------------------------------------------
# ONE CLUSTERING PER clm info CALL. THIS IS NOT A STYLE CHOICE.
#
# `clm info <graph> <cls1> <cls2> ...` is documented to accept many clusterings
# at once, but its eff and mf depend on WHICH OTHERS are in the call. Measured on
# sugarcane, same matrix, same cluster file:
#
#     cls.knone.I6 alone            -> eff=0.47281  mf=0.51576
#     cls.knone.I6 in a batch of 8  -> eff=0.37821  mf=0.46878
#
# mod and af are unaffected. Scored one at a time, eff rises monotonically across
# the whole inflation ladder and mf falls monotonically; scored in batches both
# columns develop discontinuities that look like real structure and are not.
# Costs one clm info invocation per partition. Pay it.
# ---------------------------------------------------------------------------
INFO="${SW}/info_all_methods.txt"
: > "$INFO"
n=0
for c in "${KEEP[@]}"; do
  n=$(( n + 1 ))
  [ $(( n % 10 )) -eq 0 ] && say "  scored $n/${#KEEP[@]}"
  "$BIN/clm" info "$REF" "$c" >> "$INFO" 2>/dev/null
done

# --- join clm info back onto the manifest ------------------------------------
TSV="${OUT_DIR}/cluster_methods_${STUDY}.tsv"
printf 'study\tmethod\tparam\tresource\tn_clusters\tlargest\tlargest_pct\tsingletons\teff\tmf\taf\tmod\tjury_score\tjury_word\truntime_s\n' > "$TSV"

while read -r line; do
  case "$line" in *src=*) ;; *) continue ;; esac
  f=$(awk '{for(i=1;i<=NF;i++) if($i ~ /^src=/){sub(/^src=/,"",$i); print $i; exit}}' <<<"$line")

  if [ "$f" = "$BASE" ]; then
    METHOD="baseline"; PARAM="one cluster"; RES=""; SEC="NA"
  else
    # The manifest is the source of truth for what a file IS. Nothing here parses
    # meaning out of a filename -- that is how the k-NN collector ended up with a
    # stale ERROR row it could not distinguish from a real one.
    ROW=$(awk -F'\t' -v want="$f" 'NR>1 && $4==want {print; exit}' "$MANIFEST")
    [ -n "$ROW" ] || continue
    METHOD=$(cut -f1 <<<"$ROW"); PARAM=$(cut -f2 <<<"$ROW")
    RES=$(cut -f3 <<<"$ROW");    SEC=$(cut -f7 <<<"$ROW")
  fi

  # mcl's jury grade lives in the per-cell stderr, and only mcl has one.
  JSCORE=""; JWORD=""
  ERR="${f/\/cls./\/mcl.}.stderr"
  if [ -s "$ERR" ]; then
    JURY=$(grep -m1 'jury pruning synopsis' "$ERR" 2>/dev/null | sed 's/.*synopsis: <//; s/>.*//' || true)
    JSCORE=$(printf '%s' "$JURY" | sed -n 's/^\([0-9.]*\).*/\1/p')
    JWORD=$(printf '%s' "$JURY" | sed -n 's/^[0-9.]* or \(.*\)$/\1/p')
  fi

  MAX=$(sed -n 's/.*max=\([0-9]*\).*/\1/p' <<<"$line")
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%.3f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$STUDY" "$METHOD" "$PARAM" "${RES:-}" \
    "$(sed -n 's/.*ncl=\([0-9]*\).*/\1/p' <<<"$line")" "${MAX:-NA}" \
    "$(awk -v m="${MAX:-0}" -v n="$N_NODES" 'BEGIN{print 100*m/n}')" \
    "$(sed -n 's/.*sgl=\([0-9]*\).*/\1/p' <<<"$line")" \
    "$(sed -n 's/.*eff=\([0-9.-]*\).*/\1/p' <<<"$line")" \
    "$(sed -n 's/.*mf=\([0-9.-]*\).*/\1/p'  <<<"$line")" \
    "$(sed -n 's/.*af=\([0-9.-]*\).*/\1/p'  <<<"$line")" \
    "$(sed -n 's/.*mod=\([0-9.-]*\).*/\1/p' <<<"$line")" \
    "${JSCORE:-NA}" "${JWORD:-NA}" "${SEC:-NA}" >> "$TSV"
done < "$INFO"

say "wrote $(basename "$TSV")  ($(( $(wc -l < "$TSV") - 1 )) rows)"
say ""
( head -1 "$TSV"; tail -n +2 "$TSV" | sort -t$'\t' -k9,9gr ) | column -t -s$'\t'
say ""
say "sorted by eff, the primary criterion. Note the baseline row: mf = 1.0 for a"
say "single cluster, which is why mf never decides anything on its own."
say "done: $STUDY"
