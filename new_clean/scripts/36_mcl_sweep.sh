#!/usr/bin/env bash
# =============================================================================
# 36_mcl_sweep.sh -- granularity as a surface, with MCL's own diagnostics kept.
#
# The pipeline has always run mcl at exactly one setting (-I 2, no transform) and
# thrown its stderr away. This runs a k-NN x inflation grid and keeps everything
# the tool says about it.
#
# WHY THE JURY SYNOPSIS MATTERS AND WHY IT WAS NEVER SEEN. mcl prunes each node's
# neighbour list during computation; the default -scheme 7 keeps at most ~1,200
# of them. Purple's top decile of nodes has >33,000 neighbours, so ~96% of each
# hub's edges are discarded on the fly, among weights that are nearly tied. mcl
# grades that pruning itself and prints the grade to STDERR -- which
# 05_mcl_clustering.r:67 sends to /dev/null. A low grade means the clustering
# rests on an arbitrary subset of the graph. Every cell's grade is recorded here.
#
# WHAT DECIDES. clm info, from the author's own protocol (see clmdist(1)
# EXAMPLES):
#   eff  efficiency -- balances capturing edge mass against cluster footprint
#   mf   mass fraction -- share of edge weight inside clusters
#   af   AREA FRACTION -- sum of squared cluster sizes over N^2. This is the
#        giant-module statistic: a module holding 28% of the network dominates
#        it, so it is the number to watch and it needs no bespoke metric.
#   mod  modularity, computed by clm on the same footing for every cell
# and clm dist, which says where the partition stops moving as inflation rises.
#
# RUN: through run.sh  ->  ./run.sh mclsweep sugarcane "180 240"
# =============================================================================
set -euo pipefail

STUDY="${CLEAN_STUDY:?}"; WORK="${CLEAN_WORK_DIR:?}"; BIN="${CLEAN_MCL_BIN_DIR:?}"
OUT_DIR="${CLEAN_OUT_DIR:?}"; CORES="${CLEAN_CORES:-16}"
I_LIST="${CLEAN_SWEEP_I:-1.4 2 3 4 6}"
# The UNREDUCED control gets a shorter ladder, for cost rather than principle.
# mcl on the full matrix is ~40x the work of a k-NN-reduced one (sugarcane keeps
# 2.4% of its arcs at k=180), and low inflation is the slowest end -- the -I 1.4
# cell was still raising chaos after five iterations. Three points spanning the
# range are enough for the control to answer its one question: does inflation
# alone flatten the giant module? If it does, that will be visible at 2/4/6.
I_LIST_NONE="${CLEAN_SWEEP_I_NONE:-2 4 6}"
K_LIST="${CLEAN_SWEEP_K:-}"
SCHEME="${CLEAN_MCL_SCHEME:-}"
# Extra mcl resource arguments, appended verbatim (e.g. "-S 10000").
#
# NOT A TUNING KNOB -- A CONFOUND CONTROL. mcl squares the matrix repeatedly and
# prunes every column each iteration; -scheme 7 (the default AND the highest
# preset) keeps S=1200 neighbours per node. Unpruned purple has mean degree
# 7,946, so ~85% of each list is discarded on the fly among near-tied weights,
# and mcl grades its own pruning "awful". Re-running the winning inflation at a
# higher -S is the only way to tell whether the ladder ranked inflations or
# ranked mcl's pruner. The setting is part of the cell tag so both live in one
# table.
RESOURCE="${CLEAN_MCL_RESOURCE:-}"

MCI="${WORK}/${STUDY}.mci"
SW="${WORK}/sweep_${STUDY}"; mkdir -p "$SW" "$OUT_DIR"
# Every cell, from every method, lands in one manifest so nothing downstream has
# to parse meaning out of a filename. 47 (Leiden/Louvain) appends to the same
# file, which is what lets 48 score all methods in a single clm info call.
MANIFEST="${SW}/partitions.tsv"
[ -s "$MANIFEST" ] || printf 'method\tparam\tresource\tfile\tn_clusters\tlargest\truntime_s\n' > "$MANIFEST"
# A resource setting has to reach the cell tag, or the probe run would overwrite
# the default run's clustering and the comparison would silently vanish.
RTAG=""
if [ -n "$RESOURCE" ]; then
  RTAG=".$(printf '%s' "$RESOURCE" | tr -cd 'A-Za-z0-9')"
fi
say() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }
arcs_in() { "$BIN/mcx" query -imx "$1" -t "$CORES" 2>/dev/null | awk 'NR>1{s+=$2} END{print s+0}'; }

[ -s "$MCI" ] || { echo "FATAL: no matrix at $MCI -- run ./run.sh mclload $STUDY" >&2; exit 1; }

N_NODES=$("$BIN/mcx" query -imx "$MCI" --dim 2>&1 | grep -oE '[0-9]+ x' | head -1 | grep -oE '[0-9]+')
N_ARCS=$(arcs_in "$MCI")
say "$STUDY: $N_NODES nodes, $(printf "%'d" "$N_ARCS") arcs"
say "inflation: $I_LIST   (unreduced control: $I_LIST_NONE)"
say "k-NN:      none ${K_LIST}"

# --- the grid must be legal and unambiguous BEFORE anything runs -------------
# Two failures found the hard way, both silent:
#
#  * TAG COLLISION. The cell tag strips the decimal point, so -I 1.2 and -I 12
#    both become "I12". The second is then seen as already computed, skipped,
#    and reported under the wrong inflation -- including ACROSS runs, which is
#    how it happened: an extension grid met an I12 cell left by an earlier 1.2.
#  * OUT OF RANGE. mcl accepts -I only in (0.0, 30.0]. Above that it prints
#    "float argument to -I should be in range" to STDERR and SILENTLY USES THE
#    DEFAULT 2.0. Measured: -I 40 produced a partition byte-identical to -I 2.
#    Nothing in the output says so; the row just looks like a real -I 40.
check_grid() {
  local seen="" i tag
  for i in $1; do
    awk -v v="$i" 'BEGIN{ exit !(v+0 > 0 && v+0 <= 30) }' || {
      echo "FATAL: -I $i is outside mcl's range (0.0,30.0]; mcl would silently use 2.0" >&2
      exit 1; }
    tag=$(printf '%s' "$i" | tr -d '.')
    case " $seen " in *" $tag "*)
      echo "FATAL: -I $i collides with an earlier value on cell tag 'I$tag'." >&2
      echo "  Values differing only by a decimal point (1.2 vs 12) share a cell." >&2
      exit 1;; esac
    seen="$seen $tag"
  done
}
check_grid "$I_LIST"
check_grid "$I_LIST_NONE"

TSV="${OUT_DIR}/mcl_sweep_${STUDY}.tsv"
# ADDITIVE, not truncating. This file used to be rewritten from scratch on every
# run, so an extension grid destroyed the rows the first grid had produced.
printf 'study\tknn\tknn_mode\tinflation\tn_clusters\tlargest\tlargest_pct\tmedian_size\tsingletons\tefficiency\tmass_fraction\tarea_fraction\tmodularity\tjury_score\tjury_word\truntime_s\tunderflow_vectors\n' > "${TSV}.new"
if [ -s "$TSV" ]; then
  tail -n +2 "$TSV" >> "${TSV}.new"
  say "keeping $(( $(wc -l < "$TSV") - 1 )) existing rows from $(basename "$TSV")"
fi
mv "${TSV}.new" "$TSV"

# Cluster sizes from the native cluster file: records run from the cluster index
# to a '$' terminator and may wrap across lines, so tokens are accumulated.
sizes_of() {
  awk '/^\(mclmatrix/{inm=1} !inm{next} /^begin/{go=1;next} !go{next} /^\)/{exit}
       { for(i=1;i<=NF;i++){ if($i=="$"){ print n-1; n=0 } else n++ } }' "$1" \
  | sort -n
}

for K in none $K_LIST; do
  MODE="none"
  if [ "$K" = "none" ]; then
    MAT="$MCI"
  else
    MAT="${SW}/${STUDY}.knn${K}.mci"
    if [ ! -s "$MAT" ]; then
      say "reducing with #knn($K)"
      "$BIN/mcx" alter -imx "$MCI" -tf "#knn($K)" -o "$MAT" 2>/dev/null
    fi
    # Count arcs by summing node degrees. `mcx query --dim` reports only the
    # matrix dimensions, never an entry count -- reading one out of it silently
    # yields 0 and triggers the fallback below on a perfectly good reduction.
    ENT=$(arcs_in "$MAT")
    MODE="knn"
    # #knn INTERSECTS neighbour lists, so a small k can empty the graph outright.
    # #knnj joins instead; fall back to it rather than clustering nothing.
    if [ "${ENT:-0}" -lt 1000 ]; then
      say "  #knn($K) left only ${ENT} arcs -- falling back to #knnj($K)"
      "$BIN/mcx" alter -imx "$MCI" -tf "#knnj($K)" -o "$MAT" 2>/dev/null
      MODE="knnj"
      ENT=$(arcs_in "$MAT")
    fi
    say "  k=$K ($MODE): $(printf "%'d" "$ENT") of $(printf "%'d" "$N_ARCS") arcs retained ($(awk -v a="$ENT" -v b="$N_ARCS" 'BEGIN{printf "%.1f%%", 100*a/b}'))"
  fi

  THIS_I="$I_LIST"
  [ "$K" = "none" ] && THIS_I="$I_LIST_NONE"
  CLS_LIST=()
  for I in $THIS_I; do
    TAG="k${K}.I$(printf '%s' "$I" | tr -d '.')${RTAG}"
    CLS="${SW}/cls.${TAG}"
    ERR="${SW}/mcl.${TAG}.stderr"
    # A cell left by an earlier run with a DIFFERENT inflation must never be
    # reused. check_grid cannot see across runs; this sidecar can.
    if [ -s "$CLS" ] && [ -s "${CLS}.inflation" ]; then
      PREV=$(cat "${CLS}.inflation")
      if [ "$PREV" != "$I" ]; then
        echo "FATAL: $(basename "$CLS") was produced with -I $PREV, not -I $I." >&2
        echo "  Reusing it would report one inflation's clustering under another." >&2
        exit 1
      fi
    fi
    if [ ! -s "$CLS" ]; then
      say "  mcl -I $I  (k=$K)"
      T0=$(date +%s)
      # shellcheck disable=SC2086 -- RESOURCE is a deliberate multi-word option list
      "$BIN/mcl" "$MAT" -I "$I" -te "$CORES" ${SCHEME:+-scheme "$SCHEME"} $RESOURCE -o "$CLS" 2> "$ERR" || {
        say "    mcl FAILED -- see $(basename "$ERR")"; continue; }
      echo $(( $(date +%s) - T0 )) > "${CLS}.secs"
      printf '%s' "$I" > "${CLS}.inflation"
      printf '%s\t%s\t%s\t%s\t\t\t%s\n' "mcl" "$I" "${RESOURCE:-default}" "$CLS" \
        "$(cat "${CLS}.secs")" >> "$MANIFEST"
    fi
    CLS_LIST+=("$CLS")
  done
  [ ${#CLS_LIST[@]} -gt 0 ] || { say "  no clusterings for k=$K"; continue; }

  # clm info reports every cell against the SAME matrix the clustering came from
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
  say "  clm info (one call per partition)"
  INFO="${SW}/info.k${K}.txt"
  : > "$INFO"
  for c in "${CLS_LIST[@]}"; do
    "$BIN/clm" info "$MAT" "$c" >> "$INFO" 2>/dev/null
  done

  for I in $THIS_I; do
    TAG="k${K}.I$(printf '%s' "$I" | tr -d '.')${RTAG}"
    CLS="${SW}/cls.${TAG}"; ERR="${SW}/mcl.${TAG}.stderr"
    [ -s "$CLS" ] || continue
    # Exact field match: "src=<path>" is a prefix of "src=<path>0", so a plain
    # grep for cls.knone.I4 also matches cls.knone.I40.
    LINE=$(awk -v want="src=$CLS" '{for(i=1;i<=NF;i++) if($i==want){print;exit}}' "$INFO")
    EFF=$(sed -n 's/.*eff=\([0-9.]*\).*/\1/p' <<<"$LINE")
    MOD=$(sed -n 's/.*mod=\([0-9.-]*\).*/\1/p' <<<"$LINE")
    MF=$(sed -n  's/.*mf=\([0-9.]*\).*/\1/p'  <<<"$LINE")
    AF=$(sed -n  's/.*af=\([0-9.]*\).*/\1/p'  <<<"$LINE")
    NCL=$(sed -n 's/.*ncl=\([0-9]*\).*/\1/p'  <<<"$LINE")
    MAX=$(sed -n 's/.*max=\([0-9]*\).*/\1/p'  <<<"$LINE")
    SGL=$(sed -n 's/.*sgl=\([0-9]*\).*/\1/p'  <<<"$LINE")
    # The jury synopsis is the whole reason stderr is kept.
    # mcl prints e.g. "<39.2 or deplorable>". Keep the number and the word in
    # SEPARATE columns -- the word contains spaces ("off colour"), which would
    # shift every later field of anything reading this as whitespace-delimited.
    # High inflation can underflow: mcl prints "nonpositive sum <0.000000>" and
    # returns a partition that looks plausible and is numerically meaningless.
    # Measured at -I 20 on sugarcane: a 34.8% giant with jury 12.5 "wretched".
    NDEGEN=$(grep -c 'nonpositive sum' "$ERR" 2>/dev/null || true)
    if [ "${NDEGEN:-0}" -gt 0 ]; then
      say "    WARNING: -I $I underflowed ($NDEGEN nonpositive-sum vectors); row flagged"
    fi
    JURY=$(grep -m1 'jury pruning synopsis' "$ERR" 2>/dev/null | sed 's/.*synopsis: <//; s/>.*//' || true)
    JSCORE=$(printf '%s' "$JURY" | sed -n 's/^\([0-9.]*\).*/\1/p')
    JWORD=$(printf '%s' "$JURY" | sed -n 's/^[0-9.]* or \(.*\)$/\1/p')
    MED=$(sizes_of "$CLS" | awk '{a[NR]=$1} END{print (NR? a[int((NR+1)/2)] : 0)}')
    SEC=$(cat "${CLS}.secs" 2>/dev/null || echo NA)
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%.3f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$STUDY" "$K" "$MODE" "$I" "${NCL:-NA}" "${MAX:-NA}" \
      "$(awk -v m="${MAX:-0}" -v n="$N_NODES" 'BEGIN{print 100*m/n}')" \
      "${MED:-NA}" "${SGL:-NA}" "${EFF:-NA}" "${MF:-NA}" "${AF:-NA}" "${MOD:-NA}" \
      "${JSCORE:-NA}" "${JWORD:-NA}" "$SEC" "${NDEGEN:-0}" >> "$TSV"
  done

  # Where does the partition stop moving? Distances between consecutive inflations.
  if [ ${#CLS_LIST[@]} -gt 1 ]; then
    "$BIN/clm" dist --chain "${CLS_LIST[@]}" > "${SW}/dist.k${K}.txt" 2>/dev/null || true
  fi
done

say ""
# Keep the LAST row for each (knn, inflation): a re-run supersedes its own
# earlier row rather than appearing twice.
awk -F'\t' 'NR==1{print;next}{k=$2 FS $4; row[k]=$0; if(!(k in ord)) ord[++n]=k}
     END{for(i=1;i<=n;i++) print row[ord[i]]}' "$TSV" > "${TSV}.dedup" && mv "${TSV}.dedup" "$TSV"
say "wrote $(basename "$TSV")"
column -t "$TSV" | sed -n '1,60p'
say ""
say "area fraction (af) is the giant-module statistic; jury grades the pruning."
say "done: $STUDY"
