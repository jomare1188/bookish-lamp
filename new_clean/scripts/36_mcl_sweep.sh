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

MCI="${WORK}/${STUDY}.mci"
SW="${WORK}/sweep_${STUDY}"; mkdir -p "$SW"
say() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }
arcs_in() { "$BIN/mcx" query -imx "$1" -t "$CORES" 2>/dev/null | awk 'NR>1{s+=$2} END{print s+0}'; }

[ -s "$MCI" ] || { echo "FATAL: no matrix at $MCI -- run ./run.sh mclload $STUDY" >&2; exit 1; }

N_NODES=$("$BIN/mcx" query -imx "$MCI" --dim 2>&1 | grep -oE '[0-9]+ x' | head -1 | grep -oE '[0-9]+')
N_ARCS=$(arcs_in "$MCI")
say "$STUDY: $N_NODES nodes, $(printf "%'d" "$N_ARCS") arcs"
say "inflation: $I_LIST   (unreduced control: $I_LIST_NONE)"
say "k-NN:      none ${K_LIST}"

TSV="${OUT_DIR}/mcl_sweep_${STUDY}.tsv"
printf 'study\tknn\tknn_mode\tinflation\tn_clusters\tlargest\tlargest_pct\tmedian_size\tsingletons\tefficiency\tmass_fraction\tarea_fraction\tmodularity\tjury_score\tjury_word\truntime_s\n' > "$TSV"

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
    TAG="k${K}.I$(printf '%s' "$I" | tr -d '.')"
    CLS="${SW}/cls.${TAG}"
    ERR="${SW}/mcl.${TAG}.stderr"
    if [ ! -s "$CLS" ]; then
      say "  mcl -I $I  (k=$K)"
      T0=$(date +%s)
      "$BIN/mcl" "$MAT" -I "$I" -te "$CORES" ${SCHEME:+-scheme "$SCHEME"} -o "$CLS" 2> "$ERR" || {
        say "    mcl FAILED -- see $(basename "$ERR")"; continue; }
      echo $(( $(date +%s) - T0 )) > "${CLS}.secs"
    fi
    CLS_LIST+=("$CLS")
  done
  [ ${#CLS_LIST[@]} -gt 0 ] || { say "  no clusterings for k=$K"; continue; }

  # clm info reports every cell against the SAME matrix the clustering came from
  say "  clm info"
  INFO="${SW}/info.k${K}.txt"
  "$BIN/clm" info "$MAT" "${CLS_LIST[@]}" > "$INFO" 2>/dev/null

  for I in $THIS_I; do
    TAG="k${K}.I$(printf '%s' "$I" | tr -d '.')"
    CLS="${SW}/cls.${TAG}"; ERR="${SW}/mcl.${TAG}.stderr"
    [ -s "$CLS" ] || continue
    LINE=$(grep -F "src=$CLS" "$INFO" | head -1)
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
    JURY=$(grep -m1 'jury pruning synopsis' "$ERR" 2>/dev/null | sed 's/.*synopsis: <//; s/>.*//' || true)
    JSCORE=$(printf '%s' "$JURY" | sed -n 's/^\([0-9.]*\).*/\1/p')
    JWORD=$(printf '%s' "$JURY" | sed -n 's/^[0-9.]* or \(.*\)$/\1/p')
    MED=$(sizes_of "$CLS" | awk '{a[NR]=$1} END{print (NR? a[int((NR+1)/2)] : 0)}')
    SEC=$(cat "${CLS}.secs" 2>/dev/null || echo NA)
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%.3f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$STUDY" "$K" "$MODE" "$I" "${NCL:-NA}" "${MAX:-NA}" \
      "$(awk -v m="${MAX:-0}" -v n="$N_NODES" 'BEGIN{print 100*m/n}')" \
      "${MED:-NA}" "${SGL:-NA}" "${EFF:-NA}" "${MF:-NA}" "${AF:-NA}" "${MOD:-NA}" \
      "${JSCORE:-NA}" "${JWORD:-NA}" "$SEC" >> "$TSV"
  done

  # Where does the partition stop moving? Distances between consecutive inflations.
  if [ ${#CLS_LIST[@]} -gt 1 ]; then
    "$BIN/clm" dist --chain "${CLS_LIST[@]}" > "${SW}/dist.k${K}.txt" 2>/dev/null || true
  fi
done

say ""
say "wrote $(basename "$TSV")"
column -t "$TSV" | sed -n '1,60p'
say ""
say "area fraction (af) is the giant-module statistic; jury grades the pruning."
say "done: $STUDY"
