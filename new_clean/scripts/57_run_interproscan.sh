#!/usr/bin/env bash
# =============================================================================
# 57_run_interproscan.sh -- the 12 InterPro member databases that were never run.
#
# WHAT IS MISSING. The nf-core/proteinannotator run used 5 of 17 member DBs --
# PANTHER, TIGRFAM(ncbifam), PIRSF, Hamap and SFLD. The local 5.78 install also
# carries Gene3D, SUPERFAMILY, CDD, SMART, PRINTS, PROSITE, Pfam, AntiFam, PIRSR,
# Phobius and TMHMM. Every signature that belongs to an InterPro entry yields GO
# through interpro2go, so the DBs never run are the largest untapped GO source we
# have, and all of it is local -- 49 GB of data, no downloads, no external service.
#
# WHY CHUNKS RATHER THAN ONE BIG JOB. interproscan.properties sets 6-8 embedded
# workers, so a single invocation will not use 256 cores no matter what -cpu says.
# Many concurrent small jobs will. Chunking also makes the run RESUMABLE, which
# matters when the whole thing takes many hours: a chunk that finishes is marked
# done and never repeats, so an interruption costs one chunk rather than the run.
#
# LOSSLESSNESS IS CHECKED, NOT ASSUMED. The obvious failure of a chunked run is a
# chunk silently missing from the merge. Every input sequence must land in exactly
# one chunk, and the merged output's protein count is compared against the input.
#
# GO IS NOT TAKEN FROM -goterms. It is derived later from the IPR column using our
# own pinned interpro2go, so the mapping version stays under our control. -goterms
# is still requested, as a cross-check that the two agree.
#
# RUN: through run.sh  ->  ./run.sh ipscan sugarcane
# =============================================================================
set -euo pipefail

STUDY="${CLEAN_STUDY:?}"
FASTA="${CLEAN_PROTEOME:?}"
OUTDIR="${CLEAN_OUT_DIR:?}"
IPS="${CLEAN_IPS_BIN:-/dados04/jorge/databases/InterProScan/interproscan-5.78-109.0/interproscan.sh}"
CHUNK="${CLEAN_IPS_CHUNK:-2000}"        # sequences per chunk
PAR="${CLEAN_IPS_PARALLEL:-28}"         # concurrent chunks
CPU="${CLEAN_IPS_CPU:-8}"               # -cpu per chunk (matches embedded workers)
APPS="${CLEAN_IPS_APPS:-}"              # empty = all available
TMP="${CLEAN_IPS_TMP:-/dados04/jorge/tmp/ips_${STUDY}}"

say() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }
[ -s "$FASTA" ] || { echo "FATAL: no proteome at $FASTA" >&2; exit 1; }
[ -x "$IPS" ]   || { echo "FATAL: interproscan.sh not executable at $IPS" >&2; exit 1; }
mkdir -p "$OUTDIR"/{chunks,out,done,logs} "$TMP"

N_IN=$(grep -c '^>' "$FASTA")
say "$STUDY: $(printf "%'d" "$N_IN") proteins, $CHUNK per chunk, $PAR concurrent x $CPU cpu"

# --- split once, deterministically -------------------------------------------
if [ ! -f "$OUTDIR/chunks.total" ]; then
  say "splitting"
  awk -v d="$OUTDIR/chunks" -v n="$CHUNK" '
    /^>/ { if (c % n == 0) { if (f) close(f); f = sprintf("%s/chunk_%05d.fa", d, int(c/n)) } c++ }
    { print > f } END { printf "%d\n", int((c + n - 1)/n) }' "$FASTA" > "$OUTDIR/chunks.total"
  # every sequence in exactly one chunk, or the merge below is meaningless
  N_CH=$(cat "$OUTDIR"/chunks/*.fa | grep -c '^>')
  [ "$N_CH" = "$N_IN" ] || { echo "FATAL: split holds $N_CH of $N_IN sequences" >&2; exit 1; }
  say "  $(cat "$OUTDIR/chunks.total") chunks, $(printf "%'d" "$N_CH") sequences -- lossless"
fi
TOTAL=$(cat "$OUTDIR/chunks.total")

# --- run, skipping anything already finished ---------------------------------
run_chunk() {
  local f="$1" b; b=$(basename "$f" .fa)
  [ -f "$OUTDIR/done/$b.done" ] && return 0
  if "$IPS" -i "$f" -f TSV -iprlookup -goterms -cpu "$CPU" \
       -T "$TMP/$b" -o "$OUTDIR/out/$b.tsv" ${APPS:+-appl "$APPS"} \
       > "$OUTDIR/logs/$b.log" 2>&1; then
    touch "$OUTDIR/done/$b.done"
  else
    echo "CHUNK FAILED: $b (see $OUTDIR/logs/$b.log)" >&2
  fi
  rm -rf "$TMP/$b"
}
export -f run_chunk; export OUTDIR IPS CPU TMP APPS

say "running (already done: $(ls "$OUTDIR"/done/*.done 2>/dev/null | wc -l)/$TOTAL)"
ls "$OUTDIR"/chunks/*.fa | xargs -P "$PAR" -I{} bash -c 'run_chunk "$@"' _ {}

DONE=$(ls "$OUTDIR"/done/*.done 2>/dev/null | wc -l)
say "chunks finished: $DONE / $TOTAL"
if [ "$DONE" != "$TOTAL" ]; then
  say "NOT merging -- $(( TOTAL - DONE )) chunk(s) missing. Re-run to resume."
  exit 1
fi

# --- merge, and prove nothing was lost ----------------------------------------
MERGED="$OUTDIR/${STUDY}.interproscan_full.tsv"
cat "$OUTDIR"/out/*.tsv > "$MERGED"
N_OUT=$(cut -f1 "$MERGED" | sort -u | wc -l)
say ""
say "merged: $(printf "%'d" "$(wc -l < "$MERGED")") rows, $(printf "%'d" "$N_OUT") proteins with >=1 signature"
say "  of $(printf "%'d" "$N_IN") input proteins ($(awk -v a="$N_OUT" -v b="$N_IN" 'BEGIN{printf "%.1f", 100*a/b}')%)"
say "  member DBs: $(cut -f4 "$MERGED" | sort -u | tr '\n' ' ')"
say "  with an InterPro accession: $(awk -F'\t' '$12!="-" && $12!=""{n++} END{print n+0}' "$MERGED") rows"
say "done: $STUDY"
