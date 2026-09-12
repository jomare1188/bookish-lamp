#!/usr/bin/env bash
# =============================================================================
# 59_run_pannzer.sh -- PANNZER2 GO prediction, chunked and resumable.
#
# WHAT PANNZER2 ADDS. It is sequence-similarity based but statistical rather than
# domain-transfer: it weights GO terms from homologues by taxonomic distance and
# scores each with ARGOT, returning a PPV per term. So unlike the InterPro tier it
# can annotate proteins with no recognisable domain, and unlike the curated tier it
# is thresholdable -- the PPV is what makes a prediction tier usable at all.
#
# THE CONSTRAINT THAT SHAPES THIS SCRIPT. The standalone package ships no database:
# it queries the SANSparallel and DictServer instances at Helsinki over the
# network. The compute is THEIRS, not ours, so "we have 256 cores" does not apply
# and hammering it with 30 parallel streams would be both rude and rate-limited.
# Measured throughput on our own sequences is 5.1 seq/s in a single stream
# (200 proteins in 38.9 s), so the whole job is ~24 h serial; a small number of
# concurrent streams brings that down without abusing a public academic service.
# CLEAN_PZ_PARALLEL defaults to 4 deliberately.
#
# TWO UPSTREAM BUGS had to be patched to run at all under Python 3, both recorded
# in the install notes: operators/Taxonomy.py line 50 is indented one space deeper
# than the line above it (Python 2 tolerated it), and Clustering.py uses
# numpy.int, removed in NumPy 1.24. Originals are kept as *.orig.
#
# RUN: through run.sh  ->  ./run.sh pannzer sugarcane
# =============================================================================
set -euo pipefail

STUDY="${CLEAN_STUDY:?}"
FASTA="${CLEAN_PROTEOME:?}"
OUTDIR="${CLEAN_OUT_DIR:?}"
PZDIR="${CLEAN_PZ_DIR:?}"
PY="${CLEAN_PZ_PYTHON:?}"
SPECIES="${CLEAN_PZ_SPECIES:-Saccharum officinarum}"
CHUNK="${CLEAN_PZ_CHUNK:-2000}"
PAR="${CLEAN_PZ_PARALLEL:-4}"

say() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }
[ -s "$FASTA" ] || { echo "FATAL: no proteome at $FASTA" >&2; exit 1; }
mkdir -p "$OUTDIR"/{chunks,go,done,logs}

N_IN=$(grep -c '^>' "$FASTA")
if [ ! -f "$OUTDIR/chunks.total" ]; then
  say "splitting $(printf "%'d" "$N_IN") proteins into chunks of $CHUNK"
  awk -v d="$OUTDIR/chunks" -v n="$CHUNK" '
    /^>/ { if (c % n == 0) { if (f) close(f); f = sprintf("%s/c_%05d.fa", d, int(c/n)) } c++ }
    { print > f } END { printf "%d\n", int((c + n - 1)/n) }' "$FASTA" > "$OUTDIR/chunks.total"
  N_CH=$(cat "$OUTDIR"/chunks/*.fa | grep -c '^>')
  [ "$N_CH" = "$N_IN" ] || { echo "FATAL: split holds $N_CH of $N_IN" >&2; exit 1; }
  say "  $(cat "$OUTDIR/chunks.total") chunks -- lossless"
fi
TOTAL=$(cat "$OUTDIR/chunks.total")

run_one() {
  local f="$1" b; b=$(basename "$f" .fa)
  [ -f "$OUTDIR/done/$b.done" ] && return 0
  # runsanspanz writes relative to CWD, so each chunk gets its own directory.
  local w="$OUTDIR/logs/$b.work"; mkdir -p "$w"
  if ( cd "$PZDIR" && "$PY" runsanspanz.py -R -o ",,$OUTDIR/go/$b.GO.out," \
         -s "$SPECIES" < "$f" > "$w/run.log" 2>&1 ); then
    [ -s "$OUTDIR/go/$b.GO.out" ] && touch "$OUTDIR/done/$b.done"
  else
    echo "CHUNK FAILED: $b" >&2
  fi
}
export -f run_one; export OUTDIR PZDIR PY SPECIES

say "running $PAR concurrent streams (done: $(ls "$OUTDIR"/done/*.done 2>/dev/null | wc -l)/$TOTAL)"
ls "$OUTDIR"/chunks/*.fa | xargs -P "$PAR" -I{} bash -c 'run_one "$@"' _ {}

DONE=$(ls "$OUTDIR"/done/*.done 2>/dev/null | wc -l)
say "chunks finished: $DONE / $TOTAL"
[ "$DONE" = "$TOTAL" ] || { say "NOT merging -- re-run to resume"; exit 1; }

MERGED="$OUTDIR/${STUDY}.pannzer_GO.tsv"
head -1 "$(ls "$OUTDIR"/go/*.GO.out | head -1)" > "$MERGED"
for f in "$OUTDIR"/go/*.GO.out; do tail -n +2 "$f"; done >> "$MERGED"
say "merged: $(( $(wc -l < "$MERGED") - 1 )) predictions over $(tail -n +2 "$MERGED" | cut -f1 | sort -u | wc -l) proteins"
say "done: $STUDY"
