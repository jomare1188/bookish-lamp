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
mkdir -p "$OUTDIR"/{chunks,go,go_partial,done,logs}

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

# A FAILED CHUNK LEAVES A PARTIAL FILE, AND IT MUST NOT LOOK LIKE OUTPUT. Measured
# on the 2026-09-12 sugarcane run: 21 of 195 chunks died on a ConnectTimeout to
# Helsinki, and every one had already written thousands of real lines before the
# connection dropped -- c_00008 held 14,017 lines, c_00073 1.5 MB of half-written
# records. `go/` therefore held 195 files of which only 174 were complete, with
# nothing in the filename to tell them apart. The merge below is gated on the done
# count so it could not have consumed them, but anything merging `go/*.GO.out` by
# hand would have, silently. Partials are moved to go_partial/ instead.
run_one() {
  local f="$1" b; b=$(basename "$f" .fa)
  [ -f "$OUTDIR/done/$b.done" ] && return 0
  # runsanspanz writes relative to CWD, so each chunk gets its own directory.
  local w="$OUTDIR/logs/$b.work"; mkdir -p "$w"
  if ( cd "$PZDIR" && "$PY" runsanspanz.py -R -o ",,$OUTDIR/go/$b.GO.out," \
         -s "$SPECIES" < "$f" > "$w/run.log" 2>&1 ) \
     && [ -s "$OUTDIR/go/$b.GO.out" ]; then
    touch "$OUTDIR/done/$b.done"
  else
    [ -f "$OUTDIR/go/$b.GO.out" ] &&
      mv -f "$OUTDIR/go/$b.GO.out" "$OUTDIR/go_partial/$b.GO.out.$(date +%s)"
    echo "CHUNK FAILED: $b" >&2
  fi
}
export -f run_one; export OUTDIR PZDIR PY SPECIES

# RETRY, BECAUSE THE FAILURE MODE IS TRANSIENT AND MEASURED. Every one of those 21
# failures was the same ConnectTimeout to a public academic service, and the server
# answered normally again hours later -- so a chunk that failed is not a chunk that
# cannot succeed. One pass with no retry turned a transient network blip into a
# 10.8%-of-the-proteome hole that needed a human to notice. Passes are bounded and
# back off, so a service that is genuinely down still ends the run rather than
# being hammered.
PASSES="${CLEAN_PZ_PASSES:-4}"
BACKOFF="${CLEAN_PZ_BACKOFF:-300}"
for pass in $(seq 1 "$PASSES"); do
  DONE=$(ls "$OUTDIR"/done/*.done 2>/dev/null | wc -l)
  [ "$DONE" = "$TOTAL" ] && break
  if [ "$pass" -gt 1 ]; then
    say "pass $pass/$PASSES: $(( TOTAL - DONE )) chunk(s) still missing, waiting ${BACKOFF}s"
    sleep "$BACKOFF"
  fi
  say "pass $pass/$PASSES, $PAR concurrent streams (done: $DONE/$TOTAL)"
  ls "$OUTDIR"/chunks/*.fa | xargs -P "$PAR" -I{} bash -c 'run_one "$@"' _ {}
done

DONE=$(ls "$OUTDIR"/done/*.done 2>/dev/null | wc -l)
say "chunks finished: $DONE / $TOTAL"
if [ "$DONE" != "$TOTAL" ]; then
  say "NOT merging -- $(( TOTAL - DONE )) chunk(s) failed every pass. Re-run to resume;"
  say "partial output from the failures is in go_partial/ and is NOT merged."
  exit 1
fi

MERGED="$OUTDIR/${STUDY}.pannzer_GO.tsv"
# Belt and braces: every file about to be merged must have a done marker.
for f in "$OUTDIR"/go/*.GO.out; do
  b=$(basename "$f" .GO.out)
  [ -f "$OUTDIR/done/$b.done" ] ||
    { echo "FATAL: $b has output but no done marker -- refusing to merge a partial" >&2; exit 1; }
done
head -1 "$(ls "$OUTDIR"/go/*.GO.out | head -1)" > "$MERGED"
for f in "$OUTDIR"/go/*.GO.out; do tail -n +2 "$f"; done >> "$MERGED"
say "merged: $(( $(wc -l < "$MERGED") - 1 )) predictions over $(tail -n +2 "$MERGED" | cut -f1 | sort -u | wc -l) proteins"
say "done: $STUDY"
