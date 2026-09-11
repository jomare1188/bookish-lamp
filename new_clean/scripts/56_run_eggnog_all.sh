#!/usr/bin/env bash
# =============================================================================
# 56_run_eggnog_all.sh -- re-annotate eggNOG with --go_evidence all.
#
# THE DEFECT THIS FIXES. The existing emapper.annotations has GO for 7.7% of
# proteins while 89.2% have a Description and 100% have an orthologous group.
# emapper.py:405 defaults --go_evidence to 'non-electronic', and :615 maps that to
# go_excluded = {"ND","IEA"}. The original run passed no --go_evidence, so every
# electronically-inferred term was dropped -- and for grass proteins almost all GO
# is IEA. Description at 89.2% is the realistic ceiling for what GO can reach.
#
# WHY THIS IS CHEAP. The diamond search is the expensive step and it was saved as
# emapper.seed_orthologs. --annotate_hits_table re-runs ONLY the annotation
# transfer from those hits, so this is minutes rather than hours, and the hits are
# byte-identical to the ones the original run used -- the only thing that differs
# between old and new output is the GO evidence filter.
#
# NOTHING IS OVERWRITTEN. Output goes to a new directory. The original
# emapper.annotations stays as the record of the non-electronic run, and is still
# the source for KEGG and Preferred_name.
#
# RUN: through run.sh  ->  ./run.sh eggnogall sugarcane
# =============================================================================
set -euo pipefail

STUDY="${CLEAN_STUDY:?}"
SEED="${CLEAN_SEED_ORTHOLOGS:?}"
DATA="${CLEAN_EGGNOG_DATA:?}"
OUTDIR="${CLEAN_OUT_DIR:?}"
TAXSCOPE="${CLEAN_TAX_SCOPE:-Poales}"
CPU="${CLEAN_CORES:-100}"
EMAPPER="${CLEAN_EMAPPER_BIN:-/home/genomics/miniconda3/envs/eggnog_env/bin/emapper.py}"
ORIG="${CLEAN_ORIG_ANNOT:-}"

say() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }
[ -s "$SEED" ] || { echo "FATAL: no seed orthologs at $SEED" >&2; exit 1; }
[ -d "$DATA" ] || { echo "FATAL: no eggnog data dir at $DATA" >&2; exit 1; }
mkdir -p "$OUTDIR"

say "$STUDY: re-annotating $(wc -l < "$SEED") saved hits with --go_evidence all"
say "  data_dir  $DATA"
say "  tax_scope $TAXSCOPE"

"$EMAPPER" --annotate_hits_table "$SEED" \
  -o "$STUDY" --output_dir "$OUTDIR" \
  --data_dir "$DATA" --tax_scope "$TAXSCOPE" \
  --go_evidence all --cpu "$CPU" --override

NEW="${OUTDIR}/${STUDY}.emapper.annotations"
[ -s "$NEW" ] || { echo "FATAL: emapper wrote no annotations to $NEW" >&2; exit 1; }

fill() {  # file, column index -> "count (pct)"
  awk -F'\t' -v c="$2" '!/^#/{n++; if($c!="-" && $c!="") k++}
    END{printf "%d (%.1f%%)", k+0, n?100*k/n:0}' "$1"
}
rows() { awk -F'\t' '!/^#/{n++} END{print n+0}' "$1"; }

say ""
say "$STUDY: new annotation"
say "  rows        $(rows "$NEW")"
say "  eggNOG_OGs  $(fill "$NEW" 5)"
say "  Description $(fill "$NEW" 8)"
say "  GOs         $(fill "$NEW" 10)"
say "  PFAMs       $(fill "$NEW" 21)"

# --- the check that ONLY the GO filter moved ---------------------------------
# Same hits, same tax scope, one flag different: every other column must be
# identical. If Description or PFAMs shifted, something other than --go_evidence
# changed and the comparison with the old run is not clean.
if [ -n "$ORIG" ] && [ -s "$ORIG" ]; then
  say ""
  say "against the original (non-electronic) run:"
  say "  rows        $(rows "$ORIG")  ->  $(rows "$NEW")"
  say "  eggNOG_OGs  $(fill "$ORIG" 5)  ->  $(fill "$NEW" 5)"
  say "  Description $(fill "$ORIG" 8)  ->  $(fill "$NEW" 8)"
  say "  GOs         $(fill "$ORIG" 10)  ->  $(fill "$NEW" 10)"
  say "  PFAMs       $(fill "$ORIG" 21)  ->  $(fill "$NEW" 21)"
  for pair in 5:eggNOG_OGs 8:Description 21:PFAMs; do
    c="${pair%%:*}"; nm="${pair##*:}"
    a=$(fill "$ORIG" "$c"); b=$(fill "$NEW" "$c")
    [ "$a" = "$b" ] || { echo "FATAL: $nm changed ($a -> $b); more than --go_evidence differs" >&2; exit 1; }
  done
  [ "$(rows "$ORIG")" = "$(rows "$NEW")" ] || { echo "FATAL: row count changed" >&2; exit 1; }
  say "  -> only GOs moved, as intended"
fi
touch "${OUTDIR}/.complete"
say "done: $STUDY"
