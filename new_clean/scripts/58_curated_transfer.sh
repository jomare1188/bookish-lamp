#!/usr/bin/env bash
# =============================================================================
# 58_curated_transfer.sh -- GO from CURATED plant annotation, by homology.
#
# WHY THIS TIER EXISTS. Everything else in this project's GO is domain-inferred:
# InterPro signatures and Pfam domains mapped through interpro2go/pfam2go. Those
# are IEA-equivalent -- electronically inferred, no curator, no experiment. This
# stage is the only route here to terms a human assigned after looking at data.
#
# THE SOURCE is Swiss-Prot Viridiplantae (42,106 reviewed proteins). It subsumes a
# separate TAIR download: reviewed Arabidopsis entries already carry the TAIR
# curation, and the flat file records the evidence code for every term --
#   DR   GO; GO:0005737; C:cytoplasm; IDA:TAIR.
# so experimental terms can be separated from the IEA ones that would otherwise
# smuggle the same domain-transfer evidence back in through a different door.
#
# EVIDENCE KEPT: EXP IDA IPI IMP IGI IEP and the high-throughput equivalents
# HTP HDA HMP HGI HEP. Everything else -- IEA above all, but also ISS/ISO/IBA and
# the other computational codes -- is dropped, because this tier is only worth
# having if it is genuinely different in kind from the domain tier.
#
# THE THRESHOLD IS REPORTED, NOT ASSERTED. Homology transfer is only as good as
# the homology, so coverage is printed across a grid of identity and coverage cuts
# and the chosen one is visible next to the alternatives.
#
# RUN: through run.sh  ->  ./run.sh curatedgo sugarcane
# =============================================================================
set -euo pipefail

STUDY="${CLEAN_STUDY:?}"
PROTEOME="${CLEAN_PROTEOME:?}"
REFDIR="${CLEAN_CURATED_DIR:?}"
OUT="${CLEAN_OUT_FILE:?}"
NODES="${CLEAN_NODE_METRICS:-}"
PID="${CLEAN_MIN_IDENT:-50}"
COV="${CLEAN_MIN_COV:-70}"
THREADS="${CLEAN_CORES:-200}"
TAXON="${CLEAN_CURATED_TAXON:-33090}"     # Viridiplantae

say() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }
mkdir -p "$REFDIR" "$(dirname "$OUT")"
FA="$REFDIR/sprot_${TAXON}.fasta"
TXT="$REFDIR/sprot_${TAXON}.txt"
MAP="$REFDIR/sprot_${TAXON}_curated_go.tsv"
DB="$REFDIR/sprot_${TAXON}.dmnd"

# --- reference: sequences and their curated GO -------------------------------
Q="reviewed:true+AND+taxonomy_id:${TAXON}"
if [ ! -s "$FA" ]; then
  say "fetching reviewed sequences (taxon $TAXON)"
  curl -fsS "https://rest.uniprot.org/uniprotkb/stream?query=${Q}&format=fasta" -o "$FA.part"
  mv "$FA.part" "$FA"
fi
say "reference: $(grep -c '^>' "$FA") sequences"

if [ ! -s "$TXT" ]; then
  say "fetching flat file (for evidence codes)"
  curl -fsS "https://rest.uniprot.org/uniprotkb/stream?query=${Q}&format=txt" -o "$TXT.part"
  mv "$TXT.part" "$TXT"
fi

if [ ! -s "$MAP" ]; then
  say "parsing curated GO"
  awk '
    /^AC   /{ if (!acc) { split($2,a,";"); acc=a[1] } ; next }
    /^DR   GO; GO:/{
      split($0,f,"; "); go=f[2]
      n=split(f[4],e,":"); code=e[1]; sub(/\.$/,"",code)
      if (code ~ /^(EXP|IDA|IPI|IMP|IGI|IEP|HTP|HDA|HMP|HGI|HEP)$/) print acc"\t"go"\t"code
      next }
    /^\/\//{ acc="" }
  ' "$TXT" | sort -u > "$MAP"
fi
say "curated GO: $(wc -l < "$MAP") pairs over $(cut -f1 "$MAP" | sort -u | wc -l) reference proteins"
say "  evidence: $(cut -f3 "$MAP" | sort | uniq -c | sort -rn | awk '{printf "%s=%s ",$2,$1}')"

# --- align ---------------------------------------------------------------------
[ -s "$DB" ] || { say "building diamond db"; diamond makedb --in "$FA" -d "${DB%.dmnd}" --quiet; }
HITS="$REFDIR/${STUDY}_vs_sprot.tsv"
if [ ! -s "$HITS" ]; then
  say "$STUDY: diamond blastp vs $(grep -c '^>' "$FA") curated proteins"
  diamond blastp -q "$PROTEOME" -d "$DB" -o "$HITS" -p "$THREADS" \
    --outfmt 6 qseqid sseqid pident length qcovhsp scovhsp evalue bitscore \
    --max-target-seqs 5 --evalue 1e-10 --quiet
fi
say "  $(wc -l < "$HITS") hits over $(cut -f1 "$HITS" | sort -u | wc -l) query proteins"

# --- how much survives each threshold, before one is chosen --------------------
say ""
say "coverage across thresholds (query proteins receiving >= 1 curated term):"
for p in 30 40 50 60 70; do
  row="   "
  for c in 50 70; do
    n=$(awk -F'\t' -v p="$p" -v c="$c" '
          NR==FNR { has[$1]=1; next }
          { split($2,a,"|"); acc=(length(a)>=2 ? a[2] : $2)
            if ($3>=p && $5>=c && $6>=c && (acc in has)) q[$1] }
          END { print length(q) }' "$MAP" "$HITS")
    row="$row  ident>=${p}% cov>=${c}%: $(printf "%'d" "$n")"
  done
  say "$row"
done

# --- transfer at the chosen cut -------------------------------------------------
say ""
say "transferring at ident >= ${PID}%, mutual coverage >= ${COV}%"
awk -F'\t' -v p="$PID" -v c="$COV" '
  NR==FNR { go[$1] = go[$1] "\t" $2; next }
  {
    split($2, a, "|"); acc = (length(a) >= 2 ? a[2] : $2)
    if ($3 >= p && $5 >= c && $6 >= c && (acc in go)) {
      n = split(go[acc], t, "\t")
      for (i = 2; i <= n; i++) print $1 "\t" t[i] "\tcurated"
    }
  }' "$MAP" "$HITS" | sort -u > "$OUT.body"
{ printf 'gene\tgo_id\tsource\n'; cat "$OUT.body"; } > "$OUT"
rm -f "$OUT.body"

say "wrote $(basename "$OUT"): $(( $(wc -l < "$OUT") - 1 )) pairs over $(awk -F'\t' 'NR>1{g[$1]}END{print length(g)}' "$OUT") genes"
if [ -n "$NODES" ] && [ -s "$NODES" ]; then
  awk -F'\t' 'NR>1{print $1}' "$NODES" | sort -u > /tmp/net_$$
  awk -F'\t' 'NR>1{print $1}' "$OUT" | sort -u > /tmp/ann_$$
  N=$(wc -l < /tmp/net_$$); G=$(comm -12 /tmp/net_$$ /tmp/ann_$$ | wc -l)
  say "network coverage: $G of $N genes ($(awk -v a="$G" -v b="$N" 'BEGIN{printf "%.1f",100*a/b}')%)"
  rm -f /tmp/net_$$ /tmp/ann_$$
fi
say "done: $STUDY"
