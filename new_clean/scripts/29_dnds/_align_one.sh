#!/usr/bin/env bash
# worker: one orthogroup -> protein alignment -> codon alignment
# args: <orthogroup> <dir>
set -uo pipefail
OG="$1"; D="$2"
FAA="${D}/${OG}.faa"; FNA="${D}/${OG}.fna"
ALN="${D}/${OG}.aln.faa"; PAML="${D}/${OG}.paml"

# L-INS-i: 3 short sequences, so the accurate mode costs nothing and matches
# the convention already set in 11_readouts/myb61/06_phylogeny.sh
"$MAFFT" --quiet --localpair --maxiterate 1000 --amino "$FAA" > "$ALN" 2>/dev/null || {
  echo -e "${OG}\tmafft_failed"; exit 0; }

# -nogap drops any column with a gap or an in-frame stop, so the alignment
# codeml sees is gap-free; codontable 1 is the standard code.
"$PAL2NAL" "$ALN" "$FNA" -output paml -nogap -codontable 1 > "$PAML" 2>/dev/null || {
  echo -e "${OG}\tpal2nal_failed"; rm -f "$PAML"; exit 0; }

read -r _ NCOD < <(head -1 "$PAML")
NCOD=${NCOD:-0}
if [ "$NCOD" -lt $(( MIN_CODONS * 3 )) ]; then
  echo -e "${OG}\ttoo_short_${NCOD}bp"; rm -f "$PAML"; exit 0
fi
echo -e "${OG}\tok_${NCOD}"
