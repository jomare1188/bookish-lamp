#!/usr/bin/env bash
# ============================================================================
# 05_orthogroup_check.sh — place the MYB61 orthologues on the OrthoFinder bridge
#
# The homology search (steps 02-04) identifies MYB61 in each proteome
# INDEPENDENTLY. This step asks the question that actually licenses the
# comparative claim: do the R570 and LA purple copies land in the SAME
# orthogroup(s)? If they do, "MYB61 in sugarcane" and "MYB61 in purple" are the
# same evolutionary entity and the two networks can be compared at this locus.
#
# It also cross-checks the locus published by Kiet et al. (KIET_LOCUS): if any
# of its copies shared an orthogroup with a confirmed MYB61 copy, the published
# id would be partially reconcilable. Reported either way.
#
# CAVEATS this script makes explicit rather than hides:
#   * Orthogroups.tsv has CRLF line endings — fields are stripped of \r.
#   * OrthoFinder leaves many polyploid haplotype copies UNASSIGNED (singletons).
#     A MYB61 copy missing from Orthogroups.tsv is therefore not evidence of
#     absence; each copy is reported as in_orthogroup / unassigned / not_in_input.
#
# RUN: ./05_orthogroup_check.sh
# ============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

OGDIR="$(dirname "$ORTHOGROUPS")"
UNASSIGNED="${OGDIR}/Orthogroups_UnassignedGenes.tsv"
BRIDGE_SC="${BASE}/files/fix_orthofinder/proteins/sugarcane_one_transcript.fa"
BRIDGE_PU="${BASE}/files/fix_orthofinder/proteins/one_transcript_purple_proteins.faa"

OUT="${OUTDIR}/MYB61_orthogroups.tsv"
STATUS="${OUTDIR}/MYB61_bridge_status.tsv"

# --- 1. per-orthogroup view -------------------------------------------------
awk -F'\t' -v OFS='\t' \
    -v scids="${OUTDIR}/MYB61_orthologs_sugarcane.ids" \
    -v puids="${OUTDIR}/MYB61_orthologs_purple.ids" \
    -v locus="$KIET_LOCUS" '
  function clean(s) { gsub(/\r/, "", s); gsub(/^[ \t]+|[ \t]+$/, "", s); return s }
  BEGIN {
    while ((getline l < scids) > 0) sc[clean(l)] = 1
    while ((getline l < puids) > 0) pu[clean(l)] = 1
    print "orthogroup", "n_purple_total", "n_sugarcane_total",
          "MYB61_purple", "MYB61_sugarcane", "kiet_locus_copies"
  }
  NR == 1 { next }
  {
    np = ($2 == "") ? 0 : split($2, P, /, */)
    ns = ($3 == "") ? 0 : split($3, S, /, */)
    mp = ""; ms = ""; kl = ""
    for (i = 1; i <= np; i++) { p = clean(P[i])
      if (p in pu) mp = mp (mp ? "," : "") p
      if (index(p, locus) > 0) kl = kl (kl ? "," : "") p }
    for (i = 1; i <= ns; i++) { s = clean(S[i])
      if (s in sc) ms = ms (ms ? "," : "") s }
    if (mp != "" || ms != "" || kl != "")
      print $1, np, ns, (mp ? mp : "-"), (ms ? ms : "-"), (kl ? kl : "-")
  }' "$ORTHOGROUPS" > "$OUT"

# --- 2. per-gene bridge status ----------------------------------------------
{
  printf 'species\tgene\tbridge_status\torthogroup\n'
  for spec in "sugarcane:${BRIDGE_SC}:3" "purple:${BRIDGE_PU}:2"; do
    sp="${spec%%:*}"; rest="${spec#*:}"; inp="${rest%:*}"; col="${rest##*:}"
    while read -r g; do
      og=$(awk -F'\t' -v g="$g" -v c="$col" '
             { f = $c; gsub(/\r/, "", f)
               n = split(f, A, /, */)
               for (i = 1; i <= n; i++) { x = A[i]; gsub(/^[ \t]+|[ \t]+$/, "", x)
                 if (x == g) { print $1; exit } } }' "$ORTHOGROUPS")
      if [ -n "$og" ]; then st="in_orthogroup"
      elif grep -qF "$g" "$UNASSIGNED" 2>/dev/null; then st="unassigned_singleton"; og="-"
      elif grep -q "^>${g}$" "$inp"; then st="in_input_not_reported"; og="-"
      else st="not_in_bridge_input"; og="-"; fi
      printf '%s\t%s\t%s\t%s\n' "$sp" "$g" "$st" "$og"
    done < "${OUTDIR}/MYB61_orthologs_${sp}.ids"
  done
} > "$STATUS"

# --- 3. report --------------------------------------------------------------
echo "[ok  ] orthogroups touching MYB61 or the published locus -> ${OUT##*/}"
echo
printf '%-14s %6s %6s  %s\n' ORTHOGROUP "#pu" "#sc" "content"
awk -F'\t' 'NR>1 {
  npu = ($4 == "-") ? 0 : gsub(/,/, ",", $4) + 1
  nsc = ($5 == "-") ? 0 : gsub(/,/, ",", $5) + 1
  nkl = ($6 == "-") ? 0 : gsub(/,/, ",", $6) + 1
  tag = (npu && nsc) ? "*** SHARED MYB61 ***" : (npu ? "purple-only" : (nsc ? "sugarcane-only" : ""))
  if (nkl) tag = tag (tag ? " + " : "") "kiet-locus x" nkl
  printf "%-14s %6s %6s  MYB61: %d pu / %d sc   [%s]\n", $1, $2, $3, npu, nsc, tag
}' "$OUT"

echo
echo "=== bridge status of every confirmed MYB61 copy ==="
awk -F'\t' 'NR>1 {n[$1"\t"$3]++} END {for (k in n) printf "    %-12s %-24s %d\n", substr(k,1,index(k,"\t")-1), substr(k,index(k,"\t")+1), n[k]}' "$STATUS" | sort
echo "    (full detail -> ${STATUS##*/})"

echo
echo "=== does the published locus ever share an orthogroup with a confirmed MYB61 copy? ==="
if awk -F'\t' 'NR>1 && $6 != "-" && ($4 != "-" || $5 != "-") {found=1} END{exit !found}' "$OUT"; then
  echo "    YES — the published id is partially reconcilable:"
  awk -F'\t' 'NR>1 && $6 != "-" && ($4 != "-" || $5 != "-") {print "    " $0}' "$OUT"
else
  echo "    NO — ${KIET_LOCUS} copies never co-occur with a confirmed MYB61 copy."
  echo "    The published id cannot be mapped onto MYB61 in this annotation."
fi
