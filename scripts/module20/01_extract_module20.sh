#!/usr/bin/env bash
# ============================================================================
# 01_extract_module20.sh — pull the Module-20 PROTEINS out of the Muñoz proteome
#
# `module20.csv` ids carry a trailing `.gen`; the proteome holds TransDecoder
# ORFs named `<transcript>.pN` with `type:` and `len:` in the header. So the
# join is: module20 id  -minus .gen->  transcript  ->  all <transcript>.pN.
#
# Each member has 2-4 ORFs, some of them antisense/internal artifacts, and
# "longest" and "complete" disagree for at least one member
# (SCA2_2__c65732f1p01355: .p1 = 341 aa but 5prime_partial, .p2 = 181 aa
# complete). Rather than impose a rule here, ALL ORFs are kept and the
# reciprocal filter in step 02 decides which one represents the member; the
# winning ORF is recorded.
#
# The nucleotide contigs are still extracted, but only as a QC/legacy artefact.
#
# RUN: ./01_extract_module20.sh
# ============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

[ -s "$M20_CSV" ]       || { echo "ERROR: missing $M20_CSV" >&2; exit 1; }
[ -s "$M20_PROTEOME" ]  || { echo "ERROR: missing $M20_PROTEOME" >&2; exit 1; }

# --- the published member table --------------------------------------------
tail -n +2 "$M20_CSV" | tr -d '\r' | awk -F, -v OFS='\t' '
  { id = $1; sub(/\.gen$/, "", id); print id, $2, $3, $4, $5, $6 }' \
  > "${OUTDIR}/module20_members.tsv"
sed -i '1i transcript\tmodule\tfamily\tclass\tmunoz_degree\tmunoz_betweenness' \
    "${OUTDIR}/module20_members.tsv"

cut -f1 "${OUTDIR}/module20_members.tsv" | tail -n +2 | sort -u > "${WORKDIR}/m20.ids"
echo "[ok  ] $(wc -l < "${WORKDIR}/m20.ids") Module-20 members"
awk -F'\t' 'NR>1 {n[$3]++} END {for (k in n) printf "       %-14s %d\n", k, n[k]}' \
    "${OUTDIR}/module20_members.tsv"

# --- proteins ---------------------------------------------------------------
echo "[scan] extracting ORFs from $(basename "$M20_PROTEOME") ..."
awk 'NR==FNR { want[$1]; next }
     /^>/ { id = substr($1, 2); base = id; sub(/\.p[0-9]+$/, "", base)
            keep = (base in want)
            if (keep) print ">" id "\t" base "\t" $2 "\t" $3
            next }
     keep' "${WORKDIR}/m20.ids" "$M20_PROTEOME" > "${WORKDIR}/m20_prot_raw.txt"

# split the annotated header line back into fasta + an ORF table
awk -F'\t' '/^>/ { print $1; next } { print }' "${WORKDIR}/m20_prot_raw.txt" > "${M20_FAA}.tmp"
mv "${M20_FAA}.tmp" "$M20_FAA"
{ printf 'orf\ttranscript\torf_type\torf_len\n'
  awk -F'\t' '/^>/ { t = $3; l = $4; sub(/^type:/, "", t); sub(/^len:/, "", l)
                     printf "%s\t%s\t%s\t%s\n", substr($1,2), $2, t, l }' "${WORKDIR}/m20_prot_raw.txt"
} > "${OUTDIR}/module20_orfs.tsv"

n_orf=$(grep -c '^>' "$M20_FAA")
n_tr=$(awk -F'\t' 'NR>1 {t[$2]} END {print length(t)}' "${OUTDIR}/module20_orfs.tsv")
echo "[ok  ] ${n_orf} ORFs from ${n_tr} / $(wc -l < "${WORKDIR}/m20.ids") members -> ${M20_FAA##*/}"

# members with NO predicted protein — a result, not a warning: they are
# non-coding / UTR fragments, which is why they never mapped.
comm -23 "${WORKDIR}/m20.ids" <(awk -F'\t' 'NR>1 {print $2}' "${OUTDIR}/module20_orfs.tsv" | sort -u) \
  > "${OUTDIR}/module20_no_orf.ids"
if [ -s "${OUTDIR}/module20_no_orf.ids" ]; then
  echo "[note] $(wc -l < "${OUTDIR}/module20_no_orf.ids") member(s) have NO predicted ORF (non-coding / UTR fragment):"
  sed 's/^/       /' "${OUTDIR}/module20_no_orf.ids"
fi

awk -F'\t' 'NR>1 {printf "       %-34s %-16s %s aa\n", $1, $3, $4}' "${OUTDIR}/module20_orfs.tsv"

# --- nucleotide contigs, kept only for QC ----------------------------------
if [ ! -s "$M20_FNA" ]; then
  echo "[scan] (QC) extracting nucleotide contigs ..."
  awk 'NR==FNR { want[$1]; next }
       /^>/ { id = substr($1, 2); keep = (id in want); if (keep) print ">" id; next }
       keep' "${WORKDIR}/m20.ids" "$M20_TRANSCRIPTOME" > "$M20_FNA"
fi
echo "[ok  ] contigs (QC only): $(grep -c '^>' "$M20_FNA")"
