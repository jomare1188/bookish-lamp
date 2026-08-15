#!/usr/bin/env bash
# ============================================================================
# 03_search_proteomes.sh — forward homology search in OUR two proteomes
#
# Query : the At/Sb/Os MYB61 anchors from step 02
# Target: R570 (sugarcane) and LA purple (purple) — the same proteomes used for
#         the OrthoFinder bridge and the GET_TFS TF calls, so ids line up.
#
# Deliberately permissive (E<=1e-10, up to 500 targets): the aim is to capture
# the whole MYB61 neighbourhood including the polyploid haplotype copies. The
# strict disambiguation happens in step 04 (reciprocal best hit).
# Only the query-coverage filter is applied here, to drop fragmentary matches
# that align to a single Myb repeat.
#
# NOTE: the LA purple proteome carries 542 sequence lines with '.' characters
# (annotation gap/ambiguity marks). HMMER tolerates them, DIAMOND does not, so
# each proteome is sanitised into WORKDIR first: '*' dropped, any non-standard
# residue character replaced by 'X'. Ids and lengths are unchanged.
#
# RUN: ./03_search_proteomes.sh
# ============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

BL6='qseqid sseqid pident length evalue bitscore qcovhsp scovhsp'

# make a DIAMOND-safe copy of a proteome (ids untouched)
sanitise() {
  local in="$1" out="$2"
  awk '/^>/ { print $1; next }
       { s = toupper($0); gsub(/\*/, "", s); gsub(/[^ACDEFGHIKLMNPQRSTVWYBZXUO]/, "X", s); print s }' \
      "$in" > "${out}.tmp" && mv "${out}.tmp" "$out"
}

search() {
  local sp="$1" pep="$2" strip_iso="$3"
  local clean="${WORKDIR}/${sp}.clean.faa"
  local db="${WORKDIR}/${sp}" raw="${WORKDIR}/${sp}_forward_raw.tsv"
  local hits="${OUTDIR}/${sp}_forward_hits.tsv"

  [ -s "$clean" ] || { echo "[prep] ${sp}: sanitising proteome"; sanitise "$pep" "$clean"; }
  [ -s "${db}.dmnd" ] || { echo "[db  ] ${sp} ($(grep -c '^>' "$clean") proteins)"; \
                           diamond makedb --quiet --in "$clean" -d "$db"; }

  echo "[srch] ${sp}: MYB61 anchors -> proteome"
  diamond blastp --quiet -q "$QUERY" -d "$db" -o "$raw" --outfmt 6 $BL6 \
      --evalue "$FWD_EVALUE" --max-target-seqs "$FWD_MAX_TARGETS" --threads "$THREADS"

  # best HSP per (anchor, protein) above the query-coverage floor, then collapse
  # proteins to genes keeping each gene's best-scoring anchor hit
  awk -v OFS='\t' -v mincov="$MIN_QCOV" -v strip="$strip_iso" '
    $7 >= mincov {
      gene = $2
      if (strip == 1) sub(/[.][0-9]+[.]p[0-9]*$/, "", gene)
      if (!(gene in bs) || $6 > bs[gene]) {
        bs[gene] = $6; prot[gene] = $2; anc[gene] = $1
        pid[gene] = $3; ev[gene] = $5; qc[gene] = $7; sc[gene] = $8
      }
    }
    END { for (g in bs) print g, prot[g], anc[g], pid[g], ev[g], bs[g], qc[g], sc[g] }' \
    "$raw" | sort -k6,6gr > "${hits}.body"

  { printf 'gene\tprotein\tbest_anchor\tpident\tevalue\tbitscore\tqcov\tscov\n'; cat "${hits}.body"; } > "$hits"
  rm -f "${hits}.body"
  echo "[ok  ] ${sp}: $(($(wc -l < "$hits") - 1)) genes with a MYB61-anchor hit (qcov>=${MIN_QCOV}%) -> ${hits##*/}"
  awk 'NR>1 && NR<=6 {printf "       %-34s anchor=%-22s bits=%-7s qcov=%s%%\n", $1, $3, $6, $7}' "$hits"
}

search sugarcane "$SC_PEP" 1
search purple    "$PU_PEP" 0
