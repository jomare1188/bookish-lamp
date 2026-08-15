#!/usr/bin/env bash
# ============================================================================
# 04_rbh_filter.sh — disambiguate the forward hits down to true MYB61 orthologues
#
# Two independent filters are applied to the step-03 candidates:
#
#   (1) RECIPROCAL BEST HIT — each candidate protein is searched back against
#       the whole Arabidopsis proteome; it passes only if its best Arabidopsis
#       hit is AT1G09540 again. This is what separates MYB61 from the ~125
#       other Arabidopsis R2R3-MYBs, which a forward search alone cannot do.
#
#   (2) OUR OWN TF CALL — the candidate gene must be independently called a
#       MYB-family TF by the project's GET_TFS run (Pfam GA + PlnTFDB rules).
#       Filters (1) and (2) use different evidence (homology vs domain rules),
#       so agreement between them is meaningful; disagreement is reported
#       rather than hidden, because a candidate that is a clear MYB61 homologue
#       but carries no callable Myb domain indicates a broken gene model.
#
# Output: MYB61_candidates.tsv — every step-03 hit with both verdicts, and
#         MYB61_orthologs_<species>.ids — the genes passing BOTH filters.
#
# RUN: ./04_rbh_filter.sh
# ============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

BL6='qseqid sseqid pident length evalue bitscore qcovhsp scovhsp'
AT_DB="${WORKDIR}/$(basename "${AT_PEP%.fa}")"

filter() {
  local sp="$1" strip_iso="$2" tfcalls="$3"
  local hits="${OUTDIR}/${sp}_forward_hits.tsv"
  local clean="${WORKDIR}/${sp}.clean.faa"
  local cand="${WORKDIR}/${sp}_cand.faa" bwd="${WORKDIR}/${sp}_backward.tsv"
  local out="${OUTDIR}/${sp}_MYB61_candidates.tsv"

  # candidate proteins (one representative protein per candidate gene)
  cut -f2 "$hits" | tail -n +2 | sort -u > "${WORKDIR}/${sp}_cand.ids"
  awk 'NR==FNR{w[$1];next} /^>/{k=(substr($1,2) in w)} k' \
      "${WORKDIR}/${sp}_cand.ids" "$clean" > "$cand"

  echo "[rbh ] ${sp}: $(grep -c '^>' "$cand") candidate proteins -> Arabidopsis"
  diamond blastp --quiet -q "$cand" -d "$AT_DB" -o "$bwd" --outfmt 6 $BL6 \
      --evalue 1e-5 --max-target-seqs 5 --threads "$THREADS"

  # best Arabidopsis hit per candidate protein
  awk -v OFS='\t' '{ if (!($1 in bs) || $6 > bs[$1]) { bs[$1] = $6; best[$1] = $2 } }
                    END { for (q in best) print q, best[q], bs[q] }' "$bwd" \
      > "${WORKDIR}/${sp}_besthit_at.tsv"

  # our own MYB call, collapsed protein -> gene
  awk -v OFS='\t' -v strip="$strip_iso" '
      { g = $1; if (strip == 1) sub(/[.][0-9]+[.]p[0-9]*$/, "", g)
        if (fam[g] == "") fam[g] = $2; else if (index(fam[g], $2) == 0) fam[g] = fam[g] ";" $2 }
      END { for (g in fam) print g, fam[g] }' "$tfcalls" > "${WORKDIR}/${sp}_tf_by_gene.tsv"

  awk -v OFS='\t' -v anchor="$AT_ANCHOR" '
    FNR==NR && FILENAME == ARGV[1] { at[$1] = $2; atbs[$1] = $3; next }
    FILENAME == ARGV[2] { fam[$1] = $2; next }
    FNR == 1 { print $0, "best_At_hit", "At_bitscore", "RBH_is_MYB61", "our_TF_family", "is_MYB_call", "verdict"; next }
    {
      gene = $1; prot = $2
      a  = (prot in at) ? at[prot] : "none"
      ab = (prot in at) ? atbs[prot] : "NA"
      rbh = (a == anchor) ? "yes" : "no"
      f = (gene in fam) ? fam[gene] : "not_a_TF"
      myb = (f ~ /(^|;)MYB(;|$)/) ? "yes" : "no"
      v = (rbh == "yes" && myb == "yes") ? "MYB61_ortholog" \
          : (rbh == "yes" && myb == "no") ? "MYB61_homolog_no_MYB_call" \
          : (rbh == "no"  && myb == "yes") ? "other_MYB" : "rejected"
      print $0, a, ab, rbh, f, myb, v
    }' "${WORKDIR}/${sp}_besthit_at.tsv" "${WORKDIR}/${sp}_tf_by_gene.tsv" "$hits" > "$out"

  awk -F'\t' 'NR>1 && $14=="MYB61_ortholog" {print $1}' "$out" | sort -u \
      > "${OUTDIR}/MYB61_orthologs_${sp}.ids"

  echo "[ok  ] ${sp} verdicts:"
  awk -F'\t' 'NR>1 {n[$14]++} END {for (v in n) printf "       %-28s %d\n", v, n[v]}' "$out" | sort -k2,2nr
  echo "       -> $(wc -l < "${OUTDIR}/MYB61_orthologs_${sp}.ids") confirmed orthologue genes"
}

filter sugarcane 1 "$SC_TF"
filter purple    0 "$PU_TF"
