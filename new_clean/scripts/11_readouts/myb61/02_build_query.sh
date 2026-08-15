#!/usr/bin/env bash
# ============================================================================
# 02_build_query.sh — build the MYB61 anchor query set
#
# Starts from AtMYB61 (AT1G09540) and finds its *reciprocal best hit* in
# sorghum and rice. RBH is used rather than a gene-name lookup because
# "MYB61" is not consistently annotated in monocot references, and because
# R2R3-MYB is one of the largest plant TF families (~125 genes in
# Arabidopsis) — a plain top-hit would give "some MYB", not MYB61.
#
#   AT1G09540 --forward--> Sb / Os candidates --backward--> Arabidopsis
#   keep a candidate only if its best Arabidopsis hit is AT1G09540 again.
#
# The surviving At + Sb + Os proteins become the query set for step 03, and
# their Myb repeat structure is verified with the project's own HMM library.
#
# RUN: ./02_build_query.sh
# ============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

BL6='qseqid sseqid pident length evalue bitscore qcovhsp scovhsp'

# one-time diamond databases for the reference proteomes
for ref in "$AT_PEP" "$SB_PEP" "$OS_PEP"; do
  db="${WORKDIR}/$(basename "${ref%.fa}")"
  [ -s "${db}.dmnd" ] || { echo "[db  ] $(basename "$ref")"; diamond makedb --quiet --in "$ref" -d "$db"; }
done
AT_DB="${WORKDIR}/$(basename "${AT_PEP%.fa}")"

# the anchor protein
awk -v id="$AT_ANCHOR" '/^>/{k=($0==">"id)} k' "$AT_PEP" > "${WORKDIR}/anchor.faa"
echo "[ok  ] anchor ${AT_ANCHOR}: $(awk '!/^>/{n+=length($0)}END{print n}' "${WORKDIR}/anchor.faa") aa"

# --- reciprocal best hit of the anchor in one target proteome ---------------
rbh() {
  local tag="$1" ref="$2"
  local db="${WORKDIR}/$(basename "${ref%.fa}")"
  local fwd="${WORKDIR}/${tag}_forward.tsv" bwd="${WORKDIR}/${tag}_backward.tsv"

  # forward: anchor -> target (keep a shortlist, not just the single top hit)
  diamond blastp --quiet -q "${WORKDIR}/anchor.faa" -d "$db" -o "$fwd" \
      --outfmt 6 $BL6 --evalue 1e-5 --max-target-seqs 20 --threads "$THREADS"

  # backward: each shortlisted candidate -> Arabidopsis, best hit only
  cut -f2 "$fwd" | sort -u > "${WORKDIR}/${tag}_cand.ids"
  awk 'NR==FNR{w[$1];next} /^>/{k=(substr($1,2) in w)} k' \
      "${WORKDIR}/${tag}_cand.ids" "$ref" > "${WORKDIR}/${tag}_cand.faa"
  diamond blastp --quiet -q "${WORKDIR}/${tag}_cand.faa" -d "$AT_DB" -o "$bwd" \
      --outfmt 6 $BL6 --evalue 1e-5 --max-target-seqs 1 --threads "$THREADS"

  # a candidate passes when its best Arabidopsis hit is the anchor itself
  awk -v anchor="$AT_ANCHOR" '
    { if (!($1 in best) || $6 > bs[$1]) { best[$1] = $2; bs[$1] = $6 } }
    END { for (q in best) if (best[q] == anchor) print q }' "$bwd" | sort -u
}

: > "${WORKDIR}/anchors_monocot.ids"
for spec in "Sb:${SB_PEP}" "Os:${OS_PEP}"; do
  tag="${spec%%:*}"; ref="${spec#*:}"
  echo "[rbh ] ${tag} <-> ${AT_ANCHOR}"
  hits=$(rbh "$tag" "$ref")
  if [ -z "$hits" ]; then
    echo "[warn] no reciprocal best hit for ${tag}" >&2
  else
    echo "$hits" | sed "s/^/${tag}\t/" >> "${WORKDIR}/anchors_monocot.ids"
    echo "$hits" | sed 's/^/       RBH: /'
  fi
done

# --- assemble the query set, prefixing ids with the source species ----------
{
  sed "s/^>/>At|/" "${WORKDIR}/anchor.faa"
  while IFS=$'\t' read -r tag id; do
    ref=$([ "$tag" = "Sb" ] && echo "$SB_PEP" || echo "$OS_PEP")
    awk -v id="$id" -v tag="$tag" '/^>/{k=($0==">"id); if(k) print ">" tag "|" id; next} k' "$ref"
  done < "${WORKDIR}/anchors_monocot.ids"
} > "$QUERY"
echo "[ok  ] query set: $(grep -c '^>' "$QUERY") proteins -> ${QUERY}"
grep '^>' "$QUERY" | sed 's/^/       /'

# --- verify the queries really are R2R3-MYBs (our own HMM library) ---------
HMMDB="${BASE}/GET_TFS/db/TF.db.hmm"
if command -v hmmfetch >/dev/null && command -v hmmsearch >/dev/null; then
  hmmfetch "$HMMDB" Myb_DNA-binding > "${WORKDIR}/Myb.hmm" 2>/dev/null || true
  if [ -s "${WORKDIR}/Myb.hmm" ]; then
    hmmsearch --cut_ga --domtblout "${WORKDIR}/query_myb.domtbl" -o /dev/null \
        "${WORKDIR}/Myb.hmm" "$QUERY"
    echo "[chk ] Myb_DNA-binding repeats per query (--cut_ga):"
    awk '!/^#/{n[$1]++} END{for (q in n) printf "       %-28s %d repeats\n", q, n[q]}' \
        "${WORKDIR}/query_myb.domtbl"
  fi
fi
