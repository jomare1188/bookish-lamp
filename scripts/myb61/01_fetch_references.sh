#!/usr/bin/env bash
# ============================================================================
# 01_fetch_references.sh — download the reference proteomes used as anchors
#
# Arabidopsis (TAIR10)  : source of the MYB61 anchor AT1G09540, and the target
#                         of the reciprocal-best-hit test in step 04.
# Sorghum, rice         : monocot anchors, so the search does not have to cross
#                         ~150 My of divergence in a single jump. Sorghum is the
#                         closest well-annotated relative of Saccharum (same
#                         tribe, Andropogoneae).
#
# Each proteome is reduced to ONE protein per gene (Ensembl pep headers carry
# `gene:<ID>`; we keep the longest isoform) so that "best hit" is a statement
# about genes, not about isoform counts.
#
# RUN: ./01_fetch_references.sh
# ============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

# longest isoform per Ensembl `gene:` tag, renamed to the gene id
primary_per_gene() {
  local gz="$1" out="$2"
  zcat "$gz" | awk '
    /^>/ {
      if (id != "") { seq[g] = (length(s) > length(seq[g])) ? s : seq[g] }
      s = ""; id = substr($1, 2); g = ""
      for (i = 1; i <= NF; i++) if ($i ~ /^gene:/) { g = substr($i, 6); sub(/\.[0-9]+$/, "", g) }
      if (g == "") g = id
      next
    }
    { s = s $0 }
    END {
      if (id != "") { seq[g] = (length(s) > length(seq[g])) ? s : seq[g] }
      for (k in seq) print ">" k "\n" seq[k]
    }' > "$out"
}

fetch() {
  local url="$1" out="$2" tag="$3"
  local gz="${REFDIR}/$(basename "$url")"
  if [ -s "$out" ]; then
    echo "[skip] $tag already present: $out ($(grep -c '^>' "$out") proteins)"
    return
  fi
  [ -s "$gz" ] || { echo "[get ] $tag <- $url"; curl -fsSL -o "${gz}.tmp" "$url" && mv "${gz}.tmp" "$gz"; }
  echo "[prep] $tag -> one protein per gene"
  primary_per_gene "$gz" "${out}.tmp" && mv "${out}.tmp" "$out"
  echo "[ok  ] $tag: $(grep -c '^>' "$out") genes"
}

fetch "$AT_URL" "$AT_PEP" "Arabidopsis thaliana TAIR10"
fetch "$SB_URL" "$SB_PEP" "Sorghum bicolor NCBIv3"
fetch "$OS_URL" "$OS_PEP" "Oryza sativa IRGSP-1.0"

# sanity: the anchor must exist in the Arabidopsis proteome
if grep -q "^>${AT_ANCHOR}$" "$AT_PEP"; then
  echo "[ok  ] anchor ${AT_ANCHOR} found in ${AT_PEP##*/}"
else
  echo "ERROR: anchor ${AT_ANCHOR} not found in $AT_PEP" >&2; exit 1
fi
