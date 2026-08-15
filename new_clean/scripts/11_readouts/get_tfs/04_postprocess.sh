#!/bin/bash
# ============================================================================
# 04_postprocess.sh — turn per-protein family calls into a gene-level TF set
# and intersect it with the co-expression network nodes.
#
# Produces (in $OUTDIR):
#   TF_no_Orphans.tsv     protein-level calls, Orphans removed
#   TF_genes.ids          unique TF gene ids (all proteome genes with a TF call)
#   network_genes.ids     unique gene ids present in the network
#   TF_in_network.ids     TF genes that ARE nodes of the network   <- H1 input
#   TF_in_network.tsv     gene <tab> Family <tab> Type  (network TFs only)
#   family_counts.tsv     per-family counts (all TF genes vs in-network)
#   summary.txt           headline numbers
# Usage:  ./04_postprocess.sh <species>
# ============================================================================
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${here}/config.sh" "${1:?usage: 04_postprocess.sh <species>}"

# written by 03_assign_families.sh, which is sequence-only and cached
fam="${CACHE_DIR}/family_assignment.tsv"
[ -s "$fam" ] || { echo "ERROR: ${fam} missing — run 03_assign_families.sh, or point CACHE_DIR at an existing run"; exit 1; }

# 1) drop header + Orphans (keep real TAPs: TFF / OTR)
awk -F'\t' '$1!="Gene" && $2!="" && $2!="Orphans"' "$fam" > "${OUTDIR}/TF_no_Orphans.tsv"

# 2) protein id -> gene id (species-specific), then map to family/type at gene level
paste <(cut -f1 "${OUTDIR}/TF_no_Orphans.tsv" | sed -E "$PROT_TO_GENE") \
      <(cut -f2,3 "${OUTDIR}/TF_no_Orphans.tsv") \
  | sort -u > "${OUTDIR}/TF_genes_families.tsv"

cut -f1 "${OUTDIR}/TF_genes_families.tsv" | sort -u > "${OUTDIR}/TF_genes.ids"

# 3) network node ids -> gene ids
tail -n +2 "$NODE_METRICS" | cut -f1 | sed -E "$NODE_TO_GENE" | sort -u > "${OUTDIR}/network_genes.ids"

# 4) intersection: TFs that are nodes of the co-expression network
comm -12 "${OUTDIR}/TF_genes.ids" "${OUTDIR}/network_genes.ids" > "${OUTDIR}/TF_in_network.ids"

# 5) family table restricted to network TFs
{ printf "gene\tFamily\tType\n"
  awk -F'\t' 'NR==FNR{net[$1]=1; next} ($1 in net)' \
      "${OUTDIR}/TF_in_network.ids" "${OUTDIR}/TF_genes_families.tsv"
} > "${OUTDIR}/TF_in_network.tsv"

# 6) per-family counts: all TF genes vs in-network TF genes (single awk pass)
awk -F'\t' '
  BEGIN{ OFS="\t" }
  # file1: gene<tab>Family<tab>Type  (all TF genes)
  FNR==NR { fam[$2]=$3; all[$2]++; next }
  # file2: TF_in_network.tsv  (gene<tab>Family<tab>Type, with header)
  FNR>1 { net[$2]++ }
  END{
    print "Family","Type","TF_genes","in_network"
    for (f in all) print f, fam[f], all[f], (f in net ? net[f] : 0)
  }' "${OUTDIR}/TF_genes_families.tsv" "${OUTDIR}/TF_in_network.tsv" \
  | { read -r h; echo "$h"; sort -t$'\t' -k4,4nr -k3,3nr; } > "${OUTDIR}/family_counts.tsv"

# 7) summary
{
  echo "species            : ${SPECIES}"
  echo "proteome           : ${PROTEOME}"
  echo "proteins with TAP  : $(wc -l < "${OUTDIR}/TF_no_Orphans.tsv")"
  echo "TF genes (total)   : $(wc -l < "${OUTDIR}/TF_genes.ids")"
  echo "network genes      : $(wc -l < "${OUTDIR}/network_genes.ids")"
  echo "TF genes in network: $(wc -l < "${OUTDIR}/TF_in_network.ids")"
  echo "TFF families        : $(awk -F'\t' 'NR>1 && $3=="TFF"' "${OUTDIR}/TF_in_network.tsv" | cut -f2 | sort -u | wc -l)"
} | tee "${OUTDIR}/summary.txt"

echo "[post] outputs in ${OUTDIR}"
