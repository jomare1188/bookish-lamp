#!/usr/bin/env bash
# =============================================================================
# 54_build_gene2go.sh -- derive GO from the nf-core/proteinannotator output.
#
# WHY THIS EXISTS. The GO annotation the pipeline had been using is crippled, and
# the cause is in emapper's source rather than in the data: emapper.py:405 sets
# --go_evidence default 'non-electronic', and :615 maps that to
# go_excluded = {"ND","IEA"}. The run passed no --go_evidence, so every
# electronically-inferred GO term was dropped. In that file every protein has an
# orthologous group (100%) and 89.2% have a Description, but only 7.7% have any
# GO. For grass proteins nearly all GO is IEA, so that 7.7% IS the exclusion.
#
# The consequence was a module GO analysis that could not run: the median
# responsive module had ZERO annotated genes, and 493 of 588 sugarcane and 165 of
# 182 purple responsive modules were gated out entirely.
#
# WHAT THE NEW ANNOTATION DOES AND DOES NOT GIVE. It does NOT contain GO either --
# its InterProScan TSV has 13 columns, and GO/Pathways are columns 14-15, emitted
# only under --goterms. What it gives is INTERPRO ACCESSIONS (col 12) and, from a
# separate HMMER run, PFAM DOMAINS. GO is derivable from both through the Gene
# Ontology Consortium's own mapping files, and that derivation is this script.
#
# THE TWO SOURCES ARE COMPLEMENTARY, not redundant: the InterProScan run covered
# PANTHER, TIGRFAM, PIRSF, Hamap and SFLD but NOT Pfam. Measured on network genes:
#
#              via InterPro   via Pfam   union     (was, from emapper)
#   sugarcane        33,313     46,453   56,974     8,173
#   purple           50,768     77,801   94,497    12,223
#
# THE PFAM E-VALUE THRESHOLD IS THE ONE REAL PARAMETER. pfam.resolved is NOT
# pre-filtered -- it contains hits at i_evalue 14. The default 1e-5 keeps 71% of
# rows; coverage at 1e-3 and 1e-2 is reported too so the choice can be audited
# rather than taken on trust.
#
# NO GO-DAG PROPAGATION is done here. topGO walks the ontology itself, and
# pre-propagating ancestors would double-count them.
#
# THE MAPPINGS ARE VENDORED, with their own version line recorded. A GO mapping
# that changes silently between runs makes two analyses incomparable.
#
# RUN: through run.sh  ->  ./run.sh gene2go sugarcane
# =============================================================================
set -euo pipefail

STUDY="${CLEAN_STUDY:?}"
ANNOT="${CLEAN_ANNOT_DIR:?}"          # .../sugarcane_purple/merged
MAPDIR="${CLEAN_GO_MAP_DIR:?}"
OUT="${CLEAN_OUT_FILE:?}"
NODES="${CLEAN_NODE_METRICS:-}"
EVAL="${CLEAN_PFAM_EVALUE:-1e-5}"

IPS="${CLEAN_IPS_TSV:-${ANNOT}/${STUDY}.interproscan.tsv}"
PFAM="${ANNOT}/${STUDY}.pfam.resolved.tsv.gz"
say() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }

for f in "$IPS" "$PFAM"; do
  [ -s "$f" ] || { echo "FATAL: missing $f" >&2; exit 1; }
done
mkdir -p "$MAPDIR" "$(dirname "$OUT")"

# --- vendor the mappings, once, with their version stamps --------------------
fetch_map() {
  local name="$1" url="https://current.geneontology.org/ontology/external2go/$1"
  if [ ! -s "${MAPDIR}/${name}" ]; then
    say "fetching ${name}"
    curl -fsS -o "${MAPDIR}/${name}.part" "$url" || {
      echo "FATAL: could not fetch $url" >&2; rm -f "${MAPDIR}/${name}.part"; exit 1; }
    mv "${MAPDIR}/${name}.part" "${MAPDIR}/${name}"
    { echo "source: $url"
      echo "fetched: $(date -Is)"
      grep -m1 '^!version date' "${MAPDIR}/${name}" || true; } > "${MAPDIR}/${name}.provenance"
  fi
  say "  ${name}: $(grep -m1 '^!version date' "${MAPDIR}/${name}" | sed 's/^!//')"
}
fetch_map interpro2go
fetch_map pfam2go

# --- the two mappings, as flat key -> GO tables --------------------------------
# interpro2go:  InterPro:IPR000013 desc > GO:name ; GO:0004222
# pfam2go:      Pfam:PF00001 7tm_1 > GO:name ; GO:0004930   <- join on the NAME,
#               because pfam.resolved carries the Pfam name, not the accession.
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
sed -n 's/^InterPro:\([^ ]*\) .* ; \(GO:[0-9]*\)$/\1\t\2/p' "${MAPDIR}/interpro2go" | sort -u > "$TMP/ipr2go"
sed -n 's/^Pfam:[^ ]* \([^ ]*\) .* ; \(GO:[0-9]*\)$/\1\t\2/p'  "${MAPDIR}/pfam2go"    | sort -u > "$TMP/pf2go"
say "mappings: $(wc -l < "$TMP/ipr2go") IPR->GO pairs, $(wc -l < "$TMP/pf2go") Pfam->GO pairs"

# --- InterPro route -----------------------------------------------------------
say "$STUDY: InterPro -> GO"
awk -F'\t' 'NR==FNR{m[$1]=m[$1]"\t"$2; next}
     $12 != "-" && $12 != "" && ($12 in m){ n=split(m[$12],g,"\t")
       for(i=2;i<=n;i++) print $1"\t"g[i] }' "$TMP/ipr2go" "$IPS" | sort -u > "$TMP/from_ipr"
say "  $(cut -f1 "$TMP/from_ipr" | sort -u | wc -l) proteins, $(wc -l < "$TMP/from_ipr") gene-GO pairs"

# --- Pfam route ---------------------------------------------------------------
say "$STUDY: Pfam -> GO  (i_evalue <= $EVAL)"
zcat "$PFAM" | awk -F'\t' -v ev="$EVAL" 'NR==FNR{m[$1]=m[$1]"\t"$2; next}
     FNR>1 && ($6+0)<=(ev+0) && ($3 in m){ n=split(m[$3],g,"\t")
       for(i=2;i<=n;i++) print $1"\t"g[i] }' "$TMP/pf2go" - | sort -u > "$TMP/from_pfam"
say "  $(cut -f1 "$TMP/from_pfam" | sort -u | wc -l) proteins, $(wc -l < "$TMP/from_pfam") gene-GO pairs"

# --- union, with the source kept per pair -------------------------------------
# Which source a pair came from has to stay visible: if a downstream result rests
# entirely on one of them, that is something to know, not to discover later.
awk -F'\t' 'NR==FNR{a[$0]; next}{b[$0]}
     END{ for(k in a) print k"\t" (k in b ? "both" : "ipr")
          for(k in b) if(!(k in a)) print k"\tpfam" }' "$TMP/from_ipr" "$TMP/from_pfam" \
  | sort -k1,1 -k2,2 > "$TMP/all"
{ printf 'gene\tgo_id\tsource\n'; cat "$TMP/all"; } > "$OUT"

say ""
say "wrote $(basename "$OUT"): $(wc -l < "$TMP/all") gene-GO pairs over $(cut -f1 "$TMP/all" | sort -u | wc -l) genes"
awk -F'\t' 'NR>1{s[$3]++} END{for(k in s) printf "  %-6s %d pairs\n", k, s[k]}' "$OUT"

# --- what this means for THIS network ------------------------------------------
if [ -n "$NODES" ] && [ -s "$NODES" ]; then
  awk -F'\t' 'NR>1{print $1}' "$NODES" | sort -u > "$TMP/net"
  cut -f1 "$TMP/all" | sort -u > "$TMP/annot"
  N=$(wc -l < "$TMP/net"); G=$(comm -12 "$TMP/net" "$TMP/annot" | wc -l)
  say ""
  say "network coverage: $G of $N genes carry >= 1 GO term ($(awk -v a="$G" -v b="$N" 'BEGIN{printf "%.1f", 100*a/b}')%)"
  for e in 1e-3 1e-2; do
    zcat "$PFAM" | awk -F'\t' -v ev="$e" 'NR==FNR{m[$1];next} FNR>1 && ($6+0)<=(ev+0) && ($3 in m){print $1}' "$TMP/pf2go" - \
      | sort -u | cat - "$TMP/from_ipr" | cut -f1 | sort -u > "$TMP/alt"
    say "  (if the Pfam cut were $e: $(comm -12 "$TMP/net" "$TMP/alt" | wc -l) genes)"
  done
fi
say "done: $STUDY"
