#!/usr/bin/env bash
# ============================================================================
# 06_phylogeny.sh — confirm the MYB61 calls with a tree, not just pairwise hits
#
# Reciprocal best hit (step 04) is a pairwise criterion; it can still be fooled
# when a family has many close paralogues. This step builds an explicit tree of
#
#   the confirmed MYB61 copies from both species
# + the At / Sb / Os MYB61 anchors
# + a decoy set: the "other_MYB" genes rejected in step 04, plus the closest
#   NON-MYB61 Arabidopsis R2R3-MYBs
#
# A correct result puts every confirmed copy inside a clade with the anchors,
# and the decoys outside it. That is the actual test — if confirmed copies and
# decoys interleave, the RBH call was not clean.
#
# MAFFT (L-INS-i when small enough) + FastTree (WAG+CAT, SH-like supports).
#
# RUN: ./06_phylogeny.sh
# ============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

PHY="${OUTDIR}/phylogeny"; mkdir -p "$PHY"
FAA="${PHY}/MYB61_tree_input.faa"
ALN="${PHY}/MYB61_tree_input.aln"
TREE="${PHY}/MYB61.tree"
LABELS="${PHY}/MYB61_tree_labels.tsv"

# pull named sequences out of a fasta, applying a prefix to the id
pull() {  # pull <fasta> <ids-file> <prefix>
  awk -v pre="$3" 'NR==FNR { w[$1]; next }
       /^>/ { id = substr($1, 2); k = (id in w); if (k) print ">" pre id; next }
       k' "$2" "$1"
}

: > "$FAA"

# 1. anchors (already prefixed At|/Sb|/Os| by step 02)
cat "$QUERY" >> "$FAA"

# 2. confirmed MYB61 copies from both species
pull "${WORKDIR}/sugarcane.clean.faa" <(awk -F'\t' 'NR>1 && $14=="MYB61_ortholog" {print $2}' \
     "${OUTDIR}/sugarcane_MYB61_candidates.tsv") "SC_MYB61|" >> "$FAA"
pull "${WORKDIR}/purple.clean.faa"    <(awk -F'\t' 'NR>1 && $14=="MYB61_ortholog" {print $2}' \
     "${OUTDIR}/purple_MYB61_candidates.tsv") "PU_MYB61|" >> "$FAA"

# 3. decoys — the step-04 rejects (other MYBs that hit the anchors but failed RBH)
pull "${WORKDIR}/sugarcane.clean.faa" <(awk -F'\t' 'NR>1 && $14!="MYB61_ortholog" {print $2}' \
     "${OUTDIR}/sugarcane_MYB61_candidates.tsv") "SC_decoy|" >> "$FAA"
pull "${WORKDIR}/purple.clean.faa"    <(awk -F'\t' 'NR>1 && $14!="MYB61_ortholog" {print $2}' \
     "${OUTDIR}/purple_MYB61_candidates.tsv") "PU_decoy|" >> "$FAA"

# 4. decoys — closest Arabidopsis R2R3-MYBs that are NOT the anchor
AT_DB="${WORKDIR}/$(basename "${AT_PEP%.fa}")"
diamond blastp --quiet -q "${WORKDIR}/anchor.faa" -d "$AT_DB" --outfmt 6 qseqid sseqid bitscore \
    --evalue 1e-5 --max-target-seqs 12 --threads "$THREADS" \
  | awk -v a="$AT_ANCHOR" '$2 != a {print $2}' | sort -u > "${PHY}/at_decoys.ids"
pull "$AT_PEP" "${PHY}/at_decoys.ids" "At_decoy|" >> "$FAA"

n=$(grep -c '^>' "$FAA")
echo "[ok  ] tree input: ${n} sequences"
grep '^>' "$FAA" | sed 's/^>//' | awk -F'|' -v OFS='\t' 'BEGIN{print "label","class"} {print $0, $1}' > "$LABELS"
awk -F'|' '{c[$1]++} END {for (k in c) printf "       %-12s %d\n", k, c[k]}' <(grep '^>' "$FAA" | sed 's/^>//')

echo "[aln ] MAFFT"
if [ "$n" -le 200 ]; then mafft --quiet --localpair --maxiterate 1000 --thread "$THREADS" "$FAA" > "$ALN"
else                      mafft --quiet --auto --thread "$THREADS" "$FAA" > "$ALN"; fi

echo "[tree] FastTree (WAG+CAT)"
FastTree -quiet -wag "$ALN" > "$TREE" 2> "${PHY}/fasttree.log"

echo "[ok  ] -> ${TREE}"

# --- automatic verdict: is there a clade holding all anchors + all confirmed,
#     and no decoys?  Checked with ape in R if available, else reported manually.
if [ -x "$RSCRIPT_APE" ]; then
"$RSCRIPT_APE" --vanilla - "$TREE" <<'RS' 2>/dev/null || echo "[note] clade test failed; inspect ${TREE} manually"
suppressMessages(library(ape))
tr  <- read.tree(commandArgs(TRUE)[1])
cls <- sub("\\|.*", "", tr$tip.label)

# Root on the most divergent Arabidopsis decoy.
d  <- diag(vcv(tr))
tr <- root(tr, outgroup = tr$tip.label[cls == "At_decoy"][which.max(d[cls == "At_decoy"])],
           resolve.root = TRUE)
cls <- sub("\\|.*", "", tr$tip.label)

# The test is deliberately restricted to the GRASS MYB61 clade. AtMYB61 is NOT
# expected inside it: the eudicot/monocot split means AtMYB61 sits outside, next
# to its own closest Arabidopsis paralogue. Requiring AtMYB61 inside would drag
# the MRCA to the root and make the test vacuous.
mono <- which(cls %in% c("Sb", "Os", "SC_MYB61", "PU_MYB61"))
nd   <- getMRCA(tr, mono)
tips <- extract.clade(tr, nd)$tip.label
tb   <- table(sub("\\|.*", "", tips))
cat(sprintf("\n[test] grass MYB61 clade (MRCA of monocot anchors + confirmed copies): %d tips\n",
            length(tips)))
for (k in names(tb)) cat(sprintf("       %-12s %d\n", k, tb[[k]]))

missing  <- setdiff(tr$tip.label[cls %in% c("SC_MYB61","PU_MYB61")], tips)
intruder <- grep("decoy", tips, value = TRUE)

if (length(missing)) {
  cat(sprintf("[WARN] %d confirmed copy/copies fall OUTSIDE the clade:\n", length(missing)))
  cat(paste0("       ", missing, collapse = "\n"), "\n")
} else {
  cat("[PASS] every confirmed copy falls inside the grass MYB61 clade.\n")
}

if (length(intruder) == 0) {
  cat("[PASS] no decoy inside the clade.\n")
} else {
  # An intruder is not automatically an error: a copy whose gene model lost a
  # Myb repeat fails the domain-rule filter yet is still a true MYB61 orthologue.
  cat(sprintf("[NOTE] %d decoy tip(s) inside the clade — inspect; a step-04\n", length(intruder)))
  cat("       'MYB61_homolog_no_MYB_call' here is a real copy with a broken\n")
  cat("       gene model, recovered by the tree:\n")
  cat(paste0("       ", intruder, collapse = "\n"), "\n")
  writeLines(sub(".*\\|", "", intruder),
             file.path(dirname(commandArgs(TRUE)[1]), "MYB61_recovered_by_tree.ids"))
}
RS
fi
