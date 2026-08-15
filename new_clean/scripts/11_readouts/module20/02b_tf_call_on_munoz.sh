#!/usr/bin/env bash
# ============================================================================
# 02b_tf_call_on_munoz.sh — run OUR TF identification on MUÑOZ'S OWN proteins
#
# WHY THIS SETTLES SOMETHING
# The blastx run found only ~32% of the mapped R570/purple genes to be MYB,
# against Muñoz's published 75%. Two explanations were possible:
#   (a) a METHOD difference — they classified TFs differently, or
#   (b) a MAPPING artefact — the orthologs we mapped to are partly not MYB.
#
# `sugar_cane.pep` turns out to be byte-identical to `GET_TFS/db/sugar_cane.pep`,
# the proteome the original GET_TFS run used — so their Family column very
# likely comes from the same PlnTFDB rule method this project uses. Applying
# that exact method (hmmsearch --cut_ga vs TF.db.hmm, then
# assign_family_membership.pl + RulesFull) to the Module-20 ORFs and comparing
# to `module20.csv` discriminates (a) from (b) directly.
#
# If the calls reproduce, the 32% is about the ORTHOLOGS we mapped to, not
# about how anyone classifies TFs.
#
# RUN: ./02b_tf_call_on_munoz.sh   (seconds - only ~23 proteins)
# ============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

DOM="${WORKDIR}/module20_munoz.domtbl"
ASSIGN="${WORKDIR}/module20_munoz_family.tsv"
OUT="${OUTDIR}/module20_tf_call_comparison.tsv"

[ -s "$M20_FAA" ] || { echo "ERROR: run 01_extract_module20.sh first" >&2; exit 1; }

echo "[hmm ] hmmsearch --cut_ga vs $(basename "$HMM_DB") on $(grep -c '^>' "$M20_FAA") ORFs"
hmmsearch --cpu "$THREADS" --cut_ga --domtblout "$DOM" -o /dev/null "$HMM_DB" "$M20_FAA"

echo "[rule] assign_family_membership.pl + RulesFull"
perl "$ASSIGN_PL" --pfam "$DOM" --rules "$RULES" --species module20 --out "$ASSIGN"

# our call per ORF, then per transcript (a member is family X if any of its ORFs is)
"$RSCRIPT_PLOT" --vanilla - <<'RS'
suppressMessages(library(data.table))
o <- Sys.getenv("OUTDIR"); w <- Sys.getenv("WORKDIR")

a <- tryCatch(fread(file.path(w, "module20_munoz_family.tsv"), header = FALSE),
              error = function(e) data.table(V1 = character(), V2 = character(), V3 = character()))
if (nrow(a)) {
  # drop a header row if the perl script wrote one, and drop Orphans
  a <- a[V1 != "" & V1 != "Gene" & V2 != "Orphans"]
  setnames(a, 1:3, c("orf", "our_family", "our_type"))
} else a <- data.table(orf = character(), our_family = character(), our_type = character())

orfs <- fread(file.path(o, "module20_orfs.tsv"))
mem  <- fread(file.path(o, "module20_members.tsv"))

a <- merge(a, orfs[, .(orf, transcript, orf_type, orf_len)], by = "orf", all.x = TRUE)
fwrite(a, file.path(o, "module20_our_tf_call_per_orf.tsv"), sep = "\t")

per_tr <- a[, .(our_call = paste(sort(unique(our_family)), collapse = ";"),
                our_type = paste(sort(unique(our_type)), collapse = ";"),
                calling_orf = orf[which.max(orf_len)]), by = transcript]

cmp <- merge(mem[, .(transcript, munoz_family = family, munoz_class = class)],
             per_tr, by = "transcript", all.x = TRUE)
cmp[is.na(our_call), our_call := "no_TF_call"]
# Muñoz writes MYB_related, the rule set writes MYB-related
cmp[, munoz_norm := gsub("_", "-", munoz_family)]
cmp[, agree := fifelse(munoz_norm == "no-tf", our_call == "no_TF_call",
                       mapply(function(a, b) any(strsplit(b, ";")[[1]] == a),
                              munoz_norm, our_call))]
setorder(cmp, munoz_family, transcript)
fwrite(cmp[, .(transcript, munoz_family, munoz_class, our_call, our_type, calling_orf, agree)],
       file.path(o, "module20_tf_call_comparison.tsv"), sep = "\t")

cat("\n=== our TF call on Muñoz's own proteins vs their published Family ===\n")
print(cmp[, .(transcript, munoz_family, our_call, agree)], row.names = FALSE)
cat(sprintf("\nagreement: %d / %d members\n", sum(cmp$agree, na.rm = TRUE), nrow(cmp)))
mm <- cmp[munoz_norm %in% c("MYB", "MYB-related")]
cat(sprintf("MYB-family members reproduced by our method: %d / %d\n",
            sum(mm$agree, na.rm = TRUE), nrow(mm)))
RS

echo "[ok  ] -> ${OUT}"
