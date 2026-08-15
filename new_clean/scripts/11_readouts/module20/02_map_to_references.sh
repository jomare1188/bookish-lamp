#!/usr/bin/env bash
# ============================================================================
# 02_map_to_references.sh — bridge Module-20 PROTEINS to both reference proteomes
#
# Protein-vs-protein (`diamond blastp`) against the Muñoz proteome ORFs, which
# replaces the earlier translated `blastx` of the nucleotide contigs. Better
# because: no frame shifts, no UTR in the query, and query coverage is finally
# meaningful (with blastx, qcov was measured against a contig that includes
# UTRs, which is why only subject coverage could be filtered on).
#
# Specificity still comes from the RECIPROCAL-ARABIDOPSIS filter:
#
#   member ORFs --blastp--> Arabidopsis   => at_anchor(member)   [best ORF wins]
#   member ORFs --blastp--> R570 / purple => candidate genes
#   candidate   --blastp--> Arabidopsis   => at_best(candidate)
#   keep candidate iff  at_best(candidate) == at_anchor(member)
#
# Each member has 2-4 TransDecoder ORFs, some antisense/internal artefacts, and
# "longest" and "complete" disagree for at least one member. So ALL ORFs are
# searched and the best-scoring one is attributed to the member and recorded
# (`orf` column) rather than choosing by a fixed rule.
#
# Mapping is done INDEPENDENTLY into each proteome, never sugarcane->purple, so
# the purple result is not conditioned on the sugarcane one.
#
# RUN: ./02_map_to_references.sh
# ============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

BL6='qseqid sseqid pident length evalue bitscore qcovhsp scovhsp'
[ -s "$M20_FAA" ] || { echo "ERROR: run 01_extract_module20.sh first" >&2; exit 1; }

for db in "$DB_SC:$SC_CLEAN" "$DB_PU:$PU_CLEAN" "$DB_AT:$AT_PEP"; do
  d="${db%%:*}"; f="${db#*:}"
  [ -s "${d}.dmnd" ] || { echo "[db  ] building $(basename "$d")"; diamond makedb --quiet --in "$f" -d "$d"; }
done

# --- anchor: best Arabidopsis hit per member (best ORF wins) ----------------
echo "[anch] Module-20 ORFs -> Arabidopsis"
diamond blastp --quiet -q "$M20_FAA" -d "$DB_AT" -o "${WORKDIR}/m20_vs_at.tsv" \
    --outfmt 6 $BL6 --evalue 1e-5 --max-target-seqs 5 --threads "$THREADS"

# orf -> best At hit, then collapse orf -> transcript keeping the strongest
awk -v OFS='\t' '{ if (!($1 in bs) || $6 > bs[$1]) { bs[$1] = $6; best[$1] = $2 } }
                  END { for (q in best) print q, best[q], bs[q] }' \
    "${WORKDIR}/m20_vs_at.tsv" > "${WORKDIR}/m20_orf_at.tsv"
awk -v OFS='\t' '{ tr = $1; sub(/\.p[0-9]+$/, "", tr)
                   if (!(tr in bs) || $3 > bs[tr]) { bs[tr] = $3; a[tr] = $2; o[tr] = $1 } }
                 END { for (t in a) print t, a[t], o[t], bs[t] }' \
    "${WORKDIR}/m20_orf_at.tsv" | sort > "${WORKDIR}/m20_at_anchor.tsv"

echo "[ok  ] $(wc -l < "${WORKDIR}/m20_at_anchor.tsv") / $(awk -F'\t' 'NR>1{t[$2]}END{print length(t)}' "${OUTDIR}/module20_orfs.tsv") members have an Arabidopsis anchor"
awk -F'\t' '{printf "       %-32s %-12s via %s\n", $1, $2, $3}' "${WORKDIR}/m20_at_anchor.tsv"

map_species() {
  local sp="$1" db="$2" clean="$3"
  local fwd="${WORKDIR}/m20_vs_${sp}.tsv" bwd="${WORKDIR}/${sp}_cand_vs_at.tsv"
  local out="${OUTDIR}/module20_map_${sp}.tsv"

  echo "[srch] Module-20 ORFs -> ${sp}"
  diamond blastp --quiet -q "$M20_FAA" -d "$db" -o "$fwd" --outfmt 6 $BL6 \
      --evalue "$EVALUE" --max-target-seqs "$MAX_TARGETS" --threads "$THREADS"

  # best HSP per (ORF, target protein) above the subject-coverage floor,
  # then collapse ORFs to their member transcript keeping the best-scoring ORF
  awk -v OFS='\t' -v mincov="$MIN_SCOV" '
    $8 >= mincov {
      tr = $1; sub(/\.p[0-9]+$/, "", tr)
      k = tr "\t" $2
      if (!(k in bs) || $6 > bs[k]) { bs[k] = $6; orf[k] = $1
        pid[k] = $3; ev[k] = $5; qc[k] = $7; sc[k] = $8 }
    }
    END { for (k in bs) { split(k, f, "\t")
          print f[1], f[2], orf[k], pid[k], ev[k], bs[k], qc[k], sc[k] } }' \
    "$fwd" > "${WORKDIR}/${sp}_fwd_best.tsv"

  cut -f2 "${WORKDIR}/${sp}_fwd_best.tsv" | sort -u > "${WORKDIR}/${sp}_cand.ids"
  echo "       $(wc -l < "${WORKDIR}/${sp}_cand.ids") candidate proteins (scov>=${MIN_SCOV}%)"

  awk 'NR==FNR { w[$1]; next } /^>/ { k = (substr($1,2) in w) } k' \
      "${WORKDIR}/${sp}_cand.ids" "$clean" > "${WORKDIR}/${sp}_cand.faa"

  echo "[rbh ] ${sp} candidates -> Arabidopsis"
  diamond blastp --quiet -q "${WORKDIR}/${sp}_cand.faa" -d "$DB_AT" -o "$bwd" \
      --outfmt 6 $BL6 --evalue 1e-5 --max-target-seqs 5 --threads "$THREADS"
  awk -v OFS='\t' '{ if (!($1 in bs) || $6 > bs[$1]) { bs[$1] = $6; best[$1] = $2 } }
                    END { for (q in best) print q, best[q] }' "$bwd" \
      > "${WORKDIR}/${sp}_cand_at_best.tsv"

  awk -F'\t' -v OFS='\t' -v sp="$sp" '
    FILENAME == ARGV[1] { anchor[$1] = $2; next }              # member -> At
    FILENAME == ARGV[2] { atbest[$1] = $2; next }              # protein -> At
    {
      tr = $1; prot = $2; orf = $3; gene = prot
      if (sp == "sugarcane") sub(/[.][0-9]+[.]p[0-9]*$/, "", gene)
      a = (tr   in anchor) ? anchor[tr]   : "none"
      b = (prot in atbest) ? atbest[prot] : "none"
      v = (a == "none") ? "no_at_anchor" : (a == b ? "reciprocal" : "rejected")
      print tr, gene, prot, orf, $4, $5, $6, $7, $8, a, b, v
    }' "${WORKDIR}/m20_at_anchor.tsv" "${WORKDIR}/${sp}_cand_at_best.tsv" \
       "${WORKDIR}/${sp}_fwd_best.tsv" \
  | sort -k1,1 -k7,7gr > "${out}.body"

  { printf 'transcript\tgene\tprotein\torf\tpident\tevalue\tbitscore\tqcov\tscov\tat_anchor\tat_best_hit\tverdict\n'
    cat "${out}.body"; } > "$out"
  rm -f "${out}.body"

  echo "[ok  ] ${sp} -> ${out##*/}"
  awk -F'\t' 'NR>1 {n[$12]++; if ($12=="reciprocal") {g[$2]; t[$1]}}
              END { for (v in n) printf "       %-14s %d protein hits\n", v, n[v]
                    printf "       => %d distinct genes from %d members\n", length(g), length(t) }' "$out"
}

map_species sugarcane "$DB_SC" "$SC_CLEAN"
map_species purple    "$DB_PU" "$PU_CLEAN"

# --- mapping summary + blastx/blastp comparison -----------------------------
"$RSCRIPT_PLOT" --vanilla - <<'RS'
suppressMessages(library(data.table))
o <- Sys.getenv("OUTDIR"); w <- Sys.getenv("WORKDIR")
m <- fread(file.path(o, "module20_members.tsv"))
f <- function(sp) fread(file.path(o, sprintf("module20_map_%s.tsv", sp)))[
        verdict == "reciprocal", .(n = uniqueN(gene)), by = transcript]
s <- f("sugarcane"); setnames(s, "n", "sugarcane_genes")
p <- f("purple");    setnames(p, "n", "purple_genes")
r <- merge(merge(m[, .(transcript, family, munoz_degree, munoz_betweenness)],
                 s, by = "transcript", all.x = TRUE),
           p, by = "transcript", all.x = TRUE)
r[is.na(sugarcane_genes), sugarcane_genes := 0][is.na(purple_genes), purple_genes := 0]
setorder(r, -sugarcane_genes)
print(r, row.names = FALSE)
fwrite(r, file.path(o, "module20_mapping_summary.tsv"), sep = "\t")

# how did the protein-level bridge change things vs the earlier blastx run?
old_f <- file.path(w, sprintf("module20_map_%s_blastx.tsv", c("sugarcane", "purple")))
if (all(file.exists(old_f))) {
  cat("\n=== blastp (proteome) vs blastx (contigs) ===\n")
  for (i in seq_along(c("sugarcane", "purple"))) {
    sp <- c("sugarcane", "purple")[i]
    old <- fread(old_f[i])[verdict == "reciprocal", unique(gene)]
    new <- fread(file.path(o, sprintf("module20_map_%s.tsv", sp)))[verdict == "reciprocal", unique(gene)]
    cat(sprintf("%-10s blastx %d genes | blastp %d genes | shared %d | blastp-only %d | blastx-only %d\n",
                sp, length(old), length(new), length(intersect(old, new)),
                length(setdiff(new, old)), length(setdiff(old, new))))
  }
}
RS
