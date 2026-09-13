#!/usr/bin/env bash
# =============================================================================
# go_status.sh -- live progress of the GO annotation work.
#
# Derives everything from the filesystem and from running processes, so it cannot
# drift out of date the way a hand-maintained checklist does. Safe to run at any
# time; writes GO_STATUS.md and prints the same thing.
#
#   ./scripts/go_status.sh              once
#   watch -n 60 ./scripts/go_status.sh  live
# =============================================================================
BASE=/dados04/jorge/comparative_saccharum
CLEAN="$BASE/new_clean"
A="$BASE/annotation"
OUT="$CLEAN/GO_STATUS.md"

fill() { [ -s "$1" ] && awk -F'\t' -v c="$2" '!/^#/{n++; if($c!="-" && $c!="")k++}
         END{printf "%d/%d (%.1f%%)", k+0, n, n?100*k/n:0}' "$1" || echo "-"; }
nseq()  { [ -s "$1" ] && grep -c '^>' "$1" 2>/dev/null || echo 0; }
alive() { pgrep -fa "$1" 2>/dev/null | grep -v 'go_status' | head -1; }

{
printf '# GO annotation — live status\n\n_%s_\n\n' "$(date '+%Y-%m-%d %H:%M:%S')"

printf '## Baseline to beat (current, InterPro+Pfam)\n\n'
printf '| study | GO genes / network | responsive modules testable | terms clearing BH |\n'
printf '|---|---|---|---|\n'
printf '| sugarcane | 56,974 / 101,990 (55.9%%) | 442 / 588 | 143 |\n'
printf '| purple | 94,497 / 170,135 (55.5%%) | 106 / 182 | 8 |\n\n'

printf '## Stage 1 — eggNOG `--go_evidence all`\n\n'
printf '| study | tax_scope | GOs | Description | status |\n|---|---|---|---|---|\n'
for s in sugarcane purple; do
  for v in "Poales:eggnog_go_all" "auto:eggnog_go_all_auto"; do
    sc="${v%%:*}"; d="$A/$s/${v##*:}/$s.emapper.annotations"
    st="not run"; [ -s "$d" ] && st="done (partial?)"
    [ -s "$d" ] && [ -f "$(dirname "$d")/.complete" ] && st="done"
    [ -n "$(alive "annotate_hits_table.*--output_dir $(dirname "$d")")" ] && st="**RUNNING**"
    printf '| %s | %s | %s | %s | %s |\n' "$s" "$sc" "$(fill "$d" 10)" "$(fill "$d" 8)" "$st"
  done
done
printf '\n_original non-electronic run: sugarcane GOs 12,903/168,135 (7.7%%), purple 15,572/211,596 (7.4%%)_\n\n'

printf '## Stage 2 — full InterProScan (12 member DBs never run)\n\n'
printf '| study | proteins | chunks done | merged | status |\n|---|---|---|---|---|\n'
for s in sugarcane purple; do
  cd_="$A/$s/interproscan_full"
  tot=$( [ -f "$cd_/chunks.total" ] && cat "$cd_/chunks.total" || echo "?" )
  done_=$(ls "$cd_"/done/*.done 2>/dev/null | wc -l)
  merged="$cd_/$s.interproscan_full.tsv"
  st="not run"; [ -d "$cd_" ] && st="prepared"
  [ "$done_" -gt 0 ] && st="running ($done_/$tot)"
  [ -s "$merged" ] && st="MERGED"
  printf '| %s | %s | %s / %s | %s | %s |\n' "$s" \
    "$(nseq "$BASE/files/fix_orthofinder/proteins/$( [ "$s" = sugarcane ] && echo sugarcane_one_transcript.fa || echo one_transcript_purple_proteins.faa )")" \
    "$done_" "$tot" "$( [ -s "$merged" ] && du -h "$merged" | cut -f1 || echo '-' )" "$st"
done
printf '\n'

printf '## Stage 3 — curated transfer (Swiss-Prot Viridiplantae)  [DONE]\n\n'
printf '| study | pairs | genes | network coverage |\n|---|---|---|---|\n'
for s in sugarcane purple; do
  # The non-adopted tables were suffixed .notused on 2026-09-13 so nothing can read
  # one by accident (annotation/README.md). Resolve either name: reporting "not run"
  # for a source that WAS run and then rejected would misrepresent the decision.
  f="$A/$s/gene2go_curated_$s.tsv"
  [ -s "$f" ] || f="$A/$s/gene2go_curated_$s.tsv.notused"
  if [ -s "$f" ]; then
    printf '| %s | %s | %s | %s |\n' "$s" "$(( $(wc -l < "$f") - 1 ))" \
      "$(awk -F'\t' 'NR>1{g[$1]}END{print length(g)}' "$f")" \
      "$( [ "$s" = sugarcane ] && echo '32,799 / 101,990 (32.2%)' || echo '45,233 / 170,135 (26.6%)' )"
  else printf '| %s | not run | | |\n' "$s"; fi
done

printf '\n## Stage 4 — PANNZER2  [NOT ADOPTED]\n\n'
# The output tree was DELETED on 2026-09-13, so there is nothing left to poll and the
# live probe that used to be here would report "not run" for a source that was in fact
# run and then dropped. These are the final counts, recorded before deletion.
printf '| study | chunks done | merged | status |\n|---|---|---|---|\n'
printf '| sugarcane | 174 / 195 | none | 21 chunks died on ConnectTimeout to Helsinki |\n'
printf '| purple | 242 / 242 | 5,223,314 predictions / 170,151 proteins | complete |\n'
printf '\n_Output deleted 2026-09-13 (6.1 GB). Runner kept: `59_run_pannzer.sh`, hardened in b283178._\n'

printf '\n## Adoption decisions so far\n\n'
printf '| source | coverage | coherence (fixed gene set) | verdict |\n|---|---|---|---|\n'
printf '| InterProScan full (17 DB) | 62.9%% / 64.4%% | sug +0.0748, pur +0.0361 | **ADOPTED** |\n'
printf '| eggNOG (tax_scope auto) | 46.1%% / 43.1%% of proteins | sug +0.0650, pur +0.0316 | rejected — worse on identical genes |\n'
printf '| curated (Swiss-Prot) | 32.2%% / 26.6%% | sug +0.0780, pur +0.0432 | not merged — see annotation/README.md |\n'
printf '| PANNZER2 | sugarcane incomplete | not scored | **dropped** — out by decision, and never finished |\n'

printf '\n## Stage 0/5 — the judge, and the merged table\n\n'
for s in sugarcane purple; do
  g="$A/$s/gene2go_$s.tsv"; c="$CLEAN/results/$s/go_coherence_$s.tsv"
  printf -- '- %-10s gene2go: %s | coherence: %s\n' "$s" \
    "$( [ -s "$g" ] && echo "$(( $(wc -l < "$g") - 1 )) pairs" || echo '-' )" \
    "$( [ -s "$c" ] && echo "scored" || echo 'not run' )"
done

printf '\n## Running now\n\n'
found=0
for p in annotate_hits_table interproscan runsanspanz diamond 55_go_coherence; do
  l=$(alive "$p"); [ -n "$l" ] && { printf -- '- `%s`\n' "$(echo "$l" | cut -c1-110)"; found=1; }
done
[ "$found" = 0 ] && printf -- '- nothing\n'
printf '\n_regenerate: `./scripts/go_status.sh` · live: `watch -n 60 ./scripts/go_status.sh`_\n'
} > "$OUT"
cat "$OUT"
