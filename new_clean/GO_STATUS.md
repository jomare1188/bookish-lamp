# GO annotation — live status

_2026-09-13 20:12:03_

## Baseline to beat (current, InterPro+Pfam)

| study | GO genes / network | responsive modules testable | terms clearing BH |
|---|---|---|---|
| sugarcane | 56,974 / 101,990 (55.9%) | 442 / 588 | 143 |
| purple | 94,497 / 170,135 (55.5%) | 106 / 182 | 8 |

## Stage 1 — eggNOG `--go_evidence all`

| study | tax_scope | GOs | Description | status |
|---|---|---|---|---|
| sugarcane | Poales | 22248/168135 (13.2%) | 149950/168135 (89.2%) | done (partial?) |
| sugarcane | auto | 81065/176001 (46.1%) | 156961/176001 (89.2%) | done (partial?) |
| purple | Poales | - | - | not run |
| purple | auto | 98867/229602 (43.1%) | 207932/229602 (90.6%) | done |

_original non-electronic run: sugarcane GOs 12,903/168,135 (7.7%), purple 15,572/211,596 (7.4%)_

## Stage 2 — full InterProScan (12 member DBs never run)

| study | proteins | chunks done | merged | status |
|---|---|---|---|---|
| sugarcane | 194593 | 98 / 98 | 323M | MERGED |
| purple | 241263 | 121 / 121 | 402M | MERGED |

## Stage 3 — curated transfer (Swiss-Prot Viridiplantae)  [DONE]

| study | pairs | genes | network coverage |
|---|---|---|---|
| sugarcane | 335509 | 51250 | 32,799 / 101,990 (32.2%) |
| purple | 366090 | 55267 | 45,233 / 170,135 (26.6%) |

## Stage 4 — PANNZER2  [NOT ADOPTED]

| study | chunks done | merged | status |
|---|---|---|---|
| sugarcane | 174 / 195 | none | 21 chunks died on ConnectTimeout to Helsinki |
| purple | 242 / 242 | 5,223,314 predictions / 170,151 proteins | complete |

_Output deleted 2026-09-13 (6.1 GB). Runner kept: `59_run_pannzer.sh`, hardened in b283178._

## Adoption decisions so far

| source | coverage | coherence (fixed gene set) | verdict |
|---|---|---|---|
| InterProScan full (17 DB) | 62.9% / 64.4% | sug +0.0748, pur +0.0361 | **ADOPTED** |
| eggNOG (tax_scope auto) | 46.1% / 43.1% of proteins | sug +0.0650, pur +0.0316 | rejected — worse on identical genes |
| curated (Swiss-Prot) | 32.2% / 26.6% | sug +0.0780, pur +0.0432 | not merged — see annotation/README.md |
| PANNZER2 | sugarcane incomplete | not scored | **dropped** — out by decision, and never finished |

## Stage 0/5 — the judge, and the merged table

- sugarcane  gene2go: 256099 pairs | coherence: scored
- purple     gene2go: 330621 pairs | coherence: not run

## Running now

- nothing

_regenerate: `./scripts/go_status.sh` · live: `watch -n 60 ./scripts/go_status.sh`_
