# GO annotation — live status

_2026-09-12 12:12:21_

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

## Stage 3 — curated transfer (Swiss-Prot + TAIR)

- sugarcane  not run
- purple     not run

## Stage 4 — PANNZER2

- sugarcane  not run
- purple     not run

## Stage 0/5 — the judge, and the merged table

- sugarcane  gene2go: 256099 pairs | coherence: scored
- purple     gene2go: 330621 pairs | coherence: not run

## Running now

- nothing

_regenerate: `./scripts/go_status.sh` · live: `watch -n 60 ./scripts/go_status.sh`_
