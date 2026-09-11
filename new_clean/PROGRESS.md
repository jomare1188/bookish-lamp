# Better GO for two grass genomes — progress

Branch `clustering-methods`. Four sources, tiered, **adopted only if they improve
module coherence** — coverage alone is the wrong target, since the last 7× coverage
gain produced *fewer* significant terms.

Baseline to beat: GO on 56,974/101,990 sugarcane (55.9%) and 94,497/170,135 purple
(55.5%); 442/588 and 106/182 responsive modules testable; 143 and 8 terms clearing
cross-module BH.

Legend: `[ ]` todo · `[~]` running · `[x]` done · `[!]` blocked

---

## Stage 0 — the judge (build before adopting anything)

- [x] `55_go_coherence.r` — Sorensen-Dice above a size-matched null, on GO not PFAM
- [~] sugarcane: current vs eggnog_auto, both normalised (running)
- [ ] purple

## Stage 1 — eggNOG re-annotation  [DONE]

- [x] `56_run_eggnog_all.sh` (reuses saved seed_orthologs; ~8 min per study)
- [x] verified only GOs moved: rows, eggNOG_OGs, Description, PFAMs all identical
- [x] **`--tax_scope` matters far more than `--go_evidence`** (corrects my earlier claim):

      sugarcane   Poales + non-electronic (original)   12,903 / 168,135   7.7%
                  Poales + go_evidence all             22,248 / 168,135  13.2%
                  auto   + go_evidence all             81,065 / 176,001  46.1%
      purple      auto   + go_evidence all             98,867 / 229,602  43.1%  (was 7.4%)

## Stage 1b — normalise every source to most-specific terms  [DONE]

- [x] `60_normalise_gene2go.r`. **eggNOG ships the full ancestor closure INCLUDING the
      three GO roots** (69.9 terms/gene; `biological_process` on 75,696 sugarcane genes),
      while InterPro+Pfam gives direct terms only (2.3, zero roots). Mixing them would
      (a) let the propagated source win the Dice coherence test on shared generalities
      and (b) make topGO double-count, since it propagates internally.
- [x] eggNOG 69.9 -> 16.7 terms/gene (4.3M implied pairs removed); current 2.3 -> 2.1

## Stage 2 — full InterProScan, the 12 member DBs never run (LONG POLE)

Local 5.78 install, 49 GB data. nf-core used 5 of 17: missing Gene3D, SUPERFAMILY,
CDD, SMART, PRINTS, PROSITE, Pfam.

- [x] `57_run_interproscan.sh` — chunk + queue + merge, resumable
- [x] invocation tested: 14 member DBs run (CDD, Gene3D, SUPERFAMILY, SMART, PRINTS,
      ProSite, Pfam ... all the missing ones). 62 seqs / 7m55s / 4.9 GB at -cpu 8
- [~] timing test at 600 seqs to separate fixed startup from marginal cost
- [ ] sugarcane 194,593 proteins
- [ ] purple 241,263 proteins
- [ ] chunking verified lossless (every sequence in exactly one chunk)

## Stage 3 — curated transfer (the only non-IEA evidence)

- [ ] `58_curated_transfer.sh` — diamond vs Swiss-Prot Viridiplantae + TAIR
- [ ] experimental evidence codes only; assert no IEA survives
- [ ] identity/coverage sensitivity curve reported, not a single asserted cut

## Stage 4 — PANNZER2

- [ ] install SANSPANZ / SANSparallel.3
- [ ] `59_run_pannzer.sh` (reuses the stage-2 splitter)
- [ ] score threshold chosen on the coherence curve, not the default

## Stage 5 — merge, tier, adopt

- [ ] `54_build_gene2go.sh` merges all sources with `source` + `tier`
- [ ] coherence delta per source, beside coverage
- [ ] adopt / reject each source on that evidence
- [ ] re-run modulego both studies; report coverage, testable, coherence, terms
- [ ] docs + commit

---

## Log

- 2026-09-11 — plan approved.
- 2026-09-11 — CORRECTION to my earlier diagnosis: I said emapper's `--go_evidence`
  default was THE cause of 7.7% GO. Measured: it explains 7.7 -> 13.2%. Widening
  `--tax_scope` from Poales to auto is worth another 33 points (-> 46.1%).
- 2026-09-11 — eggNOG GO is pre-propagated to the DAG roots; every source now passes
  through 60_normalise_gene2go.r first, or neither the coherence test nor topGO is valid.
