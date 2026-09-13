# Better GO for two grass genomes — progress

Branch `clustering-methods`. Four sources, tiered, **adopted only if they improve
module coherence** — coverage alone is the wrong target, since the last 7× coverage
gain produced *fewer* significant terms.

Baseline to beat: GO on 56,974/101,990 sugarcane (55.9%) and 94,497/170,135 purple
(55.5%); 442/588 and 106/182 responsive modules testable; 143 and 8 terms clearing
cross-module BH.

**CLOSED 2026-09-13. Outcome: ONE source, full InterProScan (17 DBs). No merge, no
tier.** The union scored BELOW either InterPro table alone on identical genes, in both
species, so merging was rejected by the judge this plan was built around rather than
skipped. eggNOG rejected; curated kept but not merged; PANNZER dropped and its output
deleted. See `annotation/README.md` and `docs/decisions.md` (2026-09-13).

Legend: `[ ]` todo · `[~]` running · `[x]` done · `[!]` blocked · `[-]` not adopted

---

## Stage 0 — the judge (build before adopting anything)

- [x] `55_go_coherence.r` — Sorensen-Dice above a size-matched null, on GO not PFAM
- [x] **sugarcane scored, 3 annotations on the SAME partition** (all normalised):

      annotation     genes    modules    H       null     excess
      current        54,903     2,906   0.0880   0.0260   +0.0620   <- best
      eggnog_auto    45,923     2,602   0.2531   0.1939   +0.0592
      union          66,725     3,209   0.1472   0.0899   +0.0572

      eggNOG's raw H is 3x higher and its excess is LOWER -- the null exposes that as
      term density (16.7 terms/gene), not shared function. Raw H would have picked it.
      The union buys +21% coverage and +303 scorable modules but dilutes the excess.
      **By the agreed rule eggNOG is not adopted** -- pending the InterProScan result,
      which changes the `current` baseline.
- [x] purple scored
- [x] **FLAW IN MY OWN JUDGE FOUND AND FIXED.** Scored on each annotation's own gene
      pool, EVERY coverage increase lowered the excess — because the null draws from
      the annotated pool, so a sparse annotation of well-characterised genes gets a low
      null and a big excess, while broader coverage adds generic terms that raise H and
      the null together. The metric was partly measuring specificity and would have
      rejected every expansion by construction. Fixed by scoring all candidates on the
      SAME genes and modules.
- [x] **fixed-set result, both species** (same genes, same modules):

      excess over null   nfcore5db   ipsfull   eggnog   union
      sugarcane           +0.0755    +0.0748   +0.0650  +0.0671
      purple              +0.0377    +0.0361   +0.0316  +0.0264

      InterPro-derived is equally coherent from 5 or 17 DBs (0.9% apart, within noise);
      eggNOG is worse on identical genes and drags the union down in BOTH species.

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
- [x] cost model measured: **369 s fixed startup + 1.71 s/protein** at -cpu 8
      (62 seqs 7m55s; 600 seqs 23m15s) -> 2,000/chunk = ~63 min, 218 chunks, ~8 h total
- [x] sugarcane: 98/98 chunks, 181,416/194,593 proteins with a signature (93.2%)
- [x] purple: 121/121 chunks, 233,300/241,263 (96.7%)
- [x] all 17 member DBs ran; 1.10M and 1.41M rows carry an InterPro accession
- [ ] sugarcane 194,593 proteins
- [ ] purple 241,263 proteins
- [ ] chunking verified lossless (every sequence in exactly one chunk)

## Stage 3 — curated transfer (the only non-IEA evidence)  [RUN, NOT MERGED]

- [x] `58_curated_transfer.sh` — diamond vs Swiss-Prot Viridiplantae + TAIR
- [x] experimental evidence codes only; no IEA survives
- [x] identity/coverage sensitivity curve reported, not a single asserted cut
- [-] **not merged.** Scores well per gene (+0.0780 / +0.0432) but reaches only
      32.2% / 26.6% of network genes and overlaps the adopted pairs by 6.5% / 5.4%.
      Merging it in would mean a tiered table, which the union result argues against.

## Stage 4 — PANNZER2  [DROPPED]

- [x] installed SANSPANZ / SANSparallel.3 (two upstream Python-3 bugs patched)
- [x] `59_run_pannzer.sh` (reuses the stage-2 splitter), hardened in `b283178`
- [-] **never finished:** purple 242/242 chunks (5,223,314 predictions over 170,151
      proteins); sugarcane stopped at 174/195, 21 chunks killed by `ConnectTimeout`
      to the public SANS service at Helsinki
- [-] **dropped by decision** before it was ever scored. Output deleted 2026-09-13
      (6.1 GB). No claim is made about whether it would have helped.

## Stage 5 — merge, tier, adopt

- [x] **ADOPTED: full InterProScan (ipsfull). REJECTED: eggNOG.** Same per-gene
      coherence as before, far more reach:

                    GO genes / network        testable modules   terms clearing BH
      sugarcane   56,974 -> 64,178 (62.9%)    442 -> 479 / 588     143 -> 164
      purple      94,497 -> 109,591 (64.4%)   106 -> 118 / 182       8 ->   9

      median annotated members per responsive module: sugarcane 4 -> 5, purple 3 -> 3
- [x] old 5-DB table kept as `gene2go_nfcore5db_<study>.tsv.bak`
- [-] **tier column: NOT BUILT.** The union of InterPro and eggNOG scores +0.0671
      (sugarcane) and +0.0264 (purple) against +0.0748 / +0.0361 for InterPro alone --
      merging LOWERS coherence in both species, so a tiered table would have cost a
      `source`/`tier` column on every pair and a caveat on every result to buy
      negative signal.
- [x] modulego re-run, both studies, **all three ontologies** on the adopted table:

                    responsive  testable   BP terms    MF terms    CC terms
                                           (clear BH)  (clear BH)  (clear BH)
      sugarcane        588        479         164         272          22
      purple           182        118           9          17           1

      MF is the stronger ontology in both species -- expected for a domain-derived
      annotation, which names what a protein DOES more sharply than what process it
      is in. BP stays primary (the question is a response) and is what figure 6 draws.
- [x] figure 6 rebuilt -- it was from 2026-08-23 and predated BOTH the Pearson-only
      rebuild and the ipsfull adoption
- [x] non-adopted tables suffixed `.notused`; `annotation/README.md` written
- [x] docs + commit

---

## Log

- 2026-09-11 — plan approved.
- 2026-09-11 — CORRECTION to my earlier diagnosis: I said emapper's `--go_evidence`
  default was THE cause of 7.7% GO. Measured: it explains 7.7 -> 13.2%. Widening
  `--tax_scope` from Poales to auto is worth another 33 points (-> 46.1%).
- 2026-09-11 — eggNOG GO is pre-propagated to the DAG roots; every source now passes
  through 60_normalise_gene2go.r first, or neither the coherence test nor topGO is valid.
