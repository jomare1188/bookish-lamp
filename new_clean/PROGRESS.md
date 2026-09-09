# Clustering methods benchmark — progress

Branch `clustering-methods`. Tree `results_cluster/`, work dir
`/dados04/jorge/tmp/mcl_work_cluster`. Nothing touches `results/`.

Graph: **unpruned Pearson-only**, `|r| >= 0.8`, weights rescaled to `[0.01, 1]`,
no k-NN reduction.

Legend: `[ ]` todo · `[~]` running · `[x]` done · `[!]` blocked/failed

---

## A. Scaffolding

- [x] branch `clustering-methods` off `pearson-knn-ba`
- [x] `.gitignore`: `new_clean/results_cluster/`
- [x] `config.sh`: `CLUSTER_*` knobs (inflation ladder, mcl resource probe, Leiden gamma, python)
- [x] `run.sh`: stage `pearsonmci` wired (rest as they land) — `pearsonmci leidensweep clustercompare figclustermethods`

## B. Build the matrices (`46_pearson_mci.sh`)

Streams the Pearson layer straight into `mcxload` — no 62 GB intermediate table.

- [x] `46_pearson_mci.sh` written
- [x] sugarcane `.mci` built — **2m41s**, 1.2 GB, peak RSS 2.4 GB
- [x] sugarcane verified: 101,990 nodes, **75,333,769 edges (exact match)**, 0 degree-0, mean degree 1,477
- [x] purple `.mci` built — **23m25s**, 10.8 GB, peak RSS 21.7 GB
- [x] purple verified: 170,135 nodes, **675,955,918 edges (exact match)**, 0 degree-0, mean degree 7,946

## C. The correctness gate — before any comparison is believed

Leiden partitions get scored by `clm info`, so the cluster-file writer must be
exactly right. A silent off-by-one here invalidates the whole benchmark.

- [x] **PASS** — `verify_cluster_roundtrip.py`, sugarcane `cls.knone.I12`: all 7 stats
      identical (eff 0.11458, mod 0.08018, mf 0.88508, af 0.12333, ncl 2350, max 32055, sgl 0)
- [x] negative control added: dropping the permutation *does* change the scores
      (eff 0.11458 → 0.11438), so the test is capable of failing
- [ ] re-run the gate on purple (its permutation may be more scrambled than sugarcane's 4.7%)
- [x] enforced in `47_leiden_sweep.py` (`write_cls` exits if coverage != N) and asserted in the gate

## D. MCL inflation ladder (`36_mcl_sweep.sh`, k-NN axis off)

Ladder `1.2 1.4 1.7 2 2.5 3 4 6`. The historical `-I 2` is one cell among eight.

- [x] `36_mcl_sweep.sh`: `CLEAN_MCL_RESOURCE` knob + `partitions.tsv` manifest + `mkdir -p OUT_DIR`
- [x] sugarcane ladder done (1.2 … 6) — giant module only 31.4% → 15.6%, `eff` still
      climbing at -I 6, so the ladder was extended
- [x] sugarcane extension: `-I 8` valid; `-I 12` dropped (tag collision); `-I 40` DELETED
      (mcl rejects >30 and silently used the default, giving `-I 2`'s partition byte-for-byte);
      `-I 20` kept but flagged — 35,044 of 101,990 vectors underflowed to zero
- [x] `36_mcl_sweep.sh` hardened: refuses out-of-range `-I`, refuses tag collisions
      (within a grid AND across runs via a `.inflation` sidecar), flags underflow,
      and its TSV is now additive instead of truncating. Both guards tested firing.
- [~] sugarcane `-I 9 10 11 13 15 16` running to locate the `eff` peak
- [~] purple ladder running — `logs/mclladder_purple.log`

## E. The pruning probe — is the ladder measuring inflation or MCL's pruner?

MCL keeps only S=1200 neighbours per node while computing. Purple's mean degree
is 7,946, so ~85% is discarded on the fly; MCL grades its own pruning "awful"
(20.1). Sugarcane: mean degree 1,477, "deplorable" (39.2).

- [ ] sugarcane: winning inflation re-run at `-S 4000` and `-S 10000`
- [ ] purple: same  (expect 2-5 h per cell)
- [ ] `clm dist` between the `-S 1200` and `-S 10000` partitions — reported, not buried

## F. Leiden / Louvain (`47_leiden_sweep.py`, `sbm` env)

- [x] `47_leiden_sweep.py` written (+ idempotent skip of completed cells)
- [x] sugarcane gamma scout done — clusters 19,571→76,387, largest 15.73%→4.29%,
      still moving at 0.5, so ladder set to `0.005 0.01 0.025 0.05 0.1 0.2 0.35 0.5 0.7 0.9`
- [x] sugarcane: Leiden CPM ladder (10 gammas) + Leiden modularity + Louvain
      KEY RESULT: Leiden CPM dissolves the giant module (17.2% -> 0.31%), while
      Leiden-modularity (42.9%) and Louvain (40.8%) make it WORSE than MCL's 15.6%.
      The giant module is what modularity maximisation wants, not an MCL artefact.
- [~] purple gamma scout running (load took 28 min, 35 GB resident; mean weight 0.3796)

## G. Scoring, all methods on the same graph (`48_cluster_compare.sh`)

`eff` is the primary criterion. `mf` is reported but never used alone — the
trivial one-cluster partition scores `mf = 1.0`.

- [x] `48_cluster_compare.sh` written (manifest-driven, one `clm info` over all methods)
- [ ] `cluster_methods_sugarcane.tsv`
- [ ] `cluster_methods_purple.tsv`
- [ ] one-cluster baseline included as the row that makes the `mf` trap obvious

## H. Independent biological check (`37_cluster_homogeneity.r`)

PFAM Sørensen–Dice homogeneity above a size-matched null. A setting has to win
on the graph metric *and* here.

- [x] `37` reads `partitions.tsv` and carries method/param/resource columns
- [ ] sugarcane scored
- [ ] purple scored

## I. Hierarchical scout — sugarcane only (`50_fastgreedy_scout.py`)

Classical hierarchical is not attempted (needs a dense 232 GB distance matrix
built from correlations we thresholded away). fast-greedy (CNM) is the
graph-native stand-in, under a 4 h cap.

- [ ] script written
- [ ] sugarcane run, or cap hit and recorded as the result

## J. Report

- [ ] `49_fig_cluster_methods.r` → `figure12_cluster_methods`
- [ ] `docs/results.md`: the comparison and the winner
- [ ] `docs/decisions.md`: why inflation is now chosen, and whether Leiden is adopted
- [ ] commit

---

## Log

- 2026-09-09 — branch created, config knobs added.
- 2026-09-09 — sugarcane matrix built and both assertions passed; purple building.
- 2026-09-09 17:15 — correctness gate PASSED with a negative control. Leiden and MCL
  can now be compared on the same numbers.
- 2026-09-09 17:27 — purple matrix built and verified (675,955,918 edges, exact).
- 2026-09-09 18:30 — TWO silent mcl failures caught, both now guarded: -I above 30 is
  ignored (default used instead), and -I 20 underflows 34% of the graph to zero.
- 2026-09-09 18:10 — BUG FOUND: `36_mcl_sweep.sh` cell tags strip the decimal point, so
  `-I 1.2` and `-I 12` collide on `cls.knone.I12`. The second is silently skipped and
  reported under the wrong inflation. Purple's grid has no collision; the sugarcane
  extension's `-I 12` row is discarded. Fix queued in /tmp/pending_fix.md (record the
  actual inflation beside each cell and refuse a mismatch).
- 2026-09-09 — NOTE: editing `run.sh`/`config.sh` while a job is reading them corrupts
  the running copy (bash reads scripts incrementally). It cost a spurious exit-2 on the
  purple build, after the work had completed. No edits to in-use scripts from here.
