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
- [x] gate re-run on purple (`cls.knone.I2`): **PASS**, all 7 statistics identical, negative control fires
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
- [x] sugarcane `-I 9 10 11 13` done; `-I 15` killed (pathological: 2,027 iterations,
      14,141 underflows, mcl auto-escalated inflation to 109)
- [x] **`eff` is MONOTONE for MCL across the whole valid range** (0.115 at 1.2 -> 0.496 at 13,
      still climbing) while singletons go 0 -> 15,344. Argmax on the grid boundary, so
      `eff` cannot choose an inflation — the same failure that sank the BA criterion.
      Only modularity has an interior peak, at **-I 1.7 (0.0818)**.


## E. The pruning probe — is the ladder measuring inflation or MCL's pruner?

MCL keeps only S=1200 neighbours per node while computing. Purple's mean degree
is 7,946, so ~85% is discarded on the fly; MCL grades its own pruning "awful"
(20.1). Sugarcane: mean degree 1,477, "deplorable" (39.2).

- [x] sugarcane `-I 1.7` at `-S 4000` and `-S 10000`: **identical to each other** (7,106
      clusters, largest 20.39%, eff 0.27964, mod 0.08239, jury 56.2 *acceptable*) and within
      0.7% of the default (7,064 / 21.02% / 0.27685 / 0.08184 / jury 37.1 *deplorable*).
      Pruning SATURATES at S=4000 and the converged answer is the default answer.
      **The jury grade is not a guide to whether a result is affected.**

- [x] the -S comparison is reported in docs/results.md, not buried

## F. Leiden / Louvain (`47_leiden_sweep.py`, `sbm` env)

- [x] `47_leiden_sweep.py` written (+ idempotent skip of completed cells)
- [x] sugarcane gamma scout done — clusters 19,571→76,387, largest 15.73%→4.29%,
      still moving at 0.5, so ladder set to `0.005 0.01 0.025 0.05 0.1 0.2 0.35 0.5 0.7 0.9`
- [x] sugarcane: Leiden CPM ladder (10 gammas) + Leiden modularity + Louvain
      KEY RESULT: Leiden CPM dissolves the giant module (17.2% -> 0.31%), while
      Leiden-modularity (42.9%) and Louvain (40.8%) make it WORSE than MCL's 15.6%.
      The giant module is what modularity maximisation wants, not an MCL artefact.
- [x] purple gamma scout + ladder done (8 gammas, leiden_mod, louvain)

## G. Scoring, all methods on the same graph (`48_cluster_compare.sh`)

`eff` is the primary criterion. `mf` is reported but never used alone — the
trivial one-cluster partition scores `mf = 1.0`.

- [x] `48_cluster_compare.sh` written (manifest-driven, one `clm info` over all methods)
- [x] `cluster_methods_sugarcane.tsv` — 27 partitions, all on the unpruned graph
      **`eff` DOES have an interior peak for Leiden CPM: gamma = 0.1 (0.550).**
      Best modularity: leiden_mod 0.1399 and louvain 0.1386, both by building a
      41-43% giant module. Best MCL modularity: 0.0818 at -I 1.7.
- [x] `cluster_methods_purple.tsv` — 19 partitions. **eff interior peak at gamma = 0.2
      (0.5611)**, higher than any MCL cell (best 0.4902 at -I 6). leiden_mod/louvain reach
      mod 0.2121 with only 55/50 clusters and 27% giants.
- [x] one-cluster baseline in the table: eff 0.0106, **mf 1.0000**, af 1.0000, mod -0.0000

## H. Independent biological check (`37_cluster_homogeneity.r`)

PFAM Sørensen–Dice homogeneity above a size-matched null. A setting has to win
on the graph metric *and* here.

- [x] `37` reads `partitions.tsv` and carries method/param/resource columns
- [x] sugarcane scored (26 partitions)
      **At MATCHED module counts MCL beats Leiden CPM on PFAM homogeneity**, not the
      reverse: ~19.5k modules +0.1051 (MCL -I 4) vs +0.0680 (CPM 0.01); ~31.5k
      +0.1187 (-I 9) vs +0.1009 (CPM 0.05). Leiden only scores higher by fragmenting
      past where MCL can reach, and its apparent peak at gamma 0.7 (+0.2602) is hollow —
      only 3,982 of 90,070 modules are scorable, falling to 1,174 at gamma 0.9.
      Robust result: leiden_mod (+0.0471), louvain (+0.0448) and -I 1.2 (+0.0433) are the
      three WORST, and they are exactly the three partitions with giant modules.
- [x] purple scored. **MCL beats Leiden CPM at matched granularity here too**
      (~36k modules: +0.0251 for -I 4 vs +0.0223 for gamma 0.2; ~22k: +0.0211 vs +0.0166).
      louvain +0.0002 and leiden_mod +0.0004 — NO annotation signal above a size-matched
      null, while scoring the highest modularity of anything tested.

## I. Hierarchical scout — sugarcane only (`50_fastgreedy_scout.py`)

Classical hierarchical is not attempted (needs a dense 232 GB distance matrix
built from correlations we thresholded away). fast-greedy (CNM) is the
graph-native stand-in, under a 4 h cap.

- [x] `50_fastgreedy_scout.py` written (+ 2 fixes: disconnected-dendrogram floor, and the
      parent hanging on an empty queue when the worker dies)
- [x] sugarcane: **cap hit at 4 h** and recorded. But an earlier run built the dendrogram
      in under 2 h before raising on the cut, so cost is contention-dependent and mixed.
      The decisive facts are structural: sugarcane has **958 connected components**, so no
      cut below 958 exists while CNM's `optimal_count` is **338** — its preferred partition
      is unreachable — and CNM optimises modularity, already ruled out. Purple not attempted.

## J. Report

- [x] `49_fig_cluster_methods.r` → `figure12_cluster_methods` (3 panels; underflowed cells excluded)
- [x] `docs/results.md`: full section written
- [x] `docs/decisions.md`: entry written
- [x] commit

---

## Log

- 2026-09-09 — branch created, config knobs added.
- 2026-09-09 — sugarcane matrix built and both assertions passed; purple building.
- 2026-09-09 17:15 — correctness gate PASSED with a negative control. Leiden and MCL
  can now be compared on the same numbers.
- 2026-09-09 17:27 — purple matrix built and verified (675,955,918 edges, exact).
- 2026-09-10 05:15 — all runs complete. fast-greedy hit its 4 h cap; the hierarchical
  verdict rests on the component floor (958) and the objective, not on timing.
- 2026-09-10 00:35 — my own dedup bug: `if(!(k in ord))` where `ord` is indexed by position,
  so nothing was ever deduplicated; and the key used $16 (runtime_s) instead of $17 (resource),
  so a -S probe row silently overwrote the default row. Both fixed, tables rebuilt.
- 2026-09-09 23:35 — CORRECTION to my own earlier reading: I said Leiden CPM beat MCL on
  homogeneity at matched granularity. It does not; MCL wins at every matched module count.
- 2026-09-09 20:10 — MAJOR: `clm info`'s eff and mf depend on WHICH OTHER clusterings
  share the call (I6 alone -> eff 0.47281; I6 in a batch of 8 -> eff 0.37821). mod and af
  are unaffected. Batched values are not monotone and look like real structure. All three
  scorers (36, 45, 48) now call clm info once per partition. The k-NN section of
  docs/results.md carries a correction: its mass-fraction column is affected, its
  modularity column — which carries the conclusion — is not.
- 2026-09-09 19:20 — I killed the purple ladder by accident: an awk pattern matched BOTH
  running 36_mcl_sweep.sh processes. Idempotence limited the loss to ~9 min.
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
