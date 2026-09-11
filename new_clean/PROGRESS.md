# Module analysis on the Pearson-only graphs — progress

Branch `clustering-methods`. **Overwrites `results/`** — this is the new main
analysis. Backup of everything displaced: `results/_pre_pearson_20260910/` (38 MB,
81 files, verified).

Graphs: unpruned Pearson-only, `|r| >= 0.8`, weights `[0.01, 1]`, no k-NN.
Inflation: **sugarcane 1.5**, **purple 3.5** — each species' modularity optimum
(figure 12). The module-trait test gains the experimental design.

The previous phase's checklist is kept at `docs/PROGRESS_clustering_phase.md`.

Legend: `[ ]` todo · `[~]` running · `[x]` done · `[!]` blocked/failed

---

## A. Safety

- [x] back up `results/{mcl_*,modules/,node_metrics,module_trait_*}` + figures
      → `results/_pre_pearson_20260910/`, 81 files, verified
      (keeps the OLD marginal `module_trait_*.tsv` for the blocked-vs-marginal comparison)

## B. Inputs the new analysis needs

- [x] `52_pearson_node_metrics.sh` written
- [x] sugarcane node metrics: 101,990 genes
- [x] purple node metrics: 170,135 genes

## C. The clustering contract, from the partitions already on disk

No re-clustering: `cls.knone.I15` and `cls.knone.I35` ARE the chosen settings.

- [x] `51_mcl_membership_from_cls.r` written (mirrors `27_sbm_membership.r`)
- [x] `.inflation` guard written AND tested firing (adopting I15 as "-I 15" aborts, writes nothing)
- [x] sugarcane: 5,653 modules, largest 23,439 (22.98%), 3 unassigned, Q=0.08212 — cross-checked
- [x] purple: 29,624 modules, largest 39,230 (23.06%), 15,600 unassigned, Q=0.15324 — cross-checked
- [x] every gene accounted for: 101,990 and 170,135, exactly

## D. Eigengenes

- [x] sugarcane: 3,627 eigengenes, PC1 median 80.8% of variance
- [x] purple: 7,493 eigengenes, PC1 median 88.4%

## E. The blocked module-trait test

```
sugarcane   eigengene ~ genotype + segment + N     n = 48, resid df 42
purple      eigengene ~ genotype + N               n = 18, resid df 15
```

- [x] `fit_blocked` + `verify_against_lm` + `row_midranks` in `scripts/lib/common.R`
- [x] `53_module_trait_blocked.r` written
- [x] **verified against `lm()`**: max |dt| 4.4e-16 / 1.3e-15, max |dp| 2.2e-16 — both species
- [x] model matrix rank-checked; unblocked Spearman reproduces `cor.test(exact=FALSE)`
- [x] Spearman primary (midranks); blocked Pearson beside it
- [x] marginal computed ON THESE EIGENGENES (see log — joining the old table was wrong)
- [x] sugarcane plant control: 12 plants, df 9, agreement with blocked fit **r = +0.9486**
- [x] within-block null, 1,000 draws: **0 permutations reached the observed count** in either species (empirical p <= 0.001)
- [x] sugarcane: **588 responsive** of 3,627 (16.2%), null mean 0.14
- [x] purple: **182 responsive** of 7,493 (2.4%), null mean 0.39
- [x] **blocked is a strict superset in both**: sugarcane 251 -> 588, purple 96 -> 182, none lost

## F. The rest of the chain

- [x] moduleprofile — TF enrichment: sugarcane OR 2.41 p 0.00207; purple OR 3.36 p 0.274 (1 module)
- [x] moduleheatmap — 40 per study
- [x] modulesummary — 250 x 48 (sugarcane), 182 x 18 (purple)
- [x] figmodules — figure5_modules + legend + stats

## G. GO — handed over, not run

topGO runs are the user's to launch in `topGO_env`.

- [x] inputs prepared and all four checked present per study; universe is now the Pearson-only network
- [~] exact command handed over (below)

## H. Report

- [x] `results/STALE.md` written
- [x] `docs/results.md` — new section added AND the old module section marked superseded in place
- [x] `docs/decisions.md` — entry written
- [x] confirmed: exactly 8 files replaced, merged edge tables untouched
- [ ] commit

---

## Log

- 2026-09-10 — backup taken; inflation confirmed as sugarcane 1.5 / purple 3.5
  (the values in the request were swapped relative to figure 12).
- 2026-09-10 — CAUGHT: my first version of 53 carried the marginal rho by JOINING the
  previous run's `module_trait_*.tsv` on `module`. Module names are positional
  (Module_%03d, largest first), so that joined DIFFERENT GENE SETS across two
  clusterings — old Module_001 held 19,604 genes, the new one holds 23,439. It ran
  cleanly and reported "244 marginal calls lost". The marginal fit is now computed on
  the same eigengenes, and blocked is a strict superset in both species.
