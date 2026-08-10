# TF (TAP) identification — H1, applied to our two reference proteomes

Identify transcription factors / transcriptional regulators (TAPs) in the two
reference proteomes behind our co-expression networks, using the
**PlnTFDB rule-based method** (Riaño-Pachón 2007; Pérez-Rodríguez 2010): scan
each proteome with a curated HMM library, then assign proteins to TF families
with a required/forbidden-domain rule set.

This re-runs the original `GET_TFS` pipeline on **our data**
(R570 = *sugarcane*, LA purple = *purple*), and modernises the driver: the old
`#$ -t 1-84` SGE array job (`code/run_tf.sh`) is replaced by a portable,
**resumable, locally-parallel** driver — no scheduler, no `module`, no
BioPython / GNU-parallel dependency.

> **Scripts location:** the 7 pipeline scripts live in
> **`../../scripts/get_tfs/`** (`config.sh`, `00_…`–`04_…`, `run_all.sh`),
> grouped with the rest of the project's `scripts/`. This folder keeps the
> **results** (`results/<species>/`) and this method write-up.

---

## Pipeline (5 steps)

```
proteome.fa
   │  01_split_proteome.sh     split into N_CHUNKS pieces (awk)
   ▼
chunks/chunk_*.fa
   │  02_run_hmmsearch.sh      hmmsearch --cut_ga --domtblout, MAX_JOBS in parallel, concat
   ▼
concatenated.domtbl
   │  03_assign_families.sh    assign_family_membership.pl + RulesFull
   ▼
family_assignment.tsv          protein → Family → Type (TFF / OTR / Orphans)
   │  04_postprocess.sh        drop Orphans, collapse to gene, ∩ network nodes
   ▼
TF_in_network.ids / .tsv       ← H1 input: TFs that are nodes of the network
```

Run one species end-to-end:

```bash
cd /dados04/jorge/comparative_saccharum/scripts/get_tfs
./run_all.sh sugarcane
./run_all.sh purple
```

Or step-by-step (same argument each time), e.g. `./01_split_proteome.sh purple`.
Tune parallelism via env vars (defaults suit this 256-core box):

```bash
N_CHUNKS=64 MAX_JOBS=32 CPU_PER_JOB=4 ./run_all.sh sugarcane   # 32×4 = 128 cores
```

`02_run_hmmsearch.sh` is **resumable**: finished chunks are skipped, partial
writes go to `.tmp`, so re-running continues where it stopped.

---

## ✅ Checklist — software

Run `./00_check_requirements.sh <species>` to verify all of this automatically.

| Software | Needed for | This box |
|---|---|---|
| **HMMER ≥ 3.1** (`hmmsearch`) | domain search | ✅ `/usr/bin/hmmsearch` (3.4) |
| **perl** | `assign_family_membership.pl` | ✅ `/usr/bin/perl` |
| **awk, sort, coreutils** | split + post-processing | ✅ |
| bash ≥ 4.3 | `wait -n` job throttling | ✅ |

No conda env required — all system tools. (The old `code/run_tf.sh` used
`module load Hmmer/3.2.1` + SGE; neither exists here, hence the rewrite.)

## ✅ Checklist — shared pipeline files (reused, already present)

| File | What it is |
|---|---|
| `../db/TF.db.hmm` | HMM library, 19,649 Pfam-A + PlnTFDB models, **with GA cutoffs** |
| `../mytfdb/RulesFull` | family rules: required / forbidden domains per family |
| `../mytfdb/assign_family_membership.pl` | rule engine (Riaño-Pachón), reads a hmmsearch domtblout |

## ✅ Checklist — our input data (per species)

| Species | Proteome | Network nodes |
|---|---|---|
| **sugarcane** (R570) | `files/fix_orthofinder/sugarcane/SofficinarumxspontaneumR570_771_v2.1.protein.fa` (299,731 prot) | `files/sugarcane/network_sugarcane_node_metrics.tsv` |
| **purple** (LA purple) | `files/fix_orthofinder/purple/one_transcript_purple_proteins.faa` (241,263 prot, 1/gene) | `files/purple/new/network_purple_node_metrics.tsv` |

These are the **same proteomes used for the OrthoFinder bridge**, so the TF gene
ids line up directly with the orthogroup mapping used elsewhere in the project.

### ID reconciliation (proteome → gene → network node)

| Species | Protein id | → gene | Network node id | → gene |
|---|---|---|---|---|
| sugarcane | `SoffiXsponR570.02Eg130400.1.p` | strip `.<iso>.p` | `SoffiXsponR570.02Eg130400.v2.1` | strip `.v2.1` |
| purple | `Soffic.02F0006060-2F` | identity | `Soffic.02F0006060-2F` | identity |

(Encoded as `PROT_TO_GENE` / `NODE_TO_GENE` in `config.sh`.)

---

## Outputs (`results/<species>/`)

| File | Contents |
|---|---|
| `concatenated.domtbl` | merged hmmsearch per-domain hits |
| `family_assignment.tsv` | protein → Family → Type (incl. Orphans) |
| `TF_no_Orphans.tsv` | protein-level TAP calls, Orphans removed |
| `TF_genes.ids` | all proteome genes with a TF call |
| `network_genes.ids` | genes present as network nodes |
| **`TF_in_network.ids`** | **TF genes that are network nodes — H1 input** |
| **`TF_in_network.tsv`** | gene · Family · Type (network TFs only) |
| `family_counts.tsv` | per-family counts (all TF genes vs in-network) |
| `summary.txt` | headline numbers |

---

## What changed vs the original `GET_TFS`

- **Parallelism**: SGE array (`run_tf.sh`, fixed `1-84`, `module load`) →
  portable background-job driver with a concurrency cap and `-Z` set to the
  full proteome size (comparable E-values across chunks; `--cut_ga` is the real
  filter so scoring is unaffected).
- **Resumability**: chunk skip + atomic `.tmp`→final rename.
- **No BioPython**: `split_fasta.py` → one awk one-liner.
- **Config-driven**: one `config.sh` parameterises both species (paths + id
  regexes); every step takes just `<species>`.
- **Network intersection built in**: original produced ad-hoc
  `TF_in_network_nosuffix.ids`; here it is a defined, reproducible step with a
  documented id mapping.

## Notes

- `--cut_ga` uses each model's curated Gathering threshold — the recommended,
  reproducible cutoff (no arbitrary E-value). `-Z` only rescales reported
  E-values and does **not** change which hits pass `--cut_ga`.
- Sugarcane uses **all isoforms**; calls are collapsed to gene in step 04 (a
  gene is a TF if any isoform carries the family's domain signature). Purple
  already has one transcript per gene.
- **Orphans** (proteins with a TF-associated domain but not satisfying any full
  family rule) are dropped from the final TF set, per the method's README.
