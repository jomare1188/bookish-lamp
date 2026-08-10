# Locating MYB61 in both references — H1

Kiet et al. (2025) build their nitrogen story around **ScMYB61.1**, published as
gene id `Soff.09G0002230–3D`. This pipeline finds that gene in **our** two
reference proteomes, because the published id cannot be used directly.

```bash
cd /dados04/jorge/comparative_saccharum/scripts/myb61
./run_all.sh          # everything;  ./run_all.sh 03  resumes from step 03
```

Outputs land in `GET_TFS/new/results/myb61/`.

---

## Why not the published id

Kiet et al. give no sequences, and the id does not resolve:

1. **It is structurally impossible in the LA purple annotation.** In all
   **241,263** gene ids the haplotype suffix number equals the chromosome
   number (`Soffic.09G0003220-**9**G`, `Soffic.03D0011800-**3**D`) — 241,263
   matches, **0 exceptions**. `09G0002230–3D` pairs a chr9 locus with a chr3
   suffix, so it is a chimeric/mistyped id.
2. **Neither single-field repair lands on a MYB.**
   - keep chr9 → `Soffic.09G0002230-9{A,B,C,D,G}`: no Myb domain even at a
     permissive `E<10`; the only Pfam hits are `Transposase_21`,
     `Transpos_assoc` and `DUF1218`. These gene models are TE-contaminated.
   - keep the suffix → `Soffic.03G0002230-3D`: 230 aa, no Myb signal, no
     Arabidopsis hit at all.
3. **The orthogroup bridge agrees** — copies of locus `09G0002230` never share
   an orthogroup with any confirmed MYB61 copy (step 05).

The two studies also classified TFs differently: Kiet used **PlantTFDB**
(homology to a curated database), while this project uses **Pfam GA + PlnTFDB
rules** (`scripts/get_tfs/`). Those are not directly comparable, which is a
second reason to re-derive MYB61 from sequence rather than trust a label.

---

## Method

```
AtMYB61 (AT1G09540)
   │ 02  reciprocal best hit -> sorghum + rice
   ▼
5 anchors (At, 2x Sb, 2x Os)          all verified: 2 Myb repeats under --cut_ga
   │ 03  DIAMOND blastp -> R570 + LA purple proteomes  (E<=1e-10, qcov>=50%)
   ▼
21 sugarcane / 29 purple candidate genes
   │ 04  filter 1: reciprocal best hit back to Arabidopsis must return AT1G09540
   │     filter 2: independently called MYB by our own GET_TFS run
   ▼
11 sugarcane / 16 purple confirmed
   │ 05  OrthoFinder bridge — do both species land in the same orthogroups?
   │ 06  MAFFT + FastTree with decoys — do they form a clade?
   │     -> recovers 1 extra sugarcane copy the domain rules had downgraded
   ▼
12 sugarcane / 16 purple MYB61 copies
   │ 07  network readout: degree percentile, N correlation, conserved edges
```

**Why RBH and not a top hit.** R2R3-MYB is one of the largest plant TF families
(~125 genes in Arabidopsis). A forward search alone returns "some MYB". The
reciprocal test — best Arabidopsis hit must be AT1G09540 again — is what makes
the call specific. The two filters use independent evidence (homology vs domain
rules), so their agreement is meaningful and their **disagreement is reported,
not hidden**.

**Why sorghum and rice anchors.** Sorghum is the closest well-annotated relative
of *Saccharum* (same tribe, Andropogoneae). Going straight from Arabidopsis
would cross ~150 My in one jump. Both grasses return **two** reciprocal best
hits — the expected co-orthologue pair left by the grass whole-genome
duplication, and both are recovered in our species (chr3-type and chr9-type).

---

## Scripts

| Script | Does |
|---|---|
| `config.sh` | all paths, thresholds, the `AT1G09540` anchor; sourced by every step |
| `01_fetch_references.sh` | download Arabidopsis / sorghum / rice proteomes (Ensembl Plants), reduce to one protein per gene |
| `02_build_query.sh` | RBH from AtMYB61 into sorghum + rice; verify Myb repeats with our own HMM library |
| `03_search_proteomes.sh` | DIAMOND blastp anchors → R570 + LA purple; collapse proteins to genes |
| `04_rbh_filter.sh` | reciprocal best hit + our own MYB call; emit per-gene verdicts |
| `05_orthogroup_check.sh` | place copies on the OrthoFinder bridge; cross-check the published locus |
| `06_phylogeny.sh` | MAFFT + FastTree with decoys; clade test; rescue broken gene models |
| `07_network_readout.r` | network position, N correlation, conserved-edge enrichment, figure |
| `08_myb61_expression_test.r` | per-genotype ANOVA + U-shape contrast on Kiet's 18 leaf libraries; expression heatmap |
| `run_all.sh` | runs 01→08 |

Software: `diamond` 2.1.9, `blastp` 2.12.0+, `hmmsearch`/`hmmfetch` (HMMER 3.4),
`mafft`, `FastTree`, `curl`, awk/sort — all system tools. R: `ape` from the
`comparative_network` env (step 06), `ggplot2`/`data.table`/`patchwork` from
`r_env` (step 07); both paths are set in `config.sh`.

### Two data quirks handled in code

- The **LA purple proteome contains `.` characters inside 542 sequence lines**
  (annotation gap marks). HMMER tolerates them, DIAMOND rejects the file. Step 03
  writes a sanitised copy (`*` dropped, non-standard residues → `X`); ids and
  lengths are unchanged.
- **`Orthogroups.tsv` has CRLF line endings.** Fields are stripped of `\r`
  before matching — without this every sugarcane lookup silently fails.

---

## Result

**MYB61 is present, orthologous, and shared** — 12 copies in R570, 16 in LA
purple, in two co-orthologue clades, joined by **5 shared orthogroups**. All 27
strictly-confirmed copies fall inside a single 32-tip grass MYB61 clade in the
tree, with all four monocot anchors and no true intruder.

Caveats stated rather than buried:

- **OrthoFinder leaves many copies unassigned** (6/12 sugarcane, 3/16 purple are
  singletons). This is the known polyploid haplotype-fragmentation problem, so
  absence from `Orthogroups.tsv` is *not* evidence of absence; every copy's
  status is in `MYB61_bridge_status.tsv`.
- `SoffiXsponR570.02Ag129700` is a genuine MYB61 copy whose gene model retains
  only one Myb repeat, so the strict rules called it `MYB-related`. The tree
  recovers it; it is flagged `MYB61_homolog_no_MYB_call` and counted separately.

Network readout (`MYB61_network_table.tsv`, figure `MYB61_network_overview.png`):

| | sugarcane | purple |
|---|---|---|
| copies / in network | 12 / 8 | 16 / 15 |
| median degree percentile | 45.7 | 45.1 |
| copies N-correlated (\|r\|≥0.6, padj≤0.05) | **0** | **0** |
| strongest \|r\| with nitrogen | 0.247 | 0.643 (padj 0.27, n.s.) |
| on conserved edges | 37.5% vs 37.2% bg (p=0.62) | **60% vs 24.9% bg, 2.4x, p=0.0041** |
| ...vs degree-matched null | 3 vs 2.92 exp (p=0.62) | **9 vs 3.41 exp, 2.6x, p=0.0019** |

`on conserved edges` is a **binary gene-level flag** — the gene has >=1 conserved
edge, not a fraction of its edges — and it is strongly degree-dependent (purple:
3.2% of nodes in the lowest degree decile, 43.3% in the highest). Step 07
therefore tests it twice: against all nodes, and against a **degree-matched null**
(20,000 resamples). The purple enrichment survives both.

So in our data MYB61 is **not a hub** (it sits at the median of both degree
distributions) and **not significantly nitrogen-correlated in either dataset** —
but in the purple network it *is* significantly enriched on cross-species
conserved edges. Its role here looks structural/conserved rather than
N-responsive.

## Testing Kiet's claim directly (step 08)

The project's standard `gene_trait_cor.r` uses a **linear Pearson against
treatment coded 0/2/6, pooled over both genotypes**. Kiet describe ScMYB61.1 as
**non-monotonic** (highest at both 0N and 6N, lowest at 2N) and **opposite
between genotypes** — a shape Pearson has ~zero power against, and a split that
pooling cancels. Step 08 therefore re-tests on the same 18 leaf libraries with
per-genotype ANOVA, a U-shape contrast `c(+1,-2,+1)`, a linear contrast, and a
genotype x N interaction.

- A **significant non-monotonic response does exist** in 51NG3: U-contrast
  padj<0.05 for **10 of 13** testable copies; ANOVA padj<0.05 for 7.
- It is **inverted** relative to the paper: 11 of 13 U estimates are negative,
  i.e. expression peaks at **normal** N and falls at both extremes.
- The **linear contrast is significant for 0 copies** — direct evidence that the
  earlier Pearson test was structurally blind to this, not that nothing existed.
- Genotype x N interaction is significant (`09F0003800-9F` padj 7.2e-05, four
  more below 0.01), so the "opposite between genotypes" claim holds in shape.
- **Caveat that dominates the rest**: these copies are barely expressed in leaf
  (mean TPM 0.0-1.5; only 2 of 16 above 1 TPM), and this is leaf-only whereas
  Kiet emphasise roots. Suggestive, not conclusive.
- Negative control: `Soffic.09G0002230-9G` at the published locus runs at
  **79.6 / 44.0 mean TPM** — one to two orders of magnitude above any true MYB61
  copy.

### The cloning primers settle the identification

The paper's forward primer is `GGATCC` (BamHI) + `ATGGGGAGGCATTCTTGC` → protein
starting **M-G-R-H-S-C**. All **28** copies found here begin `MGRHSCCYKQKL`;
**none** at the published locus do (`MARPPTMVIQDD`, `MDRSWIFGIKFT`,
`MVHPETMLLHEP`, `MARRRRLSASSL`, `MALSASLRISLI`, `MLLLIAHAAAAC`). So the gene
they cloned and overexpressed is one of ours — their Arabidopsis phenotype rests
on the real MYB61, while the RNA-seq signal attributed to `Soff.09G0002230-3D`
belongs to a different, much more highly expressed non-MYB locus.
