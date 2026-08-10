# Muñoz Module-20 across both networks — H1

Muñoz-Perez et al. build their nitrogen claim on **Module 20**: 12 genes, ~75%
MYB/MYB-related TAPs, high betweenness in *their* network. This pipeline asks
whether those genes keep the same characteristics — **MYB identity, high degree,
nitrogen responsiveness** — in **our** two reference-based networks.

```bash
cd /dados04/jorge/comparative_saccharum/scripts/module20
./run_all.sh          # ./run_all.sh 03 resumes mid-way
```

Outputs land in `GET_TFS/new/results/module20/`.

---

## The circularity trap, and what is actually testable

Module 20 was **defined** as MYB-rich and high-betweenness in Muñoz's own
network. "Is it MYB-rich and central?" therefore cannot be a test — it restates
the definition. The non-circular questions this pipeline asks are:

1. Do the members map to our references at all?
2. Do our **own** TF calls independently agree they are MYB?
3. Are they central **in our networks**, measured against the **MYB/MYB-related
   background** — not against all nodes, which would only re-measure "MYBs are
   well connected"?
4. Are they nitrogen-responsive in **our own re-quantification**?
5. **Does any of it survive in purple?** — the actual conservation test.

## Method

Muñoz's proteome (`raw_sugarcane/transcriptome_munoz/sugar_cane.pep`) is
**byte-identical to `GET_TFS/db/sugar_cane.pep`** — the proteome the original
GET_TFS run used. So their published `Family` column comes from the same
PlnTFDB rule method this project uses, and mapping can be done
**protein-vs-protein**.

The 12 members are TransDecoder ORFs (`<transcript>.pN`). Mapping is
`diamond blastp` into each reference proteome, with specificity from a
**reciprocal-Arabidopsis filter** — member and candidate gene must share the
same best Arabidopsis hit. Without it, a raw best hit against a family as large
as R2R3-MYB returns "some MYB", not the right one. Mapping is done
**independently** into each proteome, never sugarcane→purple.

Each member has 2–4 ORFs, some antisense/internal artefacts, and "longest" and
"complete" disagree for at least one (`SCA2_2__c65732f1p01355`: `.p1` = 341 aa
but 5prime_partial, `.p2` = 181 aa complete). So **all** ORFs are searched and
the best-scoring one is attributed to the member and recorded, rather than
imposing a fixed rule.

Earlier revisions searched the nucleotide contigs with `blastx`; those maps are
kept in `work/*_blastx.tsv` and step 02 reports the difference. **blastp is a
strict subset of blastx** (sugarcane 33 of 34 genes, purple 27 of 36) — the
protein bridge is more conservative, as expected once spurious frame hits go.

## Their classification reproduces exactly (step 02b)

Running our own pipeline (`hmmsearch --cut_ga` vs `TF.db.hmm` →
`assign_family_membership.pl` + `RulesFull`) on **their** proteins reproduces
**12 / 12** of their calls, including all **9 / 9** MYB-family members and all
3 `no_tf`. So Muñoz's "75% MYB" is correct at the source, and any disagreement
downstream is about the **orthologs we map to**, not about how TFs are called.

## The independence problem — read this before any p-value

- **2 of 12 members have no predicted ORF at all**
  (`N_R___TRINITY_DN20_c0_g1_i6`, `N_R___TRINITY_DN331_c0_g1_i3`) — non-coding /
  UTR fragments. That is why they never mapped; it is a result, not a failure.
- `SCA3__SP803280_c111248_g1_i2` / `_i4` are **isoforms of one locus** — and they
  anchor to *different* Arabidopsis genes, so that locus is likely chimeric.
- The 9 MYB members resolve to only **two MYB orthologue groups**:
  **AT5G59780 (AtMYB59)** ×5 and **AT3G46130 (AtMYB48)** ×3.

| Arabidopsis anchor | is | members |
|---|---|---|
| `AT5G59780` | **MYB59** | 5 |
| `AT3G46130` | **MYB48** | 3 |
| `AT3G55960` | uncharacterised | 1 |
| `AT4G10770` | **OPT7**, oligopeptide transporter | 1 (`no_tf`) |

So 12 published members are effectively **3–4 independent units**. Every test is
reported at gene, locus and anchor level; the **anchor level is the one to
interpret**. This is not pedantry — deduplicating haplotype copies moved the
sugarcane centrality result from p = 2.8e-05 to p = 0.118.

## Scripts

| Script | Does |
|---|---|
| `config.sh` | paths, thresholds; sourced by every step |
| `01_extract_module20.sh` | pull the member ORFs out of the Muñoz proteome |
| `02_map_to_references.sh` | `blastp` into both proteomes + reciprocal-Arabidopsis filter |
| `02b_tf_call_on_munoz.sh` | run OUR TF identification on THEIR proteins and compare |
| `03_network_readout.r` | degree percentile, our TF call, conserved edges, our modules, tests, figure |
| `04_nitrogen_response.r` | design-aware N tests in both datasets (tests only; figures live in 05) |
| `05_module20_heatmaps.r` | full per-species ComplexHeatmap figures (scico palettes) |
| `run_all.sh` | runs 01→05 |

---

## Result — the MYB component barely transfers; the N response does

**Mapping.** 10 of 12 members map in sugarcane (**33 genes, 30 in network**),
only **6** in purple (**27 genes, 25 in network**).

**What actually maps is mostly not MYB.** The members that map in volume are the
two *non-MYB* anchors; the MYB48 group maps to **nothing** that is in either
network:

| in-network genes by anchor | sugarcane | purple |
|---|---|---|
| `AT3G55960` (uncharacterised) | 11 | 12 |
| `AT4G10770` (OPT7, `no_tf`) | 10 | 11 |
| `AT5G59780` (**MYB59**) | 9 | **2** |
| `AT3G46130` (**MYB48**) | 0 | 0 |

Consequently only **30%** (sugarcane) and **8%** (purple) of the mapped
in-network genes carry a MYB call from our pipeline — even though their own
classification reproduces perfectly (step 02b). The MYB layer of Module 20 is
thin-to-absent in our reference-based networks, and in purple it is 2 genes.

**Centrality — not retained in either network.**

| | sugarcane | purple |
|---|---|---|
| median degree percentile | 63.0 | 44.1 |
| MYB/MYB-related background | 45.4 | 49.6 |
| p (gene / locus / **anchor**) | 0.118 / 0.309 / **0.444** | 0.897 / 0.716 / **0.716** |

Above the MYB background in sugarcane but **not significant at any honest level**;
in purple they sit slightly *below* it.

**Conserved edges — no enrichment** against a degree-matched null:
sugarcane 10 observed vs 12.04 expected (p = 0.85); purple 5 vs 6.01 (p = 0.77).

**Module recovery — partial in sugarcane, absent in purple.** Members scatter
over 16 sugarcane / 18 purple modules. Exactly one module in each holds ≥2
different Module-20 loci, but they are not comparable: sugarcane's is
**`Module_016` (474 genes) holding all 4 loci** — a real partial recovery
(16 vs 24.6 expected, p = 0.0003); purple's is `Module_001` (**44,050 genes**,
the giant module), and its module count is at chance (18 vs 17.0, p = 0.74).

**Nitrogen responsiveness — this is what transfers.**

| sugarcane (Muñoz's own data) | responsive RB975375 | non-responsive RB937570 |
|---|---|---|
| genes padj<0.05 (High vs Low N) | **24 / 33** | 19 / 33 |
| median \|log2FC\| | 0.85 | 0.63 |

26 of 33 genes are expressed above 1 TPM and all four loci contribute, so this
is solid — it **validates Muñoz's core claim** in our own re-quantification. The
genotype × N interaction is significant for only 3 of 33 genes (4 with leaf
segment as covariate), so the response is largely **shared between genotypes**
rather than specific to the responsive one.

**In purple, nothing survives BH correction** — 0 of 27 genes in either
genotype, by ANOVA or U-shape contrast. An earlier revision reported a
non-monotonic 51NG3 signal here; that came from genes the stricter protein-level
mapping now rejects, so it does not stand.


### Nitrogen or genotype? (step 05 variance partition)

The heatmap columns are grouped by **nitrogen only**, with genotype as an
annotation bar, so the question is visual: if the set tracks nitrogen the blocks
differ; if it tracks genotype the samples split by the Genotype bar *within*
each block. Step 05 also makes it quantitative — per gene,
`lm(expression ~ genotype + nitrogen)` and the share of the sum of squares:

| median % variance explained | nitrogen | genotype | genes with N > genotype |
|---|---|---|---|
| **sugarcane** (Muñoz) | **20.0%** | 12.1% | 18 / 33 |
| **purple** (Kiet) | **2.1%** | **40.4%** | 7 / 25 |

**In sugarcane Module 20 tracks nitrogen; in purple the same genes track
genotype almost exclusively** — a ~19x swing in the nitrogen share. In the
purple figure the columns cluster cleanly by genotype inside every nitrogen
block; in sugarcane the low-N block is uniformly warm and the high-N block
uniformly cool across all four loci, with genotypes interleaved. This is the
cleanest single statement of why the module does not transfer: it is not that
the genes are absent, it is that in LA purple their variation is genotypic, not
nutritional. Per-gene values in `module20_variance_partition_<sp>.tsv`.

**Reading.** Module 20's *nitrogen response* reproduces robustly in Muñoz's own
data, but its *MYB identity* barely reaches our references (2 MYB genes in
purple), its *centrality* is indistinguishable from ordinary MYBs in both
networks, and its *module structure* survives only partially in sugarcane and
not at all in purple. Same split the project finds everywhere: **the response is
conserved, the topology is not.**
