# dN/dS — does sequence evolution follow network conservation?

The project has shown that co-expression **edges** are conserved between the two
networks ~2.5× above a null, while the **nitrogen response** is not conserved at
all. Both are statements about expression. This stage asks whether the edge
result has a counterpart in the sequence:

> does a gene whose co-expression neighbourhood is conserved also evolve under
> stronger purifying selection?

dN/dS is measured from genomic sequence and shares no noise with the RNA-seq the
networks came from, so it can only strengthen or falsify the conservation
result — it cannot restate it.

Stage: `scripts/29_dnds/` (`./run.sh dnds`, or step by step — two steps are long).

---

## The fact that shapes the design

R570 and LA purple are **95–99% identical at the protein level**. Measured on the
test set of single-copy orthologs:

| pair | median dS | genes with dS = 0 |
|---|---|---|
| sugarcane ↔ sorghum | 0.214 | 0 |
| purple ↔ sorghum | 0.224 | 0 |
| **sugarcane ↔ LA purple** | **0.005** | **40%** |

A per-gene ω between the two *Saccharum* genomes is therefore undefined for
about two genes in five, and noise for many of the rest. A pal2nal → codeml run
on that pair alone would produce a table that looks perfectly well-formed and
means almost nothing.

So the stage does two different things with two different comparisons:

1. **Per-gene ω against an outgroup.** *Sorghum bicolor* (~8 Myr; NCBI RefSeq
   `GCF_000003195.3_Sorghum_bicolor_NCBIv3`, downloaded in full by step 00 and
   MD5-verified). dS lands around 0.21, so ω is well estimated per gene. This is
   the quantity correlated with network position.
2. **The R570 ↔ LA purple pair, only ever aggregated.** Exact NG86 counts
   (Nd, Sd, N, S) are computed per gene and **summed within bins**; one ratio is
   formed per bin. Ratios cannot be averaged when the denominator is often zero.
   Counts can be summed. This is the only honest use of a 99%-identity pair.

## Why a three-species OrthoFinder rerun, and how sorghum actually attaches

Sorghum joins the orthology by inference rather than by best-hit heuristic:
purple + sugarcane + sorghum, same flags as `files/fix_orthofinder/run_orthofinder.sh`,
into `files/orthofinder_3sp/`. OrthoFinder 3.1.3 removed the 2.x
`-b <prev> -f <new>` add-species mode (`--assign/--core` replaced it, and that is
mutually exclusive with `-f`), so an incremental add was not available. The run
took 2h17m at `-t 256` and put sorghum where it belongs:
`(sorghum,(sugarcane,purple))`.

**The two-species run is not touched.** `config.sh:ORTHOGROUPS` still points at
`Results_Jun04_2`; step 02 checksums that file before and after and fails if it
changed. It did not.

### The 1:1:1 orthogroup route does not work here, and why

The obvious unit — orthogroups that are strictly single-copy in all three species
— yields only **532** of them, against 33,982 single-copy *pairs* in the
two-species run. Adding sorghum does not just filter that set: sorghum bridges
Saccharum subfamilies that MCL had kept apart, so the clustering coarsens
(50,494 orthogroups, mean size 7.3, against 94,229 before) and almost nothing
survives as 1:1:1. 532 genes cannot carry a regression with five covariates, ten
degree deciles and five conservation bins.

> **OrthoFinder 3.1.3 caveat, found here:** its
> `Orthogroups_SingleCopyOrthologues.txt` lists 532 ids — the right *count* — but
> only **185** of them are actually 1:1:1. The others have shapes like (3,0,0)
> and (2,0,1), i.e. genes in a single species. The count is computed correctly
> and the id list is not. The file is therefore never read; where a 1:1:1 set is
> needed it is derived from `Orthogroups.tsv`.

### What is used instead

The analysis unit stays the **two-species 1:1 pair** — exactly the set every
conservation number in this repo was computed from, so the dN/dS layer and the
conservation layer are talking about the same genes by construction, with no gate
applied after the fact.

Sorghum attaches to each pair through the three-species run's **tree-based
orthologue tables** (`Orthologues/`). A pair is kept when the sugarcane gene has
exactly one sorghum orthologue, the purple gene has exactly one, and they are the
same gene:

| | pairs |
|---|---|
| two-species 1:1 pairs | 33,982 |
| dropped — no sorghum orthologue on one side | 12,455 (36.7%) |
| dropped — more than one on a side | 890 (2.6%) |
| dropped — the two sides name *different* sorghum genes | **21 (0.06%)** |
| **kept** | **20,616** |

That last row is a strong quality check obtained for free: where both sides
resolve to exactly one sorghum gene, they agree on **which** gene 99.9% of the
time.

## Pipeline

| step | what it does |
|---|---|
| `00_setup.sh` | installs pal2nal (`dnds_env`); downloads the sorghum genome, annotation, CDS and proteins in full, MD5-verified, with a `PROVENANCE.txt` |
| `01_prepare_cds.sh` | one CDS + one protein per gene per species, keyed by the gene id the networks use |
| `02_orthofinder_3sp.sh` | **long, ~3–4 h** — the three-species run |
| `02b_build_triplets.sh` | 1:1:1 orthogroups, gated on the two-species run |
| `03_align.sh` | MAFFT L-INS-i on the 3 proteins → pal2nal → codon alignment |
| `04_codeml.sh` | codeml `runmode=-2` (all three pairs at once) + exact NG86 counts |
| `04b_yn00_check.sh` | independent estimator on the same alignments |
| `05_conserved_degree.sh` | per-gene edge-conservation fraction, from the full edge tables |
| `06_join_and_test.r` | the join and the six tests |
| `07_plots.r` | the figure |

### The translation gate in 01

`translate(CDS)` must reproduce the protein. **pal2nal fails silently on a frame
mismatch**, so anything that does not translate cleanly is dropped here and
counted in `work/seqs/cds_qc_<species>.tsv`:

| species | kept | dropped |
|---|---|---|
| sugarcane | 194,593 / 194,593 (100.00%) | 0 |
| purple | 240,744 / 241,263 (99.78%) | 519 |
| sorghum | 28,153 / 28,236 (99.71%) | 83 |

The three id transforms differ. Sugarcane carries a `locus=` field that is
already the network gene id. Purple's one-transcript proteome does **not** always
keep `.t1`, so the right isoform is recovered by translating each candidate and
matching the protein string. Sorghum is RefSeq, joined CDS↔protein through
`[protein_id=XP_…]` and reduced to the longest protein per GeneID — which happens
*before* OrthoFinder, so all three proteomes share one convention.

### Why the alignments are rebuilt

OrthoFinder's `MultipleSequenceAlignments/` are not reused. FAMSA output that has
been through OrthoFinder's column trimming no longer maps 1:1 onto the CDS, which
is exactly the correspondence pal2nal requires. Three short sequences realign in
milliseconds.

## Estimator cross-check

codeml (maximum likelihood) and yn00 (counting) share the alignments and almost
nothing else, so agreement is evidence the plumbing is right — ids joined
correctly, frames intact, pal2nal mapping sound.

| set | Spearman ρ |
|---|---|
| all alignments | 0.943 |
| **dS ∈ [0.01, 2] — the analysed set** | **0.976** |

The disagreement is concentrated in saturated genes, which the dS filter removes
anyway. The filtered figure is the one that gates.

## The two directions' predictors are not equally informative

`frac_conserved` is computed per gene in both directions (step 05). They are
shaped very differently, and that governs how much the purple replicate can say:

| | sugarcane → purple | purple → sugarcane |
|---|---|---|
| genes | 103,336 | 170,736 |
| **frac_conserved = 0** | **62.0%** | **74.2%** |
| mean | 0.0776 | 0.0092 |
| 90th percentile | 0.3125 | 0.0254 |
| max | 1.0000 | 0.5588 |

Purple's is compressed for the same reason its raw conservation rate is 1.52%
against sugarcane's 10.76%: purple has 9.3× more edges, so a far smaller share of
them find a counterpart. The consequence is statistical, not biological — with
74% of purple genes tied at exactly zero, **a Spearman correlation on purple's
predictor is weak by construction**, and the binary contrast (≥1 conserved edge
vs none) is the form to read for that direction. Sugarcane's predictor has real
dynamic range across its upper 38% and is the one the graded tests are for.

Treat the purple direction as a sign-agreement check on sugarcane's result, not
as an equally powered independent replicate.

## How to read the tests

They are run in an order that makes the answer interpretable, and **test 1 is not
the result**.

1. **Marginal** — Spearman ρ between ω(sugarcane, sorghum) and the fraction of a
   gene's edges that are conserved; plus the same against purple's own
   conservation as an independent replicate.
2. **Partial — this is the result.** Degree, expression level and conservation are
   all mutually correlated, and expression level is the strongest known predictor
   of ω in every organism it has been measured in. The model is
   `log(ω) ~ frac_conserved + log(degree) + mean_log_tpm + aln_codons + GC3`,
   and the reported quantity is the **partial** coefficient on `frac_conserved`
   with its CI, plus how much adjusted R² conservation adds over a model that
   already contains the confounders. If the effect collapses once expression and
   degree are in, *that* is the finding and must be stated that way.
3. **Degree-stratified** — the same test within each degree decile. A consistent
   sign across ten bins is far stronger evidence than one coefficient.
4. **Hubs** — the positive control. Hubs are expected to be more constrained; if
   that does not reproduce, the pipeline is suspect, not the biology.
5. **Lineage-specific** — the R570 ↔ LA purple pair, NG86 counts summed within
   conservation bins, bootstrap CI over genes.
6. **Effect size over p-value.** At n ≈ 25k everything is significant. Report ρ,
   the coefficient, and the ω difference between the top and bottom conservation
   decile, and lead with those.

## Results

20,616 triplets → 20,376 codon alignments → **12,334 genes** in the analysis set
(the drop is mostly genes absent from the sugarcane network: 20,376 → 12,423).
codeml vs yn00 on the analysed set: ρ = **0.992**.

### The answer: no, at least not at the depth this can see

| test | result |
|---|---|
| 1 marginal, sugarcane | ρ = **0.016**, p = 0.083 |
| 1 binary, ≥1 conserved edge vs none | p = **0.99** (median ω 0.1841 vs 0.1808) |
| 1 replicate, purple | ρ = **−0.029**, p = 6.6e-05 — opposite sign, effect ≈ 0 |
| **2 partial coefficient on conservation** | **−0.108**, CI [−0.183, −0.033], p = 0.0047 |
| **2 adjusted R² conservation adds** | **0.0004** |
| 3 degree deciles with negative ρ | **3 of 10** |
| 6 median ω, top ÷ bottom conservation decile | **1.010** |

The partial coefficient is in the hypothesised direction and is significant, and
it means almost nothing: conservation adds **0.04%** of explained variance to a
model that already has expression, degree, length and GC3. For scale, in the same
model expression carries p = 1.8e-123 and GC3 p ≈ 0, and the whole model reaches
adjusted R² = 0.244.

### What actually drives omega here: mostly the denominator

The model in test 2 reaches adjusted R2 = 0.244, and it is worth knowing where
that comes from, because the answer is not what the literature would predict.

| predictor | R2 alone | dR2 if dropped | std. coef | omega, q10 -> q90 |
|---|---|---|---|---|
| **GC3** | **0.190** | **0.224** | **-0.543** | **x 0.33** |
| mean log2 TPM | 0.009 | 0.035 | -0.203 | x 0.68 |
| alignment length | 0.010 | 0.010 | -0.109 | — |
| log degree | 0.000 | 0.004 | -0.069 | — |
| edge conservation | 0.000 | 0.0005 | -0.024 | x 0.97 |

GC3 dominates. But it is dominating the **ratio**, not selection — it acts almost
entirely on the synonymous denominator:

| | spearman with GC3 |
|---|---|
| log dN | +0.176 |
| **log dS** | **+0.744** |
| log omega | -0.454 |

Across GC3 tertiles dN rises 0.0177 -> 0.0254 (1.4x) while dS rises 0.0726 ->
0.2198 (3x). GC-rich third positions accumulate synonymous changes faster --
mutational bias, GC-biased gene conversion, codon usage -- which inflates dS and
deflates omega. That is composition acting on the measurement, not stronger
purifying selection.

And expression, which the model makes look like the main biological driver, has
almost no marginal relationship with omega at all:

| | spearman with mean log2 TPM |
|---|---|
| log dN | -0.311 |
| log dS | -0.301 |
| **log omega** | **-0.049** |

Expression lowers dN and dS by nearly the same factor, so the ratio barely moves.
Highly expressed genes evolve more slowly overall -- a real and well-known effect,
visible here at 40% lower rates -- but that is a **rate** result, not a
selection-intensity one. Expression's apparent weight in the model is conditional
on GC3 being in it (the two correlate at -0.285).

**So the honest summary of this whole stage: nothing we measured moves omega
much.** Substitution *rates* vary strongly and predictably with expression and
base composition; the *ratio* is close to flat across expression, degree,
conservation and everything else tested. Anyone wanting to rank these genes by
selective constraint should be aware that omega here is buffered, and that dN
with an explicit rate control may be the more informative readout.

### The positive control works, so this is a real null

Hubs are more constrained — median ω **0.1774** against **0.1860**, p = 6.0e-04.
The pipeline can detect a network-position effect on ω. It simply does not find
one for edge conservation.

### The lineage-specific trend is expression, not conservation

Summing NG86 counts within conservation bins for the R570 ↔ LA purple pair gives
what looks like a clean decline — ω from **0.393** to **0.301**, with
non-overlapping bootstrap CIs. It does not survive contact with the confounders,
because the bins are not otherwise comparable:

| frac. conserved | n | median degree | **median log2 TPM** | median GC3 | ω (pair) |
|---|---|---|---|---|---|
| 0 | 5,099 | 9 | **0.45** | 0.767 | 0.393 |
| (0, 0.05] | 2,244 | 550 | **1.22** | 0.650 | 0.370 |
| (0.05, 0.1] | 975 | 166 | **1.64** | 0.590 | 0.310 |
| (0.1, 0.2] | 1,174 | 96 | **2.06** | 0.579 | 0.329 |
| > 0.2 | 2,931 | 228 | **2.16** | 0.541 | 0.301 |

Re-forming the same ratio **within expression tertiles** removes the trend:

| frac. conserved | low expr. | mid expr. | high expr. |
|---|---|---|---|
| 0 | 0.422 | 0.360 | 0.310 |
| (0, 0.05] | 0.412 | 0.373 | 0.296 |
| (0.05, 0.1] | 0.432 | 0.262 | 0.267 |
| (0.1, 0.2] | **0.485** | 0.363 | 0.261 |
| > 0.2 | 0.407 | 0.365 | 0.215 |

The low and mid columns are flat — the pooled decline does not reproduce in
them. The **high** column does decline, 0.310 → 0.215, so expression accounts for
most of the pooled trend but not all of it; a residual survives in the most
highly expressed third, which is also the third where the counts are best
measured. That residual is not clean either: GC3 falls from 0.767 to 0.541 across
the same bins and carries a coefficient of −2.15 in the per-gene model, so the
tertile split controls expression but not composition.

The per-gene analysis is the one that can hold everything at once, and it agrees:
the partial effect of conservation is real, negative, and worth 0.04% of adjusted
R². Read the two together as "there is something there, and it is far too small
to matter", not as "there is nothing".

**Conclusion: network edge conservation does not predict selective constraint, at
either timescale, once expression level is accounted for.** Panel D of the figure
is drawn to show this directly — the pooled line in grey, the stratified lines
over it.

## A second description of network position: local clustering

The analysis above tests constraint against exactly one network property,
degree. Steps 14-16 add the **local clustering coefficient** -- of all the pairs
of a gene's neighbours, what fraction are connected to each other.

### Why not betweenness or closeness

Not mainly cost. Sugarcane has mean degree **1,475** (density 1.4%) and purple
**8,265** (4.8%), so the effective diameter is **2-3 hops**. When nearly every
pair of nodes is two steps apart, closeness has almost no variance between nodes
and betweenness degenerates into a function of degree and clustering: they would
carry little information here even computed exactly and for free. (Exact
betweenness is also O(V.E) ~ 7.8e12 for sugarcane, so it is not free.)

Clustering is the measure that stays meaningful in a dense graph, because it
describes neighbourhood *shape* rather than distance.

### C(k) rises here, which was not the expectation

The plan for this step assumed C(k) ~ 1/k -- the falling curve of
preferential-attachment networks -- and predicted clustering would be largely
redundant with degree. **Measured, it is the opposite**, and the prediction was
wrong:

| degree decile | median degree | median clustering |
|---|---|---|
| 1 | 1 | 0.000 |
| 3 | 8 | 0.444 |
| 5 | 43 | 0.440 |
| 8 | 545 | 0.504 |
| **9** | **2,637** | **0.812** |
| 10 | 11,279 | 0.715 |

A thresholded correlation network is not a growth network. A hub here sits inside
a dense co-expressed module, so its neighbours are correlated with each other too
and its neighbourhood is close to complete. The two global statistics corroborate
it: mean local clustering **0.464** against global transitivity **0.690**, and
transitivity is triangle-weighted, hence dominated by exactly those hubs.

The practical consequence is that the **redundancy screen passes** for
clustering: rho(clustering, degree) = **+0.484**, and **73%** of its spread
survives within degree deciles. Clustering is a genuinely second description of
position in this network, not degree relabelled -- which is what made it worth
computing.

### Coreness was tested too, and failed the same screen

k-core number is a real centrality and cost five minutes on the same graph load
(`COMPUTE_CORENESS=1`; max core 6,668, median 43). It does not survive:

| measure | rho with degree | spread retained within degree deciles |
|---|---|---|
| clustering | +0.484 | **73%** |
| **coreness (log)** | **+0.995** | **11%** |

At rho = 0.995 coreness is degree wearing a different name, which is unsurprising
in a graph this dense -- a node's core number is bounded by its degree, and here
almost saturates it.

**It is therefore excluded from the models, not merely flagged.** That matters:
when coreness was left in, its collinearity with degree drove clustering's
coefficient to 0.0005 (dR2 2.5e-08, p = 0.98) purely through variance inflation.
Reporting that with a caveat attached would have corrupted the estimate this step
exists to make, rather than qualifying it. A measure that fails the screen is
described and set aside; only measures that pass enter a model. Its marginal
correlations are still reported (spearman with omega +0.013, with
constraint_score -0.072) as description, since those do not suffer from
collinearity.

### And it explains almost nothing about constraint

DeltaR-squared when each term is dropped -- the only currency in which the three
network terms compare:

| term | omega readout | constraint_score readout |
|---|---|---|
| **degree** | **0.0028** | **0.0029** |
| edge conservation | 0.00045 | 0.00010 |
| **clustering** | **0.00011** | 0.00045 |

Clustering's partial coefficient is -0.033 (p = 0.18, n.s.) under omega and
-0.057 (p = 0.018, Holm 0.036) under constraint_score. Its marginal correlation
even **flips sign between the two readouts** (+0.034 vs -0.063), and the
within-decile signs are inconsistent (4 of 10 negative for omega, 7 of 10 for
constraint_score). Degree remains 6-26x more explanatory than clustering, and
degree itself explains 0.3%.

Adding clustering does not disturb what was already concluded: conservation's
DeltaR-squared moves by -0.00004. Degree's moves by -0.0017, which is expected --
the two share rho = 0.48, so clustering absorbs a little of degree's share.

The `constraint_score` model is the cleanest statement available. Because that
score already has the synonymous rate, GC3, expression and length regressed out,
its R-squared is what network position explains of *residual* constraint:
**0.0064**. Network position, described three ways, accounts for well under 1%.

## The polyploid half

The result above covers 12,334 of ~103,000 sugarcane network nodes, because it
used strict 1:1 orthologs. Steps 08–13 look at the duplicated majority. The
question there is not the same one with more genes — it is sharper, and the 1:1
set structurally could not ask it.

### What the duplicates actually are

Of multi-copy families with a single sorghum anchor, **97%** have every copy on
*one chromosome across several haplotypes* — polyploid **homeologs**, the same
ancestral locus retained 2–8 times, not dispersed paralogs. Purple behaves the
same way. Each *pair* is classified too, because a homeolog family can contain a
tandem array inside one haplotype (`OG0000138` has three copies on 05C alone).

Families larger than `MAX_FAMILY_COPIES` (20) are labelled `oversized` and
excluded: R570 is ~10–12x and LA purple ~8x, so a genuine homeolog series cannot
be much bigger, and pair count grows quadratically — in purple, **136 such
families (up to 651 copies) would have contributed 70% of all copy pairs**.

| | sugarcane | purple |
|---|---|---|
| homeolog families | 12,227 | 16,778 |
| copy pairs | 93,512 | 285,811 |
| copies, median CDS identity | **0.981** | **0.985** |
| copy pairs identical in protein / in CDS | 14.3% / **5.4%** | — |

That last row is why every quantification filter here uses CDS, not protein:
salmon reads DNA, and protein identity overstates the problem nearly threefold.

### The claim that did not survive

The scoping observation was that copies of one locus differ by a median **28×**
in network degree. That number is a within-family max/min ratio, inflated by
family size. The honest pairwise figure is **6.4×** (sugarcane) and 6.6×
(purple) — and it has to be read against a baseline, which changes everything:

| | sugarcane | purple |
|---|---|---|
| homeologous copy pairs, median fold difference in degree | 6.4× | 6.6× |
| **expression-matched random gene pairs** | **15.0×** | **10.7×** |
| copies in different MCL modules | 77.6% | 63.9% |
| **random pairs in different modules** | **94.4%** | **84.0%** |

**Copies are more similar than chance, not less** — 2.3× and 1.6× closer in
degree, and markedly more likely to share a module. There is no decoupling to
report.

Nor does sequence identity predict network position: Spearman ρ = **0.048**
(sugarcane) and **0.010** (purple) against degree divergence, and the regression
coefficient on CDS identity **flips sign between the two species** (+0.75 vs
−1.03). Within expression strata the identity trend is flat; expression
divergence alone carries adjusted R² 0.026 / 0.017, and identity adds nothing
stable on top.

### And part of what is left is a mapping artefact

Homeologs 98% identical in CDS cannot be told apart by salmon, whose EM then
splits their reads on weak evidence. Counting 31-mers (salmon's default k) unique
to each copy within its family:

| | sugarcane | purple |
|---|---|---|
| copies with **< 5%** unique 31-mers | 44.0% | 68.0% |
| copies with **zero** unique 31-mers | 19.6% | **45.5%** |

Splitting the headline on that is decisive, and it points the wrong way for a
biological reading. Among near-identical pairs:

| | sugarcane | purple |
|---|---|---|
| copies salmon **can** separate | 6.31× (n=300) | **3.73×** (n=714) |
| copies that are **k-mer ambiguous** | 8.41× (n=8,262) | **7.08×** (n=52,550) |

Copies that can be quantified independently are *more* alike in the network than
copies that cannot, and the continuous version agrees at scale: the coefficient
on `frac_unique_min` is −1.53 (p = 6.8e-08) and −1.85 (p = 1.6e-31). Apparent
expression divergence rises the same way (median 0.46 → 0.80 and 0.38 → 0.85),
which is backwards for biology — sequences that are harder to distinguish should
not be *more* differently expressed — and is exactly the EM-splitting signature.

### Copy-specific selection is not measurable here, and the gate said so first

Step 11 gives every copy its own ω against the family's sorghum anchor, then
asks — before any within-family test — whether copies are distinguishable at all.
They are not:

| | sugarcane | purple |
|---|---|---|
| copies with an ω | 49,984 / 51,263 | 101,815 / 105,000 |
| **ICC(1) over families** | **0.940** | **0.939** |
| MS between families / MS within | 0.128 / 0.0020 | 0.168 / 0.0020 |
| verdict | **gate shut** | **gate shut** |

**94% of ω variance is between families and 6% within them.** ω is a property of
the gene, not of which copy of it you look at — which is what 98% identity and a
shared branch back to sorghum predict. No within-family selection test can work
on this, so none was run, and branch models were not escalated to. This is the
gate doing its job: the alternative was a well-formed table of within-family
comparisons with no power behind it.

One thing worth keeping from that table anyway. Median ω is **0.1819** across
49,599 polyploid copies, against **0.1827** across the 12,334 single-copy
orthologs — indistinguishable. The duplicated fraction of the genome is under the
same average purifying selection as the single-copy fraction, with no sign of the
relaxation duplicates are classically expected to show. It also retires the
"biased slice" caveat on the 1:1 result at the level of ω itself: the tidy slice
was not, in this respect, unrepresentative.

### The caveat this leaves for the rest of the project

This one reaches past the dN/dS stage. Among the genes scorable here:

| | sugarcane | purple |
|---|---|---|
| network nodes that are k-mer-ambiguous copies | 12.3% | **30.5%** |
| network nodes with **zero** unique 31-mers | 5.1% | **18.9%** |

Nearly a fifth of purple's network nodes carry expression values that are not
independently identifiable from the reads. They have lower degree than the rest
(median 633 vs 878), so they are not driving the hub structure, but they are
present in every module and every conservation count. This is a lower bound —
only genes in multi-copy families with a single sorghum anchor were scored
(29% and 49% of each network).

## Caveat that belongs with any result from this stage

Single-copy orthogroups are a **biased slice**: by construction they are the more
conserved, less duplicated part of a polyploid genome, and they cover a minority
of the ~103k network nodes. Anything found here generalises to *"among
single-copy orthologs"*, not to the whole network.
