# dN/dS vs network conservation

Does a gene whose co-expression neighbourhood is conserved between the two
networks also evolve under stronger purifying selection?

Full write-up, including why the design looks like this:
[`../../docs/dnds.md`](../../docs/dnds.md).

## Run it

```bash
./run_all.sh                    # every step
./run_all.sh 03                 # from step 03
DNDS_SUBSET=200 ./run_all.sh 02b   # 200-orthogroup smoke test
```

or `./run.sh dnds [from-step]` from `new_clean/`.

**Two steps are long and are meant to be launched deliberately:**

- `02_orthofinder_3sp.sh` — three-species OrthoFinder, **~3–4 h** at `-t 256`
- `04_codeml.sh` — codeml over every alignment

`05_conserved_degree.sh` defaults to sugarcane only; the purple direction is a
37 GB edge table and ~40 min (`./05_conserved_degree.sh both`).

## Headline result

**Network edge conservation does not predict selective constraint** — marginal
rho = 0.016 (p = 0.08), the partial coefficient adds 0.04% of adjusted R2, and
only 3 of 10 degree deciles even carry the right sign. The positive control does
work (hubs are more constrained, p = 6e-4), so this is a real null rather than a
dead pipeline. The apparent lineage-specific decline (omega 0.393 -> 0.301 across
conservation bins) is an **expression** effect: it vanishes within expression
tertiles. See [`../../docs/dnds.md`](../../docs/dnds.md).

## A second network measure (steps 14-16)

Local clustering coefficient, added because the constraint analysis tested only
degree. Betweenness and closeness are excluded on information grounds, not cost:
at mean degree 1,475 the diameter is 2-3 hops, so closeness has almost no
variance and betweenness collapses toward a degree function.

**C(k) RISES in this network** -- 0.00 at degree 1 to 0.81 at degree 2,637 --
the opposite of the C(k) ~ 1/k assumed when planning it. A thresholded
correlation network is not a growth network: a hub sits inside a dense module, so
its neighbourhood is near-complete. So the redundancy screen PASSES
(rho with degree = +0.48, 73% of spread retained within degree deciles) and
clustering really is a second description of position.

It still explains almost nothing. dR2 when dropped: degree 0.0028, edge
conservation 0.00045, clustering 0.00011 (omega) / 0.00045 (constraint score).
Its marginal correlation flips sign between the two readouts. With
`constraint_score` -- which already has dS, GC3, expression and length regressed
out -- network position described three ways explains **0.6%** of what is left.

Step 06 now also writes `constraint_score`: the residual of log(dN) on log(dS),
GC3, expression and length. omega assumes dN scales 1:1 with dS; the fitted slope
is 0.511 alone and 0.725 with covariates, so omega over-corrects. The two readouts
agree at rho = +0.82 and every conclusion above holds under both.

## The polyploid half (steps 08-13)

Steps 01-07 used strict 1:1 orthologs -- 12,334 of ~103k network nodes. Steps
08-13 look at the duplicated majority, where 97% of multi-copy families are
polyploid HOMEOLOGS (one chromosome, several haplotypes), ~98% identical in CDS.

**That result did not survive its controls either**, and the controls are the
output worth reading:

- Copies differ by a median **6.4x** in network degree -- but expression-matched
  random gene pairs differ by **15x**. Copies are MORE similar than chance, not
  less. (The 28x figure from scoping was a within-family max/min ratio inflated
  by family size; the pairwise median is 6.4x.)
- 77.6% of copy pairs sit in different MCL modules, against a chance rate of
  **94.4%**.
- Sequence identity does not predict network position: rho = 0.048 (sugarcane),
  0.010 (purple), and the regression coefficient flips sign between species.
- **Part of what remains is a mapping artefact.** 44% of sugarcane and 68% of
  purple copies have <5% unique 31-mers; 19.6% and 45.5% have NONE. Among
  near-identical pairs, copies salmon CAN separate diverge less (6.31x / 3.73x)
  than copies it cannot (8.41x / 7.08x) -- the opposite of a biological signal.

Step 11 also gates copy-specific selection before testing it: **ICC = 0.94 in
both species**, i.e. 94% of omega variance is between families and 6% within, so
copies of a gene are not distinguishable in omega and no within-family test was
run. Median omega is 0.1819 over 49,599 polyploid copies against 0.1827 over the
12,334 single-copy orthologs -- the duplicated genome is under the same average
constraint, with no relaxation.

**A caveat this leaves for the whole project:** 12.3% of sugarcane network nodes
and **30.5% of purple's** are copies that cannot be quantified independently of a
sibling; 18.9% of purple's nodes have zero unique 31-mers. Lower bound -- only
multi-copy families with a single sorghum anchor were scored.

## The one thing to know before reading any output

R570 and LA purple are 95–99% identical at the protein level, so **~40% of genes
have zero synonymous differences between them**. Per-gene ω for that pair does
not exist. It is used only through summed NG86 counts within bins.

Per-gene ω comes from a *Sorghum bicolor* outgroup instead (dS ≈ 0.21), which
step 00 downloads in full from NCBI and MD5-verifies. No sorghum data is read
from anywhere else on the machine.

## Two OrthoFinder facts this stage had to work around

1. **3.1.3's `Orthogroups_SingleCopyOrthologues.txt` is not trustworthy here** —
   it lists the right *number* of ids (532) but only 185 are really 1:1:1. Never
   read it; derive the set from `Orthogroups.tsv`.
2. **Its `Orthogroups.tsv` is CRLF**, so splitting on `\n` alone glues `\r` onto
   every sugarcane gene id and silently breaks every downstream lookup.

Sorghum is attached to the two-species 1:1 pairs through the three-species run's
tree-based `Orthologues/` tables, not through 1:1:1 orthogroups (only 532 exist,
because a third species coarsens the clustering). 20,616 pairs survive, and where
both sides resolve to one sorghum gene they agree on which 99.9% of the time.

## Two things this stage must never do

1. **Repoint `ORTHOGROUPS`.** Every conservation number in the repo was computed
   from the two-species run. Step 02 checksums it before and after and fails if
   it moved.
2. **Report the marginal correlation as the answer.** Degree, expression and
   conservation are mutually correlated; expression is the strongest known
   predictor of ω anywhere. The partial coefficient (test 2) is the result.
