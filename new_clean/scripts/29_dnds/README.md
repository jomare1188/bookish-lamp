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
