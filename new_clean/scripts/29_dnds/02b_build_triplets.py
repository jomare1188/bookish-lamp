#!/usr/bin/env python3
"""
02b_build_triplets.py -- attach the sorghum outgroup to the two-species pairs.

WHY NOT THE 1:1:1 ORTHOGROUPS (the obvious route, and the one first tried):

  The three-species run yields only **532** orthogroups that are strictly
  single-copy in all three species, against 33,982 single-copy PAIRS in the
  two-species run. Adding sorghum does not merely filter that set -- sorghum
  bridges Saccharum subfamilies that MCL had kept apart, so the clustering
  coarsens (50,494 orthogroups, mean size 7.3, against 94,229 before) and
  almost nothing survives as 1:1:1. 532 genes cannot carry a regression with
  five covariates, ten degree deciles and five conservation bins.

  Worse, OrthoFinder 3.1.3's `Orthogroups_SingleCopyOrthologues.txt` cannot be
  trusted here: it lists 532 ids -- the right COUNT -- but only 185 of them are
  actually 1:1:1. The rest have shapes like (3,0,0) and (2,0,1), i.e. genes in
  one species only. The count is computed correctly and the id list is not, so
  the file is not read at all; the 1:1:1 set is derived from Orthogroups.tsv
  when it is needed.

WHAT IS DONE INSTEAD:

  The analysis unit stays the **two-species 1:1 pair** -- exactly the set every
  conservation number in this repo was computed from, so the dN/dS layer and
  the conservation layer are talking about the same genes by construction
  rather than by a gate applied afterwards.

  Sorghum is attached to each pair through the three-species run's TREE-BASED
  orthologue tables (Orthologues/), not through best hits and not through
  orthogroup membership. A pair is kept when the sugarcane gene has exactly one
  sorghum orthologue, the purple gene has exactly one, and they are the same
  gene.

  That last condition is a strong, free quality check, and it passes
  emphatically: 20,616 of 20,637 pairs (99.9%) agree on which sorghum gene is
  the orthologue. The 21 that disagree are dropped and counted.
"""
import os, sys
from collections import defaultdict


def read_tsv(path):
    """OrthoFinder writes CRLF here; splitting on \\n alone leaves \\r glued to
    the last column, which silently corrupts every sugarcane gene id."""
    raw = open(path, "rb").read().decode()
    return raw.replace("\r\n", "\n").rstrip("\n").split("\n")


def load_orthologues(path, a, b):
    """gene in species a -> set of its orthologues in species b"""
    lines = read_tsv(path)
    h = lines[0].split("\t")
    ia, ib = h.index(a), h.index(b)
    m = defaultdict(set)
    for l in lines[1:]:
        p = l.split("\t")
        if len(p) <= max(ia, ib):
            continue
        A = [x for x in p[ia].split(", ") if x]
        B = [x for x in p[ib].split(", ") if x]
        for g in A:
            m[g].update(B)
    return m


def main():
    env = os.environ
    of3_out = env["OF3_OUT"]
    results = sorted(d for d in os.listdir(of3_out) if d.startswith("Results_"))
    if not results:
        sys.exit("FATAL: no Results_* under %s -- run ./02_orthofinder_3sp.sh first" % of3_out)
    res = os.path.join(of3_out, results[0])
    print("== three-species run: %s" % res)

    sc2sb = load_orthologues(
        os.path.join(res, "Orthologues", "Orthologues_sugarcane", "sugarcane__v__sorghum.tsv"),
        "sugarcane", "sorghum")
    pu2sb = load_orthologues(
        os.path.join(res, "Orthologues", "Orthologues_purple", "purple__v__sorghum.tsv"),
        "purple", "sorghum")
    print("   sugarcane genes with a sorghum orthologue: %d" % len(sc2sb))
    print("   purple    genes with a sorghum orthologue: %d" % len(pu2sb))

    # --- the analysis unit: the two-species single-copy pairs --------------
    og2 = env["ORTHOGROUPS"]
    sco_path = os.path.join(os.path.dirname(og2), "Orthogroups_SingleCopyOrthologues.txt")
    sco2 = set(l.strip() for l in open(sco_path) if l.strip())
    lines = read_tsv(og2)
    h = lines[0].split("\t")
    i_pu = 1 if "purple" in h[1] else 2
    i_sc = 3 - i_pu
    pairs = []
    for l in lines[1:]:
        p = l.split("\t")
        if p[0] in sco2 and len(p) >= 3:
            pairs.append((p[0], p[i_sc].strip(), p[i_pu].strip()))
    print("== two-species 1:1 pairs (the conservation layer's own set): %d" % len(pairs))

    kept, no_orth, multi, disagree = [], 0, 0, 0
    for og, sc, pu in pairs:
        a, b = sc2sb.get(sc), pu2sb.get(pu)
        if not a or not b:
            no_orth += 1
            continue
        if len(a) != 1 or len(b) != 1:
            multi += 1
            continue
        if a != b:
            disagree += 1
            continue
        kept.append((og, sc, pu, next(iter(a))))

    subset = int(env.get("DNDS_SUBSET", "0") or 0)
    if subset > 0:
        kept = kept[:subset]
        print("   DNDS_SUBSET=%d -- truncating to %d triplets" % (subset, len(kept)))

    out = env["TRIPLETS"]
    with open(out, "w") as fh:
        fh.write("orthogroup\tsugarcane_gene\tpurple_gene\tsorghum_gene\torthogroup_2sp\n")
        for og, sc, pu, sb in kept:
            fh.write("%s\t%s\t%s\t%s\t%s\n" % (og, sc, pu, sb, og))

    n = len(pairs)
    print("== attaching sorghum")
    print("   dropped, no sorghum orthologue on one side : %d (%.1f%%)" % (no_orth, 100.0 * no_orth / n))
    print("   dropped, more than one on a side           : %d (%.1f%%)" % (multi, 100.0 * multi / n))
    print("   dropped, the two sides name DIFFERENT genes: %d (%.2f%%)" % (disagree, 100.0 * disagree / n))
    resolved = len(kept) + disagree
    if resolved:
        print("   cross-side agreement on the outgroup       : %d/%d (%.1f%%)"
              % (len(kept), resolved, 100.0 * len(kept) / resolved))
    print("   kept                                       : %d -> %s" % (len(kept), out))


if __name__ == "__main__":
    main()
