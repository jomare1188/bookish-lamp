#!/usr/bin/env python3
"""
08_homeolog_families.py -- define the polyploid gene families, per species.

A family here is an orthogroup with EXACTLY ONE sorghum gene (the diploid
anchor) and >=2 copies in the study species. Each is classified from the copies'
gene ids:

  homeolog   every copy on ONE chromosome, across >=2 haplotypes. This is the
             polyploid signature -- the same ancestral locus retained several
             times by sugarcane's recent hybrid polyploidy -- and it is 97% of
             the parseable multi-copy families.
  dispersed  copies span chromosomes: older, relocated duplicates. Carried as a
             labelled contrast, because if the decoupling this stage tests is
             about polyploidy it should be STRONGER in homeologs.
  unplaced   ids off-grammar (scaffold/unanchored models). Kept and labelled, not
             silently dropped.
  oversized  more copies than MAX_FAMILY_COPIES. R570 is ~10-12x and LA purple
             ~8x, so a real homeolog series of one locus cannot be much bigger;
             beyond that these are lumped superfamilies. They are also
             quadratically expensive -- in purple, 136 such families (up to 651
             copies) would contribute 70% of all copy pairs and dominate every
             statistic. Labelled and excluded from the analysis, not deleted.

Id grammar, verified on both proteomes:
    sugarcane  SoffiXsponR570.<chr><HAP>g<locus>   01Ag001800    -> chr 01, hap A
    purple     Soffic.<chr>G<locus>-<chr><HAP>     01G0000010-1A -> chr 01, hap A

NOTE the locus number is NOT shared across haplotypes: 01Ag001800 and
01Bg001800 are unrelated genes. Homeology is called from orthogroup membership
plus a shared chromosome, never from the locus number.
"""
import os, re, sys
from collections import Counter

PAT = {
    "sugarcane": re.compile(r"^SoffiXsponR570\.(\d+)([A-Z])g(\d+)$"),
    "purple":    re.compile(r"^Soffic\.(\d+)[A-Z](\d+)-(\d+)([A-Z])$"),
}


def parse(study, gene):
    """-> (chromosome, haplotype) or None if the id is off-grammar."""
    m = PAT[study].match(gene)
    if not m:
        return None
    return (m.group(1), m.group(2)) if study == "sugarcane" else (m.group(3), m.group(4))


def read_tsv(path):
    """OrthoFinder writes CRLF; splitting on \\n alone glues \\r to the last column."""
    return open(path, "rb").read().decode().replace("\r\n", "\n").rstrip("\n").split("\n")


def main():
    env = os.environ
    of3 = env["OF3_OUT"]
    res = os.path.join(of3, sorted(d for d in os.listdir(of3) if d.startswith("Results_"))[0])
    lines = read_tsv(os.path.join(res, "Orthogroups", "Orthogroups.tsv"))
    h = lines[0].split("\t")
    i = {s: h.index(s) for s in ("purple", "sorghum", "sugarcane")}
    min_copies = int(env["FAMILIES_MIN_COPIES"])
    max_copies = int(env["MAX_FAMILY_COPIES"])
    subset = int(env.get("DNDS_SUBSET", "0") or 0)

    for study in ("sugarcane", "purple"):
        rows, cls = [], Counter()
        for l in lines[1:]:
            p = l.split("\t")
            p += [""] * (4 - len(p))
            sb = [x for x in p[i["sorghum"]].split(", ") if x]
            if len(sb) != 1:
                continue
            copies = [x for x in p[i[study]].split(", ") if x]
            if len(copies) < min_copies:
                continue
            if len(copies) > max_copies:
                cls["oversized"] += 1
                rows.append((p[0], study, "oversized", len(copies), "",
                             ",".join(copies), sb[0]))
                continue
            loc = [parse(study, g) for g in copies]
            if any(x is None for x in loc):
                k, chrom = "unplaced", ""
            else:
                chrs = {c for c, _ in loc}
                haps = {hp for _, hp in loc}
                if len(chrs) == 1 and len(haps) >= 2:
                    k, chrom = "homeolog", next(iter(chrs))
                elif len(chrs) == 1:
                    k, chrom = "single_haplotype", next(iter(chrs))
                else:
                    k, chrom = "dispersed", ""
            cls[k] += 1
            rows.append((p[0], study, k, len(copies), chrom, ",".join(copies), sb[0]))

        if subset > 0:
            keep = [r for r in rows if r[2] == "homeolog"][:subset]
            keep += [r for r in rows if r[2] == "dispersed"][:max(20, subset // 10)]
            rows = keep
            print("   DNDS_SUBSET=%d -- truncating to %d families" % (subset, len(rows)))

        out = os.path.join(env["OUTDIR"], "families_%s.tsv" % study)
        with open(out, "w") as fh:
            fh.write("orthogroup\tstudy\tclass\tn_copies\tchromosome\tcopies\tsorghum_anchor\n")
            for r in rows:
                fh.write("%s\t%s\t%s\t%d\t%s\t%s\t%s\n" % r)
        print("== %s: %d families with a single sorghum anchor and >=%d copies"
              % (study, len(rows), min_copies))
        for k, v in cls.most_common():
            print("   %-18s %6d" % (k, v))
        ncop = sum(r[3] for r in rows if r[2] == "homeolog")
        print("   copies in homeolog families: %d" % ncop)
        print("   -> %s" % out)


if __name__ == "__main__":
    main()
