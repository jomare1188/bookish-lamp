#!/usr/bin/env python3
"""
10_kmer_uniqueness.py -- can salmon actually tell these copies apart?

This is the control that decides whether the whole polyploid result is biology
or plumbing. Homeologs are ~98% identical; salmon assigns multireads by EM, and
where two copies share nearly all their k-mers the split between them is weakly
identified. That would manufacture expression divergence between copies, and
therefore network divergence, with nothing biological behind it.

So: for every copy, count the k-mers it does NOT share with any sibling in its
own family. k defaults to 31, salmon's default, and k-mers are canonicalised
(min of the k-mer and its reverse complement) the way salmon indexes them.

A copy with frac_unique near zero cannot be quantified independently of its
siblings, full stop. Step 12 splits its headline on this and reports both sides
rather than filtering quietly.
"""
import os, sys
from collections import defaultdict

COMP = str.maketrans("ACGTN", "TGCAN")


def canon(k):
    rc = k.translate(COMP)[::-1]
    return k if k <= rc else rc


def read_fasta(path):
    d, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    d[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        d[name] = "".join(buf)
    return d


def main():
    env = os.environ
    K = int(env["KMER_K"])
    floor = float(env["MIN_FRAC_UNIQUE"])
    for study in ("sugarcane", "purple"):
        cds = read_fasta(env["SC_CDS_CLEAN"] if study == "sugarcane" else env["PU_CDS_CLEAN"])
        fam_path = os.path.join(env["OUTDIR"], "families_%s.tsv" % study)
        out = os.path.join(env["OUTDIR"], "kmer_uniqueness_%s.tsv" % study)
        n_rows = n_below = 0
        with open(fam_path) as fh, open(out, "w") as fo:
            fh.readline()
            fo.write("gene\torthogroup\tclass\tn_kmers\tn_unique\tfrac_unique\n")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[2] not in ("homeolog", "dispersed"):
                    continue
                og, klass = f[0], f[2]
                copies = [g for g in f[5].split(",") if g in cds]
                if len(copies) < 2:
                    continue
                sets = {}
                for g in copies:
                    s = cds[g].upper()
                    sets[g] = {canon(s[i:i + K]) for i in range(len(s) - K + 1)} if len(s) >= K else set()
                # how often each k-mer occurs across the family
                cnt = defaultdict(int)
                for s in sets.values():
                    for km in s:
                        cnt[km] += 1
                for g in copies:
                    tot = len(sets[g])
                    uniq = sum(1 for km in sets[g] if cnt[km] == 1)
                    fr = uniq / tot if tot else 0.0
                    fo.write("%s\t%s\t%s\t%d\t%d\t%.5f\n" % (g, og, klass, tot, uniq, fr))
                    n_rows += 1
                    if fr < floor:
                        n_below += 1
        print("== %s: %d copies scored (k=%d)" % (study, n_rows, K))
        print("   below the frac_unique floor of %.2f: %d (%.1f%%) -- these copies"
              % (floor, n_below, 100.0 * n_below / n_rows if n_rows else 0))
        print("   cannot be quantified independently of their siblings")
        print("   -> %s" % out)


if __name__ == "__main__":
    main()
