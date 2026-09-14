#!/usr/bin/env python3
"""Independent check of 62's edge join, by a different route.

62 packs indices into int64 keys and binary-searches a sorted array. This does it
the other way: resolve gene NAMES through a fresh parse of Orthogroups.tsv, build
candidate index pairs as plain Python tuples, and confirm membership with ONE awk
pass over purple.pairs. Nothing is shared with 62 but the input files.

TWO directions, because only the second can fail informatively:
  POSITIVE  every edge 62 called conserved must have >= 1 candidate pair present
  NEGATIVE  a sample of edges 62 did NOT call conserved, whose BOTH endpoints have
            an in-network ortholog, must have ZERO candidate pairs present
A join that over-matches passes the positive test and fails the negative one.
"""
import csv, random, subprocess, sys, os

W = "/dados04/jorge/tmp/mcl_work_cluster"
OG = "/dados04/jorge/comparative_saccharum/files/fix_orthofinder/proteins/OrthoFinder/Results_Jun04_2/Orthogroups/Orthogroups.tsv"
CONS = sys.argv[1] if len(sys.argv) > 1 else \
    "/dados04/jorge/comparative_saccharum/new_clean/results/conservation/conserved_edges_sugarcane_to_purple_pearson.tsv"
N = int(sys.argv[2]) if len(sys.argv) > 2 else 4000
random.seed(20260913)

def tab(study):
    n2i, i2n = {}, {}
    for line in open(f"{W}/{study}.tab"):
        i, nm = line.rstrip("\n").split("\t", 1)
        n2i[nm] = int(i); i2n[int(i)] = nm
    return n2i, i2n

sc_n2i, sc_i2n = tab("sugarcane")
pu_n2i, pu_i2n = tab("purple")
print(f"sugarcane {len(sc_n2i):,} nodes | purple {len(pu_n2i):,} nodes")

# orthology: sugarcane gene NAME -> list of purple INDICES (fresh parse)
orth = {}
with open(OG) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        sg = [x.strip() for x in (row["sugarcane_one_transcript"] or "").split(",") if x.strip()]
        pg = [x.strip() for x in (row["one_transcript_purple_proteins"] or "").split(",") if x.strip()]
        pidx = [pu_n2i[g] for g in pg if g in pu_n2i]
        if not pidx: continue
        for g in sg:
            if g in sc_n2i: orth.setdefault(g, []).extend(pidx)
print(f"sugarcane genes with an in-network purple ortholog: {len(orth):,}")

# --- POSITIVE sample: edges 62 says are conserved ---------------------------
pos = []
with open(CONS) as fh:
    next(fh)
    for k, line in enumerate(fh):
        if len(pos) < N: pos.append(line.split("\t")[:2])
        elif random.random() < N / (k + 1): pos[random.randrange(N)] = line.split("\t")[:2]
print(f"positive sample: {len(pos):,} conserved edges")

# --- NEGATIVE sample: sugarcane edges NOT in that file, both ends mappable ---
consset = set()
with open(CONS) as fh:
    next(fh)
    for line in fh:
        a, b = line.split("\t")[:2]
        consset.add((sc_n2i[a], sc_n2i[b]) if sc_n2i[a] < sc_n2i[b] else (sc_n2i[b], sc_n2i[a]))
print(f"conserved edge set: {len(consset):,}")

neg = []
seen = 0
for line in open(f"{W}/sugarcane.pairs"):
    i, j, _ = line.split("\t")
    i, j = int(i), int(j)
    key = (i, j) if i < j else (j, i)
    if key in consset: continue
    if sc_i2n[i] not in orth or sc_i2n[j] not in orth: continue
    seen += 1
    if len(neg) < N: neg.append([sc_i2n[i], sc_i2n[j]])
    elif random.random() < N / seen: neg[random.randrange(N)] = [sc_i2n[i], sc_i2n[j]]
print(f"negative sample: {len(neg):,} of {seen:,} mappable non-conserved edges")

# --- candidate purple index pairs for both samples --------------------------
def candidates(sample):
    out = {}
    for n, (a, b) in enumerate(sample):
        s = set()
        for x in orth.get(a, []):
            for y in orth.get(b, []):
                if x != y: s.add((x, y) if x < y else (y, x))
        out[n] = s
    return out

cpos, cneg = candidates(pos), candidates(neg)
want = {}
for tag, c in (("P", cpos), ("N", cneg)):
    for n, s in c.items():
        for p in s: want.setdefault(p, []).append((tag, n))
print(f"candidate purple pairs to look up: {len(want):,}")

wf = "/tmp/claude-1004/-dados04-jorge-comparative-saccharum/44f4c7d1-8f81-46f0-968f-db386687cbda/scratchpad/_want.txt"
with open(wf, "w") as fh:
    for (x, y) in want: fh.write(f"{x}_{y}\n")

print("one awk pass over purple.pairs (14.5 GB) ...")
cmd = (f"awk -F'\\t' 'BEGIN{{while((getline l<\"{wf}\")>0) w[l]=1}} "
       f"{{k=($1<$2)?$1\"_\"$2:$2\"_\"$1; if(k in w) print k}}' {W}/purple.pairs")
found = set(subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout.split())
print(f"  {len(found):,} of {len(want):,} candidate pairs are real purple edges")

pos_ok = sum(1 for n, s in cpos.items() if any(f"{x}_{y}" in found for x, y in s))
neg_bad = sum(1 for n, s in cneg.items() if any(f"{x}_{y}" in found for x, y in s))
print()
print(f"POSITIVE: {pos_ok:,} of {len(pos):,} conserved edges confirmed "
      f"({'PASS' if pos_ok == len(pos) else 'FAIL'})")
print(f"NEGATIVE: {neg_bad:,} of {len(neg):,} non-conserved edges wrongly match "
      f"({'PASS' if neg_bad == 0 else 'FAIL'})")
sys.exit(0 if pos_ok == len(pos) and neg_bad == 0 else 1)
