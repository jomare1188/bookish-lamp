#!/usr/bin/env python3
# =============================================================================
# 06_myb59_neighbourhoods.py -- the AtMYB59 copies as NETWORK objects.
#
# Figure 7 shows that Soffic.09G0001580-9H is the only expressed AtMYB59-anchor copy
# in purple and that its nitrogen response is U-shaped. This asks what these genes
# look like in the graphs: are they hubs, and which sugarcane copy has the most
# similar NEIGHBOURHOOD to the purple one.
#
# ---------------------------------------------------------------------------
# WHAT IS ALREADY KNOWN BEFORE THIS RUNS, and must be stated before the result
# rather than fitted to it:
#
#   THEY ARE NOT HUBS. Degrees 76-2,198 in sugarcane (54.6th-84.7th percentile) and
#   635 for the purple copy -- BELOW purple's median. Munoz define Module 20 as
#   high-betweenness on THEIR network; in ours these are ordinary mid-degree genes.
#
#   ONLY ONE SUGARCANE COPY IS THE ORTHOLOG. The nine are one locus only by
#   Arabidopsis anchor; OrthoFinder splits them across >= 4 orthogroups, and exactly
#   one pairs 1:1 with the purple gene -- 06Ag012300, in OG0088735, which is also the
#   highest-degree copy. So "which copy is most alike?" has a PRIOR ANSWER, and the
#   neighbour comparison is a test of it: confirmation if 06Ag012300 also wins on
#   neighbours, and a more interesting result if it does not.
# ---------------------------------------------------------------------------
#
# THE OVERLAP NEEDS A NULL OR IT RANKS BY SIZE. Two large gene sets always share
# members. For each sugarcane copy the observed overlap with the purple gene's
# neighbourhood is therefore compared against the overlap with random purple genes of
# MATCHED DEGREE -- the question is not "do they share neighbours" but "more than
# this copy shares with any purple gene of that size".
#
# THREE SIMILARITY MEASURES, because one would mislead. The sets differ ~30-fold in
# size (76 to 2,198 neighbours), so Jaccard alone ranks mostly by size; the overlap
# coefficient |A n B| / min(|A|,|B|) does not. The raw count is carried too, because
# both ratios hide it.
#
# COPIES ARE NOT INDEPENDENT. 6 of the 9 sit in sugarcane's Module_010, so their
# neighbourhoods overlap each other heavily and nine copies are not nine
# measurements. The pairwise within-sugarcane overlap is written out so that is
# visible rather than assumed.
#
# RUN: ./run_all.sh 06     (or python3 06_myb59_neighbourhoods.py)
# =============================================================================

import csv
import os
import re
import sys
import time

import numpy as np
import pandas as pd

BASE = "/dados04/jorge/comparative_saccharum"
CLEAN = f"{BASE}/new_clean"
WORK = "/dados04/jorge/tmp/mcl_work_cluster"
RES = f"{CLEAN}/results"
OUT = f"{RES}/readouts/module20"
NBR_DIR = f"{OUT}/myb59_neighbourhoods"
OG_FILE = (f"{BASE}/files/fix_orthofinder/proteins/OrthoFinder/"
           f"Results_Jun04_2/Orthogroups/Orthogroups.tsv")
ANCHOR = "AT5G59780"
FOCUS = "Soffic.09G0001580-9H"
SP_COL = {"sugarcane": "sugarcane_one_transcript",
          "purple": "one_transcript_purple_proteins"}
N_NULL = 200
DEG_TOL = 0.20
SEED = 1188

rng = np.random.default_rng(SEED)
os.makedirs(NBR_DIR, exist_ok=True)


def say(*a):
    print("[%s]" % time.strftime("%H:%M:%S"), *a, flush=True)


_RE_P = re.compile(r"\.p[0-9]+$")
_RE_V = re.compile(r"\.v[0-9]+(\.[0-9]+)*$")
strip_version = lambda g: _RE_V.sub("", _RE_P.sub("", g))

EXPECT_EDGES = {"sugarcane": 75_333_769, "purple": 675_955_918}
EXPECT_PAIRS, EXPECT_OG = 114_681, 43_390

say("AtMYB59 copies: degree, module, neighbourhood, cross-species overlap")

# --- the gene set, from the table figure 7 uses ------------------------------
tab = [r for r in csv.DictReader(open(f"{OUT}/module20_network_table.tsv"), delimiter="\t")
       if r["at_anchor"] == ANCHOR]
TARGET = {sp: sorted({r["gene"] for r in tab if r["species"] == sp})
          for sp in ("sugarcane", "purple")}
if (len(TARGET["sugarcane"]), len(TARGET["purple"])) != (9, 3):
    sys.exit(f"FATAL: expected 9 sugarcane and 3 purple {ANCHOR} copies, got "
             f"{len(TARGET['sugarcane'])} and {len(TARGET['purple'])}")
MYB = {r["gene"]: r["our_is_MYB"].strip() in ("TRUE", "True", "1") for r in tab}
say(f"  {ANCHOR} copies: 9 sugarcane, 3 purple")


# --- index maps and node metrics ---------------------------------------------
def read_tab(study):
    names = []
    with open(f"{WORK}/{study}.tab") as fh:
        for expect, line in enumerate(fh):
            i, nm = line.rstrip("\n").split("\t", 1)
            if int(i) != expect:
                sys.exit(f"FATAL: {study}.tab is not in index order at {expect}")
            names.append(nm)
    return np.array(names, dtype=object)


NAMES, IDX, DEG = {}, {}, {}
for sp in ("sugarcane", "purple"):
    NAMES[sp] = read_tab(sp)
    IDX[sp] = {g: i for i, g in enumerate(NAMES[sp])}
    d = pd.read_csv(f"{RES}/{sp}/network_{sp}_node_metrics.tsv", sep="\t",
                    usecols=["gene", "degree"])
    dd = dict(zip(d["gene"], d["degree"].astype(int)))
    DEG[sp] = np.array([dd[g] for g in NAMES[sp]], dtype=np.int64)
    say(f"  {sp}: {len(NAMES[sp]):,} nodes, median degree {int(np.median(DEG[sp])):,}")

# A copy can be absent from the graph -- Soffic.05G0022030-5F is (0 TPM, no variance).
# It is REPORTED as absent rather than silently dropped.
ABSENT = {sp: [g for g in TARGET[sp] if g not in IDX[sp]] for sp in TARGET}
for sp, gs in ABSENT.items():
    for g in gs:
        say(f"  NOTE: {g} ({sp}) is not a node of the Pearson-only network")
PRESENT = {sp: [g for g in TARGET[sp] if g in IDX[sp]] for sp in TARGET}
if FOCUS not in PRESENT["purple"]:
    sys.exit(f"FATAL: the focus gene {FOCUS} is not in the purple network")


# --- neighbours, one streaming pass per species ------------------------------
def neighbours(study, genes):
    """Neighbour index sets for `genes`. One pass; a target can be in EITHER column,
    which is the bug that would silently halve every degree."""
    want = {IDX[study][g] for g in genes}
    want_arr = np.fromiter(want, dtype=np.int64, count=len(want))
    out = {i: [] for i in want}
    n_rows = 0
    t0 = time.time()
    for ch in pd.read_csv(f"{WORK}/{study}.pairs", sep="\t", header=None,
                          names=["i", "j", "w"],
                          dtype={"i": np.int32, "j": np.int32, "w": np.float32},
                          chunksize=5_000_000, engine="c"):
        i = ch["i"].to_numpy(); j = ch["j"].to_numpy()
        n_rows += len(i)
        # BOTH COLUMNS. dump-upper writes each edge once, in one orientation only,
        # so a target appears sometimes as i and sometimes as j. Scanning one column
        # would silently halve every degree -- which is what the assertion below
        # exists to catch, and did.
        m = np.isin(i, want_arr)
        if m.any():
            for tgt_i, nbr in zip(i[m], j[m]):
                out[int(tgt_i)].append(nbr)
        m = np.isin(j, want_arr)
        if m.any():
            for tgt_i, nbr in zip(j[m], i[m]):
                out[int(tgt_i)].append(nbr)
    if n_rows != EXPECT_EDGES[study]:
        sys.exit(f"FATAL: streamed {n_rows:,} {study} edges, expected "
                 f"{EXPECT_EDGES[study]:,}")
    say(f"  {study}: scanned {n_rows:,} edges in {time.time() - t0:.0f}s")
    return {k: np.unique(np.asarray(v, dtype=np.int64)) for k, v in out.items()}


say("")
say("extracting neighbourhoods")
NBR = {sp: neighbours(sp, PRESENT[sp]) for sp in ("sugarcane", "purple")}

# THE NEIGHBOUR COUNT MUST EQUAL THE DEGREE ALREADY ON RECORD. If the second-column
# pass above were missing, every count here would be wrong and still look plausible.
for sp in ("sugarcane", "purple"):
    for g in PRESENT[sp]:
        k = IDX[sp][g]
        if len(NBR[sp][k]) != DEG[sp][k]:
            sys.exit(f"FATAL: {g} has {len(NBR[sp][k])} neighbours but degree "
                     f"{DEG[sp][k]} -- the edge scan is incomplete")
say("  every neighbour count matches the recorded degree")

for sp in ("sugarcane", "purple"):
    for g in PRESENT[sp]:
        with open(f"{NBR_DIR}/{sp}__{g}.txt", "w") as fh:
            fh.write("\n".join(NAMES[sp][NBR[sp][IDX[sp][g]]]) + "\n")
say(f"  wrote {sum(len(v) for v in PRESENT.values())} neighbour lists to "
    f"{os.path.basename(NBR_DIR)}/")


# --- orthology, in index space, restricted to in-network pairs ----------------
# The SAME universe 61 and 62 assert, so this readout and the conservation analysis
# project orthology identically. A different universe here would make the overlap
# numbers incomparable with everything else in the paper.
say("")
say("loading orthologs")
src, tgt, n_og = [], [], 0
with open(OG_FILE) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        S = [strip_version(x.strip())
             for x in (row[SP_COL["sugarcane"]] or "").split(",") if x.strip()]
        P = [strip_version(x.strip())
             for x in (row[SP_COL["purple"]] or "").split(",") if x.strip()]
        ia = [IDX["sugarcane"][g] for g in S if g in IDX["sugarcane"]]
        ib = [IDX["purple"][g] for g in P if g in IDX["purple"]]
        if not ia or not ib:
            continue
        n_og += 1
        for x in ia:
            for y in ib:
                src.append(x); tgt.append(y)
src = np.asarray(src, dtype=np.int64); tgt = np.asarray(tgt, dtype=np.int64)
if (len(src), n_og) != (EXPECT_PAIRS, EXPECT_OG):
    sys.exit(f"FATAL: ortholog universe is {len(src):,}/{n_og:,}, but 61 and 62 use "
             f"{EXPECT_PAIRS:,}/{EXPECT_OG:,}")
say(f"  {len(src):,} pairs over {n_og:,} orthogroups (matches 61/62)")

order = np.argsort(src, kind="stable")
SRC_SORTED, TGT_SORTED = src[order], tgt[order]
rev = np.argsort(tgt, kind="stable")
TGT_R, SRC_R = tgt[rev], src[rev]


def project(idx_set, s_sorted, t_sorted):
    """Map a set of source indices through orthology to target indices."""
    lo = np.searchsorted(s_sorted, idx_set, "left")
    hi = np.searchsorted(s_sorted, idx_set, "right")
    out = [t_sorted[a:b] for a, b in zip(lo, hi) if b > a]
    return np.unique(np.concatenate(out)) if out else np.zeros(0, dtype=np.int64)


# --- the comparison -----------------------------------------------------------
def scores(a, b):
    """raw intersection, Jaccard, overlap coefficient, and the fraction of B recovered.

    FOUR NUMBERS BECAUSE NO ONE OF THEM IS SUFFICIENT HERE.
      * Jaccard is dominated by set size when the sets differ 30-fold (86 to 2,540
        projected neighbours against the focus gene's 635).
      * The overlap coefficient divides by min(|A|,|B|), and that MINIMUM FLIPS
        between rows -- it is the sugarcane set for the small copies and the focus
        set for the large ones -- so read down a column it is not one statistic.
      * frac_of_b = |A n B| / |B| keeps ONE denominator for every row: what fraction
        of the purple gene's neighbourhood this copy reproduces. That is the
        comparable one, and the one the null and the ranking use.
      * the raw count, because all three ratios hide how few genes this rests on.
    """
    if len(a) == 0 or len(b) == 0:
        return 0, 0.0, 0.0, 0.0
    inter = np.intersect1d(a, b, assume_unique=True).size
    union = len(a) + len(b) - inter
    return (inter, inter / union, inter / min(len(a), len(b)), inter / len(b))


FOCUS_I = IDX["purple"][FOCUS]
FOCUS_NBR = NBR["purple"][FOCUS_I]
FOCUS_DEG = int(DEG["purple"][FOCUS_I])
say("")
say(f"focus: {FOCUS}  degree {FOCUS_DEG:,}")

# THE DEGREE-MATCHED NULL is what licenses any ranking below. Two large sets always
# intersect; the question is whether a copy shares MORE with this particular purple
# gene than with any purple gene of the same size. Candidates exclude the focus gene
# itself and the other AtMYB59 copies.
cand = np.flatnonzero((DEG["purple"] >= FOCUS_DEG * (1 - DEG_TOL)) &
                      (DEG["purple"] <= FOCUS_DEG * (1 + DEG_TOL)))
cand = np.setdiff1d(cand, np.array([IDX["purple"][g] for g in PRESENT["purple"]]))
say(f"  degree-matched null pool: {len(cand):,} purple genes within "
    f"+/-{int(100 * DEG_TOL)}% of {FOCUS_DEG:,}")
null_idx = rng.choice(cand, size=min(N_NULL, len(cand)), replace=False)
NULL_NBR = neighbours("purple", [NAMES["purple"][i] for i in null_idx])

# --- ortholog status, stated BEFORE the neighbour result ----------------------
og_of = {}
with open(OG_FILE) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        S = {strip_version(x.strip())
             for x in (row[SP_COL["sugarcane"]] or "").split(",") if x.strip()}
        P = {strip_version(x.strip())
             for x in (row[SP_COL["purple"]] or "").split(",") if x.strip()}
        for g in (S | P) & (set(TARGET["sugarcane"]) | set(TARGET["purple"])):
            og_of[g] = row["Orthogroup"]
FOCUS_OG = og_of.get(FOCUS)
say(f"  {FOCUS} is in {FOCUS_OG}; sugarcane copies sharing it: " +
    (", ".join(g for g in TARGET["sugarcane"] if og_of.get(g) == FOCUS_OG) or "(none)"))

rows = []
for g in PRESENT["sugarcane"]:
    gi = IDX["sugarcane"][g]
    nb = NBR["sugarcane"][gi]
    proj = project(nb, SRC_SORTED, TGT_SORTED)
    inter, jac, ovl, frac = scores(proj, FOCUS_NBR)
    # The null is on frac_of_b -- one denominator for every row. Null neighbourhoods
    # are degree-matched to the focus gene, so |B| is comparable throughout.
    nulls_frac = np.array([scores(proj, NULL_NBR[i])[3] for i in null_idx])
    nulls_cnt = np.array([scores(proj, NULL_NBR[i])[0] for i in null_idx])
    sd = nulls_frac.std(ddof=1)
    z = (frac - nulls_frac.mean()) / sd if sd > 0 else np.nan
    p = (int((nulls_frac >= frac).sum()) + 1) / (len(nulls_frac) + 1)
    rows.append(dict(
        sugarcane_gene=g, orthogroup=og_of.get(g, "-"),
        is_ortholog_of_focus=(og_of.get(g) == FOCUS_OG and FOCUS_OG is not None),
        our_is_MYB=MYB.get(g), degree=int(DEG["sugarcane"][gi]),
        n_neighbours=len(nb), n_projected=len(proj),
        focus_neighbours=len(FOCUS_NBR), shared=inter,
        shared_expected=round(float(nulls_cnt.mean()), 2),
        jaccard=round(jac, 6), overlap_coef=round(ovl, 6),
        frac_of_focus=round(frac, 6),
        null_mean_frac=round(float(nulls_frac.mean()), 6),
        null_sd=round(float(sd), 6),
        z=round(float(z), 3) if np.isfinite(z) else np.nan,
        p_emp=round(p, 4), n_null=len(nulls_frac)))

df = pd.DataFrame(rows).sort_values("frac_of_focus", ascending=False)
# BH across the nine copies. Nine tests on counts this small need it, and without it
# a p of 0.015 on FOUR shared genes reads as a finding.
pv = df["p_emp"].to_numpy()
o = np.argsort(pv); ranked = pv[o]
adj = np.minimum.accumulate((ranked * len(pv) / np.arange(1, len(pv) + 1))[::-1])[::-1]
bh = np.empty_like(adj); bh[o] = np.clip(adj, 0, 1)
df["p_adj"] = np.round(bh, 4)
f_over = f"{OUT}/myb59_neighbour_overlap.tsv"
df.to_csv(f_over, sep="\t", index=False)
say("")
say(f"wrote {os.path.basename(f_over)}")

# --- the reverse projection, since orthology is many-to-many ------------------
rev_rows = []
proj_focus = project(FOCUS_NBR, TGT_R, SRC_R)
for g in PRESENT["sugarcane"]:
    nb = NBR["sugarcane"][IDX["sugarcane"][g]]
    inter, jac, ovl, frac = scores(proj_focus, nb)
    rev_rows.append(dict(sugarcane_gene=g, focus_projected=len(proj_focus),
                         n_neighbours=len(nb), shared=inter,
                         jaccard=round(jac, 6), overlap_coef=round(ovl, 6),
                         frac_of_sugarcane=round(frac, 6)))
pd.DataFrame(rev_rows).sort_values("overlap_coef", ascending=False).to_csv(
    f"{OUT}/myb59_neighbour_overlap_reverse.tsv", sep="\t", index=False)

# --- how independent are the nine copies? -------------------------------------
# 6 of 9 sit in one MCL module, so their neighbourhoods are expected to overlap each
# other. Writing it out keeps "nine copies" from being read as nine measurements.
pw = []
for a in PRESENT["sugarcane"]:
    for b in PRESENT["sugarcane"]:
        if a >= b:
            continue
        i2, j2, o2, _f = scores(NBR["sugarcane"][IDX["sugarcane"][a]],
                                NBR["sugarcane"][IDX["sugarcane"][b]])
        pw.append(dict(gene_a=a, gene_b=b, shared=i2,
                       jaccard=round(j2, 6), overlap_coef=round(o2, 6)))
pd.DataFrame(pw).sort_values("overlap_coef", ascending=False).to_csv(
    f"{OUT}/myb59_within_sugarcane_overlap.tsv", sep="\t", index=False)

print()
print("=" * 78)
print(f"CROSS-SPECIES NEIGHBOURHOOD OVERLAP with {FOCUS} ({FOCUS_DEG:,} neighbours)")
print("=" * 78)
print(f"{'sugarcane copy':30s} {'OG':>11s} {'orth':>5s} {'deg':>6s} {'proj':>6s} "
      f"{'shared':>7s} {'exp':>6s} {'%focus':>7s} {'z':>6s} {'p':>7s} {'padj':>7s}")
for _, r in df.iterrows():
    print(f"{r.sugarcane_gene:30s} {r.orthogroup:>11s} "
          f"{'YES' if r.is_ortholog_of_focus else '-':>5s} {r.degree:6,} "
          f"{r.n_projected:6,} {r.shared:7,} {r.shared_expected:6.1f} "
          f"{100 * r.frac_of_focus:6.2f}% {r.z:6.2f} {r.p_emp:7.4f} {r.p_adj:7.4f}")
print()
print("  `%focus` is the share of the purple gene's 635 neighbours this copy's")
print("  projected neighbourhood recovers -- ONE denominator for every row. `exp` is")
print("  what a purple gene of matched degree gives. The ranking means nothing without")
print("  it: two large neighbourhoods always intersect.")
print()
print("  READ THE ABSOLUTE COUNTS. The best row shares a handful of genes out of 635,")
print("  so a z of 3 rests on 4 genes against an expectation of 1. NO COPY reproduces")
print("  a meaningful part of the purple gene's neighbourhood.")
mx = pw and max(pw, key=lambda r: r["overlap_coef"])
if mx:
    print(f"\n  Within sugarcane the copies are NOT independent: the most similar pair")
    print(f"  ({mx['gene_a']} / {mx['gene_b']}) shares {mx['overlap_coef']:.1%} of the")
    print("  smaller neighbourhood. Nine copies are not nine measurements.")
say("")
say("done")
