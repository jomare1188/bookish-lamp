#!/usr/bin/env python3
# =============================================================================
# 10_myb59_annotate.py -- the AtMYB59 seeds and their neighbours, annotated.
#
# WHAT THIS IS FOR. Step 06 wrote eleven files of bare gene IDs. That is enough for
# a script and useless to a bench biologist, who needs to know what each gene IS.
# This turns each list into a table: how strongly the gene is tied to the MYB, how
# connected it is, which module it sits in, whether it is a TF, whether it responds
# to nitrogen, what GO says it does (split BP / MF / CC, as NAMES not IDs), what
# InterPro found in it, and what its counterpart is in the other species.
#
# EVERY NEIGHBOUR IS KEPT HERE. The Cytoscape step (11) draws only the |r| >= 0.85
# core, because the full graph is a hairball -- but that is a drawing decision. The
# tables below are complete, so nothing is lost and the trimmed view can always be
# checked against them.
#
# ---------------------------------------------------------------------------
# TWO JOINS CAN SILENTLY MULTIPLY ROWS, and both are aggregated BEFORE joining:
# a gene has many GO terms and many InterPro hits, so a naive merge turns one gene
# into forty rows and every downstream count is then wrong in a way that still looks
# like a table. The row count is asserted against the gene count afterwards.
# ---------------------------------------------------------------------------
#
# |r| IS RECOMPUTED, NOT STORED. The .pairs dump carries the rescaled weight, and
# both networks use the SAME scale (|r| 0.8 -> 0.01, 0.9999 -> 1, from config.sh's
# STAT_MIN/STAT_MAX), so one inverse formula serves both species and a threshold
# means the same thing in either. Spot-checked against network_<study>_edges.tsv,
# which carries `weight` and `pearson_r` side by side.
#
# RUN: ./run_all.sh 10
# =============================================================================

import csv
import glob
import os
import re
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

BASE = "/dados04/jorge/comparative_saccharum"
CLEAN = f"{BASE}/new_clean"
WORK = "/dados04/jorge/tmp/mcl_work_cluster"
RES = f"{CLEAN}/results"
M20 = f"{RES}/readouts/module20"
NBR = f"{M20}/myb59_neighbourhoods"
OUT = f"{M20}/myb59_annotated"
GO_TERMS = f"{CLEAN}/docs/data/go_terms.tsv"
OG_FILE = (f"{BASE}/files/fix_orthofinder/proteins/OrthoFinder/"
           f"Results_Jun04_2/Orthogroups/Orthogroups.tsv")
SP_COL = {"sugarcane": "sugarcane_one_transcript",
          "purple": "one_transcript_purple_proteins"}
# The purple side is the focus gene only; the sugarcane side is all nine AtMYB59
# copies. Both are what the collaborator asked to explore.
SEEDS = {"purple": ["Soffic.09G0001580-9H"], "sugarcane": None}   # None = all files
EXPECT_EDGES = {"sugarcane": 75_333_769, "purple": 675_955_918}
LO, HI = 0.8, 0.9999
CORE_MIN_R = 0.85          # what step 11 draws; also what gets cached here

os.makedirs(OUT, exist_ok=True)
say = lambda *a: print(*a, flush=True)
r_of = lambda w: LO + (np.asarray(w) - 0.01) / 0.99 * (HI - LO)
_RE_P, _RE_V = re.compile(r"\.p[0-9]+$"), re.compile(r"\.v[0-9]+(\.[0-9]+)*$")
strip_version = lambda g: _RE_V.sub("", _RE_P.sub("", g))

say("annotating the AtMYB59 seeds and their neighbours\n")

# --- GO id -> (name, ontology) ------------------------------------------------
if not os.path.exists(GO_TERMS):
    sys.exit(f"FATAL: {GO_TERMS} missing -- run step 09 first")
gt = pd.read_csv(GO_TERMS, sep="\t")
GO_NAME = dict(zip(gt["GO.ID"], gt["Term"]))
GO_ONT = dict(zip(gt["GO.ID"], gt["Ontology"]))
say(f"GO lookup: {len(GO_NAME):,} terms")


def read_tab(study):
    names = []
    with open(f"{WORK}/{study}.tab") as fh:
        for expect, line in enumerate(fh):
            i, nm = line.rstrip("\n").split("\t", 1)
            if int(i) != expect:
                sys.exit(f"FATAL: {study}.tab not in index order at {expect}")
            names.append(nm)
    return names


# --- orthology, in-network only (the universe 61/62 use) ----------------------
def ortholog_map(idx_a, idx_b, col_a, col_b):
    m = defaultdict(set)
    with open(OG_FILE) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            A = [strip_version(x.strip()) for x in (row[col_a] or "").split(",") if x.strip()]
            B = [strip_version(x.strip()) for x in (row[col_b] or "").split(",") if x.strip()]
            A = [g for g in A if g in idx_a]
            B = [g for g in B if g in idx_b]
            for g in A:
                m[g].update(B)
    return m


def aggregate_go(path, wanted):
    """gene -> {BP: [names], MF: [...], CC: [...]}, aggregated BEFORE any join."""
    out = defaultdict(lambda: defaultdict(set))
    with open(path) as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for row in rdr:
            g = row["gene"]
            if g not in wanted:
                continue
            gid = row["go_id"]
            ont = GO_ONT.get(gid)
            if ont:
                out[g][ont].add(GO_NAME.get(gid, gid))
    return out


def aggregate_ipr(path, wanted):
    """gene -> (accessions, descriptions). InterProScan TSV: col 1 protein,
    col 12 IPR accession, col 13 IPR description; '-' where a signature has none."""
    acc = defaultdict(list)
    desc = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 13:
                continue
            g = p[0]
            if g not in wanted or p[11] in ("-", ""):
                continue
            if p[11] not in acc[g]:
                acc[g].append(p[11])
                desc[g].append(p[12])
    return acc, desc


# =============================================================================
# THE SIGN OF r, AND WHY IT COSTS A 64 GB SCAN
#
# The .pairs dump carries only the rescaled WEIGHT, which is built from |r| -- the
# sign is gone. For a transcription factor that is close to the most important thing
# about a partner: a gene that falls as the MYB rises is a candidate target of
# repression, and one that rises with it is not, and both look identical in |r|.
# Measured on 06Ag012300's neighbours: 46% are NEGATIVELY correlated. Shipping a
# table that silently merged the two would misinform the reader.
#
# The sign survives only in network_<study>_edges.tsv (7.7 GB and 64.5 GB). Those
# tables describe the MERGED network, but filtering to `source != "mi"` recovers
# EXACTLY the Pearson-only edge set -- verified by arithmetic: 73,329,192 + 2,004,577
# = 75,333,769 and 625,775,930 + 50,179,988 = 675,955,918, the two network edge
# counts to the edge. So the scan is expensive but exact.
#
# ONE PASS serves both deliverables: every edge touching a seed (for the tables
# here), plus every strong edge among the neighbourhood (for the Cytoscape step),
# written to a small file so step 11 need not scan again.
# =============================================================================
def scan_signed_edges(study, seeds, union, core_min_r):
    """-> (seed_edges {seed: {gene: signed r}}, strong_edges [(g1,g2,r)])"""
    seeds, union = set(seeds), set(union)
    # .isin() wants a list-like; building it once avoids a rebuild per chunk
    seeds_l, union_l = list(seeds), list(union)
    se = defaultdict(dict)
    strong = []
    n_pearson = 0
    f = f"{RES}/{study}/network_{study}_edges.tsv"
    if not os.path.exists(f):
        sys.exit(f"FATAL: {f} is the only source of the SIGN of r, and it is missing")
    for ch in pd.read_csv(f, sep="\t", usecols=["gene1", "gene2", "pearson_r", "source"],
                          dtype={"gene1": str, "gene2": str,
                                 "pearson_r": np.float32, "source": str},
                          chunksize=4_000_000, engine="c"):
        ch = ch[ch["source"].to_numpy() != "mi"]          # Pearson-only edge set
        n_pearson += len(ch)
        # pandas .isin() is vectorised in C; the obvious
        # np.fromiter((g in union for g in g1), ...) is a Python loop over 4M
        # strings per chunk and turns a 4-minute scan into a 40-minute one.
        in1 = ch["gene1"].isin(union_l).to_numpy()
        in2 = ch["gene2"].isin(union_l).to_numpy()
        both = in1 & in2
        if not both.any():
            continue
        sub = ch[both]
        g1 = sub["gene1"].to_numpy(); g2 = sub["gene2"].to_numpy()
        rr = sub["pearson_r"].to_numpy()
        s1 = sub["gene1"].isin(seeds_l).to_numpy()
        s2 = sub["gene2"].isin(seeds_l).to_numpy()
        in1 = in2 = np.ones(len(g1), bool)
        for k in np.flatnonzero((s1 | s2) & (in1 & in2)):
            a, b, r = g1[k], g2[k], float(rr[k])
            if a in seeds:
                se[a][b] = r
            if b in seeds:
                se[b][a] = r
        keep = (in1 & in2) & (np.abs(rr) >= core_min_r)
        for k in np.flatnonzero(keep):
            strong.append((g1[k], g2[k], float(rr[k])))
    if n_pearson != EXPECT_EDGES[study]:
        sys.exit(f"FATAL: {n_pearson:,} Pearson edges in the merged table, expected "
                 f"{EXPECT_EDGES[study]:,} -- the source filter is wrong")
    return se, strong


# =============================================================================
# per species
# =============================================================================
NAMES, IDX, TABLES = {}, {}, {}
for sp in ("purple", "sugarcane"):
    NAMES[sp] = read_tab(sp)
    IDX[sp] = {g: i for i, g in enumerate(NAMES[sp])}
ORTH = {
    "sugarcane": ortholog_map(IDX["sugarcane"], IDX["purple"],
                              SP_COL["sugarcane"], SP_COL["purple"]),
    "purple": ortholog_map(IDX["purple"], IDX["sugarcane"],
                           SP_COL["purple"], SP_COL["sugarcane"]),
}

for sp in ("purple", "sugarcane"):
    names, idx = NAMES[sp], IDX[sp]
    files = sorted(glob.glob(f"{NBR}/{sp}__*.txt"))
    seeds = SEEDS[sp] or [os.path.basename(f)[:-4].split("__")[1] for f in files]
    seeds = [g for g in seeds if g in idx]
    say(f"\n=== {sp}: {len(seeds)} seed(s)")

    nbr_of = {}
    for g in seeds:
        f = f"{NBR}/{sp}__{g}.txt"
        nbr_of[g] = [x for x in open(f).read().split() if x]

    genes = sorted(set(seeds) | {x for v in nbr_of.values() for x in v})
    say(f"    {len(genes):,} genes (seeds + all neighbours)")

    # --- signed r to each seed, and the degree check ---------------------------
    say(f"    scanning {sp}'s edge table for the SIGN of r (this is the slow step)")
    se, strong = scan_signed_edges(sp, seeds, genes, CORE_MIN_R)
    for g in seeds:
        got = len(se[g])
        if got != len(nbr_of[g]):
            sys.exit(f"FATAL: {g} has {got} edges but {len(nbr_of[g])} neighbours listed")
    say(f"    every seed's edge count matches its neighbour list; "
        f"{len(strong):,} edges at |r| >= {CORE_MIN_R} cached for step 11")
    neg = sum(1 for d in se.values() for r in d.values() if r < 0)
    tot = sum(len(d) for d in se.values())
    say(f"    {neg:,} of {tot:,} seed edges are NEGATIVE correlations "
        f"({100 * neg / tot:.0f}%)")
    with open(f"{M20}/myb59_signed_edges_{sp}.tsv", "w") as fh:
        fh.write("gene1\tgene2\tr\n")
        for a, b, r in strong:
            fh.write(f"{a}\t{b}\t{r:.4f}\n")

    # --- per-gene attributes --------------------------------------------------
    nm = pd.read_csv(f"{RES}/{sp}/network_{sp}_node_metrics.tsv", sep="\t",
                     usecols=["gene", "degree", "strength"])
    deg = dict(zip(nm["gene"], nm["degree"].astype(int)))
    stg = dict(zip(nm["gene"], nm["strength"]))
    srt = np.sort(nm["degree"].to_numpy())
    pct = lambda d: 100.0 * np.searchsorted(srt, d, "left") / len(srt)

    mem = pd.read_csv(f"{RES}/{sp}/mcl_{sp}_membership.tsv", sep="\t",
                      usecols=["gene", "module_name"])
    mod = dict(zip(mem["gene"], mem["module_name"]))
    msize = mem["module_name"].value_counts().to_dict()

    tf_f = f"{RES}/readouts/get_tfs/{sp}/TF_in_network.tsv"
    tf = {}
    if os.path.exists(tf_f):
        t = pd.read_csv(tf_f, sep="\t")
        tf = dict(zip(t["gene"], t["Family"]))

    bl = pd.read_csv(f"{RES}/{sp}/gene_trait_blocked_{sp}.tsv", sep="\t",
                     usecols=["gene", "r_partial", "padj", "responsive"])
    nr = dict(zip(bl["gene"], bl["r_partial"]))
    npadj = dict(zip(bl["gene"], bl["padj"]))
    nresp = dict(zip(bl["gene"], bl["responsive"]))

    want = set(genes)
    go = aggregate_go(f"{BASE}/annotation/{sp}/gene2go_{sp}.tsv", want)
    ipr_f = f"{BASE}/annotation/{sp}/interproscan_full/{sp}.interproscan_full.tsv"
    ipr_acc, ipr_desc = ({}, {})
    if os.path.exists(ipr_f):
        ipr_acc, ipr_desc = aggregate_ipr(ipr_f, want)
        say(f"    InterPro: {len(ipr_acc):,} of {len(genes):,} genes carry an IPR hit")
    else:
        say(f"    NOTE: no InterProScan table at {ipr_f} -- IPR columns left empty")

    J = lambda xs: "; ".join(sorted(xs)) if xs else ""
    rows = []
    for g in genes:
        is_seed = g in seeds
        # which seeds is it a neighbour of, and how strongly
        of, rs = [], []
        for s in seeds:
            rv = se[s].get(g)
            if rv is not None:
                of.append(s)
                rs.append(f"{rv:.4f}")
        d = deg.get(g, 0)
        rows.append({
            "gene": g,
            "role": "seed" if is_seed else "neighbour",
            "neighbour_of": "; ".join(of),
            "n_seeds": len(of),
            "r_to_seed": "; ".join(rs),
            "r_to_seed_strongest": max((float(x) for x in rs), key=abs, default=""),
            "direction": ("" if is_seed else
                          ("positive" if rs and max((float(x) for x in rs), key=abs) > 0
                           else "negative")),
            "degree": d,
            "degree_percentile": round(float(pct(d)), 1),
            "strength": round(float(stg.get(g, float("nan"))), 3),
            "mcl_module": mod.get(g, ""),
            "module_size": msize.get(mod.get(g, ""), ""),
            "is_TF": "yes" if g in tf else "no",
            "TF_family": tf.get(g, ""),
            "nitrogen_r": round(float(nr[g]), 4) if g in nr else "",
            "nitrogen_padj": f"{npadj[g]:.3g}" if g in npadj else "",
            "nitrogen_responsive": ("yes" if str(nresp.get(g)) in ("True", "TRUE")
                                    else "no") if g in nresp else "",
            "GO_BP": J(go.get(g, {}).get("BP", set())),
            "GO_MF": J(go.get(g, {}).get("MF", set())),
            "GO_CC": J(go.get(g, {}).get("CC", set())),
            "InterPro": "; ".join(ipr_acc.get(g, [])),
            "InterPro_desc": "; ".join(ipr_desc.get(g, [])),
            "ortholog_other_species": J(ORTH[sp].get(g, set())),
        })

    df = pd.DataFrame(rows)
    # ONE ROW PER GENE. The GO and InterPro joins are the two that could multiply
    # rows; they were aggregated first, and this is the proof.
    if len(df) != len(genes) or df["gene"].duplicated().any():
        sys.exit(f"FATAL: {len(df)} rows for {len(genes)} genes -- a join multiplied rows")
    df["_abs"] = df["r_to_seed_strongest"].apply(lambda x: abs(x) if x != "" else 0)
    df = df.sort_values(["role", "_abs"], ascending=[True, False]).drop(columns="_abs")

    comb = f"{OUT}/{sp}__ALL_myb59_neighbours.tsv"
    df.to_csv(comb, sep="\t", index=False)
    say(f"    wrote {os.path.basename(comb)}  ({len(df):,} rows)")

    for s in seeds:
        sub = df[(df["gene"] == s) | (df["neighbour_of"].str.contains(re.escape(s)))].copy()
        # IN A PER-SEED FILE, r_to_seed MUST MEAN THIS SEED. A gene neighbouring two
        # of the nine copies carries two values in the combined table ("0.8094;
        # 0.9427"), and left as-is a reader cannot tell which belongs to the seed the
        # file is named after -- or parse the column as a number at all.
        sub.insert(3, "r_to_this_seed",
                   [("" if g == s else f"{se[s][g]:.4f}") for g in sub["gene"]])
        sub.insert(4, "direction_to_this_seed",
                   ["" if g == s else ("positive" if se[s][g] > 0 else "negative")
                    for g in sub["gene"]])
        sub = sub.assign(_a=[0 if g == s else abs(se[s][g]) for g in sub["gene"]]) \
                 .sort_values("_a", ascending=False).drop(columns="_a")
        f = f"{OUT}/{sp}__{s}.tsv"
        sub.to_csv(f, sep="\t", index=False)
        say(f"      {os.path.basename(f):52s} {len(sub):>6,} rows")
    TABLES[sp] = df

say("\nannotation coverage of the genes in these tables:")
for sp, df in TABLES.items():
    n = len(df)
    g = int((df["GO_BP"].astype(bool) | df["GO_MF"].astype(bool)
             | df["GO_CC"].astype(bool)).sum())
    i = int(df["InterPro"].astype(bool).sum())
    t = int((df["is_TF"] == "yes").sum())
    r = int((df["nitrogen_responsive"] == "yes").sum())
    say(f"  {sp:10s} {n:>6,} genes | GO {g:>6,} ({100*g/n:.1f}%) | "
        f"InterPro {i:>6,} ({100*i/n:.1f}%) | TF {t:>4,} | N-responsive {r:>4,}")
say("\ndone")
