#!/usr/bin/env python3
# =============================================================================
# 69_og_conservation.py -- edge conservation as a DIRECT comparison of two graphs
# on the same vertices.
#
# WHY THIS IS NOT 62. 61 and 62 exist to project an edge through orthology: a
# sugarcane edge (i,j) is conserved if SOME purple ortholog of i is joined to SOME
# purple ortholog of j. That projection is many-to-many -- it is what inflated the
# conserved-pair count to 1.25x by multiplicity -- and it is unavoidable while the
# vertex is a gene, because the two species do not have the same genes.
#
# 68 removed the problem instead of correcting for it. Both orthogroup networks
# are built on ONE shared row set (63,271 groups, same order, asserted), so the
# projection is the IDENTITY and conservation is just set intersection.
#
# THE OLD NULL HAS NOTHING LEFT TO PERMUTE, AND THAT IS THE POINT. 61 and 62
# permute the ortholog assignment; here there is no assignment. The honest
# replacement is a DEGREE-MATCHED NODE RELABELLING of the target graph: each node
# is swapped for another of similar degree, so the null keeps both degree
# sequences and destroys only which group is which. Overlap above that is
# conservation that degree alone does not explain -- which is the confound
# docs/results.md flags as unseparated at gene level ("genes on a conserved edge
# are hubs, median degree 228 vs 18, partly by arithmetic").
#
# WHAT IS REPORTED
#   overlap, Jaccard             on the full edge sets, both directions at once
#                                (the comparison is symmetric now, so "sugarcane
#                                -> purple" and the reverse are one number)
#   rate by weight decile        rank-based deciles of |r|, as 62 does, so the
#                                strength trend stays comparable
#   degree correlation           Spearman of per-group degree across species: does
#                                a well-connected orthogroup in one species tend to
#                                be well-connected in the other?
#   fold over the null           with the degree-matched relabelling
#
# RUN: through run.sh -> ./run.sh ogconservation
# =============================================================================

import json
import os
import sys

import numpy as np
import pandas as pd

STUDIES = ("sugarcane", "purple")
CHUNK = 20_000_000


def env_req(name):
    v = os.environ.get(name, "")
    if not v:
        sys.exit(f"FATAL: {name} is not set -- launch through run.sh")
    return v


say = lambda *a: print(*a, flush=True)

OG_PREFIX = {s: env_req(f"CLEAN_OG_PREFIX_{s.upper()}") for s in STUDIES}
LAYER = {s: env_req(f"CLEAN_LAYER_{s.upper()}") for s in STUDIES}
OUT_DIR = env_req("CLEAN_OUT_DIR")
NULL_REPS = int(os.environ.get("CLEAN_NULL_REPS", "20"))
N_BINS = int(os.environ.get("CLEAN_WEIGHT_BINS", "10"))
DEG_STRATA = int(os.environ.get("CLEAN_DEGREE_STRATA", "20"))
SEED = int(os.environ.get("CLEAN_SEED", "1"))

os.makedirs(OUT_DIR, exist_ok=True)
rng = np.random.default_rng(SEED)
say("edge conservation on a shared vertex set\n")

# --- the shared vertex set, asserted -----------------------------------------
nodes = {}
for s in STUDIES:
    with open(f"{OG_PREFIX[s]}.genes.txt") as fh:
        nodes[s] = [ln.strip() for ln in fh if ln.strip()]
ref = nodes[STUDIES[0]]
for s in STUDIES[1:]:
    if nodes[s] != ref:
        sys.exit(f"FATAL: {s}'s row set differs from {STUDIES[0]}'s -- this stage "
                 f"requires the shared set 68 writes; re-run ./run.sh ogmatrix")
IDX = {g: i for i, g in enumerate(ref)}
N = len(ref)
say(f"shared vertices: {N:,}")

# Undirected edge as one int64: min<<32 | max. Sorted, it makes membership a
# searchsorted and the intersection an intersect1d -- the same trick 62 uses to
# afford 676M edges.
pack = lambda a, b: (np.minimum(a, b).astype(np.int64) << 32) | np.maximum(a, b).astype(np.int64)


def load_edges(study):
    """-> (sorted int64 edge keys, |r|), with a hard completeness check.

    THE COMPLETENESS CHECK IS NOT OPTIONAL. The engine streams its edgelist for
    minutes (8.7 for purple's 118,874,423 edges), and a reader that starts while
    it is still writing gets a SILENTLY TRUNCATED graph -- which happened on
    2026-10-01, when this stage read 90,032,320 of those edges, 76% of them, and
    produced conservation rates that looked entirely plausible. summary.json is
    written only after the edgelist is closed, so its n_significant_edges is the
    authority on how many rows there should be.
    """
    keys, wts = [], []
    f = LAYER[study]
    if not os.path.exists(f):
        sys.exit(f"FATAL: {f} missing -- run  ./run.sh ognetwork {study}")
    sj = f.replace(".edgelist.tsv", ".summary.json")
    if not os.path.exists(sj):
        sys.exit(f"FATAL: {sj} missing -- the layer is unfinished or was interrupted; "
                 f"re-run  ./run.sh ognetwork {study}")
    with open(sj) as fh:
        expect = int(json.load(fh)["n_significant_edges"])
    for ch in pd.read_csv(f, sep="\t", usecols=["gene1", "gene2", "pearson"],
                          chunksize=CHUNK):
        a = ch["gene1"].map(IDX).to_numpy()
        b = ch["gene2"].map(IDX).to_numpy()
        ok = np.isfinite(a.astype(float)) & np.isfinite(b.astype(float))
        if not ok.all():
            sys.exit(f"FATAL: {study} has {int((~ok).sum()):,} edges whose endpoints "
                     f"are not in the shared vertex set")
        keys.append(pack(a.astype(np.int64), b.astype(np.int64)))
        wts.append(np.abs(ch["pearson"].to_numpy(np.float32)))
    k = np.concatenate(keys); w = np.concatenate(wts)
    if len(k) != expect:
        sys.exit(f"FATAL: {study} edgelist holds {len(k):,} edges but its "
                 f"summary.json says {expect:,}. The file is incomplete -- most "
                 f"likely still being written. Wait for  ./run.sh ognetwork "
                 f"{study}  to finish, then re-run.")
    o = np.argsort(k, kind="stable")
    return k[o], w[o]


E, W, DEG = {}, {}, {}
for s in STUDIES:
    E[s], W[s] = load_edges(s)
    deg = np.bincount(np.concatenate([(E[s] >> 32), (E[s] & 0xFFFFFFFF)]),
                      minlength=N)
    say(f"{s:10s} {len(E[s]):>12,} edges | mean degree {2 * len(E[s]) / N:8.1f} "
        f"| |r| median {np.median(W[s]):.4f}")
    DEG[s] = deg

A, B = STUDIES
say("")

# --- observed overlap --------------------------------------------------------
shared_keys = np.intersect1d(E[A], E[B], assume_unique=False)
n_shared = len(shared_keys)
union = len(E[A]) + len(E[B]) - n_shared
say(f"edges in both      {n_shared:>12,}")
say(f"jaccard            {n_shared / union:>12.6f}")
say(f"of {A:<14s} {100 * n_shared / len(E[A]):>11.3f}%")
say(f"of {B:<14s} {100 * n_shared / len(E[B]):>11.3f}%")

# --- degree correlation ------------------------------------------------------
from scipy.stats import spearmanr                                   # noqa: E402
both_present = (DEG[A] > 0) & (DEG[B] > 0)
rho, p_rho = spearmanr(DEG[A][both_present], DEG[B][both_present])
say(f"\ndegree spearman    {rho:>12.4f}  (p = {p_rho:.3g}, "
     f"{int(both_present.sum()):,} groups with degree > 0 in both)")

# --- rate by weight decile, in the source species ---------------------------
rows = []
for src, tgt in ((A, B), (B, A)):
    in_t = np.isin(E[src], E[tgt], assume_unique=False)
    q = np.quantile(W[src], np.linspace(0, 1, N_BINS + 1))
    q[0], q[-1] = -np.inf, np.inf
    binid = np.clip(np.searchsorted(q, W[src], side="right") - 1, 0, N_BINS - 1)
    for b in range(N_BINS):
        m = binid == b
        if not m.any():
            continue
        rows.append(dict(source=src, target=tgt, decile=b + 1,
                         n_edges=int(m.sum()), n_conserved=int(in_t[m].sum()),
                         rate=float(in_t[m].mean()),
                         r_min=float(W[src][m].min()), r_max=float(W[src][m].max())))
    say(f"\n{src} -> {tgt}: conserved by |r| decile")
    for d in rows[-N_BINS:]:
        say(f"  D{d['decile']:<2d} |r| {d['r_min']:.4f}-{d['r_max']:.4f}  "
            f"{d['n_edges']:>11,} edges  {100 * d['rate']:>7.3f}% conserved")
pd.DataFrame(rows).to_csv(f"{OUT_DIR}/og_conservation_by_decile.tsv",
                          sep="\t", index=False)

# --- the degree-matched null -------------------------------------------------
# Relabel B's vertices within strata of B's own degree, so the null keeps both
# degree sequences exactly and destroys only WHICH orthogroup is which.
say(f"\nnull: {NULL_REPS} degree-matched relabellings of {B} "
    f"({DEG_STRATA} degree strata)")
order = np.argsort(DEG[B], kind="stable")
strata = np.array_split(order, DEG_STRATA)
null_counts = np.empty(NULL_REPS, dtype=np.int64)
src_b = (E[B] >> 32).astype(np.int64)
dst_b = (E[B] & 0xFFFFFFFF).astype(np.int64)
for rep in range(NULL_REPS):
    perm = np.arange(N, dtype=np.int64)
    for st in strata:
        perm[st] = st[rng.permutation(len(st))]
    pk = np.sort(pack(perm[src_b], perm[dst_b]))
    null_counts[rep] = len(np.intersect1d(E[A], pk, assume_unique=False))
nm = float(null_counts.mean())
fold = n_shared / nm if nm > 0 else float("inf")
p_emp = (int((null_counts >= n_shared).sum()) + 1) / (NULL_REPS + 1)
say(f"  observed {n_shared:,} | null mean {nm:,.1f} "
    f"(sd {null_counts.std():,.1f}) | fold {fold:.3f} | p {p_emp:.4g}")

# --- write -------------------------------------------------------------------
summ = pd.DataFrame(dict(metric=[
    "n_vertices", "n_edges_" + A, "n_edges_" + B, "mean_degree_" + A,
    "mean_degree_" + B, "n_edges_shared", "jaccard",
    "pct_of_" + A, "pct_of_" + B, "degree_spearman", "degree_spearman_p",
    "n_both_degree_gt0", "null_reps", "degree_strata", "null_mean", "null_sd",
    "fold_over_null", "p_empirical"],
    value=[N, len(E[A]), len(E[B]), f"{2 * len(E[A]) / N:.2f}",
           f"{2 * len(E[B]) / N:.2f}", n_shared, f"{n_shared / union:.6f}",
           f"{100 * n_shared / len(E[A]):.4f}", f"{100 * n_shared / len(E[B]):.4f}",
           f"{rho:.4f}", f"{p_rho:.3g}", int(both_present.sum()), NULL_REPS,
           DEG_STRATA, f"{nm:.1f}", f"{null_counts.std():.1f}",
           f"{fold:.4f}", f"{p_emp:.4g}"]))
summ.to_csv(f"{OUT_DIR}/og_conservation_summary.tsv", sep="\t", index=False)
pd.DataFrame(dict(rep=np.arange(1, NULL_REPS + 1), n_conserved=null_counts)).to_csv(
    f"{OUT_DIR}/og_conservation_null.tsv", sep="\t", index=False)
pd.DataFrame(dict(og_2sp=ref, degree_sugarcane=DEG["sugarcane"],
                  degree_purple=DEG["purple"])).to_csv(
    f"{OUT_DIR}/og_degree_both.tsv", sep="\t", index=False)
say(f"\nwrote og_conservation_summary.tsv, _by_decile.tsv, _null.tsv, og_degree_both.tsv")
say("done")
