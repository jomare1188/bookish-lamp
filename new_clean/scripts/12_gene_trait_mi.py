#!/usr/bin/env python
# ============================================================================
# 12_gene_trait_mi.py — mutual information between each gene and a TRAIT
#
# WHY THIS EXISTS
#
# 07_gene_trait_cor.r correlates each gene against the trait with Pearson. For
# purple that trait is a three-level nitrogen DOSE (0/2/6 mM), and a linear
# correlation against a dose only sees a monotone, roughly linear trend. A
# saturating or threshold response -- the shape a nitrogen response is most
# likely to take -- is invisible to it. That is the same argument that motivated
# the MI layer for the network, applied to the gene-vs-trait question instead of
# gene-vs-gene.
#
# It matters here more than anywhere else in the pipeline: the node-level
# conservation result bottlenecks on purple having only 62 trait-correlated
# genes at n = 18, and the single most nitrogen-responsive-looking gene found so
# far (a MYB61 copy, |r| = 0.643) fails FDR because a linear r of that size is
# not significant at that sample size.
#
# ---------------------------------------------------------------------------
# THE ESTIMATOR
#
# The trait is DISCRETE (2 or 3 levels) and expression is continuous, so the
# KSG algorithm the network layer uses does not apply -- it assumes both
# variables are continuous. The right tool is Ross (2014, PLoS ONE 9:e87357),
# which is KSG adapted to a discrete class variable:
#
#     I = psi(N) - <psi(N_c)> + <psi(k_i)> - <psi(m_i)>
#
# where for each sample i:
#     N_c  = size of i's trait class
#     d_i  = distance to the k-th nearest neighbour OF THE SAME CLASS
#     k_i  = same-class neighbours within d_i
#     m_i  = neighbours of ANY class within d_i
#
# TIE CORRECTION, and why it is not optional here. Ross assumes continuous x, so
# distances are distinct. Expression is rank-transformed (see below), which makes
# every distance a small integer and ties are then the norm rather than the
# exception. Using a fixed k against tied distances underestimates the
# neighbourhood and the estimator goes NEGATIVE for perfectly separated classes.
# Counting the actual number of same-class points within d_i (k_i, >= k) fixes
# it: on perfectly separated classes the estimator then returns 0.657 at n = 48,
# C = 2 against the ceiling of ln 2 = 0.693, and 0.908 at n = 18, C = 3 against
# ln 3 = 1.099 -- approaching the entropy of the trait from below, as it must.
#
# ---------------------------------------------------------------------------
# WHY THE NULL IS SHARED, AND WHY THAT MAKES THIS CHEAP
#
# Each gene is rank-transformed, ties broken at random, so its expression
# becomes the positions 0..n-1. The statistic then depends ONLY on which trait
# label sits at which position -- i.e. on a length-n string over the trait's
# alphabet -- and not on the gene's actual values at all.
#
# So the null distribution is the distribution of that statistic over random
# permutations of the trait labels, and it is the SAME for every gene. One
# permutation run serves all 40,000 of them, exactly as the shared rank-
# permutation null serves all 1.5e10 gene pairs in 02_network_engine.py.
#
# It also means the null is cheap enough to resolve the threshold that actually
# matters. BH over ~40,000 genes needs p ~ 1.3e-6 for the strongest gene; 1e7
# permutations put ~13 draws past that, so the cut comes from the empirical null
# rather than an extrapolation. (The GPD tail is still fitted, for genes beyond
# the empirical range.)
#
# ---------------------------------------------------------------------------
# OUTPUT
#
#   <out>.tsv          gene trait mi pval padj pearson pearson_pval pearson_padj
#                      finding    -- which method calls the gene significant
#   <out>.null.tsv     the shared null: quantiles, mean/sd, GPD tail fit
#   <out>.summary.json counts, thresholds, timings
#
# `finding` is the point of the exercise: `mi_only` genes are those a linear
# correlation cannot see.
#
# RUN (through run.sh):  ./run.sh traitmi sugarcane
# ============================================================================

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import os
import time

import numpy as np
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "engine", os.path.join(HERE, "02_network_engine.py"))
engine = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(engine)

P_FLOOR = engine.P_FLOOR


# ---------------------------------------------------------------------------
def parse_traits(spec: str) -> dict:
    """'genotype:A=1,B=0;treatment:0N=0,2N=2' -> {trait: {level: value}}"""
    out = {}
    for block in spec.split(";"):
        block = block.strip()
        if not block:
            continue
        name, body = block.split(":", 1)
        levels = {}
        for kv in body.split(","):
            lv, val = kv.rsplit("=", 1)
            levels[lv.strip()] = float(val)
        out[name.strip()] = levels
    return out


def read_meta(path: str):
    """Sample sheet -> (list of column names lowercased, list of row dicts)."""
    import csv
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise SystemExit(f"{path} is empty")
    rows = [{(k or "").strip().lower(): v for k, v in r.items()} for r in rows]
    return rows


# ---------------------------------------------------------------------------
def ross_batch(L: torch.Tensor, class_sizes: torch.Tensor, k: int) -> torch.Tensor:
    """Ross MI with tie correction, for a batch of label sequences.

    L is (B, n) of class ids, where L[b, r] is the class of the sample sitting
    at rank position r. Distances are |r - r'|, so the whole computation is a
    function of L alone -- which is what makes the null shared.
    """
    B, n = L.shape
    dev = L.device
    pos = torch.arange(n, device=dev)
    D = (pos[:, None] - pos[None, :]).abs()                    # (n, n)
    notself = ~torch.eye(n, dtype=torch.bool, device=dev)

    same = L[:, :, None] == L[:, None, :]                      # (B, n, n)
    valid = same & notself

    BIG = n + 1
    md = torch.where(valid, D.expand(B, n, n), torch.full_like(D.expand(B, n, n), BIG))
    # k-th smallest same-class distance
    d = md.sort(dim=2).values[:, :, k - 1]                     # (B, n)

    within = D.expand(B, n, n) <= d[:, :, None]
    k_i = (within & valid).sum(2)                              # same-class within d
    m_i = (within & notself).sum(2)                            # any class within d

    # size of each sample's own class
    Nc = class_sizes.to(dev)[L]                                # (B, n)

    dg = torch.digamma
    return (dg(torch.tensor(float(n), device=dev))
            - dg(Nc.float()).mean(1)
            + dg(k_i.float()).mean(1)
            - dg(m_i.float()).mean(1))


def build_null(lab: np.ndarray, class_sizes: np.ndarray, k: int, n_perm: int,
               device, seed: int, chunk: int = 20000) -> np.ndarray:
    """Statistic under label permutation. Shared by every gene."""
    rng = np.random.default_rng(seed + 977)
    n = len(lab)
    cs = torch.from_numpy(class_sizes).to(device)
    out, done, t0 = [], 0, time.time()
    while done < n_perm:
        m = min(chunk, n_perm - done)
        # rng.permuted shuffles each row independently
        P = rng.permuted(np.broadcast_to(lab, (m, n)).copy(), axis=1)
        L = torch.from_numpy(P).to(device).long()
        out.append(ross_batch(L, cs, k).float().cpu().numpy())
        done += m
    null = np.sort(np.concatenate(out))
    print(f"[null] {len(null):,} label permutations in {time.time()-t0:.1f}s  "
          f"mean {null.mean():+.5f}  sd {null.std():.5f}  "
          f"q99.9 {np.quantile(null, 0.999):.5f}  max {null[-1]:.5f}", flush=True)
    return null


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Gene-vs-trait mutual information (Ross estimator)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--matrix", required=True, help="VST prefix from 01_export_vst.r")
    ap.add_argument("--meta", required=True, help="sample sheet CSV")
    ap.add_argument("--traits", required=True, help="trait spec, as in config.sh")
    ap.add_argument("--trait", required=True, help="which trait to test")
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--k", type=int, default=3)
    ap.add_argument("--n-perm", type=int, default=10_000_000)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--seed", type=int, default=1188)
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--gene-chunk", type=int, default=20000)
    args = ap.parse_args()

    t_start = time.time()
    device = torch.device(args.device)
    if device.type == "cuda":
        if not torch.cuda.is_available():
            print("[warn] CUDA unavailable — falling back to CPU")
            device = torch.device("cpu")
        else:
            print(f"[gpu ] {torch.cuda.get_device_name(device)}")

    # ---- data ------------------------------------------------------------
    x, genes, meta = engine.load_matrix(args.matrix)
    n_genes, n_all = x.shape
    samples = meta["samples"]
    print(f"[data] {meta.get('label','?')}: {n_genes:,} genes x {n_all} libraries")

    traits = parse_traits(args.traits)
    if args.trait not in traits:
        raise SystemExit(f"--trait {args.trait} not in the spec "
                         f"({', '.join(traits)})")
    enc = traits[args.trait]

    rows = read_meta(args.meta)
    by_sample = {r["sample"]: r for r in rows if "sample" in r}
    missing = [s for s in samples if s not in by_sample]
    if missing:
        raise SystemExit(f"samples missing from {os.path.basename(args.meta)}: "
                         f"{missing[:5]}")
    if args.trait not in rows[0]:
        raise SystemExit(f"no '{args.trait}' column in "
                         f"{os.path.basename(args.meta)}; have "
                         f"{', '.join(sorted(rows[0]))}")

    raw = [by_sample[s][args.trait] for s in samples]
    keep = np.array([v in enc for v in raw])
    if not keep.all():
        dropped = sorted({v for v, kk in zip(raw, keep) if not kk})
        print(f"[trait] dropping {int((~keep).sum())} samples with unencoded "
              f"level(s): {dropped}")
    x = x[:, keep]
    vals = np.array([enc[v] for v, kk in zip(raw, keep) if kk], dtype=np.float64)
    n = int(keep.sum())

    levels, lab = np.unique(vals, return_inverse=True)
    class_sizes = np.bincount(lab).astype(np.int64)
    print(f"[trait] {args.trait}: n = {n}, {len(levels)} levels "
          + ", ".join(f"{lv:g} (n={cnt})" for lv, cnt in zip(levels, class_sizes)))
    if class_sizes.min() <= args.k:
        raise SystemExit(f"--k {args.k} needs every class to have more than k "
                         f"members; smallest class has {class_sizes.min()}")
    H = float(-(class_sizes / n * np.log(class_sizes / n)).sum())
    print(f"[trait] H(trait) = {H:.4f} nats — the ceiling on any MI with it")

    # ---- rank transform --------------------------------------------------
    ranks, tie_rate = engine.rank_transform(x, args.seed)
    print(f"[rank] ties broken at random in {tie_rate*100:.1f}% of genes "
          f"(sampled); each gene is now a permutation of 0..{n-1}")

    # order[g, r] = index of the sample sitting at rank r for gene g
    order = np.argsort(ranks, axis=1, kind="stable")
    lab_t = torch.from_numpy(lab).to(device).long()
    cs_t = torch.from_numpy(class_sizes).to(device)

    # ---- the shared null -------------------------------------------------
    null = build_null(lab, class_sizes, args.k, args.n_perm, device, args.seed)
    gpd = engine.fit_gpd_tail(null)
    pval_fn = engine.make_pvalue_fn(null, gpd)
    print(f"[null] tail fit: {gpd['method']}  u={gpd['u']:.5f} "
          f"xi={gpd['xi']:+.4f} sigma={gpd['sigma']:.5f}")
    depth = 1.0 / (args.n_perm * (args.alpha / n_genes))
    print(f"[null] BH needs p ~ {args.alpha/n_genes:.3g} for the top gene; "
          f"{args.n_perm:,} permutations resolve it "
          + ("empirically" if depth <= 1 else f"only {1/depth:.1f}x over"))

    # ---- per-gene statistic ----------------------------------------------
    t0 = time.time()
    mi = np.empty(n_genes, dtype=np.float64)
    for s in range(0, n_genes, args.gene_chunk):
        e = min(s + args.gene_chunk, n_genes)
        L = lab_t[torch.from_numpy(order[s:e]).to(device).long()]
        mi[s:e] = ross_batch(L, cs_t, args.k).double().cpu().numpy()
    print(f"[mi  ] {n_genes:,} genes in {time.time()-t0:.1f}s")

    p_mi = pval_fn(mi)
    padj_mi = bh(p_mi)

    # How much of the signal is EXTRAPOLATED? The permutation null has a finite
    # range and the statistic has a hard ceiling at H(trait), so any gene above
    # the null's empirical maximum gets a p-value from the GPD tail rather than
    # from data -- and past the fitted endpoint it is floored outright. Those
    # p-values order genes correctly but their magnitude is not calibrated:
    # reporting padj = 1e-20 for a gene whose permutation p is merely "< 1e-7"
    # would be indefensible, so the counts are surfaced here and carried in the
    # summary and the output table.
    null_max = float(null[-1])
    endpoint = engine.gpd_endpoint(gpd)
    beyond_null = int((mi > null_max).sum())
    beyond_fit = int((mi > endpoint).sum()) if np.isfinite(endpoint) else 0
    p_resolution = 1.0 / args.n_perm
    print(f"[null] empirical max {null_max:.5f} (p = {p_resolution:.1g}); "
          f"GPD endpoint {endpoint:.5f}")
    if beyond_null:
        print(f"[warn] {beyond_null:,} genes exceed the empirical null max. Their "
              f"p-values come from the GPD tail, not from permutations: treat "
              f"them as 'p < {p_resolution:.1g}' and use the ranking, not the "
              f"magnitude.")
    if beyond_fit:
        print(f"[warn] {beyond_fit:,} of those also exceed the fitted tail's "
              f"endpoint and are floored at p = {P_FLOOR:g}.")
    extrapolated = mi > null_max

    # ---- Pearson, on the same samples, for a like-for-like comparison ----
    from scipy import stats as st
    xc = x - x.mean(1, keepdims=True)
    vc = vals - vals.mean()
    denom = np.sqrt((xc ** 2).sum(1) * (vc ** 2).sum())
    r = np.where(denom > 0, (xc * vc).sum(1) / np.maximum(denom, 1e-300), 0.0)
    df = n - 2
    tstat = np.abs(r) * np.sqrt(df / np.maximum(1 - r * r, 1e-300))
    p_r = np.clip(2 * st.t.sf(tstat, df=df), P_FLOOR, 1.0)
    padj_r = bh(p_r)

    sig_mi = padj_mi <= args.alpha
    sig_r = padj_r <= args.alpha
    finding = np.where(sig_mi & sig_r, "both",
               np.where(sig_mi, "mi_only",
                np.where(sig_r, "pearson_only", "neither")))
    p_note = np.where(extrapolated, "extrapolated", "empirical")

    # ---- write -----------------------------------------------------------
    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    fmt = "%s\t%s\t%.5f\t%.4g\t%.4g\t%.5f\t%.4g\t%.4g\t%s\t%s\n".__mod__
    with open(args.out + ".tsv", "w", buffering=1 << 22) as fh:
        fh.write("gene\ttrait\tmi\tpval\tpadj\tpearson\tpearson_pval\t"
                 "pearson_padj\tfinding\tp_source\n")
        fh.write("".join(map(fmt, zip(genes, [args.trait] * n_genes, mi,
                                      p_mi, padj_mi, r, p_r, padj_r, finding,
                                      p_note))))

    with open(args.out + ".null.tsv", "w") as fh:
        fh.write("quantile\tvalue\n")
        for q in (0.5, 0.9, 0.99, 0.999, 0.9999, 0.99999, 1.0):
            fh.write(f"{q}\t{float(np.quantile(null, q)):.6f}\n")
        fh.write(f"null_mean\t{float(null.mean()):.6f}\n")
        fh.write(f"null_sd\t{float(null.std()):.6f}\n")
        fh.write(f"n_perm\t{args.n_perm}\n")
        fh.write(f"trait_entropy\t{H:.6f}\n")
        for kk, vv in gpd.items():
            fh.write(f"gpd_{kk}\t{vv}\n")

    counts = {c: int((finding == c).sum())
              for c in ("both", "mi_only", "pearson_only", "neither")}
    print(f"\n[out ] significant at padj <= {args.alpha}:")
    print(f"       MI      {int(sig_mi.sum()):,}")
    print(f"       Pearson {int(sig_r.sum()):,}")
    for c in ("both", "mi_only", "pearson_only"):
        print(f"       {c:<13}{counts[c]:,}")
    if counts["mi_only"]:
        top = np.argsort(np.where(finding == "mi_only", p_mi, np.inf))[:5]
        print("       strongest mi_only genes (mi, |r|, p source):")
        for i in top:
            print(f"         {genes[i]:<32} {mi[i]:.4f}  |r|={abs(r[i]):.3f}  "
                  f"{p_note[i]}")

    summary = dict(
        matrix=args.matrix, label=meta.get("label"), trait=args.trait,
        n_genes=n_genes, n_samples=n, k=args.k, n_perm=args.n_perm,
        alpha=args.alpha, trait_levels=[float(v) for v in levels],
        class_sizes=[int(c) for c in class_sizes], trait_entropy=H,
        tie_rate_sampled=tie_rate,
        null_mean=float(null.mean()), null_sd=float(null.std()),
        null_max=null_max, gpd_endpoint=None if not np.isfinite(endpoint) else float(endpoint),
        p_resolution=p_resolution,
        n_beyond_empirical_null=beyond_null, n_beyond_gpd_endpoint=beyond_fit,
        n_significant_mi=int(sig_mi.sum()), n_significant_pearson=int(sig_r.sum()),
        **{f"n_{c}": v for c, v in counts.items()},
        runtime_min=(time.time() - t_start) / 60.0,
        out_file=args.out + ".tsv")
    with open(args.out + ".summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"[done] {summary['runtime_min']:.1f} min -> {args.out}.*")


def bh(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values, in the input order."""
    m = p.size
    o = np.argsort(p)
    adj = np.minimum.accumulate((m * p[o] / np.arange(1, m + 1))[::-1])[::-1]
    out = np.empty_like(adj)
    out[o] = np.clip(adj, 0.0, 1.0)
    return out


if __name__ == "__main__":
    main()
