#!/usr/bin/env python
# ============================================================================
# 10_augment_network.py — union of the Pearson and mutual-information networks
#
# WHAT THIS PRODUCES
#
# One file per study, in the SAME format general_stats.r writes
# (network_*_filtered_edges.tsv), so every downstream script consumes it
# unchanged: mcl_clustering.r, eigengene.r, module_trait_cor.r,
# comparative_networks2.r, full_matrix_edge_comparation.r and
# network_conservation_join.r all read only gene1/gene2/weight, or just the
# gene pair.
#
#   gene1  gene2  stat  pval  padj  weight  pearson_r  ksg  source
#
#   stat       effective strength on the CORRELATION scale: |r| for edges
#              Pearson found, r_eq for edges only MI found, max(...) for both.
#              general_stats.r reads column 3 positionally as "pearson", so if
#              that script is re-run on this file the column lands correctly.
#   pval/padj  the more significant of the two tests
#   weight     0.01 + (stat-lo)/(hi-lo)*0.99, identical to general_stats.r's
#              minmax with lo/hi fixed at PEARSON_MIN/PEARSON_MAX. The script
#              VERIFIES this reproduces the input file's own weight column.
#   pearson_r  the true signed r, or NA when only MI found the edge
#   ksg        mutual information in nats, or NA when only Pearson found it
#   source     pearson | mi | both
#
# WHY r_eq
#
# MI is in nats and Pearson in correlation units; the pipeline's weighting is
# defined on |r|. For a bivariate Gaussian, I = -0.5*ln(1-r^2) exactly, which
# inverts to
#
#       r_eq = sqrt(1 - exp(-2 * (I - bias)))
#
# after subtracting the estimator's small-sample null bias (KSG is biased
# UPWARD: +0.124 nats at n=48, +0.180 at n=18 — see MI_for_review.md 2.5).
# For a non-monotone edge r_eq means "as strong as a linear r of this size",
# NOT "its Pearson r is this". That is exactly why `source` exists: never
# interpret r_eq without it.
#
# THE MERGE IS ONLY VALID IF BOTH SIDES ARE EQUALLY STRINGENT.  Both must be
# cut at the same per-edge false-positive rate, i.e. Pearson |r| >= 0.8 against
# an MI edgelist filtered to the p-value that |r| = 0.8 implies at that n
# (9.02e-12 at n=48, 6.72e-05 at n=18). Use the *.p08.edgelist.tsv files.
#
# METHOD: the MI side is much the smaller of the two, so it is held in RAM as
# sorted int64 keys (gene index pair) and the Pearson file is streamed once.
# Peak RAM is ~30 bytes per MI edge (~2.4 GB for purple's 79 M).
#
# RUN:
#   python 10_augment_network.py --study sugarcane        # reads config.sh
#   python 10_augment_network.py --study purple
#   python 10_augment_network.py --study sugarcane --limit 2000000   # test
# ============================================================================

import argparse
import json
import math
import os
import re
import sys
import time

import numpy as np
import pandas as pd

# The Pearson networks carry bare locus ids; the MI gene lists (taken straight
# from the dds rownames) carry the R570 annotation version suffix, e.g.
#   genes.txt                  SoffiXsponR570.01Ag000100.v2.1
#   network_*_filtered_edges   SoffiXsponR570.01Ag000100
# Stripping it is checked for collisions before use. Purple ids
# (Soffic.01A0000290-1A) are unaffected by this pattern.
VERSION_SUFFIX = r"\.v[0-9]+(?:\.[0-9]+)*$"

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.abspath(os.path.join(HERE, "..", ".."))

# Null bias of the KSG estimator by sample size (MI_for_review.md 2.5).
NULL_BIAS = {48: 0.124, 18: 0.180}

# Per-edge p-value that |r| = PEARSON_MIN implies at each n. The MI edgelist
# must already be filtered to this; we re-check and refuse if it is not.
STUDIES = {
    "sugarcane": dict(
        n=48,
        pearson=f"{BASE}/files/sugarcane/network_sugarcane_filtered_edges.tsv",
        mi=f"{BASE}/files/sugarcane/mi/sugarcane_ksg.p08.edgelist.tsv",
        genes=f"{BASE}/files/sugarcane/mi/sugarcane.genes.txt",
        out=f"{BASE}/files/sugarcane/network_sugarcane_augmented_edges.tsv",
    ),
    "purple": dict(
        n=18,
        pearson=f"{BASE}/files/purple/new/network_purple_filtered_edges.tsv",
        mi=f"{BASE}/files/purple/new/mi/purple_ksg.p08.edgelist.tsv",
        genes=f"{BASE}/files/purple/new/mi/purple.genes.txt",
        out=f"{BASE}/files/purple/new/network_purple_augmented_edges.tsv",
    ),
}


def mi_pvalue_at(r, n):
    """Two-sided p of Pearson r at n samples — the cut the MI side must match."""
    from scipy import stats
    t = r * math.sqrt((n - 2) / (1 - r * r))
    return float(2 * stats.t.sf(abs(t), df=n - 2))


def r_equivalent(mi, bias, cap):
    """MI in nats -> equivalent |r|, via the Gaussian identity, bias-corrected."""
    d = np.maximum(mi - bias, 0.0)
    r = np.sqrt(1.0 - np.exp(-2.0 * d))
    return np.minimum(r, cap)


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ---------------------------------------------------------------------------


def load_mi(path, gidx, n_genes, bias, stat_min, stat_max, expect_p):
    """MI edgelist -> (sorted keys, ksg, pval, padj, r_eq). Aborts on drift."""
    log(f"reading MI edgelist {os.path.basename(path)}")
    df = pd.read_csv(path, sep="\t", dtype={
        "gene1": str, "gene2": str, "ksg": np.float32,
        "pval": np.float64, "padj": np.float64})
    log(f"  {len(df):,} MI edges")

    worst = df["pval"].max()
    if worst > expect_p * 1.001:
        sys.exit(
            f"ERROR: MI edgelist contains p = {worst:.4g}, looser than the "
            f"{expect_p:.4g} that |r| = {STAT_MIN_LABEL} implies at this n.\n"
            f"       The two sides would not be equally stringent and the "
            f"union would be invalid.\n"
            f"       Filter it first:  awk -F'\\t' 'NR==1 || $4+0 <= {expect_p:.4g}' "
            f"in.tsv > out.tsv")

    i = gidx.get_indexer(df["gene1"].str.replace(VERSION_SUFFIX, "", regex=True))
    j = gidx.get_indexer(df["gene2"].str.replace(VERSION_SUFFIX, "", regex=True))
    if (i < 0).any() or (j < 0).any():
        sys.exit("ERROR: MI edgelist references genes absent from genes.txt")
    i = i.astype(np.int64)
    j = j.astype(np.int64)

    lo = np.minimum(i, j)
    hi = np.maximum(i, j)
    keys = lo * n_genes + hi

    req = r_equivalent(df["ksg"].to_numpy(), bias, stat_max)
    keep = req >= stat_min
    n_drop = int((~keep).sum())
    if n_drop:
        log(f"  {n_drop:,} MI edges have r_eq < {stat_min} and are dropped "
            f"(they cannot enter a network defined on |r| >= {stat_min})")

    order = np.argsort(keys[keep], kind="stable")
    return dict(
        keys=keys[keep][order],
        ksg=df["ksg"].to_numpy()[keep][order],
        pval=df["pval"].to_numpy()[keep][order],
        padj=df["padj"].to_numpy()[keep][order],
        r_eq=req[keep][order].astype(np.float32),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--study", choices=sorted(STUDIES), required=True)
    ap.add_argument("--pearson"), ap.add_argument("--mi")
    ap.add_argument("--genes"), ap.add_argument("--out")
    ap.add_argument("--stat-min", type=float, default=0.8,
                    help="general_stats.r PEARSON_MIN")
    ap.add_argument("--stat-max", type=float, default=0.9999,
                    help="general_stats.r PEARSON_MAX")
    ap.add_argument("--null-bias", type=float, default=None,
                    help="KSG null mean at this n (default: table lookup)")
    ap.add_argument("--chunk", type=int, default=4_000_000)
    ap.add_argument("--limit", type=int, default=0,
                    help="stop after N Pearson rows (testing only — the output "
                         "is then a partial network, not a usable one)")
    args = ap.parse_args()

    cfg = STUDIES[args.study]
    pearson_f = args.pearson or cfg["pearson"]
    mi_f = args.mi or cfg["mi"]
    genes_f = args.genes or cfg["genes"]
    out_f = args.out or cfg["out"]
    n = cfg["n"]
    bias = args.null_bias if args.null_bias is not None else NULL_BIAS[n]

    global STAT_MIN_LABEL
    STAT_MIN_LABEL = args.stat_min
    expect_p = mi_pvalue_at(args.stat_min, n)

    log(f"study {args.study}  n = {n}  KSG null bias {bias:+.3f} nats")
    log(f"  matching cut: |r| >= {args.stat_min}  ==  p <= {expect_p:.4g}")
    if args.limit:
        log(f"  *** --limit {args.limit:,}: PARTIAL OUTPUT, for testing only")

    raw = [ln.strip() for ln in open(genes_f) if ln.strip()]
    genes = [re.sub(VERSION_SUFFIX, "", g) for g in raw]
    n_genes = len(genes)
    n_stripped = sum(a != b for a, b in zip(raw, genes))
    if len(set(genes)) != n_genes:
        sys.exit(f"ERROR: stripping the version suffix collapses {n_genes:,} "
                 f"gene ids into {len(set(genes)):,} unique names. Rerun with "
                 f"--keep-suffix and fix the naming instead.")
    log(f"  {n_genes:,} genes"
        + (f"  ({n_stripped:,} had a version suffix stripped to match the "
           f"Pearson network's ids)" if n_stripped else ""))
    gidx = pd.Index(genes)

    mi = load_mi(mi_f, gidx, n_genes, bias, args.stat_min, args.stat_max,
                 expect_p)
    n_mi = len(mi["keys"])
    seen = np.zeros(n_mi, dtype=bool)

    lo_w, hi_w = args.stat_min, args.stat_max
    span = hi_w - lo_w

    fmt = "%s\t%s\t%.5f\t%.4g\t%.4g\t%.6f\t%s\t%s\t%s\n".__mod__
    gname = np.asarray(genes, dtype=object)

    n_both = n_pear = 0
    weight_err = 0.0
    checked = 0
    t0 = time.time()

    log(f"streaming {os.path.basename(pearson_f)}")
    with open(out_f, "w") as out:
        out.write("gene1\tgene2\tstat\tpval\tpadj\tweight\tpearson_r\tksg\t"
                  "source\n")

        reader = pd.read_csv(
            pearson_f, sep="\t", chunksize=args.chunk,
            usecols=[0, 1, 2, 3, 4, 5],
            names=["gene1", "gene2", "pearson", "pval", "padj", "weight"],
            header=0,
            dtype={"gene1": str, "gene2": str, "pearson": np.float64,
                   "pval": np.float64, "padj": np.float64,
                   "weight": np.float64})

        done = 0
        for chunk in reader:
            i = gidx.get_indexer(chunk["gene1"])
            j = gidx.get_indexer(chunk["gene2"])
            if (i < 0).any() or (j < 0).any():
                bad = chunk["gene1"][i < 0].head(3).tolist() + \
                      chunk["gene2"][j < 0].head(3).tolist()
                sys.exit("ERROR: Pearson edges reference genes absent from "
                         f"genes.txt — the two sides are not the same gene "
                         f"set. Examples: {bad}")
            i = i.astype(np.int64)
            j = j.astype(np.int64)
            keys = np.minimum(i, j) * n_genes + np.maximum(i, j)

            r = chunk["pearson"].to_numpy()
            absr = np.abs(r)
            stat = absr.copy()
            pval = chunk["pval"].to_numpy().copy()   # pandas hands back a
            padj = chunk["padj"].to_numpy().copy()   # read-only view

            # self-check: does our weight formula reproduce general_stats.r's?
            if checked < 5_000_000:
                w_ours = 0.01 + (absr - lo_w) / span * 0.99
                weight_err = max(weight_err,
                                 float(np.abs(w_ours - chunk["weight"].to_numpy()).max()))
                checked += len(chunk)

            pos = np.searchsorted(mi["keys"], keys)
            np.clip(pos, 0, max(n_mi - 1, 0), out=pos)
            hit = (n_mi > 0) & (mi["keys"][pos] == keys)
            if hit.any():
                h = np.nonzero(hit)[0]
                p = pos[h]
                seen[p] = True
                stat[h] = np.maximum(absr[h], mi["r_eq"][p])
                pval[h] = np.minimum(pval[h], mi["pval"][p])
                padj[h] = np.minimum(padj[h], mi["padj"][p])
                n_both += len(h)

            ksg_s = np.full(len(chunk), "NA", dtype=object)
            if hit.any():
                ksg_s[h] = np.char.mod("%.5f", mi["ksg"][p]).astype(object)
            src = np.where(hit, "both", "pearson")
            n_pear += len(chunk) - int(hit.sum())

            weight = 0.01 + (stat - lo_w) / span * 0.99
            rs = np.char.mod("%.5f", r)
            out.write("".join(map(fmt, zip(
                gname[i], gname[j], stat, pval, padj, weight, rs, ksg_s, src))))

            done += len(chunk)
            if done % (args.chunk * 5) == 0:
                log(f"  {done:,} Pearson rows  ({done/(time.time()-t0):,.0f}/s)")
            if args.limit and done >= args.limit:
                log(f"  stopping at --limit {args.limit:,}")
                break

        # ---- MI-only edges ------------------------------------------------
        only = np.nonzero(~seen)[0]
        log(f"appending {len(only):,} MI-only edges")
        for s in range(0, len(only), args.chunk):
            sl = only[s:s + args.chunk]
            k = mi["keys"][sl]
            i = (k // n_genes).astype(np.int64)
            j = (k % n_genes).astype(np.int64)
            stat = mi["r_eq"][sl].astype(np.float64)
            weight = 0.01 + (stat - lo_w) / span * 0.99
            ksg_s = np.char.mod("%.5f", mi["ksg"][sl])
            out.write("".join(map(fmt, zip(
                gname[i], gname[j], stat, mi["pval"][sl], mi["padj"][sl],
                weight, ["NA"] * len(sl), ksg_s, ["mi"] * len(sl)))))
        n_mionly = len(only)

    dt = time.time() - t0
    total = n_pear + n_both + n_mionly

    log("")
    log(f"weight formula self-check: max |ours - general_stats.r| = "
        f"{weight_err:.2e} over {checked:,} rows"
        + ("   OK" if weight_err < 1e-6 else "   *** MISMATCH ***"))
    if weight_err >= 1e-6:
        log("    the input's weights were not built with minmax over "
            f"[{lo_w}, {hi_w}] — check general_stats.r's WEIGHT_METHOD")

    log("")
    log(f"{'pearson only':<16}{n_pear:>15,}  {n_pear/total:>7.2%}")
    log(f"{'both':<16}{n_both:>15,}  {n_both/total:>7.2%}")
    log(f"{'MI only':<16}{n_mionly:>15,}  {n_mionly/total:>7.2%}")
    log(f"{'TOTAL':<16}{total:>15,}")
    log(f"MI edges recovered by Pearson: "
        f"{n_both/max(n_both+n_mionly,1):.2%}  "
        f"(the rest are what Pearson could not see)")
    log(f"written to {out_f}  in {dt/60:.1f} min")

    summary = dict(
        study=args.study, n_samples=n, null_bias=bias,
        stat_min=args.stat_min, stat_max=args.stat_max,
        matched_p=expect_p, weight_max_abs_err=weight_err,
        n_pearson_only=n_pear, n_both=n_both, n_mi_only=n_mionly,
        n_total=total, partial=bool(args.limit), runtime_min=dt / 60,
        pearson_file=pearson_f, mi_file=mi_f, out_file=out_f)
    with open(out_f.replace(".tsv", ".summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)


if __name__ == "__main__":
    main()
