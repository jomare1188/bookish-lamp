#!/usr/bin/env python
# ============================================================================
# 03_merge_layers.py — union of the linear (Pearson) and non-linear (MI) layers
#
# WHAT THIS PRODUCES
#
# The network. One file per study, and the only edge table anything downstream
# reads:
#
#   gene1  gene2  stat  pval  padj  weight  pearson_r  ksg  source
#
#   stat       effective strength on the CORRELATION scale: |r| for edges the
#              Pearson layer found, r_eq for edges only MI found, the larger of
#              the two when both did. Column 3, where every consumer that reads
#              an edge statistic positionally expects it.
#   pval/padj  the more significant of the two tests
#   weight     0.01 + (stat-lo)/(hi-lo)*0.99 — minmax over [stat_min, stat_max],
#              the weighting the clustering step consumes
#   pearson_r  the true SIGNED r, or NA when only MI found the edge
#   ksg        mutual information in nats, or NA when only Pearson found it
#   source     pearson | mi | both
#
# `source` is the point of the whole exercise: it lets every downstream result
# be split by how the edge was found.
#
# WHY r_eq
#
# MI is in nats and Pearson in correlation units; the pipeline's weighting is
# defined on |r|. For a bivariate Gaussian, I = -0.5*ln(1-r^2) exactly, which
# inverts to
#
#       r_eq = sqrt(1 - exp(-2 * (I - bias)))
#
# after subtracting the estimator's small-sample null bias. KSG is biased
# UPWARD at these sample sizes, so skipping the correction would inflate every
# MI edge. The bias is read from the KSG layer's own <mi>.null.tsv (`null_mean`)
# rather than a table, so it always describes the run being merged.
#
# For a NON-MONOTONE edge, r_eq means "as strong as a linear r of this size",
# not "its Pearson r is this". That is why `source` must always be read
# alongside it.
#
# THE MERGE IS ONLY VALID IF BOTH SIDES ARE EQUALLY STRINGENT. Both layers must
# be cut at the same per-edge false-positive rate: |r| >= stat_min for the
# linear layer, and for the MI layer the MI value whose p-value equals the one
# |r| = stat_min implies at that n. This script re-derives that p and REFUSES
# to merge if the MI layer is looser — equal specificity is what licenses
# taking a union at all.
#
# METHOD: the MI side is much the smaller of the two, so it is held in RAM as
# sorted int64 gene-index pairs and the Pearson file is streamed once past it.
# Peak RAM is ~30 bytes per MI edge.
#
# RUN (normally through run.sh):
#   python 03_merge_layers.py --study sugarcane --n 48 \
#       --pearson results/sugarcane/layers/sugarcane_pearson.edgelist.tsv \
#       --mi      results/sugarcane/layers/sugarcane_ksg.edgelist.tsv \
#       --genes   results/sugarcane/vst/sugarcane.genes.txt \
#       --out     results/sugarcane/network_sugarcane_edges.tsv
# ============================================================================

import argparse
import json
import math
import os
import sys
import time

import numpy as np
import pandas as pd

# Fallback only. The real value comes from the KSG layer's null.tsv; these are
# the measured means from the 2026-08 sweeps and exist so a merge still runs if
# that file is missing.
NULL_BIAS_FALLBACK = {48: 0.124, 18: 0.180}


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def pvalue_at_r(r, n):
    """Two-sided p of Pearson r at n samples — the cut both layers must share."""
    from scipy import stats
    t = r * math.sqrt((n - 2) / (1 - r * r))
    return float(2 * stats.t.sf(abs(t), df=n - 2))


def r_equivalent(mi, bias, cap):
    """MI in nats -> equivalent |r|, via the Gaussian identity, bias-corrected."""
    d = np.maximum(mi - bias, 0.0)
    r = np.sqrt(1.0 - np.exp(-2.0 * d))
    return np.minimum(r, cap)


def read_null_bias(mi_path):
    """KSG null mean from the layer's own null.tsv, or None if unavailable."""
    null_f = mi_path.replace(".edgelist.tsv", ".null.tsv")
    if not os.path.exists(null_f):
        return None
    with open(null_f) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 2 and parts[0] == "null_mean":
                return float(parts[1])
    return None


def load_mi(path, gidx, n_genes, bias, stat_min, stat_max, expect_p, stat_col):
    """MI edgelist -> sorted keys + per-edge arrays. Aborts if not stringent."""
    log(f"reading MI layer {os.path.basename(path)}")
    df = pd.read_csv(path, sep="\t", dtype={
        "gene1": str, "gene2": str, stat_col: np.float32,
        "pval": np.float64, "padj": np.float64})
    log(f"  {len(df):,} MI edges")

    worst = float(df["pval"].max())
    if worst > expect_p * 1.001:
        sys.exit(
            f"ERROR: the MI layer contains p = {worst:.4g}, looser than the "
            f"{expect_p:.4g} that |r| = {stat_min} implies at n = {gidx.n_samples}.\n"
            f"       The two layers would not be equally stringent and the union "
            f"would be invalid. Re-run the MI layer with "
            f"--match-pearson {stat_min}.")

    i = gidx.get_indexer(df["gene1"])
    j = gidx.get_indexer(df["gene2"])
    if (i < 0).any() or (j < 0).any():
        sys.exit("ERROR: the MI layer references genes absent from genes.txt")
    i, j = i.astype(np.int64), j.astype(np.int64)
    keys = np.minimum(i, j) * n_genes + np.maximum(i, j)

    req = r_equivalent(df[stat_col].to_numpy(), bias, stat_max)
    keep = req >= stat_min
    n_drop = int((~keep).sum())
    if n_drop:
        log(f"  {n_drop:,} MI edges have r_eq < {stat_min} and are dropped "
            f"(they cannot enter a network defined on |r| >= {stat_min})")

    order = np.argsort(keys[keep], kind="stable")
    return dict(
        keys=keys[keep][order],
        ksg=df[stat_col].to_numpy()[keep][order],
        pval=df["pval"].to_numpy()[keep][order],
        padj=df["padj"].to_numpy()[keep][order],
        r_eq=req[keep][order].astype(np.float32),
    )


def main():
    ap = argparse.ArgumentParser(
        description="Merge the Pearson and MI network layers",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--study", required=True)
    ap.add_argument("--n", type=int, required=True, help="number of libraries")
    ap.add_argument("--pearson", required=True, help="pearson layer edgelist")
    ap.add_argument("--mi", required=True, help="ksg layer edgelist")
    ap.add_argument("--genes", required=True, help="<study>.genes.txt")
    ap.add_argument("--out", required=True)
    ap.add_argument("--stat-min", type=float, default=0.8)
    ap.add_argument("--stat-max", type=float, default=0.9999)
    ap.add_argument("--null-bias", type=float, default=None,
                    help="KSG null mean in nats; default reads it from the MI "
                         "layer's null.tsv")
    ap.add_argument("--chunk", type=int, default=4_000_000)
    ap.add_argument("--limit", type=int, default=0,
                    help="stop after N Pearson rows (testing only — the output "
                         "is then a partial network, not a usable one)")
    args = ap.parse_args()

    n = args.n
    bias = args.null_bias
    bias_src = "--null-bias"
    if bias is None:
        bias = read_null_bias(args.mi)
        bias_src = f"{os.path.basename(args.mi).replace('.edgelist.tsv', '.null.tsv')} null_mean"
    if bias is None:
        bias = NULL_BIAS_FALLBACK.get(n)
        bias_src = "built-in fallback table"
    if bias is None:
        sys.exit(f"ERROR: no KSG null bias available for n = {n}. Pass "
                 f"--null-bias explicitly.")

    expect_p = pvalue_at_r(args.stat_min, n)

    log(f"study {args.study}  n = {n}")
    log(f"  KSG null bias {bias:+.4f} nats   [{bias_src}]")
    log(f"  matching cut: |r| >= {args.stat_min}  ==  p <= {expect_p:.4g}")
    if args.limit:
        log(f"  *** --limit {args.limit:,}: PARTIAL OUTPUT, for testing only")

    # Gene ids were normalised once, in 01_export_vst.r, and both layers were
    # built from this exact file — so no stripping or reconciliation here. The
    # duplicate check is what would catch it if that ever stopped being true.
    genes = [ln.strip() for ln in open(args.genes) if ln.strip()]
    n_genes = len(genes)
    if len(set(genes)) != n_genes:
        sys.exit(f"ERROR: {args.genes} has duplicate gene ids "
                 f"({n_genes:,} lines, {len(set(genes)):,} unique)")
    log(f"  {n_genes:,} genes")
    gidx = pd.Index(genes)
    gidx.n_samples = n                      # for the error message in load_mi

    mi = load_mi(args.mi, gidx, n_genes, bias, args.stat_min, args.stat_max,
                 expect_p, stat_col="ksg")
    n_mi = len(mi["keys"])
    seen = np.zeros(n_mi, dtype=bool)

    lo_w, hi_w = args.stat_min, args.stat_max
    span = hi_w - lo_w

    fmt = "%s\t%s\t%.5f\t%.4g\t%.4g\t%.6f\t%s\t%s\t%s\n".__mod__
    gname = np.asarray(genes, dtype=object)

    n_both = n_pear = 0
    out_of_range = 0
    t0 = time.time()

    log(f"streaming {os.path.basename(args.pearson)}")
    with open(args.out, "w") as out:
        out.write("gene1\tgene2\tstat\tpval\tpadj\tweight\tpearson_r\tksg\t"
                  "source\n")

        reader = pd.read_csv(
            args.pearson, sep="\t", chunksize=args.chunk,
            usecols=[0, 1, 2, 3, 4],
            names=["gene1", "gene2", "pearson", "pval", "padj"],
            header=0,
            dtype={"gene1": str, "gene2": str, "pearson": np.float64,
                   "pval": np.float64, "padj": np.float64})

        done = 0
        for chunk in reader:
            i = gidx.get_indexer(chunk["gene1"])
            j = gidx.get_indexer(chunk["gene2"])
            if (i < 0).any() or (j < 0).any():
                bad = chunk["gene1"][i < 0].head(3).tolist() + \
                      chunk["gene2"][j < 0].head(3).tolist()
                sys.exit("ERROR: the Pearson layer references genes absent from "
                         f"genes.txt — the two layers are not the same gene "
                         f"set. Examples: {bad}")
            i, j = i.astype(np.int64), j.astype(np.int64)
            keys = np.minimum(i, j) * n_genes + np.maximum(i, j)

            r = chunk["pearson"].to_numpy()
            absr = np.abs(r)
            stat = absr.copy()
            pval = chunk["pval"].to_numpy().copy()   # pandas hands back a
            padj = chunk["padj"].to_numpy().copy()   # read-only view

            # The engine already applied [stat_min, stat_max]; this catches a
            # layer built with different thresholds being merged by mistake.
            out_of_range += int(((absr < lo_w) | (absr > hi_w)).sum())

            pos = np.searchsorted(mi["keys"], keys)
            np.clip(pos, 0, max(n_mi - 1, 0), out=pos)
            hit = (n_mi > 0) & (mi["keys"][pos] == keys)
            ksg_s = np.full(len(chunk), "NA", dtype=object)
            if hit.any():
                h = np.nonzero(hit)[0]
                p = pos[h]
                seen[p] = True
                stat[h] = np.maximum(absr[h], mi["r_eq"][p])
                pval[h] = np.minimum(pval[h], mi["pval"][p])
                padj[h] = np.minimum(padj[h], mi["padj"][p])
                ksg_s[h] = np.char.mod("%.5f", mi["ksg"][p]).astype(object)
                n_both += len(h)

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

    if out_of_range:
        log(f"WARNING: {out_of_range:,} Pearson edges fell outside "
            f"[{lo_w}, {hi_w}] — the layer was not built with these thresholds")

    log("")
    log(f"{'pearson only':<16}{n_pear:>15,}  {n_pear/total:>7.2%}")
    log(f"{'both':<16}{n_both:>15,}  {n_both/total:>7.2%}")
    log(f"{'MI only':<16}{n_mionly:>15,}  {n_mionly/total:>7.2%}")
    log(f"{'TOTAL':<16}{total:>15,}")
    log(f"MI edges also found by Pearson: "
        f"{n_both/max(n_both+n_mionly,1):.2%}  "
        f"(the rest are what the linear layer could not see)")
    log(f"written to {args.out}  in {dt/60:.1f} min")

    summary = dict(
        study=args.study, n_samples=n, null_bias=bias, null_bias_source=bias_src,
        stat_min=args.stat_min, stat_max=args.stat_max,
        matched_p=expect_p,
        n_pearson_only=n_pear, n_both=n_both, n_mi_only=n_mionly,
        n_total=total, partial=bool(args.limit), runtime_min=dt / 60,
        pearson_file=args.pearson, mi_file=args.mi, out_file=args.out)
    with open(args.out.replace(".tsv", ".summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)


if __name__ == "__main__":
    main()
