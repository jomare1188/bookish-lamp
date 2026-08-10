#!/usr/bin/env python
# ============================================================================
# 98_power_curve.py — what dependence strength can MI actually detect at this n?
#
# WHY THIS EXISTS
#
# Running the sweep with --match-pearson 0.7 sets the MI floor at the value with
# the same NULL TAIL PROBABILITY as |r| = 0.7. That equalises the false-positive
# rate, not the effect size, and the two are very far apart at small n: the KSG
# estimator is much noisier than Pearson, so buying the same specificity costs a
# great deal of power. On the first 20,000 sugarcane genes, 90% of the Pearson
# |r| >= 0.7 pairs did NOT clear the matched MI floor.
#
# This script makes that trade-off explicit. For bivariate-Gaussian pairs of
# known correlation it reports, at the study's own n:
#
#   * the true MI,        -0.5*ln(1-r^2)
#   * what KSG estimates (mean +/- sd) -- including its small-sample upward bias
#   * POWER: the fraction of such pairs that clear the matched floor
#   * the equivalent |r| the floor really corresponds to
#
# Read it before deciding whether to replace Pearson with MI: if power at
# |r| = 0.7 is near zero, an MI network is not a "better Pearson network", it is
# a much smaller network of much stronger relationships that happens to also
# include non-monotone ones.
#
# RUN: /home/genomics/miniconda3/envs/docling/bin/python 98_power_curve.py
# ============================================================================

import argparse
import contextlib
import importlib.util
import io
import math
import os

import numpy as np
import torch
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "cmi", os.path.join(HERE, "custom_mutual_information.py"))
cmi = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cmi)

DEV = torch.device("cuda", 0) if torch.cuda.is_available() else torch.device("cpu")


def ksg_gaussian(n, r, n_pairs, k=3, seed=0):
    """KSG estimates for n_pairs independent bivariate-Gaussian samples."""
    rng = np.random.default_rng(seed)
    z1 = rng.standard_normal((n_pairs, n))
    z2 = r * z1 + math.sqrt(max(1 - r * r, 0.0)) * rng.standard_normal((n_pairs, n))
    est = cmi.ESTIMATORS["ksg"](n=n, k=k, device=DEV, half=True)
    out = []
    for s in range(0, n_pairs, 8192):
        a = np.argsort(np.argsort(z1[s:s + 8192], axis=1), axis=1).astype(np.int32)
        b = np.argsort(np.argsort(z2[s:s + 8192], axis=1), axis=1).astype(np.int32)
        pa = est.prep(torch.from_numpy(a).to(DEV))
        pb = est.prep(torch.from_numpy(b).to(DEV))
        out.append(est.block_paired(pa, pb).float().cpu().numpy())
    return np.concatenate(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, nargs="+", default=[48, 18],
                    help="library counts to profile")
    ap.add_argument("--k", type=int, default=3)
    ap.add_argument("--match-pearson", type=float, default=0.8)
    ap.add_argument("--n-perm", type=int, default=20_000_000)
    ap.add_argument("--n-pairs", type=int, default=200_000)
    args = ap.parse_args()

    for n in args.n:
        est = cmi.ESTIMATORS["ksg"](n=n, k=args.k, device=DEV, half=True)
        with contextlib.redirect_stdout(io.StringIO()):
            null = cmi.build_null(est, n, args.n_perm, DEV, 1, chunk_pairs=131072)
        gpd = cmi.fit_gpd_tail(null)
        pf = cmi.make_pvalue_fn(null, gpd)

        r0 = args.match_pearson
        t = r0 * math.sqrt((n - 2) / (1 - r0 * r0))
        p_r = float(2 * stats.t.sf(abs(t), df=n - 2))
        lo, hi = float(null.min()) - 1, float(null.max()) + 10
        for _ in range(200):
            mid = (lo + hi) / 2
            if pf(np.array([mid]))[0] > p_r:
                lo = mid
            else:
                hi = mid
        floor = (lo + hi) / 2

        print(f"\n{'='*78}\nn = {n} libraries   k = {args.k}")
        print(f"  null: mean {null.mean():+.4f}  sd {null.std():.4f}   "
              f"(KSG's small-sample bias — the null is NOT centred on 0)")
        print(f"  |r| = {r0} at this n is p = {p_r:.3g}  ->  matched MI floor = {floor:.4f} nats")
        print(f"\n  {'|r|':>5} {'true MI':>9} {'KSG mean':>10} {'KSG sd':>8} "
              f"{'power':>8}   (power = P(KSG >= floor))")
        pow_at = {}
        for r in (0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 0.99):
            v = ksg_gaussian(n, r, args.n_pairs, args.k, seed=int(r * 1000))
            truth = -0.5 * math.log(1 - r * r)
            power = float(np.mean(v >= floor))
            pow_at[r] = power
            print(f"  {r:>5.2f} {truth:>9.4f} {v.mean():>10.4f} {v.std():>8.4f} "
                  f"{power:>8.1%}")

        # what |r| does the floor really correspond to, in power terms?
        eq = [r for r, p in sorted(pow_at.items()) if p >= 0.5]
        eq = eq[0] if eq else None
        print(f"\n  => the matched floor behaves like a Pearson cut at "
              f"|r| ~ {eq if eq else '>0.99'} (the weakest |r| detected "
              f"more than half the time).")
        if pow_at.get(r0, 0) < 0.2:
            print(f"  => power at the nominal |r| = {r0} is only {pow_at[r0]:.1%}: at "
                  f"n = {n}, MI cannot resolve a relationship of that strength from "
                  f"noise once {1.5e10:.0e} pairs are corrected for.")


if __name__ == "__main__":
    main()
