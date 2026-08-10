#!/usr/bin/env python
# ============================================================================
# 99_validate.py — correctness tests for custom_mutual_information.py
#
# The GPU KSG kernel is a fused, blocked, fp16 reimplementation of an estimator
# that is normally written as a loop over pairs. That is exactly the kind of
# code that can be fast and quietly wrong, so nothing from it should be trusted
# until these five checks pass:
#
#   1 REFERENCE   GPU KSG == a deliberately naive CPU KSG, pair by pair.
#   2 ANALYTIC    On a bivariate Gaussian with known r and large n, KSG
#                 recovers the true MI = -0.5*ln(1-r^2).
#   3 NON-LINEAR  On y = x^2, y = |x|, y = sin(3x), where Pearson collapses,
#                 KSG and xi still see the dependence. This is the entire
#                 reason for the script's existence.
#   4 NULL        Under independence the reported p-values are uniform, and the
#                 GPD tail extrapolation tracks the empirical tail. This is what
#                 makes the shared-null shortcut legitimate.
#   5 PRECISION   fp16 and fp32 distances give bit-identical MI (rank
#                 differences are small integers, so fp16 is exact).
#
# RUN: /home/genomics/miniconda3/envs/docling/bin/python 99_validate.py
# ============================================================================

import importlib.util
import os
import sys

import numpy as np
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "cmi", os.path.join(HERE, "custom_mutual_information.py"))
cmi = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cmi)

DEV = torch.device("cuda", 0) if torch.cuda.is_available() else torch.device("cpu")
FAILURES = []


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  {detail}" if detail else ""))
    if not ok:
        FAILURES.append(name)


def gpu_mi(X, k=3, half=True, estimator="ksg"):
    """Run the production kernel on a small matrix, return the full G x G result."""
    ranks, _ = cmi.rank_transform(np.asarray(X, dtype=np.float32), seed=1)
    est = cmi.ESTIMATORS[estimator](n=X.shape[1], k=k, device=DEV, half=half)
    p = est.prep(torch.from_numpy(ranks).to(DEV))
    return est.block(p, p).float().cpu().numpy()


# ---------------------------------------------------------------------------
# 1. reference implementation — written for obviousness, not speed
# ---------------------------------------------------------------------------
def cpu_ksg_pair(x, y, k=3):
    """Kraskov algorithm 1, one pair, straight from the paper."""
    from scipy.special import digamma
    n = len(x)
    dx = np.abs(x[:, None] - x[None, :])
    dy = np.abs(y[:, None] - y[None, :])
    d = np.maximum(dx, dy)
    mi = 0.0
    for i in range(n):
        eps = np.sort(d[i])[k]           # [0] is the self-distance
        nx = np.sum(dx[i] < eps) - 1
        ny = np.sum(dy[i] < eps) - 1
        mi += digamma(nx + 1) + digamma(ny + 1)
    return digamma(k) + digamma(n) - mi / n


def test_reference():
    print("\n[1] GPU kernel vs naive CPU KSG")
    rng = np.random.default_rng(11)
    n, G = 48, 40
    X = rng.standard_normal((G, n)).astype(np.float32)
    X[1] = X[0] + 0.3 * rng.standard_normal(n)          # dependent
    X[2] = X[0] ** 2 + 0.3 * rng.standard_normal(n)     # non-monotone

    ranks, _ = cmi.rank_transform(X, seed=1)
    M = gpu_mi(X, k=3)
    ref = np.array([[cpu_ksg_pair(ranks[i].astype(float), ranks[j].astype(float), 3)
                     for j in range(G)] for i in range(G)])
    off = ~np.eye(G, dtype=bool)
    err = np.abs(M[off] - ref[off]).max()
    check("max |GPU - CPU| over 1560 pairs", err < 1e-4, f"err = {err:.2e}")


# ---------------------------------------------------------------------------
# 2. analytic ground truth
# ---------------------------------------------------------------------------
def test_analytic():
    print("\n[2] bivariate Gaussian, known MI = -0.5*ln(1-r^2)")
    rng = np.random.default_rng(23)
    n = 4000
    for r in (0.3, 0.6, 0.9):
        z1 = rng.standard_normal(n)
        z2 = r * z1 + np.sqrt(1 - r * r) * rng.standard_normal(n)
        truth = -0.5 * np.log(1 - r * r)
        X = np.stack([z1, z2]).astype(np.float32)
        est = gpu_mi(X, k=5)[0, 1]
        gc = gpu_mi(X, estimator="gcmi")[0, 1]
        check(f"r = {r}: true {truth:.4f}",
              abs(est - truth) < 0.06,
              f"KSG {est:.4f} (err {est-truth:+.4f}) | GCMI {gc:.4f}")


# ---------------------------------------------------------------------------
# 3. the point of the exercise: non-linear dependence
# ---------------------------------------------------------------------------
def test_nonlinear():
    print("\n[3] non-linear dependence that Pearson cannot see  (n = 48)")
    rng = np.random.default_rng(37)
    n = 48
    x = rng.standard_normal(n)
    cases = {
        "y = x + e        (linear)":      x + 0.25 * rng.standard_normal(n),
        "y = x^2 + e      (parabola)":    x ** 2 + 0.25 * rng.standard_normal(n),
        "y = |x| + e      (V-shape)":     np.abs(x) + 0.25 * rng.standard_normal(n),
        "y = sin(3x) + e  (oscillating)": np.sin(3 * x) + 0.25 * rng.standard_normal(n),
        "y = e            (independent)": rng.standard_normal(n),
    }
    X = np.stack([x] + list(cases.values())).astype(np.float32)
    ksg = gpu_mi(X, k=3)[0, 1:]
    xi = gpu_mi(X, estimator="xi")[0, 1:]
    gc = gpu_mi(X, estimator="gcmi")[0, 1:]
    pear = np.array([abs(np.corrcoef(x, y)[0, 1]) for y in cases.values()])

    print(f"    {'relationship':32s} {'|Pearson|':>9s} {'KSG':>8s} {'GCMI':>8s} {'xi':>8s}")
    for lab, p, m, g, c in zip(cases, pear, ksg, gc, xi):
        print(f"    {lab:32s} {p:9.3f} {m:8.3f} {g:8.3f} {c:8.3f}")

    # the parabola/V/sine cases must beat the independent case on KSG and xi
    ind_k, ind_x = ksg[-1], xi[-1]
    check("KSG separates all 4 dependent cases from noise",
          bool(np.all(ksg[:4] > ind_k + 0.1)), f"noise KSG = {ind_k:.3f}")
    check("xi separates all 4 dependent cases from noise",
          bool(np.all(xi[:4] > ind_x + 0.1)), f"noise xi = {ind_x:.3f}")
    check("Pearson misses the oscillating case but KSG does not",
          pear[3] < 0.3 < ksg[3], f"|r| = {pear[3]:.3f} vs KSG = {ksg[3]:.3f}")


# ---------------------------------------------------------------------------
# 3b. the paired kernel used to build the null must equal the block kernel
# ---------------------------------------------------------------------------
def test_rectangular():
    """Non-square tiles happen at the end of every sweep, so they must give the
    same answers as the square block that contains them."""
    print("\n[3c] rectangular tiles == the corresponding sub-block")
    rng = np.random.default_rng(61)
    n, G = 48, 96
    X = rng.standard_normal((G, n)).astype(np.float32)
    X[1] = X[0] ** 2 + 0.3 * rng.standard_normal(n)
    ranks, _ = cmi.rank_transform(X, seed=1)
    for name in ("ksg", "gcmi", "xi"):
        est = cmi.ESTIMATORS[name](n=n, k=3, device=DEV, half=True)
        R = torch.from_numpy(ranks).to(DEV)
        full = est.block(est.prep(R), est.prep(R)).float().cpu().numpy()
        # a deliberately lopsided tile, like the last one of a real run
        sub = est.block(est.prep(R[:70]), est.prep(R[64:])).float().cpu().numpy()
        err = float(np.abs(sub - full[:70, 64:]).max())
        check(f"{name}: 70x32 tile matches the full matrix", err < 1e-5,
              f"max diff = {err:.2e}")


def test_paired():
    print("\n[3b] block_paired() == diagonal of block(), all three estimators")
    rng = np.random.default_rng(53)
    n, m = 48, 64
    a = np.argsort(rng.random((m, n)), axis=1).astype(np.int32)
    b = np.argsort(rng.random((m, n)), axis=1).astype(np.int32)
    for name in ("ksg", "gcmi", "xi"):
        est = cmi.ESTIMATORS[name](n=n, k=3, device=DEV, half=True)
        ta = est.prep(torch.from_numpy(a).to(DEV))
        tb = est.prep(torch.from_numpy(b).to(DEV))
        full = torch.diagonal(est.block(ta, tb)).float().cpu().numpy()
        pair = est.block_paired(ta, tb).float().cpu().numpy()
        err = float(np.abs(full - pair).max())
        check(f"{name}: paired vs block diagonal", err < 1e-5, f"max diff = {err:.2e}")


# ---------------------------------------------------------------------------
# 4. null calibration — the shared-null shortcut
# ---------------------------------------------------------------------------
def test_null():
    print("\n[4] null calibration and GPD tail")
    n, k = 48, 3
    est = cmi.ESTIMATORS["ksg"](n=n, k=k, device=DEV, half=True)
    null = cmi.build_null(est, n, 200_000, DEV, seed=5, chunk_pairs=16384)
    gpd = cmi.fit_gpd_tail(null)
    pf = cmi.make_pvalue_fn(null, gpd)

    # independent draws from a DIFFERENT seed must yield uniform p-values
    fresh = cmi.build_null(est, n, 50_000, DEV, seed=99, chunk_pairs=16384)
    p = pf(fresh)
    for a in (0.05, 0.01, 0.001):
        obs = float(np.mean(p <= a))
        check(f"P(p <= {a}) under H0 == {a}", abs(obs - a) < max(0.25 * a, 0.002),
              f"observed {obs:.5f}")

    # the GPD extrapolation must agree with the empirical tail where both exist
    q = float(np.quantile(null, 0.9999))
    emp = float(np.mean(fresh >= q))
    check("GPD tail matches empirical at the 1e-4 level",
          abs(np.log10(max(pf(np.array([q]))[0], 1e-12)) - np.log10(max(emp, 1e-12))) < 0.35,
          f"GPD {pf(np.array([q]))[0]:.2e} vs empirical {emp:.2e}")
    print(f"    tail fit: {gpd['method']}  u = {gpd['u']:.4f}  "
          f"xi = {gpd['xi']:+.4f}  sigma = {gpd['sigma']:.4f}")


# ---------------------------------------------------------------------------
# 5. fp16 == fp32
# ---------------------------------------------------------------------------
def test_precision():
    print("\n[5] fp16 vs fp32 distances")
    rng = np.random.default_rng(41)
    X = rng.standard_normal((64, 48)).astype(np.float32)
    a, b = gpu_mi(X, half=True), gpu_mi(X, half=False)
    err = float(np.abs(a - b).max())
    check("fp16 KSG is exact", err < 1e-5, f"max diff = {err:.2e}")


if __name__ == "__main__":
    print(f"device: {DEV}  torch {torch.__version__}")
    test_reference()
    test_analytic()
    test_nonlinear()
    test_rectangular()
    test_paired()
    test_null()
    test_precision()
    print("\n" + ("ALL CHECKS PASSED" if not FAILURES
                  else f"{len(FAILURES)} FAILED: {', '.join(FAILURES)}"))
    sys.exit(1 if FAILURES else 0)
