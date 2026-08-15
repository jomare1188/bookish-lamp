#!/usr/bin/env python
# ============================================================================
# validate.py — correctness tests for 02_network_engine.py
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
#   6 PEARSON     The linear layer: the GPU kernel reproduces numpy's float64
#                 corrcoef, the analytic t-null reproduces scipy, the SIGN of
#                 negative correlations survives the shard round-trip, and the
#                 rank transform is genuinely bypassed (Pearson on ranks would
#                 silently be Spearman).
#
# RUN: /home/genomics/miniconda3/envs/docling/bin/python validate.py
# ============================================================================

import importlib.util
import math
import os
import sys

import numpy as np
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "cmi", os.path.join(HERE, "02_network_engine.py"))
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


# ---------------------------------------------------------------------------
# 6. the Pearson layer
# ---------------------------------------------------------------------------
def test_pearson():
    print("\n[6] pearson estimator (the linear layer)")
    rng = np.random.default_rng(7)
    n = 48
    X = rng.standard_normal((200, n)).astype(np.float32)
    # plant some strong pairs, half of them ANTI-correlated -- the sign is the
    # part most likely to be lost, because the sweep thresholds on |r|
    for a, b, sgn in ((0, 1, +1), (2, 3, -1), (4, 5, -1), (6, 7, +1)):
        X[b] = sgn * X[a] + 0.15 * rng.standard_normal(n).astype(np.float32)

    est = cmi.ESTIMATORS["pearson"](n=n, device=DEV)

    check("uses_ranks is False", est.uses_ranks is False,
          "Pearson on ranks would be Spearman")

    p = est.prep(torch.from_numpy(X).to(DEV))
    got = est.block(p, p).float().cpu().numpy()
    ref = np.corrcoef(X.astype(np.float64))
    err = float(np.abs(got - ref).max())
    check("GPU r == numpy float64 corrcoef", err < 1e-5, f"max diff = {err:.2e}")

    # the planted anti-correlations must come back NEGATIVE
    neg = [got[2, 3], got[4, 5]]
    check("negative correlations keep their sign",
          all(v < -0.9 for v in neg), f"r = {neg[0]:.4f}, {neg[1]:.4f}")

    # magnitude() is what the sweep thresholds on
    m = est.magnitude(torch.tensor([-0.95, 0.5]))
    check("magnitude() is abs()", bool((m > 0).all()) and abs(float(m[0]) - 0.95) < 1e-6)

    # analytic null vs scipy, and the round trip that sets the network threshold
    from scipy import stats
    pval_fn, bracket = est.analytic_pval(n)
    rs = np.array([0.3, 0.5, 0.7, 0.8, 0.9, -0.85])
    t = np.abs(rs) * np.sqrt((n - 2) / (1 - rs * rs))
    exp = 2 * stats.t.sf(t, df=n - 2)
    rel = float(np.max(np.abs(pval_fn(rs) - exp) / exp))
    check("analytic p == scipy t-test", rel < 1e-10, f"max rel err = {rel:.2e}")
    check("p is computed on |r|", abs(pval_fn(np.array([-0.85]))[0]
                                     - pval_fn(np.array([0.85]))[0]) < 1e-300)
    check("bracket spans the statistic", bracket == (0.0, 1.0))

    # inverting the null at p(|r|=0.8) must land back on exactly 0.8 -- this is
    # the check run.sh relies on when it passes --match-pearson STAT_MIN
    target = float(2 * stats.t.sf(0.8 * np.sqrt((n - 2) / (1 - 0.64)), df=n - 2))
    lo, hi = bracket
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if pval_fn(np.array([mid]))[0] > target:
            lo = mid
        else:
            hi = mid
    floor = 0.5 * (lo + hi)
    check("--match-pearson 0.8 round-trips to |r| = 0.8",
          abs(floor - 0.8) < 1e-9, f"got {floor:.12f}")

    # end-to-end: the signed value must survive sweep -> shard -> assemble
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        shard_dir = os.path.join(td, "s")
        genes = [f"g{i}" for i in range(200)]
        data = torch.from_numpy(X).to(DEV)
        cmi.sweep(est, data, 200, 0.8, shard_dir, 64, False, DEV, 10**9,
                  cut_max=0.9999)
        out = os.path.join(td, "e.tsv")
        n_w, cut_val, _ = cmi.assemble(shard_dir, genes, pval_fn,
                                       200 * 199 / 2, 0.05, out, "pearson",
                                       64.0, floor=0.8, magnitude=est.magnitude)
        rows = [l.split("\t") for l in open(out).read().splitlines()[1:]]
        vals = {(r[0], r[1]): float(r[2]) for r in rows}
        check("planted +r edge survives", vals.get(("g0", "g1"), 0) > 0.9)
        check("planted -r edge survives with sign",
              vals.get(("g2", "g3"), 0) < -0.9,
              f"r = {vals.get(('g2','g3'))}")
        check("no edge below the floor was written",
              all(abs(v) >= 0.8 for v in vals.values()))
        check("no edge above max-value was written",
              all(abs(v) <= 0.9999 for v in vals.values()))
        # p-values must rank by |r|, not by r
        srt = sorted(vals.items(), key=lambda kv: -abs(kv[1]))
        pv = {(r[0], r[1]): float(r[3]) for r in rows}
        ok = all(pv[srt[i][0]] <= pv[srt[i + 1][0]] * (1 + 1e-9)
                 for i in range(len(srt) - 1))
        check("p-values are monotone in |r| (not in r)", ok)


# ---------------------------------------------------------------------------
# 7. the gene-vs-trait estimator (Ross 2014, tie-corrected)
# ---------------------------------------------------------------------------
def test_gene_trait():
    print("\n[7] gene-trait MI (Ross, discrete trait)")
    import importlib.util
    sp = importlib.util.spec_from_file_location(
        "gtm", os.path.join(HERE, "12_gene_trait_mi.py"))
    gtm = importlib.util.module_from_spec(sp); sp.loader.exec_module(gtm)

    from scipy.special import digamma as psi

    def naive(lab, k=3):
        """Deliberately obvious reference: loops, no batching, no torch."""
        lab = np.asarray(lab); n = len(lab)
        tot = 0.0; ktot = 0.0; nctot = 0.0
        for i in range(n):
            d_same = sorted(abs(j - i) for j in range(n) if j != i and lab[j] == lab[i])
            d = d_same[k - 1]
            ki = sum(1 for j in range(n) if j != i and lab[j] == lab[i] and abs(j - i) <= d)
            mi = sum(1 for j in range(n) if j != i and abs(j - i) <= d)
            ktot += psi(ki); tot += psi(mi)
            nctot += psi(np.sum(lab == lab[i]))
        return psi(n) - nctot / n + ktot / n - tot / n

    def batched(lab, k=3):
        lab = np.asarray(lab)
        cs = torch.from_numpy(np.bincount(lab).astype(np.int64)).to(DEV)
        L = torch.from_numpy(lab[None, :]).to(DEV).long()
        return float(gtm.ross_batch(L, cs, k)[0])

    rng = np.random.default_rng(11)
    # batched == naive, on random and on structured label sequences
    worst = 0.0
    for n, C in ((48, 2), (18, 3), (30, 3)):
        base = np.repeat(np.arange(C), n // C)
        for _ in range(15):
            lab = rng.permutation(base)
            worst = max(worst, abs(batched(lab) - naive(lab)))
        worst = max(worst, abs(batched(base) - naive(base)))       # perfectly separated
    check("batched == naive reference", worst < 1e-5, f"max diff = {worst:.2e}")

    # the ceiling: perfect separation must approach H(trait) FROM BELOW
    ok = True; detail = []
    for n, C in ((48, 2), (18, 3)):
        base = np.repeat(np.arange(C), n // C)
        I = batched(base); H = math.log(C)
        detail.append(f"n={n} C={C}: {I:.3f} vs ln{C}={H:.3f}")
        ok &= (0.7 * H < I <= H)
    check("perfect separation approaches H(trait) from below", ok, "; ".join(detail))

    # a trait that is constant within class but the classes interleave perfectly
    # carries no information about position -> should sit near the null
    inter = np.array([i % 2 for i in range(48)])
    check("interleaved labels give ~0", abs(batched(inter)) < 0.25,
          f"I = {batched(inter):+.4f}")

    # null calibration: p-values from the shared null must be uniform under H0
    lab = np.repeat(np.arange(2), 24)
    cs = np.bincount(lab).astype(np.int64)
    null = gtm.build_null(lab, cs, 3, 200_000, DEV, seed=5)
    gpd = cmi.fit_gpd_tail(null)
    pf = cmi.make_pvalue_fn(null, gpd)
    draw = gtm.build_null(lab, cs, 3, 100_000, DEV, seed=99)
    for a in (0.05, 0.01, 0.001):
        obs = float((pf(draw) <= a).mean())
        check(f"P(p <= {a}) under H0 == {a}", abs(obs - a) <= 0.35 * a + 2e-4,
              f"observed {obs:.5f}")

    # the null must not depend on the gene -- that is what licenses sharing it.
    # two different random "genes" ranked over the same labels give statistics
    # drawn from the same distribution, by construction; verify the statistic
    # really is a pure function of the label sequence.
    labA = rng.permutation(np.repeat(np.arange(2), 24))
    same = batched(labA) == batched(labA.copy())
    check("statistic is a pure function of the label sequence", same)


if __name__ == "__main__":
    print(f"device: {DEV}  torch {torch.__version__}")
    test_reference()
    test_analytic()
    test_nonlinear()
    test_rectangular()
    test_paired()
    test_null()
    test_precision()
    test_pearson()
    test_gene_trait()
    print("\n" + ("ALL CHECKS PASSED" if not FAILURES
                  else f"{len(FAILURES)} FAILED: {', '.join(FAILURES)}"))
    sys.exit(1 if FAILURES else 0)
