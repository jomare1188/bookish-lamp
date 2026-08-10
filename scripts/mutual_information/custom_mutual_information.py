#!/usr/bin/env python
# ============================================================================
# custom_mutual_information.py — GPU all-pairs dependence networks
#
# A drop-in alternative to scripts/pearson_cor.r + scripts/build_edgelist.r
# that measures NON-LINEAR dependence instead of linear correlation, and goes
# straight from the expression matrix to a thresholded edgelist without ever
# materialising the n_genes x n_genes matrix.
#
# ---------------------------------------------------------------------------
# WHAT IT COMPUTES  (--estimator)
#
#   ksg   Kraskov-Stogbauer-Grassberger mutual information, algorithm 1
#         (Phys Rev E 69:066138, 2004). k-nearest-neighbour estimator in the
#         joint 2-D space, Chebyshev metric. This is the actual MI, in nats,
#         and it sees ANY dependence: non-monotone, multi-modal, thresholded.
#         The expensive one.
#
#   gcmi  Gaussian-copula MI (Ince et al., Hum Brain Mapp 38:1541, 2017).
#         Rank -> normal scores, then MI = -0.5 * log(1 - r^2). It is a proven
#         LOWER BOUND on the true MI and costs one matrix multiply, but being
#         a rank correlation it only sees MONOTONE dependence. Included as a
#         fast baseline and as a sanity check on the KSG run: gcmi <= ksg
#         should hold for monotone pairs, and pairs where ksg >> gcmi are
#         exactly the non-linear relationships Pearson is missing.
#
#   xi    Chatterjee's correlation coefficient (JASA 116:2009, 2021). Not an
#         MI, but the cheapest way to detect non-monotone dependence: it -> 1
#         whenever Y is a measurable function of X, monotone or not, and has a
#         closed-form rank expression. Asymmetric, so the reported value is
#         max(xi(x->y), xi(y->x)).
#
# ---------------------------------------------------------------------------
# WHY EVERYTHING IS RANK-TRANSFORMED FIRST
#
# Each gene is replaced by its ranks 0..n-1 across libraries, ties broken
# uniformly at random. Three things follow, and the whole design rests on them:
#
#   1. MI is invariant under strictly monotone transforms of the marginals, so
#      ranking costs no information. It also removes the VST's residual
#      mean-variance trend, which otherwise leaks into a kNN estimator.
#   2. Distinct ranks mean no zero joint distances, so the KSG estimator never
#      hits its degenerate tie case.
#   3. THE NULL BECOMES DATA-INDEPENDENT. Under independence, a pair of rank
#      vectors is just a random permutation pair, so the null distribution of
#      the statistic depends ONLY on (n, k) — not on the genes. That means the
#      permutation null can be built ONCE, from --n-perm permutations, and
#      reused for all 1.5e10 pairs. Per-pair permutation testing at this scale
#      is otherwise completely out of reach.
#
# ---------------------------------------------------------------------------
# HOW IT USES THE GPU AND THE RAM
#
# GPU: the KSG kernel is blocked into (a x b) gene tiles and evaluated as one
# (a, b, n, n) tensor of joint Chebyshev distances in fp16 -- exact here, since
# every value is a small integer rank difference. --gpu-mem-gb sets the tile
# size; the loop backs off and retries at half size on CUDA OOM. Measured on an
# RTX A4500: ~1.7 M pairs/s at n = 48.
#
# RAM: the expression matrix is 33 MB and lives on the GPU for the whole run.
# The host never holds anything of order n_genes^2. Surviving edges are streamed
# to per-tile 12-byte binary shards (int32 gene_i, int32 gene_j, float32 value)
# and only assembled into a text edgelist at the very end. Shards double as
# checkpoints: --resume skips every tile that already has an .ok marker, so a
# multi-hour run survives a reboot.
#
# ---------------------------------------------------------------------------
# SAMPLE SIZE — READ THIS BEFORE TRUSTING THE OUTPUT
#
# kNN mutual information is data-hungry. At n = 48 (Munoz) the KSG estimator is
# usable but has a wide sampling distribution; at n = 18 (Kiet) it is close to
# the floor of what the estimator can do, and individual edge values should be
# read as "ranked evidence of dependence", never as calibrated nats. The shared
# permutation null keeps the FALSE-POSITIVE rate honest at any n -- that part is
# exact -- but POWER at n = 18 is low, and the non-linear relationships this is
# meant to find are precisely the ones that need the most samples. Expect the
# purple MI network to be sparser than its Pearson counterpart, and treat a
# purple-vs-sugarcane comparison of MI edge counts as confounded by n.
#
# ---------------------------------------------------------------------------
# OUTPUT
#
#   <out>.edgelist.tsv   gene1 gene2 <estimator> pval padj   -- same column
#                        layout as edgelist_*_pearson.tsv, so build_edgelist.r's
#                        consumers (comparative_networks2.r, mcl_clustering.r,
#                        network_conservation_join.r) work unchanged.
#   <out>.null.tsv       the shared null: quantiles + GPD tail fit
#   <out>.summary.json   parameters, timings, edge counts, diagnostics
#
# RUN  (normally through run_all.sh, which fills these in from config.sh)
#   python custom_mutual_information.py \
#       --matrix  files/sugarcane/mi/sugarcane \
#       --out     files/sugarcane/mi/sugarcane_ksg \
#       --estimator ksg --k 3 --match-pearson 0.8 --resume
# ============================================================================

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time

import numpy as np
import torch

EDGE_DTYPE = np.dtype([("i", np.int32), ("j", np.int32), ("v", np.float32)])
P_FLOOR = 1e-300

# Auto-sized null: aim for this many permutations to exceed the strictest
# threshold, so it is read off the empirical null rather than extrapolated.
PERM_PER_EXCEED = 30.0

# The candidate cut is held at least this many times looser than the
# effect-size floor, so BH has room below the floor and the floor itself is
# never the weakest stored value.
CAND_P_HEADROOM = 20.0


# ===========================================================================
# input
# ===========================================================================
def load_matrix(prefix: str):
    with open(prefix + ".meta.json") as fh:
        meta = json.load(fh)
    with open(prefix + ".genes.txt") as fh:
        genes = [ln.rstrip("\n") for ln in fh]
    g, n = int(meta["n_genes"]), int(meta["n_samples"])
    if len(genes) != g:
        raise SystemExit(f"gene list ({len(genes)}) disagrees with meta ({g})")
    x = np.fromfile(prefix + ".f32", dtype=np.float32)
    if x.size != g * n:
        raise SystemExit(f"{prefix}.f32 has {x.size} values, expected {g*n}")
    return x.reshape(g, n), genes, meta


def rank_transform(x: np.ndarray, seed: int):
    """Per-gene ranks 0..n-1 with ties broken uniformly at random.

    Random (not average) tie-breaking is deliberate: it guarantees every gene is
    an exact permutation of 0..n-1, which is what makes the shared permutation
    null in build_null() valid for every pair.
    """
    g, n = x.shape
    rng = np.random.default_rng(seed)

    n_tied = int(np.sum([len(np.unique(row)) < n for row in x[: min(g, 20000)]]))
    tie_rate = n_tied / min(g, 20000)

    perm = rng.permuted(np.broadcast_to(np.arange(n), (g, n)).copy(), axis=1)
    xp = np.take_along_axis(x, perm, axis=1)
    order = np.argsort(xp, axis=1, kind="stable")       # ties -> random order
    ranks = np.empty((g, n), dtype=np.int32)
    np.put_along_axis(ranks, order,
                      np.broadcast_to(np.arange(n, dtype=np.int32), (g, n)), axis=1)
    out = np.empty((g, n), dtype=np.int32)
    np.put_along_axis(out, perm, ranks, axis=1)
    return out, tie_rate


# ===========================================================================
# estimators — each takes rank tiles on the GPU and returns an (a, b) matrix
# ===========================================================================
class KSG:
    """Kraskov algorithm 1 on rank data, Chebyshev metric.

        I = psi(k) + psi(N) - < psi(nx+1) + psi(ny+1) >

    fp16 is exact for the distances (integer rank differences, |d| < 2048) and
    halves both the memory and the time of the (a, b, n, n) tensor.
    """

    name = "ksg"
    symmetric = True

    def __init__(self, n: int, k: int, device, half: bool = True):
        if k >= n:
            raise SystemExit(f"--k {k} must be < n_samples ({n})")
        self.n, self.k, self.device = n, k, device
        self.dtype = torch.float16 if half else torch.float32
        # psi_tab[m] = digamma(m + 1), so counts index it directly
        self.psi = torch.digamma(
            torch.arange(1, n + 2, dtype=torch.float32, device=device))
        self.const = float(self.psi[k - 1] + self.psi[n - 1])

    def bytes_per_pair(self) -> int:
        # (a,b,n,n): fp16 distances + topk workspace + two bool masks
        return self.n * self.n * 10

    def prep(self, ranks: torch.Tensor):
        """Pairwise within-gene |rank_i - rank_j|, shape (G, n, n)."""
        r = ranks.to(self.dtype)
        return (r[:, :, None] - r[:, None, :]).abs()

    def block(self, dx: torch.Tensor, dy: torch.Tensor) -> torch.Tensor:
        # dx (a,n,n), dy (b,n,n) -> joint distance (a,b,n,n)
        d = torch.maximum(dx[:, None], dy[None, :])
        eps = torch.topk(d, self.k + 1, dim=-1, largest=False).values[..., self.k]
        del d
        eps = eps.unsqueeze(-1)                                  # (a,b,n,1)
        # strict inequality; the self-distance 0 is always < eps, hence the -1
        nx = (dx[:, None] < eps).sum(-1) - 1                     # (a,b,n)
        ny = (dy[None, :] < eps).sum(-1) - 1
        del eps
        return self.const - (self.psi[nx] + self.psi[ny]).mean(-1)

    def block_paired(self, dx: torch.Tensor, dy: torch.Tensor) -> torch.Tensor:
        """MI of MATCHED pairs only: (dx[s], dy[s]) for each s. Shape (m,).

        Same estimator as block(), minus the outer product. Taking the diagonal
        of an m x m block instead would compute m^2 pairs to use m of them.
        """
        d = torch.maximum(dx, dy)                                # (m,n,n)
        eps = torch.topk(d, self.k + 1, dim=-1,
                         largest=False).values[..., self.k].unsqueeze(-1)
        del d
        nx = (dx < eps).sum(-1) - 1
        ny = (dy < eps).sum(-1) - 1
        del eps
        return self.const - (self.psi[nx] + self.psi[ny]).mean(-1)


class GCMI:
    """Gaussian-copula MI: normal scores, then -0.5 * log(1 - r^2)."""

    name = "gcmi"
    symmetric = True

    def __init__(self, n: int, device, **_):
        self.n, self.device = n, device
        self.dtype = torch.float32

    def bytes_per_pair(self) -> int:
        return 16

    def prep(self, ranks: torch.Tensor):
        z = torch.special.ndtri((ranks.float() + 0.5) / self.n)
        z = z - z.mean(1, keepdim=True)
        return z / z.norm(dim=1, keepdim=True).clamp_min(1e-12)   # unit rows

    def block(self, zx: torch.Tensor, zy: torch.Tensor) -> torch.Tensor:
        r = (zx @ zy.T).clamp(-0.999999, 0.999999)
        return -0.5 * torch.log1p(-r * r)

    def block_paired(self, zx: torch.Tensor, zy: torch.Tensor) -> torch.Tensor:
        r = (zx * zy).sum(1).clamp(-0.999999, 0.999999)
        return -0.5 * torch.log1p(-r * r)


class XI:
    """Chatterjee's xi, both directions, reported as the maximum.

        xi(x->y) = 1 - 3 * sum |r_{i+1} - r_i| / (n^2 - 1)

    with r the y-ranks read in x-sorted order. No ties (ranks are a permutation),
    so the simplified continuous formula applies exactly.
    """

    name = "xi"
    symmetric = False

    def __init__(self, n: int, device, **_):
        self.n, self.device = n, device
        self.denom = float(n * n - 1)
        self.dtype = torch.float32

    def bytes_per_pair(self) -> int:
        return self.n * 24

    def prep(self, ranks: torch.Tensor):
        # keep both the ranks and the order that sorts each gene
        return (ranks.float(), ranks.argsort(dim=1))

    def block(self, px, py) -> torch.Tensor:
        (rx, ox), (ry, oy) = px, py
        n = self.n

        def one_way(order_src, rank_dst):
            # shapes are taken from the arguments, NOT closed over: the second
            # call passes them the other way round, and the final tile of a
            # sweep is not square
            m1, m2 = order_src.shape[0], rank_dst.shape[0]
            idx = order_src[:, None, :].expand(m1, m2, n)
            vals = rank_dst[None, :, :].expand(m1, m2, n).gather(-1, idx)
            s = (vals[..., 1:] - vals[..., :-1]).abs().sum(-1)
            return 1.0 - 3.0 * s / self.denom

        return torch.maximum(one_way(ox, ry), one_way(oy, rx).transpose(0, 1))

    def block_paired(self, px, py) -> torch.Tensor:
        (rx, ox), (ry, oy) = px, py

        def one_way(order_src, rank_dst):
            vals = rank_dst.gather(-1, order_src)               # (m,n)
            s = (vals[:, 1:] - vals[:, :-1]).abs().sum(-1)
            return 1.0 - 3.0 * s / self.denom

        return torch.maximum(one_way(ox, ry), one_way(oy, rx))


ESTIMATORS = {"ksg": KSG, "gcmi": GCMI, "xi": XI}


# ===========================================================================
# the shared null
# ===========================================================================
def build_null(est, n: int, n_perm: int, device, seed: int, chunk_pairs: int):
    """Statistic under independence, from random rank-permutation pairs.

    Valid for every gene pair because after ranking, an independent pair IS a
    random permutation pair — see the header. Returns the sorted null sample.
    """
    rng = np.random.default_rng(seed + 977)
    out, done = [], 0
    t0 = time.time()
    while done < n_perm:
        m = min(chunk_pairs, n_perm - done)
        a = np.argsort(rng.random((m, n)), axis=1).astype(np.int32)
        b = np.argsort(rng.random((m, n)), axis=1).astype(np.int32)
        ta = est.prep(torch.from_numpy(a).to(device))
        tb = est.prep(torch.from_numpy(b).to(device))
        # only the matched diagonal of the m x m block is an independent draw
        vals = _diag_block(est, ta, tb, m)
        out.append(vals.float().cpu().numpy())
        done += m
        del ta, tb, vals
        torch.cuda.empty_cache()
    null = np.sort(np.concatenate(out))
    print(f"[null] {len(null):,} permutations in {time.time()-t0:.1f}s  "
          f"mean {null.mean():+.5f}  sd {null.std():.5f}  "
          f"q99.9 {np.quantile(null, 0.999):.5f}  max {null[-1]:.5f}", flush=True)
    return null


def _take(prepped, sl):
    return tuple(t[sl] for t in prepped) if isinstance(prepped, tuple) else prepped[sl]


def _diag_block(est, ta, tb, m: int, step: int = 65536):
    """Matched pairs (a_s, b_s) only, via each estimator's paired kernel."""
    parts = []
    for s in range(0, m, step):
        e = min(s + step, m)
        parts.append(est.block_paired(_take(ta, slice(s, e)), _take(tb, slice(s, e))))
    return torch.cat(parts)


def fit_gpd_tail(null: np.ndarray, tail_q: float = 0.99):
    """Probability-weighted-moment GPD fit to the null's upper tail.

    Needed because BH over ~1.5e10 tests rejects at p ~ 1e-12, which no
    affordable number of permutations can resolve empirically. The empirical
    ECDF is used down to p = 1 - tail_q; below that, exceedances over the
    threshold u are extrapolated with a generalised Pareto tail. Hosking &
    Wallis (1987) PWM estimator; no scipy needed.
    """
    u = float(np.quantile(null, tail_q))
    y = null[null > u] - u
    ny = y.size
    if ny < 200:
        return dict(u=u, p_u=ny / null.size, xi=0.0, sigma=float(max(y.mean(), 1e-12)),
                    n_tail=ny, method="exponential (too few exceedances)")
    ys = np.sort(y)
    a0 = ys.mean()
    a1 = float(np.mean(ys * (1.0 - (np.arange(1, ny + 1) - 0.35) / ny)))
    den = a0 - 2.0 * a1
    if abs(den) < 1e-15:
        xi, sigma = 0.0, a0
        method = "exponential (degenerate PWM)"
    else:
        k = a0 / den - 2.0                     # Hosking's k
        xi, sigma = -k, 2.0 * a0 * a1 / den    # standard GPD parametrisation
        method = "GPD-PWM"
        if sigma <= 0:
            xi, sigma, method = 0.0, a0, "exponential (non-positive scale)"
    return dict(u=u, p_u=ny / null.size, xi=float(xi), sigma=float(sigma),
                n_tail=int(ny), method=method)


def gpd_endpoint(gpd: dict) -> float:
    """Upper endpoint of the fitted tail (inf unless the shape is negative).

    A negative GPD shape means the fitted null has a hard maximum, so any edge
    above it gets p = 0, floored to P_FLOOR. That is harmless for the network --
    such edges are kept by any threshold -- but their reported p/padj are
    "smaller than the fit can resolve", not calibrated numbers. The count is
    written to the summary so it is never silently assumed otherwise.
    """
    xi, sigma, u = gpd["xi"], gpd["sigma"], gpd["u"]
    return u - sigma / xi if xi < 0 else float("inf")


def make_pvalue_fn(null: np.ndarray, gpd: dict):
    """Upper-tail p-value: empirical ECDF, GPD-extrapolated past the threshold."""
    n_null = null.size
    u, p_u, xi, sigma = gpd["u"], gpd["p_u"], gpd["xi"], gpd["sigma"]

    def pval(v: np.ndarray) -> np.ndarray:
        v = np.asarray(v, dtype=np.float64)
        # empirical: (#null >= v + 1) / (n + 1), the conservative convention
        emp = (n_null - np.searchsorted(null, v, side="left") + 1.0) / (n_null + 1.0)
        out = emp
        hi = v > u
        if np.any(hi):
            z = v[hi] - u
            if abs(xi) < 1e-8:
                tail = p_u * np.exp(-z / sigma)
            else:
                base = 1.0 + xi * z / sigma
                tail = np.where(base > 0, p_u * np.power(np.maximum(base, 1e-300),
                                                         -1.0 / xi), 0.0)
            out = out.copy()
            out[hi] = tail
        return np.clip(out, P_FLOOR, 1.0)

    return pval


# ===========================================================================
# the all-pairs sweep
# ===========================================================================
def choose_tile(est, budget_gb: float, n_genes: int) -> int:
    bpp = est.bytes_per_pair()
    t = int(math.sqrt(budget_gb * 1e9 / bpp))
    return max(32, min(t, n_genes, 8192))


def sweep(est, ranks_gpu, n_genes, cut, shard_dir, tile, resume, device, log_every):
    """Upper triangle, tile by tile, writing candidate edges to binary shards."""
    os.makedirs(shard_dir, exist_ok=True)
    row_tiles = list(range(0, n_genes, tile))
    total_pairs = n_genes * (n_genes - 1) / 2
    done_pairs, kept, t0 = 0.0, 0, time.time()

    for ti, r0 in enumerate(row_tiles):
        r1 = min(r0 + tile, n_genes)
        shard = os.path.join(shard_dir, f"tile_{ti:06d}.bin")
        okf = shard + ".ok"
        tile_pairs = sum(max(0, n_genes - max(r0, i) - 1) for i in range(r0, r1))

        if resume and os.path.exists(okf):
            done_pairs += tile_pairs
            kept += os.path.getsize(shard) // EDGE_DTYPE.itemsize
            continue

        px = est.prep(ranks_gpu[r0:r1])
        buf = []
        c0 = r0                                   # upper triangle only
        while c0 < n_genes:
            c1 = min(c0 + tile, n_genes)
            py = est.prep(ranks_gpu[c0:c1])
            v = est.block(px, py).float()
            del py

            mask = v >= cut
            if c0 == r0:                          # diagonal tile: keep j > i
                gi = torch.arange(r0, r1, device=device)[:, None]
                gj = torch.arange(c0, c1, device=device)[None, :]
                mask &= gj > gi
            idx = mask.nonzero(as_tuple=False)
            if idx.numel():
                rec = np.empty(idx.shape[0], dtype=EDGE_DTYPE)
                rec["i"] = (idx[:, 0] + r0).cpu().numpy()
                rec["j"] = (idx[:, 1] + c0).cpu().numpy()
                rec["v"] = v[mask].cpu().numpy()
                buf.append(rec)
            del v, mask, idx
            c0 = c1

        with open(shard, "wb") as fh:
            if buf:
                arr = np.concatenate(buf)
                arr.tofile(fh)
                kept += arr.size
        open(okf, "w").close()
        del px, buf
        torch.cuda.empty_cache()

        done_pairs += tile_pairs
        if (ti + 1) % log_every == 0 or r1 == n_genes:
            el = time.time() - t0
            frac = done_pairs / total_pairs
            eta = el / frac - el if frac > 0 else 0.0
            print(f"[sweep] tile {ti+1}/{len(row_tiles)}  {frac*100:5.1f}%  "
                  f"{done_pairs/max(el,1e-9)/1e6:5.2f} Mpairs/s  "
                  f"kept {kept:,}  elapsed {el/60:.1f}m  ETA {eta/60:.1f}m",
                  flush=True)
    return kept, total_pairs


# ===========================================================================
# assembly: BH over ALL pairs, then the text edgelist
# ===========================================================================
def assemble(shard_dir, genes, pval_fn, n_tests, alpha, out_tsv, stat_name,
             ram_limit_gb, floor=-np.inf):
    shards = sorted(f for f in os.listdir(shard_dir) if f.endswith(".bin"))
    n_edges = sum(os.path.getsize(os.path.join(shard_dir, f)) for f in shards) \
        // EDGE_DTYPE.itemsize
    if n_edges == 0:
        print("[bh] no candidate edges survived the value cut")
        open(out_tsv, "w").write(f"gene1\tgene2\t{stat_name}\tpval\tpadj\n")
        return 0, float("nan")

    need_gb = n_edges * 8 / 1e9
    if need_gb > ram_limit_gb:
        raise SystemExit(
            f"{n_edges:,} candidate edges need ~{need_gb:.1f} GB to FDR-correct "
            f"(limit {ram_limit_gb} GB). Re-run with a higher --min-value.")

    # --- pass A: all candidate statistics, sorted descending ----------------
    vals = np.empty(n_edges, dtype=np.float32)
    at = 0
    for f in shards:
        a = np.fromfile(os.path.join(shard_dir, f), dtype=EDGE_DTYPE)
        vals[at:at + a.size] = a["v"]
        at += a.size
    vals.sort()
    vals = vals[::-1].copy()                       # descending == p ascending

    # --- BH over ALL m pairs, not just the candidates -----------------------
    # Exact despite only holding the candidates, because every pair that was NOT
    # written to a shard has a p LARGER than every pair that was, so a
    # candidate's rank among the candidates equals its rank among all m tests.
    #
    # The one way this breaks is if BH's cut lands on the candidate boundary --
    # then edges that BH would have kept were never stored. That is exactly what
    # happens if the candidate cut is set at or below the Bonferroni value, since
    # BH always rejects a SUPERSET of Bonferroni. Hence --cand-p defaults loose
    # and the boundary is checked explicitly below.
    p_sorted = pval_fn(vals)
    ranks = np.arange(1, n_edges + 1, dtype=np.float64)
    padj = np.minimum.accumulate((n_tests * p_sorted / ranks)[::-1])[::-1]
    np.clip(padj, 0.0, 1.0, out=padj)
    n_sig = int(np.searchsorted(padj, alpha, side="right"))

    bh_val = float(vals[n_sig - 1]) if n_sig else float("inf")
    cut_val = max(bh_val, floor)          # BH and the effect-size floor, both
    n_keep = int(np.searchsorted(-vals, -np.float32(cut_val), side="right"))

    print(f"[bh] {n_edges:,} candidates | {n_tests:,.0f} tests | "
          f"BH alpha={alpha} -> {n_sig:,} edges ({stat_name} >= {bh_val:.5f})")
    if floor > -np.inf:
        print(f"[bh] effect-size floor {stat_name} >= {floor:.5f} -> "
              f"{n_keep:,} edges kept (the binding constraint is "
              f"{'the floor' if floor > bh_val else 'BH'})")
    n_floor = int(np.sum(p_sorted[:n_keep] <= P_FLOOR * 10))
    if n_floor:
        print(f"[bh] {n_floor:,} of them ({100*n_floor/max(n_keep,1):.2f}%) exceed the "
              f"fitted null's upper endpoint; their p/padj are reported at the "
              f"floor and mean 'beyond the fit', not a calibrated value.")
    # Truncation is only real if the FINAL cut reaches the candidate boundary.
    # BH saturating the candidate set is harmless whenever the effect-size floor
    # sits well above that boundary -- then the floor, not the FDR, defines the
    # network, and everything the floor admits was stored.
    if n_keep >= n_edges and n_edges:
        print("[bh] WARNING: the final cut reaches the weakest stored candidate, "
              "so edges that belong in the network were never written. The "
              "network is TRUNCATED — re-run with a larger --cand-p.")
    elif n_sig >= n_edges and n_edges:
        print(f"[bh] note: BH alone would keep every stored candidate, so the FDR "
              f"is not the binding constraint — the effect-size floor is, and it "
              f"sits {floor - float(vals[-1]):.3f} above the weakest candidate, so "
              f"the network is complete. padj is reported for reference.")

    # --- pass B: stream shards again, write only the significant edges ------
    # Row formatting, not the statistics, dominates this pass at 1e8+ edges, so
    # gene ids are looked up with one vectorised fancy-index and the rows are
    # built with "".join over a %-format map rather than per-row f-strings.
    gname = np.asarray(genes, dtype=object)
    fmt = "%s\t%s\t%.5f\t%.4g\t%.4g\n".__mod__
    CH = 1 << 20
    with open(out_tsv, "w", buffering=1 << 23) as fh:
        fh.write(f"gene1\tgene2\t{stat_name}\tpval\tpadj\n")
        written = 0
        for f in shards:
            a = np.fromfile(os.path.join(shard_dir, f), dtype=EDGE_DTYPE)
            if a.size == 0:
                continue
            a = a[a["v"] >= cut_val]
            if a.size == 0:
                continue
            pos = np.clip(np.searchsorted(-vals, -a["v"], side="left"),
                          0, n_edges - 1)
            g1, g2 = gname[a["i"]], gname[a["j"]]
            pv, pa = p_sorted[pos], padj[pos]
            for s in range(0, a.size, CH):
                e = min(s + CH, a.size)
                fh.write("".join(map(fmt, zip(g1[s:e], g2[s:e], a["v"][s:e],
                                              pv[s:e], pa[s:e]))))
            written += a.size
    print(f"[out] {written:,} edges -> {out_tsv}")
    return written, cut_val, n_floor


# ===========================================================================
def main():
    ap = argparse.ArgumentParser(
        description="GPU all-pairs mutual-information / dependence networks",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--matrix", required=True,
                    help="prefix written by 00_export_vst.r")
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--estimator", default="ksg", choices=list(ESTIMATORS))
    ap.add_argument("--k", type=int, default=3, help="KSG neighbours")
    ap.add_argument("--n-perm", type=int, default=0,
                    help="permutations for the shared null; 0 = auto-size it to "
                         "resolve the strictest threshold empirically")
    ap.add_argument("--max-perm", type=int, default=200_000_000,
                    help="ceiling for the auto-sized null (4 bytes each in RAM; "
                         "2e8 costs ~9 min at n=48 and 0.8 GB)")
    ap.add_argument("--alpha", type=float, default=0.05, help="BH FDR level")
    ap.add_argument("--cand-p", type=float, default=1e-3,
                    help="per-pair p-value below which an edge is written to a "
                         "shard. Sets how much disk the run costs and, more "
                         "importantly, the range over which BH stays exact: BH "
                         "may safely reject up to cand_p * n_tests / alpha edges")
    ap.add_argument("--match-pearson", type=float, default=0.8,
                    help="effect-size floor, given as the Pearson |r| whose "
                         "two-sided p at this n the floor should match. Mirrors "
                         "general_stats.r's PEARSON_MIN = 0.8 rule (NOT "
                         "build_edgelist.r's 0.7, which is only a pre-filter) "
                         "so the MI and Pearson networks are equally "
                         "stringent per edge. "
                         "0 disables the floor")
    ap.add_argument("--min-value", type=float, default=None,
                    help="effect-size floor given directly in the statistic's "
                         "own units; overrides --match-pearson")
    ap.add_argument("--gpu-mem-gb", type=float, default=5.0)
    ap.add_argument("--ram-limit-gb", type=float, default=24.0)
    ap.add_argument("--tile", type=int, default=0, help="0 = auto")
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--fp32", action="store_true",
                    help="fp32 KSG distances (slower; fp16 is exact for n<2048)")
    ap.add_argument("--seed", type=int, default=1188)
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--keep-shards", action="store_true")
    ap.add_argument("--limit-genes", type=int, default=0,
                    help="use only the first N genes (benchmarking)")
    ap.add_argument("--log-every", type=int, default=5, help="row tiles per log line")
    args = ap.parse_args()

    t_start = time.time()
    device = torch.device(args.device)
    if device.type == "cuda":
        if not torch.cuda.is_available():
            raise SystemExit("CUDA requested but not available")
        if device.index is None:                 # "cuda" -> "cuda:0"
            device = torch.device("cuda", torch.cuda.current_device())
        torch.cuda.set_device(device)
        print(f"[gpu ] {torch.cuda.get_device_name(device)}  "
              f"{torch.cuda.get_device_properties(device).total_memory/1e9:.1f} GB  "
              f"torch {torch.__version__}")

    # ---- data -----------------------------------------------------------
    x, genes, meta = load_matrix(args.matrix)
    if args.limit_genes:
        x, genes = x[: args.limit_genes], genes[: args.limit_genes]
    n_genes, n = x.shape
    print(f"[data] {meta.get('label','?')}: {n_genes:,} genes x {n} libraries")

    ranks, tie_rate = rank_transform(x, args.seed)
    print(f"[rank] ties broken at random in {tie_rate*100:.1f}% of genes "
          f"(sampled); marginals are now exact permutations of 0..{n-1}")
    if n < 30 and args.estimator == "ksg":
        print(f"[warn] n = {n}: the KSG estimator is at the low end of its "
              f"usable range. The permutation null keeps the false-positive "
              f"rate exact, but power is limited — see the script header.")

    ranks_gpu = torch.from_numpy(ranks).to(device)
    est = ESTIMATORS[args.estimator](n=n, k=args.k, device=device,
                                     half=not args.fp32)
    n_tests = n_genes * (n_genes - 1) / 2.0

    tile = args.tile or choose_tile(est, args.gpu_mem_gb, n_genes)
    print(f"[plan] estimator={est.name} tile={tile} "
          f"({tile*tile*est.bytes_per_pair()/1e9:.2f} GB/block) "
          f"pairs={n_tests:,.0f}")

    # ---- how deep into the tail do the thresholds actually reach? ---------
    # Both cuts are defined as p-values, and the strictest one determines how
    # many permutations the null needs. Computing it FIRST means the null can be
    # sized to resolve it empirically instead of extrapolating a long way.
    p_pearson = None
    if args.min_value is None and args.match_pearson and args.match_pearson > 0:
        from scipy import stats
        r = args.match_pearson
        t = r * math.sqrt((n - 2) / max(1.0 - r * r, 1e-12))
        p_pearson = float(2.0 * stats.t.sf(abs(t), df=n - 2))

    # The candidate cut must sit BELOW the effect-size floor with room to spare,
    # or the floor lands on the boundary and the network comes out truncated.
    # This bites at small n: |r| = 0.7 is p = 3.1e-8 at n = 48 but only 1.2e-3 at
    # n = 18, which is looser than the default candidate cut.
    cand_p = args.cand_p
    if p_pearson is not None and cand_p < CAND_P_HEADROOM * p_pearson:
        cand_p = CAND_P_HEADROOM * p_pearson
        print(f"[cut ] |r| = {args.match_pearson} at n = {n} is p = {p_pearson:.3g}; "
              f"raising --cand-p {args.cand_p:g} -> {cand_p:.3g} so the floor "
              f"clears the candidate boundary")
    target_p = min(cand_p, p_pearson) if p_pearson else cand_p

    # ---- null ------------------------------------------------------------
    # the paired kernel costs ~n^2 bytes per permutation, not n^2 per PAIR of
    # genes, so the null can use a far bigger chunk than the sweep's tile
    if args.n_perm <= 0:                          # 'auto'
        n_perm = int(min(args.max_perm, max(1_000_000, PERM_PER_EXCEED / target_p)))
        print(f"[null] auto-sizing to resolve p = {target_p:.3g}: "
              f"{n_perm:,} permutations "
              f"({n_perm * target_p:.1f} expected exceedances)")
    else:
        n_perm = args.n_perm
    null = build_null(est, n, n_perm, device, args.seed,
                      chunk_pairs=int(max(4096, min(
                          262144, args.gpu_mem_gb * 1e9 / (est.bytes_per_pair() / 8 + 1)))))
    gpd = fit_gpd_tail(null)
    pval_fn = make_pvalue_fn(null, gpd)
    print(f"[null] tail fit: {gpd['method']}  u={gpd['u']:.5f} "
          f"xi={gpd['xi']:+.4f} sigma={gpd['sigma']:.5f} "
          f"(n_exceed={gpd['n_tail']:,})")
    depth = 1.0 / (n_perm * target_p)
    if depth > 3.0:
        print(f"[warn] the strictest threshold (p = {target_p:.3g}) is {depth:.0f}x "
              f"beyond what {n_perm:,} permutations resolve, so it comes from the "
              f"GPD extrapolation, not the empirical null. Measured seed-to-seed "
              f"spread of such a threshold is a few percent of the statistic; "
              f"raise --max-perm if the edge count matters to that precision.")

    # ---- the two thresholds ----------------------------------------------
    def value_at_p(target_p: float) -> float:
        """Invert the null: the statistic value with per-pair p = target_p."""
        lo, hi = float(null.min()) - 1.0, float(null.max()) + 10.0
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if pval_fn(np.array([mid]))[0] > target_p:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    cut = value_at_p(cand_p)
    bh_headroom = cand_p * n_tests / args.alpha
    print(f"[cut ] candidate {est.name} >= {cut:.5f} (p = {cand_p:g} per pair); "
          f"BH stays exact up to {bh_headroom:,.0f} rejected edges")

    if args.min_value is not None:
        floor = args.min_value
        floor_src = "--min-value"
    elif p_pearson is not None:
        # the statistic value whose null tail probability equals that of
        # |r| = --match-pearson under the Pearson null at this n
        floor = value_at_p(p_pearson)
        floor_src = f"--match-pearson {args.match_pearson} (p = {p_pearson:.3g} at n = {n})"
    else:
        floor, floor_src = -np.inf, "disabled"
    if floor > -np.inf:
        print(f"[cut ] effect-size floor {est.name} >= {floor:.5f}  [{floor_src}]")

    with open(args.out + ".null.tsv", "w") as fh:
        fh.write("quantile\tvalue\n")
        for q in (0.5, 0.9, 0.99, 0.999, 0.9999, 0.99999, 1.0):
            fh.write(f"{q}\t{float(np.quantile(null, q)):.6f}\n")
        for kk, vv in gpd.items():
            fh.write(f"gpd_{kk}\t{vv}\n")

    # ---- sweep -----------------------------------------------------------
    shard_dir = args.out + ".shards"
    kept, _ = sweep(est, ranks_gpu, n_genes, cut, shard_dir, tile,
                    args.resume, device, args.log_every)
    t_sweep = time.time() - t_start
    print(f"[sweep] done: {kept:,} candidate edges in {t_sweep/60:.1f} min")

    # ---- FDR + edgelist ---------------------------------------------------
    n_sig, cut_val, n_floor = assemble(shard_dir, genes, pval_fn, n_tests,
                                       args.alpha, args.out + ".edgelist.tsv",
                                       est.name, args.ram_limit_gb, floor)

    if not args.keep_shards:
        for f in os.listdir(shard_dir):
            os.remove(os.path.join(shard_dir, f))
        os.rmdir(shard_dir)

    summary = dict(
        matrix=args.matrix, label=meta.get("label"), estimator=est.name,
        k=args.k if est.name == "ksg" else None, n_genes=n_genes, n_samples=n,
        n_tests=n_tests, n_perm=args.n_perm, alpha=args.alpha,
        cand_p=cand_p, cand_p_requested=args.cand_p, effect_size_floor=None if floor == -np.inf else floor,
        effect_size_floor_source=floor_src,
        tie_rate_sampled=tie_rate, candidate_cut=cut, final_cut=cut_val,
        n_candidate_edges=int(kept), n_significant_edges=int(n_sig),
        edge_density=(n_sig / n_tests) if n_tests else None,
        n_edges_beyond_null_fit=int(n_floor),
        null_upper_endpoint=gpd_endpoint(gpd),
        gpd=gpd, tile=tile, fp16=not args.fp32,
        runtime_min=(time.time() - t_start) / 60.0,
        torch=torch.__version__,
        device=torch.cuda.get_device_name(device) if device.type == "cuda" else "cpu",
    )
    with open(args.out + ".summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"[done] {summary['runtime_min']:.1f} min total -> {args.out}.*")


if __name__ == "__main__":
    main()
