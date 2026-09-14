#!/usr/bin/env python3
# =============================================================================
# 62_conserved_edges_pearson.py -- edge-level conservation on the Pearson-only
# networks, for all edges and for the nitrogen-correlated subsets.
#
# WHY A NEW SCRIPT. 06_conservation_join.r answers the same question on the MERGED
# (Pearson + MI) graph and streams network_<study>_edges.tsv. Those tables were
# DELETED -- 70 GB whose only consumer was mcxload -- so 06 cannot run on the graph
# the analysis now uses. 06 and 13 stay untouched as the record of that analysis;
# this is the answer on the current one.
#
# THE INPUT IS THE mcxdump ALREADY ON DISK. 47_leiden_sweep.py dumped each graph as
# `<study>.pairs` (mcxdump --dump-pairs --dump-upper --no-loops): one row per
# undirected edge, `idx1 idx2 weight`, indices into `<study>.tab`. Verified line for
# line against the matrices -- 75,333,769 sugarcane and 675,955,918 purple -- so
# nothing has to be regenerated, and the edge WEIGHT comes along for free.
#
# EVERYTHING RUNS IN INTEGER INDEX SPACE. Gene names are never materialised until the
# conserved edges are written. That is what makes 676M edges against a 676M-edge
# adjacency affordable: the adjacency is one sorted int64 array (key = min<<32 | max,
# 5.4 GB) and a lookup is np.searchsorted. It is built ONCE and reused by the observed
# pass and by every null replicate, which is the argument 13_conservation_null.r:34
# already makes.
#
# ---------------------------------------------------------------------------
# WHAT THE LAYER BREAKDOWN BECAME
#
# The old summary's headline was "are MI edges conserved as often as Pearson edges?"
# -- vacuous on a single-layer graph. It is replaced by conservation against EDGE
# STRENGTH, in rank-based weight deciles, nearly free from the same stream and the
# question that matters here: do the strongest co-expression edges survive across
# species? Bins are rank-based so they hold equal numbers of edges, and each bin's
# weight RANGE is reported -- the [0.01, 1] rescaling is per-study, so a decile in
# sugarcane and a decile in purple do not stand for the same |r|.
# ---------------------------------------------------------------------------
#
# THREE STRATA, ONE PASS, because re-streaming 676M edges to ask a second question is
# exactly what to avoid:
#   all           every edge in the source network
#   resp_source   both endpoints nitrogen-responsive in the SOURCE species
#   resp_both     both endpoints responsive in BOTH species
#
# WHO IS RESPONSIVE IS READ FROM 61, NEVER RECOMPUTED.
# 61_conserved_blocked_nodes.r already decided it -- the blocked test, plus purple's
# quadratic tier -- and this reads its status tables. One definition, so the node and
# edge levels cannot drift apart.
#
# THE `_pearson` SUFFIX ON EVERY OUTPUT IS NOT COSMETIC.
# conserved_genes_<study>_FULL.txt is still read by `run.sh trait`, by
# 08_conserved_cor_genes.r and by 09_go_enrichment.r -- all merged-network stages.
# Overwriting it would silently re-point three superseded stages at a different graph
# and make their outputs a mixture of two analyses. Nothing here touches a _FULL path.
#
# RUN: through run.sh  ->  ./run.sh consedges sugarcane_to_purple
# =============================================================================

import csv
import os
import re
import sys
import time

import numpy as np
import pandas as pd

T0 = time.time()


def say(*a):
    print("[%s]" % time.strftime("%H:%M:%S"), *a, flush=True)


def env_req(name):
    v = os.environ.get(name)
    if not v:
        sys.exit(f"FATAL: {name} is not set -- launch through run.sh")
    return v


DIRECTION = env_req("CLEAN_DIRECTION")
WORK = env_req("CLEAN_WORK_DIR")
OUT_DIR = env_req("CLEAN_OUT_DIR")
ORTHOGROUPS = env_req("CLEAN_ORTHOGROUPS")
RESULTS = env_req("CLEAN_RESULTS")
CHUNK = int(float(os.environ.get("CLEAN_CONS_EDGE_CHUNK", 5_000_000)))
N_BINS = int(os.environ.get("CLEAN_CONS_WEIGHT_BINS", 10))
N_REPS = int(os.environ.get("CLEAN_NULL_REPS", 20))
N_SAMPLE = int(float(os.environ.get("CLEAN_NULL_SAMPLE", 5_000_000)))
SEED = int(os.environ.get("CLEAN_SEED", 1188))
SP = {"sugarcane": os.environ.get("CLEAN_OG_SPECIES_SUGARCANE",
                                  "sugarcane_one_transcript"),
      "purple": os.environ.get("CLEAN_OG_SPECIES_PURPLE",
                               "one_transcript_purple_proteins")}

if DIRECTION == "sugarcane_to_purple":
    A, B = "sugarcane", "purple"
elif DIRECTION == "purple_to_sugarcane":
    A, B = "purple", "sugarcane"
else:
    sys.exit("FATAL: CLEAN_DIRECTION must be sugarcane_to_purple or purple_to_sugarcane")

os.makedirs(OUT_DIR, exist_ok=True)
rng = np.random.default_rng(SEED)
say(f"edge-level conservation: {A} -> {B}   (Pearson-only graphs)")

KEY_SHIFT = np.int64(32)


def pack(i, j):
    """Canonical undirected key: min<<32 | max, so (i,j) and (j,i) collide by design."""
    i = i.astype(np.int64, copy=False)
    j = j.astype(np.int64, copy=False)
    return (np.minimum(i, j) << KEY_SHIFT) | np.maximum(i, j)


# =============================================================================
# 1. index maps
# =============================================================================
def read_tab(study):
    """<study>.tab is mcl's `index<TAB>name`, in index order. The order is ASSERTED:
    every downstream array is positional, so a gap would misname every gene."""
    path = os.path.join(WORK, f"{study}.tab")
    if not os.path.exists(path):
        sys.exit(f"FATAL: no {path} -- run ./run.sh pearsonmci {study}")
    names = []
    with open(path) as fh:
        for expect, line in enumerate(fh):
            i, nm = line.rstrip("\n").split("\t", 1)
            if int(i) != expect:
                sys.exit(f"FATAL: {path} is not in index order at line {expect}")
            names.append(nm)
    return np.array(names, dtype=object)


NAMES = {s: read_tab(s) for s in (A, B)}
NN = {s: len(NAMES[s]) for s in (A, B)}
for s in (A, B):
    say(f"  {s}: {NN[s]:,} nodes")

# The id space must fit the 32-bit halves of the packed key. 170,135 does by four
# orders of magnitude, but an unchecked overflow would silently alias two different
# edges onto one key and inflate conservation.
if max(NN.values()) >= 2 ** 31:
    sys.exit("FATAL: node count exceeds the 32-bit packing this script relies on")

IDX = {s: {g: i for i, g in enumerate(NAMES[s])} for s in (A, B)}


def pairs_path(study):
    p = os.path.join(WORK, f"{study}.pairs")
    if not os.path.exists(p) or os.path.getsize(p) == 0:
        sys.exit(f"FATAL: no edge dump at {p}\n"
                 f"  regenerate:  mcxdump -imx {study}.mci --dump-pairs "
                 f"--dump-upper --no-loops")
    return p


def stream_pairs(study, chunk=CHUNK):
    """Yield (i, j, w) as int32/int32/float32. pandas' C reader is ~3x a numpy
    split-and-cast on this file (5.0 vs 1.7 M rows/s, measured)."""
    for ch in pd.read_csv(pairs_path(study), sep="\t", header=None,
                          names=["i", "j", "w"],
                          dtype={"i": np.int32, "j": np.int32, "w": np.float32},
                          chunksize=chunk, engine="c"):
        yield (ch["i"].to_numpy(), ch["j"].to_numpy(), ch["w"].to_numpy())


# =============================================================================
# 2. orthologs, in index space, restricted to nodes of BOTH networks
# =============================================================================
# A pair only exists for this analysis if both genes are nodes of their own network:
# a gene absent from the graph has no edges to conserve, whatever its orthology.
say("")
say("loading orthologs")


# The SAME two substitutions 06_conservation_join.r:57-58 applies, as regexes
# rather than string surgery, so the two routes cannot diverge on an id shape
# neither of us anticipated. The ortholog-universe assertion below is what proves
# they agree in practice: it reproduces the R path's count exactly.
_RE_P = re.compile(r"\.p[0-9]+$")
_RE_V = re.compile(r"\.v[0-9]+(\.[0-9]+)*$")


def strip_version(g):
    return _RE_V.sub("", _RE_P.sub("", g))


with open(ORTHOGROUPS) as fh:
    rdr = csv.DictReader(fh, delimiter="\t")
    col_a, col_b = SP[A], SP[B]
    for c in (col_a, col_b):
        if c not in rdr.fieldnames:
            sys.exit(f"FATAL: Orthogroups.tsv has no column '{c}'\n"
                     f"  found: {rdr.fieldnames}")
    src_idx, tgt_idx, n_og = [], [], 0
    for row in rdr:
        ga = [strip_version(x.strip()) for x in (row[col_a] or "").split(",") if x.strip()]
        gb = [strip_version(x.strip()) for x in (row[col_b] or "").split(",") if x.strip()]
        ia = [IDX[A][g] for g in ga if g in IDX[A]]
        ib = [IDX[B][g] for g in gb if g in IDX[B]]
        if not ia or not ib:
            continue
        n_og += 1
        for x in ia:
            for y in ib:
                src_idx.append(x)
                tgt_idx.append(y)

src_idx = np.asarray(src_idx, dtype=np.int32)
tgt_idx = np.asarray(tgt_idx, dtype=np.int32)
say(f"  ortholog pairs with BOTH sides in their network: {len(src_idx):,}"
    f"  over {n_og:,} orthogroups")

# THE UNIVERSE MUST MATCH 61's. The node level measured 114,681 pairs over 43,390
# orthogroups; if the edge level disagrees the two are answering different questions
# and no funnel across them means anything.
EXPECT_PAIRS, EXPECT_OG = 114_681, 43_390
if (len(src_idx), n_og) != (EXPECT_PAIRS, EXPECT_OG):
    sys.exit(f"FATAL: ortholog universe is {len(src_idx):,} pairs / {n_og:,} "
             f"orthogroups, but the node level (61) used {EXPECT_PAIRS:,} / "
             f"{EXPECT_OG:,}.\n  The two levels must share a universe.")
say(f"  universe matches the node level (61): {EXPECT_PAIRS:,} pairs / {EXPECT_OG:,} OGs")

# CSR over the SOURCE node domain: ortho of node u is indices[indptr[u]:indptr[u+1]].
order = np.argsort(src_idx, kind="stable")
csr_indices = tgt_idx[order]
counts = np.bincount(src_idx, minlength=NN[A])
csr_indptr = np.zeros(NN[A] + 1, dtype=np.int64)
np.cumsum(counts, out=csr_indptr[1:])
DEG = counts.astype(np.int32)
say(f"  source nodes with an in-network ortholog: {int((DEG > 0).sum()):,}"
    f" of {NN[A]:,}  ({100 * (DEG > 0).mean():.1f}%)")


# =============================================================================
# 3. the target adjacency -- built ONCE, reused by the observed pass and the null
# =============================================================================
say("")
say(f"building the {B} adjacency (sorted int64 keys)")
t = time.time()
n_b_edges = 0
blocks = []
for i, j, _ in stream_pairs(B):
    if i.max() >= NN[B] or j.max() >= NN[B]:
        sys.exit(f"FATAL: {B}.pairs holds an index >= {NN[B]} -- the dump does not "
                 f"match {B}.tab")
    blocks.append(pack(i, j))
    n_b_edges += len(i)
ADJ = np.concatenate(blocks)
del blocks
ADJ.sort()
say(f"  {n_b_edges:,} edges -> {ADJ.nbytes / 2**30:.1f} GiB of keys in "
    f"{time.time() - t:.0f}s")

# THE DUMP MUST BE THE GRAPH. A stale or truncated .pairs would silently shrink the
# target network and depress every conservation rate, which looks exactly like a
# biological result.
EXPECT_EDGES = {"sugarcane": 75_333_769, "purple": 675_955_918}
if n_b_edges != EXPECT_EDGES[B]:
    sys.exit(f"FATAL: {B}.pairs holds {n_b_edges:,} edges, the matrix has "
             f"{EXPECT_EDGES[B]:,} -- regenerate the dump")
# ADJ is already sorted, so duplicates are adjacent: a linear scan costs one bool
# array (0.7 GB) where np.unique would allocate and re-sort another copy of the
# 5 GB key array.
n_dup = int(np.count_nonzero(ADJ[1:] == ADJ[:-1]))
if n_dup:
    say(f"  NOTE: {n_dup:,} duplicate undirected keys collapsed (dump-upper should "
        f"give each edge once)")


def in_adj(keys):
    """Membership in the target edge set, by binary search."""
    pos = np.searchsorted(ADJ, keys)
    np.minimum(pos, len(ADJ) - 1, out=pos)
    return ADJ[pos] == keys


# =============================================================================
# 4. who is nitrogen-responsive -- read from 61, never recomputed
# =============================================================================
def responsive_sets(study):
    """(responsive in `study`, responsive in BOTH species), as boolean masks over
    the node index. Both come from 61_conserved_blocked_nodes.r's status table: its
    rows ARE the responsive genes, and status == conserved_correlated marks the ones
    whose ortholog responds too."""
    f = os.path.join(RESULTS, "conservation", f"{study}_status_blocked_nodes.tsv")
    if not os.path.exists(f):
        sys.exit(f"FATAL: missing {f}\n  run  ./run.sh consblocked 0  first -- it "
                 f"defines which genes are nitrogen-correlated")
    d = pd.read_csv(f, sep="\t", usecols=[f"{study}_gene", "status"])
    resp = np.zeros(NN[study], dtype=bool)
    both = np.zeros(NN[study], dtype=bool)
    miss = 0
    for g, st in zip(d[f"{study}_gene"], d["status"]):
        k = IDX[study].get(g)
        if k is None:
            miss += 1
            continue
        resp[k] = True
        if st == "conserved_correlated":
            both[k] = True
    if miss:
        sys.exit(f"FATAL: {miss} responsive {study} genes are not network nodes -- "
                 f"61 and this script disagree about the node universe")
    return resp, both


RESP_A, BOTH_A = responsive_sets(A)
say("")
say(f"nitrogen-correlated {A} genes (from 61): {int(RESP_A.sum()):,} responsive, "
    f"{int(BOTH_A.sum()):,} also responsive in {B}")

# =============================================================================
# 5. weight deciles, from a sample
# =============================================================================
# Rank-based bins need quantiles of 676M floats. A 5M-edge sample puts the 1%-wide
# quantile boundaries within ~0.1% of the truth, which is far finer than any effect
# here, and it avoids a whole extra pass over the stream.
say("")
say(f"weight deciles from a {N_SAMPLE:,}-edge sample")
sample_w = []
seen = 0
for _, _, w in stream_pairs(A):
    seen += len(w)
    take = rng.random(len(w)) < min(1.0, N_SAMPLE / EXPECT_EDGES[A])
    sample_w.append(w[take])
    if sum(len(x) for x in sample_w) >= N_SAMPLE:
        break
sample_w = np.concatenate(sample_w)
EDGES_BINS = np.quantile(sample_w, np.linspace(0, 1, N_BINS + 1))
EDGES_BINS[0], EDGES_BINS[-1] = -np.inf, np.inf
say("  boundaries: " + " ".join(f"{x:.4g}" for x in EDGES_BINS[1:-1]))
del sample_w


def weight_bin(w):
    return np.clip(np.searchsorted(EDGES_BINS, w, side="right") - 1, 0, N_BINS - 1)


# =============================================================================
# 6. the join
# =============================================================================
def conserved_mask(u, v, indices):
    """For each edge (u[e], v[e]), does ANY ortholog pair land on a target edge?

    Returns (mask, n_hits). Only edges whose BOTH endpoints have an in-network
    ortholog can possibly be conserved, so those are the only ones expanded --
    that alone drops ~63% of sugarcane's edges and ~81% of purple's.
    """
    m = len(u)
    out = np.zeros(m, dtype=bool)
    hits = np.zeros(m, dtype=np.int32)
    du, dv = DEG[u], DEG[v]
    live = np.flatnonzero((du > 0) & (dv > 0))
    if len(live) == 0:
        return out, hits
    lu, lv = u[live], v[live]
    cu, cv = DEG[lu].astype(np.int64), DEG[lv].astype(np.int64)
    ncand = cu * cv
    total = int(ncand.sum())
    if total == 0:
        return out, hits

    # candidate -> which edge, and its offset within that edge's cartesian product
    eid = np.repeat(np.arange(len(live), dtype=np.int64), ncand)
    starts = np.cumsum(ncand) - ncand
    off = np.arange(total, dtype=np.int64) - np.repeat(starts, ncand)
    cv_rep = np.repeat(cv, ncand)
    o1 = indices[np.repeat(csr_indptr[lu], ncand) + off // cv_rep]
    o2 = indices[np.repeat(csr_indptr[lv], ncand) + off % cv_rep]

    hit = in_adj(pack(o1, o2))
    if hit.any():
        per_edge = np.bincount(eid[hit], minlength=len(live))
        out[live] = per_edge > 0
        hits[live] = per_edge
    return out, hits


# =============================================================================
# 7. observed pass -- one stream, every stratum and every decile
# =============================================================================
STRATA = ("all", "resp_source", "resp_both")
tot = {s: np.zeros(N_BINS, dtype=np.int64) for s in STRATA}
con = {s: np.zeros(N_BINS, dtype=np.int64) for s in STRATA}
conserved_genes = np.zeros(NN[A], dtype=bool)
# THE RESPONSIVE STRATA ARE KEPT IN FULL, not sampled. See the null below: a
# Bernoulli sample sized for a 676M-edge network catches ~46 of purple's 5,692
# resp_both edges, and if none of those happens to be conserved the stratum's null
# divides by nothing and reports fold 0.0 against a true rate of 5.9%. These are
# small enough to carry whole -- 408k and 850k edges.
rs_i, rs_j, rs_w = [], [], []

out_edges = os.path.join(OUT_DIR, f"conserved_edges_{A}_to_{B}_pearson.tsv")
say("")
say(f"streaming {A} ({EXPECT_EDGES[A]:,} edges) -> {os.path.basename(out_edges)}")
t = time.time()
n_seen = 0
n_src_edges = 0
with open(out_edges, "w") as fh:
    fh.write("gene1\tgene2\tweight\tn_ortho_hits\n")
    for i, j, w in stream_pairs(A):
        n_src_edges += len(i)
        if i.max() >= NN[A] or j.max() >= NN[A]:
            sys.exit(f"FATAL: {A}.pairs holds an index >= {NN[A]}")
        mask, hits = conserved_mask(i, j, csr_indices)
        wb = weight_bin(w)

        sel = {"all": np.ones(len(i), dtype=bool),
               "resp_source": RESP_A[i] & RESP_A[j],
               "resp_both": BOTH_A[i] & BOTH_A[j]}
        for s in STRATA:
            tot[s] += np.bincount(wb[sel[s]], minlength=N_BINS)
            con[s] += np.bincount(wb[sel[s] & mask], minlength=N_BINS)

        rs = np.flatnonzero(sel["resp_source"])
        if len(rs):
            rs_i.append(i[rs]); rs_j.append(j[rs]); rs_w.append(w[rs])

        k = np.flatnonzero(mask)
        if len(k):
            conserved_genes[i[k]] = True
            conserved_genes[j[k]] = True
            block = "".join(
                f"{NAMES[A][a]}\t{NAMES[A][b]}\t{ww:.6g}\t{hh}\n"
                for a, b, ww, hh in zip(i[k], j[k], w[k], hits[k]))
            fh.write(block)

        n_seen += len(i)
        if n_seen % (CHUNK * 10) < CHUNK:
            say(f"  {n_seen:,} edges ({100 * n_seen / EXPECT_EDGES[A]:.1f}%) | "
                f"conserved {int(con['all'].sum()):,} "
                f"({100 * con['all'].sum() / max(n_seen, 1):.2f}%) | "
                f"{time.time() - t:.0f}s")

if n_src_edges != EXPECT_EDGES[A]:
    sys.exit(f"FATAL: streamed {n_src_edges:,} {A} edges, expected "
             f"{EXPECT_EDGES[A]:,}")
say(f"  done in {time.time() - t:.0f}s")


out_genes = os.path.join(OUT_DIR, f"conserved_genes_{A}_pearson.txt")
with open(out_genes, "w") as fh:
    for g in NAMES[A][conserved_genes]:
        fh.write(g + "\n")
say(f"wrote {os.path.basename(out_genes)}  ({int(conserved_genes.sum()):,} genes)")

# =============================================================================
# 8. the null -- permute the target column of the ortholog table
# =============================================================================
# 13_conservation_null.r's construction exactly. Permuting `csr_indices` while
# keeping `csr_indptr` IS "permute the target column": it preserves each source
# gene's fan-out, which source genes have any ortholog at all, and the multiset of
# target genes used -- so the only thing destroyed is WHICH target each source maps
# to. A rate that survives this is not an artifact of density or fan-out.
#
# Run on a Bernoulli sample of the source: at 5M edges the standard error on a 10%
# rate is 0.013%, far finer than the effect, and it turns a full pass into seconds so
# replicates are affordable.
say("")
say(f"null: {N_REPS} replicates")
frac = min(1.0, N_SAMPLE / EXPECT_EDGES[A])
si, sj, sw = [], [], []
for i, j, w in stream_pairs(A):
    take = rng.random(len(i)) < frac
    si.append(i[take]); sj.append(j[take]); sw.append(w[take])
si = np.concatenate(si); sj = np.concatenate(sj); sw = np.concatenate(sw)

# EACH STRATUM GETS THE EDGE SET IT CAN AFFORD.
#   all          a Bernoulli sample -- at 5M edges the SE on a 1.5% rate is 0.005%,
#                far finer than the effect, and a full pass per replicate is not
#                affordable at 676M edges.
#   resp_source  every edge, because there are only a few hundred thousand
#   resp_both    every edge, because there are a few thousand and sampling them
#                is what produced a fold of 0.0 against a real rate of 5.9%
# The basis is written into the output so a rate is never read without knowing
# which set it came from.
RS_I = np.concatenate(rs_i) if rs_i else np.zeros(0, dtype=np.int32)
RS_J = np.concatenate(rs_j) if rs_j else np.zeros(0, dtype=np.int32)
RS_W = np.concatenate(rs_w) if rs_w else np.zeros(0, dtype=np.float32)
del rs_i, rs_j, rs_w
rb = BOTH_A[RS_I] & BOTH_A[RS_J]

NULL_SET = {
    "all": (si, sj, weight_bin(sw), "sample"),
    "resp_source": (RS_I, RS_J, weight_bin(RS_W), "complete stratum"),
    "resp_both": (RS_I[rb], RS_J[rb], weight_bin(RS_W[rb]), "complete stratum"),
}
for s_ in STRATA:
    a_, _, _, basis = NULL_SET[s_]
    say(f"  {s_:<12} {len(a_):>10,} edges  ({basis})")


def rates(indices):
    r = {}
    for st in STRATA:
        a_, b_, wb_, _ = NULL_SET[st]
        if len(a_) == 0:
            r[st] = (0, 0)
            r[f"{st}__bins"] = (np.zeros(N_BINS, dtype=np.int64),
                                np.zeros(N_BINS, dtype=np.int64))
            continue
        mask, _ = conserved_mask(a_, b_, indices)
        r[st] = (len(a_), int(mask.sum()))
        r[f"{st}__bins"] = (np.bincount(wb_, minlength=N_BINS),
                            np.bincount(wb_[mask], minlength=N_BINS))
    return r


obs_s = rates(csr_indices)
say("  observed on those sets: " +
    ", ".join(f"{s_} {100 * obs_s[s_][1] / max(obs_s[s_][0], 1):.3f}%" for s_ in STRATA))

null_rate = {s_: [] for s_ in STRATA}
null_bin = {s_: [] for s_ in STRATA}
t = time.time()
for r in range(1, N_REPS + 1):
    perm = rng.permutation(csr_indices)
    nr = rates(perm)
    for s_ in STRATA:
        n, k = nr[s_]
        null_rate[s_].append(k / n if n else np.nan)
        bt, bk = nr[f"{s_}__bins"]
        null_bin[s_].append(np.where(bt > 0, bk / np.maximum(bt, 1), np.nan))
    if r % 5 == 0 or r == N_REPS:
        say(f"  rep {r}/{N_REPS}  ({(time.time() - t) / 60:.1f} min)")

# =============================================================================
# 9. outputs
# =============================================================================
rows = []
for s in STRATA:
    n_tot, n_con = int(tot[s].sum()), int(con[s].sum())
    nr = np.asarray(null_rate[s], dtype=float)
    obs_rate = obs_s[s][1] / obs_s[s][0] if obs_s[s][0] else np.nan
    rows.append(dict(
        direction=f"{A}_to_{B}", stratum=s, weight_bin="ALL",
        weight_lo="", weight_hi="",
        edges=n_tot, conserved=n_con,
        conserved_fraction=round(n_con / n_tot, 8) if n_tot else np.nan,
        # the null is measured on `null_basis`, so the rate ON THAT SET is carried
        # beside it -- comparing a full-pass rate against a null computed on a
        # different set would confound the two
        null_basis=NULL_SET[s][3],
        null_set_edges=obs_s[s][0], null_set_conserved=obs_s[s][1],
        null_set_rate=round(obs_rate, 8) if obs_s[s][0] else np.nan,
        null_mean=round(float(np.nanmean(nr)), 8) if len(nr) else np.nan,
        null_sd=round(float(np.nanstd(nr, ddof=1)), 8) if len(nr) > 1 else np.nan,
        fold_over_null=round(obs_rate / float(np.nanmean(nr)), 4)
        if len(nr) and np.nanmean(nr) > 0 else np.nan,
        p_emp=round((int((nr >= obs_rate).sum()) + 1) / (len(nr) + 1), 6)))
    for b in range(N_BINS):
        nb = np.asarray([x[b] for x in null_bin[s]], dtype=float)
        sb_t, sb_k = obs_s[f"{s}__bins"]
        sr = sb_k[b] / sb_t[b] if sb_t[b] else np.nan
        rows.append(dict(
            direction=f"{A}_to_{B}", stratum=s, weight_bin=f"D{b + 1}",
            weight_lo=round(float(EDGES_BINS[b]), 6) if b else "min",
            weight_hi=round(float(EDGES_BINS[b + 1]), 6) if b < N_BINS - 1 else "max",
            edges=int(tot[s][b]), conserved=int(con[s][b]),
            conserved_fraction=round(con[s][b] / tot[s][b], 8) if tot[s][b] else np.nan,
            null_basis=NULL_SET[s][3],
            null_set_edges=int(sb_t[b]), null_set_conserved=int(sb_k[b]),
            null_set_rate=round(float(sr), 8) if sb_t[b] else np.nan,
            null_mean=round(float(np.nanmean(nb)), 8) if len(nb) else np.nan,
            null_sd=round(float(np.nanstd(nb, ddof=1)), 8) if len(nb) > 1 else np.nan,
            fold_over_null=round(sr / float(np.nanmean(nb)), 4)
            if sb_t[b] and np.nanmean(nb) > 0 else np.nan,
            p_emp=round((int((nb >= sr).sum()) + 1) / (len(nb) + 1), 6)
            if sb_t[b] else np.nan))

summ = pd.DataFrame(rows)
f_summ = os.path.join(OUT_DIR, f"conservation_summary_{A}_to_{B}_pearson.tsv")
summ.to_csv(f_summ, sep="\t", index=False)
say(f"wrote {os.path.basename(f_summ)}  ({len(summ)} rows)")

print()
print("=" * 78)
print(f"EDGE-LEVEL CONSERVATION  {A} -> {B}   (Pearson-only)")
print("=" * 78)
hdr = summ[summ.weight_bin == "ALL"]
for _, r in hdr.iterrows():
    # The basis is printed because a fold computed on a 5M sample and a fold
    # computed on the whole stratum are not the same statistic.
    print(f"  {r.stratum:<12} {r.edges:>12,} edges  {r.conserved:>11,} conserved  "
          f"{100 * r.conserved_fraction:>7.3f}%   fold {r.fold_over_null:>6}  "
          f"p {r.p_emp:<9} [null on {r.null_basis}, {r.null_set_edges:,} edges]")
print()
print(f"  by edge strength ({A} weight deciles, stratum = all):")
for _, r in summ[(summ.stratum == "all") & (summ.weight_bin != "ALL")].iterrows():
    print(f"    {r.weight_bin:<4} w {str(r.weight_lo):>9}-{str(r.weight_hi):<9} "
          f"{r.edges:>12,} edges  {100 * r.conserved_fraction:>7.3f}%   "
          f"fold {r.fold_over_null}")
print()
print("  Each direction is read ONLY against its own null: purple's density is 3.2x")
print("  sugarcane's, so a random ortholog pair is far likelier to be an edge there")
print("  and the two rates are not comparable to each other.")
n_both = int(BOTH_A.sum())
print(f"\n  The resp_both stratum draws on {n_both:,} {A} genes. An edge needs TWO of")
print("  them, so a small count there is arithmetic before it is biology.")
say("")
say(f"done: {DIRECTION}  ({(time.time() - T0) / 60:.1f} min total)")
