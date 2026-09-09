#!/usr/bin/env python3
# =============================================================================
# 50_fastgreedy_scout.py -- the hierarchical family, bounded, sugarcane only.
#
# WHY CLASSICAL HIERARCHICAL CLUSTERING IS NOT HERE. It needs a dense pairwise
# dissimilarity matrix: 170,135^2 doubles is 232 GB for purple, ~1.4e10 pairs,
# with linkage on top at O(n^2 log n). Worse, it is not a graph method at all --
# it would have to be built from the correlation matrix this pipeline
# deliberately thresholded away at |r| >= 0.8, so it would be answering a
# different question from every other row in the comparison.
#
# WHAT THIS IS INSTEAD. fast-greedy (Clauset-Newman-Moore) is the graph-native
# member of the agglomerative family: it merges communities greedily by
# modularity gain and yields a real dendrogram, so "cut it at a different level"
# means something. It is the honest stand-in, and running it answers the
# question the plan actually needs answered -- what does the hierarchical family
# cost and give here -- without pretending the classical method was affordable.
#
# TWO THINGS TO EXPECT, both recorded rather than worked around:
#   * CNM optimises modularity, so it inherits modularity's resolution limit.
#     At 2m = 1.5e8 (sugarcane) it cannot resolve communities below ~sqrt(2m)
#     edges, and it is known to produce very lopsided communities.
#   * cost. O(m d log n) with a large constant, on 75.3M edges.
#
# THE CAP. The C call cannot be interrupted by a Python signal, so the work runs
# in a forked child and the parent kills it at CLUSTER_FASTGREEDY_CAP_S. Hitting
# the cap IS a result -- it gets written to the manifest as such, and purple is
# not attempted.
#
# RUN: through run.sh  ->  ./run.sh fastgreedy sugarcane
# =============================================================================

import multiprocessing as mp
import os
import sys
import time

import igraph


def say(*a):
    print("[%s]" % time.strftime("%H:%M:%S"), *a, flush=True)


STUDY = os.environ["CLEAN_STUDY"]
WORK = os.environ["CLEAN_WORK_DIR"]
CAP = int(os.environ.get("CLUSTER_FASTGREEDY_CAP_S", "14400"))
N_CUTS = int(os.environ.get("CLUSTER_FASTGREEDY_CUTS", "3"))

TAB = os.path.join(WORK, f"{STUDY}.tab")
PAIRS = os.path.join(WORK, f"{STUDY}.pairs")
SW = os.path.join(WORK, f"sweep_{STUDY}")
MANIFEST = os.path.join(SW, "partitions.tsv")
STATUS = os.path.join(SW, f"fastgreedy_{STUDY}.status")
N_NODES = sum(1 for _ in open(TAB))

if not os.path.exists(PAIRS):
    sys.exit(f"FATAL: no {PAIRS} -- run leidensweep first, it writes the dump")

say(f"{STUDY}: loading (same dump Leiden uses)")
G = igraph.Graph.Read_Ncol(PAIRS, names=True, weights=True, directed=False)
to_mci = [int(n) for n in G.vs["name"]]
assert sorted(to_mci) == list(range(N_NODES)), "index map is not a permutation"
say(f"  {G.vcount():,} vertices, {G.ecount():,} edges")


def write_cls(membership, path):
    groups = {}
    for v, c in enumerate(membership):
        groups.setdefault(c, []).append(to_mci[v])
    order = sorted(groups.values(), key=len, reverse=True)
    if sum(len(g) for g in order) != N_NODES:
        sys.exit("FATAL: partition does not cover every node")
    tmp = path + ".part"
    with open(tmp, "w") as fh:
        fh.write("(mclheader\nmcltype matrix\ndimensions %dx%d\n)\n(mclmatrix\nbegin\n"
                 % (N_NODES, len(order)))
        for ci, members in enumerate(order):
            members.sort()
            fh.write("%d %s $\n" % (ci, " ".join(map(str, members))))
        fh.write(")\n")
    os.replace(tmp, path)
    return len(order), len(order[0])


def worker(q):
    t0 = time.time()
    dend = G.community_fastgreedy(weights="weight")
    secs = time.time() - t0
    opt = dend.optimal_count
    # The dendrogram is the whole reason to run this: report the modularity-optimal
    # cut and two coarser ones, so "cut it elsewhere" is an option on the table
    # rather than a hypothetical.
    cuts = sorted({opt, max(2, opt // 4), max(2, opt // 16)}, reverse=True)[:N_CUTS]
    out = []
    for n in cuts:
        memb = dend.as_clustering(n).membership
        out.append((n, memb))
    q.put((secs, opt, out))


if __name__ == "__main__":
    say(f"fast-greedy (CNM), hard cap {CAP}s ({CAP / 3600:.1f} h)")
    q = mp.Queue()
    p = mp.Process(target=worker, args=(q,))
    t0 = time.time()
    p.start()
    p.join(CAP)

    if p.is_alive():
        p.terminate()
        p.join()
        msg = (f"CAP HIT: fast-greedy did not finish on {STUDY} "
               f"({G.ecount():,} edges) within {CAP}s")
        say(msg)
        say("  This is the result: the hierarchical family is not affordable here.")
        say("  purple (9x the edges) is not attempted.")
        with open(STATUS, "w") as fh:
            fh.write(msg + "\n")
        sys.exit(0)

    secs, opt, out = q.get()
    say(f"  finished in {secs:.0f}s; modularity-optimal cut = {opt:,} communities")
    if not os.path.exists(MANIFEST):
        with open(MANIFEST, "w") as fh:
            fh.write("method\tparam\tresource\tfile\tn_clusters\tlargest\truntime_s\n")
    for n, memb in out:
        path = os.path.join(SW, f"cls.fastgreedy.{n}")
        ncl, largest = write_cls(memb, path)
        with open(MANIFEST, "a") as fh:
            fh.write(f"fastgreedy\tcut{n}\t\t{path}\t{ncl}\t{largest}\t{secs:.1f}\n")
        say(f"  cut at {n:,}: {ncl:,} clusters, largest {largest:,} "
            f"({100 * largest / N_NODES:.2f}%)")
    with open(STATUS, "w") as fh:
        fh.write(f"OK: finished in {secs:.0f}s, optimal cut {opt}\n")
    say(f"done: {STUDY}")
