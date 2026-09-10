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

    # A DISCONNECTED graph cannot be cut below its component count. igraph builds
    # (n_vertices - n_components) merges, and as_clustering(k) takes
    # (n_vertices - k) steps, so any k below n_components asks for more steps than
    # the merges matrix has and igraph raises
    #   "Number of steps is greater than number of rows in merges matrix"
    # Measured here: 101,990 vertices but only 101,032 merges, i.e. 958 components,
    # and optimal_count asked for 101,652 steps. Derive the floor from the merges
    # matrix itself rather than trusting optimal_count.
    floor = G.vcount() - len(dend.merges)

    # optimal_count is typically BELOW that floor here (measured: 338 against a
    # floor of 958), so the modularity-optimal cut is not reachable at all and
    # every "coarser" cut collapses onto the floor -- which is just the connected
    # components. To get something comparable, cut instead at granularities the
    # other methods actually reach, so the row means the same thing as an MCL or
    # Leiden row at similar cluster count.
    cuts = sorted({floor, G.vcount() // 4, G.vcount() // 2}, reverse=True)[:N_CUTS]
    out = []
    for n in cuts:
        memb = dend.as_clustering(n).membership
        out.append((n, memb))
    q.put((secs, opt, floor, out))


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

    secs, opt, floor, out = q.get()
    say(f"  finished in {secs:.0f}s; modularity-optimal cut = {opt:,} communities")
    say(f"  graph has {floor:,} connected components -- no cut below that is possible")
    if opt < floor:
        say(f"  NOTE: modularity-optimal cut ({opt:,}) is BELOW the component floor,")
        say("        so CNM's own preferred partition is unreachable on this graph.")
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
