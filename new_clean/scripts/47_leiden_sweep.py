#!/usr/bin/env python3
# =============================================================================
# 47_leiden_sweep.py -- Leiden and Louvain on the unpruned graph, in MCL's own
# cluster format so clm info can score them beside MCL's partitions.
#
# WHY A SECOND METHOD. -I 2 was never compared against anything. Sweeping MCL's
# inflation says which inflation is best; it cannot say whether MCL is the right
# tool. Leiden is the fast graph-native alternative, and -- unlike MCL -- it
# prunes nothing: MCL's -scheme 7 keeps 1200 neighbours per node while computing,
# against a mean degree of 1,477 (sugarcane) and 7,946 (purple).
#
# WHY CPM RATHER THAN MODULARITY IS THE MAIN LADDER. Newman modularity has a
# resolution limit: it cannot resolve communities smaller than ~sqrt(2m) edges,
# and 2m here is 1.35e9 (purple). CPM has no such limit and its resolution gamma
# reads directly -- a community is kept while its internal weighted density
# exceeds gamma. Modularity is run once anyway, as the comparison.
#
# ---------------------------------------------------------------------------
# THE TRAP THIS SCRIPT EXISTS TO AVOID, MEASURED NOT ASSUMED
#
# igraph's Read_Ncol(names=False) does NOT treat the integers in the file as
# vertex ids. It assigns ids by ORDER OF FIRST APPEARANCE. Verified:
#
#     file: 5 3 0.5 / 3 1 0.7 / 1 0 0.9
#     Read_Ncol(names=False) -> vcount 4, edges [(0,1),(1,2),(2,3)]
#     Read_Ncol(names=True)  -> names ['5','3','1','0'], same edges
#
# So node 5 becomes vertex 0. Writing that membership out against the .mci node
# domain would assign every gene to the wrong cluster -- and clm info would
# score it happily, returning plausible numbers for a scrambled partition.
# Hence names=True and an explicit inverse permutation, asserted to be a
# permutation of 0..N-1 before anything is written.
# ---------------------------------------------------------------------------
#
# CLUSTER FILE FORMAT, verified against a real mcl output (mcxio(5)):
#
#     (mclheader
#     mcltype matrix
#     dimensions <n_nodes>x<n_clusters>
#     )
#     (mclmatrix
#     begin
#     0 0 1 2 $
#     1 3 4 5 $
#     )
#
# Columns are clusters, entries are node indices, row domain canonical.
#
# RUN: through run.sh  ->  ./run.sh leidensweep sugarcane [scout]
# =============================================================================

import os
import subprocess
import sys
import time

import igraph


def env_req(name):
    v = os.environ.get(name)
    if not v:
        sys.exit(f"FATAL: {name} is not set")
    return v


def say(*a):
    print("[%s]" % time.strftime("%H:%M:%S"), *a, flush=True)


STUDY = env_req("CLEAN_STUDY")
WORK = env_req("CLEAN_WORK_DIR")
BIN = env_req("CLEAN_MCL_BIN_DIR")
OUT_DIR = env_req("CLEAN_OUT_DIR")
GAMMA = os.environ.get("CLUSTER_LEIDEN_GAMMA", "").split()
SCOUT = os.environ.get("CLUSTER_LEIDEN_SCOUT", "0.01 0.05 0.2 0.5").split()
N_ITER = int(os.environ.get("CLUSTER_LEIDEN_ITER", "2"))
MODE = os.environ.get("CLEAN_LEIDEN_MODE", "full")  # full | scout

MCI = os.path.join(WORK, f"{STUDY}.mci")
TAB = os.path.join(WORK, f"{STUDY}.tab")
PAIRS = os.path.join(WORK, f"{STUDY}.pairs")
SW = os.path.join(WORK, f"sweep_{STUDY}")
MANIFEST = os.path.join(SW, "partitions.tsv")
os.makedirs(SW, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)

for f in (MCI, TAB):
    if not os.path.exists(f):
        sys.exit(f"FATAL: missing {f} -- run ./run.sh pearsonmci {STUDY}")

N_NODES = sum(1 for _ in open(TAB))
say(f"{STUDY}: {N_NODES:,} nodes in the tab file")


# --- the edge list, dumped once and reused -----------------------------------
# --dump-upper gives each undirected edge once; values are emitted by default.
if not os.path.exists(PAIRS) or os.path.getsize(PAIRS) == 0:
    say(f"dumping edges -> {os.path.basename(PAIRS)}")
    t0 = time.time()
    with open(PAIRS, "wb") as fh:
        subprocess.run(
            [os.path.join(BIN, "mcxdump"), "-imx", MCI,
             "--dump-pairs", "--dump-upper", "--no-loops"],
            stdout=fh, stderr=subprocess.DEVNULL, check=True)
    say(f"  dumped in {time.time() - t0:.0f}s, {os.path.getsize(PAIRS) / 2**30:.1f} GiB")
else:
    say(f"reusing {os.path.basename(PAIRS)} ({os.path.getsize(PAIRS) / 2**30:.1f} GiB)")


# --- load, and undo igraph's re-indexing -------------------------------------
say("loading into igraph (C-level reader; edges never become Python objects)")
t0 = time.time()
G = igraph.Graph.Read_Ncol(PAIRS, names=True, weights=True, directed=False)
say(f"  loaded in {time.time() - t0:.0f}s: {G.vcount():,} vertices, {G.ecount():,} edges")

if G.vcount() != N_NODES:
    sys.exit(f"FATAL: igraph sees {G.vcount()} vertices, the matrix has {N_NODES}.\n"
             "  Every node in a thresholded correlation network has degree >= 1,\n"
             "  so a shortfall means the dump is incomplete.")

# vertex i in igraph  ->  node to_mci[i] in the .mci domain
to_mci = [int(n) for n in G.vs["name"]]
if sorted(to_mci) != list(range(N_NODES)):
    sys.exit("FATAL: the igraph->mci index map is not a permutation of 0..N-1. "
             "Refusing to write a partition that would be scrambled.")
say("  index map verified as a permutation of 0..N-1")

W = G.es["weight"]
say(f"  weights: min {min(W):.4g}  max {max(W):.4g}  mean {sum(W) / len(W):.4g}")


# --- membership -> mcl cluster file ------------------------------------------
def write_cls(membership, path):
    """Emit an mclmatrix cluster file over the .mci node domain.

    membership is indexed by IGRAPH vertex id; to_mci maps it back. Clusters are
    renumbered densely from 0 and emitted in descending size order, which is what
    mcl itself does and makes the head of the file readable.
    """
    groups = {}
    for v, c in enumerate(membership):
        groups.setdefault(c, []).append(to_mci[v])
    order = sorted(groups.values(), key=len, reverse=True)

    n_assigned = sum(len(g) for g in order)
    if n_assigned != N_NODES:
        sys.exit(f"FATAL: partition covers {n_assigned} nodes, expected {N_NODES}")

    tmp = path + ".part"
    with open(tmp, "w") as fh:
        fh.write("(mclheader\nmcltype matrix\ndimensions %dx%d\n)\n"
                 "(mclmatrix\nbegin\n" % (N_NODES, len(order)))
        for ci, members in enumerate(order):
            members.sort()
            fh.write("%d %s $\n" % (ci, " ".join(map(str, members))))
        fh.write(")\n")
    os.replace(tmp, path)
    return len(order), len(order[0])


def done_already(method, param):
    """Skip cells already on disk so a re-run resumes instead of restarting.
    The scout writes real cells; the full ladder must not recompute them."""
    path = os.path.join(SW, f"cls.{method}.{param}".replace(" ", ""))
    if os.path.exists(path) and os.path.getsize(path) > 0:
        say(f"  {method} {param}: already present, skipped")
        return True
    return False


def record(method, param, membership, secs):
    tag = f"cls.{method}.{param}".replace(" ", "")
    path = os.path.join(SW, tag)
    ncl, largest = write_cls(membership, path)
    with open(MANIFEST, "a") as fh:
        fh.write(f"{method}\t{param}\t\t{path}\t{ncl}\t{largest}\t{secs:.1f}\n")
    say(f"  {method} {param}: {ncl:,} clusters, largest {largest:,} "
        f"({100 * largest / N_NODES:.2f}%), {secs:.0f}s")


if not os.path.exists(MANIFEST):
    with open(MANIFEST, "w") as fh:
        fh.write("method\tparam\tresource\tfile\tn_clusters\tlargest\truntime_s\n")


# --- the runs ----------------------------------------------------------------
gammas = SCOUT if MODE == "scout" else (GAMMA or SCOUT)
if MODE == "scout":
    say(f"SCOUT mode: gamma = {' '.join(gammas)} (calibrating the ladder, "
        "gamma is compared against edge weights so its useful range is "
        "a property of this network)")

for g in gammas:
    if done_already("leiden_cpm", g):
        continue
    say(f"Leiden CPM gamma={g}")
    t0 = time.time()
    part = G.community_leiden(objective_function="CPM", weights="weight",
                              resolution=float(g), n_iterations=N_ITER)
    record("leiden_cpm", g, part.membership, time.time() - t0)

if MODE != "scout":
    if not done_already("leiden_mod", "1.0"):
        say("Leiden modularity (resolution 1.0)")
        t0 = time.time()
        part = G.community_leiden(objective_function="modularity", weights="weight",
                                  resolution=1.0, n_iterations=N_ITER)
        record("leiden_mod", "1.0", part.membership, time.time() - t0)

    if not done_already("louvain", "1.0"):
        say("Louvain (multilevel)")
        t0 = time.time()
        part = G.community_multilevel(weights="weight")
        record("louvain", "1.0", part.membership, time.time() - t0)

say(f"wrote {MANIFEST}")
say(f"done: {STUDY}")
