#!/usr/bin/env python3
# =============================================================================
# verify_cluster_roundtrip.py -- the gate that licenses comparing Leiden to MCL.
#
# THE RISK. Leiden partitions are scored by `clm info` against the .mci, exactly
# as MCL's are. That only means anything if the cluster file 47 writes places
# each gene in the cluster Leiden actually put it in. igraph does NOT preserve
# the .mci node ids -- Read_Ncol assigns vertex ids by order of first appearance
# (measured: a file starting "5 3" makes node 5 into vertex 0) -- so 47 carries
# an explicit inverse permutation. If that permutation were wrong, or dropped,
# clm info would still return perfectly plausible eff/mf/af/mod for a scrambled
# partition and nothing else in the pipeline would notice. This is the one place
# a silent error invalidates the entire benchmark.
#
# THE TEST. Take a real MCL partition of the real graph. Push it through the
# same path a Leiden partition takes -- into igraph vertex space via the same
# permutation, then back out through 47's writer -- and require clm info to
# return the SAME numbers for the re-emitted file as for MCL's own.
#
# It is a genuine test of the permutation because the permutation is not the
# identity: mcxdump emits column by column, so vertex order follows node 0's
# neighbour list, not 0,1,2,...  The script asserts that before trusting a pass.
#
# RUN: CLEAN_STUDY=sugarcane CLEAN_WORK_DIR=... CLEAN_MCL_BIN_DIR=... \
#        python verify_cluster_roundtrip.py <a cls file>
# =============================================================================

import os
import subprocess
import sys
import time

import igraph


def say(*a):
    print("[%s]" % time.strftime("%H:%M:%S"), *a, flush=True)


STUDY = os.environ["CLEAN_STUDY"]
WORK = os.environ["CLEAN_WORK_DIR"]
BIN = os.environ["CLEAN_MCL_BIN_DIR"]
SRC = sys.argv[1]

MCI = os.path.join(WORK, f"{STUDY}.mci")
TAB = os.path.join(WORK, f"{STUDY}.tab")
PAIRS = os.path.join(WORK, f"{STUDY}.pairs")
N_NODES = sum(1 for _ in open(TAB))


def parse_cls(path):
    """mcl cluster file -> {node index: cluster index}. Records run from the
    cluster index to a '$' and may wrap across lines, so tokens are streamed."""
    memb = {}
    in_matrix = started = False
    cur = None
    head = True
    for line in open(path):
        if line.startswith("(mclmatrix"):
            in_matrix = True
            continue
        if not in_matrix:
            continue
        if line.startswith("begin"):
            started = True
            continue
        if not started:
            continue
        if line.startswith(")"):
            break
        for tok in line.split():
            if tok == "$":
                head = True
                cur = None
                continue
            if head:
                cur = int(tok)
                head = False
            else:
                memb[int(tok)] = cur
    return memb


def clm_info(cls_path):
    out = subprocess.run([os.path.join(BIN, "clm"), "info", MCI, cls_path],
                         capture_output=True, text=True, check=True).stdout
    return {k: v for k, v in
            (f.split("=", 1) for f in out.split() if "=" in f and not f.startswith("src"))}


# --- the same load 47 does ---------------------------------------------------
say("loading the graph exactly as 47 does")
G = igraph.Graph.Read_Ncol(PAIRS, names=True, weights=True, directed=False)
to_mci = [int(n) for n in G.vs["name"]]
assert sorted(to_mci) == list(range(N_NODES)), "index map is not a permutation"

identity = to_mci == list(range(N_NODES))
say(f"  {G.vcount():,} vertices; igraph->mci map is "
    f"{'THE IDENTITY -- this test proves nothing, use a real dump' if identity else 'a non-trivial permutation (good: the test has teeth)'}")
say(f"  first 10 of the map: {to_mci[:10]}")
if identity:
    sys.exit("FATAL: identity permutation, the round-trip cannot detect the bug it exists for")

# --- MCL's partition -> igraph vertex space -> 47's writer -------------------
say(f"reading {os.path.basename(SRC)}")
mcl_memb = parse_cls(SRC)
if len(mcl_memb) != N_NODES:
    sys.exit(f"FATAL: source partition covers {len(mcl_memb)} of {N_NODES} nodes")
say(f"  {len(set(mcl_memb.values())):,} clusters over {len(mcl_memb):,} nodes")

# This is the step a Leiden result takes: a membership indexed by IGRAPH vertex.
memb_igraph = [mcl_memb[to_mci[v]] for v in range(G.vcount())]

groups = {}
for v, c in enumerate(memb_igraph):
    groups.setdefault(c, []).append(to_mci[v])
order = sorted(groups.values(), key=len, reverse=True)
assert sum(len(g) for g in order) == N_NODES

RT = os.path.join(WORK, f"roundtrip.{STUDY}.cls")
with open(RT, "w") as fh:
    fh.write("(mclheader\nmcltype matrix\ndimensions %dx%d\n)\n(mclmatrix\nbegin\n"
             % (N_NODES, len(order)))
    for ci, members in enumerate(order):
        members.sort()
        fh.write("%d %s $\n" % (ci, " ".join(map(str, members))))
    fh.write(")\n")
say(f"  re-emitted through 47's writer -> {os.path.basename(RT)}")

# --- NEGATIVE CONTROL --------------------------------------------------------
# A passing round-trip is only evidence if the test could have failed. Here the
# permutation is deliberately dropped -- membership written against igraph's own
# vertex ids, as it would be if someone "simplified" 47 by using names=False --
# and clm info must return DIFFERENT numbers. Only 4.7% of sugarcane's nodes sit
# at a non-identity position, so this is not obvious a priori and is measured
# rather than argued.
groups_bad = {}
for v, c in enumerate(memb_igraph):
    groups_bad.setdefault(c, []).append(v)          # v, not to_mci[v] -- the bug
order_bad = sorted(groups_bad.values(), key=len, reverse=True)
BAD = os.path.join(WORK, f"roundtrip.{STUDY}.scrambled.cls")
with open(BAD, "w") as fh:
    fh.write("(mclheader\nmcltype matrix\ndimensions %dx%d\n)\n(mclmatrix\nbegin\n"
             % (N_NODES, len(order_bad)))
    for ci, members in enumerate(order_bad):
        members.sort()
        fh.write("%d %s $\n" % (ci, " ".join(map(str, members))))
    fh.write(")\n")

# --- the verdict -------------------------------------------------------------
a, b = clm_info(SRC), clm_info(RT)
c = clm_info(BAD)
say("")
say("%-6s %14s %14s %14s   %s" % ("", "mcl's own", "round-tripped", "no-permutation", ""))
ok = True
detected = False
for k in ("eff", "mod", "mf", "af", "ncl", "max", "sgl"):
    same = a.get(k) == b.get(k)
    ok &= same
    if a.get(k) != c.get(k):
        detected = True
    say("%-6s %14s %14s %14s   %s" % (k, a.get(k), b.get(k), c.get(k),
                                      "same" if same else "DIFFERENT"))
say("")
if not detected:
    say("FAIL (negative control) -- dropping the permutation changed NOTHING, so")
    say("     this test cannot detect the bug it exists for. A pass means little.")
    ok = False
else:
    say("negative control OK: dropping the permutation visibly changes the scores,")
    say("so the round-trip is capable of failing.")
say("")
if ok:
    say("PASS -- the writer and the permutation are correct; Leiden partitions")
    say("       can be compared to MCL's on these numbers.")
else:
    say("FAIL -- do not trust any cross-method comparison until this passes.")
sys.exit(0 if ok else 1)
