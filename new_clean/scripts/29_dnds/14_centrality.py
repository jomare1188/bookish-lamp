#!/usr/bin/env python3
"""
14_centrality.py -- a second description of network position: local clustering.

The constraint analysis so far tests exactly one network property, degree. This
adds the local clustering coefficient -- of all the pairs of a gene's neighbours,
what fraction are themselves connected. It describes the *shape* of a gene's
neighbourhood rather than its size, and unlike the path-based centralities it
stays meaningful in a graph this dense.

WHY NOT BETWEENNESS OR CLOSENESS, and it is not mainly about cost: sugarcane has
mean degree 1,475 (density 1.4%), so the effective diameter is 2-3 hops. When
almost every pair of nodes is two steps apart, closeness has nearly no variance
between nodes and betweenness collapses into a function of degree and clustering.
They would carry little information here even if they were free. (Exact
betweenness is also O(V*E) ~ 7.8e12 for this graph, so it is not free.)

WHAT THIS SCRIPT DOES NOT DO: interpret the number. In dense graphs clustering
falls off with degree roughly as C(k) ~ 1/k, so clustering and degree are
mechanically anti-correlated and a raw clustering-vs-constraint correlation would
largely restate the degree result. Step 15 screens for that before testing
anything, and reports only the partial effect after degree.

The graph is not rebuilt: the SBM stage already wrote this network in graph-tool
form. That it IS this network is verified here, not assumed.
"""
import os
import sys

import numpy as np
from graph_tool import load_graph
from graph_tool.clustering import local_clustering
from graph_tool.topology import kcore_decomposition


def die(msg):
    sys.exit("FATAL: " + msg)


def main():
    env = os.environ
    graph_path = env["GRAPH_sugarcane"]
    metrics_path = env["NODE_METRICS_sugarcane"]
    out = os.path.join(env["OUTDIR"], "centrality_sugarcane.tsv")

    print("== loading %s" % graph_path)
    g = load_graph(graph_path)
    print("   vertices %d   edges %d   directed %s"
          % (g.num_vertices(), g.num_edges(), g.is_directed()))
    if g.is_directed():
        die("graph is directed; the clustering definition below assumes undirected")

    names = np.array([g.vp["name"][v] for v in g.vertices()])
    deg = g.get_out_degrees(g.get_vertices())

    # --- is this actually the pipeline's network? -------------------------
    # The graph comes from the SBM stage, which is outside this pipeline. If it
    # is a different build -- a different threshold, a different gene set -- then
    # every number below would be quietly about the wrong network.
    print("== verifying against %s" % os.path.basename(metrics_path))
    ref = {}
    with open(metrics_path) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            ref[f[0]] = int(f[1])
    if len(ref) != g.num_vertices():
        die("vertex count %d != %d genes in node_metrics" % (g.num_vertices(), len(ref)))
    missing = [n for n in names if n not in ref]
    if missing:
        die("%d graph vertices are absent from node_metrics (e.g. %s)"
            % (len(missing), missing[0]))
    mism = [(n, int(d), ref[n]) for n, d in zip(names, deg) if ref[n] != int(d)]
    if mism:
        die("degree disagrees for %d genes (e.g. %s: graph %d vs metrics %d)"
            % (len(mism), mism[0][0], mism[0][1], mism[0][2]))
    print("   vertex count, gene set and per-gene degree all match")

    # Self-loops and parallel edges would both distort clustering. Neither should
    # exist in a thresholded correlation network, so check rather than assume.
    n_self = sum(1 for e in g.edges() if e.source() == e.target())
    if n_self:
        die("%d self-loops present; clustering would be distorted" % n_self)

    # --- the measure -------------------------------------------------------
    print("== local clustering (unweighted, C++/OpenMP)")
    # Unweighted deliberately: this is the standard definition and the one the
    # C(k) ~ 1/k literature is stated in. Edge weights here are a transformed
    # correlation, so a weighted variant (Barrat) would need its own argument
    # rather than being substituted silently.
    clust = local_clustering(g).a.copy()
    print("   mean %.4f   median %.4f   range [%.4f, %.4f]"
          % (clust.mean(), np.median(clust), clust.min(), clust.max()))
    if clust.min() < 0 or clust.max() > 1:
        die("clustering outside [0,1]; the measure is wrong")

    cores = None
    if env.get("COMPUTE_CORENESS", "0") == "1":
        print("== k-core decomposition")
        cores = kcore_decomposition(g).a.copy()
        print("   max core %d   median %d" % (cores.max(), int(np.median(cores))))

    with open(out, "w") as fh:
        fh.write("gene\tdegree\tclustering\tcoreness\n")
        for i, n in enumerate(names):
            fh.write("%s\t%d\t%.6f\t%s\n"
                     % (n, deg[i], clust[i], cores[i] if cores is not None else "NA"))
    print("== wrote %d genes -> %s" % (len(names), out))
    print("   global mean clustering %.4f  (RESULTS_2026-07-16.md reports 0.690 "
          "transitivity for this network -- they are different statistics, but a "
          "wildly different value would mean the graph or the measure is wrong)"
          % clust.mean())


if __name__ == "__main__":
    main()
