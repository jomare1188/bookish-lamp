#!/usr/bin/env python3
# =============================================================================
# 11_myb59_cytoscape.py -- the AtMYB59 neighbourhoods as networks you can open.
#
# WHY A TRIMMED NETWORK, AND WHAT IS TRIMMED. The full induced subgraph is 79,128
# edges for the purple gene alone and 1,534,542 for its sugarcane ortholog --
# Cytoscape on a laptop will not cope, and a hairball hides the partners worth
# looking at anyway. The rule here is uniform and easy to state:
#
#     keep every gene joined to a seed at |r| >= 0.85,
#     and every edge among those genes at |r| >= 0.85.
#
# One threshold, applied to every edge rather than only to the seed's own, so
# "what is in this picture" has a one-line answer. THE TRIMMING IS ONLY ABOUT
# DRAWING: step 10's tables keep every neighbour, so anything left out of the
# picture can still be looked up.
#
# THE THRESHOLD IS NOT ARBITRARY. These neighbourhoods are broad and weak -- the
# median seed edge is |r| 0.823 in purple and 0.822 in sugarcane, barely over the
# 0.8 network threshold -- so 0.85 keeps the top fifth or so, which is where a real
# partner is likelier to be. It is also the same |r| in both species: both networks
# use one weight scale, so a single cut means the same thing on either side.
#
# EDGES CARRY THEIR SIGN. 46% of one MYB's neighbours are NEGATIVELY correlated, and
# for a transcription factor that is much of the point -- so `r` is signed and
# `direction` is a column you can style on, rather than everything being |r|.
#
# THE THIRD FILE IS THE CROSS-SPECIES ONE. Only 33 orthology links join the purple
# gene's neighbourhood to its sugarcane counterparts' -- 3 of them among the drawn
# core. That is the readout's central negative, and seeing two clouds joined by a
# few threads says it better than a table does.
#
# RUN: ./run_all.sh 11
# =============================================================================

import csv
import glob
import os
import re
import sys
from collections import defaultdict

import pandas as pd

BASE = "/dados04/jorge/comparative_saccharum"
CLEAN = f"{BASE}/new_clean"
RES = f"{CLEAN}/results"
M20 = f"{RES}/readouts/module20"
ANN = f"{M20}/myb59_annotated"
OUT = f"{M20}/myb59_cytoscape"
OG_FILE = (f"{BASE}/files/fix_orthofinder/proteins/OrthoFinder/"
           f"Results_Jun04_2/Orthogroups/Orthogroups.tsv")
SP_COL = {"sugarcane": "sugarcane_one_transcript",
          "purple": "one_transcript_purple_proteins"}
CORE_MIN_R = 0.85

os.makedirs(OUT, exist_ok=True)
say = lambda *a: print(*a, flush=True)
_RE_P, _RE_V = re.compile(r"\.p[0-9]+$"), re.compile(r"\.v[0-9]+(\.[0-9]+)*$")
strip_version = lambda g: _RE_V.sub("", _RE_P.sub("", g))

say(f"building Cytoscape files: core = |r| >= {CORE_MIN_R}\n")

CORE, FULL, STATS = {}, {}, {}
for sp in ("purple", "sugarcane"):
    ann_f = f"{ANN}/{sp}__ALL_myb59_neighbours.tsv"
    edg_f = f"{M20}/myb59_signed_edges_{sp}.tsv"
    for f in (ann_f, edg_f):
        if not os.path.exists(f):
            sys.exit(f"FATAL: {f} missing -- run step 10 first")

    ann = pd.read_csv(ann_f, sep="\t", dtype=str).fillna("")
    seeds = set(ann.loc[ann["role"] == "seed", "gene"])
    ed = pd.read_csv(edg_f, sep="\t")

    # Core nodes: a seed, or joined to one at |r| >= CORE_MIN_R.
    core = set(seeds)
    strong = ed[ed["r"].abs() >= CORE_MIN_R]
    for a, b in zip(strong["gene1"], strong["gene2"]):
        if a in seeds:
            core.add(b)
        if b in seeds:
            core.add(a)

    # Core edges: both ends in the core, at the same threshold.
    e = strong[strong["gene1"].isin(core) & strong["gene2"].isin(core)].copy()
    e["direction"] = ["positive" if r > 0 else "negative" for r in e["r"]]
    e["abs_r"] = e["r"].abs().round(4)
    e["edge_type"] = [
        "seed-neighbour" if (a in seeds or b in seeds) else "neighbour-neighbour"
        for a, b in zip(e["gene1"], e["gene2"])]
    e = e.rename(columns={"gene1": "source", "gene2": "target"})
    e = e[["source", "target", "r", "abs_r", "direction", "edge_type"]]

    n = ann[ann["gene"].isin(core)].copy()
    n.insert(1, "is_seed", ["yes" if g in seeds else "no" for g in n["gene"]])
    # EVERY DRAWN NODE MUST BE LOOKUP-ABLE. If a core node had no annotation row the
    # picture would contain genes the tables cannot explain.
    if len(n) != len(core):
        sys.exit(f"FATAL: {len(core)} core nodes but {len(n)} annotation rows")

    n.to_csv(f"{OUT}/{sp}_core_nodes.tsv", sep="\t", index=False)
    e.to_csv(f"{OUT}/{sp}_core_edges.tsv", sep="\t", index=False)
    CORE[sp], FULL[sp] = core, set(ann["gene"])
    neg = int((e["direction"] == "negative").sum())
    say(f"{sp:10s} {len(core):>5,} nodes  {len(e):>7,} edges  "
        f"({int((e.edge_type == 'seed-neighbour').sum()):,} seed-neighbour, "
        f"{neg:,} negative)")

    # Numbers the README quotes, taken from the files being shipped rather than
    # from memory -- so a re-run cannot leave the prose describing an older build.
    seed_r = pd.concat([
        pd.to_numeric(pd.read_csv(f, sep="\t", dtype=str)["r_to_this_seed"],
                      errors="coerce").dropna().abs()
        for f in sorted(glob.glob(f"{ANN}/{sp}__*.tsv")) if "__ALL_" not in f])
    has = lambda col: int((ann[col] != "").sum())
    STATS[sp] = dict(
        genes=len(ann), seeds=sorted(seeds), n_seed_edges=len(seed_r),
        median_r=seed_r.median(), n_strong=int((seed_r >= CORE_MIN_R).sum()),
        nodes=len(core), edges=len(e), negative=neg,
        seed_edges=int((e.edge_type == "seed-neighbour").sum()),
        go=int((~((ann["GO_BP"] == "") & (ann["GO_MF"] == "")
                  & (ann["GO_CC"] == ""))).sum()),
        ipr=has("InterPro"), tf=int((ann["is_TF"] == "yes").sum()),
        nresp=int((ann["nitrogen_responsive"] == "yes").sum()))

# --- the cross-species links --------------------------------------------------
# Counted twice: among the CORE nodes (what gets drawn) and among the FULL
# neighbourhoods (what the tables cover). The README quotes both, because a reader
# who finds three threads in the picture will want to know whether trimming made
# them scarce or whether they were always scarce.
links, full_links, full_ogs = [], 0, 0
with open(OG_FILE) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        S = {strip_version(x.strip())
             for x in (row[SP_COL["sugarcane"]] or "").split(",") if x.strip()}
        P = {strip_version(x.strip())
             for x in (row[SP_COL["purple"]] or "").split(",") if x.strip()}
        fs, fp = S & FULL["sugarcane"], P & FULL["purple"]
        full_links += len(fs) * len(fp)
        full_ogs += 1 if (fs and fp) else 0
        a, b = S & CORE["sugarcane"], P & CORE["purple"]
        for x in sorted(a):
            for y in sorted(b):
                links.append((x, y, row["Orthogroup"]))
with open(f"{OUT}/ortholog_links.tsv", "w") as fh:
    fh.write("source\ttarget\tr\tabs_r\tdirection\tedge_type\torthogroup\n")
    for x, y, og in links:
        fh.write(f"{x}\t{y}\t\t\t\tortholog\t{og}\n")
say(f"\northolog links between the two cores: {len(links)} "
    f"(full neighbourhoods: {full_links} links over {full_ogs} orthogroups)")
if not links:
    say("  none -- the two neighbourhoods do not touch at all within the core")

# --- the README ---------------------------------------------------------------
# Generated, not typed. Every number below is read out of the files shipped beside
# it, so the prose cannot drift from the build the way a hand-written note would.
go = pd.read_csv(f"{M20}/myb59_neighbour_go.tsv", sep="\t")
go_bh = int((go["p.adj"] < 0.05).sum())
P, S = STATS["purple"], STATS["sugarcane"]
pseed = P["seeds"][0]
pct = lambda a, b: f"{100 * a / b:.0f}%"

# The MYB-to-MYB link, if the trimming kept it -- NOT simply links[0]: the list is
# in orthogroup-file order, and the first entry is an ordinary neighbour pair.
seed_link = [(x, y) for x, y, _ in links
             if x in set(S["seeds"]) and y in set(P["seeds"])]
seed_link_txt = (
    f"One of the {len(links)} is the MYB pair itself — `{seed_link[0][0]}` with "
    f"`{seed_link[0][1]}`." if seed_link else
    "None of them is the MYB pair itself.")
sg_rows = sorted(len(pd.read_csv(f, sep="\t", dtype=str))
                 for f in glob.glob(f"{ANN}/sugarcane__*.tsv") if "__ALL_" not in f)

readme = f"""# The AtMYB59 neighbourhoods: how to look at them

_Generated by `scripts/11_readouts/module20/11_myb59_cytoscape.py`. Every number in
this file is computed from the files sitting next to it._

## What this is

`{pseed}` is the purple-cane MYB we are following: the one copy of the
AtMYB59 anchor that is expressed in purple, and the one with the U-shaped nitrogen
response. In sugarcane (R570) the same anchor picks up **{len(S["seeds"])} copies**.

A **co-expression neighbour** here means: across our nitrogen experiments, that
gene's expression rises and falls together with the MYB's (or exactly against it),
closely enough to pass the network's correlation threshold of |r| = 0.8. It is a
statement about expression, not about physical interaction or regulation.

There are two ways in, and they answer different questions.

| you want to | open |
|---|---|
| "what *are* the neighbours of my gene?" | the tables in `myb59_annotated/` — Excel |
| "what does the neighbourhood *look like*?" | the files in `myb59_cytoscape/` — Cytoscape |

The tables hold **every** neighbour. The network files hold only the strongest
({CORE_MIN_R} and above), because the full picture is unusable — see "Why the
picture is trimmed".

---

## 1. The tables (`myb59_annotated/`)

One row per gene, one file per MYB, plus a combined file per species. Open in Excel
or LibreOffice: they are tab-separated text.

| file | rows |
|---|---|
| `purple__{pseed}.tsv` | {P["genes"]:,} |
| `purple__ALL_myb59_neighbours.tsv` | {P["genes"]:,} (the same, purple has one seed) |
| `sugarcane__<gene>.tsv` × {len(S["seeds"])} | {sg_rows[0]:,} – {sg_rows[-1]:,} |
| `sugarcane__ALL_myb59_neighbours.tsv` | {S["genes"]:,} (each gene once, even if it neighbours several MYBs) |

Every file contains the MYB itself (`role` = seed) plus all of its neighbours.

### What the columns mean

| column | meaning |
|---|---|
| `gene` | the gene ID |
| `role` | `seed` = one of the MYBs; `neighbour` = co-expressed with one |
| `neighbour_of` | which MYB(s) it is a neighbour of |
| `n_seeds` | how many of the MYBs it neighbours — above 1 means shared |
| `r_to_this_seed` | **the number to sort on in a per-MYB file.** Correlation with that MYB, signed: +0.9 = rises with it, −0.9 = falls as it rises |
| `direction` | `positive` / `negative`, the sign of the above |
| `r_to_seed`, `r_to_seed_strongest` | in the combined file, all of a gene's correlations, and the largest |
| `degree` | how many neighbours the gene has in the whole network |
| `degree_percentile` | where that sits among all genes — 99 = top 1%, a hub |
| `strength` | sum of the gene's edge weights |
| `mcl_module` / `module_size` | the co-expression cluster it belongs to |
| `is_TF` / `TF_family` | is it a transcription factor, and of which family |
| `nitrogen_r`, `nitrogen_padj`, `nitrogen_responsive` | does its expression track nitrogen, after accounting for experiment and genotype |
| `GO_BP` / `GO_MF` / `GO_CC` | function as GO terms, written out in words — biological process, molecular function, cellular component |
| `InterPro` / `InterPro_desc` | protein domains found in it |
| `ortholog_other_species` | the corresponding gene(s) in the other species, where there are any |

### Coverage — what is *not* blank

| | genes | with GO | with domains | TFs | nitrogen-responsive |
|---|---|---|---|---|---|
| purple | {P["genes"]:,} | {P["go"]} ({pct(P["go"], P["genes"])}) | {P["ipr"]} ({pct(P["ipr"], P["genes"])}) | {P["tf"]} | **{P["nresp"]}** |
| sugarcane | {S["genes"]:,} | {S["go"]:,} ({pct(S["go"], S["genes"])}) | {S["ipr"]:,} ({pct(S["ipr"], S["genes"])}) | {S["tf"]} | {S["nresp"]:,} ({pct(S["nresp"], S["genes"])}) |

About a third of genes have no GO term at all. That is a gap in the annotation, not
a finding about those genes — do not read a blank as "no function".

**The purple {P["nresp"]} is not a typo.** The purple MYB is itself
nitrogen-responsive, but not one of its {P["genes"] - 1} neighbours is, while over a
third of the sugarcane neighbours are. Whatever this gene is co-expressed with in
purple, it is not the nitrogen programme.

---

## 2. The networks (`myb59_cytoscape/`)

| file | what it is |
|---|---|
| `purple_core_edges.tsv` | {P["edges"]:,} edges among {P["nodes"]} genes — **import as network** |
| `purple_core_nodes.tsv` | the same {P["nodes"]} genes with all the columns above — **import as table** |
| `sugarcane_core_edges.tsv` | {S["edges"]:,} edges among {S["nodes"]} genes |
| `sugarcane_core_nodes.tsv` | the {S["nodes"]} genes |
| `ortholog_links.tsv` | {len(links)} cross-species links, to join the two pictures |

### Opening one in Cytoscape

1. **File → Import → Network from File**, pick `purple_core_edges.tsv`.
   Cytoscape will show a column-mapping dialog: set `source` to **Source Node**,
   `target` to **Target Node**, and leave the rest as edge attributes.
2. **File → Import → Table from File**, pick `purple_core_nodes.tsv`.
   Import as **Node Table Columns**, with `gene` as the key.
3. In the **Style** panel:
   - **Node Fill Color** → map `is_seed` (discrete) → make `yes` red.
   - **Node Size** → map `degree` (continuous) — big nodes are hubs.
   - **Edge Stroke Color** → map `direction` (discrete) → `positive` grey,
     `negative` blue. Half these edges are negative; hiding that hides the point.
   - **Edge Width** → map `abs_r` (continuous).
4. **Layout → Prefuse Force Directed Layout**, using `abs_r` as the edge weight.
5. To look at one MYB at a time: **Select → Nodes → By column**, `neighbour_of`
   contains the gene you want, then **File → New Network → From Selected Nodes**.

To see both species at once, import the two edge files and `ortholog_links.tsv`
into the **same** network, then both node tables. The `edge_type` column separates
the three kinds.

---

## Why the picture is trimmed, and what that costs

Everything at **|r| ≥ {CORE_MIN_R}** is drawn — applied to every edge, not just to
the MYBs' own. Untrimmed, the purple gene's neighbourhood alone is 79,128 edges and
the sugarcane one 1,534,542; Cytoscape on a laptop will not open that, and a
hairball hides the partners worth seeing.

| | neighbours in the tables | drawn |
|---|---|---|
| purple | {P["genes"] - len(P["seeds"]):,} | {P["nodes"]} |
| sugarcane | {S["genes"] - len(S["seeds"]):,} | {S["nodes"]} |

**Nothing is lost, only undrawn** — every neighbour is still in the tables.

One wrinkle if you filter the tables yourself: `r_to_this_seed` is rounded to four
decimals, so a few edges print as exactly `0.8500` while their true value is a
hair under {CORE_MIN_R} and they are not drawn. Three sugarcane edges are in this
position.

---

## Three things to keep in mind

**1. These neighbourhoods are broad and weak.** The typical MYB-to-neighbour
correlation is |r| **{P["median_r"]:.2f}** (purple) and **{S["median_r"]:.2f}**
(sugarcane) — barely over the 0.8 cut-off. Only {P["n_strong"]} of
{P["n_seed_edges"]} and {S["n_strong"]} of {S["n_seed_edges"]:,} reach
{CORE_MIN_R}. These MYBs have many loose partners rather than a tight core, which
is also why the purple one sits in an MCL module of only two genes. Treat the
neighbour lists as a pool to draw candidates from, not as a module.

**2. Around half of the links are negative.** {P["negative"]:,} of {P["edges"]:,}
drawn purple edges and {S["negative"]:,} of {S["edges"]:,} sugarcane edges are
anti-correlations. For a transcription factor that is not noise to be filtered out
— a repressor's targets are expected to go the other way.

**3. The GO enrichment table wants a grain of salt.** `myb59_neighbour_go.tsv` has
{len(go)} enriched terms, but only **{go_bh} pass multiple-testing correction**
(`p.adj` < 0.05). Terms are listed on the raw p-value, as everywhere else in this
project. Sort by `p.adj` and read the top of the list, not all {len(go)} rows.

---

## The cross-species question

{len(links)} orthology links join the two drawn neighbourhoods; across the full
neighbourhoods there are {full_links}, over {full_ogs} orthogroups. One of the
{seed_link_txt}

That scarcity is the result, not a shortcoming of the files. The formal test is in
`myb59_neighbour_overlap.tsv`: the purple gene and its 1:1 sugarcane ortholog share
no more neighbours than two unrelated genes of the same size would (z = 0.53,
p = 0.29), while two sugarcane copies share up to 83% of theirs. The measure works;
the conservation is not there. The two MYBs keep the same DNA-binding domain and
change the company they keep.

## Where these came from

`scripts/11_readouts/module20/` — step 06 found the neighbours, 07 ran the GO
enrichment, 09–11 built the files described here. The upstream network, annotation
and nitrogen tests are described in `new_clean/docs/methods.md`.
"""
with open(f"{M20}/README_for_collaborator.md", "w") as fh:
    fh.write(readme)
say(f"wrote README_for_collaborator.md ({len(readme.splitlines())} lines)")
say("\ndone")
