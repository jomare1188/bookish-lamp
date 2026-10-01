#!/usr/bin/env python3
# =============================================================================
# 64_og_crosswalk.py -- bridge the two orthogroup namespaces, and annotate the
# 2-species orthogroups with what the dN/dS stage knows about them.
#
# WHY THIS EXISTS. There are TWO OrthoFinder runs in this project and their
# orthogroup ids are not comparable:
#
#   2-species  files/fix_orthofinder/.../Results_Jun04_2   94,229 OGs
#              backs EVERY conservation result (06, 08, 13, 61, 62) and the
#              hardcoded 114,681-pair / 43,390-orthogroup contract
#   3-species  files/orthofinder_3sp/out/Results_Aug31     50,494 OGs
#              backs families_*, kmer_uniqueness_*, percopy_omega_*, triplets
#
# `OG0000017` exists in both and means DIFFERENT GENE SETS. Joining on the id
# string is silently wrong, so this joins on GENE MEMBERSHIP instead.
#
# THE TWO RUNS ARE NOT ALTERNATIVE LABELLINGS -- one is finer. Measured here:
#   2sp -> 3sp   95.8% of mappable 2sp OGs land in exactly ONE 3sp OG
#   3sp -> 2sp   59.1% -- adding sorghum MERGED groups the 2-species run kept apart
# So the 2-species run is the finer partition and is the unit of analysis; this
# table is the lookup that lets a 2sp orthogroup inherit 3sp attributes.
#
# WHAT NEEDS THE CROSSWALK AND WHAT DOES NOT. Only GROUP-level attributes do --
# the sorghum anchor and the family class, which are properties of a 3sp
# orthogroup. GENE-level attributes (frac_unique, per-copy omega) are keyed on
# the gene, so they join straight onto the 2sp orthogroup with NO ambiguity at
# all. That distinction is why the ambiguous 4.2% costs so little: it withholds
# an anchor, never a mappability or an omega.
#
# AMBIGUITY IS FLAGGED, NEVER RESOLVED. A 2sp OG spanning several 3sp OGs gets
# ambiguous=TRUE, og_3sp="" and no anchor. Picking the largest overlap would
# invent a fact.
#
# RUN: through run.sh -> ./run.sh ogcrosswalk
# =============================================================================

import csv
import os
import statistics as st
import sys
from collections import defaultdict

STUDIES = ("sugarcane", "purple")

# Measured 2026-10-01 on Results_Jun04_2 + Results_Aug31. These are the shape of
# the join, not of the biology: a change means an input moved, which must be
# noticed rather than absorbed.
EXPECT = {
    "og_2sp": 94_229,
    "both_species": 71_587,
    "mappable": 88_593,
    "one_to_one": 84_852,
}


def env_req(name):
    v = os.environ.get(name, "")
    if not v:
        sys.exit(f"FATAL: {name} is not set -- launch through run.sh")
    return v


def read_tsv(path):
    """OrthoFinder writes CRLF; splitting on \\n alone glues \\r to the last column."""
    return (open(path, "rb").read().decode()
            .replace("\r\n", "\n").rstrip("\n").split("\n"))


def cells(s):
    return [x.strip() for x in (s or "").split(",") if x.strip()]


def say(*a):
    print(*a, flush=True)


# --- inputs -------------------------------------------------------------------
OG2_FILE = env_req("CLEAN_ORTHOGROUPS")
OG3_FILE = env_req("CLEAN_ORTHOGROUPS_3SP")
DNDS_DIR = env_req("CLEAN_DNDS_DIR")
OUT_FILE = env_req("CLEAN_OUT_FILE")
SP_COL = {"sugarcane": os.environ.get("CLEAN_OG_SPECIES_SUGARCANE",
                                      "sugarcane_one_transcript"),
          "purple": os.environ.get("CLEAN_OG_SPECIES_PURPLE",
                                   "one_transcript_purple_proteins")}

say("building the 2sp <-> 3sp orthogroup crosswalk\n")

# --- the 2-species table: the unit of analysis --------------------------------
lines = read_tsv(OG2_FILE)
h2 = lines[0].split("\t")
i2 = {c: h2.index(SP_COL[c]) for c in STUDIES}
og2_genes = {}                      # og_2sp -> {study: [genes]}
gene_og2 = {}                       # gene -> og_2sp
for ln in lines[1:]:
    p = ln.split("\t")
    og = p[0]
    g = {c: cells(p[i2[c]]) if i2[c] < len(p) else [] for c in STUDIES}
    og2_genes[og] = g
    for c in STUDIES:
        for x in g[c]:
            gene_og2[x] = og
if len(og2_genes) != EXPECT["og_2sp"]:
    sys.exit(f"FATAL: {len(og2_genes):,} orthogroups in the 2-species table, "
             f"expected {EXPECT['og_2sp']:,} -- did ORTHOGROUPS move?")
both = sum(1 for g in og2_genes.values() if g["sugarcane"] and g["purple"])
if both != EXPECT["both_species"]:
    sys.exit(f"FATAL: {both:,} orthogroups span both species, "
             f"expected {EXPECT['both_species']:,}")
say(f"2sp: {len(og2_genes):,} orthogroups | {both:,} span both species")

# --- the 3-species table: where the anchors and the family classes live -------
lines = read_tsv(OG3_FILE)
h3 = lines[0].split("\t")
c3 = {k: next(i for i, x in enumerate(h3) if k in x.lower())
      for k in ("sugarcane", "purple", "sorghum")}
gene_og3 = {}
og3_anchor = {}
for ln in lines[1:]:
    p = ln.split("\t")
    og = p[0]
    sb = cells(p[c3["sorghum"]]) if c3["sorghum"] < len(p) else []
    # One sorghum gene is the families stage's own criterion for a usable anchor
    # (29_dnds/08_homeolog_families.py:72); carry it only when it is unambiguous.
    og3_anchor[og] = sb[0] if len(sb) == 1 else ""
    for k in ("sugarcane", "purple"):
        if c3[k] < len(p):
            for x in cells(p[c3[k]]):
                gene_og3[x] = og
say(f"3sp: {len(lines) - 1:,} orthogroups | {len(gene_og3):,} Saccharum genes placed")

# --- the crosswalk, by gene membership ---------------------------------------
fan = defaultdict(set)              # og_2sp -> {og_3sp}
shared = defaultdict(int)           # og_2sp -> genes with a 3sp placement
for gene, og2 in gene_og2.items():
    og3 = gene_og3.get(gene)
    if og3 is not None:
        fan[og2].add(og3)
        shared[og2] += 1
one_to_one = sum(1 for v in fan.values() if len(v) == 1)
if (len(fan), one_to_one) != (EXPECT["mappable"], EXPECT["one_to_one"]):
    sys.exit(f"FATAL: crosswalk shape changed -- {len(fan):,} mappable / "
             f"{one_to_one:,} one-to-one, expected "
             f"{EXPECT['mappable']:,} / {EXPECT['one_to_one']:,}")
amb = len(fan) - one_to_one
say(f"crosswalk: {len(fan):,} of {len(og2_genes):,} orthogroups reach the 3sp run | "
    f"{one_to_one:,} unambiguous ({100 * one_to_one / len(fan):.1f}%) | "
    f"{amb:,} ambiguous ({100 * amb / len(fan):.1f}%)")

# --- group-level attributes: family class, via the crosswalk -----------------
fam_class = {s: {} for s in STUDIES}        # og_3sp -> class
for s in STUDIES:
    f = f"{DNDS_DIR}/families_{s}.tsv"
    if not os.path.exists(f):
        sys.exit(f"FATAL: {f} missing -- run the dN/dS stage's step 08 first")
    with open(f) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            fam_class[s][row["orthogroup"]] = row["class"]
    say(f"  families_{s}.tsv: {len(fam_class[s]):,} 3sp families")

# --- gene-level attributes: NO crosswalk needed, these key on the gene -------
uniq = {s: {} for s in STUDIES}             # gene -> frac_unique
omega = {s: {} for s in STUDIES}            # gene -> omega
for s in STUDIES:
    f = f"{DNDS_DIR}/kmer_uniqueness_{s}.tsv"
    if os.path.exists(f):
        with open(f) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                try:
                    uniq[s][row["gene"]] = float(row["frac_unique"])
                except (TypeError, ValueError):
                    pass
    f = f"{DNDS_DIR}/percopy_omega_{s}.tsv"
    if os.path.exists(f):
        with open(f) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                try:
                    omega[s][row["gene"]] = float(row["omega"])
                except (TypeError, ValueError):
                    pass
    say(f"  {s}: frac_unique on {len(uniq[s]):,} genes | "
        f"omega on {len(omega[s]):,} genes")

# --- write --------------------------------------------------------------------
os.makedirs(os.path.dirname(OUT_FILE), exist_ok=True)
COLS = ["og_2sp", "n_sugarcane", "n_purple", "both_species",
        "og_3sp", "n_3sp_touched", "ambiguous", "n_genes_crosswalked",
        "sorghum_anchor", "class_sugarcane", "class_purple",
        "frac_unique_min_sugarcane", "frac_unique_min_purple",
        "n_zero_unique_sugarcane", "n_zero_unique_purple",
        "omega_median_sugarcane", "omega_median_purple"]

fmt = lambda v, d=4: "" if v is None else f"{v:.{d}f}"
n_anchor = n_amb = 0
with open(OUT_FILE, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(COLS)
    for og2 in sorted(og2_genes):
        g = og2_genes[og2]
        hits = fan.get(og2, set())
        is_amb = len(hits) > 1
        og3 = next(iter(hits)) if len(hits) == 1 else ""
        n_amb += is_amb
        anchor = og3_anchor.get(og3, "") if og3 else ""
        n_anchor += bool(anchor)
        row = [og2, len(g["sugarcane"]), len(g["purple"]),
               "TRUE" if (g["sugarcane"] and g["purple"]) else "FALSE",
               og3, len(hits), "TRUE" if is_amb else "FALSE", shared.get(og2, 0),
               anchor,
               fam_class["sugarcane"].get(og3, "") if og3 else "",
               fam_class["purple"].get(og3, "") if og3 else ""]
        for s in STUDIES:
            vals = [uniq[s][x] for x in g[s] if x in uniq[s]]
            row.append(fmt(min(vals)) if vals else "")
        for s in STUDIES:
            row.append(sum(1 for x in g[s]
                           if uniq[s].get(x, 1.0) == 0.0))
        for s in STUDIES:
            vals = [omega[s][x] for x in g[s] if x in omega[s]]
            row.append(fmt(st.median(vals)) if vals else "")
        w.writerow(row)

say(f"\nwrote {os.path.basename(OUT_FILE)}  ({len(og2_genes):,} rows)")
say(f"  {n_anchor:,} orthogroups carry a sorghum anchor")
say(f"  {n_amb:,} flagged ambiguous -- no 3sp id, no anchor, gene-level columns intact")
say("\ndone")
