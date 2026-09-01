#!/usr/bin/env python3
"""
04b_yn00_check.py -- independent estimator on the same alignments.

yn00 (Yang & Nielsen 2000) is a counting method; codeml -2 is maximum
likelihood. They share the alignments but almost nothing else, so agreement is
evidence the plumbing is right -- ids joined correctly, frames intact, pal2nal
mapping sound. Disagreement points at the alignment, not at biology.

Run on a subset: this is a control, not a result.
"""
import os, re, random, subprocess, tempfile, shutil, sys
from concurrent.futures import ProcessPoolExecutor

CTL = """      seqfile = in.paml
      outfile = out
      verbose = 0
        icode = 0
    weighting = 0
   commonf3x4 = 0
"""
# yn00's table:  seq. seq.  S  N  t  kappa  omega  dN +- SE  dS +- SE
RE_ROW = re.compile(
    r"^\s*(\d+)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+"
    r"(-?[\d.]+)\s+(-?[\d.]+)\s*\+-\s*[\d.]+\s+(-?[\d.]+)", re.M)


def one(args):
    og, d = args
    paml = os.path.join(d, og + ".paml")
    tmp = tempfile.mkdtemp(prefix="yn_" + og + "_", dir=os.environ["DNDS_TMPDIR"])
    try:
        shutil.copy(paml, os.path.join(tmp, "in.paml"))
        open(os.path.join(tmp, "yn00.ctl"), "w").write(CTL)
        subprocess.run([os.environ["YN00"], "yn00.ctl"], cwd=tmp,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=300)
        outp = os.path.join(tmp, "out")
        if not os.path.exists(outp):
            return og, None
        names = [l.split()[0] for l in open(paml).read().split("\n")[1:] if l.strip()
                 and not l[0].isdigit() and not set(l.strip()) <= set("ACGTN-")]
        txt = open(outp).read()
        # sequence order in the file is sugarcane, purple, sorghum (03_explode)
        idx = {1: "sugarcane", 2: "purple", 3: "sorghum"}
        res = {}
        for m in RE_ROW.finditer(txt):
            a, b = idx.get(int(m.group(1))), idx.get(int(m.group(2)))
            if a and b:
                res[frozenset((a, b))] = float(m.group(7))
        return og, res.get(frozenset(("sugarcane", "sorghum")))
    except Exception:
        return og, None
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main():
    env = os.environ
    n_want = int(env.get("YN00_N", "2000"))
    tasks = []
    for line in open(os.path.join(env["WORKDIR"], "og_list.txt")):
        og, d = line.rstrip("\n").split("\t")
        if os.path.exists(os.path.join(d, og + ".paml")):
            tasks.append((og, d))
    random.Random(0).shuffle(tasks)
    tasks = tasks[:n_want]
    print("== yn00 cross-check on %d alignments" % len(tasks))

    ml, ds = {}, {}
    with open(os.path.join(env["OUTDIR"], "dnds_raw.tsv")) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        i_og, i_om, i_ds = h.index("orthogroup"), h.index("sc_sb_omega"), h.index("sc_sb_dS")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            ml[f[i_og]] = float(f[i_om])
            ds[f[i_og]] = float(f[i_ds])

    yn = {}
    with ProcessPoolExecutor(max_workers=int(env["DNDS_THREADS"])) as ex:
        for og, om in ex.map(one, tasks, chunksize=8):
            if om is not None:
                yn[og] = om

    def comparable(keep):
        xs, ys = [], []
        for og, om in yn.items():
            if og in ml and keep(og) and 0 < om < 5 and 0 < ml[og] < 5:
                xs.append(ml[og]); ys.append(om)
        return xs, ys

    pairs = list(zip(*comparable(lambda og: True))) or []
    if len(yn) < 10:
        sys.exit("FATAL: yn00 returned too few comparable values (%d)" % len(yn))

    def spearman(xs, ys):
        def rank(v):
            order = sorted(range(len(v)), key=lambda i: v[i])
            r = [0.0] * len(v)
            for pos, i in enumerate(order):
                r[i] = pos
            return r
        rx, ry = rank(xs), rank(ys)
        n = len(xs)
        mx, my = sum(rx) / n, sum(ry) / n
        num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
        den = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5
        return num / den if den else float("nan")

    # Report BOTH: every alignment, and the dS window the analysis actually
    # uses. The two estimators are expected to drift apart on saturated genes
    # (yn00's counting correction degrades as dS grows), and those genes are
    # filtered out in 06 -- so the filtered rho is the one that gates.
    lo, hi = float(env["DNDS_DS_MIN"]), float(env["DNDS_DS_MAX"])
    out = os.path.join(env["OUTDIR"], "yn00_crosscheck.tsv")
    rows = []
    for label, keep in (("all", lambda og: True),
                        ("ds_filtered", lambda og: lo <= ds[og] <= hi)):
        xs, ys = comparable(keep)
        rho = spearman(xs, ys) if len(xs) >= 10 else float("nan")
        rows.append((label, len(xs), rho))
        print("   %-12s n = %4d   Spearman rho = %.4f" % (label, len(xs), rho))
    with open(out, "w") as fh:
        fh.write("set\tn\tspearman_rho_codeml_vs_yn00\n")
        for label, n, rho in rows:
            fh.write("%s\t%d\t%.4f\n" % (label, n, rho))
    gate = rows[1][2]
    print("   -> %s" % out)
    print("   %s" % ("PASS: the two estimators agree on the analysed set (rho > 0.95)"
                     if gate > 0.95 else
                     "WARNING: rho <= 0.95 on the ANALYSED set -- inspect alignments "
                     "before trusting omega"))


if __name__ == "__main__":
    main()
