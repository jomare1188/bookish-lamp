#!/usr/bin/env python3
"""
11_percopy_omega.py -- omega for every copy against its sorghum anchor, then a
POWER GATE on whether copy-specific selection is measurable at all.

Each copy gets its own omega from a 2-sequence codeml run against the family's
diploid sorghum anchor. That much is straightforward and reuses the same
MAFFT -> pal2nal -> codeml chain steps 03/04 use, collapsed into one worker so
165k pairs do not each leave four files on disk.

The gate is the point. Copies of a family are ~98% identical and share almost
all of their branch back to sorghum, so their omegas are expected to be nearly
the same by construction, not by biology. Before any within-family test of
"does the copy with the better network position evolve under stronger
constraint", measure whether copies are distinguishable in omega at all:

    ICC = (MSB - MSW) / (MSB + (k-1) MSW)      one-way, families as groups

ICC near 1 means omega is a property of the FAMILY and the copies within it are
interchangeable -- in which case no within-family test can work, and the honest
output is to say so. Only an ICC leaving real within-family variance justifies
escalating to branch models, which is otherwise out of scope.
"""
import os, re, subprocess, sys, tempfile, shutil
from concurrent.futures import ProcessPoolExecutor

CTL = """      seqfile = in.paml
      outfile = out
        noisy = 0
      verbose = 0
      runmode = -2
      seqtype = 1
    CodonFreq = 2
        model = 0
      NSsites = 0
        icode = 0
    fix_kappa = 0
        kappa = 2
    fix_omega = 0
        omega = 0.4
    cleandata = 1
"""
RE_VALS = re.compile(
    r"t=\s*([\d.eE+-]+)\s+S=\s*([\d.eE+-]+)\s+N=\s*([\d.eE+-]+)\s+"
    r"dN/dS=\s*([\d.eE+-]+)\s+dN\s*=\s*([\d.eE+-]+)\s+dS\s*=\s*([\d.eE+-]+)")

PROT = CDS = SBP = SBC = None


def read_fasta(path):
    d, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    d[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        d[name] = "".join(buf)
    return d


def _init(pp, pc, sp, sc):
    global PROT, CDS, SBP, SBC
    PROT, CDS = read_fasta(pp), read_fasta(pc)
    SBP, SBC = read_fasta(sp), read_fasta(sc)


def one(args):
    og, gene, anchor = args
    if gene not in PROT or anchor not in SBP:
        return None
    tmp = tempfile.mkdtemp(prefix="pc_", dir=os.environ["DNDS_TMPDIR"])
    try:
        faa = os.path.join(tmp, "in.faa")
        fna = os.path.join(tmp, "in.fna")
        with open(faa, "w") as f1, open(fna, "w") as f2:
            f1.write(">copy\n%s\n>anchor\n%s\n" % (PROT[gene], SBP[anchor]))
            f2.write(">copy\n%s\n>anchor\n%s\n" % (CDS[gene], SBC[anchor]))
        aln = os.path.join(tmp, "aln.faa")
        with open(aln, "w") as fo:
            r = subprocess.run([os.environ["MAFFT"], "--quiet", "--localpair",
                                "--maxiterate", "1000", "--amino", faa],
                               stdout=fo, stderr=subprocess.DEVNULL, timeout=600)
        if r.returncode != 0:
            return None
        paml = os.path.join(tmp, "in.paml")
        with open(paml, "w") as fo:
            r = subprocess.run([os.environ["PAL2NAL"], aln, fna, "-output", "paml",
                                "-nogap", "-codontable", "1"],
                               stdout=fo, stderr=subprocess.DEVNULL, timeout=600)
        if r.returncode != 0 or os.path.getsize(paml) == 0:
            return None
        with open(paml) as fh:
            head = fh.readline().split()
        if len(head) < 2 or int(head[1]) < 3 * int(os.environ["MIN_CODONS"]):
            return None
        open(os.path.join(tmp, "codeml.ctl"), "w").write(CTL)
        r = subprocess.run([os.environ["CODEML"], "codeml.ctl"], cwd=tmp,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=600)
        outp = os.path.join(tmp, "out")
        if r.returncode != 0 or not os.path.exists(outp):
            return None
        m = RE_VALS.search(open(outp).read())
        if not m:
            return None
        return (og, gene, anchor, int(head[1]) // 3, float(m.group(1)), float(m.group(2)),
                float(m.group(3)), float(m.group(5)), float(m.group(6)), float(m.group(4)))
    except Exception:
        return None
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def icc_oneway(groups):
    """one-way ICC over families; groups = list of lists of omega."""
    groups = [g for g in groups if len(g) >= 2]
    if len(groups) < 10:
        return None
    n = sum(len(g) for g in groups)
    k = len(groups)
    grand = sum(sum(g) for g in groups) / n
    msb = sum(len(g) * (sum(g) / len(g) - grand) ** 2 for g in groups) / (k - 1)
    msw = sum(sum((x - sum(g) / len(g)) ** 2 for x in g) for g in groups) / (n - k)
    # average group size, Shrout-Fleiss correction for unequal sizes
    k0 = (n - sum(len(g) ** 2 for g in groups) / n) / (k - 1)
    denom = msb + (k0 - 1) * msw
    return (msb - msw) / denom if denom > 0 else None, msb, msw, k0, k, n


def main():
    env = os.environ
    env["PAL2NAL"] = open(os.path.join(env["WORKDIR"], "pal2nal.path")).read().strip()
    os.makedirs(env["DNDS_TMPDIR"], exist_ok=True)
    ds_min, ds_max = float(env["DNDS_DS_MIN"]), float(env["DNDS_DS_MAX"])

    for study in ("sugarcane", "purple"):
        tasks = []
        with open(os.path.join(env["OUTDIR"], "families_%s.tsv" % study)) as fh:
            fh.readline()
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[2] not in ("homeolog", "dispersed"):
                    continue
                for g in f[5].split(","):
                    tasks.append((f[0], g, f[6]))
        print("== %s: omega for %d copies against their sorghum anchors" % (study, len(tasks)))
        pp = env["SC_PEP_CLEAN"] if study == "sugarcane" else env["PU_PEP_CLEAN"]
        pc = env["SC_CDS_CLEAN"] if study == "sugarcane" else env["PU_CDS_CLEAN"]
        out = os.path.join(env["OUTDIR"], "percopy_omega_%s.tsv" % study)
        rows = []
        with open(out, "w") as fo:
            fo.write("orthogroup\tgene\tsorghum_anchor\taln_codons\tt\tS\tN\tdN\tdS\tomega\n")
            with ProcessPoolExecutor(max_workers=int(env["DNDS_THREADS"]),
                                     initializer=_init,
                                     initargs=(pp, pc, env["SB_PEP_CLEAN"],
                                               env["SB_CDS_CLEAN"])) as ex:
                for r in ex.map(one, tasks, chunksize=16):
                    if r is None:
                        continue
                    fo.write("%s\t%s\t%s\t%d\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\n" % r)
                    rows.append(r)
        print("   %d / %d copies produced an omega" % (len(rows), len(tasks)))

        # --- the gate --------------------------------------------------------
        from collections import defaultdict
        fam = defaultdict(list)
        for r in rows:
            if ds_min <= r[8] <= ds_max and 0 < r[9] <= float(env["DNDS_OMEGA_FLAG"]):
                fam[r[0]].append(r[9])
        res = icc_oneway(list(fam.values()))
        gate = os.path.join(env["OUTDIR"], "percopy_omega_icc_%s.tsv" % study)
        if res is None:
            print("   ICC: not computable")
            open(gate, "w").write("metric\tvalue\nicc\tNA\n")
            continue
        icc, msb, msw, k0, k, n = res
        med = sorted(x for v in fam.values() for x in v)
        print("   median omega %.4f over %d copies in %d families (dS-filtered)"
              % (med[len(med) // 2], n, k))
        print("   ICC(1) = %.4f   [MS_between %.4f, MS_within %.4f, mean family size %.2f]"
              % (icc, msb, msw, k0))
        verdict = ("GATE SHUT: omega is essentially a property of the family; copies "
                   "are not distinguishable, so no within-family selection test can work"
                   if icc >= 0.75 else
                   "GATE OPEN: there is real within-family omega variance to explain")
        print("   %s" % verdict)
        with open(gate, "w") as fo:
            fo.write("metric\tvalue\n")
            for kk, vv in (("icc", "%.6f" % icc), ("ms_between", "%.6f" % msb),
                           ("ms_within", "%.6f" % msw), ("mean_family_size", "%.4f" % k0),
                           ("n_families", k), ("n_copies", n),
                           ("median_omega", "%.6f" % med[len(med) // 2]),
                           ("gate", "shut" if icc >= 0.75 else "open")):
                fo.write("%s\t%s\n" % (kk, vv))


if __name__ == "__main__":
    main()
