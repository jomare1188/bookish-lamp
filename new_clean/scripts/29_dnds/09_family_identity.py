#!/usr/bin/env python3
"""
09_family_identity.py -- how similar are the copies within each family?

Two identities per copy pair, and the distinction matters:

  pid_prot  amino-acid identity. What you would quote as "these are the same
            gene".
  pid_cds   nucleotide identity over the SAME alignment. This is the one that
            governs whether salmon can tell the copies apart, because salmon
            reads DNA. Copies can be 100% identical in protein and comfortably
            distinct in CDS, so the protein figure alone would overstate the
            quantification problem.

The CDS alignment is derived from the protein alignment by codon mapping (3 nt
per aligned residue) rather than by aligning nucleotides independently. Step 01
already guaranteed translate(CDS) == protein for every gene, so the mapping is
exact and the two identities are measured over identical columns.

Each PAIR is also classified, which the family-level label cannot do:
  homeolog  different haplotypes of the same chromosome -- a polyploid copy
  tandem    the same haplotype -- a local duplication sitting inside a
            homeolog family (e.g. OG0000138 has three copies on 05C alone)
"""
import os, re, subprocess, sys, tempfile, shutil
from concurrent.futures import ProcessPoolExecutor

PAT = {
    "sugarcane": re.compile(r"^SoffiXsponR570\.(\d+)([A-Z])g(\d+)$"),
    "purple":    re.compile(r"^Soffic\.(\d+)[A-Z](\d+)-(\d+)([A-Z])$"),
}


def hap(study, gene):
    m = PAT[study].match(gene)
    if not m:
        return None
    return (m.group(1), m.group(2)) if study == "sugarcane" else (m.group(3), m.group(4))


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


PROT = CDS = SBP = SBC = None
STUDY = None


def _init(study, pp, pc, sp, sc):
    global PROT, CDS, SBP, SBC, STUDY
    STUDY = study
    PROT, CDS = read_fasta(pp), read_fasta(pc)
    SBP, SBC = read_fasta(sp), read_fasta(sc)


def codonise(aln_prot, cds):
    """aligned protein -> aligned CDS, 3 nt per residue, '---' per gap."""
    out, k = [], 0
    for aa in aln_prot:
        if aa == "-":
            out.append("---")
        else:
            out.append(cds[k:k + 3] if k + 3 <= len(cds) else "NNN")
            k += 3
    return "".join(out)


def ident(a, b, step):
    """identity over positions/columns where NEITHER is a gap."""
    n = same = 0
    for i in range(0, len(a), step):
        x, y = a[i:i + step], b[i:i + step]
        if "-" in x or "-" in y:
            continue
        n += 1
        if x == y:
            same += 1
    return (same / n, n) if n else (None, 0)


def one(args):
    og, copies, anchor = args
    seqs = {g: PROT[g] for g in copies if g in PROT}
    if len(seqs) < 2 or anchor not in SBP:
        return og, [], []
    seqs["ANCHOR"] = SBP[anchor]
    tmp = tempfile.mkdtemp(prefix="fam_" + og + "_", dir=os.environ["DNDS_TMPDIR"])
    try:
        fa = os.path.join(tmp, "in.faa")
        with open(fa, "w") as fh:
            for k, v in seqs.items():
                fh.write(">%s\n%s\n" % (k, v))
        # --auto, not L-INS-i: families reach 57 copies and this is a distance
        # measurement, not a phylogeny.
        r = subprocess.run([os.environ["MAFFT"], "--quiet", "--auto", "--amino", fa],
                           capture_output=True, text=True, timeout=1800)
        if r.returncode != 0:
            return og, [], []
        aln, name, buf = {}, None, []
        for line in r.stdout.split("\n"):
            if line.startswith(">"):
                if name:
                    aln[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            elif line.strip():
                buf.append(line.strip())
        if name:
            aln[name] = "".join(buf)

        acds = {}
        for k in aln:
            src = SBC[anchor] if k == "ANCHOR" else CDS.get(k)
            if src is None:
                return og, [], []
            acds[k] = codonise(aln[k], src)

        keys = [k for k in aln if k != "ANCHOR"]
        pairs, anchors = [], []
        for i in range(len(keys)):
            a = keys[i]
            pa, na = ident(aln[a], aln["ANCHOR"], 1)
            ca, _ = ident(acds[a], acds["ANCHOR"], 3)
            if pa is not None:
                anchors.append((og, a, "%.5f" % pa, "%.5f" % (ca if ca is not None else 0), na))
            for j in range(i + 1, len(keys)):
                b = keys[j]
                pp, np_ = ident(aln[a], aln[b], 1)
                cc, _ = ident(acds[a], acds[b], 3)
                if pp is None or cc is None:
                    continue
                ha, hb = hap(STUDY, a), hap(STUDY, b)
                if ha and hb and ha[0] == hb[0]:
                    cl = "tandem" if ha[1] == hb[1] else "homeolog"
                elif ha and hb:
                    cl = "dispersed"
                else:
                    cl = "unplaced"
                pairs.append((og, a, b, cl, "%.5f" % pp, "%.5f" % cc, np_))
        return og, pairs, anchors
    except Exception:
        return og, [], []
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main():
    env = os.environ
    os.makedirs(env["DNDS_TMPDIR"], exist_ok=True)
    for study in ("sugarcane", "purple"):
        fam_path = os.path.join(env["OUTDIR"], "families_%s.tsv" % study)
        tasks = []
        with open(fam_path) as fh:
            fh.readline()
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[2] in ("homeolog", "dispersed"):
                    tasks.append((f[0], f[5].split(","), f[6]))
        print("== %s: %d families" % (study, len(tasks)))
        pp = env["SC_PEP_CLEAN"] if study == "sugarcane" else env["PU_PEP_CLEAN"]
        pc = env["SC_CDS_CLEAN"] if study == "sugarcane" else env["PU_CDS_CLEAN"]
        outp = os.path.join(env["OUTDIR"], "family_pairs_%s.tsv" % study)
        outa = os.path.join(env["OUTDIR"], "family_anchor_%s.tsv" % study)
        # Write to .part and rename only on success. A half-written table under
        # the final name is worse than no table at all: step 12 would read it as
        # complete and quietly analyse a fraction of the data.
        npair = nfail = 0
        with open(outp + ".part", "w") as fp, open(outa + ".part", "w") as fa:
            fp.write("orthogroup\tgene_a\tgene_b\tpair_class\tpid_prot\tpid_cds\taln_cols\n")
            fa.write("orthogroup\tgene\tpid_prot_anchor\tpid_cds_anchor\taln_cols\n")
            with ProcessPoolExecutor(max_workers=int(env["DNDS_THREADS"]),
                                     initializer=_init,
                                     initargs=(study, pp, pc, env["SB_PEP_CLEAN"],
                                               env["SB_CDS_CLEAN"])) as ex:
                for og, pairs, anchors in ex.map(one, tasks, chunksize=8):
                    if not pairs:
                        nfail += 1
                    for r in pairs:
                        fp.write("%s\t%s\t%s\t%s\t%s\t%s\t%d\n" % r)
                        npair += 1
                    for r in anchors:
                        fa.write("%s\t%s\t%s\t%s\t%d\n" % r)
        os.replace(outp + ".part", outp)
        os.replace(outa + ".part", outa)
        print("   %d copy pairs, %d families produced nothing" % (npair, nfail))
        print("   -> %s" % outp)


if __name__ == "__main__":
    main()
