#!/usr/bin/env python3
"""
04_codeml.py -- pairwise dN/dS for every codon alignment.

TWO estimators, each used where it is the right tool:

  codeml runmode=-2 (pairwise ML) on the 3-sequence alignment returns ALL
  THREE pairs in one run -- sugarcane/sorghum, purple/sorghum AND
  sugarcane/purple. The two outgroup comparisons are the per-gene omega the
  analysis is built on: sorghum is ~8 Myr away, dS lands around 0.15, and the
  estimate is stable.

  NG86 counting, implemented here, gives EXACT INTEGER counts of synonymous
  and non-synonymous differences (Sd, Nd) and site counts (S, N) for the
  sugarcane/purple pair. codeml does not report those integers, and they are
  precisely what the aggregation in 06 needs: R570 and LA purple are 95-99%
  identical, so per-gene omega between them is undefined for a large share of
  genes, and the only honest route is to SUM counts within a bin and form one
  ratio per bin. Ratios cannot be summed; counts can.
"""
import os, re, subprocess, sys, tempfile, shutil, glob
from concurrent.futures import ProcessPoolExecutor

CODONS = {}
_B = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, b1 in enumerate(_B):
    for j, b2 in enumerate(_B):
        for k, b3 in enumerate(_B):
            CODONS[b1 + b2 + b3] = _AA[i * 16 + j * 4 + k]

# --- NG86 site counts, precomputed once per codon --------------------------
_SITES = {}
for c, aa in CODONS.items():
    if aa == "*":
        continue
    syn = 0.0
    for pos in range(3):
        for b in "TCAG":
            if b == c[pos]:
                continue
            alt = c[:pos] + b + c[pos + 1:]
            if CODONS[alt] == aa:
                syn += 1.0 / 3.0
    _SITES[c] = (syn, 3.0 - syn)          # (S sites, N sites)


def _path_counts(c1, c2, diffs):
    """Average syn/nonsyn steps over every mutational pathway between codons."""
    if not diffs:
        return 0.0, 0.0
    tot_s = tot_n = 0.0
    n_paths = 0
    def walk(cur, remaining, s, n):
        nonlocal tot_s, tot_n, n_paths
        if not remaining:
            tot_s += s; tot_n += n; n_paths += 1
            return
        for k, pos in enumerate(remaining):
            nxt = cur[:pos] + c2[pos] + cur[pos + 1:]
            if CODONS.get(nxt, "*") == "*":
                continue                   # never route through a stop codon
            if CODONS[nxt] == CODONS[cur]:
                walk(nxt, remaining[:k] + remaining[k + 1:], s + 1, n)
            else:
                walk(nxt, remaining[:k] + remaining[k + 1:], s, n + 1)
    walk(c1, diffs, 0.0, 0.0)
    if n_paths == 0:                       # every route passes a stop
        return 0.0, float(len(diffs))
    return tot_s / n_paths, tot_n / n_paths


def ng86(a, b):
    """-> (S, N, Sd, Nd) for two aligned, gap-free, in-frame sequences."""
    S = N = Sd = Nd = 0.0
    for i in range(0, len(a), 3):
        c1, c2 = a[i:i + 3], b[i:i + 3]
        if c1 not in _SITES or c2 not in _SITES:
            continue                       # stop codon or ambiguity: skip site
        s1, n1 = _SITES[c1]; s2, n2 = _SITES[c2]
        S += (s1 + s2) / 2.0
        N += (n1 + n2) / 2.0
        diffs = [p for p in range(3) if c1[p] != c2[p]]
        if diffs:
            ds, dn = _path_counts(c1, c2, diffs)
            Sd += ds; Nd += dn
    return S, N, Sd, Nd


def read_paml(path):
    toks = open(path).read().split()
    n, ln = int(toks[0]), int(toks[1])
    i, seqs = 2, {}
    while len(seqs) < n:
        name = toks[i]; i += 1
        buf = ""
        while len(buf) < ln:
            buf += toks[i]; i += 1
        seqs[name] = buf
    return ln, seqs


CTL = """      seqfile = {seq}
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

# "2 (purple) ... 1 (sugarcane)"  then  "t= .. S= .. N= .. dN/dS= .. dN = .. dS = .."
RE_PAIR = re.compile(r"^\s*\d+\s+\((\S+)\)\s+\.\.\.\s+\d+\s+\((\S+)\)", re.M)
RE_VALS = re.compile(
    r"t=\s*([\d.eE+-]+)\s+S=\s*([\d.eE+-]+)\s+N=\s*([\d.eE+-]+)\s+"
    r"dN/dS=\s*([\d.eE+-]+)\s+dN\s*=\s*([\d.eE+-]+)\s+dS\s*=\s*([\d.eE+-]+)")


def parse_codeml(text):
    out = {}
    pos = [(m.start(), m.group(1), m.group(2)) for m in RE_PAIR.finditer(text)]
    for k, (start, a, b) in enumerate(pos):
        end = pos[k + 1][0] if k + 1 < len(pos) else len(text)
        m = RE_VALS.search(text, start, end)
        if m:
            out[frozenset((a, b))] = dict(
                t=float(m.group(1)), S=float(m.group(2)), N=float(m.group(3)),
                omega=float(m.group(4)), dN=float(m.group(5)), dS=float(m.group(6)))
    return out


def one(args):
    og, d = args
    paml = os.path.join(d, og + ".paml")
    tmp = tempfile.mkdtemp(prefix=og + "_", dir=os.environ["DNDS_TMPDIR"])
    try:
        shutil.copy(paml, os.path.join(tmp, "in.paml"))
        with open(os.path.join(tmp, "codeml.ctl"), "w") as fh:
            fh.write(CTL.format(seq="in.paml"))
        r = subprocess.run([os.environ["CODEML"], "codeml.ctl"], cwd=tmp,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                           timeout=600)
        outp = os.path.join(tmp, "out")
        if r.returncode != 0 or not os.path.exists(outp):
            return og, None, "codeml_failed"
        pairs = parse_codeml(open(outp).read())
        ln, seqs = read_paml(paml)
        S, N, Sd, Nd = ng86(seqs["sugarcane"], seqs["purple"])
        row = {"orthogroup": og, "aln_codons": ln // 3,
               "ng86_S": S, "ng86_N": N, "ng86_Sd": Sd, "ng86_Nd": Nd,
               "gc3": gc3(seqs["sugarcane"])}
        for tag, (a, b) in (("sc_sb", ("sugarcane", "sorghum")),
                            ("pu_sb", ("purple", "sorghum")),
                            ("sc_pu", ("sugarcane", "purple"))):
            v = pairs.get(frozenset((a, b)))
            if v is None:
                return og, None, "missing_pair_" + tag
            for k in ("t", "S", "N", "dN", "dS", "omega"):
                row["%s_%s" % (tag, k)] = v[k]
        return og, row, "ok"
    except subprocess.TimeoutExpired:
        return og, None, "timeout"
    except Exception as e:
        return og, None, "error_%s" % type(e).__name__
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def gc3(seq):
    third = seq[2::3]
    if not third:
        return float("nan")
    return sum(1 for c in third if c in "GC") / len(third)


COLS = ["orthogroup", "aln_codons", "gc3",
        "sc_sb_t", "sc_sb_S", "sc_sb_N", "sc_sb_dN", "sc_sb_dS", "sc_sb_omega",
        "pu_sb_t", "pu_sb_S", "pu_sb_N", "pu_sb_dN", "pu_sb_dS", "pu_sb_omega",
        "sc_pu_t", "sc_pu_S", "sc_pu_N", "sc_pu_dN", "sc_pu_dS", "sc_pu_omega",
        "ng86_S", "ng86_N", "ng86_Sd", "ng86_Nd"]


def main():
    env = os.environ
    os.makedirs(env["DNDS_TMPDIR"], exist_ok=True)
    tasks = []
    for line in open(os.path.join(env["WORKDIR"], "og_list.txt")):
        og, d = line.rstrip("\n").split("\t")
        if os.path.exists(os.path.join(d, og + ".paml")):
            tasks.append((og, d))
    print("== codeml: %d alignments, %s workers" % (len(tasks), env["DNDS_THREADS"]))

    out = os.path.join(env["OUTDIR"], "dnds_raw.tsv")
    status = {}
    n_ok = 0
    with open(out, "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        with ProcessPoolExecutor(max_workers=int(env["DNDS_THREADS"])) as ex:
            for og, row, st in ex.map(one, tasks, chunksize=8):
                status[st] = status.get(st, 0) + 1
                if row is None:
                    continue
                fh.write("\t".join(
                    ("%.6g" % row[c]) if isinstance(row[c], float) else str(row[c])
                    for c in COLS) + "\n")
                n_ok += 1
    print("== status")
    for k, v in sorted(status.items(), key=lambda x: -x[1]):
        print("   %-22s %d" % (k, v))
    print("== wrote %d rows -> %s" % (n_ok, out))
    if n_ok == 0:
        sys.exit("FATAL: codeml produced nothing")


if __name__ == "__main__":
    main()
