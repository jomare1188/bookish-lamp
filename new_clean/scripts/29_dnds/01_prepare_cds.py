#!/usr/bin/env python3
"""
01_prepare_cds.py -- normalise CDS and protein FASTAs to ONE record per gene,
headered with the gene id the networks and OrthoFinder use.

Three species, three different id problems:

  sugarcane  Phytozome primaryTranscriptOnly. Transcript-headered, but every
             record carries `locus=`, which IS the proteome/network gene id.
             Counts already match (194,593), so this is a pure rename.

  purple     255,706 CDS records with `.tN` isoform suffixes against a
             241,263-gene proteome. The isoform the one-transcript proteome
             kept is NOT always `.t1`, so the right isoform is recovered by
             translating each candidate and matching the protein string.

  sorghum    RefSeq, many isoforms per gene. Keyed on `[db_xref=GeneID:N]`,
             joined to the protein FASTA through `[protein_id=XP_...]`, and
             reduced to the longest protein per gene -- which must happen HERE,
             before OrthoFinder, so the three proteomes share one convention.

Every species passes the same gate: translate(CDS) must reproduce the protein.
pal2nal fails silently on a frame mismatch, so anything that does not translate
cleanly is dropped here and counted, never carried forward.
"""
import gzip, os, re, sys
from collections import defaultdict

CODONS = {}
_B = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, b1 in enumerate(_B):
    for j, b2 in enumerate(_B):
        for k, b3 in enumerate(_B):
            CODONS[b1 + b2 + b3] = _AA[i * 16 + j * 4 + k]


def translate(seq):
    seq = seq.upper().replace("U", "T")
    return "".join(CODONS.get(seq[i:i + 3], "X") for i in range(0, len(seq) - len(seq) % 3, 3))


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_fasta(path):
    name, buf = None, []
    with opener(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name, buf = line[1:].rstrip("\n"), []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


def prot_match(cds, prot):
    """Does this CDS translate to this protein? Returns (ok, identity)."""
    t = translate(cds).rstrip("*")
    p = prot.rstrip("*")
    if not p:
        return False, 0.0
    if t == p:
        return True, 1.0
    n = min(len(t), len(p))
    if n == 0 or abs(len(t) - len(p)) > 1:
        return False, 0.0
    same = sum(1 for a, b in zip(t[:n], p[:n]) if a == b or a == "X" or b == "X")
    ident = same / max(len(t), len(p))
    return ident >= 0.99, ident


def write_fasta(path, records, width=60):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(">%s\n" % name)
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


def qc_report(tag, kept, dropped, reasons, qc_path):
    with open(qc_path, "w") as fh:
        fh.write("gene\treason\tdetail\n")
        for g, r, d in reasons:
            fh.write("%s\t%s\t%s\n" % (g, r, d))
    total = kept + dropped
    pct = 100.0 * kept / total if total else 0.0
    print("   %-10s kept %7d / %7d  (%.2f%%)  dropped %d -> %s"
          % (tag, kept, total, pct, dropped, os.path.basename(qc_path)))
    return pct


def do_sugarcane(cds_path, pep_path, out_cds, out_pep, qc_path):
    prots = {n.split()[0]: s for n, s in read_fasta(pep_path)}
    locus_re = re.compile(r"locus=(\S+)")
    kept_c, kept_p, reasons, seen = [], [], [], set()
    dropped = 0
    for name, seq in read_fasta(cds_path):
        m = locus_re.search(name)
        if not m:
            dropped += 1
            reasons.append((name.split()[0], "no_locus_field", ""))
            continue
        gene = m.group(1)
        if gene in seen:
            dropped += 1
            reasons.append((gene, "duplicate_locus", ""))
            continue
        prot = prots.get(gene)
        if prot is None:
            dropped += 1
            reasons.append((gene, "no_protein", ""))
            continue
        ok, ident = prot_match(seq, prot)
        if not ok:
            dropped += 1
            reasons.append((gene, "translation_mismatch", "%.3f" % ident))
            continue
        seen.add(gene)
        kept_c.append((gene, seq))
        kept_p.append((gene, prot.rstrip("*")))
    write_fasta(out_cds, kept_c)
    write_fasta(out_pep, kept_p)
    return qc_report("sugarcane", len(kept_c), dropped, reasons, qc_path)


def do_purple(cds_path, pep_path, out_cds, out_pep, qc_path):
    prots = {n.split()[0]: s for n, s in read_fasta(pep_path)}
    isoforms = defaultdict(list)
    for name, seq in read_fasta(cds_path):
        tid = name.split()[0]
        gene = re.sub(r"\.t\d+$", "", tid)
        isoforms[gene].append((tid, seq))
    kept_c, kept_p, reasons = [], [], []
    dropped = 0
    for gene, cands in isoforms.items():
        prot = prots.get(gene)
        if prot is None:
            # gene absent from the one-transcript proteome: not an error, that
            # proteome is the OrthoFinder input and defines the gene universe
            continue
        pick, best = None, (0.0, None)
        for tid, seq in cands:
            ok, ident = prot_match(seq, prot)
            if ok:
                pick = (tid, seq)
                break
            if ident > best[0]:
                best = (ident, tid)
        if pick is None:
            dropped += 1
            reasons.append((gene, "no_isoform_matches_protein",
                            "best=%s ident=%.3f n_iso=%d" % (best[1], best[0], len(cands))))
            continue
        kept_c.append((gene, pick[1]))
        kept_p.append((gene, prot.rstrip("*")))
    write_fasta(out_cds, kept_c)
    write_fasta(out_pep, kept_p)
    return qc_report("purple", len(kept_c), dropped, reasons, qc_path)


def do_sorghum(cds_path, pep_path, out_cds, out_pep, qc_path):
    prots = {n.split()[0]: s for n, s in read_fasta(pep_path)}
    gid_re = re.compile(r"\[db_xref=GeneID:(\d+)\]")
    pid_re = re.compile(r"\[protein_id=([^\]]+)\]")
    sym_re = re.compile(r"\[gene=([^\]]+)\]")
    pseudo_re = re.compile(r"\[pseudo=true\]")

    by_gene = defaultdict(list)
    sym_of = {}
    for name, seq in read_fasta(cds_path):
        if pseudo_re.search(name):
            continue
        g, p = gid_re.search(name), pid_re.search(name)
        if not g or not p:
            continue
        gid, pid = g.group(1), p.group(1)
        s = sym_re.search(name)
        if s:
            sym_of.setdefault(gid, s.group(1))
        by_gene[gid].append((pid, seq))

    # a gene symbol is only usable as an id if it names exactly one GeneID
    sym_count = defaultdict(int)
    for gid, s in sym_of.items():
        sym_count[s] += 1

    kept_c, kept_p, reasons = [], [], []
    dropped = 0
    for gid, cands in by_gene.items():
        best = None
        for pid, seq in cands:
            prot = prots.get(pid)
            if prot is None:
                continue
            ok, ident = prot_match(seq, prot)
            if not ok:
                continue
            if best is None or len(prot) > len(best[2]):
                best = (pid, seq, prot)
        if best is None:
            dropped += 1
            reasons.append((gid, "no_isoform_translates", "n_iso=%d" % len(cands)))
            continue
        sym = sym_of.get(gid)
        gene = sym if (sym and sym_count[sym] == 1) else "GeneID%s" % gid
        kept_c.append((gene, best[1]))
        kept_p.append((gene, best[2].rstrip("*")))
    write_fasta(out_cds, kept_c)
    write_fasta(out_pep, kept_p)
    return qc_report("sorghum", len(kept_c), dropped, reasons, qc_path)


if __name__ == "__main__":
    env = os.environ
    seqdir = env["SEQDIR"]
    os.makedirs(seqdir, exist_ok=True)
    sorg = env["SORGHUM_DIR"]
    acc = env["SORGHUM_ACC"]

    print("== normalising CDS + protein to one record per gene")
    pcts = []
    pcts.append(do_sugarcane(env["CDS_sugarcane"], env["SC_PEP"],
                             env["SC_CDS_CLEAN"], env["SC_PEP_CLEAN"],
                             os.path.join(seqdir, "cds_qc_sugarcane.tsv")))
    pcts.append(do_purple(env["CDS_purple"], env["PU_PEP"],
                          env["PU_CDS_CLEAN"], env["PU_PEP_CLEAN"],
                          os.path.join(seqdir, "cds_qc_purple.tsv")))
    pcts.append(do_sorghum(os.path.join(sorg, acc + "_cds_from_genomic.fna.gz"),
                           os.path.join(sorg, acc + "_protein.faa.gz"),
                           env["SB_CDS_CLEAN"], env["SB_PEP_CLEAN"],
                           os.path.join(seqdir, "cds_qc_sorghum.tsv")))
    if min(pcts) < 95.0:
        print("\nWARNING: a species kept < 95%% of its genes -- inspect cds_qc_*.tsv "
              "before trusting anything downstream.", file=sys.stderr)
