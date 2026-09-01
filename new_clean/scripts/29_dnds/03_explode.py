#!/usr/bin/env python3
"""
03_explode.py -- one protein FASTA and one CDS FASTA per orthogroup.

Done in a single pass with the three sequence sets held in memory, because the
alternative (a `seqkit grep` per orthogroup) is tens of thousands of process
launches over 300 MB files. Files are sharded 100 ways so no directory holds
tens of thousands of entries.

The protein and CDS files are written in the SAME order with the SAME headers,
which is what pal2nal needs to map an amino-acid alignment back onto codons.
"""
import os, sys


def read_fasta(path):
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


def load(path):
    return dict(read_fasta(path))


def shard(og):
    return "%02d" % (hash(og) % 100) if False else og[-2:]


def main():
    env = os.environ
    ogdir = os.path.join(env["WORKDIR"], "og")
    os.makedirs(ogdir, exist_ok=True)

    prot = {"sugarcane": load(env["SC_PEP_CLEAN"]),
            "purple":    load(env["PU_PEP_CLEAN"]),
            "sorghum":   load(env["SB_PEP_CLEAN"])}
    cds = {"sugarcane": load(env["SC_CDS_CLEAN"]),
           "purple":    load(env["PU_CDS_CLEAN"]),
           "sorghum":   load(env["SB_CDS_CLEAN"])}
    print("== loaded  sugarcane %d  purple %d  sorghum %d"
          % (len(prot["sugarcane"]), len(prot["purple"]), len(prot["sorghum"])))

    written, missing = 0, 0
    todo = open(os.path.join(env["WORKDIR"], "og_list.txt"), "w")
    with open(env["TRIPLETS"]) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            og, genes = f[0], {"sugarcane": f[1], "purple": f[2], "sorghum": f[3]}
            if any(genes[s] not in prot[s] or genes[s] not in cds[s] for s in genes):
                missing += 1
                continue
            d = os.path.join(ogdir, shard(og))
            os.makedirs(d, exist_ok=True)
            # headers are the SPECIES, not the gene id: codeml's pairwise output
            # is keyed by sequence name, and a fixed three-name vocabulary makes
            # the parse in 04 trivial and species-unambiguous. The gene ids stay
            # in triplets.tsv, which is the join key everywhere else.
            with open(os.path.join(d, og + ".faa"), "w") as pf, \
                 open(os.path.join(d, og + ".fna"), "w") as cf:
                for s in ("sugarcane", "purple", "sorghum"):
                    pf.write(">%s\n%s\n" % (s, prot[s][genes[s]]))
                    cf.write(">%s\n%s\n" % (s, cds[s][genes[s]]))
            todo.write("%s\t%s\n" % (og, d))
            written += 1
    todo.close()
    print("== exploded %d orthogroups (%d skipped, sequence missing)" % (written, missing))
    if written == 0:
        sys.exit("FATAL: nothing to align")


if __name__ == "__main__":
    main()
