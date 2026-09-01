# ============================================================================
# config.sh — dN/dS as a sequence-evolution layer on the comparative networks
#
# Question: does a gene whose co-expression neighbourhood is CONSERVED between
# the two networks also evolve under stronger purifying selection?
#
# The design fact that shapes everything here: R570 and LA purple are 95-99%
# identical at the protein level, so a per-gene omega between THEM is undefined
# or noise (dS is frequently exactly 0). So:
#
#   per-gene omega  -> measured against a Sorghum bicolor OUTGROUP (dS ~0.15)
#   the R570/purple pair -> kept, but only ever AGGREGATED (sum Nd,Sd,N,S
#                           within bins), never used per gene.
#
# Sourced by every step:  source config.sh
# NOTE: no `set -e` here (this file is *sourced*); step scripts set their own.
# ============================================================================

BASE="/dados04/jorge/comparative_saccharum"
CLEAN="${BASE}/new_clean"

OUTDIR="${CLEAN}/results/dnds"
WORKDIR="${OUTDIR}/work"
mkdir -p "$OUTDIR" "$WORKDIR"

# --- inputs: proteomes (the SAME files OrthoFinder was given) ---------------
SC_PEP="${BASE}/files/fix_orthofinder/proteins/sugarcane_one_transcript.fa"
PU_PEP="${BASE}/files/fix_orthofinder/proteins/one_transcript_purple_proteins.faa"

# --- inputs: CDS ------------------------------------------------------------
# sugarcane CDS lives OUTSIDE this repo, in the Phytozome reference tree. It is
# the primaryTranscriptOnly set: 194,593 records, exactly the count of
# sugarcane_one_transcript.fa, and it carries a `locus=` field that maps a
# transcript straight onto the gene id the proteome and the networks use.
CDS_sugarcane="/dados04/jorge/rnaseq_diatraea/reference_genomes/sugarcane/annotation/SofficinarumxspontaneumR570_771_v2.1.cds_primaryTranscriptOnly.fa.gz"
# purple CDS is in the repo: 255,706 records with `.tN` isoform suffixes, so it
# has MORE records than the 241,263-gene proteome. The isoform kept by the
# one-transcript file is not always `.t1`, so 01 recovers it by translating and
# matching the protein string rather than by assuming a suffix.
CDS_purple="${BASE}/files/fix_orthofinder/purple/purple_cds.fna"

# --- the outgroup: downloaded in full by 00, used from nowhere else ---------
SORGHUM_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/195/GCF_000003195.3_Sorghum_bicolor_NCBIv3"
SORGHUM_ACC="GCF_000003195.3_Sorghum_bicolor_NCBIv3"
SORGHUM_DIR="${BASE}/files/sorghum_ncbi"

# --- normalised sequence sets written by 01 (gene id as the ONLY header) ----
SEQDIR="${WORKDIR}/seqs"
SC_CDS_CLEAN="${SEQDIR}/sugarcane_cds.fna"
PU_CDS_CLEAN="${SEQDIR}/purple_cds.fna"
SB_CDS_CLEAN="${SEQDIR}/sorghum_cds.fna"
SC_PEP_CLEAN="${SEQDIR}/sugarcane_prot.faa"
PU_PEP_CLEAN="${SEQDIR}/purple_prot.faa"
SB_PEP_CLEAN="${SEQDIR}/sorghum_prot.faa"      # 1 protein per GeneID -> OrthoFinder input

# --- orthology --------------------------------------------------------------
# The EXISTING two-species run. This is the orthology every conservation number
# in the repo was computed from, and 02 must leave it byte-identical. It is used
# here ONLY as the consistency gate on the new triplets. NEVER repoint it.
ORTHOGROUPS="${BASE}/files/fix_orthofinder/proteins/OrthoFinder/Results_Jun04_2/Orthogroups/Orthogroups.tsv"

# The NEW three-species run (purple + sugarcane + sorghum), launched by 02 into
# its own tree. OrthoFinder 3.1.3 dropped the 2.x `-b prev -f new` add-species
# mode, so this is a clean rerun rather than an incremental one.
ORTHOFINDER_3SP="${BASE}/files/orthofinder_3sp"
OF3_IN="${ORTHOFINDER_3SP}/proteins3"
OF3_OUT="${ORTHOFINDER_3SP}/out"
OF3_ENV="othofinder3.1.3"     # MUST be `conda activate`d: the orthofinder
                              # launcher is #!/usr/bin/env python3, so calling
                              # it by absolute path picks up system python and
                              # dies on `import scipy`.
OF3_THREADS="${OF3_THREADS:-256}"
OF3_ANALYSIS="${OF3_ANALYSIS:-128}"

TRIPLETS="${OUTDIR}/triplets.tsv"

# --- alignment / codeml -----------------------------------------------------
MAFFT="${MAFFT:-/usr/bin/mafft}"
PAL2NAL="${PAL2NAL:-}"        # resolved by 00 (conda env, else bin/pal2nal.pl)
CODEML="${CODEML:-/home/genomics/miniconda3/envs/feelnc/bin/codeml}"
YN00="${YN00:-/home/genomics/miniconda3/envs/feelnc/bin/yn00}"
DNDS_ENV="dnds_env"

MIN_CODONS="${MIN_CODONS:-100}"   # after pal2nal -nogap
DNDS_THREADS="${DNDS_THREADS:-60}"
# codeml writes fixed filenames, so every orthogroup needs its own directory.
# Do not let tens of thousands of them land on a small system /tmp.
DNDS_TMPDIR="${DNDS_TMPDIR:-/dados04/jorge/tmp/dnds}"

# --- filters, applied in 06 and always REPORTED as counts, never silently ---
DNDS_DS_MIN="${DNDS_DS_MIN:-0.01}"    # below this omega is not estimable
DNDS_DS_MAX="${DNDS_DS_MAX:-2.0}"     # above this synonymous sites saturate
DNDS_OMEGA_FLAG="${DNDS_OMEGA_FLAG:-10}"   # flag as outlier, do not drop quietly

# --- network tables the dN/dS values are joined onto ------------------------
NODE_METRICS_sugarcane="${CLEAN}/results/sugarcane/network_sugarcane_node_metrics.tsv"
NODE_METRICS_purple="${CLEAN}/results/purple/network_purple_node_metrics.tsv"
HUBS_sugarcane="${CLEAN}/results/sugarcane/mcl_sugarcane_hub_genes.tsv"
HUBS_purple="${CLEAN}/results/purple/mcl_purple_hub_genes.tsv"
CONSDIR="${CLEAN}/results/conservation"
CONS_EDGES_sugarcane="${CONSDIR}/conserved_edges_sugarcane_to_purple_FULL.tsv"
CONS_EDGES_purple="${CONSDIR}/conserved_edges_purple_to_sugarcane_FULL.tsv"
TPM_sugarcane="${BASE}/run1/salmon/salmon.merged.gene_tpm.tsv"
TPM_purple="${BASE}/china/run2_onlyL/salmon/salmon.merged.gene_tpm.tsv"

# --- R interpreters (same envs the rest of the pipeline pins) ---------------
RSCRIPT_NET="${RSCRIPT_NET:-/home/genomics/miniconda3/envs/r_net_env/bin/Rscript}"
RSCRIPT_PLOT="${RSCRIPT_PLOT:-/home/genomics/miniconda3/envs/r_env/bin/Rscript}"

# --- subset switch: every step honours it, so the whole chain is testable ---
# DNDS_SUBSET=200 ./run_all.sh   -> 200 orthogroups end to end, minutes not hours
DNDS_SUBSET="${DNDS_SUBSET:-0}"

export BASE CLEAN OUTDIR WORKDIR SEQDIR \
       SC_PEP PU_PEP CDS_sugarcane CDS_purple \
       SORGHUM_URL SORGHUM_ACC SORGHUM_DIR \
       SC_CDS_CLEAN PU_CDS_CLEAN SB_CDS_CLEAN \
       SC_PEP_CLEAN PU_PEP_CLEAN SB_PEP_CLEAN \
       ORTHOGROUPS ORTHOFINDER_3SP OF3_IN OF3_OUT OF3_ENV OF3_THREADS OF3_ANALYSIS \
       TRIPLETS MAFFT PAL2NAL CODEML YN00 DNDS_ENV \
       MIN_CODONS DNDS_THREADS DNDS_TMPDIR \
       DNDS_DS_MIN DNDS_DS_MAX DNDS_OMEGA_FLAG \
       NODE_METRICS_sugarcane NODE_METRICS_purple HUBS_sugarcane HUBS_purple \
       CONSDIR CONS_EDGES_sugarcane CONS_EDGES_purple TPM_sugarcane TPM_purple \
       RSCRIPT_NET RSCRIPT_PLOT DNDS_SUBSET
