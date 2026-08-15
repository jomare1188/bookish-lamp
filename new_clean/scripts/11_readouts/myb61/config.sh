# ============================================================================
# config.sh — shared configuration for the MYB61 orthologue search (H1)
#
# Goal: locate the sugarcane orthologue(s) of MYB61 in BOTH of our reference
# proteomes by homology to a curated reference anchor set, because the gene id
# published by Kiet et al. (`Soff.09G0002230-3D`) cannot be reconciled with our
# LA purple annotation (see README.md, "Why not the published id").
#
# Sourced by every step:  source config.sh
# NOTE: no `set -e` here (this file is *sourced*); step scripts set their own.
# ============================================================================

BASE="/dados04/jorge/comparative_saccharum"

# --- cache vs output --------------------------------------------------------
# CACHE_DIR is the ORIGINAL GET_TFS tree, read-only here. Steps 01-06 are sequence/orthology/tree
# work and are network-independent, so their output (work/, refs/, phylogeny/,
# *_MYB61_candidates.tsv) is reused rather than recomputed.
# OUTDIR is under new_clean/, so a re-run never overwrites the previous results.
CACHE_DIR="${BASE}/GET_TFS/new/results/myb61"
OUTDIR="/dados04/jorge/comparative_saccharum/new_clean/results/readouts/myb61"
mkdir -p "$OUTDIR"
CACHE_REFDIR="${CACHE_DIR}/refs"
REFDIR="${OUTDIR}/refs"          # downloaded reference proteomes
WORKDIR="${OUTDIR}/work"         # diamond dbs + raw hit tables
mkdir -p "$REFDIR" "$WORKDIR"

# --- the anchor -------------------------------------------------------------
# AtMYB61 = AT1G09540 (R2R3-MYB; stomatal aperture, lignin, seed-coat mucilage)
AT_ANCHOR="AT1G09540"

# --- Ensembl Plants reference proteomes (open FTP, stable filenames) ---------
ENSEMBL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/fasta"
AT_URL="${ENSEMBL}/arabidopsis_thaliana/pep/Arabidopsis_thaliana.TAIR10.pep.all.fa.gz"
SB_URL="${ENSEMBL}/sorghum_bicolor/pep/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.pep.all.fa.gz"
OS_URL="${ENSEMBL}/oryza_sativa/pep/Oryza_sativa.IRGSP-1.0.pep.all.fa.gz"

AT_PEP="${REFDIR}/Arabidopsis_thaliana.TAIR10.pep.primary.fa"   # 1 protein / gene
SB_PEP="${REFDIR}/Sorghum_bicolor.pep.primary.fa"
OS_PEP="${REFDIR}/Oryza_sativa.pep.primary.fa"
QUERY="${OUTDIR}/MYB61_query.faa"                               # At + Sb + Os anchors

# --- our two proteomes (same files as the OrthoFinder bridge & GET_TFS) -----
SC_PEP="${BASE}/files/fix_orthofinder/sugarcane/SofficinarumxspontaneumR570_771_v2.1.protein.fa"
PU_PEP="${BASE}/files/fix_orthofinder/purple/one_transcript_purple_proteins.faa"

# protein id -> gene id (identical to scripts/get_tfs/config.sh)
SC_PROT_TO_GENE='s/\.[0-9]+\.p[0-9]*$//'      # SoffiXsponR570.02Eg130400.1.p -> ...130400
PU_PROT_TO_GENE='s/$//'                        # Soffic.02F0006060-2F == gene

# --- our own TF calls (proteome-wide, from the GET_TFS rerun) ---------------
SC_TF="/dados04/jorge/comparative_saccharum/new_clean/results/readouts/get_tfs/sugarcane/TF_no_Orphans.tsv"
PU_TF="/dados04/jorge/comparative_saccharum/new_clean/results/readouts/get_tfs/purple/TF_no_Orphans.tsv"

# --- OrthoFinder bridge -----------------------------------------------------
ORTHOGROUPS="${BASE}/files/fix_orthofinder/proteins/OrthoFinder/Results_Jun04_2/Orthogroups/Orthogroups.tsv"

# --- the locus Kiet et al. published, kept as a labelled cross-check --------
KIET_LOCUS="09G0002230"

# --- R interpreters ---------------------------------------------------------
# `ape` (tree handling, step 06) lives in the comparative_network env;
# ggplot2/data.table/patchwork (step 07 figure) live in r_env.
RSCRIPT_APE="${RSCRIPT_APE:-/home/genomics/miniconda3/envs/comparative_network/bin/Rscript}"
RSCRIPT_PLOT="${RSCRIPT_PLOT:-/home/genomics/miniconda3/envs/r_env/bin/Rscript}"

# --- search thresholds ------------------------------------------------------
# Permissive on the forward search (we want the whole MYB61 neighbourhood),
# strict on the disambiguation (reciprocal best hit must return the anchor).
FWD_EVALUE="${FWD_EVALUE:-1e-10}"
FWD_MAX_TARGETS="${FWD_MAX_TARGETS:-500}"
MIN_QCOV="${MIN_QCOV:-50}"       # % of the reference query covered
THREADS="${THREADS:-32}"

export BASE OUTDIR REFDIR WORKDIR AT_ANCHOR \
       AT_URL SB_URL OS_URL AT_PEP SB_PEP OS_PEP QUERY \
       SC_PEP PU_PEP SC_PROT_TO_GENE PU_PROT_TO_GENE SC_TF PU_TF \
       ORTHOGROUPS KIET_LOCUS FWD_EVALUE FWD_MAX_TARGETS MIN_QCOV THREADS \
       RSCRIPT_APE RSCRIPT_PLOT
