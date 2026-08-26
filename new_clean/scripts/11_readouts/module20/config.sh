# ============================================================================
# config.sh — shared configuration for the Muñoz Module-20 cross-network test
#
# Muñoz-Perez et al. define Module 20 (M20) on THEIR OWN network, built from a
# de novo pan-transcriptome. M20 is ~75% MYB/MYB-related TAPs and is their
# central nitrogen claim. This pipeline asks whether those genes keep the same
# characteristics — MYB identity, high degree, nitrogen responsiveness — in OUR
# two reference-based networks (R570 = sugarcane, LA purple = purple).
#
# NOTE the circularity trap: M20 was DEFINED as MYB-rich and high-betweenness in
# their network, so "are they MYBs with high degree?" is not a real test. The
# non-circular questions are (1) do they map, (2) do they stand out against the
# MYB background IN OUR networks, (3) are they N-responsive in our own
# re-quantification, (4) does any of it survive in purple.
#
# Sourced by every step:  source config.sh
# ============================================================================

BASE="/dados04/jorge/comparative_saccharum"
MUNOZ="${BASE}/raw_sugarcane/transcriptome_munoz"
M20_CSV="${MUNOZ}/module20.csv"
M20_TRANSCRIPTOME="${MUNOZ}/Projects_Sugarcane_NitrogenResponsiveGenotypes_Salmon_All_masked__100p__unmasked.fasta"
# Muñoz PROTEOME (TransDecoder ORFs, ids <transcript>.pN). Byte-identical to
# GET_TFS/db/sugar_cane.pep — i.e. the very proteome the original GET_TFS TF
# run used, so Muñoz's own Family column comes from the same PlnTFDB rule
# method this project uses. Mapping is done protein-vs-protein against it.
M20_PROTEOME="${MUNOZ}/sugar_cane.pep"

# Steps 01-02b are sequence-only and cached in the original tree; 03-05 are
# network-dependent and write under new_clean/.
CACHE_DIR="${BASE}/GET_TFS/new/results/module20"
OUTDIR="/dados04/jorge/comparative_saccharum/new_clean/results/readouts/module20"
mkdir -p "$OUTDIR"
WORKDIR="${OUTDIR}/work"
mkdir -p "$OUTDIR" "$WORKDIR"

M20_FNA="${OUTDIR}/module20.fna"          # the 12 contigs (nucleotide, legacy/QC)
M20_FAA="${OUTDIR}/module20.faa"          # the Module-20 ORFs (protein) - primary input

# --- DIAMOND databases (built by the myb61 pipeline; rebuilt here if absent) -
MYB61_WORK="${BASE}/GET_TFS/new/results/myb61/work"
DB_SC="${MYB61_WORK}/sugarcane"                                   # R570 proteome
DB_PU="${MYB61_WORK}/purple"                                      # LA purple proteome
DB_AT="${MYB61_WORK}/Arabidopsis_thaliana.TAIR10.pep.primary"     # reciprocal target
AT_PEP="${BASE}/GET_TFS/new/results/myb61/refs/Arabidopsis_thaliana.TAIR10.pep.primary.fa"
SC_CLEAN="${MYB61_WORK}/sugarcane.clean.faa"
PU_CLEAN="${MYB61_WORK}/purple.clean.faa"

# --- our networks -----------------------------------------------------------
SC_NODES="${BASE}/new_clean/results/sugarcane/network_sugarcane_node_metrics.tsv"
PU_NODES="${BASE}/new_clean/results/purple/network_purple_node_metrics.tsv"
SC_NODE_STRIP='\.v[0-9.]+$'      # node ids carry .v2.1; gene ids do not
PU_NODE_STRIP=''

# NOTE: 03_network_readout.r does not read these two -- it carries its own
# defaults (overridable via CLEAN_M20_MODS_SUGARCANE / _PURPLE). Change both
# together or they will disagree about which clustering the readout used.
SC_MODULES="${BASE}/new_clean/results/sugarcane/mcl_sugarcane_membership.tsv"
PU_MODULES="${BASE}/new_clean/results/purple/mcl_purple_membership.tsv"

SC_TF_NET="/dados04/jorge/comparative_saccharum/new_clean/results/readouts/get_tfs/sugarcane/TF_in_network.tsv"
PU_TF_NET="/dados04/jorge/comparative_saccharum/new_clean/results/readouts/get_tfs/purple/TF_in_network.tsv"

SC_CONSERVED="${BASE}/new_clean/results/conservation/conserved_genes_sugarcane_FULL.txt"
PU_CONSERVED="${BASE}/new_clean/results/conservation/conserved_genes_purple_FULL.txt"

# --- expression (design-aware nitrogen tests, step 04) ----------------------
# sugarcane: 48 libs = 2 genotypes x 2 N levels x 3 leaf segments x reps.
#   The three leaf segments (Leaf Apex / Leaf Base / Leaf) are treated as ONE
#   leaf tissue, per the project's decision; segment is kept only as a
#   sensitivity covariate.
SC_TPM="${BASE}/run1/salmon/salmon.merged.gene_tpm.tsv"
SC_META="${BASE}/samplesheet.csv"
# purple: 18 leaf libs = 2 genotypes x (0N,2N,6N) x 3 reps
PU_TPM="${BASE}/china/run2_onlyL/salmon/salmon.merged.gene_tpm.tsv"
PU_META="${BASE}/china/samplesheet_china.csv"

ORTHOGROUPS="${BASE}/files/fix_orthofinder/proteins/OrthoFinder/Results_Jun04_2/Orthogroups/Orthogroups.tsv"

# --- our own TF identification, applied to the Module-20 proteins (step 02b) -
HMM_DB="${BASE}/GET_TFS/db/TF.db.hmm"
RULES="${BASE}/GET_TFS/mytfdb/RulesFull"
ASSIGN_PL="${BASE}/GET_TFS/mytfdb/assign_family_membership.pl"

# --- search thresholds ------------------------------------------------------
# Permissive forward search (polyploid: one transcript legitimately hits many
# haplotype copies); specificity comes from the reciprocal-Arabidopsis filter.
EVALUE="${EVALUE:-1e-10}"
MAX_TARGETS="${MAX_TARGETS:-500}"
MIN_SCOV="${MIN_SCOV:-50}"      # % of the subject protein covered
THREADS="${THREADS:-32}"

RSCRIPT_PLOT="${RSCRIPT_PLOT:-/home/genomics/miniconda3/envs/r_env/bin/Rscript}"

export BASE MUNOZ M20_CSV M20_TRANSCRIPTOME M20_PROTEOME OUTDIR WORKDIR M20_FNA M20_FAA \
       HMM_DB RULES ASSIGN_PL \
       MYB61_WORK DB_SC DB_PU DB_AT AT_PEP SC_CLEAN PU_CLEAN \
       SC_NODES PU_NODES SC_NODE_STRIP PU_NODE_STRIP SC_MODULES PU_MODULES \
       SC_TF_NET PU_TF_NET SC_CONSERVED PU_CONSERVED \
       SC_TPM SC_META PU_TPM PU_META ORTHOGROUPS \
       EVALUE MAX_TARGETS MIN_SCOV THREADS RSCRIPT_PLOT
