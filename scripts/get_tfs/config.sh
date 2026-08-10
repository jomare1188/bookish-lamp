# ============================================================================
# config.sh — per-species configuration for the TF (TAP) identification pipeline
#
# Method: Riaño-Pachón 2007 / Pérez-Rodríguez 2010 (PlnTFDB-style rules).
# Sourced by every step script:  source config.sh <species>
#   <species> is one of:  sugarcane | purple
#
# NOTE: intentionally NO `set -e` here (this file is *sourced*); the step
# scripts set their own shell options.
# ============================================================================

SPECIES="${1:-}"

# --- shared resources (reused from the original GET_TFS run) ----------------
GETTFS_DIR="/dados04/jorge/comparative_saccharum/GET_TFS"
HMM_DB="${GETTFS_DIR}/db/TF.db.hmm"                        # 19,649 Pfam-A + PlnTFDB models, with GA cutoffs
RULES="${GETTFS_DIR}/mytfdb/RulesFull"                     # family-assignment rules
ASSIGN_PL="${GETTFS_DIR}/mytfdb/assign_family_membership.pl"
BASE="${GETTFS_DIR}/new"

# --- local parallelism (replaces the old SGE `#$ -t 1-84` array job) --------
# Keep  MAX_JOBS * CPU_PER_JOB  <=  nproc  (this box has 256 cores).
N_CHUNKS="${N_CHUNKS:-64}"        # split the proteome into this many pieces
MAX_JOBS="${MAX_JOBS:-32}"        # concurrent hmmsearch processes
CPU_PER_JOB="${CPU_PER_JOB:-4}"   # threads per hmmsearch

case "$SPECIES" in
  sugarcane)
    # R570 reference proteome (all isoforms)
    PROTEOME="/dados04/jorge/comparative_saccharum/files/fix_orthofinder/sugarcane/SofficinarumxspontaneumR570_771_v2.1.protein.fa"
    NODE_METRICS="/dados04/jorge/comparative_saccharum/files/sugarcane/network_sugarcane_node_metrics.tsv"
    # protein id  SoffiXsponR570.02Eg130400.1.p    -> gene  SoffiXsponR570.02Eg130400
    PROT_TO_GENE='s/\.[0-9]+\.p[0-9]*$//'
    # node id     SoffiXsponR570.02Eg130400.v2.1   -> gene  SoffiXsponR570.02Eg130400
    NODE_TO_GENE='s/\.v[0-9.]+$//'
    ;;
  purple)
    # LA purple reference proteome (one transcript per gene = same ids as the network nodes)
    PROTEOME="/dados04/jorge/comparative_saccharum/files/fix_orthofinder/purple/one_transcript_purple_proteins.faa"
    NODE_METRICS="/dados04/jorge/comparative_saccharum/files/purple/new/network_purple_node_metrics.tsv"
    # protein id  Soffic.02F0006060-2F   == gene (identity)
    PROT_TO_GENE='s/$//'
    # node id     Soffic.02F0006060-2F   == gene (identity)
    NODE_TO_GENE='s/$//'
    ;;
  *)
    echo "ERROR: pass a species -> sugarcane | purple" >&2
    return 1 2>/dev/null || exit 1
    ;;
esac

OUTDIR="${BASE}/results/${SPECIES}"
CHUNKDIR="${OUTDIR}/chunks"
DOMTBLDIR="${OUTDIR}/domtbl"

export SPECIES GETTFS_DIR HMM_DB RULES ASSIGN_PL BASE \
       N_CHUNKS MAX_JOBS CPU_PER_JOB \
       PROTEOME NODE_METRICS PROT_TO_GENE NODE_TO_GENE \
       OUTDIR CHUNKDIR DOMTBLDIR
