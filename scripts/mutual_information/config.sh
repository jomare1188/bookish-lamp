# ============================================================================
# config.sh — shared paths and defaults for the mutual-information pipeline
#
# Sourced by run_all.sh and 00_export_vst.r (via environment variables).
# ============================================================================

BASE="/dados04/jorge/comparative_saccharum"
MIDIR="${BASE}/scripts/mutual_information"

# --- interpreters -----------------------------------------------------------
# torch 2.11 + cu130 lives in the `docling` env; it is the only env on this
# machine with a CUDA build. See README for how to make a dedicated env.
PYTORCH="/home/genomics/miniconda3/envs/docling/bin/python"
RSCRIPT_DESEQ="/home/genomics/miniconda3/envs/cor_env/bin/Rscript"

# --- inputs: the SAME DESeq2 objects the Pearson networks were built from ---
#
# NOTE ON THE PURPLE INPUT: the purple Pearson network came from china/run1 --
# 36 libraries, leaf and root -- narrowed to the 18 leaf libraries by
# Group1 == "L". That is NOT the same as china/run2_onlyL, which the module20
# and myb61 expression scripts use: the two quantifications disagree about which
# genes are flat, and the CV>=15 gene sets differ by 68 genes. Reproducing the
# Pearson gene set exactly requires run1 + the Group1=L filter, so that is what
# is used here. 00_export_vst.r verifies this and aborts if it ever drifts.
DDS_SUGARCANE="${BASE}/run1/salmon/deseq2_qc/deseq2.dds.RData"
COLS_SUGARCANE=""                 # all 48 libraries
DDS_PURPLE="${BASE}/china/run1/salmon/deseq2_qc/deseq2.dds.RData"
COLS_PURPLE="Group1=L"            # the 18 leaf libraries

# Existing Pearson matrices — only their HEADER is read, to verify that the
# gene set going into MI is byte-identical to the one that went into Pearson.
# Leave empty to skip the check.
PEARSON_SUGARCANE="${BASE}/files/sugarcane/matrix_sugarcane_pearson.tsv"
PEARSON_PURPLE="${BASE}/files/purple/new/matrix_purple_pearson.tsv"

# --- outputs ----------------------------------------------------------------
OUT_SUGARCANE="${BASE}/files/sugarcane/mi"
OUT_PURPLE="${BASE}/files/purple/new/mi"

# --- parameters -------------------------------------------------------------
MIN_CV=15          # identical to pearson_cor.r — do not change without reason
ESTIMATOR=ksg      # ksg | gcmi | xi
KSG_K=3            # Kraskov k; 3 is the usual small-sample choice
ALPHA=0.05         # BH FDR level
GPU_MEM_GB=5       # working-set budget on the GPU (A4500 has 20 GB); measured
                   # to be saturated past ~3 GB, so more buys nothing

# Effect-size floor, expressed as the Pearson |r| whose p-value at this n the
# floor should match. This is what makes the MI and Pearson networks equally
# stringent per edge, which in turn is what licenses taking their union.
#
# 0.8, NOT 0.7. build_edgelist.r's PEARSON_THRESHOLD <- 0.7 is only a pre-filter
# for the intermediate edgelist_*_pearson.tsv; the networks actually analysed
# come from general_stats.r, whose PEARSON_MIN <- 0.8 is the real cut. Verified
# against the data: network_{sugarcane,purple}_filtered_edges.tsv both bottom
# out at |r| = 0.8 (p = 9.06e-12 and 6.74e-05 respectively).
MATCH_PEARSON=0.8

# Per-pair p below which an edge is written to a shard. Also the range over
# which BH stays exact (up to CAND_P * n_tests / ALPHA rejections). At 1e-3 the
# sugarcane sweep stores ~4.5e8 candidates, about 5.4 GB of shards.
CAND_P=1e-3

# Ceiling on the auto-sized permutation null. At n=48 the |r|=0.8 floor sits at
# p = 9.0e-12, far past what any affordable permutation count resolves, so the
# GPD tail does the work; 2e8 permutations pin the fit well enough that the
# measured seed-to-seed spread of the floor is ~0.01 nats.
MAX_PERM=200000000

# Host RAM the FDR step may use (8 bytes per candidate edge). Purple is the
# demanding one: ~1.9e9 candidates, about 15 GB. The machine has ~440 GB free.
RAM_LIMIT_GB=64

export BASE MIDIR PYTORCH RSCRIPT_DESEQ MIN_CV
