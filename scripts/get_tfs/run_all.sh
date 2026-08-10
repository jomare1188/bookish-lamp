#!/bin/bash
# ============================================================================
# run_all.sh — end-to-end TF identification for one species.
# Usage:  ./run_all.sh <species>            (sugarcane | purple)
#         N_CHUNKS=32 MAX_JOBS=16 ./run_all.sh purple   # override parallelism
# ============================================================================
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
sp="${1:?usage: run_all.sh <species>   (sugarcane | purple)}"

echo "############ TF pipeline :: ${sp} ############"
"${here}/00_check_requirements.sh" "$sp"
"${here}/01_split_proteome.sh"     "$sp"
"${here}/02_run_hmmsearch.sh"      "$sp"    # the long step
"${here}/03_assign_families.sh"    "$sp"
"${here}/04_postprocess.sh"        "$sp"
echo "############ DONE :: ${sp} ############"
