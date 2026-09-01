#!/usr/bin/env bash
# ============================================================================
# 02b_build_triplets.sh -- 1:1:1 orthogroups, gated on the two-species run
# ============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh
python3 ./02b_build_triplets.py
