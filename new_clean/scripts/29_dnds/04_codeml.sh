#!/usr/bin/env bash
# ============================================================================
# 04_codeml.sh -- pairwise dN/dS (codeml -2) + exact NG86 counts
#
#   *** the full pass is ~30 min on 60 cores. LAUNCH IT YOURSELF. ***
# ============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh
python3 ./04_codeml.py
