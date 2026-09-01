#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh
python3 ./04b_yn00_check.py
