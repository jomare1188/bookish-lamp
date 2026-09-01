#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh
python3 ./10_kmer_uniqueness.py
