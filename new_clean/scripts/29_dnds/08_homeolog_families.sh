#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh
python3 ./08_homeolog_families.py
