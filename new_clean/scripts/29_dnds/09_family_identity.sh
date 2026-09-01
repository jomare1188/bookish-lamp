#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh
python3 ./09_family_identity.py
