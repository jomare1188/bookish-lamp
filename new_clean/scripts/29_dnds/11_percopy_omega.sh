#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh
export PAL2NAL="$(cat "${WORKDIR}/pal2nal.path")"
python3 ./11_percopy_omega.py
