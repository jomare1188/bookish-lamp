#!/usr/bin/env bash
# ============================================================================
# 00_setup.sh — tooling + the Sorghum bicolor outgroup, downloaded in full
#
# Everything sorghum in this analysis comes from ONE NCBI release fetched here.
# Nothing is picked up from elsewhere on the machine, so the outgroup has a
# single, checkable provenance.
#
# Re-runnable: files already present whose MD5 matches are not re-fetched.
# ============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

mkdir -p "$SORGHUM_DIR" "$SEQDIR" "$DNDS_TMPDIR"

# ---------------------------------------------------------------------------
# 1. pal2nal — the one tool that is genuinely missing on this machine
# ---------------------------------------------------------------------------
echo "== pal2nal"
resolve_pal2nal() {
  # already configured and working?
  if [ -n "${PAL2NAL:-}" ] && [ -x "${PAL2NAL}" ]; then echo "$PAL2NAL"; return; fi
  # a dedicated conda env, if it exists
  local envp="/home/genomics/miniconda3/envs/${DNDS_ENV}/bin/pal2nal.pl"
  if [ -x "$envp" ]; then echo "$envp"; return; fi
  # anything on PATH
  if command -v pal2nal.pl >/dev/null 2>&1; then command -v pal2nal.pl; return; fi
  # the vendored copy
  if [ -x "./bin/pal2nal.pl" ]; then echo "$(pwd)/bin/pal2nal.pl"; return; fi
  echo ""
}

P2N="$(resolve_pal2nal)"
if [ -z "$P2N" ]; then
  echo "   not found -- creating conda env '${DNDS_ENV}'"
  if conda create -y -n "$DNDS_ENV" -c bioconda -c conda-forge pal2nal paml mafft; then
    P2N="/home/genomics/miniconda3/envs/${DNDS_ENV}/bin/pal2nal.pl"
  else
    echo "   conda route failed -- fetching the standalone script into ./bin/"
    wget -q -O bin/pal2nal.tar.gz \
      "http://www.bork.embl.de/pal2nal/distribution/pal2nal.v14.tar.gz"
    tar -xzf bin/pal2nal.tar.gz -C bin --strip-components=1 pal2nal.v14/pal2nal.pl
    chmod +x bin/pal2nal.pl
    P2N="$(pwd)/bin/pal2nal.pl"
  fi
fi
[ -x "$P2N" ] || { echo "FATAL: pal2nal.pl still unavailable at '$P2N'"; exit 1; }
echo "   pal2nal.pl -> $P2N"
echo "$P2N" > "${WORKDIR}/pal2nal.path"

# ---------------------------------------------------------------------------
# 2. the rest of the toolchain must already be here -- fail loudly if not
# ---------------------------------------------------------------------------
echo "== toolchain"
for t in "$MAFFT" "$CODEML" "$YN00"; do
  [ -x "$t" ] || { echo "FATAL: missing $t"; exit 1; }
  echo "   ok  $t"
done
command -v seqkit >/dev/null || { echo "FATAL: seqkit not on PATH"; exit 1; }
command -v parallel >/dev/null || echo "   WARNING: GNU parallel absent; 03/04 fall back to xargs -P"

# ---------------------------------------------------------------------------
# 3. the outgroup: genome + annotation + sequences, one release
# ---------------------------------------------------------------------------
echo "== sorghum outgroup: ${SORGHUM_ACC}"
FILES=(
  "${SORGHUM_ACC}_genomic.fna.gz"           # assembly
  "${SORGHUM_ACC}_genomic.gff.gz"           # annotation
  "${SORGHUM_ACC}_genomic.gtf.gz"
  "${SORGHUM_ACC}_cds_from_genomic.fna.gz"  # used by 01
  "${SORGHUM_ACC}_protein.faa.gz"           # used by 01 -> 02
  "${SORGHUM_ACC}_translated_cds.faa.gz"    # frame cross-check
  "${SORGHUM_ACC}_assembly_report.txt"
)

wget -q -O "${SORGHUM_DIR}/md5checksums.txt" "${SORGHUM_URL}/md5checksums.txt"

for f in "${FILES[@]}"; do
  want="$(awk -v f="./$f" '$2==f {print $1}' "${SORGHUM_DIR}/md5checksums.txt")"
  if [ -z "$want" ]; then echo "FATAL: $f absent from md5checksums.txt"; exit 1; fi
  if [ -s "${SORGHUM_DIR}/${f}" ]; then
    have="$(md5sum "${SORGHUM_DIR}/${f}" | cut -d' ' -f1)"
    if [ "$have" = "$want" ]; then echo "   cached  $f"; continue; fi
    echo "   stale   $f -- refetching"
  fi
  echo "   fetch   $f"
  wget -q -c -O "${SORGHUM_DIR}/${f}" "${SORGHUM_URL}/${f}"
  have="$(md5sum "${SORGHUM_DIR}/${f}" | cut -d' ' -f1)"
  [ "$have" = "$want" ] || { echo "FATAL: MD5 mismatch on $f ($have != $want)"; exit 1; }
done
echo "   all files match md5checksums.txt"

cat > "${SORGHUM_DIR}/PROVENANCE.txt" <<PROV
Sorghum bicolor outgroup for the dN/dS analysis.

accession   ${SORGHUM_ACC}
source      ${SORGHUM_URL}
downloaded  $(date -u +%Y-%m-%dT%H:%M:%SZ) UTC
by          new_clean/scripts/29_dnds/00_setup.sh
verified    every file against md5checksums.txt from the same directory

This is the ONLY sorghum data used anywhere in the analysis. Other sorghum
files exist elsewhere on this machine (Phytozome v3.1.1 CDS, an Ensembl
proteome under the MYB61 readout); none of them are read by this stage.
PROV

echo
echo "00 done."
echo "  pal2nal : $P2N"
echo "  sorghum : ${SORGHUM_DIR} ($(ls -1 "${SORGHUM_DIR}" | wc -l) files)"
