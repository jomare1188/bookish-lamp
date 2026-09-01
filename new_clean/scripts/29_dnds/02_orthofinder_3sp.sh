#!/usr/bin/env bash
# ============================================================================
# 02_orthofinder_3sp.sh -- add sorghum to the orthology, by inference
#
#   *** THIS IS THE LONG STEP. ~3-4 h at -t 256. LAUNCH IT YOURSELF. ***
#
# OrthoFinder 3.1.3 removed the 2.x `-b <prev> -f <new>` add-species mode
# (`--assign/--core` replaced it, and that is mutually exclusive with -f), so
# sorghum joins through a clean three-proteome run into its OWN output tree.
#
# WHAT THIS MUST NOT DO: touch the existing two-species run. config.sh's
# ORTHOGROUPS points at Results_Jun04_2, and every conservation number in the
# repo -- the 2.5x fold over null, the responsive-ortholog-pair counts, all of
# it -- was computed from that file. It is checksummed before and after here,
# and the run is sent elsewhere with -o.
#
# Inputs are the CLEANED proteomes from step 01, not the raw ones: those are
# exactly the genes with a verified, in-frame CDS, so every single-copy
# orthogroup this run finds is guaranteed usable downstream. The cost is 519
# purple genes (0.2%) that the raw proteome had and whose CDS does not
# translate -- they could never have yielded a dN/dS value anyway.
# ============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

[ -s "$SC_PEP_CLEAN" ] && [ -s "$PU_PEP_CLEAN" ] && [ -s "$SB_PEP_CLEAN" ] || {
  echo "FATAL: run ./01_prepare_cds.sh first"; exit 1; }

# --- guard the two-species run ---------------------------------------------
GUARD="${WORKDIR}/orthogroups_2sp.md5"
md5sum "$ORTHOGROUPS" > "${GUARD}.before"
echo "== guarding two-species orthology"
cat "${GUARD}.before" | sed 's/^/   /'

# --- inputs -----------------------------------------------------------------
rm -rf "$OF3_IN"
mkdir -p "$OF3_IN"
ln -sf "$PU_PEP_CLEAN" "${OF3_IN}/purple.faa"
ln -sf "$SC_PEP_CLEAN" "${OF3_IN}/sugarcane.faa"
ln -sf "$SB_PEP_CLEAN" "${OF3_IN}/sorghum.faa"
echo "== inputs"
for f in "${OF3_IN}"/*.faa; do
  printf "   %-14s %8d proteins\n" "$(basename "$f")" "$(grep -c '^>' "$f")"
done

# OrthoFinder refuses to write into an existing -o directory. Do NOT silently
# delete one: it may hold a finished 3-4 h run. Refuse and make the caller say so.
if [ -e "$OF3_OUT" ]; then
  if [ "${FORCE:-0}" = "1" ]; then
    echo "== FORCE=1: discarding the previous run at ${OF3_OUT}"
    rm -rf "$OF3_OUT"
  else
    echo "FATAL: ${OF3_OUT} already exists."
    echo "       That may be a completed 3-4 h run. To replace it deliberately:"
    echo "         FORCE=1 ./02_orthofinder_3sp.sh"
    echo "       To keep it and just rebuild the triplets: ./02b_build_triplets.sh"
    exit 1
  fi
fi

# --- run --------------------------------------------------------------------
# Same flags as files/fix_orthofinder/run_orthofinder.sh, so the two runs are
# methodologically comparable. The env MUST be activated: the orthofinder
# launcher is #!/usr/bin/env python3, so invoking it by absolute path picks up
# system python and dies on `import scipy`.
echo "== OrthoFinder ${OF3_ENV}: 3 species, -t ${OF3_THREADS} -a ${OF3_ANALYSIS}"
echo "   started $(date -u +%Y-%m-%dT%H:%M:%SZ)"
# This env's activate hook (aster_activate.sh) dereferences LD_LIBRARY_PATH
# unconditionally, so it aborts under `set -u`. Give it a default and drop -u
# across the activation only -- the usual conda idiom -- rather than weakening
# the whole script.
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
set +u
# shellcheck disable=SC1091
source /home/genomics/miniconda3/etc/profile.d/conda.sh
conda activate "$OF3_ENV"
set -u

orthofinder -f "$OF3_IN" -o "$OF3_OUT" \
            -S diamond -M msa -A famsa -T fasttree -X \
            -t "$OF3_THREADS" -a "$OF3_ANALYSIS"

set +u
conda deactivate
set -u
echo "   finished $(date -u +%Y-%m-%dT%H:%M:%SZ)"

# --- the two-species run must be untouched ---------------------------------
md5sum "$ORTHOGROUPS" > "${GUARD}.after"
if ! diff -q "${GUARD}.before" "${GUARD}.after" >/dev/null; then
  echo "FATAL: the two-species Orthogroups.tsv CHANGED. Every conservation"
  echo "       result in the repo depends on it. Restore it before going on."
  exit 1
fi
echo "== two-species orthology verified byte-identical"

RES="$(ls -d "${OF3_OUT}"/Results_* 2>/dev/null | head -1)"
echo
echo "02 done. Results: ${RES}"
echo "  species tree : ${RES}/Species_Tree/SpeciesTree_rooted.txt"
echo "  1:1:1 set    : ${RES}/Orthogroups/Orthogroups_SingleCopyOrthologues.txt"
echo "Next: ./02b_build_triplets.sh"
