#!/usr/bin/env bash
# =============================================================================
# Execute every verification/phase*.ipynb notebook in place using nbconvert.
# Exits non-zero on the first failure. Run this AFTER `setup_and_run.sh` has
# produced all the data/jw_480km_*/output.nc files.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
MPAS_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
NB_DIR="${MPAS_DIR}/verification"

if ! command -v jupyter >/dev/null 2>&1; then
    echo "ERROR: jupyter not found on PATH." >&2
    echo "  Activate the venv:  source .venv/bin/activate" >&2
    echo "  Or install:         pip install -r verification/requirements.txt" >&2
    exit 1
fi

# Generous per-cell timeout: the Phase 8 wet-vs-dry comparison loads two
# 24-h netCDF files and computes whole-domain area-weighted means.
TIMEOUT="${NBCONVERT_TIMEOUT:-600}"

shopt -s nullglob
notebooks=("${NB_DIR}"/phase*.ipynb)
if [ "${#notebooks[@]}" -eq 0 ]; then
    echo "ERROR: No phase*.ipynb files under ${NB_DIR}" >&2
    exit 1
fi

failed=()
for nb in "${notebooks[@]}"; do
    echo ""
    echo "====== Executing $(basename "${nb}") ======"
    if jupyter nbconvert --to notebook --execute --inplace \
        --ExecutePreprocessor.timeout="${TIMEOUT}" "${nb}"; then
        echo "  PASS  $(basename "${nb}")"
    else
        echo "  FAIL  $(basename "${nb}")"
        failed+=("$(basename "${nb}")")
    fi
done

echo ""
echo "====== Notebook execution summary ======"
echo "  Total:  ${#notebooks[@]}"
echo "  Passed: $(( ${#notebooks[@]} - ${#failed[@]} ))"
echo "  Failed: ${#failed[@]}"
if [ "${#failed[@]}" -gt 0 ]; then
    printf '    - %s\n' "${failed[@]}"
    exit 1
fi
