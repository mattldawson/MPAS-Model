#!/usr/bin/env bash
# =============================================================================
# CheMPAS-A: one-command setup.
#
# Builds the dev container, fetches physics externals, compiles MPAS-A with
# MUSICA, downloads the 480-km mesh, and runs the JW baroclinic wave test.
#
# Idempotent — always rebuilds from clean to avoid stale-artifact issues.
# Only the container image and mesh download are cached (slow to fetch).
#
# Usage:  bash scripts/setup_and_run.sh [NPROCS]
#   NPROCS: MPI ranks for the JW test (default: 1)
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
MPAS_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
NPROCS="${1:-1}"
CONTAINER_RT="${CONTAINER_RT:-podman}"   # or "docker"
IMAGE="localhost/chempas-dev"

cd "${MPAS_DIR}"

# ---------- helpers ----------------------------------------------------------
step()  { echo ""; echo "====== $* ======"; }

# ---------- Step 1: Build container image ------------------------------------
step "Step 1/6: Build chempas-dev container image"
if ${CONTAINER_RT} image exists "${IMAGE}" 2>/dev/null; then
    echo "  (image exists — skipping build)"
else
    ${CONTAINER_RT} build -f docker/Containerfile --target dev -t chempas-dev .
fi

# ---------- Step 2: Fetch physics externals ----------------------------------
step "Step 2/6: Fetch physics externals"

PHYS="src/core_atmosphere/physics"

rm -rf "${PHYS}/physics_mmm"
git clone --depth 1 --branch 20250616-MPASv8.3 \
    https://github.com/NCAR/MMM-physics.git "${PHYS}/physics_mmm"

rm -rf "${PHYS}/physics_noaa/UGWP"
mkdir -p "${PHYS}/physics_noaa"
git clone --depth 1 --branch MPAS_20241223 \
    https://github.com/NOAA-GSL/UGWP.git "${PHYS}/physics_noaa/UGWP"

rm -rf "${PHYS}/physics_wrf/files"
mkdir -p "${PHYS}/physics_wrf/files"
wget -q --show-progress -O /tmp/mpas-data-v8.2.tar.gz \
    https://github.com/MPAS-Dev/MPAS-Data/archive/refs/tags/v8.2.tar.gz
tar xzf /tmp/mpas-data-v8.2.tar.gz -C "${PHYS}/physics_wrf/files" \
    --strip-components=4 "MPAS-Data-8.2/atmosphere/physics_wrf/files"
rm -f /tmp/mpas-data-v8.2.tar.gz

# ---------- Step 3: Build MPAS -----------------------------------------------
step "Step 3/6: Build MPAS-A (atmosphere + init_atmosphere)"
rm -f atmosphere_model init_atmosphere_model
${CONTAINER_RT} run --rm -v "${MPAS_DIR}:/mpas:Z" -w /mpas "${IMAGE}" bash -c '
    # No-op the external-checkout scripts (we already fetched everything)
    printf "#!/bin/sh\nexec true\n" > src/core_atmosphere/tools/manage_externals/checkout_externals
    chmod +x src/core_atmosphere/tools/manage_externals/checkout_externals
    printf "#!/bin/sh\nexec true\n" > src/core_atmosphere/physics/checkout_data_files.sh
    chmod +x src/core_atmosphere/physics/checkout_data_files.sh

    # Build atmosphere core first (needs MUSICA flags)
    make -j$(nproc) gnu CORE=atmosphere USE_PIO2=false \
        MPAS_EXTERNAL_LIBS="$(pkg-config --libs musica-fortran) -lstdc++" \
        MPAS_EXTERNAL_INCLUDES="$(pkg-config --cflags musica-fortran)"
    # Save atmosphere_model — AUTOCLEAN below will remove it
    cp atmosphere_model /tmp/atmosphere_model
    # Build init_atmosphere core (AUTOCLEAN re-compiles the shared
    # framework that was built with different options above)
    make -j$(nproc) gnu CORE=init_atmosphere USE_PIO2=false AUTOCLEAN=true
    # Restore atmosphere_model
    cp /tmp/atmosphere_model atmosphere_model
'

# ---------- Step 4: Download mesh data ---------------------------------------
step "Step 4/6: Download 480-km mesh"
if [ -f data/x1.2562.grid.nc ]; then
    echo "  (mesh exists — skipping download)"
else
    bash scripts/download_data.sh data
fi

# ---------- Step 5: Run JW tests ---------------------------------------------
step "Step 5/6: Run JW baroclinic wave test — chapman (${NPROCS} MPI ranks)"
rm -rf data/jw_480km_chapman
${CONTAINER_RT} run --rm -v "${MPAS_DIR}:/mpas:Z" -w /mpas "${IMAGE}" \
    bash scripts/run_jw_test.sh "${NPROCS}" chapman

step "Step 6/6: Run JW baroclinic wave test — analytical (${NPROCS} MPI ranks)"
rm -rf data/jw_480km_analytical
${CONTAINER_RT} run --rm -v "${MPAS_DIR}:/mpas:Z" -w /mpas "${IMAGE}" \
    bash scripts/run_jw_test.sh "${NPROCS}" analytical

# Backward-compatible symlink for phases 0-2 (which expect data/jw_480km/)
ln -sfn jw_480km_chapman data/jw_480km

# ---------- Done -------------------------------------------------------------
step "Done"
echo "  Output:"
ls -lh data/jw_480km_chapman/output.nc data/jw_480km_analytical/output.nc
echo ""
echo "  To view results:  jupyter notebook verification/"
