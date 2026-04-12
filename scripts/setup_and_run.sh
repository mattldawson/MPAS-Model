#!/usr/bin/env bash
# =============================================================================
# CheMPAS-A: one-command setup from a fresh clone.
#
# Builds the dev container, fetches physics externals, compiles MPAS-A with
# MUSICA, downloads the 480-km mesh, and runs the JW baroclinic wave test.
#
# Usage:  bash scripts/setup_and_run.sh [NPROCS]
#   NPROCS: MPI ranks for the JW test (default: 1)
#
# Safe to re-run — every step is skipped if its output already exists.
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
skip()  { echo "  (already done — skipping)"; }

# ---------- Step 1: Build container image ------------------------------------
step "Step 1/5: Build chempas-dev container image"
if ${CONTAINER_RT} image exists "${IMAGE}" 2>/dev/null; then
    skip
else
    ${CONTAINER_RT} build -f docker/Containerfile --target dev -t chempas-dev .
fi

# ---------- Step 2: Fetch physics externals ----------------------------------
step "Step 2/5: Fetch physics externals"

PHYS="src/core_atmosphere/physics"

if [ -d "${PHYS}/physics_mmm/.git" ]; then
    echo "  physics_mmm already present"
else
    rm -rf "${PHYS}/physics_mmm"
    git clone --depth 1 --branch 20250616-MPASv8.3 \
        https://github.com/NCAR/MMM-physics.git "${PHYS}/physics_mmm"
fi

if [ -d "${PHYS}/physics_noaa/UGWP/.git" ]; then
    echo "  physics_noaa/UGWP already present"
else
    rm -rf "${PHYS}/physics_noaa/UGWP"
    mkdir -p "${PHYS}/physics_noaa"
    git clone --depth 1 --branch MPAS_20241223 \
        https://github.com/NOAA-GSL/UGWP.git "${PHYS}/physics_noaa/UGWP"
fi

if [ -f "${PHYS}/physics_wrf/files/COMPATIBILITY" ]; then
    echo "  physics_wrf lookup tables already present"
else
    mkdir -p "${PHYS}/physics_wrf/files"
    echo "  Downloading WRF physics lookup tables..."
    wget -q --show-progress -O /tmp/mpas-data-v8.2.tar.gz \
        https://github.com/MPAS-Dev/MPAS-Data/archive/refs/tags/v8.2.tar.gz
    tar xzf /tmp/mpas-data-v8.2.tar.gz -C "${PHYS}/physics_wrf/files" \
        --strip-components=4 "MPAS-Data-8.2/atmosphere/physics_wrf/files"
    rm -f /tmp/mpas-data-v8.2.tar.gz
fi

# ---------- Step 3: Build MPAS -----------------------------------------------
step "Step 3/5: Build MPAS-A (atmosphere + init_atmosphere)"
if [ -x atmosphere_model ] && [ -x init_atmosphere_model ]; then
    skip
else
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
        # Build init_atmosphere core (AUTOCLEAN lets it re-compile the
        # shared framework that was built with different options above)
        make -j$(nproc) gnu CORE=init_atmosphere USE_PIO2=false AUTOCLEAN=true
    '
fi

# ---------- Step 4: Download mesh data ---------------------------------------
step "Step 4/5: Download 480-km mesh"
if [ -f data/x1.2562.grid.nc ]; then
    skip
else
    bash scripts/download_data.sh data
fi

# ---------- Step 5: Run JW test ----------------------------------------------
step "Step 5/5: Run JW baroclinic wave test (${NPROCS} MPI ranks)"
if [ -f data/jw_480km/output.nc ]; then
    echo "  output.nc already exists — delete data/jw_480km/ to re-run"
    skip
else
    ${CONTAINER_RT} run --rm -v "${MPAS_DIR}:/mpas:Z" -w /mpas "${IMAGE}" \
        bash scripts/run_jw_test.sh "${NPROCS}"
fi

# ---------- Done -------------------------------------------------------------
step "Done"
echo "  Output: data/jw_480km/output.nc"
ls -lh data/jw_480km/output.nc
echo ""
echo "  To view results:  jupyter notebook verification/"
