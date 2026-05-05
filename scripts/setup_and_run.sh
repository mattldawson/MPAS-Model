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
# Usage:  bash scripts/setup_and_run.sh [--force-rebuild] [--quick] [NPROCS]
#   --force-rebuild: delete and rebuild the container image
#   --quick: run short simulations for fast performance/debug iterations
#   NPROCS: MPI ranks for the JW test (default: 1)
# =============================================================================
set -euo pipefail

FORCE_REBUILD=false
QUICK_MODE=false

while [ "$#" -gt 0 ]; do
    case "$1" in
        --force-rebuild)
            FORCE_REBUILD=true
            shift
            ;;
        --quick)
            QUICK_MODE=true
            shift
            ;;
        *)
            break
            ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
MPAS_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
NPROCS="${1:-$(nproc)}"
CONTAINER_RT="${CONTAINER_RT:-podman}"   # or "docker"
IMAGE="localhost/chempas-dev"
MUSICA_GIT_REPOSITORY="${MUSICA_GIT_REPOSITORY:-https://github.com/mattldawson/musica.git}"
MUSICA_GIT_TAG="${MUSICA_GIT_TAG:-develop-something-ambitious}"
JW_RUN_DURATION="${JW_RUN_DURATION:-1_00:00:00}"
JW_OUTPUT_INTERVAL="${JW_OUTPUT_INTERVAL:-01:00:00}"

if [ "${QUICK_MODE}" = true ]; then
    JW_RUN_DURATION="00:03:00"
    JW_OUTPUT_INTERVAL="00:03:00"
fi

cd "${MPAS_DIR}"

# ---------- helpers ----------------------------------------------------------
step()  { echo ""; echo "====== $* ======"; }

build_image() {
    ${CONTAINER_RT} build -f docker/Containerfile --target dev -t "${IMAGE}" \
        --build-arg MUSICA_GIT_REPOSITORY="${MUSICA_GIT_REPOSITORY}" \
        --build-arg MUSICA_GIT_TAG="${MUSICA_GIT_TAG}" .
}

has_required_musica_api() {
        ${CONTAINER_RT} run --rm "${IMAGE}" bash -lc '
cat > /tmp/musica_api_probe.f90 <<"EOF"
program musica_api_probe
    use musica_micm, only : RosenbrockDAE4StandardOrder, rosenbrock_solver_parameters_t
    implicit none
    integer :: solver_type
    type(rosenbrock_solver_parameters_t) :: params
    solver_type = RosenbrockDAE4StandardOrder
    params%constraint_init_max_iterations = 100
    params%constraint_init_tolerance = 1.0d-8
    print *, solver_type
end program musica_api_probe
EOF
mpif90 -c /tmp/musica_api_probe.f90 -I/usr/local/include/musica/fortran -J/tmp
'
}

# ---------- Step 1: Build container image ------------------------------------
step "Step 1/10: Build chempas-dev container image"
if [ "${FORCE_REBUILD}" = true ]; then
    echo "  (--force-rebuild: removing existing image)"
    ${CONTAINER_RT} rmi -f "${IMAGE}" 2>/dev/null || true
fi
if ${CONTAINER_RT} image exists "${IMAGE}" 2>/dev/null; then
    echo "  (image exists — skipping build)"
else
    build_image
fi

# Validate the container has required TUV-x data
if ! ${CONTAINER_RT} run --rm "${IMAGE}" test -d /usr/local/share/musica/tuvx_data/quantum_yields; then
    echo "  Container image is stale (missing TUV-x data). Rebuilding..."
    ${CONTAINER_RT} rmi -f "${IMAGE}" 2>/dev/null || true
    build_image
fi

# Validate the container has the MUSICA API expected by MPAS chemistry.
if ! has_required_musica_api; then
    echo "  Container image is stale (MUSICA API mismatch). Rebuilding..."
    ${CONTAINER_RT} rmi -f "${IMAGE}" 2>/dev/null || true
    build_image
    if ! has_required_musica_api; then
        echo "ERROR: Rebuilt container still has an incompatible MUSICA API." >&2
        exit 1
    fi
fi

# ---------- Step 2: Fetch physics externals ----------------------------------
step "Step 2/10: Fetch physics externals"

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
step "Step 3/10: Build MPAS-A (atmosphere + init_atmosphere)"
rm -f atmosphere_model init_atmosphere_model
${CONTAINER_RT} run --rm -v "${MPAS_DIR}:/mpas:Z" -w /mpas "${IMAGE}" bash -c '
    set -euo pipefail

    # No-op the external-checkout scripts (we already fetched everything)
    printf "#!/bin/sh\nexec true\n" > src/core_atmosphere/tools/manage_externals/checkout_externals
    chmod +x src/core_atmosphere/tools/manage_externals/checkout_externals
    printf "#!/bin/sh\nexec true\n" > src/core_atmosphere/physics/checkout_data_files.sh
    chmod +x src/core_atmosphere/physics/checkout_data_files.sh

    # Clean everything first to avoid stale object files
    make clean CORE=atmosphere 2>/dev/null || true

    # MPAS Makefile has incomplete dependency tracking; high parallelism
    # can cause race conditions.  Cap at 8 and retry once on failure.
    JLEVEL=$(( $(nproc) > 8 ? 8 : $(nproc) ))

    # Build atmosphere core (needs MUSICA flags)
    make -j${JLEVEL} gnu CORE=atmosphere USE_PIO2=false \
        MPAS_EXTERNAL_LIBS="$(pkg-config --libs musica-fortran) -lstdc++" \
        MPAS_EXTERNAL_INCLUDES="$(pkg-config --cflags musica-fortran)" \
    || make -j${JLEVEL} gnu CORE=atmosphere USE_PIO2=false \
        MPAS_EXTERNAL_LIBS="$(pkg-config --libs musica-fortran) -lstdc++" \
        MPAS_EXTERNAL_INCLUDES="$(pkg-config --cflags musica-fortran)"
    # Save atmosphere_model — AUTOCLEAN below will remove it
    cp atmosphere_model /tmp/atmosphere_model
    # Build init_atmosphere core (AUTOCLEAN re-compiles the shared
    # framework that was built with different options above)
    make -j${JLEVEL} gnu CORE=init_atmosphere USE_PIO2=false AUTOCLEAN=true \
    || make -j${JLEVEL} gnu CORE=init_atmosphere USE_PIO2=false AUTOCLEAN=true
    # Restore atmosphere_model
    cp /tmp/atmosphere_model atmosphere_model
'

# ---------- Step 4: Download mesh data ---------------------------------------
step "Step 4/10: Download 480-km mesh"
if [ -f data/x1.2562.grid.nc ]; then
    echo "  (mesh exists — skipping download)"
else
    bash scripts/download_data.sh data
fi

# ---------- Step 5: Run JW tests ---------------------------------------------
echo "  JW runtime controls: run_duration=${JW_RUN_DURATION}, output_interval=${JW_OUTPUT_INTERVAL}"

step "Step 5/10: Run JW baroclinic wave test — chapman (${NPROCS} MPI ranks)"
rm -rf data/jw_480km_chapman
${CONTAINER_RT} run --rm -e JW_RUN_DURATION -e JW_OUTPUT_INTERVAL -v "${MPAS_DIR}:/mpas:Z" -w /mpas "${IMAGE}" \
    bash scripts/run_jw_test.sh "${NPROCS}" chapman

step "Step 6/10: Run JW baroclinic wave test — analytical (${NPROCS} MPI ranks)"
rm -rf data/jw_480km_analytical
${CONTAINER_RT} run --rm -e JW_RUN_DURATION -e JW_OUTPUT_INTERVAL -v "${MPAS_DIR}:/mpas:Z" -w /mpas "${IMAGE}" \
    bash scripts/run_jw_test.sh "${NPROCS}" analytical

step "Step 7/10: Run JW baroclinic wave test — chapman_emis_dep (${NPROCS} MPI ranks)"
rm -rf data/jw_480km_chapman_emis_dep
${CONTAINER_RT} run --rm -e JW_RUN_DURATION -e JW_OUTPUT_INTERVAL -v "${MPAS_DIR}:/mpas:Z" -w /mpas "${IMAGE}" \
    bash scripts/run_jw_test.sh "${NPROCS}" chapman_emis_dep

step "Step 8/10: Run JW baroclinic wave test — ts1 (${NPROCS} MPI ranks)"
rm -rf data/jw_480km_ts1
${CONTAINER_RT} run --rm -e JW_RUN_DURATION -e JW_OUTPUT_INTERVAL -v "${MPAS_DIR}:/mpas:Z" -w /mpas "${IMAGE}" \
    bash scripts/run_jw_test.sh "${NPROCS}" ts1

step "Step 9/10: Run JW baroclinic wave test — ts1_cloud (${NPROCS} MPI ranks)"
rm -rf data/jw_480km_ts1_cloud
${CONTAINER_RT} run --rm -e JW_RUN_DURATION -e JW_OUTPUT_INTERVAL -v "${MPAS_DIR}:/mpas:Z" -w /mpas "${IMAGE}" \
    bash scripts/run_jw_test.sh "${NPROCS}" ts1_cloud

step "Step 10/10: Run JW baroclinic wave test — ts1_cloud_dry (${NPROCS} MPI ranks)"
rm -rf data/jw_480km_ts1_cloud_dry
${CONTAINER_RT} run --rm -e JW_RUN_DURATION -e JW_OUTPUT_INTERVAL -v "${MPAS_DIR}:/mpas:Z" -w /mpas "${IMAGE}" \
    bash scripts/run_jw_test.sh "${NPROCS}" ts1_cloud_dry

# Backward-compatible symlink for phases 0-2 (which expect data/jw_480km/)
rm -rf data/jw_480km
ln -sfn jw_480km_chapman data/jw_480km

# ---------- Done -------------------------------------------------------------
step "Done"
echo "  Output:"
ls -lh data/jw_480km_chapman/output.nc data/jw_480km_analytical/output.nc \
       data/jw_480km_chapman_emis_dep/output.nc data/jw_480km_ts1/output.nc \
       data/jw_480km_ts1_cloud/output.nc data/jw_480km_ts1_cloud_dry/output.nc
echo ""
echo "  To view results:  jupyter notebook verification/"
