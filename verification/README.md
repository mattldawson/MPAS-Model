# CheMPAS-A Verification Notebooks

Visual verification notebooks for each development phase. Each notebook loads
MPAS output from `data/` (produced by a container run) and plots key fields so
you can see at a glance that the phase is working.

## Quick Start

```bash
# 1. Build the container (from the MPAS-Model directory)
cd MPAS-Model
podman build -f docker/Containerfile --target deps -t chempas-deps .
podman build -f docker/Containerfile --target dev  -t chempas-dev  .

# 2. Fetch physics externals (only needed once after a fresh clone)
#    This replaces the fragile manage_externals tool with direct git clones.
PHYS=src/core_atmosphere/physics
git clone --depth 1 --branch 20250616-MPASv8.3 \
    https://github.com/NCAR/MMM-physics.git $PHYS/physics_mmm
mkdir -p $PHYS/physics_noaa
git clone --depth 1 --branch MPAS_20241223 \
    https://github.com/NOAA-GSL/UGWP.git $PHYS/physics_noaa/UGWP
mkdir -p $PHYS/physics_wrf/files && cd $PHYS/physics_wrf/files
wget -q https://github.com/MPAS-Dev/MPAS-Data/archive/refs/tags/v8.2.tar.gz
tar xzf v8.2.tar.gz --strip-components=4 "MPAS-Data-8.2/atmosphere/physics_wrf/files"
rm v8.2.tar.gz
cd -

# 3. Build MPAS inside the container
podman run --rm -v "$(pwd):/mpas:Z" -w /mpas localhost/chempas-dev bash -c '
  printf "#!/bin/sh\nexec true\n" > src/core_atmosphere/tools/manage_externals/checkout_externals
  chmod +x src/core_atmosphere/tools/manage_externals/checkout_externals
  printf "#!/bin/sh\nexec true\n" > src/core_atmosphere/physics/checkout_data_files.sh
  chmod +x src/core_atmosphere/physics/checkout_data_files.sh
  make -j$(nproc) gnu CORE=atmosphere USE_PIO2=false \
    MPAS_EXTERNAL_LIBS="$(pkg-config --libs musica-fortran) -lstdc++" \
    MPAS_EXTERNAL_INCLUDES="$(pkg-config --cflags musica-fortran)"
'

# 4. Run the JW test (produces data/jw_480km/output.nc)
podman run --rm \
    -v "$(pwd):/mpas:Z" \
    -w /mpas localhost/chempas-dev \
    bash scripts/run_jw_test.sh 1

# 5. Set up a Python virtual environment
python3 -m venv .venv
source .venv/bin/activate
pip install -r verification/requirements.txt

# 6. Open the notebooks
jupyter notebook verification/
# — or open in VS Code —
```

## Notebooks

| Notebook | Phase | What it shows |
|----------|-------|---------------|
| `phase00_jw_baseline.ipynb` | 0 | Surface pressure, potential temperature, zonal wind from the JW baroclinic wave |
| `phase01_tracers.ipynb` | 1 | Passive tracer mass conservation and spatial distribution |
| `phase02_chapman.ipynb` | 2 | O3 diurnal cycle, comparison with tutorial 10 |
| `phase03_runtime_species.ipynb` | 3 | Mechanism-switching: Chapman vs analytical species in output |
| `phase04_emissions_deposition.ipynb` | 4 | Mass budget closure: emissions − deposition ± chemistry |
| `phase05_tuvx_photolysis.ipynb` | 5 | Photolysis rate profiles, regression against Phase 2 |
| `phase06_ts1.ipynb` | 6 | TS1 O3/NO/NO2/CO diurnal patterns, negative species check |
| `phase07_miam_bindings.ipynb` | 7 | MIAM Fortran unit test output vs tutorial 14 |
| `phase08_aerosol_config.ipynb` | 8 | Aerosol fields from MIAM representation config |
| `phase09_cloud_chemistry.ipynb` | 9 | Sulfate production in cloudy cells, DAE convergence |
| `phase10_regression.ipynb` | 10 | Full-stack regression, multi-resolution comparison |

## Data Layout

The notebooks expect output files under `data/` (gitignored). The JW test
script writes output to a single directory used by all current phases:

```
data/
  jw_480km/output.nc   ← Phases 0, 1, and 2 (includes tracers + chemistry species)
```

The output contains all registered scalars (`tracer_1`, `tracer_2`, `tracer_3`,
`o3`, `o`, `o1d`) plus standard dynamical fields (`theta`, `uReconstructZonal`,
`uReconstructMeridional`, `surface_pressure`, etc.).

## Dependencies

All dependencies are pure Python wheels — no system libraries required:

- `netCDF4` — reads MPAS NetCDF output
- `matplotlib` — plotting
- `numpy` — array operations
- `jupyter` — notebook runtime
