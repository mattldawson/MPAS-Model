# CheMPAS-A Verification Notebooks

Visual verification notebooks for each development phase. Each notebook loads
MPAS output from `data/` (produced by a container run) and plots key fields so
you can see at a glance that the phase is working.

## Quick Start

```bash
# 1. Build the container and run the JW test (produces data/jw_480km/output.nc)
#    See docker/Containerfile and scripts/ for details.
bash scripts/download_data.sh data
podman build -f docker/Containerfile --target build -t chempas-build .
podman run --rm \
    -v "$(pwd)/data:/mpas/data:Z" \
    -v "$(pwd)/scripts:/mpas/scripts:Z" \
    -w /mpas chempas-build \
    bash scripts/run_jw_test.sh 2

# 2. Set up a Python virtual environment
python3 -m venv .venv
source .venv/bin/activate
pip install -r verification/requirements.txt

# 3. Open the notebooks
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

The notebooks expect output files under `data/` (gitignored). Each phase's
container run script places output in a phase-specific subdirectory:

```
data/
  jw_480km/output.nc          ← Phase 0 (and base for later phases)
  jw_480km_tracers/output.nc  ← Phase 1
  jw_480km_chapman/output.nc  ← Phase 2
  ...
```

## Dependencies

All dependencies are pure Python wheels — no system libraries required:

- `netCDF4` — reads MPAS NetCDF output
- `matplotlib` — plotting
- `numpy` — array operations
- `jupyter` — notebook runtime
