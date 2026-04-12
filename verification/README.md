# CheMPAS-A Verification Notebooks

Visual verification notebooks for each development phase. Each notebook loads
MPAS output from `data/` (produced by a container run) and plots key fields so
you can see at a glance that the phase is working.

## Quick Start

```bash
cd MPAS-Model
bash scripts/setup_and_run.sh      # builds container, compiles MPAS, runs JW test

# Set up Python for the notebooks
python3 -m venv .venv
source .venv/bin/activate
pip install -r verification/requirements.txt

# Open the notebooks
jupyter notebook verification/
# — or open in VS Code —
```

The setup script is safe to re-run — every step is skipped if its output
already exists.  Delete `data/jw_480km/` to force a fresh run, or
`atmosphere_model` to force a rebuild.  Set `CONTAINER_RT=docker` if you
use Docker instead of Podman.

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
  jw_480km_chapman/output.nc    ← Phases 0, 1, 2, and 3 (Chapman chemistry)
  jw_480km_analytical/output.nc ← Phase 3 (analytical A→B→C mechanism)
  jw_480km -> jw_480km_chapman  ← symlink for backward compatibility
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
