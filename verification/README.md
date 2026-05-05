# CheMPAS-A Verification Notebooks

Visual verification notebooks for each development phase. Each notebook loads
MPAS output from `data/` (produced by a container run) and plots key fields so
you can see at a glance that the phase is working.

## Quick Start

```bash
cd MPAS-Model
bash scripts/setup_and_run.sh      # builds container, compiles MPAS, runs JW test

# Fast debug loop (short MPAS runtime, much quicker turnaround)
# bash scripts/setup_and_run.sh --quick 1

# Set up Python for the notebooks
python3 -m venv .venv
source .venv/bin/activate
pip install -r verification/requirements.txt

# Open the notebooks
jupyter notebook verification/
# — or open in VS Code —
```

The setup script validates the cached container image before building/running.
If the image is stale (for example missing TUV-x data or using an incompatible
MUSICA API), it is rebuilt automatically. The container builds MUSICA from the
latest commit on the `develop-something-ambitious` branch from
`https://github.com/mattldawson/musica.git`. Set `CONTAINER_RT=docker` if you
use Docker instead of Podman. You can override the source explicitly via
`MUSICA_GIT_REPOSITORY` and `MUSICA_GIT_TAG`.

For performance debugging, set a shorter simulation horizon with
`--quick` (3-minute run) or override explicitly using
`JW_RUN_DURATION` and `JW_OUTPUT_INTERVAL`.

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
| `phase08_cloud_chemistry.ipynb` | 8 | Cloud chemistry (`ts1_cloud` vs `ts1_cloud_dry` control): in-cloud signal localized to the prescribed 700–850 hPa layer |
| `phase09_regression.ipynb` | 9 | Full-stack regression, multi-resolution comparison |

## Data Layout

The notebooks expect output files under `data/` (gitignored). The JW test
script writes output to a single directory used by all current phases:

```
data/
  jw_480km_chapman/output.nc          ← Phases 0, 1, 2, and 3 (Chapman chemistry)
  jw_480km_analytical/output.nc       ← Phase 3 (analytical A→B→C mechanism)
  jw_480km_chapman_emis_dep/output.nc ← Phase 4 (emissions/deposition stubs)
  jw_480km_ts1/output.nc              ← Phases 5–6 (TS1 gas-phase mechanism)
  jw_480km_ts1_cloud/output.nc        ← Phase 8 (TS1 + CAM Cloud Chemistry, LWC=3e-4)
  jw_480km_ts1_cloud_dry/output.nc    ← Phase 8 dry control (same mechanism, LWC=0)
  jw_480km -> jw_480km_chapman        ← symlink for backward compatibility
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

## Run all notebooks headlessly

After `setup_and_run.sh` produces all six `data/jw_480km_*/output.nc` files,
you can execute every notebook in place (figures get embedded as cell outputs)
with:

```bash
bash scripts/run_all_notebooks.sh
```

This loops over `verification/phase*.ipynb` with `jupyter nbconvert --execute
--inplace` and exits non-zero if any notebook raises. Notebook outputs are
**not** committed to the repository — running locally is the way to see the
figures.
