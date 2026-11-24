# QUACS Development Instructions

This document describes how to get up-and-running as a QUACS developer. It covers:
* Setting up your VSCode to develop in the CheMPAS-A container environment
* Building and running CheMPAS-A in your development container
* Alternative command-line option for building and running CHeMPAS-A in a container

---

## Prerequisites

* [Docker](https://docs.docker.com/get-docker/) or [Podman](https://podman.io/getting-started/installation) (with the `docker` CLI alias)
* [VS Code](https://code.visualstudio.com/) with the [Dev Containers](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers) extension
* Git

---

## Setup
### Option 1: Developing inside the VS Code Dev Container (recommended)

The `.devcontainer/` configuration mounts the repository source at `/mpas` inside a pre-built container that has all build dependencies (MPI, NetCDF, PnetCDF, MUSICA, QUACS, etc.) already installed.

1. **Open the repo in VS Code** and, when prompted, click **Reopen in Container**.  
   Alternatively, open the Command Palette (`Ctrl+Shift+P`) and run **Dev Containers: Reopen in Container**.

2. VS Code will pull the `dev` image, start the container, and open a terminal at `/mpas`.

3. **Build the init_atmosphere core**:
   ```bash
   make -j$(nproc) gnu CORE=init_atmosphere AUTOCLEAN=true MUSICA=true
   ```

4. **Build the atmosphere core**:
   ```bash
   make -j$(nproc) gnu CORE=atmosphere AUTOCLEAN=true MUSICA=true QUACS=true
   ```
   Both executables (`init_atmosphere_model`, `atmosphere_model`) are placed in `/mpas`.

---

### Option 2: Building the container image locally (alternative approach)

0. This is only needed if you want to build the `dev` Docker image from scratch (e.g., to modify the `Containerfile`).

   ```bash
   docker build --target dev -t chempas-deps -f docker/Containerfile .
   ```

1. To build the full runtime image with MPAS executables baked in (the `build` stage clones physics externals and downloads WRF lookup tables automatically):
   ```bash
   docker build -t chempas --target build -f docker/Containerfile .
   ```

2. You can then run the container and proceed to the next section:
   ```bash
   docker run -it chempas bash
   ```
---

## Downloading test data and running the JW baroclinic wave test

These steps assume the executables have already been built (see above).

1. **Download the 480-km mesh and test case** (into the `data/` directory):
   ```bash
   quacs/download_data.sh data
   ```
   Pass `--all` as a second argument to also fetch the 240-km mesh:
   ```bash
   quacs/download_data.sh data --all
   ```

2. **Run the Jablonowski-Williamson baroclinic wave test**:
   ```bash
   quacs/run_jw_test.sh [NPROCS]
   ```
   `NPROCS` is the number of MPI ranks (default: 1). For example, with 2 ranks:
   ```bash
   quacs/run_jw_test.sh 2
   ```
   Output is written to `data/jw_480km/output.nc`.

   The test runs the QUACS conformance checks in `quacs/conformance_config.py` at the
   start of the run and in every time step, and fails unless every check passes on
   every MPI rank. The results are at the end of `data/jw_480km/log.atmosphere.0000.out`.

---

## Running QUACS parameterizations

QUACS runs when MPAS is built with `QUACS=true` and the `&quacs` namelist record names a
QUACS configuration file:

```
&quacs
    config_quacs_config_file = 'my_config.py'
/
```

The configuration file is a Python file that sets `parameterizations` and `hooks`
(see `embed/README.md` in the QUACS repository). It must list both of the hooks that
MPAS runs, even if no parameterizations are run in one of them:

| Hook | When MPAS runs it |
|------|-------------------|
| `initial_conditions` | Once, at the start of a non-restart run, before the first time step. |
| `before_dynamics` | Once per time step, before the dynamics |

### How MPAS provides variables to QUACS

- A variable that a QUACS parameterization owns is created by MPAS. If it is transported, it
  becomes a constituent of `scalars`. Otherwise, it becomes a field that can be added to an
  output stream by name.
- Any other variable must be an MPAS field of the same name, or a constituent of `scalars`
  such as `qv`.
- Parameterizations work on every cell of a block, including halo cells.

The container image installs QUACS in the Python virtual environment `/usr/local/quacs-venv`.
Set the `QUACS_PYTHON` environment variable to the `python` of another virtual environment
where QUACS is installed to use that QUACS instead.

## CI/CD

The GitHub Actions workflow (`.github/workflows/chempas.yml`) runs on every push and pull request to the `proto-quacs` branch. It:

1. **Builds the `dev` dependency image** and pushes it to the GitHub Container Registry (`ghcr.io`), using a SHA-tagged image for reproducibility and a `latest` tag on the default branch.
2. **Builds both MPAS cores** (`init_atmosphere` and `atmosphere`) inside that image.
3. **Runs the JW baroclinic wave test** on the 480-km mesh with 2 MPI ranks and verifies that `output.nc` is produced.

