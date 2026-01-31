# Installation Guide

Complete installation instructions for the QM NMR Calculator, from system dependencies to environment validation.

**Target audience:** Academic researchers and developers comfortable with command line basics.

## Prerequisites Overview

Before installing the QM NMR Calculator, you'll need:

- **Linux operating system** (Ubuntu 20.04+, Fedora 35+, or similar)
- **Python 3.11 or later**
- **NWChem 7.0.2+** quantum chemistry package
- **OpenMPI** for parallel DFT calculations
- **uv package manager** for Python dependency management
- **Git** for cloning the repository

Optional for high-accuracy conformer ensembles:
- **CREST** conformer search tool
- **xTB** semi-empirical calculator

## System Dependencies

### NWChem Quantum Chemistry Package

NWChem is the computational engine that performs DFT geometry optimizations and NMR shielding tensor calculations.

**Ubuntu/Debian installation:**

```bash
sudo apt update
sudo apt install nwchem
```

**Fedora/RHEL installation:**

```bash
sudo dnf install nwchem
```

**Building from source:**

If NWChem is not available in your distribution's repositories, or you need a specific version:

1. Visit [NWChem Documentation](https://nwchemgit.github.io/)
2. Follow the build instructions for your platform
3. Ensure the `nwchem` binary is added to your `PATH`

**Verify installation:**

```bash
nwchem --version
# Or if your installation doesn't support --version flag:
which nwchem
# Should return: /usr/bin/nwchem (or similar path)
```

**Required version:** NWChem 7.0.2 or later

> **Note:** NWChem is essential for this application. The calculator cannot function without a working NWChem installation.

### OpenMPI for Parallel Execution

OpenMPI enables multi-core parallelization for DFT calculations, significantly reducing computation time.

**Ubuntu/Debian installation:**

```bash
sudo apt install openmpi-bin libopenmpi-dev
```

**Fedora/RHEL installation:**

```bash
sudo dnf install openmpi openmpi-devel
```

**Why OpenMPI matters:**

- Draft preset calculations: 5 minutes (single core) → 2 minutes (4 cores)
- Production preset calculations: 30 minutes (single core) → 10 minutes (4 cores)

> **Note:** The application works without MPI but runs slower. For production use, OpenMPI is strongly recommended.

### Python 3.11 or Later

The application requires Python 3.11+ for modern type hints and performance optimizations.

**Check current version:**

```bash
python3 --version
# Should show: Python 3.11.x or later
```

**Ubuntu/Debian installation (if needed):**

```bash
sudo apt update
sudo apt install python3.11 python3.11-venv python3.11-dev
```

**Fedora/RHEL installation:**

```bash
sudo dnf install python3.11 python3.11-devel
```

> **Note:** Python 3.10 and earlier are not supported due to dependency requirements.

## Project Setup with uv

### What is uv?

[uv](https://github.com/astral-sh/uv) is a fast Python package manager from Astral (creators of Ruff). It's significantly faster than pip and handles virtual environment creation automatically.

Key features:
- 10-100x faster than pip for dependency resolution
- Automatic virtual environment management
- Lock file support for reproducible installs
- Drop-in pip replacement

### Install uv

**Recommended method (shell script):**

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

**Alternative (via pip):**

```bash
pip install uv
```

**Verify installation:**

```bash
uv --version
# Should show: uv 0.x.x or later
```

### Clone and Install the Project

**Clone the repository:**

```bash
git clone https://github.com/steinbeck/qm-nmr-calc.git
cd qm-nmr-calc
```

**Install all dependencies:**

```bash
uv sync
```

**What `uv sync` does:**

1. Creates a virtual environment at `.venv/` (if it doesn't exist)
2. Reads `pyproject.toml` for dependency requirements
3. Installs all runtime and development dependencies
4. Creates a lock file for reproducible builds

The virtual environment includes:
- FastAPI for the web server
- Huey for async task queue
- RDKit for molecular structure handling
- NumPy, SciPy, Pandas for numerical processing
- Matplotlib for spectrum visualization
- All other dependencies from `pyproject.toml`

**Activate the virtual environment (optional):**

```bash
source .venv/bin/activate
```

> **Note:** You don't need to activate the environment manually. You can use `uv run` to execute commands in the virtual environment directly (recommended).

**Test the installation:**

```bash
# Using uv run (no activation needed)
uv run python -c "import rdkit; print('RDKit version:', rdkit.__version__)"

# Or with activated environment
source .venv/bin/activate
python -c "import rdkit; print('RDKit version:', rdkit.__version__)"
```

If this prints the RDKit version without errors, your basic installation is complete.

