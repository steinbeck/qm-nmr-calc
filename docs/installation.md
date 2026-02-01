# Installation Guide

Complete installation instructions for the QM NMR Calculator, from system dependencies to environment validation.

**Target audience:** Academic researchers and developers comfortable with command line basics.

## Quick Start

**Prerequisites:** Linux, Python 3.11+, NWChem, uv

```bash
# Clone and install
git clone https://github.com/steinbeck/qm-nmr-calc.git
cd qm-nmr-calc
uv sync

# Start services (two terminals)
uv run python scripts/run_consumer.py  # Terminal 1: Job processor
uv run python scripts/run_api.py       # Terminal 2: API server

# Open http://localhost:8000
```

For detailed setup including optional CREST/xTB for ensemble calculations, continue reading below.

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

## Optional: CREST and xTB Setup

### Purpose

CREST (Conformer-Rotamer Ensemble Sampling Tool) and xTB (extended Tight-Binding) enable high-accuracy conformer generation for ensemble calculations.

**What they do:**

- **xTB**: Semi-empirical quantum chemistry method (GFN2-xTB force field)
- **CREST**: Automated conformer search using metadynamics
- **Together**: Generate comprehensive conformer ensembles for flexible molecules

**When to use:**

- Flexible molecules with multiple rotatable bonds
- Molecules where conformer populations affect NMR shifts
- Publication-quality results requiring ensemble averaging

**What happens without them:**

The application works fully without CREST/xTB. It will use RDKit's Knowledge-based Distance Geometry (KDG) method instead, which is:
- Faster (seconds vs minutes)
- Good for rigid molecules
- Less exhaustive for highly flexible systems

> **Note:** Both binaries (`crest` and `xtb`) must be available on your PATH for CREST mode to be enabled.

### Installation via Conda (Recommended)

The easiest way to install both tools:

```bash
# Install Miniconda if you don't have conda
# See: https://docs.conda.io/en/latest/miniconda.html

# Create a new environment (optional but recommended)
conda create -n crest-env
conda activate crest-env

# Install both tools
conda install -c conda-forge crest xtb
```

**Verify installation:**

```bash
which crest && crest --version
which xtb && xtb --version
```

Both commands should return paths and version information.

### Manual Installation (Advanced)

If you prefer to install binaries manually:

**Download binaries:**

- **CREST**: https://github.com/crest-lab/crest/releases
- **xTB**: https://github.com/grimme-lab/xtb/releases

**Extract and add to PATH:**

```bash
# Example for CREST
tar -xf crest-3.0.tar.gz
sudo mv crest /usr/local/bin/

# Example for xTB
tar -xf xtb-6.7.1.tar.gz
cd xtb-6.7.1/bin
sudo mv xtb /usr/local/bin/
```

**Verify PATH configuration:**

```bash
which crest  # Should return: /usr/local/bin/crest
which xtb    # Should return: /usr/local/bin/xtb
```

### How the Application Detects CREST

The application uses `detect_crest_available()` to check if both binaries are on PATH:

```python
# From src/qm_nmr_calc/conformers/crest_generator.py
def detect_crest_available() -> bool:
    """Detect if CREST and xTB binaries are available on PATH."""
    crest_found = shutil.which("crest") is not None
    xtb_found = shutil.which("xtb") is not None
    return crest_found and xtb_found
```

You can verify detection via the health endpoint:

```bash
curl http://localhost:8000/api/v1/health/ready | jq '.crest_available'
# Should return: true (if installed) or false (if not found)
```

### GFN2-xTB Force Field

The application uses the GFN2-xTB (Geometry, Frequency, Noncovalent interactions, version 2) force field for:

1. **Conformer ranking** during ensemble generation
2. **Energy estimation** for Boltzmann weighting
3. **Structure optimization** within CREST workflow

GFN2-xTB is built into the xTB binary - no additional configuration needed.

## Environment Validation

After installation, validate that all components work correctly.

### Step-by-Step Validation

**1. Check Python version:**

```bash
python3 --version
# Should show: Python 3.11.x or later
```

**2. Check NWChem installation:**

```bash
which nwchem
# Should return a path like: /usr/bin/nwchem

# Try running NWChem (should not error immediately)
echo "start test" | nwchem /dev/stdin 2>&1 | head -5
# Should show NWChem banner or version info
```

**3. Check OpenMPI (optional but recommended):**

```bash
which mpirun
# Should return: /usr/bin/mpirun or similar
```

**4. Start the API server:**

```bash
# Start in background
uv run python scripts/run_api.py &
# Wait a few seconds for startup
sleep 3
```

**5. Check health endpoint:**

```bash
curl http://localhost:8000/api/v1/health/ready
```

**Expected response:**

```json
{
  "status": "ready",
  "checks": {
    "data_directory": "ok",
    "task_queue": "ok",
    "crest_available": true
  },
  "crest_available": true,
  "timestamp": "2026-01-31T22:30:00Z"
}
```

Notes:
- `status` should be `"ready"`
- `data_directory` and `task_queue` should be `"ok"`
- `crest_available` reflects whether CREST/xTB are installed (true or false)

**6. Stop the API server:**

```bash
# If running in background with &
kill %1

# Or use Ctrl+C if running in foreground
```

**7. Optional: Check CREST/xTB:**

```bash
which crest && which xtb
# Should return both paths if installed
# Returns nothing if not installed (which is fine - app still works)

# Check versions
crest --version
xtb --version
```

### Next Steps

If validation succeeds, you're ready to use the calculator:

1. See [Usage Guide](usage.md) for web UI workflow
2. See [API Documentation](usage.md#api-reference) for REST API usage
3. Run your first calculation with the web interface at http://localhost:8000

## Troubleshooting Common Issues

### NWChem not found

**Symptom:**

```
FileNotFoundError: [Errno 2] No such file or directory: 'nwchem'
```

Or job fails at the DFT optimization step with "command not found" errors.

**Fix:**

1. Verify NWChem is installed:
   ```bash
   which nwchem
   ```

2. If not found, install NWChem:
   ```bash
   sudo apt install nwchem  # Ubuntu/Debian
   sudo dnf install nwchem  # Fedora
   ```

3. If installed but not on PATH, add it:
   ```bash
   export PATH="/path/to/nwchem/bin:$PATH"
   # Add to ~/.bashrc or ~/.profile to make permanent
   ```

### RDKit import errors

**Symptom:**

```
ModuleNotFoundError: No module named 'rdkit'
```

**Fix:**

This means you're not using the project's virtual environment.

**Solution 1 (recommended):** Use `uv run`

```bash
uv run python scripts/run_api.py
```

**Solution 2:** Activate the virtual environment

```bash
source .venv/bin/activate
python scripts/run_api.py
```

**Solution 3:** Reinstall dependencies

```bash
cd /path/to/qm-nmr-calc
uv sync
```

### MPI/OpenMPI issues

**Symptom:**

NWChem hangs indefinitely, or fails with MPI-related errors:

```
MPI_ABORT was invoked on rank 0
```

**Fix:**

1. Install OpenMPI:
   ```bash
   sudo apt install openmpi-bin libopenmpi-dev  # Ubuntu
   sudo dnf install openmpi openmpi-devel       # Fedora
   ```

2. Check MPI is working:
   ```bash
   mpirun --version
   ```

3. Verify NWChem can find MPI:
   ```bash
   which mpirun
   ldd $(which nwchem) | grep mpi
   ```

> **Note:** NWChem works without MPI in single-core mode, but is significantly slower. If you consistently see MPI errors, you can run calculations without parallelization - they'll just take longer.

### Permission denied on scratch/data directory

**Symptom:**

```
PermissionError: [Errno 13] Permission denied: './data/jobs/...'
```

**Fix:**

1. Check directory permissions:
   ```bash
   ls -ld ./data/jobs
   ```

2. Ensure the directory is writable:
   ```bash
   mkdir -p ./data/jobs
   chmod 755 ./data/jobs
   ```

3. Check disk space:
   ```bash
   df -h .
   ```

The application creates job directories at `./data/jobs/{job_id}/` relative to the working directory where you start the server.

> **Note:** Each job requires ~10-50 MB of scratch space for intermediate files.

### CREST timeout on large molecules

**Symptom:**

CREST ensemble generation times out after 2 hours:

```
TimeoutError: CREST ensemble generation exceeded 7200 seconds
```

**Why this happens:**

CREST performs extensive conformational searches. For highly flexible molecules (many rotatable bonds), this can take hours.

**Fix options:**

1. **Use RDKit mode instead** (faster, good for most molecules):
   - In the web UI: select "RDKit" conformer method
   - Via API: set `conformer_method: "rdkit_kdg"`

2. **Reduce conformer search space** (not currently configurable):
   - The default timeout is 7200 seconds (2 hours)
   - This is hard-coded to prevent runaway processes

3. **Use single conformer mode** for very large/flexible molecules:
   - Select "Single Conformer" mode in the UI
   - Uses lowest-energy RDKit conformer only

> **Note:** For molecules with >10 rotatable bonds, consider using single conformer mode or RDKit ensemble instead of CREST.

### Memory issues

**Symptom:**

Job fails with:
```
MemoryError
```

Or NWChem crashes with "insufficient memory" messages.

**Memory requirements:**

- **Minimum:** 4 GB RAM for small molecules (<20 atoms)
- **Recommended:** 8 GB RAM for typical organic molecules
- **Production presets:** May need 16 GB for large molecules (>50 atoms)

**Fix:**

1. Use draft preset for testing (less memory intensive):
   - Draft uses 6-31G* basis set
   - Production uses 6-311+G(2d,p) basis set (larger)

2. Close other applications to free memory

3. Monitor memory usage:
   ```bash
   free -h
   htop  # or top
   ```

4. For large molecules, use a machine with more RAM

### Health endpoint shows crest_available: false

**Symptom:**

```json
{
  "crest_available": false
}
```

**This is not an error** - it means CREST/xTB are not installed.

**What this means:**

- The application works fully without CREST
- Ensemble mode will use RDKit instead
- You can still run all calculations

**To enable CREST:**

Follow the [Optional: CREST and xTB Setup](#optional-crest-and-xtb-setup) section above.

### API server won't start - port already in use

**Symptom:**

```
OSError: [Errno 98] Address already in use
```

**Fix:**

1. Check what's using port 8000:
   ```bash
   lsof -i :8000
   ```

2. Kill the existing process:
   ```bash
   kill <PID>
   ```

3. Or use a different port:
   ```bash
   uv run uvicorn qm_nmr_calc.api.main:app --host 0.0.0.0 --port 8001
   ```

## Getting Help

### Health Diagnostics

The health endpoint provides diagnostic information:

```bash
curl http://localhost:8000/api/v1/health/ready | jq
```

Include this output when reporting issues.

### Reporting Issues

If you encounter problems not covered here:

1. Check the [GitHub Issues](https://github.com/steinbeck/qm-nmr-calc/issues)
2. Search for similar issues
3. Create a new issue with:
   - Error message
   - Output from health endpoint
   - Steps to reproduce
   - System info (OS, Python version, NWChem version)

## Related Documentation

- [Usage Guide](usage.md) - How to run calculations
- [Architecture](architecture.md) - System design and data flow
- [NMR Methodology](science.md) - Scientific methodology and DP4+

