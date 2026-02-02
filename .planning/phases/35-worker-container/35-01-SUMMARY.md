# Phase 35-01: Worker Container - Summary

**Completed:** 2026-02-02
**Duration:** ~45 minutes
**Status:** Complete

## Objective

Create and validate the worker Docker container with NWChem, CREST, xTB, and Huey consumer.

## Deliverables

| File | Purpose |
|------|---------|
| `Dockerfile.worker` | Multi-stage Docker build for worker container |
| `scripts/validate-worker.sh` | Validation script for all components |
| `scripts/test-nwchem-container.sh` | Full NWChem DFT calculation test |

## Technical Details

### Base Image
- `ghcr.io/nwchemgit/nwchem-dev/amd64:latest` (Ubuntu 20.04)
- Pre-built NWChem with OpenMPI and OpenBLAS

### Python Environment
- Miniconda Python 3.11 (glibc-compatible with Ubuntu 20.04)
- Virtual environment with all dependencies from pyproject.toml
- qm-nmr-calc package installed via pip

### Scientific Tools
| Tool | Version | Source |
|------|---------|--------|
| NWChem | 7.2.0 | Base image |
| xTB | 6.7.1 | GitHub release binary |
| CREST | 3.0.2 | GitHub release binary |
| RDKit | 2025.9.4 | PyPI |

### Critical Environment Variables
```dockerfile
# OpenMP for CREST/xTB (prevents stack overflow)
OMP_STACKSIZE="2G"
OMP_NUM_THREADS="4"

# Allow MPI to run as root in container
OMPI_ALLOW_RUN_AS_ROOT="1"
OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1"

# xTB data path
XTBPATH="/opt/xtb/share/xtb"

# Python environment
VIRTUAL_ENV="/app/.venv"
```

### Critical Runtime Requirement
**MPI requires shared memory:** Must run with `--shm-size=512m` (minimum) for NWChem calculations.

### Image Size
- Final image: ~2.1GB
- Base NWChem: ~600MB
- Miniconda + packages: ~600MB
- xTB/CREST binaries: ~100MB

## Verification Results

### Build Test
```bash
docker build -f Dockerfile.worker -t qm-nmr-calc-worker:test .
# Success
```

### Component Validation
```bash
docker run --rm qm-nmr-calc-worker:test /app/scripts/validate-worker.sh
# All validations passed: NWChem, xTB, CREST, Python, Huey, RDKit
```

### NWChem DFT Test
```bash
docker run --rm --shm-size=512m qm-nmr-calc-worker:test /app/scripts/test-nwchem-container.sh
# Water molecule optimization completed
# Total DFT energy = -75.844986299864
# Geometry optimization converged
```

### Huey Consumer Test
```bash
docker run --rm -d --name test qm-nmr-calc-worker:test
# Huey consumer started with 1 process
# Tasks registered: run_optimization_task, run_nmr_task, run_ensemble_nmr_task
```

## Challenges Resolved

1. **glibc Incompatibility:** UV's Python 3.12 (Debian 12) was incompatible with NWChem's Ubuntu 20.04 base. Solved by using Miniconda Python 3.11.

2. **NWChem Entrypoint:** Base image has NWChem as entrypoint. Reset with `ENTRYPOINT []`.

3. **MPI Root Restriction:** MPI refuses to run as root by default. Added `OMPI_ALLOW_RUN_AS_ROOT` env vars.

4. **RDKit X11 Libraries:** RDKit drawing requires libXrender1 and libxext6.

5. **CREST Path:** Archive extracts to nested directory. Path is `/opt/crest/crest/crest`.

## Requirements Coverage

| Requirement | Status |
|-------------|--------|
| DOCK-02: Worker container includes NWChem, CREST, and xTB | ✅ Complete |

## Next Steps

Phase 36: API Container - Create Dockerfile.api for FastAPI server with Python dependencies.
