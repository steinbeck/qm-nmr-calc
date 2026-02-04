---
phase: 41-arm64-dockerfile
plan: 01
subsystem: docker
tags: [arm64, dockerfile, micromamba, conda-forge]

dependency-graph:
  requires: []
  provides: [arm64-dockerfile, arm64-env-yaml, arm64-validation-script]
  affects: [42-build-publish]

tech-stack:
  added:
    - mambaorg/micromamba:2.5.0-debian12-slim
  patterns:
    - conda-forge-only-packages
    - micromamba-activation-scripts
    - architecture-aware-validation

key-files:
  created:
    - env-worker-arm64.yaml
    - Dockerfile.worker.arm64
    - scripts/validate-worker-arm64.sh
  modified: []

decisions:
  - id: arm64-base-image
    choice: mambaorg/micromamba:2.5.0-debian12-slim
    rationale: All ARM64 packages from conda-forge, no pre-built binaries available
  - id: conda-env-name
    choice: "name: base"
    rationale: micromamba activates base by default, simplifies Dockerfile
  - id: rdkit-source
    choice: conda-forge (not pip)
    rationale: conda handles C++ dependencies correctly for ARM64
  - id: arch-check-behavior
    choice: warning not error
    rationale: allows testing on x86 with QEMU emulation

metrics:
  duration: 3 min
  completed: 2026-02-04
---

# Phase 41 Plan 01: ARM64 Dockerfile Creation Summary

**ARM64 Docker build files using micromamba base and conda-forge packages**

## One-liner

ARM64 Dockerfile with micromamba + conda-forge for NWChem/xTB/CREST/RDKit, replacing x86 binary downloads.

## What Was Built

### 1. Conda Environment File (env-worker-arm64.yaml)

Defines all packages for ARM64 worker container:

```yaml
channels:
  - conda-forge
dependencies:
  - python=3.11
  - nwchem=7.3.1
  - xtb=6.7.1
  - crest=3.0.2
  - rdkit>=2025.9
  - xorg-libxrender  # X11 for RDKit drawing
  - xorg-libxext
  - pip:
    - fastapi, huey, uvicorn, etc. (from pyproject.toml)
```

Key difference from x86: RDKit comes from conda-forge (not pip) because conda handles the C++ dependencies correctly for ARM64.

### 2. ARM64 Dockerfile (Dockerfile.worker.arm64)

Uses micromamba base image instead of NWChem base:

```dockerfile
FROM mambaorg/micromamba:2.5.0-debian12-slim
COPY env-worker-arm64.yaml /tmp/
RUN micromamba install -y -n base -f /tmp/env-worker-arm64.yaml
ARG MAMBA_DOCKERFILE_ACTIVATE=1  # CRITICAL for pip/python
```

Environment variables set:
- `OMP_STACKSIZE="2G"` - Prevent stack overflow on large molecules
- `OMP_NUM_THREADS="4"` - OpenMP threads
- `OPENBLAS_NUM_THREADS="4"` - BLAS threading
- `GFORTRAN_UNBUFFERED_ALL="1"` - I/O buffering
- `NWCHEM_NPROC="4"` - MPI processes

**NOT set** (conda activation handles automatically):
- `NWCHEM_BASIS_LIBRARY` - Set by conda activate script
- `XTBPATH` - Set by conda activate script
- `OMPI_ALLOW_RUN_AS_ROOT` - Not needed (micromamba runs as non-root)

### 3. Validation Script (scripts/validate-worker-arm64.sh)

Architecture-aware validation:

```bash
# Architecture check (warning, not error)
if [ "$(uname -m)" != "aarch64" ]; then
    echo "WARNING: Not running on ARM64"
fi

# NWCHEM_BASIS_LIBRARY check (set by conda)
if [ -n "$NWCHEM_BASIS_LIBRARY" ] && [ -d "$NWCHEM_BASIS_LIBRARY" ]; then
    echo "PASS: Basis library exists"
fi
```

Validates:
- Architecture (aarch64 expected, warning if different)
- NWChem binary and basis library
- xTB binary (--version)
- CREST binary (--version)
- Python 3.11
- RDKit import and SMILES parsing
- qm_nmr_calc.queue.huey import
- qm_nmr_calc.tasks.run_nmr_task import
- Environment variables

## Key Differences from x86 Dockerfile

| Aspect | x86 (Dockerfile.worker) | ARM64 (Dockerfile.worker.arm64) |
|--------|-------------------------|----------------------------------|
| Base image | ghcr.io/nwchemgit/nwchem-dev/amd64 | mambaorg/micromamba:2.5.0-debian12-slim |
| NWChem | Pre-installed in base | conda install nwchem |
| xTB | wget GitHub release binary | conda install xtb |
| CREST | wget GitHub release binary | conda install crest |
| RDKit | pip install | conda install rdkit |
| NWCHEM_BASIS_LIBRARY | Hardcoded path | Conda activation script |
| USER root | Required for apt/wget | Not needed |
| ENTRYPOINT reset | Required | Not needed |

## Commits

| Hash | Description |
|------|-------------|
| 685c80c | feat(41-01): create ARM64 conda environment file |
| cf3aa44 | feat(41-01): create ARM64 Dockerfile using micromamba |
| 92bbd63 | feat(41-01): create ARM64 validation script |

## Files Created

```
env-worker-arm64.yaml          # Conda environment specification
Dockerfile.worker.arm64        # ARM64 container definition
scripts/validate-worker-arm64.sh  # Validation script
```

## Deviations from Plan

None - plan executed exactly as written.

## Verification Checklist

- [x] env-worker-arm64.yaml exists with nwchem, xtb, crest, rdkit from conda-forge
- [x] Dockerfile.worker.arm64 uses mambaorg/micromamba base image
- [x] Dockerfile.worker.arm64 includes MAMBA_DOCKERFILE_ACTIVATE=1 ARG
- [x] Dockerfile.worker.arm64 sets OMP_STACKSIZE, OMP_NUM_THREADS, OPENBLAS_NUM_THREADS
- [x] scripts/validate-worker-arm64.sh includes aarch64 architecture check
- [x] All files pass syntax validation

## Next Phase Readiness

Phase 42 (Build and Publish) can proceed:
- Dockerfile.worker.arm64 ready for building
- Validation script ready to test built image
- Phase 42 will need ARM64 runner (GitHub Actions ubuntu-24.04-arm or Mac builder)
