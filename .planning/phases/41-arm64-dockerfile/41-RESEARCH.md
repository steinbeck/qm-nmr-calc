# Phase 41: ARM64 Dockerfile Creation - Research

**Researched:** 2026-02-04
**Domain:** Docker containerization for ARM64 computational chemistry (micromamba + conda-forge)
**Confidence:** HIGH

## Summary

This phase creates a Dockerfile for ARM64 (linux-aarch64) that replaces the x86-specific approach of downloading pre-built binaries. The key insight is that ALL required computational chemistry packages (NWChem, xTB, CREST, RDKit) are now available on conda-forge for linux-aarch64, making the `mambaorg/micromamba` base image the ideal foundation.

The micromamba approach differs fundamentally from the x86 worker: instead of layering binaries on an NWChem base image, we use a clean micromamba image and install everything via conda-forge. This is actually simpler and more maintainable than the x86 approach because environment variables are automatically configured by conda activation scripts.

**Primary recommendation:** Use `mambaorg/micromamba:2.5.0-debian12-slim` as base, create an `environment.yaml` file specifying all packages from conda-forge (nwchem, xtb, crest, rdkit, python, plus pip dependencies), and leverage conda's activation scripts for automatic environment variable configuration.

## Standard Stack

### Core

| Package | Version | Channel | Purpose | ARM64 Support |
|---------|---------|---------|---------|---------------|
| micromamba | 2.5.0 | base image | Conda package manager | Multi-arch (amd64/arm64) |
| nwchem | 7.3.1 | conda-forge | DFT calculations | linux-aarch64 verified |
| xtb | 6.7.1 | conda-forge | Semiempirical for CREST | linux-aarch64 verified |
| crest | 3.0.2 | conda-forge | Conformer search | linux-aarch64 verified |
| rdkit | 2025.09.5 | conda-forge | Cheminformatics | linux-aarch64 verified |
| python | 3.11 | conda-forge | Runtime | linux-aarch64 verified |

### Supporting (via conda-forge or pip)

| Package | Source | Purpose | Notes |
|---------|--------|---------|-------|
| openblas | conda-forge | Linear algebra (NWChem dep) | Auto-installed with nwchem |
| openmpi | conda-forge | MPI parallelism (NWChem dep) | Auto-installed with nwchem |
| fastapi | pip | API framework | From pyproject.toml |
| huey | pip | Task queue | From pyproject.toml |
| matplotlib | conda-forge/pip | Plotting | May need X11 libs |

### Alternatives Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| micromamba debian12-slim | debian13-slim (latest) | debian12 more stable, matches existing infrastructure |
| conda-forge packages | Build from source | Source builds take hours, complex dependencies |
| Single environment.yaml | Multi-stage build | Simpler approach, micromamba handles caching well |

**Installation (environment.yaml):**
```yaml
name: qm-nmr-calc
channels:
  - conda-forge
dependencies:
  - python=3.11
  - nwchem=7.3.1
  - xtb=6.7.1
  - crest=3.0.2
  - rdkit>=2025.9
  - pip
  - pip:
    - qm-nmr-calc  # or install from source
```

## Architecture Patterns

### Dockerfile Structure

```dockerfile
# ARM64 Worker Container for qm-nmr-calc
# Uses conda-forge packages (no x86 binaries)

FROM mambaorg/micromamba:2.5.0-debian12-slim

# Labels
LABEL org.opencontainers.image.source="https://github.com/steinbeck/qm-nmr-calc"
LABEL org.opencontainers.image.description="QM NMR Calculator - ARM64 Worker"

# Copy environment definition
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yaml /tmp/environment.yaml

# Install conda packages (includes nwchem, xtb, crest, rdkit)
RUN micromamba install -y -n base -f /tmp/environment.yaml && \
    micromamba clean --all --yes

# Enable conda activation for subsequent RUN commands
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Copy and install application
WORKDIR /app
COPY --chown=$MAMBA_USER:$MAMBA_USER pyproject.toml README.md ./
COPY --chown=$MAMBA_USER:$MAMBA_USER src/ ./src/
RUN pip install --no-cache-dir .

# Environment variables (supplements conda activation)
ENV OMP_STACKSIZE="2G"
ENV OMP_NUM_THREADS="4"
ENV OPENBLAS_NUM_THREADS="4"
ENV GFORTRAN_UNBUFFERED_ALL="1"
ENV NWCHEM_NPROC="4"

# Data directory
RUN mkdir -p /app/data

# Scripts
COPY --chown=$MAMBA_USER:$MAMBA_USER scripts/validate-worker.sh /app/scripts/
RUN chmod +x /app/scripts/validate-worker.sh

WORKDIR /app
CMD ["huey_consumer", "qm_nmr_calc.queue.huey", "-w", "1", "-k", "process"]
```

### Environment Variables (Required)

| Variable | Value | Source | Purpose |
|----------|-------|--------|---------|
| `NWCHEM_BASIS_LIBRARY` | `$CONDA_PREFIX/share/nwchem/libraries/` | Conda activation | Basis set location |
| `NWCHEM_NWPW_LIBRARY` | `$CONDA_PREFIX/share/nwchem/libraryps/` | Conda activation | Pseudopotential location |
| `XTBPATH` | `$CONDA_PREFIX/share/xtb/` | Conda activation | xTB parameter files |
| `OMP_STACKSIZE` | `2G` | Dockerfile ENV | Prevent stack overflow |
| `OMP_NUM_THREADS` | `4` | Dockerfile ENV | OpenMP threads |
| `OPENBLAS_NUM_THREADS` | `4` | Dockerfile ENV | BLAS threading |
| `GFORTRAN_UNBUFFERED_ALL` | `1` | Dockerfile ENV | I/O buffering |
| `NWCHEM_NPROC` | `4` | Dockerfile ENV | MPI processes |

**Critical insight:** The conda-forge nwchem package includes activation scripts that automatically set `NWCHEM_BASIS_LIBRARY` and `NWCHEM_NWPW_LIBRARY` when the environment is activated. This is handled by micromamba's entrypoint, so manual ENV statements for these are NOT needed (unlike the x86 Dockerfile).

### Anti-Patterns to Avoid

- **Downloading x86 binaries:** GitHub release binaries are x86_64 only. Always use conda-forge for ARM64.
- **Using ghcr.io/nwchemgit base image:** The NWChem Docker images are x86-only. Use micromamba + conda-forge.
- **Skipping MAMBA_DOCKERFILE_ACTIVATE:** Required for pip and python commands in RUN statements.
- **Setting NWCHEM_BASIS_LIBRARY manually:** Conda activation handles this; manual paths may conflict.
- **Using openmpi from apt:** The conda-forge openmpi is tuned for the conda nwchem; mixing can cause issues.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| NWChem for ARM64 | Cross-compile or emulate | `conda install nwchem` | 30-60 min compile, complex MPI/BLAS deps |
| xTB for ARM64 | Build from Fortran source | `conda install xtb` | Requires Fortran toolchain |
| CREST for ARM64 | Build from source | `conda install crest` | Complex dependencies on xTB |
| RDKit for ARM64 | Build from C++ source | `conda install rdkit` | 60+ minute build, boost deps |
| OpenBLAS tuning | Manual configuration | Conda defaults | Conda packages are pre-tuned |

**Key insight:** Conda-forge has done the hard work of building and testing all scientific packages for ARM64. The packages are tested together for compatibility.

## Common Pitfalls

### Pitfall 1: Emulating x86 During Build

**What goes wrong:** Docker build hangs or takes hours when building ARM64 image on x86 host with QEMU emulation
**Why it happens:** micromamba install with QEMU emulation has a ~20% hang rate due to QEMU bugs
**How to avoid:** Build on native ARM64 host (Mac M1/M2, AWS Graviton, ARM64 CI runner)
**Warning signs:** Build hangs at "Solving environment", extremely slow package downloads

### Pitfall 2: Missing MAMBA_DOCKERFILE_ACTIVATE

**What goes wrong:** `pip install` or `python` commands fail with "command not found"
**Why it happens:** Conda environment not activated during build by default
**How to avoid:** Add `ARG MAMBA_DOCKERFILE_ACTIVATE=1` before RUN commands that need conda
**Warning signs:** "pip: command not found", "python: command not found" during build

### Pitfall 3: Wrong NWCHEM_BASIS_LIBRARY Path

**What goes wrong:** NWChem errors: "could not find basis set library"
**Why it happens:** Hardcoding x86 paths or build-time paths instead of conda paths
**How to avoid:** Rely on conda activation scripts, or use `$CONDA_PREFIX/share/nwchem/libraries/`
**Warning signs:** "Error reading basis set" in NWChem output

### Pitfall 4: OpenBLAS Threading Conflict

**What goes wrong:** NWChem hangs or uses unexpected CPU count
**Why it happens:** OpenBLAS compiled with OpenMP conflicts with NWChem's MPI
**How to avoid:** Set `OPENBLAS_NUM_THREADS=4` (or match NWCHEM_NPROC); conda-forge openblas uses pthreads by default on Linux
**Warning signs:** 100% CPU on single core, or unexpected multi-threading

### Pitfall 5: Stack Overflow on Large Molecules

**What goes wrong:** CREST/xTB crashes with SIGSEGV on molecules with 50+ atoms
**Why it happens:** Default stack size (8MB) insufficient for deep recursion
**How to avoid:** Set `OMP_STACKSIZE=2G` in Dockerfile ENV
**Warning signs:** Works on small molecules, crashes on large ones

### Pitfall 6: X11 Libraries for RDKit Drawing

**What goes wrong:** RDKit mol.Draw() fails with X11 errors
**Why it happens:** Headless container missing libxrender, libxext
**How to avoid:** Install X11 libs: `apt-get install libxrender1 libxext6`
**Warning signs:** "cannot open display", Xlib errors in logs

## Code Examples

### environment.yaml (Complete)

```yaml
# Source: ARM64 worker conda environment
name: base
channels:
  - conda-forge
dependencies:
  # Core computational chemistry
  - python=3.11
  - nwchem=7.3.1
  - xtb=6.7.1
  - crest=3.0.2
  - rdkit>=2025.9

  # X11 for RDKit drawing (headless)
  - xorg-libxrender
  - xorg-libxext

  # pip for remaining dependencies
  - pip
  - pip:
    - fastapi>=0.128.0
    - huey>=2.6.0
    - uvicorn>=0.40.0
    - aiosmtplib>=5.0.0
    - email-validator>=2.3.0
    - jinja2>=3.1.0
    - matplotlib>=3.10.0
    - orjson>=3.11.5
    - pandas>=2.3.3
    - pydantic>=2.12.5
    - python-multipart>=0.0.21
    - scipy>=1.17.0
    - statsmodels>=0.14.0
    - tqdm>=4.67.1
```

### Validation Script

```bash
#!/bin/bash
# validate-worker-arm64.sh
# Verify ARM64 worker container has all components

set -e

echo "=== ARM64 Worker Validation ==="

# Check architecture
echo "Architecture: $(uname -m)"
if [ "$(uname -m)" != "aarch64" ]; then
    echo "WARNING: Not running on ARM64"
fi

# Check NWChem
echo -n "NWChem: "
which nwchem && nwchem --version 2>&1 | head -1 || echo "NWChem returns exit 1 (normal)"
echo "NWCHEM_BASIS_LIBRARY: $NWCHEM_BASIS_LIBRARY"
ls "$NWCHEM_BASIS_LIBRARY" | head -3

# Check xTB
echo -n "xTB: "
xtb --version 2>&1 | head -1

# Check CREST
echo -n "CREST: "
crest --version 2>&1 | head -1

# Check Python and packages
echo -n "Python: "
python --version
python -c "from rdkit import Chem; print('RDKit:', Chem.rdBase.rdkitVersion)"
python -c "from qm_nmr_calc.queue import huey; print('Huey import: OK')"

echo "=== Validation Complete ==="
```

### NWChem Test (Verify Environment)

```python
# Test that NWChem basis sets are accessible
import os
import subprocess

# Conda activation should set this
basis_lib = os.environ.get("NWCHEM_BASIS_LIBRARY")
print(f"NWCHEM_BASIS_LIBRARY: {basis_lib}")

# Verify library exists
if basis_lib and os.path.isdir(basis_lib):
    files = os.listdir(basis_lib)
    print(f"Basis sets available: {len(files)}")
    assert "6-31g*" in files or "6-31g_st_" in files
else:
    raise RuntimeError("NWCHEM_BASIS_LIBRARY not configured")
```

## State of the Art

| Old Approach (x86) | New Approach (ARM64) | Why Different |
|--------------------|----------------------|---------------|
| NWChem official Docker image | micromamba + conda-forge nwchem | No ARM64 official image |
| xTB GitHub release binary | conda-forge xtb | No ARM64 binary releases |
| CREST GitHub release binary | conda-forge crest | No ARM64 binary releases |
| Miniconda + pip for Python | micromamba + conda + pip | Faster, smaller, better caching |
| Manual ENV for NWCHEM_BASIS_LIBRARY | Conda activation scripts | Automatic, less error-prone |

**Key ecosystem change:** As of 2024-2025, conda-forge has comprehensive ARM64 (linux-aarch64) support for computational chemistry. This eliminates the need for source compilation or architecture-specific workarounds.

## Open Questions

1. **MPI Configuration for ARM64**
   - What we know: conda-forge nwchem includes openmpi, `--bind-to none` should work
   - What's unclear: Whether any ARM64-specific MPI tuning is needed
   - Recommendation: Test with existing mpirun flags first, tune if issues arise

2. **shm_size for ARM64**
   - What we know: x86 worker needs `shm_size: 512m` for MPI
   - What's unclear: Whether ARM64 has same requirements
   - Recommendation: Keep same setting in docker-compose, test to verify

3. **Performance Comparison**
   - What we know: ARM64 (Apple M-series, Graviton) is competitive for scientific workloads
   - What's unclear: Exact performance vs x86 for NWChem DFT calculations
   - Recommendation: Document baseline performance once container is running

## Sources

### Primary (HIGH confidence)
- [conda-forge nwchem package](https://anaconda.org/conda-forge/nwchem/files) - Verified 7.3.1 for linux-aarch64
- [conda-forge xtb package](https://anaconda.org/conda-forge/xtb/files) - Verified 6.7.1 for linux-aarch64
- [conda-forge crest package](https://anaconda.org/conda-forge/crest/files) - Verified 3.0.2 for linux-aarch64
- [conda-forge rdkit package](https://anaconda.org/conda-forge/rdkit/files) - Verified 2025.09.5 for linux-aarch64
- [micromamba-docker documentation](https://micromamba-docker.readthedocs.io/en/latest/quick_start.html) - Dockerfile patterns
- [nwchem-feedstock NWCHEM_BASIS_LIBRARY issue](https://github.com/conda-forge/nwchem-feedstock/issues/6) - Activation scripts solution

### Secondary (MEDIUM confidence)
- [OpenBLAS GitHub USAGE.md](https://github.com/OpenMathLib/OpenBLAS/blob/develop/USAGE.md) - Threading environment variables
- [NWChem Compiling docs](https://nwchemgit.github.io/Compiling-NWChem.html) - OpenBLAS configuration
- [mamba-org/micromamba-docker releases](https://github.com/mamba-org/micromamba-docker/releases) - v2.5.0 ARM64 support confirmed

### Codebase (HIGH confidence - verified)
- `Dockerfile.worker` - Existing x86 patterns to adapt
- `docker-compose.yml` - shm_size, environment variables
- `pyproject.toml` - Python dependencies to mirror
- `src/qm_nmr_calc/conformers/crest_generator.py` - OMP_STACKSIZE requirements

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - All packages verified on conda-forge for linux-aarch64
- Architecture: HIGH - Based on micromamba official docs and existing x86 patterns
- Pitfalls: MEDIUM - ARM64-specific issues extrapolated from x86 experience and docs

**Research date:** 2026-02-04
**Valid until:** 2026-03-04 (30 days - stable domain, conda-forge packages update regularly)
