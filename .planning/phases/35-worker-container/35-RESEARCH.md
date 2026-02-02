# Phase 35: Worker Container - Research

**Researched:** 2026-02-02
**Domain:** Docker containerization for scientific computing (NWChem, CREST, xTB)
**Confidence:** HIGH

## Summary

This phase focuses on containerizing the Huey worker with all computational dependencies: NWChem (DFT calculations), CREST (conformer search), and xTB (semiempirical methods). The research confirmed that the official NWChem Docker image (`ghcr.io/nwchemgit/nwchem-dev/amd64`) is the correct base, providing a pre-built NWChem with OpenMPI. CREST and xTB can be installed from pre-built GitHub release binaries. Python and the application are installed on top using a multi-stage build pattern with uv.

The critical insight is that NWChem in containers requires special MPI configuration: `shm_size: 512m` (minimum 256m) to avoid segfaults, and the codebase already uses `--bind-to none` in mpirun which is container-compatible. CREST/xTB require `OMP_STACKSIZE=2G` (already set in the codebase) to handle large molecules. The worker container will be approximately 3-4GB due to NWChem and scientific libraries.

**Primary recommendation:** Use `ghcr.io/nwchemgit/nwchem-dev/amd64:latest` as base, install CREST/xTB from GitHub release binaries, and Python/app via uv. Set all environment variables in Dockerfile. Test with validation commands for each component.

## Standard Stack

### Core

| Component | Version | Purpose | Why Standard |
|-----------|---------|---------|--------------|
| `ghcr.io/nwchemgit/nwchem-dev/amd64` | latest | NWChem base image | Official image with OpenMPI, BLAS, tested configuration |
| xTB | 6.7.1 | Semiempirical calculations for CREST | Latest stable release, pre-built Linux binary |
| CREST | latest | Conformer search | Continuous release with GFN2-xTB support |
| Python | 3.12 | Application runtime | Matches project pyproject.toml |
| uv | latest | Python package manager | Fast installs, good Docker caching |

### Supporting

| Component | Purpose | When to Use |
|-----------|---------|-------------|
| OpenMPI | MPI parallelism for NWChem | Pre-installed in NWChem image |
| OpenBLAS | Linear algebra | Pre-installed in NWChem image |
| peewee | SQLite ORM for Huey | Bundled with application |

### Alternatives Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| NWChem dev image | Building from source | Source build takes 30-60 min, complex dependencies |
| xTB binary | conda-forge | Conda adds 400MB+, binary is self-contained |
| Python system packages | pip only | uv provides faster, more reliable installs |

**Installation commands:**
```bash
# xTB (in Dockerfile)
wget https://github.com/grimme-lab/xtb/releases/download/v6.7.1/xtb-6.7.1-linux-x86_64.tar.xz
tar xf xtb-6.7.1-linux-x86_64.tar.xz -C /opt

# CREST (in Dockerfile)
wget https://github.com/grimme-lab/crest/releases/download/latest/crest-gnu-12-ubuntu-latest.tar.xz
tar xf crest-gnu-12-ubuntu-latest.tar.xz -C /opt/crest
```

## Architecture Patterns

### Dockerfile Structure

```dockerfile
# Stage 1: Build Python dependencies
FROM ghcr.io/astral-sh/uv:python3.12-bookworm AS python-builder
# Install app dependencies into virtual environment
WORKDIR /app
COPY pyproject.toml uv.lock ./
RUN uv sync --frozen --no-dev

# Stage 2: Final image with NWChem + CREST + xTB + Python
FROM ghcr.io/nwchemgit/nwchem-dev/amd64:latest
# Install xTB and CREST binaries
# Copy Python environment from builder
# Set all environment variables
# Run huey_consumer
```

### Environment Variables (Required)

```dockerfile
# NWChem environment (from base image, verify these exist)
ENV NWCHEM_TOP="/opt/nwchem"
ENV NWCHEM_BASIS_LIBRARY="${NWCHEM_TOP}/src/basis/libraries/"

# xTB environment
ENV XTBPATH="/opt/xtb/share/xtb"
ENV PATH="/opt/xtb/bin:/opt/crest:${PATH}"

# OpenMP for CREST/xTB (CRITICAL - prevents stack overflow)
ENV OMP_STACKSIZE="2G"
ENV OMP_NUM_THREADS="4"
ENV GFORTRAN_UNBUFFERED_ALL="1"

# Python environment
ENV VIRTUAL_ENV="/app/.venv"
ENV PATH="${VIRTUAL_ENV}/bin:${PATH}"

# Application paths
ENV NWCHEM_NPROC="4"
WORKDIR /app
```

### Anti-Patterns to Avoid

- **Building NWChem from source:** Takes 30-60 minutes, complex MPI/BLAS dependencies. Use official image.
- **Installing xTB/CREST via conda:** Adds 400MB+, version conflicts. Use binaries.
- **Running as root:** Security risk. Create non-root user.
- **Hardcoded paths:** Use environment variables for all data paths.
- **Missing shm_size in compose:** NWChem will segfault. Must set in docker-compose.yaml.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| NWChem installation | Source compilation | Official Docker image | MPI/BLAS configuration is complex |
| xTB/CREST installation | Source compilation | Pre-built binaries | Static binaries work out of the box |
| Python environment isolation | Manual venv in Dockerfile | uv multi-stage build | Better caching, faster builds |
| MPI process management | Custom scripts | mpirun with --bind-to none | Already tested in codebase |

**Key insight:** Scientific software compilation is extremely fragile due to Fortran/MPI/BLAS interdependencies. Pre-built binaries and official images encapsulate years of debugging.

## Common Pitfalls

### Pitfall 1: Insufficient Shared Memory (shm_size)
**What goes wrong:** NWChem segfaults immediately or MPI hangs
**Why it happens:** Docker default /dev/shm is 64MB; NWChem MPI needs 256-512MB minimum
**How to avoid:** Set `shm_size: '512m'` in docker-compose.yaml worker service
**Warning signs:** "Bus error", "Segmentation fault", MPI hangs at startup

### Pitfall 2: Stack Overflow in CREST/xTB
**What goes wrong:** SIGSEGV when processing large molecules (50+ atoms)
**Why it happens:** Default Linux stack (8MB) is insufficient for deep recursion
**How to avoid:** Set `OMP_STACKSIZE=2G` in Dockerfile ENV
**Warning signs:** Crashes only on large molecules, works on small ones

### Pitfall 3: MPI CPU Binding Fails in Container
**What goes wrong:** "not enough slots available" or wrong CPU count detection
**Why it happens:** OpenMPI tries to bind to host CPUs visible through /sys/devices
**How to avoid:** Use `--bind-to none` in mpirun command (already in codebase)
**Warning signs:** Works on host, fails in container; CPU count mismatch errors

### Pitfall 4: Missing Environment Variables
**What goes wrong:** "nwchem not found", "xtb not found", wrong basis set path
**Why it happens:** PATH and library paths not set in Dockerfile
**How to avoid:** Explicitly set all ENV variables in Dockerfile
**Warning signs:** "command not found", "basis set not found" errors

### Pitfall 5: Wrong Working Directory for Scratch
**What goes wrong:** "Permission denied" or "No space left on device"
**Why it happens:** NWChem writes scratch to current directory by default
**How to avoid:** Use job_dir/scratch paths (already in codebase), ensure volume mounted
**Warning signs:** Permission errors in logs, disk space warnings

### Pitfall 6: Huey Database Path Mismatch
**What goes wrong:** Worker can't find queued jobs, "database is locked"
**Why it happens:** Huey uses `./data/huey.db` which changes based on WORKDIR
**How to avoid:** Mount data volume to /app/data, ensure WORKDIR is /app
**Warning signs:** Empty queue, sqlite lock errors

## Code Examples

### NWChem Invocation (Existing Code - runner.py line 77)
```python
# Source: src/qm_nmr_calc/nwchem/runner.py
cmd = f"unset DISPLAY; mpirun --bind-to none -n {processes} nwchem {infile} > {outfile} 2> {logfile}"
```
The `--bind-to none` is critical for container compatibility.

### CREST Invocation (Existing Code - crest_generator.py line 210-227)
```python
# Source: src/qm_nmr_calc/conformers/crest_generator.py
cmd = [
    "crest",
    str(input_xyz),
    "--gfn2",
    "--alpb",
    solvent,
    "--chrg",
    str(charge),
    "--ewin",
    str(ewin_kcal),
    "-T",
    str(num_threads),
]

env = os.environ.copy()
env["OMP_STACKSIZE"] = "2G"
env["GFORTRAN_UNBUFFERED_ALL"] = "1"
```
Environment variables already set in code, but should also be in Dockerfile.

### Huey Consumer Startup (queue.py)
```python
# Source: src/qm_nmr_calc/queue.py
huey = SqliteHuey(
    'qm-nmr-calc',
    filename='./data/huey.db',
    fsync=True
)
```
Consumer command: `huey_consumer qm_nmr_calc.tasks.huey`

### Validation Commands for Container
```bash
# Verify NWChem
nwchem --version 2>&1 | head -5 || echo "NWChem returns exit 1 on --version but works"
which nwchem && echo "NWChem found"

# Verify xTB
xtb --version

# Verify CREST
crest --version

# Verify Python environment
python --version
python -c "from qm_nmr_calc.tasks import huey; print('Huey import OK')"
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Build NWChem in Dockerfile | Use official ghcr.io image | 2023 | 30-60 min build time saved |
| Conda for xTB/CREST | GitHub release binaries | 2024 | 400MB smaller, no version conflicts |
| pip install | uv sync | 2024 | 10x faster Python installs |
| Single-stage Dockerfile | Multi-stage build | Standard | Smaller final image |

**Current versions (as of research date):**
- NWChem: 7.3.1 (latest release), dev image tracks master
- xTB: 6.7.1 (July 2024)
- CREST: Continuous release (November 2025)

## Open Questions

1. **arm64 Support**
   - What we know: NWChem has arm64 image (`ghcr.io/nwchemgit/nwchem-dev/arm64`)
   - What's unclear: Whether CREST/xTB have arm64 binaries
   - Recommendation: Focus on amd64 first, defer arm64 to v2.5

2. **Optimal NWCHEM_NPROC Default**
   - What we know: Codebase uses `processes=4` default
   - What's unclear: Whether this should be configurable via environment
   - Recommendation: Allow override via `NWCHEM_NPROC` env var, default 4

3. **Image Size Optimization**
   - What we know: NWChem base is ~2GB, final will be 3-4GB
   - What's unclear: Whether further reduction is possible
   - Recommendation: Accept size for v2.4, optimize later if needed

## Sources

### Primary (HIGH confidence)
- [NWChem Containers Documentation](https://nwchemgit.github.io/Containers.html) - Official Docker usage, shm_size requirements
- [NWChem Dockerfiles Repository](https://github.com/nwchemgit/nwchem-dockerfiles) - Dockerfile patterns, environment variables
- [xTB GitHub Releases](https://github.com/grimme-lab/xtb/releases) - v6.7.1 binary downloads
- [CREST GitHub Releases](https://github.com/grimme-lab/crest/releases) - Latest binary downloads
- [xTB Documentation - Setup](https://xtb-docs.readthedocs.io/en/latest/setup.html) - Environment variables
- [Huey Consumer Documentation](https://huey.readthedocs.io/en/latest/consumer.html) - Consumer command options

### Secondary (MEDIUM confidence)
- [CREST Installation Docs](https://crest-lab.github.io/crest-docs/page/installation/install_basic.html) - Installation steps, PATH setup
- [uv Docker Guide](https://docs.astral.sh/uv/guides/integration/docker/) - Multi-stage build patterns

### Codebase (HIGH confidence - verified)
- `src/qm_nmr_calc/nwchem/runner.py` - NWChem invocation with --bind-to none
- `src/qm_nmr_calc/conformers/crest_generator.py` - CREST invocation, environment setup
- `src/qm_nmr_calc/queue.py` - Huey configuration, database path
- `src/qm_nmr_calc/tasks.py` - Task definitions, consumer entry point

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - Official images and releases verified
- Architecture: HIGH - Based on existing codebase patterns and official docs
- Pitfalls: HIGH - Verified from NWChem docs, xTB issues, and codebase

**Research date:** 2026-02-02
**Valid until:** 2026-03-02 (30 days - stable domain, versions may update)
