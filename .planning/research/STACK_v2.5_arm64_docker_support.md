# Technology Stack: ARM64 Worker Container

**Project:** qm-nmr-calc v2.5 - ARM64 Docker Support
**Researched:** 2026-02-04
**Confidence:** HIGH (verified with conda-forge, Docker Hub, and GitHub documentation)

## Executive Summary

ARM64 worker container is **fully feasible** using conda-forge packages. All required computational chemistry packages (NWChem, xTB, CREST, RDKit) are available for linux-aarch64. The recommended approach is a **separate ARM64-native Dockerfile** using `mambaorg/micromamba` as the base image, with native GitHub Actions runners for CI/CD builds.

**Key insight:** The current Dockerfile.worker comment "arm64 not supported - CREST/xTB lack arm64 binaries" is now **outdated**. conda-forge has provided linux-aarch64 builds for all required packages since late 2025.

## Recommended Stack

### Base Image

| Technology | Version | Platform Support | Why |
|------------|---------|------------------|-----|
| mambaorg/micromamba | 2.5.0 | linux/amd64, linux/arm64, linux/ppc64le | Fast package resolution (2x faster than miniconda), multi-arch support, smaller image size, actively maintained |

**Specific tag recommendation:** `mambaorg/micromamba:2.5.0-debian12-slim`

**Verified platform support (via `docker manifest inspect`):**
```
- linux/amd64
- linux/arm64
- linux/ppc64le
```

Rationale:
- Debian base has broader compatibility than Alpine for scientific packages
- `-slim` variant reduces image size while retaining glibc
- Pinned version (2.5.0) ensures reproducible builds
- Multi-arch manifest means same tag works for both architectures
- Alpine has known issues: musl libc incompatibility with some conda-forge packages; QEMU bug can hang `micromamba install` when emulating different CPU

### Core Computational Packages (conda-forge)

| Package | Version | ARM64 Status | Verified Date | Source |
|---------|---------|--------------|---------------|--------|
| nwchem | 7.3.1 | AVAILABLE (linux-aarch64) | 2026-02-04 | [conda-forge](https://anaconda.org/conda-forge/nwchem) |
| xtb | 6.7.1 | AVAILABLE (hfc05e43_4, 5.05 MB) | 2026-02-04 | [conda-forge](https://anaconda.org/conda-forge/xtb) |
| crest | 3.0.2 | AVAILABLE (Aug 5, 2025 build) | 2026-02-04 | [conda-forge](https://anaconda.org/conda-forge/crest) |
| rdkit | 2025.09.5 | AVAILABLE (py3.10-3.14) | 2026-02-04 | [conda-forge](https://anaconda.org/conda-forge/rdkit) |
| python | 3.11.x | AVAILABLE | N/A | Standard |

**All packages verified on [conda-forge](https://anaconda.org/conda-forge) as of 2026-02-04.**

### Python Dependencies Strategy

| Source | Packages | Rationale |
|--------|----------|-----------|
| conda-forge | rdkit, matplotlib, pandas, scipy, statsmodels | Heavy scientific packages with native ARM64 builds |
| pip (in conda env) | huey, jinja2, orjson, pydantic, tqdm, aiosmtplib, email-validator | Pure-Python packages, no native compilation needed |

**Recommendation:** Install heavy scientific packages via conda-forge for native ARM64 binaries. Install pure-Python packages via pip within the conda environment.

### CI/CD Infrastructure

| Component | Specification | Why |
|-----------|---------------|-----|
| x86_64 Runner | `ubuntu-latest` | Standard, fast |
| ARM64 Runner | `ubuntu-24.04-arm` | Native ARM64, FREE for public repos (GA since Aug 2025) |
| Build Tool | `docker buildx` | Multi-platform support |
| Registry | GHCR (ghcr.io) | Already in use, supports multi-arch manifests |

**GitHub Actions ARM64 runner status:**
- Generally Available since August 7, 2025
- Free for public repositories
- Runner image: `ubuntu-24.04-arm`
- Hardware: Cobalt 100 processors (Arm Neoverse N2), 4 vCPUs

## Build Strategy Recommendation

### RECOMMENDED: Separate Dockerfile (Option A)

Create `Dockerfile.worker.arm64` that uses conda-forge packages instead of pre-compiled binaries.

**Advantages:**
- Clean separation of concerns
- No build-arg complexity
- Each Dockerfile optimized for its architecture
- Easier to maintain and debug
- Proven x86_64 approach unchanged

**Disadvantages:**
- Two files to maintain
- Some duplication (application code installation)

### NOT RECOMMENDED: Single Dockerfile with Build Args (Option B)

Use `ARG TARGETPLATFORM` to conditionally install packages.

**Why not:**
- Complex conditional logic (`if [ "$TARGETPLATFORM" = "linux/arm64" ]; then...`)
- Harder to test each path
- Mixed installation methods (binary vs conda)
- Fragile cross-compilation

### NOT RECOMMENDED: Unified Conda-Only Dockerfile (Option C)

Both architectures use conda-forge packages.

**Why not:**
- Larger x86_64 image (conda overhead vs pre-compiled binaries)
- Changes proven x86_64 approach unnecessarily
- Risk of introducing regressions in production x86_64

## GitHub Actions Strategy

### Current Workflow (publish-images.yml)

```yaml
# Current: Worker builds x86_64 only via QEMU-capable runner
- image: qm-nmr-calc-worker
  dockerfile: Dockerfile.worker
  platforms: linux/amd64  # arm64 not supported
```

### Recommended Workflow Update

```yaml
jobs:
  build-api:
    runs-on: ubuntu-latest
    # API image: Small Python image, QEMU acceptable
    # Continue using linux/amd64,linux/arm64 with QEMU

  build-worker-amd64:
    runs-on: ubuntu-latest
    steps:
      - uses: docker/build-push-action@v6
        with:
          file: Dockerfile.worker
          platforms: linux/amd64
          tags: ghcr.io/${{ github.repository_owner }}/qm-nmr-calc-worker:${{ github.ref_name }}-amd64
          # Push architecture-specific tag

  build-worker-arm64:
    runs-on: ubuntu-24.04-arm  # Native ARM64 runner
    steps:
      - uses: docker/build-push-action@v6
        with:
          file: Dockerfile.worker.arm64
          platforms: linux/arm64
          tags: ghcr.io/${{ github.repository_owner }}/qm-nmr-calc-worker:${{ github.ref_name }}-arm64
          # Push architecture-specific tag

  create-worker-manifest:
    needs: [build-worker-amd64, build-worker-arm64]
    runs-on: ubuntu-latest
    steps:
      # Combine into multi-arch manifest
      - run: |
          docker manifest create \
            ghcr.io/${{ github.repository_owner }}/qm-nmr-calc-worker:${{ github.ref_name }} \
            ghcr.io/${{ github.repository_owner }}/qm-nmr-calc-worker:${{ github.ref_name }}-amd64 \
            ghcr.io/${{ github.repository_owner }}/qm-nmr-calc-worker:${{ github.ref_name }}-arm64
          docker manifest push ghcr.io/${{ github.repository_owner }}/qm-nmr-calc-worker:${{ github.ref_name }}
```

**Why this approach:**
- Native ARM64 runner avoids 22x QEMU slowdown
- Each architecture builds with its optimized Dockerfile
- `docker manifest create` combines into single multi-arch tag
- Users pull single tag, Docker selects correct architecture

### Build Time Comparison

| Strategy | x86_64 Build | ARM64 Build | Total | Notes |
|----------|--------------|-------------|-------|-------|
| QEMU (current approach for API) | 3-5 min | 30-60 min | 35-65 min | ARM64 via emulation |
| Native runners (recommended) | 3-5 min | 3-5 min | 8-12 min | Parallel builds |
| Single QEMU build | N/A | N/A | 30-60 min | Blocked on ARM64 |

Reference: [actuated.com](https://actuated.com/blog/native-arm64-for-github-actions) reports 22x improvement (33.5 min to 1.5 min) for native vs QEMU builds.

## What NOT to Use

### DO NOT: Use Pre-compiled Binaries for ARM64

**Current x86_64 approach (Dockerfile.worker lines 30-41):**
```dockerfile
RUN wget https://github.com/grimme-lab/xtb/releases/download/v6.7.1/xtb-6.7.1-linux-x86_64.tar.xz
RUN wget https://github.com/crest-lab/crest/releases/download/v3.0.2/crest-gnu-12-ubuntu-latest.tar.xz
```

**Why this fails for ARM64:**
- grimme-lab/xtb releases do NOT include linux-aarch64 binaries
- crest-lab/crest releases do NOT include linux-aarch64 binaries
- NWChem base image is explicitly amd64 only (`ghcr.io/nwchemgit/nwchem-dev/amd64:latest`)

**What to use instead:** conda-forge packages provide all ARM64 builds.

### DO NOT: Use ghcr.io/nwchemgit/nwchem-dev/amd64:latest as Base

**Why:**
- Explicitly x86_64 only (it's in the image name: `/amd64:`)
- No ARM64 variant available from NWChem team
- Ties you to NWChem team's build schedule

**What to use instead:** `mambaorg/micromamba:2.5.0-debian12-slim` + conda-forge nwchem package.

### DO NOT: Use QEMU Emulation for Worker Builds

**Why:**
- 22x slower than native builds (measured)
- Can hang during conda operations (known bug with different CPU emulation)
- Wastes CI minutes
- Risk of subtle emulation bugs in floating-point calculations

**What to use instead:** `ubuntu-24.04-arm` native runners (free for public repos since August 2025).

### DO NOT: Use Alpine Base for Scientific Packages

**Why:**
- Alpine uses musl libc, not glibc
- Many conda-forge packages require glibc
- OpenMP threading can have issues on musl
- Known bug: QEMU + Alpine + micromamba can hang

**What to use instead:** Debian-based images (`mambaorg/micromamba:2.5.0-debian12-slim`).

### DO NOT: Mix Conda and pip Recklessly

**Why:**
- Package conflicts between pip-installed and conda-installed packages
- Harder to reproduce builds
- Can break conda's solver

**What to do instead:**
- Use conda/mamba for packages with native code (rdkit, scipy, matplotlib, pandas)
- Use pip (within the conda env) for pure-Python packages (huey, jinja2, pydantic)
- Define everything in a single `env.yaml` file with pip section

### DO NOT: Create Architecture-Specific Tags

**Bad:**
```
ghcr.io/steinbeck/qm-nmr-calc-worker:v2.5.0-arm64
ghcr.io/steinbeck/qm-nmr-calc-worker:v2.5.0-amd64
```

**Why:**
- Forces users to know their architecture
- Breaks `docker compose up` on different machines
- Defeats purpose of multi-arch images

**What to use instead:** Single multi-arch manifest tag (`:v2.5.0`, `:latest`). Docker automatically selects correct variant.

## Environment Configuration

### env-worker-arm64.yaml

```yaml
# Conda environment for ARM64 worker
name: base
channels:
  - conda-forge
dependencies:
  - python=3.11
  # Computational chemistry (native ARM64 builds)
  - nwchem=7.3.1
  - xtb=6.7.1
  - crest=3.0.2
  # Cheminformatics
  - rdkit=2025.09.5
  # Scientific computing (native ARM64 builds)
  - matplotlib>=3.10.0
  - pandas>=2.3.3
  - scipy>=1.17.0
  - statsmodels>=0.14.0
  # Pip for pure-Python packages
  - pip
  - pip:
    - huey>=2.6.0
    - jinja2>=3.1.0
    - orjson>=3.11.5
    - pydantic>=2.12.5
    - tqdm>=4.67.1
    - aiosmtplib>=5.0.0
    - email-validator>=2.3.0
    - python-multipart>=0.0.21
```

### Dockerfile.worker.arm64 (Skeleton)

```dockerfile
# ARM64 Worker container for qm-nmr-calc
# Uses conda-forge packages instead of pre-compiled binaries

FROM mambaorg/micromamba:2.5.0-debian12-slim

# OCI labels for GHCR repository linking
LABEL org.opencontainers.image.source="https://github.com/steinbeck/qm-nmr-calc"
LABEL org.opencontainers.image.description="QM NMR Calculator - ARM64 Worker with NWChem, CREST, xTB"
LABEL org.opencontainers.image.licenses="MIT"

# Copy environment definition
COPY --chown=$MAMBA_USER:$MAMBA_USER env-worker-arm64.yaml /tmp/env.yaml

# Install all dependencies via micromamba
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

# Enable environment activation in RUN commands
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# OpenMP for CREST/xTB (CRITICAL - prevents stack overflow on large molecules)
ENV OMP_STACKSIZE="2G"
ENV OMP_NUM_THREADS="4"
ENV GFORTRAN_UNBUFFERED_ALL="1"

# MPI settings for NWChem (conda-forge build uses OpenMPI)
ENV OMPI_ALLOW_RUN_AS_ROOT="1"
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1"

# NWChem process count (configurable)
ENV NWCHEM_NPROC="4"

# Copy application code
WORKDIR /app
COPY src/ /app/src/
COPY pyproject.toml README.md ./

# Install application via pip (within conda env)
RUN pip install --no-cache-dir .

# Create data directory
RUN mkdir -p /app/data

# Copy validation scripts
COPY scripts/validate-worker.sh /app/scripts/
COPY scripts/test-nwchem-container.sh /app/scripts/
RUN chmod +x /app/scripts/validate-worker.sh /app/scripts/test-nwchem-container.sh

# Default command: start Huey consumer
CMD ["huey_consumer", "qm_nmr_calc.queue.huey", "-w", "1", "-k", "process"]
```

### Key Differences from x86_64 Dockerfile

| Aspect | x86_64 (Current) | ARM64 (Proposed) |
|--------|------------------|------------------|
| Base image | ghcr.io/nwchemgit/nwchem-dev/amd64:latest | mambaorg/micromamba:2.5.0-debian12-slim |
| NWChem install | Pre-installed in base image | conda-forge package |
| xTB install | GitHub release binary | conda-forge package |
| CREST install | GitHub release binary | conda-forge package |
| RDKit install | pip (from PyPI) | conda-forge package |
| Python env | venv + pip | micromamba base env |
| Image size | ~2-3GB (estimate) | ~3-4GB (conda overhead) |

## Version Compatibility Matrix

| Component | x86_64 (Current) | ARM64 (Proposed) | Parity |
|-----------|------------------|------------------|--------|
| NWChem | 7.x (from base image) | 7.3.1 (conda-forge) | YES |
| xTB | 6.7.1 (binary) | 6.7.1 (conda-forge) | YES |
| CREST | 3.0.2 (binary) | 3.0.2 (conda-forge) | YES |
| RDKit | 2025.09.3 (pip) | 2025.09.5 (conda) | YES (newer) |
| Python | 3.11 | 3.11 | YES |

## Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| Conda package version mismatch with current worker | LOW | LOW | Versions align closely; test suite validates |
| Performance difference between architectures | MEDIUM | LOW | Benchmark during validation; document expected differences |
| Missing NWChem basis sets in conda package | LOW | HIGH | Verify NWCHEM_BASIS_LIBRARY path in conda package |
| OpenMPI configuration differences | MEDIUM | MEDIUM | Test MPI parallelism in container before release |
| Image size larger than x86_64 | HIGH | LOW | Conda overhead expected; functionality > size |
| CI build time increase | LOW | LOW | Native runners + parallel builds mitigate |

## Migration Path

### Phase 1: Create ARM64 Dockerfile
1. Create `Dockerfile.worker.arm64` with micromamba base
2. Create `env-worker-arm64.yaml` with conda-forge packages
3. Build locally on Apple Silicon for testing

### Phase 2: Validate Functionality
1. Run validation scripts (`validate-worker.sh`)
2. Execute test NWChem calculation
3. Compare results with x86_64 (< 0.01 ppm difference)
4. Test all calculation modes (single, CREST, RDKit-only)

### Phase 3: Update CI/CD
1. Add `ubuntu-24.04-arm` job for ARM64 worker build
2. Configure separate architecture builds
3. Add manifest creation job
4. Update matrix to remove API QEMU builds (optional optimization)

### Phase 4: Documentation & Release
1. Update README with ARM64 support notes
2. Document any performance differences
3. Tag v2.5.0 release
4. Verify multi-arch image in GHCR

## Sources

### conda-forge Package Verification (2026-02-04)
- [NWChem on conda-forge](https://anaconda.org/conda-forge/nwchem) - Version 7.3.1, linux-aarch64 available
- [xTB on conda-forge](https://anaconda.org/conda-forge/xtb) - Version 6.7.1, linux-aarch64 available
- [CREST on conda-forge](https://anaconda.org/conda-forge/crest) - Version 3.0.2, linux-aarch64 available
- [RDKit on conda-forge](https://anaconda.org/conda-forge/rdkit) - Version 2025.09.5, linux-aarch64 available

### Docker Base Image
- [mambaorg/micromamba on Docker Hub](https://hub.docker.com/r/mambaorg/micromamba) - Multi-arch image
- [micromamba-docker GitHub](https://github.com/mamba-org/micromamba-docker) - Source repository
- [micromamba-docker Documentation](https://micromamba-docker.readthedocs.io/en/latest/) - Official documentation
- [micromamba-docker Quick Start](https://micromamba-docker.readthedocs.io/en/latest/quick_start.html) - Dockerfile patterns

### Docker Multi-Architecture Builds
- [Docker Multi-Platform Builds](https://docs.docker.com/build/building/multi-platform/) - Official documentation
- [Docker Multi-Arch Build Blog](https://www.docker.com/blog/multi-arch-build-and-images-the-simple-way/) - Best practices
- [Multi-Arch Docker Images Guide](https://oneuptime.com/blog/post/2026-01-06-docker-multi-architecture-images/view) - 2026 guide
- [Multi-Arch Performance Analysis](https://skyworkz.nl/blog/multi-arch-docker-image-10x-faster/) - Manifest vs buildx performance

### GitHub Actions ARM64 Runners
- [GitHub ARM64 Runners GA Announcement](https://github.blog/changelog/2025-08-07-arm64-hosted-runners-for-public-repositories-are-now-generally-available/) - August 2025
- [GitHub Actions Multi-Platform Docker](https://docs.docker.com/build/ci/github-actions/multi-platform/) - Official guide
- [Native ARM64 Build Performance](https://actuated.com/blog/native-arm64-for-github-actions) - 22x speedup data
- [Multi-Arch Workflow Without QEMU](https://github.com/sredevopsorg/multi-arch-docker-github-workflow) - Reference implementation
