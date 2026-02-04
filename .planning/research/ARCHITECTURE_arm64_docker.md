# Architecture: ARM64 Docker Worker Integration

**Project:** qm-nmr-calc ARM64 Worker Support
**Researched:** 2026-02-04
**Confidence:** HIGH (verified with official sources)

## Executive Summary

ARM64 support for the qm-nmr-calc worker container is **feasible** and follows established Docker multi-arch patterns. The key insight is that all upstream dependencies (NWChem, xTB, CREST) now support ARM64 via either official images or conda-forge packages, eliminating the previous blockers.

**Architecture change:** Replace pre-built x86_64 binaries with conda-forge packages for xTB and CREST, enabling cross-platform support while using NWChem's official ARM64 image as the base.

---

## Current Architecture (x86_64 Only)

```
+------------------------------------------------------------------+
|                    Current Worker Image                           |
+------------------------------------------------------------------+
|  BASE: ghcr.io/nwchemgit/nwchem-dev/amd64:latest                 |
|        (x86_64 ONLY)                                              |
+------------------------------------------------------------------+
|  LAYER 1: Miniconda3-py311_24.11.1-0-Linux-x86_64.sh             |
|           (x86_64 binary download)                                |
+------------------------------------------------------------------+
|  LAYER 2: xTB v6.7.1 - xtb-6.7.1-linux-x86_64.tar.xz             |
|           (x86_64 binary download from GitHub releases)           |
+------------------------------------------------------------------+
|  LAYER 3: CREST v3.0.2 - crest-gnu-12-ubuntu-latest.tar.xz       |
|           (x86_64 binary download from GitHub releases)           |
+------------------------------------------------------------------+
|  LAYER 4: Python venv + application dependencies                  |
+------------------------------------------------------------------+
```

**Problem:** Three hardcoded x86_64 dependencies:
1. NWChem base image uses `/amd64` tag
2. xTB binary downloaded from GitHub releases (x86_64 only)
3. CREST binary downloaded from GitHub releases (x86_64 only)

---

## Proposed ARM64 Architecture

### Option A: Separate Dockerfiles (Recommended)

Create architecture-specific Dockerfiles that share common layers where possible.

```
+------------------------------------------------------------------+
|                    ARM64 Worker Image                             |
+------------------------------------------------------------------+
|  BASE: ghcr.io/nwchemgit/nwchem-dev/arm64:latest                 |
|        (ARM64 native)                                             |
+------------------------------------------------------------------+
|  LAYER 1: Miniconda3-latest-Linux-aarch64.sh                     |
|           OR Miniforge3-Linux-aarch64.sh                          |
+------------------------------------------------------------------+
|  LAYER 2: xTB 6.7.1 from conda-forge (linux-aarch64 package)     |
|           conda install -c conda-forge xtb                        |
+------------------------------------------------------------------+
|  LAYER 3: CREST 3.x from conda-forge (linux-aarch64 package)     |
|           conda install -c conda-forge crest                      |
+------------------------------------------------------------------+
|  LAYER 4: Python venv + application dependencies                  |
+------------------------------------------------------------------+
```

**Implementation:**
- `Dockerfile.worker.amd64` - Current approach with pre-built binaries (faster builds)
- `Dockerfile.worker.arm64` - Uses conda-forge for xTB/CREST

### Option B: Single Dockerfile with ARG-based Selection

Use Docker build arguments and conditional logic:

```dockerfile
# Build-time architecture detection
ARG TARGETARCH
ARG TARGETPLATFORM

# Select NWChem base based on architecture
FROM ghcr.io/nwchemgit/nwchem-dev/${TARGETARCH}:latest

# Architecture-specific Miniconda installer
ARG MINICONDA_ARCH=${TARGETARCH}
RUN if [ "${MINICONDA_ARCH}" = "amd64" ]; then \
        MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-py311_24.11.1-0-Linux-x86_64.sh"; \
    else \
        MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh"; \
    fi && \
    wget -q ${MINICONDA_URL} -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh

# Use conda-forge for architecture-portable xTB/CREST installation
RUN /opt/conda/bin/conda install -y -c conda-forge xtb crest && \
    /opt/conda/bin/conda clean -afy
```

**Trade-off:** Single Dockerfile but conda-forge adds ~500MB vs pre-built binaries.

### Recommendation: Option A (Separate Dockerfiles)

**Rationale:**
1. **Build speed:** Pre-built binaries for x86_64 are faster to download than conda package resolution
2. **Image size:** x86_64 can stay lean (~3GB), ARM64 slightly larger due to conda (~3.5GB)
3. **Maintenance:** Changes to one architecture don't risk breaking the other
4. **CI clarity:** Clear matrix separation in GitHub Actions

---

## Component Integration Points

### 1. NWChem Base Image

| Architecture | Image Path | Status |
|--------------|-----------|--------|
| x86_64 (amd64) | `ghcr.io/nwchemgit/nwchem-dev/amd64:latest` | Currently used |
| ARM64 (aarch64) | `ghcr.io/nwchemgit/nwchem-dev/arm64:latest` | Available, verified |
| ppc64le | `ghcr.io/nwchemgit/nwchem-dev/ppc64le:latest` | Available (not needed) |

**Source:** [NWChem Containers Documentation](https://nwchemgit.github.io/Containers.html)

The ARM64 image uses the same OpenMPI configuration as x86_64, ensuring computational parity.

### 2. Miniconda/Miniforge

| Architecture | Installer | Source |
|--------------|-----------|--------|
| x86_64 | `Miniconda3-py311_24.11.1-0-Linux-x86_64.sh` | repo.anaconda.com |
| ARM64 | `Miniconda3-latest-Linux-aarch64.sh` | repo.anaconda.com |
| ARM64 (alt) | `Miniforge3-Linux-aarch64.sh` | github.com/conda-forge/miniforge |

**Recommendation:** Use Miniforge for ARM64 builds - better tested on aarch64, conda-forge default channel.

### 3. xTB (Extended Tight-Binding)

| Source | x86_64 | ARM64 | Notes |
|--------|--------|-------|-------|
| GitHub Releases | Yes (v6.7.1) | **NO** | Only x86_64 binaries published |
| conda-forge | Yes (v6.7.1) | **YES** (v6.7.1) | linux-aarch64 package available |

**Source:** [xTB conda-forge](https://anaconda.org/conda-forge/xtb)

### 4. CREST (Conformer-Rotamer Ensemble Sampling)

| Source | x86_64 | ARM64 | Notes |
|--------|--------|-------|-------|
| GitHub Releases | Yes (v3.0.2) | **NO** | Only x86_64 binaries published |
| conda-forge | Yes | **YES** | linux-aarch64 package available |

**Source:** [CREST conda-forge](https://anaconda.org/conda-forge/crest)

---

## Multi-Arch Build Strategy

### Docker Buildx Manifest Lists

The goal is a single image tag that works on both architectures:

```
ghcr.io/steinbeck/qm-nmr-calc-worker:v2.5.0
    |
    +-- linux/amd64 manifest -> built from Dockerfile.worker.amd64
    |
    +-- linux/arm64 manifest -> built from Dockerfile.worker.arm64
```

When users pull `ghcr.io/steinbeck/qm-nmr-calc-worker:v2.5.0`, Docker automatically selects the correct architecture.

### Build Approach: Separate Build + Manifest Merge

Since xTB/CREST installation differs significantly between architectures, use the "build separately, merge manifests" pattern:

```yaml
# Phase 1: Build architecture-specific images
jobs:
  build-amd64:
    runs-on: ubuntu-latest
    steps:
      - uses: docker/build-push-action@v6
        with:
          file: Dockerfile.worker.amd64
          platforms: linux/amd64
          tags: ghcr.io/steinbeck/qm-nmr-calc-worker:${{ github.sha }}-amd64
          push: true

  build-arm64:
    runs-on: ubuntu-24.04-arm  # Native ARM runner (faster)
    # OR: ubuntu-latest with QEMU (slower but works)
    steps:
      - uses: docker/build-push-action@v6
        with:
          file: Dockerfile.worker.arm64
          platforms: linux/arm64
          tags: ghcr.io/steinbeck/qm-nmr-calc-worker:${{ github.sha }}-arm64
          push: true

  # Phase 2: Create multi-arch manifest
  create-manifest:
    needs: [build-amd64, build-arm64]
    steps:
      - uses: docker/login-action@v3
        with:
          registry: ghcr.io
      - run: |
          docker buildx imagetools create -t ghcr.io/steinbeck/qm-nmr-calc-worker:latest \
            ghcr.io/steinbeck/qm-nmr-calc-worker:${{ github.sha }}-amd64 \
            ghcr.io/steinbeck/qm-nmr-calc-worker:${{ github.sha }}-arm64
```

**Source:** [Docker Multi-Platform Docs](https://docs.docker.com/build/building/multi-platform/)

### Alternative: Single-Job QEMU Build

Simpler but slower (~3x build time for ARM64 under emulation):

```yaml
jobs:
  build:
    runs-on: ubuntu-latest
    steps:
      - uses: docker/setup-qemu-action@v3
      - uses: docker/setup-buildx-action@v3
      - uses: docker/build-push-action@v6
        with:
          file: Dockerfile.worker  # Must handle both architectures
          platforms: linux/amd64,linux/arm64
          push: true
          tags: ghcr.io/steinbeck/qm-nmr-calc-worker:latest
```

**Recommendation:** Use separate builds with native ARM runner (if available) or accept QEMU overhead for simpler workflow.

---

## CI/CD Changes Required

### Current Workflow: `.github/workflows/publish-images.yml`

```yaml
# Current (x86_64 only)
matrix:
  include:
    - image: qm-nmr-calc-worker
      dockerfile: Dockerfile.worker
      platforms: linux/amd64  # ARM64 blocked by CREST/xTB binaries
```

### Proposed Workflow Update

```yaml
# Updated for multi-arch
matrix:
  include:
    - image: qm-nmr-calc-api
      dockerfile: Dockerfile.api
      platforms: linux/amd64,linux/arm64  # Already works

    # Worker needs architecture-specific Dockerfiles
    - image: qm-nmr-calc-worker
      dockerfile: Dockerfile.worker.amd64
      platforms: linux/amd64

    - image: qm-nmr-calc-worker
      dockerfile: Dockerfile.worker.arm64
      platforms: linux/arm64

# Add manifest merge job after builds complete
merge-worker-manifests:
  needs: build-and-push
  # ... create multi-arch manifest
```

### Build Time Estimates

| Build | Method | Estimated Time |
|-------|--------|----------------|
| Worker amd64 | Native x86_64 runner | ~5 minutes |
| Worker arm64 | Native ARM runner | ~8 minutes (conda package resolution) |
| Worker arm64 | QEMU emulation | ~15-20 minutes |
| API (both) | Native + QEMU | ~3 minutes |

---

## Dockerfile Structure

### Dockerfile.worker.arm64 (New)

```dockerfile
# Worker container for qm-nmr-calc (ARM64/aarch64)
# Uses conda-forge for xTB and CREST (no pre-built ARM64 binaries available)

FROM ghcr.io/nwchemgit/nwchem-dev/arm64:latest

LABEL org.opencontainers.image.source="https://github.com/steinbeck/qm-nmr-calc"
LABEL org.opencontainers.image.description="QM NMR Calculator - Worker (ARM64)"
LABEL org.opencontainers.image.licenses="MIT"

USER root

# Install required tools and X11 libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    ca-certificates \
    libxrender1 \
    libxext6 \
    && rm -rf /var/lib/apt/lists/*

# Install Miniforge (better ARM64 support than Miniconda)
ENV CONDA_DIR=/opt/conda
RUN wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-aarch64.sh -O /tmp/miniforge.sh \
    && bash /tmp/miniforge.sh -b -p $CONDA_DIR \
    && rm /tmp/miniforge.sh
ENV PATH="$CONDA_DIR/bin:$PATH"

# Install xTB and CREST from conda-forge (ARM64 packages available)
RUN conda install -y -c conda-forge xtb crest \
    && conda clean -afy

# [Rest of Dockerfile identical to amd64 version]
# Create venv, install Python dependencies, set environment variables...
```

### Dockerfile.worker.amd64 (Renamed from Dockerfile.worker)

Keep the current `Dockerfile.worker` approach with pre-built binaries for better build speed on x86_64.

---

## Docker Compose Compatibility

No changes required to `docker-compose.yml`. The worker service definition works with both architectures:

```yaml
worker:
  build:
    context: .
    dockerfile: Dockerfile.worker  # Or use image: ghcr.io/.../worker:latest
  # All other settings work identically on both architectures
```

When using pre-built images:
```yaml
worker:
  image: ghcr.io/steinbeck/qm-nmr-calc-worker:latest
  # Docker automatically pulls correct architecture
```

---

## Testing Strategy

### Local Testing on ARM64 Mac (Apple Silicon)

```bash
# Build ARM64 image locally
docker build -f Dockerfile.worker.arm64 -t qm-nmr-calc-worker:arm64 .

# Run validation script
docker run --rm qm-nmr-calc-worker:arm64 /app/scripts/validate-worker.sh
```

### CI Testing Matrix

```yaml
test:
  strategy:
    matrix:
      include:
        - arch: amd64
          runner: ubuntu-latest
        - arch: arm64
          runner: ubuntu-24.04-arm  # Or use QEMU
  runs-on: ${{ matrix.runner }}
  steps:
    - run: docker run --rm ghcr.io/.../worker:${{ matrix.arch }} /app/scripts/validate-worker.sh
```

---

## Suggested Implementation Order

### Phase 1: Prepare (Low Risk)

1. Create `Dockerfile.worker.arm64` with conda-forge approach
2. Test locally on ARM64 Mac or via QEMU
3. Verify NWChem, xTB, CREST all function correctly

### Phase 2: CI Integration (Medium Risk)

4. Rename `Dockerfile.worker` to `Dockerfile.worker.amd64`
5. Update `.github/workflows/publish-images.yml`:
   - Add ARM64 build job
   - Add manifest merge job
6. Test workflow on a branch before merging

### Phase 3: Validation (Low Risk)

7. Verify multi-arch manifest works:
   ```bash
   docker buildx imagetools inspect ghcr.io/steinbeck/qm-nmr-calc-worker:latest
   ```
8. Test pull on both x86_64 and ARM64 machines
9. Run full calculation test on ARM64

---

## Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| conda-forge xTB/CREST versions lag behind | LOW | MEDIUM | Pin versions, monitor feedstocks |
| ARM64 performance differs from x86_64 | MEDIUM | LOW | Benchmark; document differences |
| QEMU build timeout in CI | MEDIUM | LOW | Use native ARM runner or increase timeout |
| NWChem ARM64 bugs | LOW | HIGH | Test thoroughly before release |

---

## Sources

### Authoritative (HIGH Confidence)

- [NWChem Container Documentation](https://nwchemgit.github.io/Containers.html) - Confirms ARM64 image availability
- [Docker Multi-Platform Documentation](https://docs.docker.com/build/building/multi-platform/) - Buildx patterns
- [xTB conda-forge](https://anaconda.org/conda-forge/xtb) - ARM64 package availability
- [CREST conda-forge](https://anaconda.org/conda-forge/crest) - ARM64 package availability
- [Miniforge Releases](https://github.com/conda-forge/miniforge) - ARM64 installer

### Supporting (MEDIUM Confidence)

- [GitHub Actions Multi-Arch Guide](https://docs.docker.com/build/ci/github-actions/multi-platform/) - CI workflow patterns
- [xTB GitHub Releases](https://github.com/grimme-lab/xtb/releases) - Confirms no ARM64 binaries
- [CREST GitHub Releases](https://github.com/crest-lab/crest/releases) - Confirms no ARM64 binaries
- [Miniconda System Requirements](https://www.anaconda.com/docs/getting-started/miniconda/system-requirements) - aarch64 support

---

## Summary

ARM64 worker support is achievable with:

1. **NWChem:** Use official `ghcr.io/nwchemgit/nwchem-dev/arm64` base image
2. **xTB/CREST:** Switch from pre-built binaries to conda-forge packages
3. **Build Strategy:** Separate Dockerfiles per architecture, merged into single multi-arch manifest
4. **CI Changes:** Add ARM64 build job + manifest merge step

The main trade-off is ARM64 image size (~500MB larger due to conda) and slightly longer build times for conda package resolution. This is acceptable given the benefit of native ARM64 support for Apple Silicon Macs and AWS Graviton instances.
