# Phase 39: CI/CD + GHCR Publishing - Research

**Researched:** 2026-02-03
**Domain:** GitHub Actions, Docker multi-platform builds, GHCR publishing
**Confidence:** HIGH

## Summary

This phase implements automated Docker image building and publishing to GitHub Container Registry (GHCR) using GitHub Actions. The standard approach uses Docker's official GitHub Actions (`docker/build-push-action`, `docker/metadata-action`, etc.) with multi-architecture support via buildx.

The critical finding concerns **arm64 support for scientific binaries**: NWChem provides official arm64 Docker images, but CREST and xTB do NOT provide pre-built arm64 binaries. For arm64, these must be installed via conda-forge. This requires modifying Dockerfile.worker to use conditional architecture-based installation, significantly complicating the arm64 build or requiring an amd64-only approach for the worker image.

**Primary recommendation:** Build both images for amd64 initially. For arm64 support of the API image (which has no architecture-specific binaries), include it in multi-arch builds. Defer full arm64 worker support or implement conditional Dockerfile logic using `TARGETARCH`.

## Standard Stack

The established tools for Docker CI/CD with GitHub Actions:

### Core Actions
| Action | Version | Purpose | Why Standard |
|--------|---------|---------|--------------|
| `docker/setup-qemu-action` | v3 | Cross-platform emulation | Required for multi-arch builds on single runner |
| `docker/setup-buildx-action` | v3 | Enable buildx builder | Required for multi-platform builds |
| `docker/login-action` | v3 | Registry authentication | Official, supports GHCR with GITHUB_TOKEN |
| `docker/build-push-action` | v6 | Build and push images | Most maintained, full feature set |
| `docker/metadata-action` | v5 | Generate tags/labels | Automatic semver handling |

### Supporting
| Tool | Purpose | When to Use |
|------|---------|-------------|
| `actions/checkout@v4` | Clone repository | Every workflow |
| GITHUB_TOKEN | Authentication | Built-in, no PAT needed |
| GitHub Actions cache | Speed up builds | Large images (use `type=gha`) |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| QEMU multi-arch | Native ARM64 runners | Faster but costs extra (GitHub charges for ARM runners) |
| Single workflow | Matrix per image | Single workflow simpler for 2 images |
| `type=gha` cache | Registry cache | Registry cache persists across branches but uses registry storage |

## Architecture Patterns

### Recommended Project Structure
```
.github/
└── workflows/
    └── publish-images.yml    # Single workflow for both images

# Image naming on GHCR:
ghcr.io/[owner]/qm-nmr-calc-api:v2.4.0
ghcr.io/[owner]/qm-nmr-calc-api:latest
ghcr.io/[owner]/qm-nmr-calc-worker:v2.4.0
ghcr.io/[owner]/qm-nmr-calc-worker:latest
```

### Pattern 1: Release Tag Trigger
**What:** Workflow triggers on semantic version tags (v*)
**When to use:** Release-based publishing (matches requirements)
**Example:**
```yaml
# Source: https://docs.github.com/actions/using-workflows/triggering-a-workflow
on:
  push:
    tags:
      - 'v*'  # Triggers on v1.0.0, v2.4.0, etc.
```

### Pattern 2: Multi-Image Build with Matrix
**What:** Single workflow builds multiple images using matrix strategy
**When to use:** Multiple Dockerfiles in same repo
**Example:**
```yaml
# Source: https://docs.docker.com/build/ci/github-actions/
strategy:
  matrix:
    include:
      - image: qm-nmr-calc-api
        dockerfile: Dockerfile.api
        platforms: linux/amd64,linux/arm64
      - image: qm-nmr-calc-worker
        dockerfile: Dockerfile.worker
        platforms: linux/amd64  # arm64 requires conda-forge for CREST/xTB
```

### Pattern 3: Automatic Tag Generation from Git Tag
**What:** metadata-action generates image tags from git tag
**When to use:** Always - ensures consistent versioning
**Example:**
```yaml
# Source: https://github.com/docker/metadata-action
- name: Docker meta
  id: meta
  uses: docker/metadata-action@v5
  with:
    images: ghcr.io/${{ github.repository_owner }}/${{ matrix.image }}
    tags: |
      type=semver,pattern={{version}}
      type=semver,pattern={{major}}.{{minor}}
      type=raw,value=latest,enable={{is_default_branch}}
```

### Anti-Patterns to Avoid
- **Hardcoding versions in workflow:** Use metadata-action instead - prevents tag mismatch
- **Using PAT when GITHUB_TOKEN works:** GITHUB_TOKEN is more secure, auto-rotates
- **Building arm64 via QEMU for large images:** Worker image is 2.1GB, QEMU build will timeout or take 30+ minutes
- **Separate workflows per image:** Harder to maintain, duplicate code

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Tag generation | Shell scripts parsing git tags | `docker/metadata-action` | Handles edge cases, OCI labels |
| Multi-arch manifest | Manual `docker manifest create` | `docker/build-push-action` with platforms | Automatic manifest list creation |
| GHCR auth | Manual docker login scripts | `docker/login-action` | Handles token securely |
| Build caching | Manual layer caching | `cache-from: type=gha` | Integrates with GH cache API |
| Conditional build | Complex shell conditionals | Matrix strategy with includes | Cleaner, parallelizable |

**Key insight:** Docker's official Actions handle all edge cases (manifest lists, OCI annotations, cache invalidation) that are easy to get wrong manually.

## Common Pitfalls

### Pitfall 1: QEMU Timeout on Large Images
**What goes wrong:** arm64 builds via QEMU take 10-30x longer than native, causing timeouts
**Why it happens:** QEMU emulates CPU instructions, extremely slow for compilation-heavy builds
**How to avoid:**
- For worker image (2.1GB with compiled binaries): Build amd64-only OR use native ARM runners
- For API image (733MB, no compilation): QEMU is acceptable
**Warning signs:** Build time > 30 minutes, workflow timeout at 6 hours

### Pitfall 2: Missing GITHUB_TOKEN Permissions
**What goes wrong:** Push to GHCR fails with 403 Forbidden
**Why it happens:** Default GITHUB_TOKEN lacks `packages: write` permission
**How to avoid:** Explicitly set permissions in workflow:
```yaml
permissions:
  contents: read
  packages: write
```
**Warning signs:** "denied: permission_denied" in logs

### Pitfall 3: Cache Size Limits
**What goes wrong:** Cache eviction causes full rebuilds
**Why it happens:** GitHub Actions cache limited to 10GB per repo, shared across all workflows
**How to avoid:**
- Use `mode=max` only when needed
- Use scope parameter for multiple images: `cache-to: type=gha,mode=max,scope=api`
**Warning signs:** "Cache not found" despite previous successful builds

### Pitfall 4: arm64 Binary Availability (CRITICAL)
**What goes wrong:** Worker container fails to build for arm64
**Why it happens:** CREST and xTB only provide x86_64 pre-built binaries
**How to avoid:**
- Option A: Build worker for amd64 only
- Option B: Modify Dockerfile.worker to use conda-forge for arm64:
  ```dockerfile
  ARG TARGETARCH
  RUN if [ "$TARGETARCH" = "amd64" ]; then \
        # Use pre-built binaries (current approach)
        wget xtb-6.7.1-linux-x86_64.tar.xz; \
      else \
        # Use conda-forge for arm64
        conda install -c conda-forge xtb crest; \
      fi
  ```
**Warning signs:** 404 errors downloading binaries during arm64 build

### Pitfall 5: NWChem Base Image Architecture
**What goes wrong:** Using wrong NWChem base image tag for architecture
**Why it happens:** NWChem uses `/amd64` and `/arm64` suffixes, not standard multi-arch manifest
**How to avoid:** Use ARG to select correct base:
```dockerfile
ARG TARGETARCH
FROM ghcr.io/nwchemgit/nwchem-dev/${TARGETARCH}:latest
```
**Warning signs:** "exec format error" when running container

### Pitfall 6: Cache API Version Mismatch
**What goes wrong:** Build fails with "legacy service is shutting down"
**Why it happens:** GitHub Cache API v1 deprecated April 2025
**How to avoid:** Use Docker Buildx >= v0.21.0 (docker/setup-buildx-action@v3 handles this)
**Warning signs:** Error mentioning "April 15, 2025" or "legacy service"

## Code Examples

Verified patterns from official sources:

### Complete Workflow for Multi-Image GHCR Publishing
```yaml
# Source: https://docs.docker.com/build/ci/github-actions/
# File: .github/workflows/publish-images.yml
name: Publish Docker Images

on:
  push:
    tags:
      - 'v*'

env:
  REGISTRY: ghcr.io

jobs:
  build-and-push:
    runs-on: ubuntu-latest
    permissions:
      contents: read
      packages: write

    strategy:
      fail-fast: false
      matrix:
        include:
          - image: qm-nmr-calc-api
            dockerfile: Dockerfile.api
            platforms: linux/amd64,linux/arm64
          - image: qm-nmr-calc-worker
            dockerfile: Dockerfile.worker
            platforms: linux/amd64  # arm64 requires Dockerfile changes

    steps:
      - name: Checkout
        uses: actions/checkout@v4

      - name: Set up QEMU
        uses: docker/setup-qemu-action@v3

      - name: Set up Docker Buildx
        uses: docker/setup-buildx-action@v3

      - name: Log in to GHCR
        uses: docker/login-action@v3
        with:
          registry: ${{ env.REGISTRY }}
          username: ${{ github.actor }}
          password: ${{ secrets.GITHUB_TOKEN }}

      - name: Extract metadata
        id: meta
        uses: docker/metadata-action@v5
        with:
          images: ${{ env.REGISTRY }}/${{ github.repository_owner }}/${{ matrix.image }}
          tags: |
            type=semver,pattern={{version}}
            type=semver,pattern={{major}}.{{minor}}
            type=raw,value=latest

      - name: Build and push
        uses: docker/build-push-action@v6
        with:
          context: .
          file: ${{ matrix.dockerfile }}
          platforms: ${{ matrix.platforms }}
          push: true
          tags: ${{ steps.meta.outputs.tags }}
          labels: ${{ steps.meta.outputs.labels }}
          cache-from: type=gha,scope=${{ matrix.image }}
          cache-to: type=gha,mode=max,scope=${{ matrix.image }}
```

### Dockerfile with Architecture-Conditional Installation
```dockerfile
# Source: https://docs.docker.com/build/building/multi-platform/
# Pattern for arm64-compatible worker (if arm64 support added later)

ARG TARGETARCH

# Select NWChem base image by architecture
FROM ghcr.io/nwchemgit/nwchem-dev/${TARGETARCH}:latest

ARG TARGETARCH

# Miniconda: use architecture-specific installer
RUN if [ "$TARGETARCH" = "amd64" ]; then \
      CONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-py311_24.11.1-0-Linux-x86_64.sh"; \
    else \
      CONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh"; \
    fi && \
    wget -q "$CONDA_URL" -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh

# xTB/CREST: pre-built for amd64, conda-forge for arm64
RUN if [ "$TARGETARCH" = "amd64" ]; then \
      # Pre-built binaries (faster, smaller)
      wget -q https://github.com/grimme-lab/xtb/releases/download/v6.7.1/xtb-6.7.1-linux-x86_64.tar.xz && \
      tar xf xtb-6.7.1-linux-x86_64.tar.xz -C /opt && \
      rm xtb-6.7.1-linux-x86_64.tar.xz; \
    else \
      # conda-forge for arm64
      /opt/conda/bin/conda install -y -c conda-forge xtb crest; \
    fi
```

### OCI Labels for Repository Linking
```dockerfile
# Source: https://docs.github.com/en/packages/working-with-a-github-packages-registry/working-with-the-container-registry
# Add to Dockerfiles to link images to repository

LABEL org.opencontainers.image.source="https://github.com/OWNER/qm-nmr-calc"
LABEL org.opencontainers.image.description="QM NMR Calculator - [api|worker] container"
LABEL org.opencontainers.image.licenses="MIT"
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `docker build` + `docker push` | `docker/build-push-action` | 2020+ | Single action, better caching |
| Manual tag parsing | `docker/metadata-action` | 2021+ | Automatic OCI annotations |
| Docker Hub default | GHCR for GitHub projects | 2020+ | Free for public repos, integrated |
| Single arch builds | Multi-arch with buildx | 2020+ | ARM support for Apple Silicon, Graviton |
| Cache API v1 | Cache API v2 | April 2025 | Must use Buildx >= v0.21.0 |

**Deprecated/outdated:**
- `docker/build-push-action@v5`: Update to v6
- `docker/metadata-action@v4`: Update to v5
- Manual `docker manifest` commands: Replaced by buildx multi-platform
- GitHub Cache API v1: Deprecated April 2025

## Open Questions

Things that couldn't be fully resolved:

1. **arm64 Worker Image Strategy**
   - What we know: CREST/xTB lack arm64 binaries; conda-forge provides them
   - What's unclear: Whether conda-forge versions match pre-built versions exactly
   - Recommendation: Start with amd64-only for worker; document arm64 as future enhancement

2. **NWChem arm64 Build Time**
   - What we know: NWChem provides arm64 base image
   - What's unclear: Build time for arm64 with QEMU emulation (2.1GB image)
   - Recommendation: If arm64 worker attempted, monitor build time closely; may need native runners

3. **GitHub Actions Cache Limits**
   - What we know: 10GB limit per repo, shared across all workflows
   - What's unclear: Whether 2.1GB worker image will cause cache thrashing
   - Recommendation: Use scope parameter to isolate caches; monitor cache eviction

## Sources

### Primary (HIGH confidence)
- [Docker Official Docs - Multi-platform GitHub Actions](https://docs.docker.com/build/ci/github-actions/multi-platform/) - Workflow patterns
- [Docker Official Docs - Cache Management](https://docs.docker.com/build/ci/github-actions/cache/) - Cache configuration
- [GitHub Docs - Container Registry](https://docs.github.com/en/packages/working-with-a-github-packages-registry/working-with-the-container-registry) - GHCR authentication
- [docker/metadata-action](https://github.com/docker/metadata-action) - Tag generation syntax
- [docker/build-push-action](https://github.com/docker/build-push-action) - Build action usage

### Secondary (MEDIUM confidence)
- [NWChem Dockerfiles](https://github.com/nwchemgit/nwchem-dockerfiles) - Confirmed arm64 image availability
- [xTB Releases](https://github.com/grimme-lab/xtb/releases) - Confirmed no arm64 binaries
- [CREST Releases](https://github.com/crest-lab/crest/releases) - Confirmed no arm64 binaries
- [conda-forge xtb](https://anaconda.org/conda-forge/xtb) - Confirmed linux-aarch64 availability
- [Miniconda Downloads](https://repo.anaconda.com/miniconda/) - Confirmed aarch64 installer exists

### Tertiary (LOW confidence)
- Various community blog posts on QEMU performance - Used for timeout estimates only

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - Official Docker documentation and Actions
- Architecture: HIGH - Verified with official docs
- Pitfalls: MEDIUM - arm64 binary availability verified; build times estimated from community reports

**Research date:** 2026-02-03
**Valid until:** 2026-05-03 (GitHub Actions ecosystem relatively stable)
