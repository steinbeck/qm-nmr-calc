# Technology Stack: Docker Deployment

**Project:** qm-nmr-calc v2.4 Docker Containerization
**Researched:** 2026-02-02
**Focus:** Containerization and deployment stack only (existing app stack validated)

## Executive Summary

This stack research focuses exclusively on containerizing and deploying the existing qm-nmr-calc application. The core application stack (FastAPI, Huey, NWChem, RDKit, CREST/xTB) is already validated and working. This research identifies the optimal tools for containerization, orchestration, and production deployment with auto-HTTPS.

**Key Decision:** Use official NWChem Docker images from GHCR rather than building from source. Building NWChem from source requires extensive compilation time (30+ minutes) and complex dependency management. The official images provide tested, optimized builds with MPI support.

---

## Recommended Stack

### Container Runtime

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| Docker Engine | 27.x | Container runtime | Industry standard, required for target deployment |
| Docker Compose | v2.x | Multi-container orchestration | Declarative, version-controlled deployment |
| BuildKit | Built-in | Image building | Cache mounts, multi-stage builds |

**Rationale:** Docker Compose v2 (integrated with Docker CLI) provides declarative multi-container deployment. The `compose.yaml` format is the current standard (supersedes `docker-compose.yml`).

### Base Images

| Image | Tag | Purpose | Size | Source |
|-------|-----|---------|------|--------|
| `ghcr.io/nwchemgit/nwchem-dev` | `amd64` | NWChem calculations | ~2GB | [NWChem GHCR](https://github.com/nwchemgit/nwchem-dockerfiles) |
| `ghcr.io/astral-sh/uv` | `python3.12-bookworm-slim` | Python/uv base | ~200MB | [Astral GHCR](https://github.com/astral-sh/uv) |
| `caddy` | `2.10-alpine` | Reverse proxy | ~50MB | [Docker Hub Official](https://hub.docker.com/_/caddy) |

### NWChem Container Strategy

**CRITICAL DECISION: Use official NWChem image as computation base, NOT multi-stage from scratch.**

**Why:**
1. NWChem compilation takes 30-60 minutes and requires ~10GB build space
2. Requires OpenMPI, BLAS, LAPACK, and Fortran compiler chains
3. Official images are tested with specific MPI configurations
4. The `ghcr.io/nwchemgit/nwchem-dev/amd64` image includes optimized OpenMPI builds

**Container Architecture:**

```
+-------------------------------------------------------------+
|                    Docker Compose Stack                      |
+-------------------------------------------------------------+
|  +---------+  +--------------+  +-------------------------+ |
|  |  Caddy  |  |   FastAPI    |  |     Huey Worker         | |
|  | :80/443 |--|    :8000     |  |  (NWChem + CREST/xTB)   | |
|  |         |  |              |  |                         | |
|  | alpine  |  | python-slim  |  |  nwchem-dev + python    | |
|  +---------+  +--------------+  +-------------------------+ |
|       |              |                     |                 |
|       +--------------+---------------------+                 |
|                      |                                       |
|              +-------+-------+                              |
|              |    Volumes    |                              |
|              |  jobs_data    |                              |
|              |  huey_db      |                              |
|              |  caddy_data   |                              |
|              +---------------+                              |
+-------------------------------------------------------------+
```

### Reverse Proxy

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| Caddy | 2.10.x | Reverse proxy, auto-HTTPS | Zero-config TLS, automatic Let's Encrypt |

**Why Caddy over Traefik:**
- **Simplicity:** Caddyfile is human-readable, ~10 lines for full HTTPS setup
- **Auto-HTTPS:** Works out of box, no ACME configuration needed
- **Single binary:** No external dependencies
- **HTTP/3:** Enabled by default
- **Sufficient for single-node:** qm-nmr-calc is single-worker, doesn't need Traefik's dynamic service discovery

**When Traefik would be better:**
- Kubernetes deployment
- Multiple dynamic services
- Complex routing rules with middlewares

**Recommendation:** Caddy is the right choice for this single-application deployment.

### CI/CD & Registry

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| GitHub Actions | N/A | CI/CD pipeline | Native to GitHub, free for public repos |
| GHCR (ghcr.io) | N/A | Container registry | Co-located with source, free, no Docker Hub limits |
| docker/build-push-action | v6 | Build & publish | Official, multi-arch support, attestations |
| docker/setup-buildx-action | v3 | BuildKit setup | Required for cache mounts |

**GHCR Publishing Pattern:**
```yaml
permissions:
  contents: read
  packages: write
  id-token: write  # For attestations

steps:
  - uses: docker/setup-buildx-action@v3
  - uses: docker/login-action@v3
    with:
      registry: ghcr.io
      username: ${{ github.actor }}
      password: ${{ secrets.GITHUB_TOKEN }}
  - uses: docker/build-push-action@v6
    with:
      push: true
      tags: ghcr.io/${{ github.repository }}:latest
      cache-from: type=gha
      cache-to: type=gha,mode=max
```

---

## Container Build Strategy

### Image 1: API Server (`qm-nmr-calc-api`)

**Base:** `ghcr.io/astral-sh/uv:python3.12-bookworm-slim`

**Build pattern:** Multi-stage with uv

```dockerfile
# Stage 1: Build dependencies
FROM ghcr.io/astral-sh/uv:python3.12-bookworm-slim AS builder
WORKDIR /app
COPY pyproject.toml uv.lock ./
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --locked --no-install-project --no-dev

# Stage 2: Runtime
FROM python:3.12-slim-bookworm
COPY --from=builder /app/.venv /app/.venv
ENV PATH="/app/.venv/bin:$PATH"
COPY src/ /app/src/
EXPOSE 8000
CMD ["uvicorn", "qm_nmr_calc.api.main:app", "--host", "0.0.0.0"]
```

**Key optimizations:**
- `UV_COMPILE_BYTECODE=1` for faster startup
- `--no-dev` excludes test dependencies
- Cache mounts persist across builds
- Separate layers for dependencies vs code

### Image 2: Huey Worker (`qm-nmr-calc-worker`)

**Base:** `ghcr.io/nwchemgit/nwchem-dev/amd64` (extended with Python)

**Build pattern:** Layer Python on NWChem base

```dockerfile
FROM ghcr.io/nwchemgit/nwchem-dev/amd64 AS nwchem-base

# Install Python and uv on top of NWChem
FROM nwchem-base
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3.12 python3.12-venv python3-pip curl && \
    rm -rf /var/lib/apt/lists/*

# Install uv
COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv

# Install CREST and xTB from conda-forge (if not using conda base)
# Alternative: Copy pre-built binaries
COPY --from=crest-builder /usr/local/bin/crest /usr/local/bin/
COPY --from=crest-builder /usr/local/bin/xtb /usr/local/bin/

WORKDIR /app
COPY pyproject.toml uv.lock ./
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --locked --no-dev

COPY src/ /app/src/
ENV PATH="/app/.venv/bin:$PATH"
CMD ["python", "-m", "huey.bin.huey_consumer", "qm_nmr_calc.queue.huey"]
```

**CREST/xTB options:**
1. **Pre-built binaries:** Download from GitHub releases, copy into image
2. **Conda-forge:** Use `micromamba` to install from conda-forge (adds ~500MB)
3. **Skip in base image:** Mark as optional, document manual addition

**Recommendation:** Pre-built binaries from GitHub releases. CREST and xTB publish static binaries that can be directly copied. This avoids conda complexity.

### What NOT to Containerize

| Component | Reason | Alternative |
|-----------|--------|-------------|
| Job data filesystem | Needs persistence, may be large | Docker volume (named or bind) |
| Huey SQLite database | Persistence required | Docker volume |
| Caddy certificates | Must persist across restarts | Named volume (`caddy_data`) |
| Configuration files | May need runtime editing | Bind mounts or ConfigMaps |

---

## Production Configuration

### Health Checks

**API Server:**
```yaml
healthcheck:
  test: ["CMD", "curl", "-f", "http://localhost:8000/api/v1/health/ready"]
  interval: 30s
  timeout: 10s
  start_period: 10s
  retries: 3
```

**Worker:**
```yaml
healthcheck:
  test: ["CMD", "python", "-c", "import qm_nmr_calc.queue; print('ok')"]
  interval: 60s
  timeout: 30s
  start_period: 30s
  retries: 3
```

### Restart Policies

| Service | Policy | Rationale |
|---------|--------|-----------|
| caddy | `unless-stopped` | Critical path, should always run |
| api | `unless-stopped` | Web server, should always run |
| worker | `on-failure:5` | May fail on bad jobs, limit restart loops |

### Resource Limits

```yaml
services:
  worker:
    deploy:
      resources:
        limits:
          memory: 8G  # NWChem DFT needs RAM
        reservations:
          memory: 4G
    # NWChem shared memory requirement
    shm_size: 256m
```

**CRITICAL:** The `shm_size: 256m` is required for NWChem MPI communication. Without it, parallel calculations will fail.

### Volumes

```yaml
volumes:
  jobs_data:      # /app/data/jobs - calculation results
  huey_db:        # /app/data/huey.db - task queue state
  caddy_data:     # /data - Caddy certificates and state
  caddy_config:   # /config - Caddy configuration
```

---

## Caddy Configuration

### Minimal Caddyfile for Auto-HTTPS

```caddyfile
{
    email admin@example.com
}

example.com {
    reverse_proxy api:8000
}
```

That's it. Caddy automatically:
- Obtains Let's Encrypt certificate
- Redirects HTTP to HTTPS
- Renews certificates before expiry
- Enables HTTP/3

### Development (localhost)

```caddyfile
:80 {
    reverse_proxy api:8000
}
```

For localhost, Caddy uses self-signed certificates from its internal CA.

---

## GHCR Publishing Workflow

### Complete Workflow Example

```yaml
name: Build and Push

on:
  push:
    branches: [main]
    tags: ['v*']

env:
  REGISTRY: ghcr.io
  API_IMAGE: ghcr.io/${{ github.repository }}-api
  WORKER_IMAGE: ghcr.io/${{ github.repository }}-worker

jobs:
  build:
    runs-on: ubuntu-latest
    permissions:
      contents: read
      packages: write
      id-token: write

    steps:
      - uses: actions/checkout@v4

      - uses: docker/setup-buildx-action@v3

      - uses: docker/login-action@v3
        with:
          registry: ${{ env.REGISTRY }}
          username: ${{ github.actor }}
          password: ${{ secrets.GITHUB_TOKEN }}

      - uses: docker/metadata-action@v5
        id: meta-api
        with:
          images: ${{ env.API_IMAGE }}
          tags: |
            type=ref,event=branch
            type=semver,pattern={{version}}
            type=sha,prefix=

      - uses: docker/build-push-action@v6
        with:
          context: .
          file: ./docker/Dockerfile.api
          push: true
          tags: ${{ steps.meta-api.outputs.tags }}
          cache-from: type=gha
          cache-to: type=gha,mode=max
```

### Image Tagging Strategy

| Trigger | Tags Generated |
|---------|----------------|
| Push to main | `main`, `sha-abc1234` |
| Tag v1.2.3 | `1.2.3`, `1.2`, `1`, `latest` |

---

## Version Summary

| Component | Version | Confidence | Source |
|-----------|---------|------------|--------|
| Docker Compose | v2.x | HIGH | [Docker Docs](https://docs.docker.com/compose/) |
| Caddy | 2.10.x | HIGH | [Caddy endoflife.date](https://endoflife.date/caddy) |
| NWChem image | dev/amd64 | HIGH | [NWChem Containers](https://nwchemgit.github.io/Containers.html) |
| uv | 0.9.x | HIGH | [Astral uv Docs](https://docs.astral.sh/uv/guides/integration/docker/) |
| Python | 3.12 | HIGH | Matches pyproject.toml requirement |
| build-push-action | v6 | HIGH | [GitHub docker/build-push-action](https://github.com/docker/build-push-action) |

---

## Alternatives Considered

| Decision | Chosen | Alternative | Why Not Alternative |
|----------|--------|-------------|---------------------|
| Reverse proxy | Caddy | Traefik | Overkill for single-app deployment |
| Registry | GHCR | Docker Hub | Rate limits, requires separate account |
| NWChem source | Official image | Build from source | 30-60 min build time, complex deps |
| CREST/xTB | Pre-built binaries | Conda-forge | Adds ~500MB, conda complexity |
| Base OS | Debian (bookworm) | Alpine | RDKit and scientific packages better supported |
| Task queue | Keep Huey | Redis/Celery | Already works, no benefit from switching |

---

## Sources

### Authoritative (HIGH confidence)
- [NWChem Container Documentation](https://nwchemgit.github.io/Containers.html)
- [uv Docker Integration Guide](https://docs.astral.sh/uv/guides/integration/docker/)
- [Caddy Official Docker Image](https://hub.docker.com/_/caddy)
- [GitHub Container Registry Docs](https://docs.github.com/en/packages/working-with-a-github-packages-registry/working-with-the-container-registry)
- [docker/build-push-action](https://github.com/docker/build-push-action)

### Supporting (MEDIUM confidence)
- [Caddy vs Traefik Comparison](https://www.programonaut.com/reverse-proxies-compared-traefik-vs-caddy-vs-nginx-docker/)
- [Docker Compose Health Checks Guide](https://last9.io/blog/docker-compose-health-checks/)
- [Multi-Stage Docker Builds for Python with uv](https://hynek.me/articles/docker-uv/)
- [NWChem Dockerfiles Repository](https://github.com/nwchemgit/nwchem-dockerfiles)
- [Publishing Docker Images to GHCR](https://docs.github.com/actions/guides/publishing-docker-images)

### CREST/xTB Sources
- [CREST conda-forge](https://anaconda.org/conda-forge/crest)
- [xTB conda-forge](https://anaconda.org/conda-forge/xtb)
- [CREST GitHub Releases](https://github.com/crest-lab/crest/releases)
