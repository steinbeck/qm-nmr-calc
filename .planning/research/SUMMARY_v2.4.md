# Research Summary: v2.4 Docker Deployment

**Project:** qm-nmr-calc v2.4 Docker Containerization
**Domain:** Containerized scientific computing deployment
**Researched:** 2026-02-02
**Confidence:** HIGH

## Executive Summary

Docker deployment for qm-nmr-calc requires balancing simplicity for academic researchers (who want `docker compose up` to work immediately) with production robustness (auto-HTTPS, health checks, persistent data). The research conclusively points to a three-container architecture: Caddy reverse proxy, FastAPI API server, and Huey worker with NWChem/CREST/xTB. The critical insight is to use official NWChem Docker images as the worker base rather than building from source, which saves 30-60 minutes of compilation and avoids complex MPI/BLAS dependency management.

The recommended approach uses Caddy for automatic HTTPS (zero configuration Let's Encrypt), pre-built images published to GHCR, and named Docker volumes for data persistence. This architecture supports three user personas: individual chemists running locally, research group admins with shared lab servers, and cloud deployers on VPS instances. The key differentiation from manual installation is the bundled HTTPS and production defaults that "just work" with a domain name in `.env`.

Primary risks center on NWChem containerization: shared memory requirements (`shm_size: 512m`), MPI process binding in containers (`--bind-to none`), and scratch disk exhaustion for large DFT calculations. Secondary risks include SQLite queue locking under multi-container load (mitigated by WAL mode and busy timeout; Redis recommended for scale). All critical pitfalls have documented prevention strategies that must be implemented in Phase 1.

## Key Findings

### Recommended Stack

The deployment stack builds on proven, official images to minimize complexity and maximize reliability.

**Core technologies:**
- **Docker Compose v2:** Declarative multi-container orchestration with `compose.yaml` format
- **Caddy 2.10:** Auto-HTTPS reverse proxy with zero-config Let's Encrypt and HTTP/3
- **ghcr.io/nwchemgit/nwchem-dev:amd64:** Official NWChem image with optimized OpenMPI builds
- **ghcr.io/astral-sh/uv:python3.12-bookworm-slim:** Fast Python packaging for API container
- **GHCR (ghcr.io):** Container registry co-located with source, no Docker Hub rate limits
- **GitHub Actions + docker/build-push-action v6:** CI/CD with cache mounts and multi-arch support

**Why Caddy over Traefik:** Caddy wins for single-application deployments. The Caddyfile is human-readable (~10 lines for full HTTPS), auto-HTTPS works without ACME configuration, and it's sufficient for single-node deployment. Traefik would be better for Kubernetes or multiple dynamic services.

**Why GHCR over Docker Hub:** Native GitHub integration, no rate limits for public repos, no separate account required.

### Expected Features

**Must have (table stakes):**
- Single `docker compose up` command for deployment
- Pre-built images on GHCR (users should not need to build)
- All dependencies bundled (NWChem, CREST, xTB in worker image)
- Volume persistence for job data and queue state
- Health checks and restart policies
- Environment configuration via `.env` file
- Quick start documentation (5-minute deployment)
- `.env.example` with all variables documented
- Non-root container user (security baseline)

**Should have (differentiators):**
- Auto-HTTPS via Caddy with Let's Encrypt
- Graceful shutdown handling (no orphaned calculations)
- Container resource limits (memory caps for runaway calculations)
- Multi-architecture images (amd64 + arm64)
- Backup-friendly named volumes
- Worker parallelism configuration via `NWCHEM_NPROC`

**Defer to v2.5+:**
- Kubernetes manifests (overkill for target users)
- Redis-backed queue (keep SQLite for single-node MVP)
- Built-in monitoring stack (Grafana/Prometheus)
- Multi-user authentication
- Rate limiting in application
- Horizontal autoscaling

**Anti-features (deliberately exclude):**
- Docker Swarm mode
- Custom orchestration
- Service mesh
- Blue-green deployment
- Secrets manager integration (`.env` file is sufficient)

### Architecture Approach

The architecture maps naturally to three containers with a shared data volume. The API container is lightweight (Python + FastAPI), the worker container is heavy but single-instance (NWChem + CREST + xTB + Huey), and Caddy handles external traffic.

**Container split:**

```
Internet --> Caddy:80/443 --> API:8000 <-- shared volume --> Worker
                                  |
                            [qm-data volume]
                            - /data/huey.db
                            - /data/jobs/
```

**Major components:**
1. **Caddy service:** TLS termination, reverse proxy, security headers, HTTP/2/3
2. **API service:** FastAPI server, web UI, job submission, status queries
3. **Worker service:** Huey consumer, NWChem/CREST/xTB calculations, result processing

**Shared state handling:**
- SQLite Huey queue with WAL mode and 5-second busy timeout
- Single named volume (`qm-data`) for both queue DB and job directories
- Both containers run as same UID for permission consistency

**Critical container configuration:**
```yaml
worker:
  shm_size: '512m'      # Required for NWChem MPI
  restart: unless-stopped
  deploy:
    resources:
      limits:
        memory: 8G      # DFT calculations need RAM
```

### Critical Pitfalls

The top pitfalls that can cause complete deployment failure or data loss:

1. **NWChem shared memory exhaustion** -- Docker default `/dev/shm` is 64MB; NWChem MPI needs 256-512MB minimum. **Fix:** Set `shm_size: '512m'` in compose. Without this, calculations segfault immediately.

2. **MPI process binding failures** -- OpenMPI 5.x has container CPU detection issues causing hangs or "not enough slots" errors. **Fix:** Always use `--bind-to none` in mpirun commands (already in codebase, must preserve in container).

3. **CREST/xTB stack overflow** -- Default 8MB stack is insufficient for 50+ atom molecules. **Fix:** Set `OMP_STACKSIZE=2G` and `ulimit -s unlimited` in container entrypoint.

4. **Scratch directory disk exhaustion** -- Single DFT calculation can generate 5GB scratch files. **Fix:** Mount dedicated scratch volume sized for `max_jobs * 5GB`, implement cleanup.

5. **SQLite locking under load** -- Multi-container access without WAL mode causes "database is locked" errors. **Fix:** Configure WAL mode with `busy_timeout: 5000`. For production scale, migrate to Redis.

6. **Hardcoded paths break in container** -- Relative paths assume working directory context. **Fix:** Use environment variables for all data paths, set `WORKDIR /app` in Dockerfile.

7. **Container image size explosion** -- Naive Dockerfiles can reach 10GB. **Fix:** Multi-stage builds, clean apt lists in same layer, use `.dockerignore`.

## Implications for Roadmap

Based on research, Docker deployment should proceed in six phases, ordered by dependency and complexity.

### Phase 1: Worker Container (Foundation)

**Rationale:** The worker is the most complex container and blocks everything else. NWChem, CREST, and xTB containerization must be proven before any other work.

**Delivers:**
- Dockerfile.worker based on `ghcr.io/nwchemgit/nwchem-dev:amd64`
- CREST and xTB binaries installed (pre-built from GitHub releases)
- Python + Huey consumer configuration
- All environment variables set (OMP_STACKSIZE, NWCHEM paths, MPI settings)
- Validation that calculations run correctly in container

**Addresses features:** Dependencies bundled, worker parallelism config
**Avoids pitfalls:** #1 shm_size, #2 MPI binding, #3 stack overflow, #8 binary not found, #10 env conflicts

### Phase 2: API Container

**Rationale:** Simpler than worker, but must be validated separately before orchestration.

**Delivers:**
- Dockerfile.api based on `ghcr.io/astral-sh/uv:python3.12-bookworm-slim`
- Multi-stage build for minimal image size
- Health check endpoint validation
- Non-root user configuration

**Addresses features:** Non-root user, minimal images
**Avoids pitfalls:** #9 image size, #11 root user

### Phase 3: Docker Compose Integration

**Rationale:** With both containers working, wire them together with shared state.

**Delivers:**
- `compose.yaml` with api, worker services
- Shared `qm-data` volume configuration
- SQLite WAL mode verification
- Environment variable handling via `.env`
- `.env.example` with all options documented
- Health checks and restart policies

**Addresses features:** Single command deployment, volume persistence, restart policies, health checks, .env configuration
**Avoids pitfalls:** #1 SQLite locking, #4 scratch exhaustion, #6 orphaned jobs, #7 hardcoded paths

### Phase 4: Caddy + HTTPS

**Rationale:** Production readiness requires HTTPS. Caddy is simple enough to add after compose is working.

**Delivers:**
- Caddyfile with auto-HTTPS and reverse proxy
- Caddy service in compose with `caddy_data` volume
- Domain configuration via `DOMAIN` env var
- HTTP-to-HTTPS redirect
- Security headers (X-Content-Type-Options, X-Frame-Options)

**Addresses features:** Auto-HTTPS, production-ready defaults
**Avoids pitfalls:** None specific (Caddy is well-documented)

### Phase 5: CI/CD + GHCR Publishing

**Rationale:** Pre-built images are table stakes for user adoption. Must be last because it publishes the completed containers.

**Delivers:**
- GitHub Actions workflow for build and push
- GHCR authentication with `GITHUB_TOKEN`
- Image tagging strategy (semver + sha + branch)
- Cache-from/cache-to for build performance
- Multi-architecture builds (amd64 + arm64)

**Addresses features:** Pre-built images on GHCR, multi-arch support
**Avoids pitfalls:** None specific

### Phase 6: Documentation + Polish

**Rationale:** Documentation ensures adoption. Must be written after implementation is stable.

**Delivers:**
- Quick start section in README (5-minute deployment)
- Full deployment guide for production
- Troubleshooting section with common issues
- Update and backup procedures
- Resource sizing recommendations

**Addresses features:** Quick start docs, deployment guide, troubleshooting
**Avoids pitfalls:** Documents all pitfall prevention strategies

### Phase Ordering Rationale

1. **Worker first:** All other phases depend on calculations working in containers. This is the riskiest and most complex component.
2. **API second:** Validates the simpler container before combining them.
3. **Compose third:** Integration testing with shared state happens here.
4. **Caddy fourth:** HTTPS is important but not blocking for local development testing.
5. **CI/CD fifth:** Can only publish images after containers are stable.
6. **Docs sixth:** Documentation reflects final implementation.

### Research Flags

**Phases likely needing deeper research during planning:**
- **Phase 1 (Worker Container):** NWChem MPI tuning in containers, CREST binary compatibility, optimal memory settings for different molecule sizes
- **Phase 5 (CI/CD):** Multi-arch build configuration, GHCR authentication edge cases

**Phases with standard patterns (skip research-phase):**
- **Phase 2 (API Container):** Well-documented FastAPI Docker patterns
- **Phase 3 (Docker Compose):** Standard orchestration, existing codebase patterns
- **Phase 4 (Caddy):** Extensive documentation, simple configuration
- **Phase 6 (Documentation):** Standard technical writing

## Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| Stack | HIGH | Official Docker docs, NWChem containers wiki, uv Docker guide |
| Features | HIGH | Docker best practices, scientific computing patterns established |
| Architecture | HIGH | Standard microservices pattern, matches existing process model |
| Pitfalls | MEDIUM-HIGH | Real issues from GitHub issues, some may have been fixed in recent versions |

**Overall confidence:** HIGH

### Gaps to Address

1. **OpenMPI version compatibility:** Research indicates OpenMPI 5.x has container issues; may need to test 4.x fallback. Validate during Phase 1.

2. **Multi-arch builds for scientific binaries:** NWChem and CREST may not have arm64 builds. Verify availability or mark arm64 as unsupported.

3. **SQLite lock contention under real load:** WAL mode should suffice for single-node, but stress test during Phase 3 integration testing.

4. **CREST/xTB installation method:** Research recommends pre-built binaries over conda-forge. Validate binary compatibility with container base image.

## Sources

### Primary (HIGH confidence)
- [NWChem Containers Documentation](https://nwchemgit.github.io/Containers.html) -- Official container guidance, shm_size requirements
- [uv Docker Integration Guide](https://docs.astral.sh/uv/guides/integration/docker/) -- Multi-stage builds, cache mounts
- [Caddy Official Documentation](https://caddyserver.com/docs/) -- Reverse proxy, auto-HTTPS
- [GitHub Container Registry Docs](https://docs.github.com/en/packages/working-with-a-github-packages-registry/working-with-the-container-registry) -- GHCR publishing
- [docker/build-push-action](https://github.com/docker/build-push-action) -- CI/CD patterns

### Secondary (MEDIUM confidence)
- [SQLite WAL Mode](https://sqlite.org/wal.html) -- Concurrency patterns
- [NWChem Dockerfiles Repository](https://github.com/nwchemgit/nwchem-dockerfiles) -- Build patterns
- [xTB Setup Documentation](https://xtb-docs.readthedocs.io/en/latest/setup.html) -- OpenMP configuration
- [Huey Documentation](https://huey.readthedocs.io/) -- Storage backend options

### Issue Trackers (Real-World Problems)
- [Huey SqliteHuey locking #445](https://github.com/coleifer/huey/issues/445) -- Validates locking concerns
- [OpenMPI Docker issues #12431](https://github.com/open-mpi/ompi/issues/12431) -- Container binding problems
- [xTB Stack overflow #191](https://github.com/grimme-lab/xtb/issues/191) -- Stack size requirements

---
*Research completed: 2026-02-02*
*Ready for roadmap: yes*
