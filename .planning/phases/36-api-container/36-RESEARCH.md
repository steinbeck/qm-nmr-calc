# Phase 36: API Container - Research

**Researched:** 2026-02-02
**Domain:** Docker containerization for FastAPI web application
**Confidence:** HIGH

## Summary

This phase focuses on containerizing the FastAPI API server (not the computational worker) for the qm-nmr-calc application. The API container handles web requests, serves static files (CSS, templates), and queues jobs to Huey for the worker container to process. Unlike the worker container (Phase 35) which requires NWChem/CREST/xTB, the API container only needs Python with the application dependencies.

The research confirms that `ghcr.io/astral-sh/uv:python3.11-bookworm-slim` is the optimal base image, providing uv pre-installed with a minimal Debian base. A multi-stage build pattern separates dependency installation from the final runtime, producing images under 200MB. Health check endpoints already exist in the codebase (`/health` for liveness, `/health/ready` for readiness) and should be used with Docker HEALTHCHECK or Kubernetes probes.

**Primary recommendation:** Use the official uv Python image, implement multi-stage build with non-root user, run uvicorn with `--workers 1` (Kubernetes scaling) or `--workers 4` (single server), and configure HEALTHCHECK using the existing `/health` endpoint.

## Standard Stack

### Core

| Component | Version | Purpose | Why Standard |
|-----------|---------|---------|--------------|
| `ghcr.io/astral-sh/uv:python3.11-bookworm-slim` | latest | Base image with uv + Python 3.11 | Official Astral image, matches pyproject.toml Python version |
| uvicorn | 0.40.0+ | ASGI server | Built-in worker management, no Gunicorn needed |
| python:3.11-slim-bookworm | 3.11 | Runtime-only stage | 51MB download, minimal attack surface |

### Supporting

| Component | Purpose | When to Use |
|-----------|---------|-------------|
| curl | Health check probe | HEALTHCHECK command in Dockerfile |
| tini | Init process | If PID 1 signal handling needed (optional) |

### Alternatives Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| uv image | python:3.11-slim + pip | uv 10x faster installs, better Docker caching |
| python:slim-bookworm | python:alpine | Alpine musl libc causes issues with some Python packages |
| uvicorn --workers | gunicorn + uvicorn | Gunicorn adds complexity, uvicorn now handles worker management natively |
| COPY static files | nginx sidecar | Overkill for small static files, FastAPI StaticFiles sufficient |

**Installation:**
```bash
# Base image already has uv, just need:
uv sync --locked --no-dev
```

## Architecture Patterns

### Recommended Dockerfile Structure

```dockerfile
# Stage 1: Build dependencies
FROM ghcr.io/astral-sh/uv:python3.11-bookworm-slim AS builder

WORKDIR /app

# Environment for uv
ENV UV_COMPILE_BYTECODE=1
ENV UV_LINK_MODE=copy
ENV UV_NO_DEV=1

# Install dependencies first (cache layer)
RUN --mount=type=cache,target=/root/.cache/uv \
    --mount=type=bind,source=uv.lock,target=uv.lock \
    --mount=type=bind,source=pyproject.toml,target=pyproject.toml \
    uv sync --locked --no-install-project

# Copy source and install project
COPY . /app
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --locked

# Stage 2: Runtime (minimal)
FROM python:3.11-slim-bookworm

# Install curl for health checks
RUN apt-get update && apt-get install -y --no-install-recommends curl \
    && rm -rf /var/lib/apt/lists/*

# Create non-root user
RUN groupadd --system --gid 999 appgroup \
    && useradd --system --uid 999 --gid appgroup --create-home appuser

WORKDIR /app

# Copy virtual environment from builder
COPY --from=builder --chown=appuser:appgroup /app/.venv /app/.venv

# Copy application source (includes templates, static files)
COPY --from=builder --chown=appuser:appgroup /app/src /app/src

# Add venv to PATH
ENV PATH="/app/.venv/bin:$PATH"

# Switch to non-root user
USER appuser

# Health check using existing endpoint
HEALTHCHECK --interval=30s --timeout=5s --start-period=10s --retries=3 \
    CMD curl --fail http://localhost:8000/health || exit 1

EXPOSE 8000

# Run uvicorn
CMD ["uvicorn", "qm_nmr_calc.api.app:app", "--host", "0.0.0.0", "--port", "8000"]
```

### Project Structure in Container

```
/app/
├── .venv/              # Python virtual environment (from builder)
│   └── bin/
│       └── uvicorn     # ASGI server
├── src/
│   └── qm_nmr_calc/
│       ├── api/
│       │   ├── app.py          # FastAPI application
│       │   ├── routers/        # API routes
│       │   ├── templates/      # Jinja2 HTML templates
│       │   └── static/         # CSS files
│       └── ...                 # Other modules
└── data/               # Volume mount point
    ├── jobs/           # Job directories (shared with worker)
    └── huey.db         # SQLite queue database (shared with worker)
```

### Pattern 1: Environment Variable Configuration

**What:** Configure the application via environment variables
**When to use:** Always in Docker deployments
**Example:**
```dockerfile
# Uvicorn settings
ENV UVICORN_HOST=0.0.0.0
ENV UVICORN_PORT=8000
ENV UVICORN_WORKERS=1

# Application settings (used by existing code)
ENV DATA_DIR=/app/data
```

### Pattern 2: Volume Mounts for Shared State

**What:** API and worker containers share data via volume mounts
**When to use:** Docker Compose deployments
**Example:**
```yaml
# docker-compose.yaml
services:
  api:
    volumes:
      - qm_data:/app/data
  worker:
    volumes:
      - qm_data:/app/data

volumes:
  qm_data:
```

### Anti-Patterns to Avoid

- **Running as root:** Security risk, always use non-root user
- **Using `--reload` in production:** Disables workers, only for development
- **Including dev dependencies:** Use `UV_NO_DEV=1` or `--no-dev` flag
- **Hardcoding ports:** Use environment variables for flexibility
- **Missing health checks:** Container orchestrators need them to manage containers properly
- **Using shell form for CMD:** Use exec form `["uvicorn", ...]` for proper signal handling

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Process management | Supervisor/systemd | uvicorn --workers | Uvicorn now handles worker restarts natively |
| Static file serving | Custom middleware | FastAPI StaticFiles | Already implemented in codebase, production-ready |
| Health endpoints | New endpoints | Existing /health, /health/ready | Already implemented with proper checks |
| Python package management | pip freeze/requirements.txt | uv sync --locked | Deterministic builds, faster installs |
| Signal handling | Custom handlers | Exec form CMD | Proper SIGTERM propagation built-in |

**Key insight:** The codebase already has health endpoints with data directory and Huey queue checks. Leverage existing code rather than creating Docker-specific health checks.

## Common Pitfalls

### Pitfall 1: Static Files Not Found

**What goes wrong:** 404 errors for CSS/templates after containerization
**Why it happens:** Paths are relative to module location, WORKDIR mismatch
**How to avoid:**
- The codebase uses `Path(__file__).resolve().parent` which is correct
- Ensure COPY includes `src/qm_nmr_calc/api/static/` and `src/qm_nmr_calc/api/templates/`
- Verify with `curl http://localhost:8000/static/css/base.css`
**Warning signs:** Web UI loads but looks unstyled

### Pitfall 2: Huey Database Inaccessible

**What goes wrong:** "database is locked" or jobs don't appear in queue
**Why it happens:** SQLite database path mismatch between API and worker
**How to avoid:**
- Mount same volume to `/app/data` in both containers
- The codebase uses `./data/huey.db` relative path
- Set `WORKDIR /app` so relative path resolves to `/app/data/huey.db`
**Warning signs:** Jobs submitted but worker never picks them up

### Pitfall 3: Permission Denied on Data Directory

**What goes wrong:** 403 or IOError when creating job directories
**Why it happens:** Non-root user can't write to volume mount
**How to avoid:**
- In docker-compose, use named volumes (Docker manages permissions)
- For bind mounts, ensure host directory has appropriate permissions
- Alternative: Create data directory in Dockerfile with correct ownership
**Warning signs:** Job submission fails immediately

### Pitfall 4: Container Doesn't Stop Gracefully

**What goes wrong:** 10-second delay on `docker stop`, SIGKILL instead of SIGTERM
**Why it happens:** Shell form CMD doesn't forward signals properly
**How to avoid:**
- Use exec form: `CMD ["uvicorn", "..."]` not `CMD uvicorn ...`
- Uvicorn's `--timeout-graceful-shutdown` default is reasonable (30s)
**Warning signs:** Slow container restarts, abrupt connection drops

### Pitfall 5: Workers > 1 with SQLite

**What goes wrong:** "database is locked" errors under load
**Why it happens:** Multiple uvicorn workers + SQLite concurrency limits
**How to avoid:**
- For single container: `--workers 1` (scale with replicas instead)
- For multiple workers: Consider PostgreSQL for Huey (future enhancement)
**Warning signs:** Sporadic 500 errors under concurrent requests

### Pitfall 6: Health Check Fails During Startup

**What goes wrong:** Container marked unhealthy before app starts
**Why it happens:** Health check runs before uvicorn is ready
**How to avoid:**
- Use `--start-period=10s` in HEALTHCHECK (grace period)
- Alternative: Kubernetes startupProbe separate from livenessProbe
**Warning signs:** Container restart loop, "unhealthy" status

## Code Examples

### Existing Health Endpoints (Source: src/qm_nmr_calc/api/routers/health.py)

```python
# Liveness probe - minimal, always returns 200 if process is running
@router.get("/health")
async def liveness():
    return {"status": "alive"}

# Readiness probe - checks dependencies
@router.get("/health/ready")
async def readiness():
    checks = {}
    # Checks data directory writable
    # Checks Huey queue importable
    # Reports CREST availability (informational)
    return {"status": "ready", "checks": checks, ...}
```

These endpoints are already production-ready. Use `/health` for Docker HEALTHCHECK (fast, no side effects) and `/health/ready` for Kubernetes readiness probes.

### Docker HEALTHCHECK Command

```dockerfile
# Use curl with fail flag - returns non-zero on HTTP errors
HEALTHCHECK --interval=30s --timeout=5s --start-period=10s --retries=3 \
    CMD curl --fail http://localhost:8000/health || exit 1
```

### Uvicorn Production Command

```dockerfile
# Single worker (recommended for Kubernetes, scale with replicas)
CMD ["uvicorn", "qm_nmr_calc.api.app:app", "--host", "0.0.0.0", "--port", "8000"]

# Multiple workers (for single VM deployment)
CMD ["uvicorn", "qm_nmr_calc.api.app:app", "--host", "0.0.0.0", "--port", "8000", "--workers", "4"]

# With proxy headers (behind nginx/traefik)
CMD ["uvicorn", "qm_nmr_calc.api.app:app", "--host", "0.0.0.0", "--port", "8000", "--proxy-headers"]
```

### Non-Root User Setup

```dockerfile
# Create system user with specific UID (reproducible)
RUN groupadd --system --gid 999 appgroup \
    && useradd --system --uid 999 --gid appgroup --create-home appuser

# Set ownership when copying files
COPY --from=builder --chown=appuser:appgroup /app/.venv /app/.venv
COPY --from=builder --chown=appuser:appgroup /app/src /app/src

# Switch before CMD
USER appuser
```

### Environment Variables for Production

```dockerfile
# uv build settings (builder stage)
ENV UV_COMPILE_BYTECODE=1     # Faster startup
ENV UV_LINK_MODE=copy         # Required for cache mounts
ENV UV_NO_DEV=1               # Exclude dev dependencies

# Runtime settings
ENV PATH="/app/.venv/bin:$PATH"
ENV PYTHONUNBUFFERED=1        # Real-time logging
ENV PYTHONDONTWRITEBYTECODE=1 # Prevent .pyc in container
```

### Kubernetes Probe Configuration

```yaml
# For Kubernetes deployments
livenessProbe:
  httpGet:
    path: /health
    port: 8000
  initialDelaySeconds: 5
  periodSeconds: 10
  timeoutSeconds: 3
  failureThreshold: 3

readinessProbe:
  httpGet:
    path: /health/ready
    port: 8000
  initialDelaySeconds: 10
  periodSeconds: 10
  timeoutSeconds: 5
  failureThreshold: 3
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Gunicorn + Uvicorn workers | uvicorn --workers | 2024 | Simpler configuration, fewer dependencies |
| pip install requirements.txt | uv sync --locked | 2024 | 10x faster builds, deterministic |
| Single-stage Dockerfile | Multi-stage build | Standard practice | 50-70% smaller images |
| tiangolo/uvicorn-gunicorn-fastapi | Custom Dockerfile | 2024 | Image deprecated, better control |
| Root user in container | Non-root user | Standard practice | Required for security compliance |

**Deprecated/outdated:**
- `tiangolo/uvicorn-gunicorn-fastapi` image: Deprecated, Uvicorn now handles workers natively
- `gunicorn -k uvicorn.workers.UvicornWorker`: Unnecessary complexity for most deployments
- Alpine-based Python images: musl libc causes issues with scientific Python packages

## Open Questions

1. **Optimal Worker Count for Single VM**
   - What we know: Uvicorn can use multiple workers with `--workers N`
   - What's unclear: Whether SQLite Huey handles concurrent access from multiple workers
   - Recommendation: Start with `--workers 1`, scale via container replicas; if needed, increase workers and monitor for SQLite lock errors

2. **Graceful Shutdown Timeout**
   - What we know: Uvicorn default is 30 seconds for graceful shutdown
   - What's unclear: Whether long-running requests exist in API (they shouldn't, jobs are async)
   - Recommendation: Use default timeout, adjust if needed

3. **TLS Termination**
   - What we know: Application expects to run behind a reverse proxy (Traefik/nginx)
   - What's unclear: Whether to add `--proxy-headers` by default
   - Recommendation: Add `--proxy-headers` flag, harmless if not behind proxy

## Sources

### Primary (HIGH confidence)
- [uv Docker Guide](https://docs.astral.sh/uv/guides/integration/docker/) - Official multi-stage patterns, environment variables
- [uv-docker-example repository](https://github.com/astral-sh/uv-docker-example) - Official example Dockerfile with non-root user
- [Uvicorn Settings](https://uvicorn.dev/settings/) - All configuration options
- [FastAPI Docker deployment](https://fastapi.tiangolo.com/deployment/docker/) - Official FastAPI Docker guide

### Secondary (MEDIUM confidence)
- [Better Stack FastAPI Docker Guide](https://betterstack.com/community/guides/scaling-python/fastapi-docker-best-practices/) - Best practices compilation
- [Docker USER instruction](https://www.docker.com/blog/understanding-the-docker-user-instruction/) - Non-root user patterns
- [Python Docker Security (Snyk)](https://snyk.io/advisor/docker/python) - Image vulnerability analysis

### Codebase (HIGH confidence - verified)
- `src/qm_nmr_calc/api/app.py` - FastAPI application entry point, static files mounting
- `src/qm_nmr_calc/api/routers/health.py` - Existing health check endpoints
- `src/qm_nmr_calc/storage.py` - Data directory path configuration (`DATA_DIR = Path("./data/jobs")`)
- `src/qm_nmr_calc/queue.py` - Huey database path (`'./data/huey.db'`)

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - Official Astral images and uvicorn documentation verified
- Architecture: HIGH - Based on existing codebase patterns and official Docker guides
- Pitfalls: HIGH - Verified from uvicorn docs, codebase paths, and SQLite behavior

**Research date:** 2026-02-02
**Valid until:** 2026-03-02 (30 days - stable domain)
