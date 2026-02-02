# Architecture Research: Docker Compose Deployment

**Project:** qm-nmr-calc v2.4 Docker Deployment
**Researched:** 2026-02-02
**Confidence:** HIGH

## Executive Summary

The qm-nmr-calc architecture naturally maps to a three-container Docker Compose stack: API server, background worker, and reverse proxy. The key architectural decisions center on handling shared state (SQLite queue database and job filesystem) and properly containerizing the NWChem/CREST/xTB computational dependencies.

**Key findings:**
- SQLite with WAL mode supports concurrent access from API and worker containers via shared volume
- Single shared volume handles both huey.db and job directories cleanly
- NWChem requires MPI configuration and shared memory allocation in container
- Caddy provides automatic HTTPS with minimal configuration
- Multi-stage builds keep image sizes reasonable despite scientific computing dependencies

## Recommended Container Architecture

### Container Split

```
                    Internet
                        |
                        v
                +---------------+
                |     Caddy     |  <- Auto-HTTPS, reverse proxy
                |   (caddy)     |
                +---------------+
                        |
            +-----------+-----------+
            |                       |
            v                       v
    +---------------+       +---------------+
    |   FastAPI     |       |    Huey       |
    |    (api)      |       |   (worker)    |
    +---------------+       +---------------+
            |                       |
            +-----------+-----------+
                        |
                        v
                +---------------+
                |  Shared Data  |  <- Named volume: qm-data
                |   Volume      |
                +---------------+
                        |
            +-----------+-----------+
            |                       |
            v                       v
        ./data/                 ./data/
        huey.db                 jobs/
```

### Service Definitions

| Service | Container | Purpose | Base Image | Build Priority |
|---------|-----------|---------|------------|----------------|
| `api` | qm-nmr-calc-api | FastAPI HTTP server | python:3.11-slim | 2 (after base) |
| `worker` | qm-nmr-calc-worker | Huey consumer + NWChem/CREST | python:3.11-slim + conda | 1 (complex deps) |
| `caddy` | caddy:2-alpine | Reverse proxy, auto-HTTPS | caddy:2-alpine | N/A (official image) |

### Why This Split

**API container (lightweight):**
- Only needs Python + FastAPI + RDKit (via conda-forge)
- No NWChem/CREST/xTB required (only worker runs calculations)
- Fast startup, minimal memory footprint
- Scales independently (multiple replicas possible)

**Worker container (heavy):**
- Contains NWChem + OpenMPI + CREST + xTB
- Large base image (~2GB) but only one instance needed
- CPU-intensive, benefits from all available cores
- Long-running processes (minutes to hours per job)

**Caddy container (external):**
- Official image, no build required
- Automatic Let's Encrypt certificates
- Zero-config HTTPS for production
- Development can bypass with direct port mapping

## Shared State Handling

### SQLite Huey Queue Database

The Huey SQLite database (`./data/huey.db`) must be shared between API (task producer) and worker (task consumer).

**SQLite concurrency in containers:**
- SQLite supports one writer + multiple readers simultaneously with WAL mode
- Containers are just processes; SQLite handles cross-process coordination via file locks
- Shared volume must support proper file locking (local filesystem works; NFS may not)

**Configuration for reliable sharing:**

```python
# queue.py - enable WAL mode for better concurrency
huey = SqliteHuey(
    'qm-nmr-calc',
    filename='./data/huey.db',
    fsync=True,
    pragmas={
        'journal_mode': 'wal',      # Enable WAL mode
        'synchronous': 'normal',     # Good performance + safety
        'busy_timeout': 5000,        # Wait 5s for locks
    }
)
```

**Volume requirements:**
- Must be a Docker named volume or bind mount to local filesystem
- MUST NOT use NFS/CIFS without proper lock support
- Both containers mount same volume at same path (`/app/data`)

### Job Directory Structure

Jobs are stored in `./data/jobs/{job_id}/`:

```
./data/
    huey.db              # Shared queue database
    huey.db-wal          # WAL file (auto-created)
    huey.db-shm          # Shared memory file (auto-created)
    jobs/
        abc123def456/
            status.json       # Job metadata
            output/           # Results
            scratch/          # NWChem working files
            logs/             # NWChem output logs
```

**Access patterns:**
- API creates job directories, writes initial status.json, reads completed results
- Worker reads job input, writes scratch/output/logs, updates status.json

**File ownership:**
- Both containers run as same UID (non-root recommended)
- Volume permissions allow read/write from both services

## Volume Mounts

### Named Volume Configuration

```yaml
# docker-compose.yml
volumes:
  qm-data:
    name: qm-nmr-calc-data

services:
  api:
    volumes:
      - qm-data:/app/data
    working_dir: /app

  worker:
    volumes:
      - qm-data:/app/data
    working_dir: /app
```

### Bind Mount Alternative (Development)

```yaml
# docker-compose.dev.yml
services:
  api:
    volumes:
      - ./data:/app/data
      - ./src:/app/src:ro    # Code changes without rebuild

  worker:
    volumes:
      - ./data:/app/data
      - ./src:/app/src:ro
```

### Volume Considerations

| Concern | Named Volume | Bind Mount |
|---------|--------------|------------|
| Data persistence | Survives `docker compose down` | Lives in host filesystem |
| Portability | Managed by Docker | Tied to host paths |
| Performance | Native Docker storage | Host filesystem speed |
| Inspection | `docker volume inspect` | Direct filesystem access |
| Backup | `docker cp` or volume plugins | Standard filesystem tools |

**Recommendation:** Named volume for production, bind mount for development.

## Network Configuration

### Internal Network

```yaml
networks:
  qm-internal:
    driver: bridge

services:
  api:
    networks:
      - qm-internal
    # No ports exposed directly

  worker:
    networks:
      - qm-internal
    # No ports needed (pulls from Huey queue)

  caddy:
    networks:
      - qm-internal
    ports:
      - "80:80"
      - "443:443"
```

### Service Discovery

- Caddy reaches API via `http://api:8000` (Docker DNS)
- Worker has no network requirements (communicates via shared volume/SQLite)
- Only Caddy exposes ports to host

### Caddyfile Configuration

```caddyfile
{
    email {$ACME_EMAIL}
}

{$DOMAIN:localhost} {
    reverse_proxy api:8000

    # Health check endpoint (no auth)
    handle /health* {
        reverse_proxy api:8000
    }

    # Compress responses
    encode gzip

    # Security headers
    header {
        X-Content-Type-Options nosniff
        X-Frame-Options DENY
        Referrer-Policy strict-origin-when-cross-origin
    }
}
```

## Resource Limits

### CPU and Memory Guidelines

```yaml
services:
  api:
    deploy:
      resources:
        limits:
          cpus: '2'
          memory: 1G
        reservations:
          cpus: '0.5'
          memory: 256M

  worker:
    deploy:
      resources:
        limits:
          cpus: '8'      # NWChem uses MPI processes
          memory: 16G    # DFT calculations need RAM
        reservations:
          cpus: '2'
          memory: 4G
    # Shared memory for MPI
    shm_size: 512m

  caddy:
    deploy:
      resources:
        limits:
          cpus: '0.5'
          memory: 128M
```

### NWChem-Specific Requirements

**Shared memory:**
- NWChem MPI communication requires shared memory
- Docker default `/dev/shm` is 64MB (too small)
- Set `shm_size: 256m` or higher for MPI

**MPI process count:**
- Configure via `NWCHEM_NPROC` environment variable
- Should match worker CPU limit
- Current code uses `processes=4` hardcoded; make configurable

**Memory per process:**
- DFT calculations: ~1-2GB per MPI process typical
- Large molecules: May need 4GB+ per process
- Safe default: `cpus * 2GB` total memory limit

## Health Checks

### API Health Check

```yaml
services:
  api:
    healthcheck:
      test: ["CMD", "curl", "-f", "http://localhost:8000/health"]
      interval: 30s
      timeout: 10s
      retries: 3
      start_period: 10s
```

Leverages existing `/health` endpoint which returns `{"status": "alive"}`.

### Worker Health Check

Worker doesn't expose HTTP. Options:

**Option 1: Process-based check**
```yaml
services:
  worker:
    healthcheck:
      test: ["CMD", "pgrep", "-f", "huey_consumer"]
      interval: 60s
      timeout: 10s
      retries: 3
```

**Option 2: Custom health file (recommended)**
```yaml
services:
  worker:
    healthcheck:
      test: ["CMD", "test", "-f", "/app/data/.worker_heartbeat"]
      interval: 60s
      timeout: 5s
      retries: 3
```

With worker code that touches heartbeat file periodically:
```python
# In Huey signal handler
@huey.periodic_task(crontab(minute='*'))
def worker_heartbeat():
    Path('/app/data/.worker_heartbeat').touch()
```

### Restart Policies

```yaml
services:
  api:
    restart: unless-stopped

  worker:
    restart: unless-stopped

  caddy:
    restart: unless-stopped
```

**Why `unless-stopped`:**
- Automatic restart on crash
- Survives host reboot (if Docker starts on boot)
- Manual `docker compose stop` is respected

## Image Build Strategy

### Multi-Stage Build for Worker

The worker needs NWChem, CREST, xTB plus Python dependencies. Multi-stage build reduces final image size.

```dockerfile
# Dockerfile.worker

# Stage 1: Build NWChem from source (or use pre-built)
FROM python:3.11-slim AS nwchem-builder
# Install build deps, compile NWChem
# (Alternative: copy from official nwchem image)

# Stage 2: Conda environment with RDKit + Python deps
FROM continuumio/miniconda3 AS conda-builder
COPY environment.yml .
RUN conda env create -f environment.yml

# Stage 3: Final runtime image
FROM python:3.11-slim AS runtime

# Copy NWChem binaries
COPY --from=nwchem-builder /usr/local/bin/nwchem /usr/local/bin/
COPY --from=nwchem-builder /usr/local/share/nwchem /usr/local/share/nwchem

# Copy conda environment
COPY --from=conda-builder /opt/conda/envs/qm-nmr-calc /opt/conda/envs/qm-nmr-calc

# Install MPI runtime (not full dev tools)
RUN apt-get update && apt-get install -y --no-install-recommends \
    openmpi-bin \
    libopenmpi3 \
    && rm -rf /var/lib/apt/lists/*

# Copy application code
COPY src /app/src
COPY pyproject.toml /app/

WORKDIR /app
ENV PATH="/opt/conda/envs/qm-nmr-calc/bin:$PATH"

CMD ["huey_consumer", "qm_nmr_calc.queue.huey", "-w", "1", "-k", "process"]
```

### Simpler Alternative: Use Official NWChem Docker Image

```dockerfile
# Dockerfile.worker (simpler approach)
FROM ghcr.io/nwchemgit/nwchem-dev/amd64 AS nwchem-base

# Install Python and app dependencies on top
RUN apt-get update && apt-get install -y \
    python3.11 \
    python3.11-venv \
    python3-pip

# Continue with app installation...
```

**Trade-offs:**
- Official image: Easier, but larger (includes build tools)
- Multi-stage: Smaller final image, more complex build

### API Dockerfile (Simpler)

```dockerfile
# Dockerfile.api
FROM python:3.11-slim

# Install system deps for RDKit
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
COPY pyproject.toml .
RUN pip install --no-cache-dir .

# Copy application code
COPY src /app/src

WORKDIR /app
EXPOSE 8000

# Health check
HEALTHCHECK --interval=30s --timeout=10s --retries=3 \
    CMD curl -f http://localhost:8000/health || exit 1

CMD ["uvicorn", "qm_nmr_calc.api.app:app", "--host", "0.0.0.0", "--port", "8000"]
```

**Note:** API image may also need RDKit for SMILES validation. Two options:
1. Install RDKit via conda in API container too
2. Move SMILES validation to worker (simplifies API image)

## Build Order for Dockerfiles

### Recommended Development Order

1. **Dockerfile.worker** (most complex)
   - NWChem installation/configuration
   - CREST/xTB installation
   - RDKit via conda
   - MPI configuration
   - Test with simple NWChem calculation

2. **Dockerfile.api** (simpler)
   - Python + FastAPI + uvicorn
   - RDKit for validation (or defer to worker)
   - Health check endpoint verification

3. **Caddyfile** (configuration only)
   - Uses official Caddy image
   - Just needs configuration file

4. **docker-compose.yml** (orchestration)
   - Wire services together
   - Volume mounts
   - Network configuration
   - Environment variables

### Layer Caching Strategy

**Worker Dockerfile layers (slowest to fastest changing):**
```dockerfile
# Layer 1: Base image + system deps (rarely changes)
FROM python:3.11-slim
RUN apt-get update && apt-get install -y openmpi-bin ...

# Layer 2: NWChem/CREST binaries (rarely changes)
COPY --from=nwchem-builder ...

# Layer 3: Python dependencies (changes with requirements)
COPY pyproject.toml .
RUN pip install .

# Layer 4: Application code (changes frequently)
COPY src /app/src
```

**Benefits:**
- Dependency changes don't rebuild NWChem
- Code changes don't reinstall dependencies
- Faster iteration during development

## Environment Configuration

### Production .env File

```bash
# .env.production

# Domain and HTTPS
DOMAIN=nmr.example.com
ACME_EMAIL=admin@example.com

# Worker configuration
NWCHEM_NPROC=4
HUEY_WORKERS=1

# Email notifications (optional)
SMTP_HOST=smtp.example.com
SMTP_PORT=587
SMTP_USER=notifications@example.com
SMTP_PASSWORD=secret

# Feature flags
CREST_ENABLED=true
```

### Development .env File

```bash
# .env.development

# Local development (no HTTPS)
DOMAIN=localhost

# Reduced resources for dev
NWCHEM_NPROC=2
HUEY_WORKERS=1

# Skip email
SMTP_HOST=

# Debug mode
LOG_LEVEL=debug
```

### Docker Compose Environment Handling

```yaml
services:
  api:
    env_file:
      - .env
    environment:
      - LOG_LEVEL=${LOG_LEVEL:-info}

  worker:
    env_file:
      - .env
    environment:
      - NWCHEM_NPROC=${NWCHEM_NPROC:-4}
      - OMP_STACKSIZE=2G
      - GFORTRAN_UNBUFFERED_ALL=1
```

## Deployment Workflow

### Initial Deployment

```bash
# 1. Clone repository
git clone https://github.com/user/qm-nmr-calc.git
cd qm-nmr-calc

# 2. Create environment file
cp .env.example .env
# Edit .env with domain, email, etc.

# 3. Build and start
docker compose build
docker compose up -d

# 4. Verify
docker compose ps
docker compose logs -f
curl https://your-domain.com/health
```

### Updates

```bash
# Pull new code
git pull

# Rebuild and restart
docker compose build
docker compose up -d

# Or for zero-downtime (if using replicas)
docker compose up -d --no-deps --build api
```

### Backup and Restore

```bash
# Backup data volume
docker run --rm -v qm-nmr-calc-data:/data -v $(pwd):/backup \
    alpine tar czf /backup/qm-data-backup.tar.gz /data

# Restore
docker run --rm -v qm-nmr-calc-data:/data -v $(pwd):/backup \
    alpine tar xzf /backup/qm-data-backup.tar.gz -C /
```

## Anti-Patterns to Avoid

### 1. Running as Root

**Don't:**
```dockerfile
# Runs as root by default
CMD ["uvicorn", "app:app"]
```

**Do:**
```dockerfile
RUN useradd -m appuser
USER appuser
CMD ["uvicorn", "app:app"]
```

**Why:** Security. Container escape vulnerabilities are less severe with non-root user.

### 2. Storing Secrets in Image

**Don't:**
```dockerfile
ENV SMTP_PASSWORD=mysecret
```

**Do:**
```yaml
# docker-compose.yml
services:
  api:
    env_file:
      - .env  # Not committed to git
```

**Why:** Images are often pushed to registries. Secrets in images are exposed.

### 3. Ignoring Health Checks

**Don't:**
```yaml
services:
  api:
    # No healthcheck defined
    restart: always
```

**Do:**
```yaml
services:
  api:
    healthcheck:
      test: ["CMD", "curl", "-f", "http://localhost:8000/health"]
      interval: 30s
    restart: unless-stopped
```

**Why:** Without health checks, Docker can't distinguish crashed from starting. `restart: always` may restart-loop forever.

### 4. Hardcoding Paths

**Don't:**
```python
huey = SqliteHuey(filename='/home/chris/dev/qm-nmr-calc/data/huey.db')
```

**Do:**
```python
huey = SqliteHuey(filename='./data/huey.db')
# Working directory set by WORKDIR in Dockerfile
```

**Why:** Paths differ between development and container environments.

### 5. Single Monolithic Container

**Don't:**
```dockerfile
# One container with supervisor running uvicorn + huey + caddy
RUN apt-get install supervisor
COPY supervisord.conf /etc/
CMD ["supervisord"]
```

**Do:**
```yaml
# Separate services
services:
  api:
    command: uvicorn ...
  worker:
    command: huey_consumer ...
  caddy:
    image: caddy:2-alpine
```

**Why:**
- Harder to debug (logs mixed)
- Can't scale independently
- Restart affects all processes
- Defeats Docker's process model

## Roadmap Implications

### Suggested Phase Structure

Based on this architecture research, the Docker deployment milestone should proceed in this order:

**Phase 1: Worker Container (Most Complex)**
- Dockerfile.worker with NWChem + MPI
- CREST/xTB installation
- Verify calculations work in container
- Establish base image for reuse

**Phase 2: API Container**
- Dockerfile.api with FastAPI
- RDKit installation decision
- Health check verification

**Phase 3: Shared Volume Setup**
- SQLite WAL mode configuration
- Volume mount testing
- File permissions verification

**Phase 4: Docker Compose Integration**
- Service orchestration
- Network configuration
- Environment variable handling

**Phase 5: Caddy + HTTPS**
- Caddyfile configuration
- Let's Encrypt integration
- Production domain setup

**Phase 6: CI/CD + Registry**
- GitHub Actions workflow
- GHCR image publishing
- Version tagging strategy

### Research Flags for Phases

| Phase | Needs Deeper Research? | Notes |
|-------|------------------------|-------|
| Worker Container | YES | NWChem Dockerfile specifics, MPI tuning |
| API Container | NO | Standard FastAPI Docker patterns |
| Shared Volume | MAYBE | SQLite locking edge cases under load |
| Docker Compose | NO | Standard orchestration patterns |
| Caddy + HTTPS | NO | Well-documented, simple config |
| CI/CD | MAYBE | GHCR authentication, multi-arch builds |

## Confidence Assessment

| Area | Confidence | Reason |
|------|------------|--------|
| Container split | HIGH | Standard microservices pattern, matches existing process model |
| SQLite sharing | HIGH | Verified SQLite supports this; WAL mode well-documented |
| Volume mounts | HIGH | Standard Docker volume patterns |
| NWChem container | MEDIUM | Official images exist, but MPI config needs verification |
| Health checks | HIGH | Existing endpoints, standard Docker patterns |
| Caddy config | HIGH | Well-documented, simple reverse proxy case |

## Sources

### Docker Compose and FastAPI
- [FastAPI in Containers - Docker - FastAPI](https://fastapi.tiangolo.com/deployment/docker/)
- [Docker Compose Volumes Guide - devopsroles.com](https://www.devopsroles.com/docker-compose-volumes-a-comprehensive-guide/)
- [Better Stack FastAPI Docker Best Practices](https://betterstack.com/community/guides/scaling-python/fastapi-docker-best-practices/)

### SQLite Concurrency
- [SQLite Write-Ahead Logging](https://sqlite.org/wal.html)
- [SQLite WAL mode with multiple processes - SQLite Forum](https://sqlite.org/forum/forumpost/c4dbf6ca17)
- [SQLite performance tuning - phiresky](https://phiresky.github.io/blog/2020/sqlite-performance-tuning/)
- [How SQLite Scales Read Concurrency - Fly.io](https://fly.io/blog/sqlite-internals-wal/)

### NWChem Docker
- [NWChem Containers Wiki](https://github.com/nwchemgit/nwchem/wiki/Containers)
- [nwchem-dockerfiles Repository](https://github.com/nwchemgit/nwchem-dockerfiles)
- [marcindulak/nwchem-openmpi Docker Hub](https://hub.docker.com/r/marcindulak/nwchem-openmpi)

### Huey Task Queue
- [Huey GitHub Repository](https://github.com/coleifer/huey)
- [Huey Documentation](https://huey.readthedocs.io/)

### Health Checks
- [Implementing Health Checks for FastAPI with Docker](https://medium.com/@ntjegadeesh/implementing-health-checks-and-auto-restarts-for-fastapi-applications-using-docker-and-4245aab27ece)
- [FastAPI Health Check Endpoint Example](https://www.index.dev/blog/how-to-implement-health-check-in-python)
