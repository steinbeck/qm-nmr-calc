# Phase 37: Docker Compose Integration - Research

**Researched:** 2026-02-02
**Domain:** Docker Compose, multi-container orchestration, persistent storage, graceful shutdown
**Confidence:** HIGH

## Summary

This phase integrates the existing Dockerfile.worker and Dockerfile.api into a complete Docker Compose deployment. The research reveals a critical signal handling requirement: Huey uses SIGINT for graceful shutdown (allowing tasks to complete), while Docker sends SIGTERM by default. The compose file must configure `stop_signal: SIGINT` for the worker service.

The standard approach is a single `docker-compose.yml` with two services (api and worker) sharing a named volume for the `/app/data` directory. Both services use the same codebase but different Dockerfiles. SQLite (Huey's queue database) works correctly across multiple containers sharing a local Docker volume - the kernel's file locking works properly on local filesystems.

Key requirements addressed: named volumes for persistence, `unless-stopped` restart policy, health checks (existing HEALTHCHECK directives preserved), `.env` configuration, extended `stop_grace_period` for long-running NMR calculations, and `shm_size: 512m` for MPI shared memory.

**Primary recommendation:** Create docker-compose.yml with `stop_signal: SIGINT` and `stop_grace_period: 300s` for worker to enable graceful task completion.

## Standard Stack

The established tools for this domain:

### Core
| Tool | Version | Purpose | Why Standard |
|------|---------|---------|--------------|
| Docker Compose | v2+ | Multi-container orchestration | Native Docker tooling, Compose Specification standard |
| Named Volumes | - | Persistent data storage | Survives container restarts, managed by Docker |
| .env files | - | Environment configuration | Standard Docker Compose pattern, keeps secrets out of compose file |

### Supporting
| Tool | Version | Purpose | When to Use |
|------|---------|---------|-------------|
| Compose healthchecks | - | Service dependency ordering | `depends_on: condition: service_healthy` |
| stop_signal | - | Custom shutdown signal | SIGINT for Huey graceful shutdown |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Named volumes | Bind mounts | Bind mounts good for dev, named volumes better for production portability |
| SQLite Huey | Redis Huey | Redis adds complexity; SQLite works fine for single-host deployment |
| unless-stopped | always | unless-stopped respects manual stops; always restarts even after `docker compose stop` |

## Architecture Patterns

### Recommended Project Structure
```
project/
├── docker-compose.yml      # Main compose file
├── .env                    # Local environment overrides (gitignored)
├── .env.example            # Documented configuration template
├── Dockerfile.api          # API container (existing)
├── Dockerfile.worker       # Worker container (existing)
└── data/                   # Volume mount point (gitignored)
    ├── huey.db             # Huey queue SQLite database
    └── jobs/               # Job directories
```

### Pattern 1: Shared Volume for API + Worker
**What:** Both containers mount the same named volume at `/app/data`
**When to use:** Always - API writes jobs, Worker reads/writes results, both need Huey DB
**Example:**
```yaml
# Source: Docker Compose official docs
services:
  api:
    volumes:
      - app-data:/app/data

  worker:
    volumes:
      - app-data:/app/data

volumes:
  app-data:
```

### Pattern 2: Service Health Dependencies
**What:** Worker waits for API to be healthy before starting
**When to use:** Ensures data directory initialized before worker attempts access
**Example:**
```yaml
# Source: Docker Compose official docs
services:
  worker:
    depends_on:
      api:
        condition: service_healthy
```

### Pattern 3: Graceful Shutdown with Custom Signal
**What:** Configure stop_signal and stop_grace_period for task completion
**When to use:** Worker containers with long-running tasks
**Example:**
```yaml
# Source: Docker Compose + Huey documentation
services:
  worker:
    stop_signal: SIGINT       # Huey graceful shutdown
    stop_grace_period: 300s   # 5 minutes for task completion
```

### Anti-Patterns to Avoid
- **SIGTERM for Huey:** SIGTERM immediately kills Huey, losing in-progress tasks. Use SIGINT instead.
- **Short stop_grace_period:** NMR calculations take minutes. 10s default causes task loss.
- **Hardcoded configuration:** Use .env files for all configurable values (ports, worker count, etc.)
- **Anonymous volumes:** Named volumes are easier to manage, backup, and restore.
- **Root ownership conflicts:** API runs as non-root (UID 999), must handle volume permissions.

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Service restart on failure | Custom supervisor | `restart: unless-stopped` | Docker handles it natively |
| Health check logic | Custom HTTP endpoints | Docker HEALTHCHECK | Already in Dockerfiles |
| Signal forwarding | Wrapper scripts | `init: true` or proper signal handling | Huey handles SIGINT natively |
| Volume permissions | Manual chown in entrypoint | Proper Dockerfile USER setup | Already done in API Dockerfile |
| Environment config | Custom config loader | Docker Compose .env support | Standard pattern |

**Key insight:** Docker Compose provides all orchestration primitives needed. The complexity is in correct configuration, not custom code.

## Common Pitfalls

### Pitfall 1: Wrong Shutdown Signal for Huey
**What goes wrong:** Task in progress is killed, job left in 'running' state forever
**Why it happens:** Docker sends SIGTERM by default, Huey interprets SIGTERM as "stop immediately"
**How to avoid:** Configure `stop_signal: SIGINT` in docker-compose.yml
**Warning signs:** Jobs stuck in 'running' status after container restart
**Source:** [Huey Consumer Documentation](https://huey.readthedocs.io/en/latest/consumer.html)

### Pitfall 2: Insufficient stop_grace_period
**What goes wrong:** Docker sends SIGKILL before task completes
**Why it happens:** Default is 10 seconds, NMR calculations take 5-30 minutes
**How to avoid:** Set `stop_grace_period: 300s` or longer based on max expected task duration
**Warning signs:** Container logs show "Received SIGKILL" before task completion

### Pitfall 3: Volume Permission Mismatch
**What goes wrong:** API container (UID 999) can't write to volume initially owned by root
**Why it happens:** Named volumes default to root ownership
**How to avoid:** API Dockerfile already creates data dir with correct ownership; worker runs as root
**Warning signs:** "Permission denied" errors on first write
**Source:** [Docker Volume Permissions](https://denibertovic.com/posts/handling-permissions-with-docker-volumes/)

### Pitfall 4: MPI Shared Memory Too Small
**What goes wrong:** NWChem crashes with SIGBUS or "not enough shared memory"
**Why it happens:** Docker default /dev/shm is 64MB, MPI needs more
**How to avoid:** Configure `shm_size: 512m` for worker container
**Warning signs:** SIGBUS errors in NWChem output
**Source:** [Intel MPI Docker Limitation](https://www.intel.com/content/www/us/en/docs/mpi-library/developer-guide-linux/2021-6/problem-mpi-limitation-for-docker.html)

### Pitfall 5: SQLite on Network Filesystem
**What goes wrong:** Database corruption, "database is locked" errors
**Why it happens:** SQLite byte-range locking doesn't work on NFS/CIFS
**How to avoid:** Use local Docker volumes only (not NFS-backed)
**Warning signs:** Intermittent database errors, corrupted huey.db
**Note:** Local Docker named volumes work fine for SQLite concurrent access

### Pitfall 6: Missing .env File
**What goes wrong:** Compose fails or uses wrong defaults
**Why it happens:** .env is gitignored, not present after fresh clone
**How to avoid:** Provide .env.example with all variables documented; fail fast with clear error
**Warning signs:** Unexpected default values, missing configuration

## Code Examples

Verified patterns from official sources:

### Complete docker-compose.yml Structure
```yaml
# Docker Compose Specification (v2+)
# Source: https://docs.docker.com/compose/compose-file/

services:
  api:
    build:
      context: .
      dockerfile: Dockerfile.api
    ports:
      - "${API_PORT:-8000}:8000"
    volumes:
      - app-data:/app/data
    environment:
      - PYTHONUNBUFFERED=1
    env_file:
      - path: .env
        required: false
    restart: unless-stopped
    # HEALTHCHECK already in Dockerfile.api
    healthcheck:
      test: ["CMD", "curl", "-f", "http://localhost:8000/health"]
      interval: 30s
      timeout: 5s
      start_period: 10s
      retries: 3

  worker:
    build:
      context: .
      dockerfile: Dockerfile.worker
    volumes:
      - app-data:/app/data
    environment:
      - PYTHONUNBUFFERED=1
      - NWCHEM_NPROC=${NWCHEM_NPROC:-4}
      - OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
    env_file:
      - path: .env
        required: false
    depends_on:
      api:
        condition: service_healthy
    restart: unless-stopped
    # CRITICAL: Huey uses SIGINT for graceful shutdown
    stop_signal: SIGINT
    stop_grace_period: 300s  # 5 minutes for long calculations
    shm_size: 512m           # MPI shared memory requirement
    # Worker healthcheck via process presence
    healthcheck:
      test: ["CMD", "pgrep", "-f", "huey_consumer"]
      interval: 30s
      timeout: 5s
      start_period: 30s
      retries: 3

volumes:
  app-data:
    name: qm-nmr-calc-data
```

### .env.example Template
```bash
# Source: Docker Compose env_file documentation
# https://docs.docker.com/compose/how-tos/environment-variables/

# ============================================================================
# API Configuration
# ============================================================================
# Port to expose API on host (default: 8000)
API_PORT=8000

# ============================================================================
# Worker Configuration
# ============================================================================
# Number of MPI processes for NWChem calculations (default: 4)
# Increase for faster calculations on multi-core systems
NWCHEM_NPROC=4

# Number of OpenMP threads for CREST/xTB (default: 4)
OMP_NUM_THREADS=4

# ============================================================================
# Email Notifications (optional)
# ============================================================================
# SMTP server for sending job completion notifications
# Leave empty to disable email notifications
SMTP_HOST=
SMTP_PORT=587
SMTP_USER=
SMTP_PASS=
SMTP_FROM=

# ============================================================================
# Advanced Options
# ============================================================================
# Stop grace period in seconds (how long to wait for jobs to complete)
# Default is 300 (5 minutes). Increase for very long calculations.
# STOP_GRACE_PERIOD=300
```

### Health Check for Non-HTTP Worker
```yaml
# Source: Docker healthcheck best practices
# For containers without HTTP endpoint, check process presence
healthcheck:
  test: ["CMD", "pgrep", "-f", "huey_consumer"]
  interval: 30s
  timeout: 5s
  start_period: 30s
  retries: 3
```

### depends_on with Health Condition
```yaml
# Source: Docker Compose official docs
# Worker waits for API to be healthy
services:
  worker:
    depends_on:
      api:
        condition: service_healthy
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `version: "3.x"` in compose file | No version required | Compose Specification (2020+) | Omit version for modern compose |
| `links:` for service discovery | Service names as hostnames | Docker Compose v2+ | Services communicate by name |
| `mem_limit` | `deploy.resources.limits.memory` | Compose Specification | Legacy syntax still works |
| Build-time health checks | Runtime HEALTHCHECK directive | Docker 1.12+ | Health checks run during container lifetime |

**Current:**
- Compose Specification is the standard (no `version:` key needed)
- `depends_on` with `condition: service_healthy` is the correct way to order service startup
- Named volumes declared at top level survive `docker compose down`
- `init: true` for proper signal handling (optional, Huey handles signals well)

## Open Questions

Things that couldn't be fully resolved:

1. **Email notification configuration in containers**
   - What we know: SMTP settings need to be passed via environment
   - What's unclear: Whether current notification code reads from environment or hardcoded
   - Recommendation: Document in .env.example, verify notification.py reads from env

2. **Optimal NWCHEM_NPROC value**
   - What we know: Default is 4, should match available CPU cores
   - What's unclear: Performance impact of over-provisioning in container
   - Recommendation: Document that users should set based on host CPU cores

3. **Volume backup strategy**
   - What we know: Named volumes can be backed up with `docker run --rm -v qm-nmr-calc-data:/data -v $(pwd):/backup alpine tar cvf /backup/backup.tar /data`
   - What's unclear: Whether automated backup is needed for this phase
   - Recommendation: Document manual backup in README, defer automation

## Sources

### Primary (HIGH confidence)
- [Docker Compose Services Specification](https://docs.docker.com/compose/compose-file/05-services/) - healthcheck, restart, shm_size, stop_signal, stop_grace_period, depends_on, volumes syntax
- [Huey Consumer Documentation](https://huey.readthedocs.io/en/latest/consumer.html) - Signal handling: SIGINT=graceful, SIGTERM=immediate, SIGHUP=restart
- [Docker Compose Volumes](https://docs.docker.com/reference/compose-file/volumes/) - Named volume declaration and usage

### Secondary (MEDIUM confidence)
- [Docker Compose Environment Variables Best Practices](https://docs.docker.com/compose/how-tos/environment-variables/best-practices/) - .env file patterns
- [SQLite with Docker Containers](https://rbranson.medium.com/sharing-sqlite-databases-across-containers-is-surprisingly-brilliant-bacb8d753054) - Confirms SQLite works across containers on same host
- [Intel MPI Docker Limitations](https://www.intel.com/content/www/us/en/docs/mpi-library/developer-guide-linux/2021-6/problem-mpi-limitation-for-docker.html) - shm_size requirements for MPI

### Tertiary (LOW confidence)
- Community patterns for FastAPI + task queue deployments - General architecture confirmed but Huey-specific examples limited

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - Docker Compose is mature, well-documented
- Architecture: HIGH - Patterns verified with official documentation
- Signal handling: HIGH - Verified directly from Huey docs, critical for OPS-01
- Pitfalls: HIGH - Based on official docs and multiple sources

**Research date:** 2026-02-02
**Valid until:** 60 days (Docker Compose is stable, Huey signals unlikely to change)

---

## Appendix: Signal Handling Summary

This is critical for requirement OPS-01 (graceful shutdown):

| Signal | Huey Behavior | Docker Default | Required Config |
|--------|---------------|----------------|-----------------|
| SIGINT | Graceful - finish current task | Not sent | `stop_signal: SIGINT` |
| SIGTERM | Immediate - interrupt task | Sent on stop | Override with SIGINT |
| SIGHUP | Clean restart | Not used | N/A |
| SIGKILL | Force kill | Sent after grace period | Extend grace period |

**The project already handles SIGNAL_INTERRUPTED** in `queue.py` to mark jobs as failed if interrupted. However, using SIGINT prevents interruption entirely by letting tasks complete.
