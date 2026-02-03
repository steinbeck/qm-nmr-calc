---
phase: 37
plan: 01
subsystem: deployment
tags: [docker-compose, orchestration, graceful-shutdown, health-checks]

dependency_graph:
  requires: [35-worker-container, 36-api-container]
  provides: [single-command-deployment, shared-volume-orchestration, graceful-worker-shutdown]
  affects: [38-production-hardening, 39-ci-cd]

tech_stack:
  added: []
  patterns: [compose-specification, named-volumes, service-health-dependencies]

key_files:
  created:
    - docker-compose.yml
    - .env.example
  modified:
    - .gitignore

decisions:
  - id: SIGINT-shutdown
    choice: "stop_signal: SIGINT for worker"
    reason: "Huey uses SIGINT for graceful task completion; SIGTERM causes immediate kill"
  - id: grace-period
    choice: "stop_grace_period: 300s"
    reason: "NMR calculations can take 5-30 minutes; allow completion before SIGKILL"
  - id: shm-size
    choice: "shm_size: 512m"
    reason: "MPI requires larger shared memory than Docker's 64MB default"

metrics:
  duration: 2m
  completed: 2026-02-03
---

# Phase 37 Plan 01: Docker Compose Configuration Summary

**One-liner:** Multi-service orchestration with SIGINT shutdown, 5-min grace period, and shared volume for API/worker communication.

## What Was Built

Created production-ready Docker Compose configuration for single-command deployment of the qm-nmr-calc stack.

### docker-compose.yml

Two-service orchestration:

1. **api** service:
   - Builds from Dockerfile.api
   - Configurable port via API_PORT env var
   - Health check on /health endpoint
   - Restart policy: unless-stopped

2. **worker** service:
   - Builds from Dockerfile.worker
   - SIGINT signal for graceful Huey shutdown
   - 5-minute stop_grace_period for long calculations
   - 512MB shm_size for MPI shared memory
   - Configurable memory limit (default 8GB)
   - Depends on API being healthy before starting

Shared infrastructure:
- Named volume `qm-nmr-calc-data` for /app/data (Huey DB + job files)
- Optional .env file for configuration overrides

### .env.example

Documented configuration template with:
- API_PORT for custom port mapping
- NWCHEM_NPROC for MPI process count
- OMP_NUM_THREADS for CREST/xTB parallelization
- WORKER_MEMORY_LIMIT for container resources
- SMTP settings for future email notifications
- Clear recommendations for optimal settings

### .gitignore Update

Added `.env` to prevent accidental commit of secrets while keeping `.env.example` tracked.

## Decisions Made

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Shutdown signal | SIGINT | Huey graceful shutdown; SIGTERM kills immediately |
| Grace period | 300s | NMR calculations take minutes; allow completion |
| Shared memory | 512MB | MPI requires more than Docker's 64MB default |
| Volume naming | qm-nmr-calc-data | Explicit naming for backup/management clarity |
| Restart policy | unless-stopped | Respects manual stops; survives host reboot |

## Deviations from Plan

None - plan executed exactly as written.

## Commits

| Commit | Type | Description |
|--------|------|-------------|
| eb7db12 | feat | docker-compose.yml with full orchestration config |
| 275f045 | docs | .env.example configuration template |
| 7c2d8f2 | chore | .gitignore update for .env |

## Testing Notes

- `docker compose config` validates syntax successfully
- All required fields verified present (stop_signal, stop_grace_period, shm_size, healthchecks)
- Ready for integration testing with full stack

## Next Phase Readiness

**Ready for Phase 38 (Production Hardening):**
- Base compose configuration complete
- Health checks in place for monitoring
- Resource limits configurable
- Graceful shutdown properly configured

No blockers identified.
