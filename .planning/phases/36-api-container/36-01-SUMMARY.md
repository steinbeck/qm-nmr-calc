---
phase: 36-api-container
plan: 01
subsystem: infra
tags: [docker, fastapi, uvicorn, uv, multi-stage-build]

# Dependency graph
requires:
  - phase: 35-worker-container
    provides: Docker containerization pattern with validation scripts
provides:
  - Dockerfile.api for FastAPI server container
  - scripts/validate-api.sh for container validation
  - Docker HEALTHCHECK configuration
affects: [37-compose, 38-monitoring, 39-cicd, 40-staging]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - Multi-stage uv build pattern for API container
    - Non-root user (UID 999) for security
    - Docker HEALTHCHECK with curl

key-files:
  created:
    - Dockerfile.api
    - scripts/validate-api.sh
  modified: []

key-decisions:
  - "Image size ~733MB due to RDKit and scientific dependencies (numpy, scipy, pandas)"
  - "Added libxrender1, libxext6, libexpat1 for RDKit drawing support"
  - "Used exec form CMD for proper SIGTERM handling"

patterns-established:
  - "Validation script pattern: Python, uvicorn, app import, static files, templates, data dir, user, health endpoint"
  - "HEALTHCHECK pattern: curl --fail http://localhost:8000/health"

# Metrics
duration: 4min
completed: 2026-02-02
---

# Phase 36 Plan 01: API Container Summary

**Multi-stage Docker build for FastAPI with uv dependency management, non-root user (UID 999), and Docker HEALTHCHECK**

## Performance

- **Duration:** 4 min
- **Started:** 2026-02-02T19:01:25Z
- **Completed:** 2026-02-02T19:05:46Z
- **Tasks:** 3
- **Files modified:** 2

## Accomplishments
- Multi-stage Docker build using `ghcr.io/astral-sh/uv:python3.11-bookworm-slim` builder and `python:3.11-slim-bookworm` runtime
- Non-root user (appuser, UID 999) for container security
- Docker HEALTHCHECK hitting `/health` endpoint with curl
- Comprehensive validation script testing all container components

## Task Commits

Each task was committed atomically:

1. **Task 1 + Task 2: Dockerfile.api and validation script** - `8784f7f` (feat)
   - Tasks combined because Dockerfile requires validation script to exist

**Plan metadata:** (pending)

## Files Created/Modified
- `Dockerfile.api` - Multi-stage build for FastAPI API server
- `scripts/validate-api.sh` - Container validation script testing Python, uvicorn, app import, static files, templates, data directory, non-root user, health endpoint

## Decisions Made

1. **Image size ~733MB** - RDKit and scientific dependencies (numpy, scipy, pandas, statsmodels, matplotlib) are large. Original target of ~200MB was based on research not accounting for full dependency chain. This is acceptable for an NMR calculation application.

2. **Added X11 libraries** - RDKit drawing requires libxrender1, libxext6, libexpat1 for molecule visualization.

3. **exec form CMD** - Using `CMD ["uvicorn", ...]` instead of shell form ensures proper signal handling for container orchestration.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Added missing X11 libraries for RDKit**
- **Found during:** Task 2 (validation script execution)
- **Issue:** RDKit drawing imports failed with missing libXrender.so.1 and libexpat.so.1
- **Fix:** Added libxrender1, libxext6, libexpat1 to apt-get install
- **Files modified:** Dockerfile.api
- **Verification:** App imports successfully, validation script passes
- **Committed in:** 8784f7f

---

**Total deviations:** 1 auto-fixed (blocking)
**Impact on plan:** Essential for RDKit functionality. No scope creep.

## Issues Encountered
- Port 8000 and 8001 were already in use during health check verification; used port 8888 instead

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- API container ready for docker-compose integration
- Health endpoint verified working at `/health`
- Static files and templates accessible
- Ready for Phase 37: Docker Compose

---
*Phase: 36-api-container*
*Completed: 2026-02-02*
