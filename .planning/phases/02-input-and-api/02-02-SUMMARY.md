---
phase: 02-input-and-api
plan: 02
subsystem: api
tags: [fastapi, routers, rest-api, health-checks, job-submission]

# Dependency graph
requires:
  - phase: 02-input-and-api
    plan: 01
    provides: Validation functions and API schemas
provides:
  - Health check endpoints (/health, /health/ready)
  - Job submission endpoints (POST /jobs, POST /jobs/upload)
  - Job status endpoint (GET /jobs/{job_id})
  - Routers ready for mounting in FastAPI app
affects: [02-03 (app.py will mount these routers)]

# Tech tracking
tech-stack:
  added: []
  patterns: [FastAPI router pattern, HTTP 202 Accepted with Location header]

key-files:
  created:
    - src/qm_nmr_calc/api/routers/__init__.py
    - src/qm_nmr_calc/api/routers/health.py
    - src/qm_nmr_calc/api/routers/jobs.py
  modified: []

key-decisions:
  - "Health liveness probe is minimal (just return alive status)"
  - "Readiness probe checks data directory writable and Huey importable"
  - "Jobs router prefix is /jobs, gets /api/v1 prefix when mounted"
  - "Use JSONResponse with explicit headers for 202 responses"

patterns-established:
  - "job_status_to_response() helper for JobStatus to API response conversion"
  - "HTTPException with detail dict for RFC 7807 error format"

# Metrics
duration: 2min
completed: 2026-01-19
---

# Phase 2 Plan 2: Routers Summary

**FastAPI routers for health checks (liveness/readiness) and job endpoints (SMILES submission, file upload, status polling)**

## Performance

- **Duration:** 2 min
- **Started:** 2026-01-19T17:15:19Z
- **Completed:** 2026-01-19T17:16:56Z
- **Tasks:** 2
- **Files created:** 3

## Accomplishments
- Created health router with /health (liveness) and /health/ready (readiness) endpoints
- Created jobs router with POST /jobs (SMILES), POST /jobs/upload (MOL/SDF), GET /jobs/{job_id}
- Implemented job_status_to_response() helper to flatten JobStatus to API response format
- Error responses use RFC 7807 ProblemDetail format via HTTPException detail dict

## Task Commits

Each task was committed atomically:

1. **Task 1: Create health router** - `ac446b6` (feat)
2. **Task 2: Create jobs router** - `746182a` (feat)

## Files Created/Modified
- `src/qm_nmr_calc/api/routers/__init__.py` - Package initialization
- `src/qm_nmr_calc/api/routers/health.py` - Liveness and readiness endpoints
- `src/qm_nmr_calc/api/routers/jobs.py` - Job submission and status endpoints

## Decisions Made
- Health liveness just returns {"status": "alive"} - minimal for load balancer probes
- Readiness checks data directory writable and Huey queue importable
- Use JSONResponse for 202 responses to add Location and Retry-After headers
- Convert helper function flattens JobStatus input.smiles/name to input_smiles/input_name

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- Health and jobs routers ready to be mounted in FastAPI app
- Ready for 02-03 (app.py creation and integration)

---
*Phase: 02-input-and-api*
*Completed: 2026-01-19*
