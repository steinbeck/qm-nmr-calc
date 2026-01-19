---
phase: 04-results-delivery
plan: 01
subsystem: api
tags: [fastapi, downloads, xyz, sdf, zip, results]

# Dependency graph
requires:
  - phase: 03-nmr-calculations
    provides: NMR results stored in job status, XYZ geometry files
provides:
  - Results retrieval endpoint returning NMR shifts as JSON
  - Geometry download as XYZ file
  - Geometry download as SDF file with coordinates
  - Raw NWChem output download as ZIP archive
affects: [04-notifications, testing, documentation]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - 409 Conflict for incomplete job state (vs 404 for missing)
    - FileResponse for static file downloads
    - In-memory ZIP creation with io.BytesIO

key-files:
  created: []
  modified:
    - src/qm_nmr_calc/storage.py
    - src/qm_nmr_calc/api/routers/jobs.py

key-decisions:
  - "Return 409 Conflict for incomplete jobs (vs 404) - job exists but not ready"
  - "XYZ uses chemical/x-xyz MIME type, SDF uses chemical/x-mdl-sdfile"
  - "SDF generated on-the-fly from SMILES + XYZ coordinates using RDKit"
  - "Output ZIP includes only .out and .nw files from scratch directory"

patterns-established:
  - "Job state validation pattern: 404 not found, 409 not complete, then proceed"
  - "File path helpers in storage.py for output file location"

# Metrics
duration: 3min
completed: 2026-01-19
---

# Phase 4 Plan 1: Results and Download Endpoints Summary

**Four GET endpoints for retrieving NMR results as JSON, downloading optimized geometry as XYZ/SDF, and raw NWChem output as ZIP**

## Performance

- **Duration:** 3 min
- **Started:** 2026-01-19T21:30:44Z
- **Completed:** 2026-01-19T21:33:34Z
- **Tasks:** 3
- **Files modified:** 2

## Accomplishments
- Added `get_geometry_file()` and `get_output_files()` helpers to storage.py
- Created GET /jobs/{job_id}/results endpoint returning NMRResultsResponse JSON
- Created GET /jobs/{job_id}/geometry endpoint for XYZ file download
- Created GET /jobs/{job_id}/geometry.sdf endpoint generating SDF from SMILES + coordinates
- Created GET /jobs/{job_id}/output endpoint for ZIP archive of NWChem files
- All endpoints properly return 409 Conflict for incomplete jobs

## Task Commits

Each task was committed atomically:

1. **Task 1: Add file path helpers to storage.py** - `9c4624e` (feat)
2. **Task 2: Add results and download endpoints to jobs.py** - `ce6528f` (feat)
3. **Task 3: Test endpoints with curl** - verification only, no commit

**Plan metadata:** (pending)

## Files Created/Modified
- `src/qm_nmr_calc/storage.py` - Added get_geometry_file() and get_output_files() helpers
- `src/qm_nmr_calc/api/routers/jobs.py` - Added four new GET endpoints for results and downloads

## Decisions Made
- Use 409 Conflict status for jobs that exist but are not complete (clearer than 404)
- Generate SDF on-the-fly from original SMILES with XYZ coordinates (no separate SDF storage)
- ZIP only includes .out and .nw files (main NWChem outputs) from scratch directory
- Use chemical/x-xyz and chemical/x-mdl-sdfile MIME types for proper file associations

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None - all verification tests passed on first attempt.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- Results retrieval endpoints complete and functional
- Ready for 04-02 (email notifications) which will notify users when these results become available
- All endpoints documented in OpenAPI at /api/v1/openapi.json

---
*Phase: 04-results-delivery*
*Completed: 2026-01-19*
