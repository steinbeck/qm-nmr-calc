---
phase: 06-web-ui
plan: 02
subsystem: ui
tags: [jinja2, forms, polling, javascript, fastapi]

# Dependency graph
requires:
  - phase: 06-01
    provides: Template infrastructure with Jinja2, Pico CSS, base.html
provides:
  - Submission form with SMILES/file input and validation
  - Status page with auto-polling and redirect
  - Web routes for /submit, /status/{job_id}
affects: [06-03, 06-04]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - Form data preservation on validation error
    - JavaScript mutual exclusion for inputs
    - Polling with setTimeout recursion
    - 303 redirect after POST

key-files:
  created:
    - src/qm_nmr_calc/api/templates/status.html
  modified:
    - src/qm_nmr_calc/api/templates/submit.html
    - src/qm_nmr_calc/api/routers/web.py

key-decisions:
  - "Form preserves values on validation error via form_data context"
  - "JavaScript disables opposite input when one has value"
  - "3-second polling interval with backoff on errors"
  - "303 See Other redirect after successful POST"

patterns-established:
  - "_get_form_context() helper for consistent form data"
  - "render_error() nested function for validation error responses"
  - "Polling pattern: async fetch -> update UI -> setTimeout recurse"

# Metrics
duration: 4min
completed: 2026-01-20
---

# Phase 6 Plan 2: Submission and Status Pages Summary

**Form submission with SMILES/file validation, status page with 3-second polling and auto-redirect on completion**

## Performance

- **Duration:** 4 min
- **Started:** 2026-01-20T12:13:18Z
- **Completed:** 2026-01-20T12:17:37Z
- **Tasks:** 3
- **Files modified:** 3

## Accomplishments

- Full submission form with SMILES text input and MOL/SDF file upload
- JavaScript mutual exclusion preventing both inputs from being used
- Form validation with error display and preserved input values
- Status page with live elapsed time counter
- Polling API every 3 seconds with auto-redirect on completion
- Error display on job failure with no further polling

## Task Commits

Each task was committed atomically:

1. **Task 1: Create submission form template** - `a5421b9` (feat)
2. **Task 2: Create status page template with polling** - `c91372e` (feat)
3. **Task 3: Add submit and status routes to web router** - `9ebe623` (feat)

## Files Created/Modified

- `src/qm_nmr_calc/api/templates/submit.html` - Full form with SMILES, file, solvent, preset, name, email fields
- `src/qm_nmr_calc/api/templates/status.html` - Status display with polling JavaScript
- `src/qm_nmr_calc/api/routers/web.py` - POST /submit and GET /status/{job_id} routes

## Decisions Made

- **Form data preservation:** On validation error, form values are passed back to template to avoid user re-entry
- **Input mutual exclusion:** JavaScript disables file input when SMILES has value and vice versa
- **Poll interval:** 3 seconds chosen as balance between responsiveness and server load
- **Error backoff:** On fetch failure, retry interval doubles (6 seconds) to handle temporary issues
- **303 redirect:** Using 303 See Other after POST ensures browser won't resubmit on refresh

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- Submission form fully functional at /
- Status page with polling at /status/{job_id}
- Results page placeholder ready for Plan 06-03
- All validation errors displayed with preserved form values

---
*Phase: 06-web-ui*
*Completed: 2026-01-20*
