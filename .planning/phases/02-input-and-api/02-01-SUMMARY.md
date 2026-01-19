---
phase: 02-input-and-api
plan: 01
subsystem: api
tags: [fastapi, pydantic, rdkit, validation, schemas]

# Dependency graph
requires:
  - phase: 01-foundation
    provides: Job storage model (JobStatus, JobInput) and directory management
provides:
  - FastAPI + uvicorn + python-multipart dependencies installed
  - Molecule validation functions (SMILES and MOL/SDF)
  - API request/response Pydantic schemas
  - JobInput model extended with name field
affects: [02-02 (routers will use validation and schemas), 02-03 (health endpoints)]

# Tech tracking
tech-stack:
  added: [fastapi, uvicorn, python-multipart]
  patterns: [RFC 7807 error responses, tuple return (success, error) pattern]

key-files:
  created:
    - src/qm_nmr_calc/validation.py
    - src/qm_nmr_calc/api/__init__.py
    - src/qm_nmr_calc/api/schemas.py
    - README.md
  modified:
    - pyproject.toml
    - uv.lock
    - src/qm_nmr_calc/models.py
    - src/qm_nmr_calc/storage.py

key-decisions:
  - "Tuple return pattern for validation: (mol, None) success, (None, error) failure"
  - "Detect SDF vs MOL by $$$$ delimiter presence in content"
  - "Flattened input fields in JobStatusResponse (input_smiles, input_name) for simpler API"
  - "RFC 7807 ProblemDetail schema for standardized error responses"

patterns-established:
  - "Validation functions return (success_value, None) or (None, error_message)"
  - "API schemas use Field() with descriptions and examples for OpenAPI docs"

# Metrics
duration: 4min
completed: 2026-01-19
---

# Phase 2 Plan 1: Validation and Schemas Summary

**FastAPI dependencies installed, RDKit-based SMILES/MOL/SDF validation module, and Pydantic API schemas with RFC 7807 error model**

## Performance

- **Duration:** 4 min
- **Started:** 2026-01-19T17:10:00Z
- **Completed:** 2026-01-19T17:14:00Z
- **Tasks:** 3
- **Files modified:** 8

## Accomplishments
- Installed FastAPI 0.128.0, uvicorn 0.40.0, python-multipart 0.0.21
- Created validation module with validate_smiles() and validate_mol_file() functions
- Created API schemas: JobSubmitRequest, JobStatusResponse, ProblemDetail
- Extended JobInput model with optional name field for molecule labels

## Task Commits

Each task was committed atomically:

1. **Task 1: Install FastAPI dependencies** - `ee12331` (chore)
2. **Task 2: Create molecule validation module** - `2ca99b8` (feat)
3. **Task 3: Create API schemas** - `b14f829` (feat)

## Files Created/Modified
- `src/qm_nmr_calc/validation.py` - SMILES and MOL/SDF validation using RDKit
- `src/qm_nmr_calc/api/__init__.py` - API package initialization
- `src/qm_nmr_calc/api/schemas.py` - Request/response Pydantic models
- `src/qm_nmr_calc/models.py` - Added name field to JobInput
- `src/qm_nmr_calc/storage.py` - Updated create_job_directory() to accept name
- `pyproject.toml` - Added fastapi, uvicorn, python-multipart dependencies
- `uv.lock` - Updated lockfile
- `README.md` - Created to fix hatchling build error

## Decisions Made
- Used tuple return pattern (mol, None) / (None, error) for validation functions - clear success/failure without exceptions
- Detect SDF vs MOL files by checking for $$$$ delimiter in content
- Flattened input fields in JobStatusResponse (input_smiles, input_name) for simpler REST API response
- Used RFC 7807 ProblemDetail schema for standardized machine-readable errors

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Created README.md to fix hatchling build**
- **Found during:** Task 1 (Install FastAPI dependencies)
- **Issue:** pyproject.toml references README.md which didn't exist, causing build failure
- **Fix:** Created minimal README.md with project description and basic usage
- **Files modified:** README.md
- **Verification:** `uv add` command succeeded after creating file
- **Committed in:** ee12331 (Task 1 commit)

---

**Total deviations:** 1 auto-fixed (1 blocking)
**Impact on plan:** Necessary to unblock build. No scope creep.

## Issues Encountered
- RDKit stderr capture doesn't work for C-level output (RDKit writes errors at C level, not Python sys.stderr). Fallback generic error message "Invalid SMILES string" used when specific error unavailable. This is documented as open question in RESEARCH.md.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- Validation functions ready for router integration
- API schemas ready for endpoint implementation
- FastAPI installed and importable
- Ready for 02-02 (job submission and status endpoints)

---
*Phase: 02-input-and-api*
*Completed: 2026-01-19*
