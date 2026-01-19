# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-01-19)

**Core value:** Reliable async NMR predictions with full control over calculation parameters -- submit a molecule, get back accurate 1H/13C shifts without babysitting long-running calculations.
**Current focus:** Phase 3 - Calculation Pipeline (Phase 2 complete)

## Current Position

Phase: 2 of 6 (Input and API) - COMPLETE
Plan: 3 of 3 in phase (complete)
Status: Phase 2 complete, ready for Phase 3
Last activity: 2026-01-19 -- Completed 02-03-PLAN.md (Application Assembly)

Progress: [██████░░░░] ~30%

## Performance Metrics

**Velocity:**
- Total plans completed: 6
- Average duration: 2.8 min
- Total execution time: 17 min

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 01-foundation | 3 | 9 min | 3.0 min |
| 02-input-and-api | 3 | 8 min | 2.7 min |

**Recent Trend:**
- Last 5 plans: 01-03 (2 min), 02-01 (4 min), 02-02 (2 min), 02-03 (2 min)
- Trend: Stable

*Updated after each plan completion*

## Accumulated Context

### Decisions

Decisions are logged in PROJECT.md Key Decisions table.
Recent decisions affecting current work:

- Use hatchling build system for src layout (uv compatible)
- 12-character hex job IDs from uuid4 for URL-safe identifiers
- orjson with OPT_INDENT_2 for human-readable status.json files
- ISiCLE installed as editable from local fork
- SqliteHuey with fsync=True for crash-safe job queue
- Job ID as first argument convention for all tasks (signal handler extraction)
- Scratch directory inside job directory for cleanup
- Let exceptions propagate to Huey for consistent error handling
- validate_environment() exits process on failure (fail-fast)
- recover_interrupted_jobs() scans for 'running' jobs at startup
- Quick test for CI, full test requires running consumer
- Tuple return pattern for validation: (mol, None) success, (None, error) failure
- Detect SDF vs MOL by $$$$ delimiter presence in content
- Flattened input fields in JobStatusResponse for simpler API
- RFC 7807 ProblemDetail schema for standardized error responses
- Health liveness is minimal (just return alive status)
- Readiness checks data directory writable and Huey importable
- Use JSONResponse with explicit headers for 202 Accepted responses
- OpenAPI JSON served at /api/v1/openapi.json (versioned)
- Health endpoints at root (no /api/v1 prefix)
- TestClient module-level client for test efficiency

### Pending Todos

None yet.

### Blockers/Concerns

- RDKit stderr capture doesn't work for C-level output (known limitation, fallback error messages used)

## Session Continuity

Last session: 2026-01-19T17:20:00Z
Stopped at: Completed 02-03-PLAN.md (Application Assembly) - Phase 2 Complete
Resume file: None
