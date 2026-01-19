# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-01-19)

**Core value:** Reliable async NMR predictions with full control over calculation parameters -- submit a molecule, get back accurate 1H/13C shifts without babysitting long-running calculations.
**Current focus:** Phase 2 - Input and API

## Current Position

Phase: 2 of 6 (Input and API)
Plan: 2 of 3 in phase (complete)
Status: In progress
Last activity: 2026-01-19 -- Completed 02-02-PLAN.md (Routers)

Progress: [█████░░░░░] ~25%

## Performance Metrics

**Velocity:**
- Total plans completed: 5
- Average duration: 3.0 min
- Total execution time: 15 min

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 01-foundation | 3 | 9 min | 3.0 min |
| 02-input-and-api | 2 | 6 min | 3.0 min |

**Recent Trend:**
- Last 5 plans: 01-02 (3 min), 01-03 (2 min), 02-01 (4 min), 02-02 (2 min)
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

### Pending Todos

None yet.

### Blockers/Concerns

- RDKit stderr capture doesn't work for C-level output (known limitation, fallback error messages used)

## Session Continuity

Last session: 2026-01-19T17:16:56Z
Stopped at: Completed 02-02-PLAN.md (Routers)
Resume file: None
