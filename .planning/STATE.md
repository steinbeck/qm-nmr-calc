# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-01-19)

**Core value:** Reliable async NMR predictions with full control over calculation parameters -- submit a molecule, get back accurate 1H/13C shifts without babysitting long-running calculations.
**Current focus:** Phase 5 - Visualization (Plan 1 complete, Plan 2 remaining)

## Current Position

Phase: 5 of 6 (Visualization)
Plan: 1 of 2 in phase (05-01-PLAN.md complete)
Status: In progress
Last activity: 2026-01-20 -- Completed 05-01-PLAN.md (Visualization Module)

Progress: [█████████████░░] ~71%

## Performance Metrics

**Velocity:**
- Total plans completed: 12
- Average duration: 3.0 min
- Total execution time: 37 min

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 01-foundation | 3 | 9 min | 3.0 min |
| 02-input-and-api | 3 | 8 min | 2.7 min |
| 03-nmr-calculations | 3 | 11 min | 3.7 min |
| 04-results-delivery | 2 | 7 min | 3.5 min |
| 05-visualization | 1 | 2 min | 2.0 min |

**Recent Trend:**
- Last 5 plans: 03-03 (5 min), 04-01 (3 min), 04-02 (4 min), 05-01 (2 min)
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
- Production preset as default (reliability over speed)
- TypedDict for preset config (not Pydantic - just config data)
- TMS reference scaling factors from Pierens et al. for B3LYP/6-311+G(2d,p)
- 11 common NMR solvents supported initially
- Solvent is required field (no default) - users must explicitly choose
- AtomShift stores both shielding and shift; API returns only shift
- NMRResults includes calculation metadata (functional, basis_set, solvent)
- Two-step DFT: geometry optimization then NMR shielding
- run_nmr_task updates job status with NMR results
- Validate solvent at API level before job creation
- Return 409 Conflict for incomplete jobs (vs 404 for missing)
- SDF generated on-the-fly from SMILES + XYZ coordinates
- Output ZIP includes only .out and .nw files from scratch directory
- aiosmtplib for async SMTP, email-validator for Pydantic EmailStr
- Best-effort email delivery (logs errors, never fails jobs)
- Environment variables for all SMTP config (no hardcoded credentials)
- Async/sync wrapper pattern for Huey signal handlers
- Agg backend set before pyplot import for headless rendering
- NWChem 1-based to RDKit 0-based index conversion for atomNote
- 2400x1800 PNG for 300 DPI equivalent at 8x6 inches
- atomNote set before PrepareMolForDrawing()
- plt.close(fig) after savefig to prevent memory leaks

### Pending Todos

None yet.

### Blockers/Concerns

- RDKit stderr capture doesn't work for C-level output (known limitation, fallback error messages used)

## Session Continuity

Last session: 2026-01-20T10:29:17Z
Stopped at: Completed 05-01-PLAN.md (Visualization Module)
Resume file: None
