# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-01-21)

**Core value:** Reliable async NMR predictions with full control over calculation parameters -- submit a molecule, get back accurate 1H/13C shifts without babysitting long-running calculations.
**Current focus:** Milestone v1.1 — Accurate Chemical Shifts

## Current Position

Phase: 8.1 of 11 (DELTA50 Data Viewer)
Plan: 1 of 1 in phase (COMPLETE)
Status: Phase complete
Last activity: 2026-01-21 — Completed 08.1-01-PLAN.md

Progress: [████████████████████████████████░░░░░░░░░] 59% (24/41 plans)

## Performance Metrics

**Velocity:**
- Total plans completed: 24
- Average duration: 3.3 min
- Total execution time: 80 min

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 01-foundation | 3 | 9 min | 3.0 min |
| 02-input-and-api | 3 | 8 min | 2.7 min |
| 03-nmr-calculations | 3 | 11 min | 3.7 min |
| 04-results-delivery | 2 | 7 min | 3.5 min |
| 05-visualization | 2 | 5 min | 2.5 min |
| 06-web-ui | 3 | 6 min | 2.0 min |
| 07-nwchem-integration | 4 | 15 min | 3.75 min |
| 08-delta50-setup | 2 | 12 min | 6.0 min |
| 08.1-delta50-viewer | 1 | 4 min | 4.0 min |

**Recent Trend:**
- Last 5 plans: 07-04 (6 min), 08-01 (9 min), 08-02 (3 min), 08.1-01 (4 min)
- Trend: Stable

*Updated after each plan completion*

## Accumulated Context

### Decisions

Decisions are logged in PROJECT.md Key Decisions table.
Recent decisions affecting current work:

- Use hatchling build system for src layout (uv compatible)
- 12-character hex job IDs from uuid4 for URL-safe identifiers
- orjson with OPT_INDENT_2 for human-readable status.json files
- ISiCLE removed, direct NWChem integration via nwchem module
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
- CHCl3 and DMSO solvents supported (those with COSMO dielectric values)
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
- Visualizations generated before job status update (ready when complete)
- Helper function _get_visualization() for common endpoint logic
- Pico CSS blue theme via CDN for minimal/clean scientific aesthetic
- Static files mounted at /static before routers to avoid route conflicts
- Web router at root (no /api prefix) for browser-friendly URLs
- Native HTML dialog over JavaScript modal library (zero dependencies)
- 303 See Other redirect for incomplete jobs to status page
- Re-use submit.html template for error display (404, 500)
- Quoted basis set names in NWChem input to handle special characters
- Case-insensitive solvent validation with lowercase lookup
- COSMO_DIELECTRIC dict pattern for solvent dielectric constants
- Flexible regex patterns with fallbacks for NWChem version variations
- Shielding data format: {index: [], atom: [], shielding: []} for shifts.py compatibility
- Descriptive RuntimeError messages with expected NWChem section headers
- ETKDGv3 as default conformer generation method with deterministic seeding (0xF00D)
- XYZ bond determination requires explicit charge parameter (default=0)
- run_calculation() as single entry point for NMR calculations
- COSMO solvation applied to BOTH geometry optimization and NMR shielding
- Keep isicle_version field in models with "N/A" for backwards compatibility
- DELTA50: compound_XX.xyz naming convention for benchmark molecules
- DELTA50: Both unique shift lists and atom-index mappings in JSON
- DELTA50: Separate BENCHMARK_PRESETS with WP04 for 1H-optimized calculations
- DELTA50: Resume via shifts.json marker file existence check
- 3Dmol.js from CDN for minimal build complexity
- Full XYZ format reconstruction for 3Dmol.js compatibility

### Roadmap Evolution

- Phase 8.1 inserted after Phase 8: DELTA50 Data Viewer - verify extracted structures and experimental shifts before benchmark execution (user request)

### Pending Todos

None yet.

### Blockers/Concerns

- RDKit stderr capture doesn't work for C-level output (known limitation, fallback error messages used)

## Session Continuity

Last session: 2026-01-21T20:52:02Z
Stopped at: Completed Phase 8.1 (DELTA50 Data Viewer)
Resume file: None
