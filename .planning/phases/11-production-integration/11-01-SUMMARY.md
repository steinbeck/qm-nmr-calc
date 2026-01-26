---
phase: 11-production-integration
plan: 01
subsystem: nmr-calculations
tags: [delta50, scaling-factors, regression, nwchem, shifts]

# Dependency graph
requires:
  - phase: 10-scaling-factors
    provides: "DELTA50-derived regression factors in data/benchmark/delta50/scaling_factors.json"
provides:
  - "Production calculations use DELTA50 regression factors instead of CHESHIRE TMS references"
  - "Scaling factors bundled in package via importlib.resources"
  - "Factor lookup by composite key: functional/basis_set/nucleus/solvent"
  - "ISiCLE dependency completely removed"
affects: [production-calculations, api-responses, benchmarking]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Lazy loading with @cache decorator for JSON data"
    - "importlib.resources for package data access"
    - "Explicit ValueError on missing scaling factors"

key-files:
  created:
    - "src/qm_nmr_calc/data/scaling_factors.json"
  modified:
    - "src/qm_nmr_calc/shifts.py"
    - "src/qm_nmr_calc/tasks.py"
    - "src/qm_nmr_calc/benchmark/runner.py"
    - "pyproject.toml"

key-decisions:
  - "Replace preset-based scaling with direct functional/basis_set/solvent parameters"
  - "Fail explicitly with ValueError when no scaling factor exists (no fallback)"
  - "Use @cache decorator for one-time JSON loading at first call"
  - "Add jinja2, scipy, statsmodels as direct dependencies after ISiCLE removal"

patterns-established:
  - "Package data bundled via hatchling force-include configuration"
  - "Composite keys for scaling factor lookup: functional/basis/nucleus/solvent"

# Metrics
duration: 5min
completed: 2026-01-23
---

# Phase 11 Plan 01: Production Integration Summary

**DELTA50 regression factors integrated into production via lazy-loaded JSON with explicit functional/basis_set/solvent parameters; ISiCLE dependency removed**

## Performance

- **Duration:** 5 min
- **Started:** 2026-01-23T15:49:29Z
- **Completed:** 2026-01-23T15:54:41Z
- **Tasks:** 3
- **Files modified:** 5

## Accomplishments

- Bundled DELTA50 scaling factors in package data directory with hatchling force-include
- Rewrote shifts.py to load factors from JSON using importlib.resources with @cache decorator
- Updated tasks.py and benchmark/runner.py to call shielding_to_shift with new signature
- Removed ISiCLE from dependencies entirely (replaced with direct jinja2, scipy, statsmodels)

## Task Commits

Each task was committed atomically:

1. **Task 1: Copy scaling factors and configure hatchling build** - `426d7a6` (chore)
2. **Task 2: Rewrite shifts.py to use DELTA50 regression factors** - `3003218` (refactor)
3. **Task 3: Update tasks.py and runner.py to use new shielding_to_shift signature** - `8ad1716` (feat)

**Dependency fix:** `ca0838d` (fix: add jinja2, scipy, statsmodels as direct dependencies)

## Files Created/Modified

- `src/qm_nmr_calc/data/scaling_factors.json` - DELTA50 regression factors (4 entries: B3LYP/CHCl3+DMSO/1H+13C)
- `src/qm_nmr_calc/shifts.py` - Load factors from JSON, lookup by composite key, apply regression formula
- `src/qm_nmr_calc/tasks.py` - Pass functional/basis_set/solvent from preset config
- `src/qm_nmr_calc/benchmark/runner.py` - Pass functional/basis_set/solvent from task params
- `pyproject.toml` - Removed isicle, added force-include for data/, added jinja2/scipy/statsmodels

## Decisions Made

1. **Composite key format:** Used `functional/basis_set/nucleus/solvent` for factor lookup (e.g., "B3LYP/6-311+G(2d,p)/1H/CHCl3")
2. **No fallback behavior:** Raise explicit ValueError if requested combination not supported (fail-fast)
3. **Lazy loading with cache:** Use @cache decorator so JSON loads once on first call, not at import time
4. **Direct dependencies:** Added jinja2, scipy, statsmodels as direct deps after ISiCLE removal exposed missing transitives

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Added jinja2 as direct dependency**
- **Found during:** Task 1 verification
- **Issue:** ISiCLE removal uninstalled jinja2, breaking FastAPI's Jinja2Templates
- **Fix:** Added `jinja2>=3.1.0` to pyproject.toml dependencies
- **Files modified:** pyproject.toml
- **Verification:** uv sync succeeded, templates import without error
- **Committed in:** `ca0838d` (fix commit)

**2. [Rule 3 - Blocking] Added scipy and statsmodels as direct dependencies**
- **Found during:** Task 1 verification
- **Issue:** ISiCLE removal uninstalled scipy and statsmodels, breaking benchmark/analysis.py imports
- **Fix:** Added `scipy>=1.17.0` and `statsmodels>=0.14.0` to pyproject.toml
- **Files modified:** pyproject.toml
- **Verification:** uv sync succeeded, analysis.py imports without error
- **Committed in:** `ca0838d` (fix commit)

---

**Total deviations:** 2 auto-fixed (2 blocking issues)
**Impact on plan:** Both fixes essential for unblocking imports after ISiCLE removal. No scope creep - restoring functionality that ISiCLE provided transitively.

## Issues Encountered

None beyond the expected dependency cleanup after ISiCLE removal.

## Next Phase Readiness

- Production calculations now use DELTA50-derived factors with explicit R²/MAE metadata available
- Factor lookup fails explicitly for unsupported solvent/functional combinations
- ISiCLE completely removed from project
- Ready for API/UI updates to display factor metadata and accuracy estimates
- WP04 functional deferred (not supported by NWChem without custom compilation)

---
*Phase: 11-production-integration*
*Completed: 2026-01-23*
