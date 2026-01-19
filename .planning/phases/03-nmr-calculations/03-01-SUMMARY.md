---
phase: 03-nmr-calculations
plan: 01
subsystem: calculations
tags: [nmr, dft, presets, solvents, shielding]

# Dependency graph
requires:
  - phase: 02-input-and-api
    provides: FastAPI application structure, job models
provides:
  - PresetName enum and PRESETS configuration dict
  - Shielding-to-shift conversion with TMS reference
  - COSMO solvent validation for NWChem
affects: [03-02, 03-03, calculation pipeline]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - TypedDict for configuration data (vs BaseModel for API validation)
    - Linear regression for shielding-to-shift conversion

key-files:
  created:
    - src/qm_nmr_calc/presets.py
    - src/qm_nmr_calc/shifts.py
    - src/qm_nmr_calc/solvents.py
  modified: []

key-decisions:
  - "Production preset as default (reliability over speed)"
  - "TypedDict for preset config (not Pydantic - just config data)"
  - "TMS reference scaling factors from Pierens et al. for B3LYP/6-311+G(2d,p)"
  - "11 common NMR solvents supported initially"

patterns-established:
  - "Preset enum pattern: str + Enum base for JSON serialization"
  - "Shift conversion: m * shielding + b with nucleus-specific factors"

# Metrics
duration: 2min
completed: 2026-01-19
---

# Phase 3 Plan 1: Calculation Support Modules Summary

**Preset configurations (draft/production), shielding-to-shift conversion with TMS reference, and NWChem COSMO solvent validation**

## Performance

- **Duration:** 2 min
- **Started:** 2026-01-19T18:52:51Z
- **Completed:** 2026-01-19T18:55:09Z
- **Tasks:** 3
- **Files modified:** 3

## Accomplishments
- Created presets module with draft (6-31G*) and production (6-311+G(2d,p)) configurations
- Implemented shielding-to-shift conversion using linear regression (shift = m * shielding + b)
- Added solvent validation supporting 11 common NMR solvents with case-insensitive matching

## Task Commits

Each task was committed atomically:

1. **Task 1: Create presets module** - `05996f5` (feat)
2. **Task 2: Create shifts module** - `75553e0` (feat)
3. **Task 3: Create solvents module** - `a75e8c1` (feat)

## Files Created/Modified
- `src/qm_nmr_calc/presets.py` - PresetName enum, PRESETS dict, DEFAULT_PRESET constant
- `src/qm_nmr_calc/shifts.py` - SCALING_FACTORS, shielding_to_shift() with input validation
- `src/qm_nmr_calc/solvents.py` - SUPPORTED_SOLVENTS dict, validate_solvent(), get_supported_solvents()

## Decisions Made
- **Production as default:** Users get reliable results by default, opt-in to draft for speed
- **TypedDict over BaseModel:** Presets are pure configuration data, not API validation
- **TMS reference values:** Using published scaling factors (m=-1.0, b=31.8 for H, b=182.5 for C) for B3LYP/6-311+G(2d,p)
- **11 solvents initially:** Common NMR solvents (chcl3, dmso, methanol, acetone, benzene, h2o, acetntrl, dcm, thf, pyridine, toluene)

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- Presets module ready for integration into calculation tasks
- Shifts module ready to convert NMR shielding output to chemical shifts
- Solvents module ready for job input validation
- Next plan (03-02) can extend isicle_wrapper.py and tasks.py to use these modules

---
*Phase: 03-nmr-calculations*
*Completed: 2026-01-19*
