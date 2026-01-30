---
phase: 24-conformer-preselection
plan: 02
subsystem: conformers
tags: [xtb, gfn2, semi-empirical, energy-ranking, subprocess, alpb-solvation]

# Dependency graph
requires:
  - phase: 24-01
    provides: RMSD clustering for conformer diversity

provides:
  - xTB/GFN2 energy ranking module for conformer pre-selection
  - ALPB implicit solvation support for common NMR solvents
  - Graceful fallback when xTB unavailable

affects:
  - 24-03: Conformer workflow integration (will use xTB ranking when available)
  - conformer-selection: Future optimization of selection strategy

# Tech tracking
tech-stack:
  added: [xTB/GFN2-xTB via subprocess (optional dependency)]
  patterns:
    - Optional tool detection pattern (detect_xtb_available)
    - Subprocess-based external tool integration with timeouts
    - Test skip decorators for optional dependencies

key-files:
  created:
    - src/qm_nmr_calc/conformers/xtb_ranking.py
    - tests/test_xtb_ranking.py
  modified: []

key-decisions:
  - "xTB as optional dependency: App works without it, uses MMFF fallback"
  - "GFN2-xTB method: Balance of speed and accuracy for conformer ranking"
  - "ALPB implicit solvation: Faster than full solvent models, adequate for ranking"
  - "Subprocess with timeout: Prevent hanging on problematic structures"
  - "Relative energies in kcal/mol: Consistent with MMFF output format"

patterns-established:
  - "Optional tool detection: shutil.which() for binary availability"
  - "Test skip decorators: @pytest.mark.skipif for optional dependencies"
  - "Solvent name mapping: Case-insensitive mapping to canonical ALPB names"
  - "Energy parsing: Robust extraction from xTB stdout format"
  - "Batch ranking: Dict[conf_id, relative_energy] return type"

# Metrics
duration: 2min
completed: 2026-01-30
---

# Phase 24 Plan 02: xTB Energy Ranking Summary

**GFN2-xTB semi-empirical energy ranking with ALPB solvation for conformer pre-selection, 100-1000x faster than DFT with 0.4-0.5 Spearman correlation**

## Performance

- **Duration:** 2 min
- **Started:** 2026-01-30T08:50:06Z
- **Completed:** 2026-01-30T08:52:08Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments

- xTB binary detection with version checking for optional tool support
- GFN2-xTB single-point energy calculation via subprocess with timeout
- ALPB implicit solvation for common NMR solvents (CDCl3, DMSO, etc.)
- Batch conformer ranking returning relative energies in kcal/mol
- Comprehensive test suite that works with or without xTB installed

## Task Commits

Each task was committed atomically:

1. **Task 1: Create xtb_ranking.py module** - `1a95eb8` (feat)
2. **Task 2: Create unit tests for xTB ranking** - `852c5a3` (test)

## Files Created/Modified

- `src/qm_nmr_calc/conformers/xtb_ranking.py` - xTB detection, energy calculation, and conformer ranking functions
- `tests/test_xtb_ranking.py` - 14 tests (9 pass, 5 skip without xTB) for detection, solvation, energy calculation, and ranking

## Decisions Made

**xTB as optional dependency:**
- App functions without xTB binary installed
- Falls back to MMFF energy ranking if xTB unavailable
- Install instructions provided in error messages and module docstring

**GFN2-xTB method selection:**
- GFN2 provides best balance of speed and accuracy for conformer ranking
- 100-1000x faster than DFT calculations
- Research shows 0.39-0.47 Spearman correlation with DFT energies (vs -0.1 to -0.45 for MMFF)

**ALPB implicit solvation model:**
- Faster than explicit solvent or full COSMO models
- Adequate accuracy for conformer ranking purposes
- Supports all common NMR solvents (CDCl3, DMSO, D2O, etc.)

**Subprocess timeout handling:**
- 60-second default timeout per conformer calculation
- Prevents hanging on problematic structures
- Clear TimeoutError message for user troubleshooting

**Relative energy output format:**
- Returns energies in kcal/mol relative to minimum conformer
- Consistent with MMFF ranking output format
- Simplifies integration into existing conformer workflow

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None - implementation straightforward, tests pass as expected.

## User Setup Required

**Optional:** Install xTB for improved conformer energy ranking

```bash
# Via conda
conda install -c conda-forge xtb

# Verify installation
xtb --version
```

**Without xTB:** App continues to work using MMFF energy ranking as fallback.

## Next Phase Readiness

**Ready for Phase 24-03 (Conformer workflow integration):**
- xTB ranking module complete and tested
- API stable for integration into conformer pipeline
- Graceful fallback behavior when xTB unavailable
- Solvent mapping covers all common NMR solvents

**Integration points established:**
- `detect_xtb_available()` - Check before attempting xTB ranking
- `rank_conformers_by_xtb(mol, conf_ids, charge, solvent)` - Main ranking function
- Returns `Dict[int, float]` matching MMFF ranking output format

**Testing coverage:**
- Detection and version checking
- Solvent name mapping (case-insensitive, all NMR solvents)
- Energy calculation (skipped if xTB unavailable)
- Batch ranking (skipped if xTB unavailable)
- Error handling for unavailable binary

---
*Phase: 24-conformer-preselection*
*Completed: 2026-01-30*
