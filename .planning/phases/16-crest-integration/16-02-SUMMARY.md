---
phase: 16-crest-integration
plan: 02
type: tdd
subsystem: conformers
completed: 2026-01-28
duration: 17min
tags: [crest, subprocess, ensemble-generation, energy-conversion, tdd]

requires:
  - "16-01: CREST binary detection and XYZ parsing"
  - "13-03: ConformerEnsemble model and sequential conformer IDs"
  - "13-02: Energy window filtering"
  - "14-01: HARTREE_TO_KCAL constant"

provides:
  - "run_crest() subprocess runner with timeout handling"
  - "generate_crest_ensemble() complete SMILES-to-ConformerEnsemble pipeline"
  - "Energy conversion from CREST Hartree to relative kcal/mol"
  - "XYZ file writing per conformer"

affects:
  - "16-03: Will integrate into conformer_mode selection logic"
  - "17: API will expose CREST as conformer_method option"

tech-stack:
  added:
    - "subprocess module for CREST execution"
    - "OMP_STACKSIZE and GFORTRAN_UNBUFFERED_ALL env vars"
  patterns:
    - "Subprocess timeout handling for long-running CREST jobs"
    - "Error chaining with helpful user-facing messages"
    - "Pre-DFT energy filtering in kcal/mol for consistency with RDKit path"

key-files:
  created: []
  modified:
    - "src/qm_nmr_calc/conformers/crest_generator.py: Added run_crest() and generate_crest_ensemble()"
    - "tests/test_crest_generator.py: Added 12 new tests (27 total)"

decisions:
  - id: "subprocess-timeout-default"
    choice: "Default 7200 seconds (2 hours) for CREST timeout"
    rationale: "Macrocycles can take hours; 2hr balances completeness vs responsiveness"
    alternatives: ["3600s (1hr): too short for macrocycles", "no timeout: job can hang forever"]

  - id: "crest-timeout-message"
    choice: "Suggest RDKit fallback in timeout error message"
    rationale: "User-facing guidance for when CREST is too slow"
    alternatives: ["Generic error: less helpful", "Automatic fallback: loses explicit control"]

  - id: "energy-units-crest"
    choice: "Store relative kcal/mol in ConformerData, not absolute Hartree"
    rationale: "Consistency with RDKit path; pre-DFT filtering uses same scale"
    alternatives: ["Hartree: mismatch with RDKit path", "Convert on-demand: extra complexity"]

  - id: "skip-rmsd-dedup-crest"
    choice: "No RMSD deduplication after CREST"
    rationale: "CREST handles deduplication internally during search"
    alternatives: ["Duplicate dedup: unnecessary computational cost", "Trust CREST completely"]

  - id: "directory-creation-safety"
    choice: "mkdir(parents=True, exist_ok=True) before XYZ write"
    rationale: "Defensive programming; handles edge cases where create_conformer_directories mock doesn't create dirs"
    alternatives: ["Assume dirs exist: brittle in tests", "Check exists before mkdir: extra syscall"]
---

# Phase 16 Plan 02: CREST Runner and Ensemble Builder Summary

**One-liner:** Subprocess execution with timeout handling and complete SMILES-to-ConformerEnsemble pipeline for CREST conformer generation

## What Was Built

Implemented two critical functions for CREST integration:

1. **run_crest()**: Subprocess runner
   - Builds CREST command with GFN2-xTB, ALPB solvent, charge, energy window, thread count
   - Sets OpenMP environment variables (OMP_STACKSIZE=2G, GFORTRAN_UNBUFFERED_ALL=1)
   - Timeout handling with user-friendly error suggesting RDKit fallback
   - Non-zero exit code handling with stderr capture
   - Output file verification (crest_conformers.xyz must exist)

2. **generate_crest_ensemble()**: Complete pipeline
   - Generate initial 3D geometry from SMILES using RDKit
   - Write input XYZ for CREST
   - Run CREST conformational search
   - Parse multi-structure XYZ output
   - Convert energies from Hartree to relative kcal/mol
   - Apply energy window filter (6.0 kcal/mol default)
   - Skip RMSD deduplication (CREST handles internally)
   - Write individual XYZ files per conformer
   - Build ConformerEnsemble with method="crest"

## TDD Commits

| Phase | Commit | Description |
|-------|--------|-------------|
| RED   | 60cb696 | Added 12 failing tests for run_crest() and generate_crest_ensemble() |
| GREEN | 550e495 | Implemented both functions, all tests passing |

## Test Coverage

27 tests total in test_crest_generator.py:
- 15 from Plan 01 (binary detection, ALPB mapping, XYZ parsing)
- 12 new for Plan 02:
  - 7 tests for run_crest() (success, timeout, error codes, command args, env vars, threads)
  - 5 tests for generate_crest_ensemble() (success, energy conversion, filtering, XYZ files, ConformerData fields)

## Technical Decisions

### Energy Unit Consistency
**Challenge:** CREST outputs absolute Hartree energies; RDKit path uses relative kcal/mol.

**Solution:** Convert to relative kcal/mol during generate_crest_ensemble():
- Subtract minimum Hartree energy
- Multiply by HARTREE_TO_KCAL (627.5095)
- Store in ConformerData.energy with energy_unit="kcal_mol"

**Benefit:** Pre-DFT filtering uses same scale as RDKit path. Downstream DFT processing sees identical data structure.

### Timeout Handling
**Challenge:** CREST can hang on macrocycles (hours or never complete).

**Solution:**
- subprocess.run() with timeout parameter (default: 7200s = 2 hours)
- Catch subprocess.TimeoutExpired
- Raise RuntimeError with message: "CREST timeout... Consider using RDKit conformer generation mode instead."

**Benefit:** Jobs don't hang forever; user gets actionable guidance.

### Environment Variables
**Challenge:** CREST/xTB need specific OpenMP and Fortran settings for stability.

**Solution:**
- Copy os.environ
- Set OMP_STACKSIZE=2G (prevents stack overflow)
- Set GFORTRAN_UNBUFFERED_ALL=1 (immediate output for logging)
- Pass env dict to subprocess.run()

**Benefit:** Stable CREST execution; output captured correctly.

### No RMSD Deduplication
**Challenge:** Should we run RMSD dedup after CREST like RDKit path?

**Solution:** Skip it. CREST performs internal deduplication during conformational search.

**Benefit:** Avoid redundant computation; trust CREST's internal algorithm.

## Deviations from Plan

None - plan executed exactly as written.

## Integration Points

### Upstream Dependencies
- parse_crest_ensemble() (Plan 01) for XYZ parsing
- get_alpb_solvent() (Plan 01) for solvent mapping (caller validates)
- filter_by_energy_window() (Plan 13-02) for pre-DFT filtering
- HARTREE_TO_KCAL (Plan 14-01) for unit conversion
- ConformerData, ConformerEnsemble (Plan 12-01) for data model
- smiles_to_xyz() (Phase 3) for initial geometry generation

### Downstream Usage
Plan 16-03 will:
- Call generate_crest_ensemble() when conformer_method="crest"
- Handle CREST unavailable fallback to RDKit
- Integrate into unified generate_ensemble() function

Phase 17 API will:
- Expose conformer_method="crest" option in job submission
- Document CREST timeout behavior
- Provide status updates during long CREST runs

## Files Modified

**src/qm_nmr_calc/conformers/crest_generator.py** (225 lines added)
- Added imports: os, subprocess, HARTREE_TO_KCAL, filter_by_energy_window, ConformerData, ConformerEnsemble
- Implemented run_crest() with subprocess.run(), timeout handling, env vars
- Implemented generate_crest_ensemble() with full pipeline orchestration
- Directory creation safety: mkdir(parents=True, exist_ok=True)

**tests/test_crest_generator.py** (361 lines added)
- TestRunCrest: 7 tests covering success, timeout, errors, command building, env vars
- TestGenerateCrestEnsemble: 5 tests covering pipeline, energy conversion, filtering, files, data fields
- All tests use mocked subprocess.run() - no real CREST execution

## Next Phase Readiness

Phase 16 Plan 03 can proceed:
- run_crest() and generate_crest_ensemble() tested and working
- Error handling provides clear user feedback
- Energy units match RDKit path for downstream consistency

Ready for integration into conformer_mode selection logic.

## Lessons Learned

1. **Subprocess mocking:** Mock at module level (qm_nmr_calc.conformers.crest_generator.subprocess.run) not subprocess directly
2. **Test filesystem isolation:** mkdir(parents=True, exist_ok=True) makes code resilient to mock variations
3. **Timeout messaging:** Include actionable alternative (RDKit mode) in error messages for better UX
4. **Energy unit consistency:** Decide early whether to store absolute or relative energies; changes propagate through entire pipeline

## Metrics

- **Duration:** 17 minutes
- **LOC Added:** 586 (225 implementation + 361 tests)
- **Tests Added:** 12 (27 total in file)
- **Functions Added:** 2 (run_crest, generate_crest_ensemble)
- **Dependencies Added:** 0 (used existing stdlib subprocess)
- **Test Coverage:** 100% of new functions via TDD
