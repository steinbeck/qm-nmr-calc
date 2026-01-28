---
phase: 16-crest-integration
plan: 01
subsystem: conformer-generation
tags: [crest, xtb, alpb, xyz-parsing, tdd]

# Dependency graph
requires:
  - phase: 13-rdkit-conformer-generation
    provides: Conformer generation pipeline structure
  - phase: 12-conformer-data-model
    provides: ConformerEnsemble and ConformerData models
provides:
  - CREST/xTB binary detection with caching
  - ALPB solvent model mapping for supported solvents
  - Multi-structure XYZ file parser for CREST ensemble output
  - CRESTConformer data structure
affects: [16-02-crest-runner, 16-03-crest-pipeline-dispatch]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "lru_cache for binary detection (result never changes during process)"
    - "Case-insensitive solvent mapping via lowercase normalization"
    - "Sequential 1-based conformer IDs (conf_001, conf_002, ...)"
    - "NamedTuple for simple data containers"

key-files:
  created:
    - src/qm_nmr_calc/conformers/crest_generator.py
    - tests/test_crest_generator.py
  modified: []

key-decisions:
  - "lru_cache(maxsize=1) for detect_crest_available - binary availability immutable during process lifetime"
  - "Both crest and xtb required - CREST depends on xTB, return False if either missing"
  - "ALPB solvent map hardcoded to chcl3/dmso - only supported solvents, extensible later"
  - "Case-insensitive solvent matching via lowercase normalization"
  - "parse_crest_ensemble extracts energy from first token of comment line (Hartree units)"
  - "Sequential 1-based conformer IDs for consistency with RDKit pipeline"

patterns-established:
  - "TDD with RED-GREEN-REFACTOR workflow: test → implementation → cleanup"
  - "lru_cache.cache_clear() in tests to isolate cache state between test cases"
  - "Patch at module level (qm_nmr_calc.conformers.crest_generator.shutil.which) for cached functions"

# Metrics
duration: 12min
completed: 2026-01-28
---

# Phase 16 Plan 01: CREST Foundation Utilities Summary

**CREST/xTB binary detection, ALPB solvent mapping (CHCl3/DMSO), and multi-structure XYZ parser with energy extraction**

## Performance

- **Duration:** 12 min
- **Started:** 2026-01-28T10:29:00Z
- **Completed:** 2026-01-28T10:41:40Z
- **Tasks:** 1 (TDD feature with 3 phases)
- **Files modified:** 2 created
- **Tests:** 15 new, 244 total passing

## Accomplishments

- `detect_crest_available()` with lru_cache for efficient binary detection (requires both crest and xtb)
- `get_alpb_solvent()` maps job solvent names to ALPB identifiers (CHCl3→chcl3, DMSO→dmso, vacuum→None)
- `parse_crest_ensemble()` parses concatenated multi-structure XYZ files with energy extraction
- `CRESTConformer` NamedTuple for conformer data (id, energy_hartree, xyz_block)
- Full test coverage with 15 new tests covering all behaviors

## Task Commits

Each TDD phase was committed atomically:

1. **RED Phase: Failing tests** - `b5ff287` (test)
   - Created tests/test_crest_generator.py with 15 test functions
   - Tests for binary detection (4 cases), solvent mapping (6 cases), XYZ parsing (3 cases)
   - Tests fail due to missing module

2. **GREEN Phase: Implementation** - `362bfe7` (feat)
   - Created src/qm_nmr_calc/conformers/crest_generator.py
   - Implemented all three utilities plus CRESTConformer and ALPB_SOLVENT_MAP
   - Fixed test cache isolation (patch at module level, clear cache between tests)
   - All 244 tests passing (229 existing + 15 new)

3. **REFACTOR Phase:** No changes needed - code clean and follows project patterns

## Files Created/Modified

- `src/qm_nmr_calc/conformers/crest_generator.py` - CREST utilities module with detection, mapping, and parsing functions
- `tests/test_crest_generator.py` - Comprehensive test suite with 15 test functions across 5 test classes

## Decisions Made

1. **lru_cache(maxsize=1) for binary detection**: Binary availability is immutable during process lifetime, cache eliminates redundant PATH lookups
2. **Require both binaries**: CREST depends on xTB, so `detect_crest_available()` returns False unless both present
3. **Hardcoded ALPB solvent map**: Only chcl3 and dmso currently supported, but map is extensible for future additions
4. **Case-insensitive matching**: Normalize input to lowercase before lookup to handle CHCl3, chcl3, CHCL3 variants
5. **Energy from comment line first token**: CREST XYZ format puts Hartree energy as first token on comment line
6. **Sequential 1-based IDs**: conf_001, conf_002, ... matches RDKit pipeline naming convention (established in 13-03)

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Fixed test cache isolation for lru_cache**
- **Found during:** GREEN phase (test execution)
- **Issue:** First test cached result, subsequent tests with different mocks failed because cache returned stale value
- **Fix:**
  - Added `detect_crest_available.cache_clear()` before each test
  - Changed patch target from `shutil.which` to `qm_nmr_calc.conformers.crest_generator.shutil.which` (module-level)
- **Files modified:** tests/test_crest_generator.py
- **Verification:** All 4 binary detection tests pass with different mock configurations
- **Committed in:** 362bfe7 (GREEN phase commit)

---

**Total deviations:** 1 auto-fixed (blocking test issue)
**Impact on plan:** Essential fix for test isolation with cached functions. No scope change, proper testing practice.

## Issues Encountered

None - TDD workflow proceeded smoothly through RED-GREEN-REFACTOR.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

**Ready for Plan 02 (CREST Runner):**
- Binary detection available for runtime capability checks
- ALPB solvent mapping ready for --alpb flag construction
- XYZ parser ready for crest_conformers.xyz output processing
- CRESTConformer structure defined for data passing

**Foundation complete:**
- All three pure utility functions fully tested
- Clear inputs/outputs for CREST runner integration
- Pattern established: check availability → map solvent → run CREST → parse ensemble

**Next:**
- Plan 02 will build `run_crest_conformer_search()` using these utilities
- Plan 03 will integrate CREST into pipeline dispatch logic with RDKit fallback

---
*Phase: 16-crest-integration*
*Completed: 2026-01-28*
