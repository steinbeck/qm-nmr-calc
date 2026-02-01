---
phase: 32-core-nmredata-module
plan: 01
subsystem: nmr-export
tags: [nmredata, sdf, rdkit, chemical-shifts, data-export]

# Dependency graph
requires:
  - phase: 31-documentation-polish
    provides: Complete API documentation and user guides
provides:
  - NMReData SDF generation module with all required tags
  - Solvent name mapping (NWChem → NMReData format)
  - 1-indexed atom assignment formatting
  - Ensemble metadata support (Boltzmann averaging provenance)
affects: [33-api-ui-integration, nmr-export, data-interoperability]

# Tech tracking
tech-stack:
  added: []  # No new dependencies (uses existing RDKit)
  patterns:
    - "NMReData tag formatting with proper SDF structure"
    - "1-based atom indexing preserved from NWChem output"
    - "Solvent mapping table for deuterated NMR solvents"

key-files:
  created:
    - src/qm_nmr_calc/nmredata.py (270 lines)
    - tests/test_nmredata.py (434 lines)
  modified: []

key-decisions:
  - "Use lowercase atom labels (h5, c3) for cleaner NMReData format"
  - "Preserve 1-based atom indices from NWChem (no conversion needed)"
  - "Map NWChem solvents to deuterated forms: chcl3→CDCl3, dmso→(CD3)2SO"
  - "Format shifts with 4 decimal precision per NMReData spec"
  - "Generate SDF on-demand (no pre-generation) following existing /geometry.sdf pattern"

patterns-established:
  - "format_sdf_tag() helper for consistent SDF tag formatting"
  - "NMREDATA_SEP constant for ', ' separator throughout"
  - "Comprehensive test coverage with RDKit round-trip parsing"

# Metrics
duration: 10min
completed: 2026-02-01
---

# Phase 32 Plan 01: Core NMReData Generation Module Summary

**Complete NMReData SDF export with 8 required tags, 1-indexed atom assignments, and Boltzmann ensemble metadata**

## Performance

- **Duration:** 10 min
- **Started:** 2026-02-01T13:15:44Z
- **Completed:** 2026-02-01T13:25:12Z
- **Tasks:** 3
- **Files modified:** 2 (created)

## Accomplishments

- Created `nmredata.py` module with `generate_nmredata_sdf()` function producing NMReData-compliant SDF files
- Implemented all 8 required NMReData tags (VERSION, LEVEL, SOLVENT, TEMPERATURE, ASSIGNMENT, FORMULA, SMILES, ID)
- Verified 1-indexed atom assignments match NWChem convention (no off-by-one errors)
- Added ensemble mode support with Boltzmann averaging metadata in provenance tag
- Comprehensive test suite with 36 tests including RDKit round-trip parsing validation

## Task Commits

Each task was committed atomically:

1. **Task 1: Create nmredata.py module with tag infrastructure and solvent mapping** - `1712886` (feat)
   - Module constants: NMREDATA_VERSION, NMREDATA_LEVEL, NMREDATA_SEP
   - Solvent mapping for chcl3, dmso, vacuum
   - Atom label formatting (lowercase h/c)
   - 14 unit tests

2. **Task 2: Implement ASSIGNMENT tag formatting with 1-indexed atoms** - `528e8a4` (feat)
   - format_assignment_tag() with multiline format
   - Preserves 1-based indices from NWChem
   - 4 decimal precision, ", " separator
   - 8 unit tests verifying indexing and format

3. **Task 3: Implement generate_nmredata_sdf() with all tags and ensemble support** - `0ed7201` (feat)
   - Main generation function with XYZ parsing
   - RDKit Mol creation with 3D coordinates
   - All 8 NMReData tags injection
   - Ensemble provenance metadata
   - 14 comprehensive unit tests

## Files Created/Modified

**Created:**
- `src/qm_nmr_calc/nmredata.py` (270 lines) - NMReData SDF generation module
  - `NMREDATA_VERSION = "1.1"` - Format version constant
  - `NMREDATA_LEVEL = "0"` - Predicted data level
  - `NMREDATA_SEP = ", "` - Separator constant
  - `map_solvent_to_nmredata(solvent)` - Solvent name mapping
  - `format_atom_label(atom, index)` - Lowercase atom labels
  - `format_assignment_tag(h1_shifts, c13_shifts)` - ASSIGNMENT tag formatting
  - `format_sdf_tag(name, value)` - SDF tag helper
  - `generate_nmredata_sdf(...)` - Main generation function

- `tests/test_nmredata.py` (434 lines) - Comprehensive unit tests
  - TestConstants (3 tests) - Verify format constants
  - TestSolventMapping (6 tests) - Solvent name mapping
  - TestAtomLabel (5 tests) - Atom label formatting
  - TestAssignmentTag (8 tests) - ASSIGNMENT tag with 1-indexed atoms
  - TestSDFTagFormat (2 tests) - SDF tag formatting helper
  - TestGenerateNMReDataSDF (12 tests) - Full SDF generation with RDKit round-trip

**Modified:** None

## Decisions Made

1. **Lowercase atom labels:** Use `h5`, `c3` format instead of `H5`, `C3` for cleaner, more readable NMReData output
2. **1-based indexing preserved:** NWChem outputs are already 1-based, so no conversion needed (avoiding potential off-by-one errors)
3. **Solvent mapping table:** Explicit mapping from NWChem names to NMReData deuterated forms:
   - `chcl3` → `CDCl3`
   - `dmso` → `(CD3)2SO` (NMReData convention)
   - `vacuum` → `vacuum` (gas phase)
4. **4 decimal precision:** Format shifts as `{shift:.4f}` per NMReData specification
5. **Separator constant:** Use `NMREDATA_SEP = ", "` throughout for consistency and easy updates
6. **Ensemble metadata in provenance:** Include Boltzmann averaging info in NMREDATA_ID tag when `is_ensemble=True`

## Deviations from Plan

**None - plan executed exactly as written**

All functionality implemented as specified in plan. No auto-fixes needed, no architectural changes required.

## Issues Encountered

**RDKit SDMolSupplier hydrogen handling:** Initial round-trip test failed because RDKit's default `SDMolSupplier()` removes explicit hydrogens. Fixed by adding `removeHs=False` parameter to preserve full atom count during parsing validation.

## Test Coverage

**36 tests total, all passing:**
- 3 constant verification tests
- 6 solvent mapping tests (all 3 solvents + error cases)
- 5 atom label formatting tests
- 8 ASSIGNMENT tag tests (format, precision, indexing, separator)
- 2 SDF tag format tests
- 12 full SDF generation tests (all tags, RDKit parsing, ensemble mode, error handling)

**Coverage highlights:**
- 1-indexed atom verification (no off-by-one errors)
- RDKit round-trip parsing (export → parse → verify atom count)
- Ensemble mode provenance metadata
- Error handling (invalid SMILES, atom count mismatch)
- Tag format compliance (all 8 required NMReData tags)

## Next Phase Readiness

**Ready for Phase 33 (API and UI Integration):**
- ✅ `generate_nmredata_sdf()` function ready for API endpoint consumption
- ✅ Handles both single-conformer and ensemble jobs
- ✅ Comprehensive test coverage ensures reliability
- ✅ No external dependencies added (uses existing RDKit)

**API integration can follow existing `/geometry.sdf` pattern:**
- Load job status
- Read optimized.xyz geometry
- Extract h1_shifts and c13_shifts from nmr_results
- Call `generate_nmredata_sdf()`
- Return as `chemical/x-mdl-sdfile` response

**No blockers or concerns.**

---
*Phase: 32-core-nmredata-module*
*Completed: 2026-02-01*
