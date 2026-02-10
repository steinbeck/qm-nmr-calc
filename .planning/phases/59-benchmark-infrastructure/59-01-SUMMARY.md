---
phase: 59-benchmark-infrastructure
plan: 01
subsystem: benchmark-infrastructure
tags: [nwchem, cosmo, solvents, benchmark-cli, input-generation]

requires:
  - phase-58: Expanded solvent support with 6 new solvents (CHCl3, DMSO, Methanol, Water, Acetone, Benzene)

provides:
  - 6 new solvents in NWChem input generation (pyridine, thf, toluene, dcm, acetonitrile, dmf)
  - COSMO name mapping layer for acetonitrile -> acetntrl
  - Benchmark CLI accepts 12 solvents via --solvents argument
  - Opt-in access to new solvents (not in default SOLVENTS list)

affects:
  - phase-60: Dispatch mechanics can now use new solvents
  - phase-61: Will add NWChem calculation execution for new solvents
  - phase-62: Analysis will process results from new solvents
  - phase-63: Scaling factor derivation will compute factors for new solvents

tech-stack:
  added: []
  patterns:
    - name-mapping-layer
    - validation-with-special-cases

key-files:
  created: []
  modified:
    - src/qm_nmr_calc/nwchem/input_gen.py
    - src/qm_nmr_calc/benchmark/__main__.py
    - tests/test_nwchem_input.py

decisions:
  - id: COSMO-NAME-MAPPING
    title: Map user-friendly solvent names to NWChem COSMO names
    rationale: NWChem COSMO recognizes "acetntrl" not "acetonitrile" - need mapping layer to accept user-friendly names while generating correct NWChem input
    alternatives: [hardcode "acetntrl" in all docs and CLI, reject "acetonitrile" entirely]
    tradeoffs: Adds one indirection layer but improves UX
  - id: OPT-IN-NEW-SOLVENTS
    title: New solvents are opt-in only via --solvents flag
    rationale: Scaling factors don't exist yet (Phase 63 derives them) - can't add to default SOLVENTS list until then
    alternatives: [add to defaults immediately with placeholder factors, block CLI access until factors exist]
    tradeoffs: Users must explicitly request new solvents but avoids shipping invalid defaults

metrics:
  duration: 17 minutes
  completed: 2026-02-10
---

# Phase 59 Plan 01: Add New Solvents with COSMO Mapping Summary

**One-liner:** Add 6 new solvents (pyridine, thf, toluene, dcm, acetonitrile, dmf) to NWChem input generation with acetonitrile->acetntrl COSMO name mapping.

## What Was Built

### COSMO Name Mapping Layer

Added a mapping layer between user-friendly solvent names and NWChem COSMO internal names:

- **COSMO_NAME_MAP**: Dictionary mapping `"acetonitrile"` to `"acetntrl"`
- **_get_cosmo_solvent_name()**: Helper function that applies mapping, passes through other solvents unchanged
- Applied in both `generate_optimization_input()` and `generate_shielding_input()`

NWChem COSMO parameter tables use abbreviated names for some solvents. Acetonitrile is the only one in our 13-solvent set that requires mapping.

### Extended Solvent Support

**Input Generation (input_gen.py):**
- Expanded SUPPORTED_SOLVENTS from 7 to 13 solvents
- Added: pyridine, thf, toluene, dcm, acetonitrile, dmf
- All new solvents validated and generate correct COSMO blocks
- Acetonitrile produces "solvent acetntrl" in COSMO block
- Others produce "solvent {name}" directly

**Benchmark CLI (__main__.py):**
- Expanded --solvents choices from 6 to 12 options
- Added: Pyridine, THF, Toluene, DCM, Acetonitrile, DMF (Title-case)
- Help text unchanged (still shows 6 default solvents)
- Users opt-in via `python -m qm_nmr_calc.benchmark run --solvents Pyridine`

**Runner defaults (runner.py):**
- SOLVENTS list UNCHANGED - still 6 entries (CHCl3, DMSO, Methanol, Water, Acetone, Benzene)
- New solvents not added to defaults because scaling factors don't exist yet
- Phase 63 will derive scaling factors, then defaults can be updated

### Test Coverage

**TestCOSMONameMapping (7 new tests):**
1. Acetonitrile maps to acetntrl in optimization input
2. Acetonitrile maps to acetntrl in shielding input
3. Case-insensitive acetonitrile handling
4. Direct-name solvents (pyridine, thf, toluene, dcm, dmf) use names unchanged
5. _get_cosmo_solvent_name() unit test
6. COSMO_NAME_MAP only contains acetonitrile

**Updated test_all_supported_solvents_accepted:**
- Now uses _get_cosmo_solvent_name() to check expected COSMO name
- Parametrized across all 12 non-vacuum solvents
- All pass including acetonitrile special case

**Test results:**
- 40 NWChem input tests pass (7 new)
- 72 focused tests pass (NWChem + benchmark modules)
- No regressions in existing functionality

## Decisions Made

### COSMO Name Mapping
**Decision:** Implement user-friendly name → NWChem COSMO name mapping layer

**Context:** NWChem COSMO parameter tables use abbreviated internal names. "acetonitrile" is recognized as "acetntrl". User-facing documentation and CLI use full names for clarity.

**Options considered:**
1. **Map at input generation** (chosen) - Accept "acetonitrile" everywhere, map to "acetntrl" when writing COSMO block
2. Hardcode "acetntrl" in docs/CLI - Users must know NWChem internal names
3. Reject "acetonitrile" - Only accept "acetntrl"

**Rationale:** Option 1 provides best UX. Users see consistent naming across CLI/docs/code. Mapping layer is simple (1 function, 1 dict) and testable. Future solvents with mapping needs can be added to COSMO_NAME_MAP.

**Impact:** One indirection layer (minimal overhead). Clear separation between external names (user-facing) and internal names (NWChem COSMO).

### Opt-In New Solvents
**Decision:** New solvents accessible via --solvents flag only, not in default SOLVENTS list

**Context:** runner.py SOLVENTS constant defines default solvents for benchmark runs. Scaling factors (slope/intercept for computed → experimental shift mapping) don't exist for new solvents until Phase 63 derives them from benchmark data. Calculations in Phase 60-62 will generate raw computed shifts but can't produce accurate final predictions without scaling factors.

**Options considered:**
1. **Keep out of defaults** (chosen) - Users opt-in via `--solvents Pyridine`
2. Add to defaults with placeholder factors - Would produce invalid predictions
3. Block CLI access entirely - Would prevent Phase 60-62 data collection

**Rationale:** Option 1 allows Phase 60-62 to collect benchmark data (needed to derive factors in Phase 63) while preventing users from getting invalid predictions if they run defaults. After Phase 63 completes, we can add new solvents to SOLVENTS with correct scaling factors.

**Impact:** Documentation must explain opt-in access. Phase 60-62 will use `--solvents Pyridine THF Toluene DCM Acetonitrile DMF --functionals B3LYP WP04` explicitly.

## Deviations from Plan

None - plan executed exactly as written.

## Next Phase Readiness

**Phase 60 (Dispatch Mechanics):** Ready
- Benchmark CLI accepts all 6 new solvents via --solvents
- Input generation produces valid NWChem COSMO input for all 6
- No changes needed in runner.py for Phase 60-62 (users opt-in via CLI)

**Phase 61 (NWChem Execution):** Ready
- COSMO blocks for new solvents validated in tests
- Acetonitrile mapping ensures NWChem will recognize solvent name
- Phase 61 will run actual NWChem calculations to verify COSMO parameters exist

**Phase 62 (Result Parsing):** Ready
- No changes needed - parsing logic is solvent-agnostic

**Phase 63 (Scaling Factor Derivation):** Ready
- Will process results from all 6 new solvents
- Will derive slope/intercept for each (functional, solvent) pair
- After Phase 63, new solvents can be added to runner.py SOLVENTS defaults

**Blockers/Concerns:** None

## Implementation Notes

### COSMO Name Mapping Architecture

```python
# User code / CLI
solvent = "Acetonitrile"  # Title-case from argparse

# Benchmark runner (runner.py line 234)
solvent_lower = solvent.lower()  # "acetonitrile"

# Input generation (input_gen.py)
solvent_name = _validate_solvent(solvent_lower)  # "acetonitrile" in SUPPORTED_SOLVENTS ✓
cosmo_name = _get_cosmo_solvent_name(solvent_name)  # "acetntrl" from COSMO_NAME_MAP

# NWChem input
f"solvent {cosmo_name}"  # "solvent acetntrl"
```

Chain works: CLI Title-case → runner.py lowercase → validation → mapping → NWChem COSMO block

### Why Only Acetonitrile Needs Mapping

NWChem COSMO parameter database (hardcoded in NWChem source) uses these names:
- pyridine → pyridine ✓
- thf → thf ✓
- toluene → toluene ✓
- dcm → dcm ✓
- acetonitrile → **acetntrl** (abbreviated)
- dmf → dmf ✓

Only acetonitrile diverges from common name. If future solvents need mapping, add to COSMO_NAME_MAP dict.

### Test Strategy

**Unit tests:** _get_cosmo_solvent_name() returns correct mapping/passthrough
**Integration tests:** Full input generation produces correct COSMO block text
**Parametrized tests:** All 12 non-vacuum solvents generate valid input
**Coverage:** Both optimization and shielding input paths tested

### Future Solvent Additions

To add another solvent (e.g., dichloroethane):

1. Check NWChem COSMO parameter name (may differ from common name)
2. Add to SUPPORTED_SOLVENTS in input_gen.py
3. If COSMO name differs, add to COSMO_NAME_MAP
4. Add to CLI --solvents choices in __main__.py
5. Run Phase 60-63 benchmark cycle to derive scaling factors
6. Add to runner.py SOLVENTS defaults once factors exist

## Files Changed

**src/qm_nmr_calc/nwchem/input_gen.py** (feat(59-01) ad28ea1):
- SUPPORTED_SOLVENTS: 7 → 13 solvents
- Added COSMO_NAME_MAP and _get_cosmo_solvent_name()
- Applied mapping in generate_optimization_input()
- Applied mapping in generate_shielding_input()

**src/qm_nmr_calc/benchmark/__main__.py** (feat(59-01) ad28ea1):
- CLI --solvents choices: 6 → 12 solvents (added Pyridine, THF, Toluene, DCM, Acetonitrile, DMF)

**tests/test_nwchem_input.py** (test(59-01) 20d3148):
- Added imports: COSMO_NAME_MAP, _get_cosmo_solvent_name
- Added TestCOSMONameMapping class with 7 tests
- Updated test_all_supported_solvents_accepted to use _get_cosmo_solvent_name()

## Verification

All success criteria met:

✅ Benchmark CLI accepts pyridine, thf, toluene, dcm, acetonitrile, dmf as valid --solvents values
✅ NWChem input for acetonitrile contains "solvent acetntrl" (COSMO mapping applied)
✅ NWChem input for pyridine, thf, toluene, dcm, dmf contains "solvent {name}" (direct names)
✅ Parametrized test covers all 12 non-vacuum solvents without failure
✅ All existing tests still pass (no regressions - 72 focused tests, 40 NWChem tests)
✅ runner.py SOLVENTS unchanged (new solvents are opt-in only)

## Task Summary

| Task | Description | Commit | Files Modified |
|------|-------------|--------|----------------|
| 1 | Add 6 new solvents to input_gen.py with COSMO name mapping | ad28ea1 | input_gen.py, __main__.py |
| 2 | Add tests for COSMO name mapping and new solvent acceptance | 20d3148 | test_nwchem_input.py |

**Total execution time:** 17 minutes
**Commits:** 2 (feat + test)
**Tests added:** 7 new tests
**Tests modified:** 1 parametrized test updated
**Tests passing:** 40 NWChem input tests, 72 focused tests
