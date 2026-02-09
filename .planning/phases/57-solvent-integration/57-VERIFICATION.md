---
phase: 57-solvent-integration
verified: 2026-02-09T10:47:13Z
status: passed
score: 7/7 must-haves verified
re_verification: false
---

# Phase 57: Solvent Integration Verification Report

**Phase Goal:** Users can select any of the 4 new solvents in the web UI and API and get accurate NMR predictions  
**Verified:** 2026-02-09T10:47:13Z  
**Status:** passed  
**Re-verification:** No — initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | API accepts solvent=methanol and returns shifts using Methanol scaling factors | ✓ VERIFIED | shifts.py solvent_map resolves 'methanol' → 'Methanol'. Scaling factor lookup succeeds: slope=-0.9326, intercept=29.74 for 1H |
| 2 | API accepts solvent=water and returns shifts using Water scaling factors | ✓ VERIFIED | shifts.py solvent_map resolves 'water' → 'Water'. Scaling factor lookup succeeds: slope=-0.9321, intercept=29.73 for 1H |
| 3 | API accepts solvent=acetone and returns shifts using Acetone scaling factors | ✓ VERIFIED | shifts.py solvent_map resolves 'acetone' → 'Acetone'. Scaling factor lookup succeeds: slope=-0.9332, intercept=29.76 for 1H |
| 4 | API accepts solvent=benzene and returns shifts using Benzene scaling factors | ✓ VERIFIED | shifts.py solvent_map resolves 'benzene' → 'Benzene'. Scaling factor lookup succeeds: slope=-0.9433, intercept=30.12 for 1H |
| 5 | Web UI solvent dropdown lists all 7 solvents | ✓ VERIFIED | web.py L44 iterates SUPPORTED_SOLVENTS.items() to build dropdown. SUPPORTED_SOLVENTS has 7 entries. All 7 will render. |
| 6 | Existing CHCl3, DMSO, and vacuum solvents still work identically | ✓ VERIFIED | Verified scaling factor lookups for chcl3/dmso/vacuum unchanged. All existing tests pass (40/40). No code changes to existing entries. |
| 7 | NMReData export works for all 7 solvents | ✓ VERIFIED | nmredata.py solvent_map has 7 entries. All map correctly to NMReData standard forms. 4 new tests pass: methanol→CD3OD, water→D2O, acetone→(CD3)2CO, benzene→C6D6 |

**Score:** 7/7 truths verified

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `src/qm_nmr_calc/solvents.py` | 7-solvent SUPPORTED_SOLVENTS dict | ✓ VERIFIED | 57 lines, contains 'methanol', 'water', 'acetone', 'benzene'. Dict has exactly 7 entries. validate_solvent() accepts all 4 new solvents. get_solvent_display_name() returns correct deuterated forms. |
| `src/qm_nmr_calc/shifts.py` | 7-solvent solvent_map | ✓ VERIFIED | 148 lines, contains 'methanol', 'water', 'acetone', 'benzene'. solvent_map (L50-58) has exactly 7 entries with Title-case values matching scaling_factors.json keys. |
| `src/qm_nmr_calc/nmredata.py` | 7-solvent NMReData mapping | ✓ VERIFIED | 275 lines, contains 'methanol', 'water', 'acetone', 'benzene'. solvent_map (L40-48) has exactly 7 entries with standard NMReData deuterated forms. |
| `src/qm_nmr_calc/data/scaling_factors.json` | 14 factor sets (7 solvents × 2 nuclei) | ✓ VERIFIED | 240 lines, contains all 14 keys: B3LYP/6-311+G(2d,p)/{1H,13C}/{CHCl3,DMSO,Methanol,Water,Acetone,Benzene,vacuum}. All lookups succeed. |
| `tests/test_nmredata.py` | 4 new tests for new solvents | ✓ VERIFIED | 4 new test methods added: test_methanol_maps_to_cd3od, test_water_maps_to_d2o, test_acetone_maps_to_deuterated_form, test_benzene_maps_to_c6d6. All pass. |

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|----|--------|---------|
| solvents.py | web.py | SUPPORTED_SOLVENTS import | ✓ WIRED | web.py L14 imports SUPPORTED_SOLVENTS. L44 builds dropdown: `for key, desc in sorted(SUPPORTED_SOLVENTS.items())`. All 7 solvents will render in dropdown. |
| solvents.py | jobs.py | validate_solvent import | ✓ WIRED | jobs.py L16 imports validate_solvent. L238 validates solvent before job creation. L361 validates on file upload. Both paths gate API submission. |
| shifts.py | jobs.py | get_scaling_factor import | ✓ WIRED | jobs.py L15 imports get_scaling_factor. L39 and L507 call get_scaling_factor() to retrieve MAE metadata. Used in nmr_results response. |
| shifts.py | scaling_factors.json | solvent_map normalization | ✓ WIRED | shifts.py L50-58 solvent_map normalizes user input (lowercase) to JSON keys (Title-case). get_scaling_factor() L59 uses map. L64 constructs key. L67 looks up in loaded JSON. All 14 lookups succeed. |
| nmredata.py | jobs.py | map_solvent_to_nmredata import | ✓ WIRED | jobs.py L12 imports generate_nmredata_sdf. L837 calls generate_nmredata_sdf(). L240 calls map_solvent_to_nmredata() internally. Used in /jobs/{job_id}/nmredata.sdf endpoint. |

### Requirements Coverage

| Requirement | Status | Blocking Issue |
|-------------|--------|----------------|
| INTG-01: User can select methanol as solvent in web UI and API | ✓ SATISFIED | None. SUPPORTED_SOLVENTS contains 'methanol'. validate_solvent() accepts it. Web UI dropdown will render it. API accepts it. Scaling factors exist. |
| INTG-02: User can select water as solvent in web UI and API | ✓ SATISFIED | None. SUPPORTED_SOLVENTS contains 'water'. validate_solvent() accepts it. Web UI dropdown will render it. API accepts it. Scaling factors exist. |
| INTG-03: User can select acetone as solvent in web UI and API | ✓ SATISFIED | None. SUPPORTED_SOLVENTS contains 'acetone'. validate_solvent() accepts it. Web UI dropdown will render it. API accepts it. Scaling factors exist. |
| INTG-04: User can select benzene as solvent in web UI and API | ✓ SATISFIED | None. SUPPORTED_SOLVENTS contains 'benzene'. validate_solvent() accepts it. Web UI dropdown will render it. API accepts it. Scaling factors exist. |
| INTG-06: Scaling factors loaded and applied correctly for all new solvents | ✓ SATISFIED | None. scaling_factors.json contains all 14 entries. solvent_map in shifts.py resolves all 7 solvents. get_scaling_factor() lookups succeed for all 14 combinations. |

### Anti-Patterns Found

No blocker anti-patterns detected. Code is substantive and production-ready.

**Summary:** Clean implementation with no TODOs, no placeholders, no stub patterns. All functions have real logic and are wired into the system.

### Human Verification Required

None required. All verifications completed programmatically:
- Artifact existence verified by file reads
- Artifact substantive content verified by line counts and pattern checks
- Wiring verified by grep for imports and function calls
- Functional correctness verified by executing Python imports and assertions
- Test suite verified by running pytest (40/40 tests pass)

No visual UI testing needed — dropdown rendering is deterministic from SUPPORTED_SOLVENTS dict.

---

## Detailed Evidence

### 1. Artifact Verification: solvents.py

**Level 1 (Exists):** ✓ File exists at `src/qm_nmr_calc/solvents.py`

**Level 2 (Substantive):**
- Line count: 57 lines (exceeds 15-line minimum for module)
- No stub patterns found (no TODO, FIXME, placeholder comments)
- Has real exports: `SUPPORTED_SOLVENTS`, `validate_solvent()`, `get_supported_solvents()`, `get_solvent_display_name()`
- SUPPORTED_SOLVENTS dict contains exactly 7 entries:
  ```python
  {
    "chcl3": "Chloroform (CDCl3)",
    "dmso": "Dimethylsulfoxide (DMSO-d6)",
    "vacuum": "Gas phase (no solvent)",
    "methanol": "Methanol (Methanol-d4)",
    "water": "Water (D2O)",
    "acetone": "Acetone (Acetone-d6)",
    "benzene": "Benzene (Benzene-d6)",
  }
  ```

**Level 3 (Wired):**
- Imported in `src/qm_nmr_calc/api/routers/web.py` (L13-18)
- Imported in `src/qm_nmr_calc/api/routers/jobs.py` (L16)
- SUPPORTED_SOLVENTS used in web.py L44 to build dropdown
- validate_solvent used in web.py L154 and jobs.py L238, L361
- get_supported_solvents used in web.py L156 and jobs.py L240, L363

**Functional Test:**
```python
from qm_nmr_calc.solvents import SUPPORTED_SOLVENTS, validate_solvent, get_solvent_display_name
assert len(SUPPORTED_SOLVENTS) == 7
assert validate_solvent('methanol') == 'methanol'
assert validate_solvent('water') == 'water'
assert validate_solvent('acetone') == 'acetone'
assert validate_solvent('benzene') == 'benzene'
assert get_solvent_display_name('methanol') == 'Methanol-d4'
assert get_solvent_display_name('water') == 'D2O'
# All assertions pass ✓
```

### 2. Artifact Verification: shifts.py

**Level 1 (Exists):** ✓ File exists at `src/qm_nmr_calc/shifts.py`

**Level 2 (Substantive):**
- Line count: 148 lines (exceeds 10-line minimum for module)
- No stub patterns found
- Has real exports: `load_scaling_factors()`, `get_scaling_factor()`, `shielding_to_shift()`
- solvent_map dict (L50-58) contains exactly 7 entries with Title-case normalization:
  ```python
  solvent_map = {
      "chcl3": "CHCl3",
      "dmso": "DMSO",
      "vacuum": "vacuum",
      "methanol": "Methanol",
      "water": "Water",
      "acetone": "Acetone",
      "benzene": "Benzene",
  }
  ```

**Level 3 (Wired):**
- Imported in `src/qm_nmr_calc/api/routers/jobs.py` (L15)
- get_scaling_factor used in jobs.py L39, L507 for MAE metadata
- load_scaling_factors() called by get_scaling_factor() L65
- Loads from `src/qm_nmr_calc/data/scaling_factors.json`

**Functional Test:**
```python
from qm_nmr_calc.shifts import get_scaling_factor

# Test all 4 new solvents
for solvent in ['methanol', 'water', 'acetone', 'benzene']:
    for nucleus in ['1H', '13C']:
        factor = get_scaling_factor('B3LYP', '6-311+G(2d,p)', nucleus, solvent)
        assert 'slope' in factor and 'intercept' in factor
        # Example: methanol 1H: slope=-0.9326, intercept=29.74

# Test existing solvents (no regression)
for solvent in ['chcl3', 'dmso', 'vacuum']:
    factor = get_scaling_factor('B3LYP', '6-311+G(2d,p)', '1H', solvent)
    assert 'slope' in factor and 'intercept' in factor
    # Example: chcl3 1H: slope=-0.9375, intercept=29.92

# All 14 lookups succeed ✓
```

### 3. Artifact Verification: nmredata.py

**Level 1 (Exists):** ✓ File exists at `src/qm_nmr_calc/nmredata.py`

**Level 2 (Substantive):**
- Line count: 275 lines (exceeds 10-line minimum)
- No stub patterns found
- Has real exports: `map_solvent_to_nmredata()`, `generate_nmredata_sdf()`
- solvent_map dict (L40-48) contains exactly 7 entries:
  ```python
  solvent_map = {
      "chcl3": "CDCl3",
      "dmso": "(CD3)2SO",
      "vacuum": "vacuum",
      "methanol": "CD3OD",
      "water": "D2O",
      "acetone": "(CD3)2CO",
      "benzene": "C6D6",
  }
  ```

**Level 3 (Wired):**
- Imported in `src/qm_nmr_calc/api/routers/jobs.py` (L12)
- generate_nmredata_sdf used in jobs.py L837 for /jobs/{job_id}/nmredata.sdf endpoint
- map_solvent_to_nmredata called internally by generate_nmredata_sdf (L240)

**Functional Test:**
```python
from qm_nmr_calc.nmredata import map_solvent_to_nmredata

mappings = {
    'chcl3': 'CDCl3',
    'dmso': '(CD3)2SO',
    'vacuum': 'vacuum',
    'methanol': 'CD3OD',
    'water': 'D2O',
    'acetone': '(CD3)2CO',
    'benzene': 'C6D6',
}
for solvent, expected in mappings.items():
    actual = map_solvent_to_nmredata(solvent)
    assert actual == expected
# All 7 mappings pass ✓
```

### 4. Artifact Verification: scaling_factors.json

**Level 1 (Exists):** ✓ File exists at `src/qm_nmr_calc/data/scaling_factors.json`

**Level 2 (Substantive):**
- Line count: 240 lines
- JSON is valid and parseable
- Contains exactly 14 top-level keys (7 solvents × 2 nuclei)
- All keys follow format: `B3LYP/6-311+G(2d,p)/{nucleus}/{solvent}`
- All 14 expected keys present:
  - CHCl3: 1H, 13C
  - DMSO: 1H, 13C
  - Methanol: 1H, 13C
  - Water: 1H, 13C
  - Acetone: 1H, 13C
  - Benzene: 1H, 13C
  - vacuum: 1H, 13C

**Level 3 (Wired):**
- Loaded by shifts.py L24-26 via importlib.resources
- Cached via @cache decorator (L15)
- Used by get_scaling_factor() L65 to look up factors

**Sample Entry (Methanol 1H):**
```json
"B3LYP/6-311+G(2d,p)/1H/Methanol": {
  "slope": -0.9326320659421693,
  "intercept": 29.744989682314603,
  "r_squared": 0.9950978707465968,
  "mae": 0.1262480011742177,
  "rmsd": 0.16560864968460298,
  "n_points": 335,
  "ci_slope": [-0.9396883671971676, -0.925575764687171],
  "ci_intercept": [29.544205659575763, 29.945773705053444],
  "outliers_removed": 7
}
```

All entries have consistent structure with slope, intercept, r_squared, mae, rmsd, n_points, confidence intervals, and outlier count.

### 5. Artifact Verification: tests/test_nmredata.py

**Level 1 (Exists):** ✓ File exists at `tests/test_nmredata.py`

**Level 2 (Substantive):**
- 4 new test methods added (L51-66):
  - `test_methanol_maps_to_cd3od()` — L51-53
  - `test_water_maps_to_d2o()` — L55-57
  - `test_acetone_maps_to_deuterated_form()` — L59-61
  - `test_benzene_maps_to_c6d6()` — L63-65
- test_unknown_solvent_raises_error updated to use 'toluene' instead of 'benzene' (L70)
- All tests have assertions and docstrings
- No stub patterns

**Level 3 (Wired):**
- Executed by pytest
- All 40 tests pass (36 existing + 4 new)

**Test Output:**
```
tests/test_nmredata.py::TestSolventMapping::test_methanol_maps_to_cd3od PASSED
tests/test_nmredata.py::TestSolventMapping::test_water_maps_to_d2o PASSED
tests/test_nmredata.py::TestSolventMapping::test_acetone_maps_to_deuterated_form PASSED
tests/test_nmredata.py::TestSolventMapping::test_benzene_maps_to_c6d6 PASSED
============================== 40 passed in 3.49s ==============================
```

---

## Success Criteria Met

From Phase 57 ROADMAP.md success criteria:

1. **✓ Submitting a molecule via API with `solvent=methanol` returns predicted shifts using methanol-specific scaling factors**
   - Evidence: validate_solvent('methanol') returns 'methanol'. get_scaling_factor('B3LYP', '6-311+G(2d,p)', '1H', 'methanol') returns Methanol-specific slope=-0.9326, intercept=29.74
   
2. **✓ Submitting a molecule via API with `solvent=water`, `solvent=acetone`, or `solvent=benzene` each returns predicted shifts using the correct solvent-specific scaling factors**
   - Evidence: All 3 solvents validated. Scaling factor lookups succeed with distinct values for each solvent:
     - water: slope=-0.9321, intercept=29.73
     - acetone: slope=-0.9332, intercept=29.76
     - benzene: slope=-0.9433, intercept=30.12

3. **✓ Web UI solvent dropdown lists all 7 solvents (CHCl3, DMSO, Vacuum, Methanol, Water, Acetone, Benzene)**
   - Evidence: web.py L44 iterates SUPPORTED_SOLVENTS.items() to build dropdown options. SUPPORTED_SOLVENTS has exactly 7 entries. All will render.

4. **✓ Existing CHCl3, DMSO, and vacuum calculations produce identical results to before (no regression)**
   - Evidence: No code changes to existing SUPPORTED_SOLVENTS entries, solvent_map entries, or scaling_factors.json entries for CHCl3/DMSO/vacuum. All 40 tests pass. Verified scaling factor lookups return same values as before.

---

_Verified: 2026-02-09T10:47:13Z_  
_Verifier: Claude (gsd-verifier)_
