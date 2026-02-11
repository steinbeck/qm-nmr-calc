---
phase: 64-solvent-integration
verified: 2026-02-11T09:11:20Z
status: passed
score: 9/9 must-haves verified
---

# Phase 64: Solvent Integration Verification Report

**Phase Goal:** Add 6 new solvents (pyridine, THF, toluene, DCM, acetonitrile, DMF) to the 3 gatekeeper modules so users can select them in the web UI and API for NMR predictions.

**Verified:** 2026-02-11T09:11:20Z
**Status:** PASSED
**Re-verification:** No — initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | API accepts solvent=pyridine and returns shifts using Pyridine scaling factors | ✓ VERIFIED | `get_scaling_factor()` successfully resolves pyridine → Pyridine → B3LYP/6-311+G(2d,p)/1H/Pyridine (slope=-0.9340) |
| 2 | API accepts solvent=thf and returns shifts using THF scaling factors | ✓ VERIFIED | `get_scaling_factor()` successfully resolves thf → THF → B3LYP/6-311+G(2d,p)/1H/THF (slope=-0.9356) |
| 3 | API accepts solvent=toluene and returns shifts using Toluene scaling factors | ✓ VERIFIED | `get_scaling_factor()` successfully resolves toluene → Toluene → B3LYP/6-311+G(2d,p)/1H/Toluene (slope=-0.9429) |
| 4 | API accepts solvent=dcm and returns shifts using DCM scaling factors | ✓ VERIFIED | `get_scaling_factor()` successfully resolves dcm → DCM → B3LYP/6-311+G(2d,p)/1H/DCM (slope=-0.9350) |
| 5 | API accepts solvent=acetonitrile and returns shifts using Acetonitrile scaling factors | ✓ VERIFIED | `get_scaling_factor()` successfully resolves acetonitrile → Acetonitrile → B3LYP/6-311+G(2d,p)/1H/Acetonitrile (slope=-0.9326) |
| 6 | API accepts solvent=dmf and returns shifts using DMF scaling factors | ✓ VERIFIED | `get_scaling_factor()` successfully resolves dmf → DMF → B3LYP/6-311+G(2d,p)/1H/DMF (slope=-0.9325) |
| 7 | Web UI solvent dropdown lists all 13 solvents (12 visible + vacuum) | ✓ VERIFIED | `SUPPORTED_SOLVENTS` has 13 entries, web.py sorts and displays all (12 visible + vacuum) |
| 8 | Existing CHCl3, DMSO, vacuum, methanol, water, acetone, benzene solvents still work identically | ✓ VERIFIED | All 7 legacy solvents pass validation, scaling factor lookup, and NMReData mapping (regression test passed) |
| 9 | NMReData export works for all 13 solvents with correct deuterated formulas | ✓ VERIFIED | `map_solvent_to_nmredata()` successfully maps all 13 solvents to correct deuterated formulas (e.g., pyridine → C5D5N) |

**Score:** 9/9 truths verified (100%)

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `src/qm_nmr_calc/solvents.py` | 13-solvent SUPPORTED_SOLVENTS dict containing "pyridine" | ✓ VERIFIED | EXISTS (62 lines), SUBSTANTIVE (13 entries: acetone, acetonitrile, benzene, chcl3, dcm, dmf, dmso, methanol, pyridine, thf, toluene, vacuum, water), WIRED (imported by web.py, jobs.py, input_gen.py) |
| `src/qm_nmr_calc/shifts.py` | 13-solvent solvent_map for scaling factor lookup containing "pyridine" | ✓ VERIFIED | EXISTS (153 lines), SUBSTANTIVE (13 entries in solvent_map dict, all map to Title-case keys matching scaling_factors.json), WIRED (used by jobs.py for factor lookup) |
| `src/qm_nmr_calc/nmredata.py` | 13-solvent NMReData solvent mapping containing "pyridine" | ✓ VERIFIED | EXISTS (280 lines), SUBSTANTIVE (13 entries in solvent_map with correct deuterated formulas), WIRED (imported by jobs.py for SDF export) |
| `tests/test_nmredata.py` | 6 new solvent mapping tests + updated unknown solvent test containing "test_pyridine_maps_to_c5d5n" | ✓ VERIFIED | EXISTS (474 lines), SUBSTANTIVE (46 tests total including 6 new: test_pyridine_maps_to_c5d5n, test_thf_maps_to_c4d8o, test_toluene_maps_to_c7d8, test_dcm_maps_to_cd2cl2, test_acetonitrile_maps_to_cd3cn, test_dmf_maps_to_deuterated_form + updated test_unknown_solvent_raises_error using "ethanol"), WIRED (all 46 tests pass) |

All artifacts verified at 3 levels: EXISTS ✓, SUBSTANTIVE ✓, WIRED ✓

### Key Link Verification

| From | To | Via | Status | Details |
|------|-----|-----|--------|---------|
| `src/qm_nmr_calc/solvents.py` | `src/qm_nmr_calc/api/routers/web.py` | SUPPORTED_SOLVENTS import drives dropdown rendering | ✓ WIRED | web.py line 14 imports SUPPORTED_SOLVENTS, line 44 uses it to build dropdown with sorted(SUPPORTED_SOLVENTS.items()). Verified: 13 solvents displayed. |
| `src/qm_nmr_calc/shifts.py` | `src/qm_nmr_calc/data/scaling_factors.json` | solvent_map normalizes keys to match JSON keys | ✓ WIRED | solvent_map (line 50-64) maps lowercase solvent names to Title-case keys. Verified: All 26 keys exist in scaling_factors.json (13 solvents x 2 nuclei). |
| `src/qm_nmr_calc/solvents.py` | `src/qm_nmr_calc/api/routers/jobs.py` | validate_solvent gates API submission | ✓ WIRED | jobs.py line 16 imports validate_solvent. Verified: All 13 solvents pass validation (returns solvent name). |
| `src/qm_nmr_calc/nmredata.py` | `src/qm_nmr_calc/api/routers/jobs.py` | NMReData export for SDF generation | ✓ WIRED | jobs.py line 12 imports generate_nmredata_sdf. Verified: map_solvent_to_nmredata() handles all 13 solvents. |
| `src/qm_nmr_calc/solvents.py` | `src/qm_nmr_calc/nwchem/input_gen.py` | COSMO solvent parameter generation | ✓ WIRED | input_gen.py line 8 has SUPPORTED_SOLVENTS set with all 13 solvents including special acetonitrile → acetntrl mapping. |

All key links verified as WIRED ✓

### Requirements Coverage

No explicit requirements file mapped to this phase. The phase goal encompasses all integration requirements.

### Anti-Patterns Found

**No anti-patterns detected.**

Scanned files:
- `src/qm_nmr_calc/solvents.py` (62 lines)
- `src/qm_nmr_calc/shifts.py` (153 lines)
- `src/qm_nmr_calc/nmredata.py` (280 lines)
- `tests/test_nmredata.py` (474 lines)

Checks performed:
- No TODO/FIXME/PLACEHOLDER comments
- No stub implementations (empty returns, console.log only)
- No hardcoded placeholders
- All functions have substantive implementations
- All new solvents have real data (not placeholders)

### Detailed Verification Results

#### 1. SUPPORTED_SOLVENTS Count and Content

```bash
$ python -c "from qm_nmr_calc.solvents import SUPPORTED_SOLVENTS; print(len(SUPPORTED_SOLVENTS))"
13

$ python -c "from qm_nmr_calc.solvents import SUPPORTED_SOLVENTS; print(sorted(SUPPORTED_SOLVENTS.keys()))"
['acetone', 'acetonitrile', 'benzene', 'chcl3', 'dcm', 'dmf', 'dmso', 'methanol', 'pyridine', 'thf', 'toluene', 'vacuum', 'water']
```

**Result:** ✓ All 13 solvents present (7 original + 6 new)

#### 2. Solvent Validation (API Gatekeeper)

```bash
$ python -c "from qm_nmr_calc.solvents import validate_solvent; new_solvents = ['pyridine', 'thf', 'toluene', 'dcm', 'acetonitrile', 'dmf']; results = [(s, validate_solvent(s)) for s in new_solvents]; [print(f'{s}: {r}') for s, r in results]"
pyridine: pyridine
thf: thf
toluene: toluene
dcm: dcm
acetonitrile: acetonitrile
dmf: dmf
```

**Result:** ✓ All 6 new solvents pass validation

#### 3. Display Name Extraction

```bash
$ python -c "from qm_nmr_calc.solvents import get_solvent_display_name; new_solvents = ['pyridine', 'thf', 'toluene', 'dcm', 'acetonitrile', 'dmf']; results = [(s, get_solvent_display_name(s)) for s in new_solvents]; [print(f'{s} -> {r}') for s, r in results]"
pyridine -> Pyridine-d5
thf -> THF-d8
toluene -> Toluene-d8
dcm -> CD2Cl2
acetonitrile -> CD3CN
dmf -> DMF-d7
```

**Result:** ✓ Correct deuterated display names extracted from parentheses

#### 4. Scaling Factor Lookup (26 Total)

```bash
$ python -c "from qm_nmr_calc.shifts import get_scaling_factor; [get_scaling_factor('B3LYP', '6-311+G(2d,p)', n, s) for n in ['1H', '13C'] for s in ['chcl3', 'dmso', 'vacuum', 'methanol', 'water', 'acetone', 'benzene', 'pyridine', 'thf', 'toluene', 'dcm', 'acetonitrile', 'dmf']]; print('All 26 factor lookups succeed')"
All 26 factor lookups succeed
```

**Detailed check:** All 26 keys found in scaling_factors.json:
- B3LYP/6-311+G(2d,p)/1H/{CHCl3, DMSO, vacuum, Methanol, Water, Acetone, Benzene, Pyridine, THF, Toluene, DCM, Acetonitrile, DMF}
- B3LYP/6-311+G(2d,p)/13C/{same 13 solvents}

**Sample slopes for new solvents:**
- pyridine: slope=-0.9340
- thf: slope=-0.9356
- toluene: slope=-0.9429
- dcm: slope=-0.9350
- acetonitrile: slope=-0.9326
- dmf: slope=-0.9325

**Result:** ✓ All scaling factors resolve correctly

#### 5. NMReData Solvent Mapping

```bash
$ python -c "from qm_nmr_calc.nmredata import map_solvent_to_nmredata; [map_solvent_to_nmredata(s) for s in ['chcl3', 'dmso', 'vacuum', 'methanol', 'water', 'acetone', 'benzene', 'pyridine', 'thf', 'toluene', 'dcm', 'acetonitrile', 'dmf']]; print('All 13 NMReData mappings work')"
All 13 NMReData mappings work
```

**Mappings for new solvents:**
- pyridine → C5D5N
- thf → C4D8O
- toluene → C7D8
- dcm → CD2Cl2
- acetonitrile → CD3CN
- dmf → (CD3)2NCDO

**Result:** ✓ All NMReData deuterated formulas correct

#### 6. Web UI Dropdown Content

Simulated web.py dropdown rendering:
```
Web UI will display 13 solvents in dropdown (sorted by label):
 1. acetone         -> Acetone (Acetone-d6)
 2. acetonitrile    -> Acetonitrile (CD3CN)
 3. benzene         -> Benzene (Benzene-d6)
 4. chcl3           -> Chloroform (CDCl3)
 5. dcm             -> Dichloromethane (CD2Cl2)
 6. dmso            -> Dimethylsulfoxide (DMSO-d6)
 7. vacuum          -> Gas phase (no solvent)
 8. methanol        -> Methanol (Methanol-d4)
 9. dmf             -> N,N-Dimethylformamide (DMF-d7)
10. pyridine        -> Pyridine (Pyridine-d5)
11. thf             -> Tetrahydrofuran (THF-d8)
12. toluene         -> Toluene (Toluene-d8)
13. water           -> Water (D2O)
```

**Result:** ✓ All 13 solvents (12 visible + vacuum) displayed correctly

#### 7. Regression Testing (Legacy Solvents)

**Validation:**
```
chcl3: chcl3 ✓
dmso: dmso ✓
vacuum: vacuum ✓
methanol: methanol ✓
water: water ✓
acetone: acetone ✓
benzene: benzene ✓
```

**Scaling factors (1H):**
```
chcl3: slope=-0.9375 ✓
dmso: slope=-0.9323 ✓
vacuum: slope=-0.9554 ✓
methanol: slope=-0.9326 ✓
water: slope=-0.9321 ✓
acetone: slope=-0.9332 ✓
benzene: slope=-0.9433 ✓
```

**NMReData mappings:**
```
chcl3 -> CDCl3 ✓
dmso -> (CD3)2SO ✓
vacuum -> vacuum ✓
methanol -> CD3OD ✓
water -> D2O ✓
acetone -> (CD3)2CO ✓
benzene -> C6D6 ✓
```

**Result:** ✓ All 7 legacy solvents work identically (no regressions)

#### 8. Test Suite Results

**nmredata tests:**
```
tests/test_nmredata.py::TestSolventMapping::test_pyridine_maps_to_c5d5n PASSED
tests/test_nmredata.py::TestSolventMapping::test_thf_maps_to_c4d8o PASSED
tests/test_nmredata.py::TestSolventMapping::test_toluene_maps_to_c7d8 PASSED
tests/test_nmredata.py::TestSolventMapping::test_dcm_maps_to_cd2cl2 PASSED
tests/test_nmredata.py::TestSolventMapping::test_acetonitrile_maps_to_cd3cn PASSED
tests/test_nmredata.py::TestSolventMapping::test_dmf_maps_to_deuterated_form PASSED
tests/test_nmredata.py::TestSolventMapping::test_unknown_solvent_raises_error PASSED
========================= 46 tests passed in 4.01s ==========================
```

**Key changes:**
- Added 6 new test methods (lines 67-89 in test_nmredata.py)
- Updated test_unknown_solvent_raises_error to use "ethanol" instead of "toluene" (line 94)
- Total test count: 46 (40 existing + 6 new)

**Result:** ✓ All tests pass

#### 9. NWChem COSMO Integration

Verified `src/qm_nmr_calc/nwchem/input_gen.py` line 8:
```python
SUPPORTED_SOLVENTS = {"chcl3", "dmso", "water", "acetone", "methanol", "benzene", "vacuum", "pyridine", "thf", "toluene", "dcm", "acetonitrile", "dmf"}
```

Special mapping for NWChem COSMO (line 12-14):
```python
COSMO_NAME_MAP: dict[str, str] = {
    "acetonitrile": "acetntrl",  # NWChem uses abbreviated name
}
```

**Result:** ✓ All 13 solvents supported by NWChem input generator

### Conclusion

**Phase goal ACHIEVED.**

All 9 must-have truths verified through programmatic checks and test execution. The 6 new solvents (pyridine, THF, toluene, DCM, acetonitrile, DMF) are now fully integrated into:

1. **Gatekeeper Module 1 (solvents.py):** SUPPORTED_SOLVENTS dict with 13 entries, validate_solvent() accepts all new solvents
2. **Gatekeeper Module 2 (shifts.py):** solvent_map with 13 entries, all resolve to correct scaling_factors.json keys
3. **Gatekeeper Module 3 (nmredata.py):** solvent_map with 13 entries, all map to correct deuterated formulas

All key links verified as wired:
- Web UI dropdown will display all 13 solvents (sorted alphabetically by label)
- API validation accepts all 13 solvents
- Scaling factor lookup succeeds for all 26 combinations (13 solvents x 2 nuclei)
- NMReData export works for all 13 solvents
- NWChem COSMO input generation supports all 13 solvents

No regressions detected: All 7 legacy solvents work identically.

No anti-patterns found: All implementations are substantive, no stubs or placeholders.

**Ready to proceed to next phase.**

---

_Verified: 2026-02-11T09:11:20Z_
_Verifier: Claude (gsd-verifier)_
