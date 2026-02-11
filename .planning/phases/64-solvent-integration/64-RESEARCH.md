# Phase 64: Solvent Integration - Research

**Researched:** 2026-02-11
**Domain:** NMR solvent integration pipeline (gatekeeper modules)
**Confidence:** HIGH

## Summary

Phase 64 follows the EXACT same pattern as Phase 57 (v2.8 solvent integration), which successfully wired 4 new solvents (Methanol, Water, Acetone, Benzene) through the production pipeline. This phase wires 6 new solvents (Pyridine, THF, Toluene, DCM, Acetonitrile, DMF) through the same pipeline.

The integration involves updating 3 "gatekeeper" modules that control solvent access across the entire application stack: `solvents.py` (validation/display), `shifts.py` (scaling factor lookup), and `nmredata.py` (export format). The UI and API automatically pick up changes from these modules via imports, requiring no direct modification.

Phase 63 (just completed) already added all 26 scaling factors (13 solvents × 2 nuclei) to `scaling_factors.json`. Phase 59 already added NWChem COSMO support for all 6 solvents to `input_gen.py`. This phase is purely wiring—connecting existing infrastructure.

**Primary recommendation:** Follow Phase 57 plan structure exactly, updating the same 3 modules with the same entry format. Add 6 solvent entries instead of 4. Update test that used "toluene" as unknown solvent (now valid).

## Standard Stack

Phase 57 established the integration pattern. No new libraries or tools required.

### Core Components
| Module | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| solvents.py | N/A | Solvent validation & display names | Central registry, imported by API & UI |
| shifts.py | N/A | Scaling factor lookup | Maps solvent names to JSON keys |
| nmredata.py | N/A | NMReData export format | Maps solvent names to deuterated forms |

### Dependencies (Already Present)
| Component | Location | Purpose | Integration Point |
|-----------|----------|---------|-------------------|
| NWChem COSMO | input_gen.py | Input generation | Already has all 6 (Phase 59) |
| Scaling factors | scaling_factors.json | Shift calculation | Already has 26 (Phase 63) |
| UI dropdown | submit.html | User selection | Auto-imports from solvents.py |
| API validation | routers/jobs.py | Request gating | Auto-imports from solvents.py |

**Key insight:** The 3 gatekeeper modules are the ONLY bottleneck. Everything else is already wired.

## Architecture Patterns

### Pattern 1: Central Registry with Auto-Propagation
**What:** Single source of truth (`SUPPORTED_SOLVENTS` dict in `solvents.py`) that propagates to all consumers via imports.

**How it works:**
```python
# solvents.py - Single source of truth
SUPPORTED_SOLVENTS: dict[str, str] = {
    "chcl3": "Chloroform (CDCl3)",
    "pyridine": "Pyridine (Pyridine-d5)",  # NEW
    # ...
}

# web.py - UI auto-picks up changes (line 44)
solvents = [
    {"value": key, "label": desc}
    for key, desc in sorted(SUPPORTED_SOLVENTS.items(), key=lambda x: x[1])
]

# jobs.py - API auto-validates via validate_solvent() import
# No changes needed - validation function reads SUPPORTED_SOLVENTS
```

**Why this works:** No hardcoded lists anywhere else. Add to `SUPPORTED_SOLVENTS` → UI dropdown updates, API accepts, validation passes.

### Pattern 2: Solvent Name Normalization Chain
**What:** 3-step mapping from user input → internal key → external format.

**The chain:**
1. **User input** → `validate_solvent()` → **normalized key** (lowercase)
2. **Normalized key** → `shifts.py:solvent_map` → **Title-case** (for JSON lookup)
3. **Normalized key** → `nmredata.py:solvent_map` → **Deuterated formula** (for export)

**Example for pyridine:**
```
User: "Pyridine" or "pyridine" or "PYRIDINE"
  ↓ validate_solvent()
Internal: "pyridine" (stored in job, passed to all functions)
  ↓ shifts.py:solvent_map
JSON key: "Pyridine" (matches "B3LYP/6-311+G(2d,p)/1H/Pyridine")
  ↓ nmredata.py:solvent_map
NMReData: "C5D5N" (deuterated chemical formula)
```

**Critical invariants:**
- `solvents.py` keys: lowercase, match NWChem names (except acetonitrile)
- `shifts.py` values: Title-case, match JSON keys EXACTLY
- `nmredata.py` values: Chemical formulas with deuterium (D)

### Pattern 3: Display Name Extraction
**What:** `SUPPORTED_SOLVENTS` values encode BOTH full name AND display name using parentheses.

**Format:** `"Full Name (Display Name)"`

**Examples:**
```python
"chcl3": "Chloroform (CDCl3)"         # Display: "CDCl3"
"pyridine": "Pyridine (Pyridine-d5)"  # Display: "Pyridine-d5"
"dmf": "N,N-Dimethylformamide (DMF-d7)"  # Display: "DMF-d7"
```

**Extraction:** `get_solvent_display_name()` parses text in parentheses (lines 42-56 of solvents.py).

**Where used:**
- UI dropdown shows full description: "Pyridine (Pyridine-d5)"
- Job results show short form: "Pyridine-d5"
- Both derived from same source → guaranteed consistency

## Don't Hand-Roll

Phase 57 established the pattern. DO NOT deviate.

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Solvent validation | Custom validator per endpoint | `validate_solvent()` from solvents.py | Single source prevents drift |
| UI dropdown options | Hardcoded `<option>` list | Loop over `SUPPORTED_SOLVENTS` | Auto-updates when module changes |
| Display names | String manipulation in templates | `get_solvent_display_name()` | Consistent extraction logic |
| Scaling factor lookup | Direct JSON key construction | `get_scaling_factor()` with solvent_map | Handles case normalization |
| NMReData export | Custom deuterated name logic | `map_solvent_to_nmredata()` | Centralized mapping, tested |

**Key insight:** The gatekeeper pattern prevents fragmentation. Every access point goes through the same 3 modules. Adding a solvent updates ALL access points simultaneously.

## Common Pitfalls

### Pitfall 1: Mismatched Capitalization in JSON Keys
**What goes wrong:** Scaling factor lookup fails with "No scaling factor for B3LYP/6-311+G(2d,p)/1H/pyridine" even though file has "B3LYP/6-311+G(2d,p)/1H/Pyridine".

**Why it happens:** JSON keys use Title-case ("Pyridine"), but internal code uses lowercase ("pyridine"). The `solvent_map` in `shifts.py` performs this normalization.

**How to avoid:**
1. Check `scaling_factors.json` for exact key format: `grep "1H/" src/qm_nmr_calc/data/scaling_factors.json`
2. Add `solvent_map` entry with exact capitalization match
3. Verify with: `get_scaling_factor('B3LYP', '6-311+G(2d,p)', '1H', 'pyridine')`

**Warning signs:** KeyError or ValueError mentioning scaling factor lookup during shifts calculation.

### Pitfall 2: Display Name Format Inconsistency
**What goes wrong:** UI shows "Pyridine" instead of "Pyridine-d5" because parentheses are missing or malformed.

**Why it happens:** `get_solvent_display_name()` extracts text between `(` and `)`. If format is wrong, extraction fails silently and returns the key.

**How to avoid:**
1. Follow exact pattern: `"Full Name (Display Name)"`
2. Ensure closing parenthesis: `"Pyridine (Pyridine-d5)"`
3. Test extraction: `get_solvent_display_name('pyridine')` should return `"Pyridine-d5"`

**Warning signs:** Job results show lowercase internal names instead of formatted display names.

### Pitfall 3: NMReData Chemical Formula Errors
**What goes wrong:** NMReData validation tools reject SDF file because solvent name doesn't match chemical formula convention.

**Why it happens:** NMReData expects specific deuterated chemical formulas (e.g., "C5D5N" for pyridine-d5, not "Pyridine-d5").

**How to avoid:**
1. Use chemical formulas with deuterium (D): "C5D5N", "C4D8O", "C7D8"
2. For complex molecules, use parentheses for groups: "(CD3)2CO", "(CD3)2NCDO"
3. Verify against existing patterns in `nmredata.py` (lines 40-48)
4. Cross-reference with NMR literature for standard deuterated formulas

**Warning signs:** NMReData export generates but fails external validation; users report SDF incompatibility.

### Pitfall 4: Test Breakage from "Unknown" Solvent Becoming Valid
**What goes wrong:** `test_unknown_solvent_raises_error` in `test_nmredata.py` fails because it tests "toluene" which is now valid (line 70).

**Why it happens:** Tests use real solvent names as examples of "unknown" solvents. When that solvent gets added, the test premise breaks.

**How to avoid:**
1. Update test to use genuinely unknown solvent (e.g., "ethanol", "hexane")
2. Check ALL test files for uses of new solvent names as test data
3. Run full test suite BEFORE committing

**Warning signs:** `test_unknown_solvent_raises_error` fails with "Expected ValueError but got 'Toluene'".

### Pitfall 5: Forgetting Vacuum in Solvent Count
**What goes wrong:** Documentation says "14 solvents" but UI shows 13. User confusion.

**Why it happens:** "vacuum" (gas phase) is internal but not exposed in UI dropdown (it's valid but hidden from users).

**How to avoid:**
1. When counting solvents, specify "13 visible in UI, 14 total including vacuum"
2. Document why vacuum is special (gas phase, no COSMO, internal use)
3. UI tests should verify 13 dropdown options (excluding vacuum)

**Warning signs:** Off-by-one count in documentation, UI, or tests.

## Code Examples

Verified patterns from Phase 57 and current codebase:

### Adding 6 Solvents to solvents.py
```python
# Source: Phase 57 plan, lines 81-91
# Pattern: "key": "Full Name (Display Name)"
SUPPORTED_SOLVENTS: dict[str, str] = {
    "chcl3": "Chloroform (CDCl3)",
    "dmso": "Dimethylsulfoxide (DMSO-d6)",
    "vacuum": "Gas phase (no solvent)",
    "methanol": "Methanol (Methanol-d4)",
    "water": "Water (D2O)",
    "acetone": "Acetone (Acetone-d6)",
    "benzene": "Benzene (Benzene-d6)",
    # NEW: Phase 64 additions
    "pyridine": "Pyridine (Pyridine-d5)",
    "thf": "Tetrahydrofuran (THF-d8)",
    "toluene": "Toluene (Toluene-d8)",
    "dcm": "Dichloromethane (CD2Cl2)",
    "acetonitrile": "Acetonitrile (CD3CN)",
    "dmf": "N,N-Dimethylformamide (DMF-d7)",
}
```

**Key points:**
- Keys: lowercase, match NWChem names
- Display names in parentheses: match NMR community conventions
- Order: keep existing entries first, append new ones

### Adding 6 Solvents to shifts.py
```python
# Source: Phase 57 plan, lines 96-112
# Pattern: "lowercase": "Title-case"
# Values MUST match scaling_factors.json keys EXACTLY
solvent_map = {
    "chcl3": "CHCl3",
    "dmso": "DMSO",
    "vacuum": "vacuum",
    "methanol": "Methanol",
    "water": "Water",
    "acetone": "Acetone",
    "benzene": "Benzene",
    # NEW: Phase 64 additions
    "pyridine": "Pyridine",
    "thf": "THF",
    "toluene": "Toluene",
    "dcm": "DCM",
    "acetonitrile": "Acetonitrile",
    "dmf": "DMF",
}
```

**Critical:** These values construct JSON keys like "B3LYP/6-311+G(2d,p)/1H/Pyridine". Capitalization MUST match JSON exactly (verified via grep above).

### Adding 6 Solvents to nmredata.py
```python
# Source: Phase 57 plan, lines 138-148; current nmredata.py lines 40-48
# Pattern: "lowercase": "Chemical Formula with D"
solvent_map = {
    "chcl3": "CDCl3",
    "dmso": "(CD3)2SO",
    "vacuum": "vacuum",
    "methanol": "CD3OD",
    "water": "D2O",
    "acetone": "(CD3)2CO",
    "benzene": "C6D6",
    # NEW: Phase 64 additions
    "pyridine": "C5D5N",           # Pyridine: C5H5N → C5D5N
    "thf": "C4D8O",                # THF: C4H8O → C4D8O
    "toluene": "C7D8",             # Toluene: C7H8 → C7D8
    "dcm": "CD2Cl2",               # DCM: CH2Cl2 → CD2Cl2
    "acetonitrile": "CD3CN",       # Acetonitrile: CH3CN → CD3CN
    "dmf": "(CD3)2NCDO",           # DMF: (CH3)2NCHO → (CD3)2NCDO
}
```

**Key points:**
- Use chemical formulas, not common names
- Replace all H with D (deuterium)
- Use parentheses for molecular groups: "(CD3)2NCDO"
- Verified against NMR solvent data charts and WebSearch findings

### Test Update for test_nmredata.py
```python
# Source: Phase 57 plan, lines 152-157
# OLD: test used "benzene" as unknown solvent (now valid)
# NEW: use genuinely unknown solvent
def test_unknown_solvent_raises_error(self):
    """Unknown solvent should raise ValueError."""
    with pytest.raises(ValueError, match="Unknown solvent"):
        map_solvent_to_nmredata("ethanol")  # Changed from "toluene"
```

**Why:** "toluene" is now valid (Phase 64), so test needs different unknown solvent. "ethanol" not in any current phase roadmap.

### Adding Tests for New Solvents
```python
# Source: Phase 57 plan, lines 159-176
# Add 6 new test methods to TestSolventMapping class

def test_pyridine_maps_to_c5d5n(self):
    """Pyridine maps to deuterated form C5D5N."""
    assert map_solvent_to_nmredata("pyridine") == "C5D5N"

def test_thf_maps_to_c4d8o(self):
    """THF maps to deuterated form C4D8O."""
    assert map_solvent_to_nmredata("thf") == "C4D8O"

def test_toluene_maps_to_c7d8(self):
    """Toluene maps to deuterated form C7D8."""
    assert map_solvent_to_nmredata("toluene") == "C7D8"

def test_dcm_maps_to_cd2cl2(self):
    """DCM maps to deuterated form CD2Cl2."""
    assert map_solvent_to_nmredata("dcm") == "CD2Cl2"

def test_acetonitrile_maps_to_cd3cn(self):
    """Acetonitrile maps to deuterated form CD3CN."""
    assert map_solvent_to_nmredata("acetonitrile") == "CD3CN"

def test_dmf_maps_to_deuterated_form(self):
    """DMF maps to NMReData convention (CD3)2NCDO."""
    assert map_solvent_to_nmredata("dmf") == "(CD3)2NCDO"
```

**Insert location:** After existing solvent tests, before `test_unknown_solvent_raises_error` (currently line 67).

## State of the Art

| Component | Current State (Phase 63) | After Phase 64 | Notes |
|-----------|-------------------------|----------------|-------|
| NWChem COSMO | 13 solvents in input_gen.py | No change | Phase 59 completed this |
| Scaling factors | 26 factors in JSON | No change | Phase 63 completed this |
| API validation | 7 solvents in SUPPORTED_SOLVENTS | 13 solvents | This phase |
| UI dropdown | 7 options (6 visible) | 13 options (12 visible) | This phase |
| NMReData export | 7 mappings | 13 mappings | This phase |
| Test coverage | 7 solvents tested | 13 solvents tested | This phase |

**Key insight:** Infrastructure layers (NWChem, scaling factors) were completed in earlier phases. This phase is the final wiring to expose them to users.

**Deprecated/outdated:**
- None. Pattern from Phase 57 remains current best practice.

## Open Questions

**Q1: Should "vacuum" be exposed in UI dropdown?**
- **What we know:** "vacuum" exists in `SUPPORTED_SOLVENTS` but historically hidden from UI
- **What's unclear:** Whether users should be able to select gas-phase calculations via UI
- **Recommendation:** Keep hidden (status quo). Gas phase is for advanced users; API still accepts it. Document in API docs only.

**Q2: Do NMReData tools validate deuterated solvent names?**
- **What we know:** NMReData spec v1.1 defines NMREDATA_SOLVENT tag format as free text
- **What's unclear:** Whether external validators check chemical formula correctness
- **Recommendation:** Follow existing patterns (CDCl3, D2O, etc.) which are widely accepted. Test export with real NMReData tools if available.

**Q3: Should we add integration tests for all 13 solvents?**
- **What we know:** Current tests are unit tests on individual modules
- **What's unclear:** Whether we need end-to-end API tests for each solvent
- **Recommendation:** Not in this phase. Existing test pattern (unit tests per module) is sufficient. Phase 57 didn't add integration tests and worked fine.

## Sources

### Primary (HIGH confidence)
- Phase 57 plan: `.planning/phases/57-solvent-integration/57-01-PLAN.md`
- Current codebase: `src/qm_nmr_calc/solvents.py`, `shifts.py`, `nmredata.py`, `nwchem/input_gen.py`
- Scaling factors JSON: `src/qm_nmr_calc/data/scaling_factors.json` (verified via grep)
- Test file: `tests/test_nmredata.py`

### Secondary (MEDIUM confidence)
- [NWChem COSMO Solvation Model](https://nwchemgit.github.io/COSMO-Solvation-Model.html) - Official NWChem documentation listing supported solvent names
- [NMR Solvent Data Chart (Cambridge Isotope Laboratories)](https://chem.washington.edu/sites/chem/files/documents/facilities/nmrsolventschart_001.pdf) - Standard deuterated solvent naming conventions

### Tertiary (LOW confidence - WebSearch verified patterns)
- Multiple NMR solvent charts confirming deuterated formulas: Pyridine-d5 (C5D5N), THF-d8 (C4D8O), Toluene-d8 (C7D8), DCM-d2 (CD2Cl2), Acetonitrile-d3 (CD3CN), DMF-d7 (C3D7NO)
- Pattern from existing code: Use chemical formulas with D, not commercial names like "Pyridine-d5"

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - Established by Phase 57, no changes
- Architecture: HIGH - Current codebase implements these patterns consistently
- Code examples: HIGH - Directly from Phase 57 plan and current working code
- NWChem COSMO names: HIGH - Verified from official NWChem documentation
- Deuterated formulas: MEDIUM - Verified from multiple NMR sources and existing code patterns, but NMReData spec validation uncertain
- Pitfalls: HIGH - Identified from Phase 57 success and potential failure modes

**Research date:** 2026-02-11
**Valid until:** 2026-06-30 (120 days - stable codebase patterns, NMR conventions don't change rapidly)
