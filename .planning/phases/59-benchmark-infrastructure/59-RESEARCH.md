# Phase 59: Benchmark Infrastructure - Research

**Researched:** 2026-02-10
**Domain:** NWChem COSMO solvent extension, Python CLI validation, test infrastructure
**Confidence:** HIGH

## Summary

Phase 59 extends the benchmark CLI to support 6 additional solvents (pyridine, THF, toluene, DCM, acetonitrile, DMF) by adding them to the existing 3-module solvent gatekeeper pattern. The task is straightforward: add solvent names to validation lists and implement a COSMO name mapping layer to handle the one exception (acetonitrile → acetntrl).

The project already has a well-established pattern for solvent validation across three modules (solvents.py, shifts.py, nmredata.py), comprehensive test coverage using pytest with parametrized tests, and a working CLI validation system using argparse choices. The implementation is a mechanical extension of existing patterns rather than architectural design.

**Primary recommendation:** Follow the established 3-module pattern exactly. Add the 6 new solvents to all three validation dictionaries simultaneously to maintain consistency. Create a mapping layer in input_gen.py to translate user-friendly names (acetonitrile) to NWChem COSMO names (acetntrl) before input file generation.

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| Python | 3.11+ | Language runtime | Project standard (pyproject.toml requires-python) |
| pytest | 9.0.2+ | Testing framework | Used across entire test suite |
| argparse | stdlib | CLI parsing | Python standard library, already in use for benchmark CLI |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| pytest.mark.parametrize | - | Parameterized tests | Testing all solvents systematically |
| orjson | 3.11.5+ | JSON handling | Already used for scaling_factors.json |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| argparse choices | Custom validation | argparse choices provides automatic error messages and tab completion |
| Direct NWChem names | User-friendly names | User-friendly names improve DX but require mapping layer |

**Installation:**
No new dependencies required - all tools already in project.

## Architecture Patterns

### Recommended Project Structure
```
src/qm_nmr_calc/
├── solvents.py              # Validation layer (SUPPORTED_SOLVENTS dict)
├── shifts.py                # Scaling factor lookup (solvent_map dict)
├── nmredata.py              # NMReData export (solvent_map dict)
└── nwchem/
    └── input_gen.py         # NWChem input generation (SUPPORTED_SOLVENTS set + mapping)
```

### Pattern 1: Three-Module Solvent Gatekeeper

**What:** Solvents are validated in three separate modules, each with its own dictionary/set. Each module serves a different purpose but must stay synchronized.

**When to use:** When extending solvent support - all three must be updated together.

**Current implementation:**
```python
# solvents.py - Validation layer
SUPPORTED_SOLVENTS: dict[str, str] = {
    "chcl3": "Chloroform (CDCl3)",
    "dmso": "Dimethylsulfoxide (DMSO-d6)",
    # ... 7 total solvents
}

# shifts.py - Scaling factor lookup (Title-case normalization)
solvent_map = {
    "chcl3": "CHCl3",
    "dmso": "DMSO",
    # ... maps to scaling_factors.json keys
}

# nmredata.py - NMReData export (deuterated forms)
solvent_map = {
    "chcl3": "CDCl3",
    "dmso": "(CD3)2SO",
    # ... maps to NMReData solvent names
}

# nwchem/input_gen.py - NWChem input generation
SUPPORTED_SOLVENTS = {"chcl3", "dmso", "water", ...}
```

**Extension strategy for Phase 59:**
Add 6 new solvents to all three modules simultaneously:
- solvents.py: Add with descriptive names (e.g., "pyridine": "Pyridine")
- shifts.py: Add with Title-case (e.g., "pyridine": "Pyridine")
- nmredata.py: Add with deuterated forms (e.g., "pyridine": "C5D5N")
- input_gen.py: Add to set AND create mapping function for acetonitrile

### Pattern 2: COSMO Name Mapping Layer

**What:** NWChem COSMO uses specific internal names that don't always match user expectations. Most are intuitive (pyridine → pyridine), but acetonitrile requires special handling (acetonitrile → acetntrl).

**When to use:** In input_gen.py after validation, before writing to NWChem input file.

**Recommended implementation:**
```python
# nwchem/input_gen.py

# Map user-friendly names to NWChem COSMO names
COSMO_NAME_MAP = {
    "acetonitrile": "acetntrl",  # NWChem uses unusual abbreviation
    # All other solvents map to themselves
}

def _get_cosmo_solvent_name(solvent: str) -> str:
    """Map user-friendly solvent name to NWChem COSMO name.

    Most solvents use the same name, but some (acetonitrile) require mapping.
    """
    return COSMO_NAME_MAP.get(solvent, solvent)

def generate_optimization_input(..., solvent: str):
    solvent_name = _validate_solvent(solvent)  # Existing validation
    cosmo_name = _get_cosmo_solvent_name(solvent_name)  # NEW: mapping
    # Use cosmo_name in COSMO block
```

**Why this pattern:**
- Isolates NWChem naming quirks to one place
- Allows user-friendly names in CLI and API
- Easy to extend if more mappings are needed
- Explicit about the one special case

### Pattern 3: CLI Validation with argparse choices

**What:** Use argparse choices parameter to provide automatic validation and helpful error messages.

**Current implementation:**
```python
# benchmark/__main__.py
run_parser.add_argument(
    "--solvents",
    nargs="+",
    choices=["CHCl3", "DMSO", "Methanol", "Water", "Acetone", "Benzene"],
    help="Solvents to test (default: CHCl3 DMSO Methanol Water Acetone Benzene)",
)
```

**Extension for Phase 59:**
```python
run_parser.add_argument(
    "--solvents",
    nargs="+",
    choices=[
        "CHCl3", "DMSO", "Methanol", "Water", "Acetone", "Benzene",
        "Pyridine", "THF", "Toluene", "DCM", "Acetonitrile", "DMF"  # NEW
    ],
    help="Solvents to test (default: all 13 solvents)",
)
```

**Note:** CLI uses Title-case (Acetonitrile), but internal code uses lowercase (acetonitrile). The .lower() normalization happens in runner.py line 234.

### Pattern 4: Parametrized Test Coverage

**What:** Use pytest.mark.parametrize to test all solvents systematically without code duplication.

**Existing pattern:**
```python
# tests/test_nwchem_input.py
@pytest.mark.parametrize(
    "solvent",
    sorted(SUPPORTED_SOLVENTS - {"vacuum"}),
)
def test_all_supported_solvents_accepted(self, solvent):
    """All supported solvents should generate valid optimization input."""
    result = generate_optimization_input(
        geometry_xyz="C 0.0 0.0 0.0",
        functional="b3lyp",
        basis_set="6-31G*",
        solvent=solvent,
    )
    assert "cosmo" in result.lower()
    assert f"solvent {solvent}" in result
```

**Extension strategy:** Add new solvents to SUPPORTED_SOLVENTS and parametrized test automatically covers them. No test code changes needed beyond adding one test for the acetonitrile → acetntrl mapping.

### Anti-Patterns to Avoid

- **Inconsistent solvent lists:** If one module supports a solvent but another doesn't, subtle bugs occur. Always update all three modules together.
- **Hardcoded solvent lists:** The CLI hardcodes the list in choices. This is acceptable since the benchmark runner is explicitly tied to specific solvents (scaling factors must exist).
- **Case sensitivity bugs:** Always normalize to lowercase early in the pipeline (already done via .lower() in validation functions).

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| CLI argument validation | Custom string parsing with if/elif chains | argparse choices parameter | Automatic error messages, tab completion, self-documenting |
| Solvent name validation | Regex patterns or string matching | Dictionary membership test (in SUPPORTED_SOLVENTS) | Simpler, more maintainable, explicit list |
| Test data iteration | Multiple copy-paste test methods | pytest.mark.parametrize | Automatic coverage of all cases, less code duplication |

**Key insight:** The project already has validation patterns. This phase is purely additive - no new patterns needed.

## Common Pitfalls

### Pitfall 1: Forgetting to Update All Three Modules

**What goes wrong:** Adding a solvent to input_gen.py but forgetting shifts.py or nmredata.py. NWChem calculations succeed but shift conversion or export fails with cryptic errors.

**Why it happens:** The three modules are in different files and serve different purposes. Easy to miss one.

**How to avoid:**
- Create a checklist: solvents.py ✓, shifts.py ✓, nmredata.py ✓, input_gen.py ✓
- Add all solvents in one commit to ensure atomicity
- Run full test suite (pytest will catch missing mappings)

**Warning signs:**
- ValueError during shifts.py conversion: "No scaling factor for X"
- ValueError during nmredata.py export: "Unknown solvent"
- Tests pass for input generation but fail for shifts/export

### Pitfall 2: Using Wrong COSMO Name (acetonitrile vs acetntrl)

**What goes wrong:** Using "acetonitrile" in the COSMO block instead of "acetntrl". NWChem will either fail to recognize the solvent or use wrong parameters.

**Why it happens:** NWChem uses an unusual abbreviation "acetntrl" (not "acetonitrile" or "acn"). This is documented in NWChem COSMO tables but easy to miss.

**How to avoid:**
- Verify NWChem COSMO name in official documentation before adding solvent
- Create explicit mapping function (_get_cosmo_solvent_name) for the translation
- Add integration test that runs actual NWChem calculation with acetonitrile

**Warning signs:**
- NWChem error: "Unknown solvent"
- NWChem uses default dielectric constant instead of 35.688

### Pitfall 3: Case Sensitivity in Solvent Names

**What goes wrong:** CLI accepts "Pyridine" but code expects "pyridine". Validation fails or lookups return None.

**Why it happens:** Different layers use different conventions (CLI uses Title-case for readability, internal uses lowercase for NWChem compatibility).

**How to avoid:**
- Always normalize to lowercase in validation functions (already done via .lower().strip())
- Keep internal code lowercase-only
- Use Title-case only in CLI help text and choices

**Warning signs:**
- ValueError: "Unknown solvent 'Pyridine'" even though it's in the list
- Dictionary lookup returns None for valid solvent

### Pitfall 4: Missing Scaling Factors

**What goes wrong:** Adding solvent to validation but forgetting it won't have scaling factors yet. Calculations run but shift conversion fails.

**Why it happens:** Phase 59 adds solvents to CLI but doesn't run benchmark or derive factors. Those are future phases.

**How to avoid:**
- Document in code comments: "Scaling factors for these solvents will be derived in Phase 60"
- Handle missing scaling factors gracefully in runner.py (already done - see lines 242-256)
- Don't add new solvents to SOLVENTS constant in runner.py (line 24) yet - they're not ready for default runs

**Warning signs:**
- ValueError: "No scaling factor for B3LYP/6-311+G(2d,p)/1H/Pyridine"
- This is expected until factors are derived

## Code Examples

Verified patterns from existing codebase:

### Adding Solvent to solvents.py
```python
# Source: /home/chris/develop/qm-nmr-calc/src/qm_nmr_calc/solvents.py
SUPPORTED_SOLVENTS: dict[str, str] = {
    "chcl3": "Chloroform (CDCl3)",
    "dmso": "Dimethylsulfoxide (DMSO-d6)",
    "vacuum": "Gas phase (no solvent)",
    "methanol": "Methanol (Methanol-d4)",
    "water": "Water (D2O)",
    "acetone": "Acetone (Acetone-d6)",
    "benzene": "Benzene (Benzene-d6)",
    # NEW: Add these 6 solvents
    "pyridine": "Pyridine (Pyridine-d5)",
    "thf": "Tetrahydrofuran (THF-d8)",
    "toluene": "Toluene (Toluene-d8)",
    "dcm": "Dichloromethane (DCM-d2)",
    "acetonitrile": "Acetonitrile (CD3CN)",
    "dmf": "N,N-Dimethylformamide (DMF-d7)",
}
```

### Adding Solvent to shifts.py
```python
# Source: /home/chris/develop/qm-nmr-calc/src/qm_nmr_calc/shifts.py
solvent_map = {
    "chcl3": "CHCl3",
    "dmso": "DMSO",
    "vacuum": "vacuum",
    "methanol": "Methanol",
    "water": "Water",
    "acetone": "Acetone",
    "benzene": "Benzene",
    # NEW: Add these 6 solvents (Title-case for scaling_factors.json keys)
    "pyridine": "Pyridine",
    "thf": "THF",
    "toluene": "Toluene",
    "dcm": "DCM",
    "acetonitrile": "Acetonitrile",
    "dmf": "DMF",
}
```

### Adding Solvent to nmredata.py
```python
# Source: /home/chris/develop/qm-nmr-calc/src/qm_nmr_calc/nmredata.py
solvent_map = {
    "chcl3": "CDCl3",
    "dmso": "(CD3)2SO",
    "vacuum": "vacuum",
    "methanol": "CD3OD",
    "water": "D2O",
    "acetone": "(CD3)2CO",
    "benzene": "C6D6",
    # NEW: Add these 6 solvents (deuterated forms for NMReData)
    "pyridine": "C5D5N",  # Pyridine-d5
    "thf": "C4D8O",  # THF-d8
    "toluene": "C7D8",  # Toluene-d8
    "dcm": "CD2Cl2",  # DCM-d2
    "acetonitrile": "CD3CN",  # Acetonitrile-d3
    "dmf": "C3D7NO",  # DMF-d7
}
```

### Adding Solvent to input_gen.py with COSMO Mapping
```python
# Source: /home/chris/develop/qm-nmr-calc/src/qm_nmr_calc/nwchem/input_gen.py
SUPPORTED_SOLVENTS = {
    "chcl3", "dmso", "water", "acetone", "methanol", "benzene", "vacuum",
    # NEW: Add these 6 solvents
    "pyridine", "thf", "toluene", "dcm", "acetonitrile", "dmf"
}

# NEW: COSMO name mapping for special cases
COSMO_NAME_MAP = {
    "acetonitrile": "acetntrl",  # NWChem uses this abbreviation
}

def _get_cosmo_solvent_name(solvent: str) -> str:
    """Map user-friendly solvent name to NWChem COSMO name.

    Most solvents use the same name, but acetonitrile maps to 'acetntrl'.

    Args:
        solvent: Normalized solvent name (lowercase)

    Returns:
        NWChem COSMO solvent name
    """
    return COSMO_NAME_MAP.get(solvent, solvent)

def generate_optimization_input(...):
    solvent_name = _validate_solvent(solvent)

    if solvent_name == "vacuum":
        cosmo_block = ""
    else:
        cosmo_name = _get_cosmo_solvent_name(solvent_name)  # NEW
        cosmo_block = f"""
cosmo
  do_gasphase False
  solvent {cosmo_name}
end
"""
    # ... rest unchanged
```

### Test for COSMO Name Mapping
```python
# NEW test in tests/test_nwchem_input.py
class TestCOSMONameMapping:
    """Tests for NWChem COSMO solvent name mapping."""

    def test_acetonitrile_maps_to_acetntrl(self):
        """Acetonitrile should use 'acetntrl' in COSMO block."""
        result = generate_optimization_input(
            geometry_xyz="C 0.0 0.0 0.0",
            functional="b3lyp",
            basis_set="6-31G*",
            solvent="acetonitrile",
        )

        assert "cosmo" in result.lower()
        assert "solvent acetntrl" in result  # NWChem name
        assert "acetonitrile" not in result.lower()  # User name not in output

    def test_other_solvents_use_direct_names(self):
        """Most solvents use their name directly (no mapping)."""
        for solvent in ["pyridine", "thf", "toluene", "dcm", "dmf"]:
            result = generate_optimization_input(
                geometry_xyz="C 0.0 0.0 0.0",
                functional="b3lyp",
                basis_set="6-31G*",
                solvent=solvent,
            )
            assert f"solvent {solvent}" in result
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| 7 solvents (v2.8) | 13 solvents (v2.9 after Phase 59) | Phase 59 | Extended solvent coverage for broader applicability |
| Direct NWChem names in user input | User-friendly names + mapping layer | Phase 59 | Better UX, hides NWChem naming quirks |

**Note:** This is an incremental extension, not a paradigm shift. The architecture remains the same.

## Open Questions

1. **NMReData deuterated forms**
   - What we know: Standard forms exist (C5D5N for pyridine-d5, etc.)
   - What's unclear: Whether these match NMReData convention exactly
   - Recommendation: Verify with NMReData specification or check existing NMR database exports. Fall back to chemical formulas if unclear.

2. **Scaling factors for new solvents**
   - What we know: Phase 59 adds CLI support, but scaling factors don't exist yet
   - What's unclear: Whether to add solvents to default SOLVENTS list in runner.py now
   - Recommendation: Add to SUPPORTED_SOLVENTS but NOT to default SOLVENTS constant. Users can opt-in with --solvents flag. Default list gets updated after Phase 60 derives factors.

3. **Integration test feasibility**
   - What we know: Full NWChem calculations are expensive
   - What's unclear: Whether to add integration tests for new solvents in Phase 59
   - Recommendation: Add unit tests (input generation, name mapping) now. Integration tests come later when running benchmark calculations.

## Sources

### Primary (HIGH confidence)
- [NWChem Solvation Models Documentation](https://nwchemgit.github.io/Solvation-Models.html) - Official COSMO solvent list with dielectric constants
- [NWChem COSMO Solvent Discussion](https://groups.google.com/g/nwchem-forum/c/GST-BYQc8Mg/m/pFtlb9H-BgAJ) - Confirmation of acetonitrile → acetntrl mapping
- Project source code:
  - /home/chris/develop/qm-nmr-calc/src/qm_nmr_calc/solvents.py
  - /home/chris/develop/qm-nmr-calc/src/qm_nmr_calc/shifts.py
  - /home/chris/develop/qm-nmr-calc/src/qm_nmr_calc/nmredata.py
  - /home/chris/develop/qm-nmr-calc/src/qm_nmr_calc/nwchem/input_gen.py
  - /home/chris/develop/qm-nmr-calc/src/qm_nmr_calc/benchmark/__main__.py
  - /home/chris/develop/qm-nmr-calc/tests/test_nwchem_input.py
  - /home/chris/develop/qm-nmr-calc/tests/test_nmredata.py

### Secondary (MEDIUM confidence)
- NWChem COSMO dielectric constants from documentation (verified values):
  - pyridine: ε=12.978
  - THF: ε=7.4257
  - toluene: ε=2.3741
  - DCM: ε=8.930
  - acetonitrile: ε=35.688
  - DMF: ε=37.219

### Tertiary (LOW confidence)
None - all findings verified with official documentation or project source code.

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - All tools already in project, no new dependencies
- Architecture: HIGH - Following established 3-module pattern exactly
- Pitfalls: HIGH - Identified from code review and test suite analysis

**Research date:** 2026-02-10
**Valid until:** 2026-04-10 (60 days - stable domain, NWChem COSMO list rarely changes)
