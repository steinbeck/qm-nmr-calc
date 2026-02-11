# Phase 63: Scaling Factor Derivation - Research

**Researched:** 2026-02-11
**Domain:** Linear regression for NMR scaling factor derivation
**Confidence:** HIGH

## Summary

Phase 63 derives OLS scaling factors for 6 new solvents (Pyridine, THF, Toluene, DCM, Acetonitrile, DMF) using the EXACT same tooling and patterns established in Phase 56 for the previous 4 solvents (Methanol, Water, Acetone, Benzene). The analysis module is already production-ready and requires only a simple default solvents list update to process the new benchmark data.

**Key finding:** All 300 benchmark calculations (50 compounds × 6 solvents) are complete and verified in `data/benchmark/results/`. The existing `python -m qm_nmr_calc.benchmark analyze` command will derive the 12 new factor sets when the default solvents list in `analysis.py` is updated to include the 6 new solvents.

**Primary recommendation:** Update `derive_all_factors()` default solvents parameter from the current 6 solvents to all 12 solvents, run analysis, merge 12 new factors with existing 14 to get 26 total.

## Standard Stack

The stack is unchanged from Phase 56 - all tools are already in place:

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| statsmodels | 0.14+ | OLS regression with confidence intervals | Industry standard for statistical regression in Python |
| pandas | 2.0+ | Data aggregation and grouping | Standard for tabular data manipulation |
| numpy | 1.24+ | Numerical operations and outlier detection | Foundation for scientific computing |
| matplotlib | 3.7+ | Regression and residual plots | Standard visualization for scientific plots |
| orjson | 3.9+ | Fast JSON serialization | Used for package data and analysis output |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| pathlib | stdlib | File path handling | All file operations |

**Installation:**
All dependencies already installed in project environment. No new packages needed.

## Architecture Patterns

### Project Structure (Existing)
```
src/qm_nmr_calc/benchmark/
├── analysis.py           # OLS regression and report generation
├── __main__.py          # CLI entry point (analyze command)
├── models.py            # ScalingFactor pydantic model
├── data_loader.py       # Load experimental shifts
└── runner.py            # Benchmark presets and results dir

data/benchmark/
├── results/             # Raw calculation results (300 shifts.json files)
│   └── compound_XX/B3LYP_{Solvent}/shifts.json
└── delta50/             # Analysis output directory
    ├── scaling_factors.json
    ├── SCALING-FACTORS.md
    └── plots/
```

### Pattern 1: Default Solvents Update
**What:** Update the `derive_all_factors()` function's default solvents parameter to include all 12 solvents
**When to use:** When new solvent benchmark data becomes available
**Example:**
```python
# File: src/qm_nmr_calc/benchmark/analysis.py
# Line: ~213 (in derive_all_factors function signature)

def derive_all_factors(
    functionals: list[str] | None = None,
    solvents: list[str] | None = None,
) -> dict[str, ScalingFactor]:
    """Derive scaling factors for all functional/solvent/nucleus combinations."""
    if functionals is None:
        functionals = ["B3LYP"]  # WP04 not yet complete
    if solvents is None:
        # BEFORE (Phase 56): 6 solvents
        # solvents = ["CHCl3", "DMSO", "Methanol", "Water", "Acetone", "Benzene"]

        # AFTER (Phase 63): 12 solvents (add 6 new)
        solvents = [
            "CHCl3", "DMSO", "Methanol", "Water", "Acetone", "Benzene",
            "Pyridine", "THF", "Toluene", "DCM", "Acetonitrile", "DMF"
        ]
```

### Pattern 2: Analysis Execution
**What:** Run CLI analyze command to derive factors and generate report
**When to use:** After updating default solvents list
**Example:**
```bash
# Derive all factors and generate full report
python -m qm_nmr_calc.benchmark analyze --output-dir data/benchmark/delta50

# Or just print factors without plots (for testing)
python -m qm_nmr_calc.benchmark analyze --factors-only
```

### Pattern 3: Factor Merging Strategy
**What:** Merge newly derived factors with existing factors, preserving vacuum entries
**When to use:** After analysis completes successfully
**Merge logic:**
```python
import orjson
from pathlib import Path

# Load newly derived factors (24 entries: 12 existing + 12 new)
new_factors = orjson.loads(Path("data/benchmark/delta50/scaling_factors.json").read_bytes())

# Load existing package factors to get vacuum entries (2 entries)
pkg_factors = orjson.loads(Path("src/qm_nmr_calc/data/scaling_factors.json").read_bytes())

# Add vacuum entries (not in benchmark results, from Phase 11.2)
for key, value in pkg_factors.items():
    if "vacuum" in key:
        new_factors[key] = value

# Should have 26 entries total (24 solvent + 2 vacuum)
assert len(new_factors) == 26, f"Expected 26, got {len(new_factors)}"

# Write merged file
Path("src/qm_nmr_calc/data/scaling_factors.json").write_bytes(
    orjson.dumps(new_factors, option=orjson.OPT_INDENT_2)
)
```

## Don't Hand-Roll

Problems with existing solutions from Phase 56:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| OLS regression with CI | Custom least squares | statsmodels.api.OLS | Handles confidence intervals, diagnostics correctly |
| Outlier detection | Manual threshold logic | 3-sigma rule with refit | Standard statistical approach, proven in Phase 56 |
| Regression plots | Custom matplotlib | Existing plot_regression() | Already handles scatter, line, residuals, outliers |
| Factor serialization | Custom JSON logic | Existing save_factors_json() | Handles pydantic model conversion correctly |
| Report generation | Custom markdown builder | Existing generate_report() | Comprehensive, proven format |

**Key insight:** ALL tooling exists and is production-tested from Phase 56. This phase is a pure data expansion, not a code development phase.

## Common Pitfalls

### Pitfall 1: Solvent Name Case Sensitivity
**What goes wrong:** Analysis fails to find results because solvent name casing doesn't match directory names
**Why it happens:** Benchmark results directories use exact casing (e.g., `B3LYP_Pyridine`, not `B3LYP_pyridine`)
**How to avoid:** Use exact casing from benchmark results directories in the solvents list
**Warning signs:** Empty DataFrame from `aggregate_regression_data()`, warnings about "No data for {functional}/{solvent}"
**Verified safe names:**
- `Pyridine` (not pyridine)
- `THF` (all caps)
- `Toluene` (capital T)
- `DCM` (all caps)
- `Acetonitrile` (capital A)
- `DMF` (all caps)

### Pitfall 2: Acetonitrile COSMO Name Confusion
**What goes wrong:** Developer confuses user-facing name "Acetonitrile" with COSMO internal name "acetntrl"
**Why it happens:** NWChem uses abbreviated "acetntrl" in COSMO parameter tables, but results are stored as "Acetonitrile"
**How to avoid:** Analysis uses directory names ("Acetonitrile"), not COSMO names. The COSMO mapping is only in `input_gen.py` for calculation time.
**Warning signs:** Looking for "acetntrl" in results directories (doesn't exist)
**Verification:** `ls data/benchmark/results/compound_01/` shows `B3LYP_Acetonitrile`, not `B3LYP_acetntrl`

### Pitfall 3: Vacuum Factor Overwrite
**What goes wrong:** Vacuum factors get lost during merge because they're not in the newly derived factors
**Why it happens:** Vacuum factors came from Phase 11.2's separate vacuum benchmark, not from the DELTA50 solvent benchmark
**How to avoid:** Explicitly preserve vacuum entries from existing package file during merge (see Pattern 3)
**Warning signs:** Final package file has 24 entries instead of 26
**Verification:** `jq 'keys[] | select(contains("vacuum"))' src/qm_nmr_calc/data/scaling_factors.json` returns 2 entries

### Pitfall 4: Forgetting to Update shifts.py Solvent Map
**What goes wrong:** Runtime `get_scaling_factor()` fails for new solvents even though factors exist
**Why it happens:** `shifts.py` has a `solvent_map` dict for normalization that doesn't include new solvents
**How to avoid:** After merging factors, update the `solvent_map` in `shifts.py` to include lowercase mappings for new solvents
**Warning signs:** `ValueError: No scaling factor for B3LYP/6-311+G(2d,p)/1H/Pyridine`
**Fix required:** Add entries to `solvent_map` in `shifts.py` (deferred to Phase 64 per ROADMAP)

## Code Examples

### Example 1: Verify Benchmark Data Exists
```bash
# Count results per new solvent (should be 50 each)
for solvent in Pyridine THF Toluene DCM Acetonitrile DMF; do
    count=$(find data/benchmark/results/ -name "shifts.json" -path "*B3LYP_${solvent}*" | wc -l)
    echo "$solvent: $count"
done

# Expected output:
# Pyridine: 50
# THF: 50
# Toluene: 50
# DCM: 50
# Acetonitrile: 50
# DMF: 50
```

### Example 2: Test Analysis with --factors-only
```bash
# Quick test without generating plots (fast feedback)
python -m qm_nmr_calc.benchmark analyze --factors-only

# Should print 24 factor sets (12 existing + 12 new)
# If it prints 12, the solvents list wasn't updated yet
```

### Example 3: Validate Quality Gates
```python
import orjson
from pathlib import Path

factors = orjson.loads(Path("data/benchmark/delta50/scaling_factors.json").read_bytes())

# Check R² > 0.99 for all factors
failures = []
for key, factor in factors.items():
    if factor["r_squared"] < 0.99:
        failures.append(f"{key}: R²={factor['r_squared']:.4f}")

if failures:
    print("R² GATE FAILURES:")
    for f in failures:
        print(f"  {f}")
else:
    print("✓ All factors pass R² > 0.99 gate")

# Check 1H MAE < 0.2 ppm
h1_failures = []
for key, factor in factors.items():
    if "1H" in key and factor["mae"] >= 0.2:
        h1_failures.append(f"{key}: MAE={factor['mae']:.3f}")

if h1_failures:
    print("1H MAE GATE FAILURES:")
    for f in h1_failures:
        print(f"  {f}")
else:
    print("✓ All 1H factors pass MAE < 0.2 ppm gate")

# Check 13C MAE < 3.0 ppm
c13_failures = []
for key, factor in factors.items():
    if "13C" in key and factor["mae"] >= 3.0:
        c13_failures.append(f"{key}: MAE={factor['mae']:.3f}")

if c13_failures:
    print("13C MAE GATE FAILURES:")
    for f in c13_failures:
        print(f"  {f}")
else:
    print("✓ All 13C factors pass MAE < 3.0 ppm gate")
```

### Example 4: Verify Final Package Data
```python
from qm_nmr_calc.shifts import load_scaling_factors

factors = load_scaling_factors()

print(f"Total factors: {len(factors)}")  # Should be 26

# Count by solvent
solvents = set()
for key in factors.keys():
    solvent = key.split("/")[-1]
    solvents.add(solvent)

print(f"Unique solvents: {sorted(solvents)}")
# Should include: Acetone, Acetonitrile, Benzene, CHCl3, DCM, DMF, DMSO,
#                 Methanol, Pyridine, THF, Toluene, Water, vacuum

# Verify each solvent has both 1H and 13C
for solvent in sorted(solvents):
    h1_key = f"B3LYP/6-311+G(2d,p)/1H/{solvent}"
    c13_key = f"B3LYP/6-311+G(2d,p)/13C/{solvent}"
    h1_present = h1_key in factors
    c13_present = c13_key in factors
    print(f"{solvent:15} 1H:{h1_present}  13C:{c13_present}")
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Manual factor derivation | CLI analyze command | Phase 8-11 (v1.1) | Automated, reproducible |
| No outlier removal | 3-sigma outlier removal with refit | Phase 11 (v1.1) | Better R² and MAE |
| 2 solvents (CHCl3, DMSO) | 14 solvents after Phase 63 | Phase 56 (v2.8), Phase 63 (v2.9) | Broader applicability |
| Vacuum factors only | Solvent-specific factors | Phase 11 (v1.1) | Accurate solvent-dependent shifts |

**No deprecated patterns:** Phase 56 established the current methodology, which remains state-of-the-art for this project.

## Open Questions

**None.** All patterns are proven from Phase 56. The only unknown is whether the new solvents will pass quality gates, but this will be determined during execution.

## Sources

### Primary (HIGH confidence)
- `/home/chris/develop/qm-nmr-calc/src/qm_nmr_calc/benchmark/analysis.py` - Lines 197-243 (`derive_all_factors` function with current default solvents list)
- `/home/chris/develop/qm-nmr-calc/src/qm_nmr_calc/benchmark/__main__.py` - Lines 201-221 (`cmd_analyze` function showing CLI usage)
- `/home/chris/develop/qm-nmr-calc/src/qm_nmr_calc/data/scaling_factors.json` - Current package data with 14 entries (verified via `jq`)
- `.planning/phases/56-scaling-factor-derivation/56-01-PLAN.md` - Complete reference implementation
- `.planning/phases/56-scaling-factor-derivation/56-01-SUMMARY.md` - Verified outcomes and quality metrics

### Secondary (HIGH confidence)
- `/home/chris/develop/qm-nmr-calc/data/benchmark/results/` - 300 verified shifts.json files (50 × 6 solvents)
- `/home/chris/develop/qm-nmr-calc/src/qm_nmr_calc/shifts.py` - Lines 49-59 (solvent_map for runtime loading)
- `/home/chris/develop/qm-nmr-calc/src/qm_nmr_calc/nwchem/input_gen.py` - Lines 1-30 (COSMO_NAME_MAP showing acetonitrile→acetntrl)

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - All tools verified installed and working from Phase 56
- Architecture: HIGH - Patterns proven in production, no new code needed
- Pitfalls: HIGH - All identified from Phase 56 execution and current codebase inspection

**Research date:** 2026-02-11
**Valid until:** 2026-03-11 (30 days - stable domain, proven tooling)

**Benchmark data verified:**
- Pyridine: 50 calculations ✓
- THF: 50 calculations ✓
- Toluene: 50 calculations ✓
- DCM: 50 calculations ✓
- Acetonitrile: 50 calculations ✓
- DMF: 50 calculations ✓
- **Total: 300 calculations ready for analysis**

**Key insight:** This is a data expansion phase, not a tooling development phase. The entire infrastructure exists and is production-tested. Execution should be straightforward: update solvents list, run analysis, merge factors, validate.
