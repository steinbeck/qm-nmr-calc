# Phase 11: Production Integration - Research

**Researched:** 2026-01-23
**Domain:** Python production refactoring - replacing TMS references with regression factors and removing legacy dependencies
**Confidence:** HIGH

## Summary

This phase integrates NWChem-derived scaling factors from Phase 10 into production calculations, replacing the legacy TMS reference approach from CHESHIRE. The work involves three main components: (1) loading scaling factors from JSON at runtime using modern Python package data patterns, (2) removing the ISiCLE dependency completely from the codebase, and (3) updating API responses to include factor metadata for transparency and reproducibility.

The standard approach for this type of production integration in Python scientific packages is to bundle small data files (< 100 kB) within the package using `importlib.resources` for reliable access across different installation types. The existing scaling_factors.json file (1.7 KB) fits well within this threshold. For dependency removal, the key is ensuring clean breaks - remove all imports, update version tracking, and handle backward compatibility explicitly rather than assuming graceful degradation.

The current codebase already has "N/A" placeholders for isicle_version in the API endpoints, indicating Phase 7 partially addressed ISiCLE removal. This phase completes the cleanup by removing the model field entirely and eliminating the pyproject.toml dependency. The shifts.py module currently uses hardcoded SCALING_FACTORS dictionaries with CHESHIRE data; these will be replaced with runtime JSON loading from the Phase 10 derived factors.

**Primary recommendation:** Use importlib.resources.files() to load scaling_factors.json from a new data/ subdirectory within the package. Build a factor lookup service that fails fast with clear error messages when no factor exists for a given functional/solvent/nucleus combination. Remove isicle_version field from models using Pydantic schema evolution patterns.

## Standard Stack

The established libraries/tools for this domain:

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| importlib.resources | stdlib (3.11+) | Package data access | Official Python stdlib for resource loading, no performance overhead vs pkg_resources |
| orjson | 3.11.5 | JSON parsing | Already in dependencies, 33% faster than stdlib json (280ms vs 420ms), returns bytes for throughput |
| Pydantic | 2.12.5 | Schema evolution | Already in use, supports model_config for field deprecation and migration |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| hatchling | current | Build backend | Already configured, handles package_data via [tool.hatch.build] |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| importlib.resources | pathlib.Path(__file__).parent | __file__ breaks in zip distributions, not portable |
| orjson | stdlib json | orjson already in use, no reason to switch |
| Pooch | Manual download | Overkill for 1.7 KB static data |

**Installation:**
No new dependencies required - using only stdlib and existing packages.

## Architecture Patterns

### Recommended Project Structure
```
src/qm_nmr_calc/
├── data/
│   └── scaling_factors.json    # Moved from data/benchmark/delta50/
├── shifts.py                   # Updated to load JSON
├── models.py                   # Remove isicle_version field
└── api/
    └── schemas.py              # Update response schemas
```

### Pattern 1: Package Data Loading with importlib.resources
**What:** Load JSON data from package using stdlib resource API
**When to use:** Any time loading static data files bundled with package
**Example:**
```python
# Source: https://docs.python.org/3/library/importlib.resources.html
from importlib.resources import files
import orjson

def load_scaling_factors() -> dict:
    """Load scaling factors from package data."""
    data_dir = files("qm_nmr_calc").joinpath("data")
    json_bytes = data_dir.joinpath("scaling_factors.json").read_bytes()
    return orjson.loads(json_bytes)
```

### Pattern 2: Factor Lookup with Composite Keys
**What:** Build lookup key from calculation parameters (functional/basis/nucleus/solvent)
**When to use:** When factors depend on multiple parameters
**Example:**
```python
# Source: Project CONTEXT.md decisions
def get_scaling_factor(functional: str, basis_set: str, nucleus: str, solvent: str) -> dict:
    """Get scaling factor for specific calculation parameters.

    Args:
        functional: DFT functional (e.g., 'B3LYP')
        basis_set: Basis set (e.g., '6-311+G(2d,p)')
        nucleus: Nucleus type ('1H' or '13C')
        solvent: Solvent name (e.g., 'CHCl3', 'DMSO')

    Returns:
        Dict with keys: slope, intercept, r_squared, mae, rmsd, n_points

    Raises:
        ValueError: If no factor exists for this combination
    """
    key = f"{functional}/{basis_set}/{nucleus}/{solvent}"
    factors = load_scaling_factors()  # Consider caching this

    if key not in factors:
        raise ValueError(
            f"No scaling factor available for {key}. "
            f"Available combinations: {list(factors.keys())}"
        )

    return factors[key]
```

### Pattern 3: Schema Evolution with Pydantic
**What:** Remove deprecated fields from models while maintaining API compatibility
**When to use:** When removing fields from Pydantic models
**Example:**
```python
# Source: https://docs.pydantic.dev/latest/concepts/models/
# OLD (keep for reference):
# class JobStatus(BaseModel):
#     isicle_version: str
#     nwchem_version: str

# NEW:
class JobStatus(BaseModel):
    nwchem_version: str
    # isicle_version removed - ISiCLE no longer used
```

### Pattern 4: Factor Metadata in API Responses
**What:** Include provenance metadata with calculation results
**When to use:** Scientific APIs where reproducibility matters
**Example:**
```python
# Source: Phase 11 CONTEXT.md decisions
class NMRResultsResponse(BaseModel):
    h1_shifts: list[AtomShiftResponse]
    c13_shifts: list[AtomShiftResponse]
    functional: str
    basis_set: str
    solvent: str
    # New fields:
    scaling_factor_source: str = "DELTA50"  # Factor source
    h1_mae: float  # Expected MAE for 1H in ppm
    c13_mae: float  # Expected MAE for 13C in ppm
```

### Anti-Patterns to Avoid
- **Silent fallbacks:** Don't fall back to TMS if factor missing - fail explicitly so users know their results are unreliable
- **Global imports:** Don't load JSON at module level - lazy load on first use to avoid import-time I/O
- **Ignoring errors:** Don't catch FileNotFoundError and continue - if scaling_factors.json is missing, the package is broken
- **Version field removal without cleanup:** Don't just stop using isicle_version - remove it from models, API schemas, and all call sites

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| JSON parsing | Custom parser | orjson (already in deps) | Handles edge cases, 33% faster, well-tested |
| Package data paths | __file__ + relative paths | importlib.resources | Works in zip files, handles namespace packages |
| Schema migration | Manual field filtering | Pydantic model_config | Type-safe, automatic validation |
| Error messages | Generic ValueError | Specific message with available options | User sees what went wrong and how to fix it |

**Key insight:** Python stdlib and existing dependencies already handle these concerns. Focus on integration logic, not reimplementing resource loading or JSON parsing.

## Common Pitfalls

### Pitfall 1: Loading Data at Module Import Time
**What goes wrong:** Loading JSON in global scope causes import-time I/O, slowing imports and making errors harder to debug
**Why it happens:** Convenient to have SCALING_FACTORS as module constant like old code
**How to avoid:** Lazy load on first use with @functools.cache decorator
**Warning signs:** "ModuleNotFoundError during import" or slow import times
**Example:**
```python
# BAD - loads at import time
SCALING_FACTORS = orjson.loads(
    files("qm_nmr_calc").joinpath("data/scaling_factors.json").read_bytes()
)

# GOOD - lazy load with cache
from functools import cache

@cache
def load_scaling_factors() -> dict:
    """Load scaling factors from package data (cached after first call)."""
    data_dir = files("qm_nmr_calc").joinpath("data")
    json_bytes = data_dir.joinpath("scaling_factors.json").read_bytes()
    return orjson.loads(json_bytes)
```

### Pitfall 2: Graceful Degradation for Missing Factors
**What goes wrong:** Code falls back to TMS reference or skips factor application when factor is missing, producing inaccurate results
**Why it happens:** Developer wants to "be robust" and not crash on missing data
**How to avoid:** Fail fast with clear error message listing available combinations
**Warning signs:** Results differ from Phase 10 analysis, users confused why accuracy varies
**Example:**
```python
# BAD - silent fallback
try:
    factor = get_scaling_factor(functional, basis_set, nucleus, solvent)
except KeyError:
    factor = TMS_REFERENCE_FALLBACK  # Wrong! Produces inaccurate shifts

# GOOD - explicit failure
factor = get_scaling_factor(functional, basis_set, nucleus, solvent)
# Let ValueError propagate with clear message about what's missing
```

### Pitfall 3: Incomplete Dependency Removal
**What goes wrong:** ISiCLE removed from pyproject.toml but still imported in scripts/ or __pycache__ references remain
**Why it happens:** Forgetting to check non-package code (scripts/) and compiled bytecode
**How to avoid:** Systematic search for all references: grep -r "isicle" including scripts/ and exclude __pycache__ results
**Warning signs:** Import errors when running scripts, confusing error messages about missing isicle
**Example:**
```bash
# Check for all ISiCLE references
grep -r "isicle" src/ scripts/ --include="*.py" | grep -v __pycache__

# Remove ISiCLE from pyproject.toml dependencies
# Update scripts/ to not import from isicle_wrapper
```

### Pitfall 4: Forgetting to Configure Package Data
**What goes wrong:** scaling_factors.json not included in wheel distribution, works in development but fails after pip install
**Why it happens:** hatchling doesn't include non-Python files by default unless configured
**How to avoid:** Add package data configuration to pyproject.toml and verify with built wheel
**Warning signs:** "FileNotFoundError: scaling_factors.json" only in installed package, not in development
**Example:**
```toml
# pyproject.toml
[tool.hatch.build.targets.wheel]
packages = ["src/qm_nmr_calc"]
# If using include_package_data pattern, ensure MANIFEST.in exists
# or explicitly list in artifacts

[tool.hatch.build.targets.wheel.force-include]
"data/benchmark/delta50/scaling_factors.json" = "qm_nmr_calc/data/scaling_factors.json"
```

### Pitfall 5: Breaking Existing Job Status Files
**What goes wrong:** Old job status.json files have isicle_version field, new code crashes when loading them
**Why it happens:** Pydantic strict mode requires all fields in JSON match model
**How to avoid:** Set model_config strict=False temporarily or use model_validate with from_attributes
**Warning signs:** "ValidationError: extra fields not permitted" when loading old jobs
**Example:**
```python
# models.py - handle legacy data
class JobStatus(BaseModel):
    model_config = ConfigDict(
        strict=False,  # Allow extra fields from old JSON
        # Or use extra="ignore" in Pydantic v2
    )

    nwchem_version: str
    # isicle_version removed - old files may still have it, that's OK
```

## Code Examples

Verified patterns from official sources:

### Loading Scaling Factors (importlib.resources)
```python
# Source: https://docs.python.org/3/library/importlib.resources.html
from functools import cache
from importlib.resources import files
import orjson

@cache
def load_scaling_factors() -> dict:
    """Load scaling factors from package data.

    Returns:
        Dict mapping factor keys to dicts with slope, intercept, r_squared, etc.
        Factor key format: "B3LYP/6-311+G(2d,p)/1H/CHCl3"
    """
    data_dir = files("qm_nmr_calc").joinpath("data")
    json_bytes = data_dir.joinpath("scaling_factors.json").read_bytes()
    return orjson.loads(json_bytes)
```

### Factor Lookup with Clear Errors
```python
# Source: Phase 11 CONTEXT.md and error handling best practices
def get_scaling_factor(
    functional: str,
    basis_set: str,
    nucleus: str,
    solvent: str
) -> dict:
    """Get scaling factor for specific calculation parameters.

    Raises:
        ValueError: If no factor exists for this combination
    """
    key = f"{functional}/{basis_set}/{nucleus}/{solvent}"
    factors = load_scaling_factors()

    if key not in factors:
        available = [k for k in factors.keys() if k.split('/')[-1] == solvent]
        raise ValueError(
            f"No scaling factor for {key}.\n"
            f"This solvent/functional combination is not supported.\n"
            f"Available factors for {solvent}: {available}"
        )

    return factors[key]
```

### Applying Regression Factors
```python
# Source: Phase 10 CONTEXT.md - regression approach
def shielding_to_shift(
    shielding_data: dict,
    functional: str,
    basis_set: str,
    solvent: str,
) -> dict[str, list[dict]]:
    """Convert shielding values to chemical shifts using regression factors.

    Uses linear regression: shift = slope * shielding + intercept
    Factors are derived from DELTA50 benchmark data (Phase 10).

    Args:
        shielding_data: Dict with 'index', 'atom', 'shielding' keys
        functional: DFT functional used (e.g., 'B3LYP')
        basis_set: NMR basis set used (e.g., '6-311+G(2d,p)')
        solvent: Solvent used (e.g., 'CHCl3', 'DMSO')

    Returns:
        Dict with '1H' and '13C' keys, each containing list of shift dicts
    """
    shifts = []

    for idx, atom, shield in zip(
        shielding_data["index"],
        shielding_data["atom"],
        shielding_data["shielding"],
    ):
        # Map atom symbol to nucleus notation
        nucleus = "1H" if atom == "H" else "13C"

        # Get factor for this nucleus/solvent
        factor = get_scaling_factor(functional, basis_set, nucleus, solvent)

        # Apply linear regression: shift = m * shielding + b
        shift = factor["slope"] * shield + factor["intercept"]

        shifts.append({
            "index": idx,
            "atom": atom,
            "shielding": shield,
            "shift": round(shift, 2),
        })

    # Separate and sort by nucleus
    h_shifts = sorted([s for s in shifts if s["atom"] == "H"],
                     key=lambda x: x["shift"], reverse=True)
    c_shifts = sorted([s for s in shifts if s["atom"] == "C"],
                     key=lambda x: x["shift"], reverse=True)

    return {"1H": h_shifts, "13C": c_shifts}
```

### Enhanced API Response with Factor Metadata
```python
# Source: Phase 11 CONTEXT.md - include factor metadata
def build_nmr_results_response(
    nmr_results: NMRResults,
    functional: str,
    basis_set: str,
    solvent: str,
) -> NMRResultsResponse:
    """Build API response with factor metadata."""
    # Get factor metadata for both nuclei
    h1_factor = get_scaling_factor(functional, basis_set, "1H", solvent)
    c13_factor = get_scaling_factor(functional, basis_set, "13C", solvent)

    return NMRResultsResponse(
        h1_shifts=[
            AtomShiftResponse(index=s.index, atom=s.atom, shift=s.shift)
            for s in nmr_results.h1_shifts
        ],
        c13_shifts=[
            AtomShiftResponse(index=s.index, atom=s.atom, shift=s.shift)
            for s in nmr_results.c13_shifts
        ],
        functional=nmr_results.functional,
        basis_set=nmr_results.basis_set,
        solvent=nmr_results.solvent,
        # Factor metadata (new fields)
        scaling_factor_source="DELTA50",
        h1_expected_mae=round(h1_factor["mae"], 2),
        c13_expected_mae=round(c13_factor["mae"], 2),
        h1_r_squared=round(h1_factor["r_squared"], 4),
        c13_r_squared=round(c13_factor["r_squared"], 4),
    )
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| TMS reference (CHESHIRE) | Regression factors (DELTA50) | Phase 11 (2026) | More accurate shifts, solvent-specific factors |
| ISiCLE wrapper | Direct NWChem | Phase 7 (2026) | Removed unnecessary abstraction layer |
| pkg_resources | importlib.resources | Python 3.9+ (2020) | Faster, stdlib, works in zip imports |
| Hardcoded factors | JSON data file | Phase 11 (2026) | Extensible to new functionals/solvents |
| Generic error messages | Explicit factor availability | Phase 11 (2026) | Users know exactly what's supported |

**Deprecated/outdated:**
- **pkg_resources:** Slower than importlib.resources, adds packaging dependency overhead
- **ujson:** Maintainers put in maintenance-only mode, orjson now preferred
- **TMS reference approach:** Less accurate than regression factors, doesn't account for DFT systematic errors
- **"isicle_version" field:** ISiCLE removed in Phase 7, field persisted for compatibility, now being removed completely

## Open Questions

Things that couldn't be fully resolved:

1. **Should factor loading be thread-safe?**
   - What we know: @functools.cache is thread-safe in CPython via GIL
   - What's unclear: Whether Huey workers might cause race conditions on first load
   - Recommendation: Test with concurrent worker startup; if issues arise, use threading.Lock

2. **How to handle legacy job status.json files with isicle_version?**
   - What we know: Pydantic v2 has model_config(extra="ignore") for extra fields
   - What's unclear: Whether ALL old jobs should remain loadable or just warn/migrate
   - Recommendation: Use extra="ignore" to allow loading old files without breaking; document that isicle_version is ignored

3. **Should scaling_factors.json stay in data/benchmark/delta50/ or move to src/qm_nmr_calc/data/?**
   - What we know: Best practice is data files inside package for importlib.resources
   - What's unclear: Whether to keep benchmark data separate or consolidate
   - Recommendation: Copy to src/qm_nmr_calc/data/ for production use; keep original in data/benchmark/delta50/ for reproducibility

## Sources

### Primary (HIGH confidence)
- [Python importlib.resources documentation](https://docs.python.org/3/library/importlib.resources.html) - Official stdlib docs (updated 2026-01-22)
- [Scientific Python data files guide](https://learn.scientific-python.org/development/patterns/data-files/) - Authoritative patterns for package data
- [Python Packaging Best Practices 2026](https://dasroot.net/posts/2026/01/python-packaging-best-practices-setuptools-poetry-hatch/) - Current build system patterns
- Project codebase analysis - shifts.py, models.py, tasks.py examined directly

### Secondary (MEDIUM confidence)
- [orjson performance comparison](https://en.chr.fan/2026/01/07/python-json/) - 280ms vs 420ms benchmark (Jan 2026)
- [Managing API Changes 2026](https://www.theneo.io/blog/managing-api-changes-strategies) - 70% reduction in update-related incidents
- [Python dependency management best practices](https://www.geeksforgeeks.org/python/best-practices-for-managing-python-dependencies/) - Auditing and removal strategies
- [Setuptools 78.0.1 breaking changes](https://pydevtools.com/blog/setuptools-78-0-1-breaking-package-installation/) - Recent real-world example (Jan 2026)

### Tertiary (LOW confidence)
- [Graceful degradation patterns](https://medium.com/@RampantLions/robust-error-handling-in-python-tracebacks-graceful-degradation-and-suppression-11f7a140720b) - General patterns (Jul 2025)
- [Python versioning semantic](https://fonzi.ai/blog/python-versioning-semantic) - Background on SemVer

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - Using only stdlib and existing dependencies, verified in docs
- Architecture: HIGH - importlib.resources is official pattern, factor lookup is straightforward
- Pitfalls: MEDIUM - Some derived from general best practices, others from codebase analysis
- Code examples: HIGH - Based on official docs and project CONTEXT.md decisions

**Research date:** 2026-01-23
**Valid until:** 2026-04-23 (90 days - Python ecosystem stable, importlib.resources is mature stdlib)
