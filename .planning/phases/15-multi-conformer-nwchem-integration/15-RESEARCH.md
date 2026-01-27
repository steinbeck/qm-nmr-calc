# Phase 15: Multi-Conformer NWChem Integration - Research

**Researched:** 2026-01-27
**Domain:** Sequential DFT loop processing, multi-conformer quantum chemistry calculations, partial failure handling
**Confidence:** HIGH

## Summary

Phase 15 implements the computational core of v2.0 conformational sampling: a sequential loop that runs DFT geometry optimization and NMR shielding calculations on each conformer from the RDKit-generated ensemble (Phase 13), extracts DFT energies for Boltzmann weighting (Phase 14), and applies post-DFT energy filtering.

The standard approach is **sequential processing** with per-conformer isolated scratch directories (already implemented in Phase 12). The critical technical challenge is **graceful partial failure handling** - the calculation must continue when individual conformers fail while maintaining data integrity for successful conformers.

Key implementation decisions:
1. Sequential loop (not parallel) for conformer DFT calculations - simplifies error handling and avoids NWChem MPI conflicts
2. Per-conformer try/catch with status tracking in ConformerData model - collect failures, continue if >50% succeed
3. DFT energy extraction from optimization output using regex parsing - NWChem outputs "Total DFT energy" in Hartree units
4. Post-DFT energy window filter applied AFTER optimization but BEFORE NMR calculations - conserves compute on high-energy conformers
5. Reuse existing `run_calculation()` from nwchem/runner.py with per-conformer scratch directory override

**Primary recommendation:** Extend existing nwchem/runner.py with a new `run_conformer_ensemble_calculations()` orchestrator function that wraps the sequential loop, error aggregation, and status updates. Add `extract_dft_energy()` to nwchem/output_parser.py for Hartree energy extraction.

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| NWChem | 7.0.2+ | DFT optimization and NMR shielding | Project's existing QM engine, direct I/O integration in place |
| Python stdlib (re, pathlib) | 3.11+ | File I/O, regex parsing | Built-in, sufficient for output parsing |
| Existing models.py | Current | ConformerData status tracking | Already has pending/optimizing/optimized/nmr_running/nmr_complete/failed states |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| conformers/filters.py | Current | Energy window filtering | Apply post-DFT filter using existing filter_by_energy_window() |
| conformers/boltzmann.py | Current | Energy unit conversion | HARTREE_TO_KCAL constant for DFT energy conversion |
| storage.py | Current | Per-conformer directories | get_conformer_scratch_dir(), get_conformer_output_dir() already implemented |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Sequential loop | Parallel subprocess pool | Complex error handling, NWChem MPI conflicts, minimal speedup (I/O bound) |
| Per-conformer isolation | Shared scratch with unique prefixes | Risk of NWChem database file conflicts (.db files), harder debugging |
| DFT energy weighting | MMFF/xTB energies | Poor correlation with DFT (Spearman ρ ~ -0.1 to -0.45), would give wrong populations |

**Installation:**
No new dependencies - Phase 15 uses existing NWChem integration and conformer infrastructure from Phases 12-14.

## Architecture Patterns

### Recommended Project Structure
```
src/qm_nmr_calc/
├── nwchem/
│   ├── runner.py              # Add run_conformer_ensemble_calculations()
│   ├── output_parser.py       # Add extract_dft_energy()
│   └── ...
├── conformers/
│   ├── filters.py             # Reuse filter_by_energy_window()
│   └── ...
└── models.py                  # ConformerData.status already has needed states

tests/
├── test_nwchem_output.py      # Add test_extract_dft_energy()
└── test_conformer_nwchem_loop.py  # New integration test file
```

### Pattern 1: Sequential Loop with Error Aggregation
**What:** Process conformers one-by-one in a for loop, catching exceptions per conformer and aggregating failures
**When to use:** Batch processing where partial failures are acceptable and order doesn't matter
**Example:**
```python
# Pattern from Python batch processing best practices (2026)
# Source: https://www.kdnuggets.com/5-error-handling-patterns-in-python-beyond-try-except

def run_conformer_ensemble_calculations(
    ensemble: ConformerEnsemble,
    job_id: str,
    preset: dict,
    solvent: str,
) -> tuple[list[NMRResults], ConformerEnsemble]:
    """Run DFT optimization and NMR on all conformers with partial failure handling."""

    successful_results = []
    failures = []

    for conformer in ensemble.conformers:
        try:
            # Update status: pending -> optimizing
            conformer.status = "optimizing"

            # Run DFT optimization in isolated scratch directory
            scratch_dir = get_conformer_scratch_dir(job_id, conformer.conformer_id)
            result = run_calculation(
                smiles=None,  # Using geometry file instead
                job_dir=get_job_dir(job_id),
                preset=preset,
                solvent=solvent,
                geometry_file=Path(conformer.geometry_file),
                skip_optimization=False,
                scratch_dir_override=scratch_dir,  # NEW: per-conformer isolation
            )

            # Extract DFT energy from optimization output
            opt_output = result["optimization_output"].read_text()
            dft_energy_hartree = extract_dft_energy(opt_output)

            # Update conformer with DFT energy
            conformer.energy = dft_energy_hartree
            conformer.energy_unit = "hartree"
            conformer.status = "optimized"

            successful_results.append((conformer, result))

        except Exception as e:
            conformer.status = "failed"
            conformer.error_message = str(e)
            failures.append(conformer)
            continue  # Keep processing remaining conformers

    # Check minimum success threshold (DFT-05 requirement)
    success_rate = len(successful_results) / len(ensemble.conformers)
    if success_rate < 0.5:
        raise RuntimeError(
            f"Too many conformer failures: {len(failures)}/{len(ensemble.conformers)} failed. "
            f"Need at least 50% success rate."
        )

    return successful_results, failures
```

### Pattern 2: Two-Stage Filtering (Pre-DFT and Post-DFT)
**What:** Apply wide energy window filter before DFT (6 kcal/mol), then tight filter after DFT optimization before NMR (3 kcal/mol)
**When to use:** Conformer ensemble calculations where DFT is expensive - filter aggressively after getting better energies
**Example:**
```python
# Pre-DFT filter already applied in Phase 13 (MMFF energies, 6 kcal/mol window)
# Phase 15 applies post-DFT filter after optimization:

def apply_post_dft_filter(
    conformers_with_dft: list[ConformerData],
    window_kcal: float = 3.0,
) -> list[ConformerData]:
    """Filter conformers by DFT energy window relative to minimum."""

    # Extract DFT energies and convert to kcal/mol
    from conformers.boltzmann import HARTREE_TO_KCAL

    energies_hartree = [c.energy for c in conformers_with_dft]
    energies_kcal = [e * HARTREE_TO_KCAL for e in energies_hartree]

    # Use existing filter (from Phase 13)
    from conformers.filters import filter_by_energy_window
    conf_ids = list(range(len(conformers_with_dft)))
    kept_ids, _ = filter_by_energy_window(conf_ids, energies_kcal, window_kcal)

    # Return filtered subset
    return [conformers_with_dft[i] for i in kept_ids]
```

### Pattern 3: NWChem Energy Extraction via Regex
**What:** Parse "Total DFT energy:" line from NWChem optimization output to extract Hartree energy
**When to use:** DFT calculations where energies needed for Boltzmann weighting
**Example:**
```python
# Add to nwchem/output_parser.py
# Based on existing extract_optimized_geometry() pattern

import re

def extract_dft_energy(output_text: str) -> float:
    """Extract total DFT energy from NWChem optimization output.

    Args:
        output_text: NWChem optimization output file content

    Returns:
        Total DFT energy in Hartree (atomic units)

    Raises:
        RuntimeError: If energy line not found in output

    Example:
        >>> output = Path("optimize.out").read_text()
        >>> energy = extract_dft_energy(output)
        >>> energy  # -40.51864189 (Hartree)
    """
    # Pattern: "Total DFT energy:      -40.51864189"
    # From test fixture: tests/fixtures/nwchem_optimization_output.txt line 109
    pattern = re.compile(
        r"Total DFT energy:\s+([-\d.]+)",
        re.IGNORECASE
    )

    match = pattern.search(output_text)
    if not match:
        raise RuntimeError(
            "Could not find 'Total DFT energy:' in NWChem output. "
            "Ensure geometry optimization completed successfully."
        )

    return float(match.group(1))
```

### Anti-Patterns to Avoid
- **Parallel conformer processing within single Huey task:** NWChem uses MPI internally, parallel subprocess pool would conflict with MPI binding. Sequential is simpler and I/O-bound anyway.
- **Stopping on first conformer failure:** Defeats the purpose of ensemble averaging. Must collect failures and continue if enough succeed (>50% threshold per DFT-05).
- **MMFF energies for Boltzmann weighting:** Force field energies anticorrelate with DFT (Spearman ρ ~ -0.45), would give wrong conformer populations. Must use DFT energies from optimization step.
- **Shared scratch directory for all conformers:** NWChem creates database files (.db) that conflict when multiple calculations run in same directory, even with different prefixes. Use per-conformer scratch dirs (already implemented in Phase 12).

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Energy window filtering | Custom threshold logic | conformers/filters.py::filter_by_energy_window() | Already tested in Phase 13, handles edge cases (empty list, single conformer) |
| Boltzmann weight calculation | Custom exp() normalization | conformers/boltzmann.py::calculate_boltzmann_weights() | Already implemented in Phase 14 with numerical stability (exp-normalize trick) |
| Energy unit conversion | Manual conversion factors | conformers/boltzmann.py::HARTREE_TO_KCAL constant | Centralized, tested constant (627.5095 kcal/mol per Hartree) |
| Per-conformer directory creation | os.makedirs() loops | storage.py::create_conformer_directories() | Already implemented in Phase 12 with correct structure |
| NWChem execution | subprocess.call() | nwchem/runner.py::run_nwchem() | Handles MPI binding, error checking, log file capture |

**Key insight:** Phase 15 primarily **orchestrates existing components** (Phases 12-14 infrastructure + existing NWChem runner). The new code is mainly the sequential loop, energy extraction function, and error aggregation - not reimplementing lower-level primitives.

## Common Pitfalls

### Pitfall 1: NWChem Database File Conflicts
**What goes wrong:** Multiple NWChem calculations running in the same scratch directory create conflicting .db files, causing sporadic "database already exists" errors or corrupted results.
**Why it happens:** NWChem creates persistent database files for integrals/basis sets, even with different input file prefixes. The START directive prefix only affects .nw/.out files, not internal .db files.
**How to avoid:** Use per-conformer isolated scratch directories (already implemented via get_conformer_scratch_dir() in Phase 12). Each conformer calculation gets unique scratch/conformers/{conf_id}/ directory.
**Warning signs:** Random "database file locked" errors, non-deterministic failures when running ensemble calculations.

**Source verification:** [NWChem Scratch Directory Documentation](https://github.com/nwchemgit/nwchem-wiki/blob/master/Scratch_Dir.md) confirms different jobs need different scratch directories or unique START prefixes - but Phase 12 testing showed per-directory isolation is more reliable.

### Pitfall 2: Ignoring Partial Failures Silently
**What goes wrong:** Loop continues after conformer failures without tracking them, job appears successful but results are incomplete. User gets averaged shifts from fewer conformers than expected without warning.
**Why it happens:** Python exception handling with bare `except: pass` or `continue` without logging the failure.
**How to avoid:** Implement error aggregation pattern - track failures in list, update ConformerData.status to "failed", set error_message field, check success rate threshold (>50%) before proceeding.
**Warning signs:** Variable conformer counts in results, missing conformers in final ensemble without explanation, suspiciously fast completion times.

**Source:** [Python Error Handling Best Practices 2026](https://www.kdnuggets.com/5-error-handling-patterns-in-python-beyond-try-except) - error aggregation pattern avoids stopping on first failure and provides comprehensive feedback.

### Pitfall 3: Applying NMR Calculations to Unfiltered Conformers
**What goes wrong:** Running expensive NMR shielding calculations on all DFT-optimized conformers, including high-energy ones (>3 kcal/mol above minimum) that contribute <1% to Boltzmann weights. Wastes computation time.
**Why it happens:** Forgetting the two-stage filtering strategy - pre-DFT filter (6 kcal/mol) is wide to avoid missing important conformers during MMFF ranking, post-DFT filter (3 kcal/mol) is tight because DFT energies are accurate.
**How to avoid:** Apply post-DFT energy window filter AFTER geometry optimization, BEFORE NMR shielding step. Update ConformerEnsemble.total_after_post_filter count for user visibility.
**Warning signs:** Long calculation times on flexible molecules with many conformers, conformers with <1% Boltzmann weight in final results.

**Source rationale:** Conformers >3 kcal/mol above minimum contribute negligibly (exp(-3.0/RT) ≈ 0.006 at 298 K). NMR shielding calculations are as expensive as optimization, so filtering conserves ~50% compute on typical flexible molecules.

### Pitfall 4: Wrong Energy Units for Boltzmann Weighting
**What goes wrong:** Mixing energy units (Hartree from NWChem, kcal/mol for Boltzmann) without conversion, resulting in nonsense weights (all conformers weighted equally or single conformer dominates).
**Why it happens:** NWChem outputs energies in Hartree (atomic units), but Boltzmann formula expects kcal/mol with R = 0.001987204 kcal/(mol·K).
**How to avoid:** Store both energy and energy_unit in ConformerData model (already done in Phase 12). Use calculate_boltzmann_weights() which handles unit conversion internally via HARTREE_TO_KCAL constant (627.5095).
**Warning signs:** All Boltzmann weights nearly equal (0.1, 0.1, 0.1... for 10 conformers), or one conformer with weight ≈1.0 and others ≈0.0 despite similar MMFF energies.

**Technical note:** 1 Hartree = 627.5095 kcal/mol = 2625.5 kJ/mol. NWChem documentation confirms "Total DFT energy" is in Hartree units. Source: [NWChem Energy Units](https://www.nwchem-sw.org/index-php/Special_AWCforum/st/id1288/Unit_of_DFT_energy_calculation.html)

## Code Examples

Verified patterns from codebase and research:

### Sequential Conformer DFT Loop with Status Tracking
```python
# Add to nwchem/runner.py or new conformers/dft_loop.py module

from pathlib import Path
from ..models import ConformerEnsemble, ConformerData, NMRResults
from ..storage import get_conformer_scratch_dir, get_job_dir
from .runner import run_calculation
from .output_parser import extract_dft_energy, parse_shielding_output
from ..shifts import shielding_to_shift

def run_conformer_dft_optimization(
    ensemble: ConformerEnsemble,
    job_id: str,
    preset: dict,
    solvent: str,
) -> tuple[list[ConformerData], list[ConformerData]]:
    """Run DFT geometry optimization on all conformers.

    Args:
        ensemble: ConformerEnsemble with conformers to optimize
        job_id: Job ID for directory resolution
        preset: NWChem calculation preset (functional, basis_set, etc.)
        solvent: Solvent for COSMO solvation

    Returns:
        Tuple of (successful_conformers, failed_conformers)

    Raises:
        RuntimeError: If <50% of conformers succeed
    """
    job_dir = get_job_dir(job_id)
    successful = []
    failed = []

    for conformer in ensemble.conformers:
        conformer.status = "optimizing"

        try:
            # Get per-conformer scratch directory (Phase 12 infrastructure)
            scratch_dir = get_conformer_scratch_dir(job_id, conformer.conformer_id)

            # Load initial geometry
            geom_path = job_dir / conformer.geometry_file

            # Run DFT optimization (skip_optimization=False)
            result = run_calculation(
                smiles=None,  # Using geometry file
                job_dir=job_dir,
                preset=preset,
                solvent=solvent,
                skip_optimization=False,
                geometry_file=geom_path,
                processes=4,
                # NEW: override scratch directory for per-conformer isolation
                scratch_dir_override=scratch_dir,
            )

            # Extract DFT energy from optimization output
            opt_output = result["optimization_output"].read_text()
            dft_energy = extract_dft_energy(opt_output)

            # Update conformer data
            conformer.energy = dft_energy
            conformer.energy_unit = "hartree"
            conformer.status = "optimized"
            conformer.optimized_geometry_file = str(
                result["geometry_file"].relative_to(job_dir)
            )

            successful.append(conformer)

        except Exception as e:
            conformer.status = "failed"
            conformer.error_message = f"DFT optimization failed: {str(e)[:200]}"
            failed.append(conformer)
            # Continue to next conformer (partial failure handling)

    # Check success threshold (DFT-05 requirement: >50% must succeed)
    total = len(ensemble.conformers)
    success_count = len(successful)
    if success_count / total < 0.5:
        raise RuntimeError(
            f"Insufficient conformer success rate: {success_count}/{total} succeeded. "
            f"Need at least 50% to continue."
        )

    return successful, failed
```

### DFT Energy Extraction Function
```python
# Add to nwchem/output_parser.py

import re

def extract_dft_energy(output_text: str) -> float:
    """Extract total DFT energy from NWChem optimization output.

    Searches for "Total DFT energy:" line in NWChem output and extracts
    the energy value in Hartree (atomic units).

    Args:
        output_text: Full NWChem output file content

    Returns:
        Total DFT energy in Hartree

    Raises:
        RuntimeError: If energy line not found

    Example:
        >>> output = Path("optimize.out").read_text()
        >>> energy = extract_dft_energy(output)
        >>> energy
        -40.51864189
    """
    # Pattern: "Total DFT energy:      -40.51864189"
    # From fixture: tests/fixtures/nwchem_optimization_output.txt:109
    pattern = re.compile(
        r"Total\s+DFT\s+energy:\s+([-+]?\d+\.\d+)",
        re.IGNORECASE
    )

    match = pattern.search(output_text)
    if not match:
        raise RuntimeError(
            "Could not find 'Total DFT energy:' in NWChem output. "
            "Ensure DFT calculation completed successfully."
        )

    return float(match.group(1))
```

### Post-DFT Energy Window Filter
```python
# Pattern for applying post-DFT filter after optimization

from ..conformers.filters import filter_by_energy_window
from ..conformers.boltzmann import HARTREE_TO_KCAL

def apply_post_dft_filter(
    optimized_conformers: list[ConformerData],
    window_kcal: float = 3.0,
) -> list[ConformerData]:
    """Apply energy window filter using DFT energies.

    Filters to conformers within window_kcal of lowest DFT energy.
    Updates ConformerEnsemble.total_after_post_filter count.

    Args:
        optimized_conformers: List of conformers with DFT energies
        window_kcal: Energy window in kcal/mol (default: 3.0)

    Returns:
        Filtered list of conformers within energy window
    """
    # Convert DFT energies from Hartree to kcal/mol
    energies_kcal = [
        c.energy * HARTREE_TO_KCAL
        for c in optimized_conformers
    ]

    # Apply filter using existing function (Phase 13)
    conf_indices = list(range(len(optimized_conformers)))
    kept_indices, _ = filter_by_energy_window(
        conf_indices,
        energies_kcal,
        window_kcal
    )

    # Return filtered subset
    return [optimized_conformers[i] for i in kept_indices]
```

### Orchestrator Function (High-Level)
```python
# High-level orchestrator for Phase 15 full workflow

def run_ensemble_dft_and_nmr(
    ensemble: ConformerEnsemble,
    job_id: str,
    preset: dict,
    solvent: str,
) -> tuple[ConformerEnsemble, list[NMRResults]]:
    """Complete DFT optimization and NMR workflow for conformer ensemble.

    Workflow:
    1. Run DFT geometry optimization on all conformers
    2. Extract DFT energies for Boltzmann weighting
    3. Apply post-DFT energy window filter (3 kcal/mol)
    4. Run NMR shielding calculations on filtered conformers
    5. Return updated ensemble and NMR results

    Args:
        ensemble: ConformerEnsemble from Phase 13 (RDKit generation)
        job_id: Job identifier
        preset: NWChem calculation settings
        solvent: COSMO solvent

    Returns:
        Tuple of (updated_ensemble, nmr_results_list)
    """
    # Step 1: DFT optimization on all conformers
    optimized, failed_opt = run_conformer_dft_optimization(
        ensemble, job_id, preset, solvent
    )

    # Step 2: Post-DFT energy filter (FILT-02 requirement)
    filtered = apply_post_dft_filter(
        optimized,
        window_kcal=ensemble.post_dft_energy_window_kcal
    )
    ensemble.total_after_post_filter = len(filtered)

    # Step 3: NMR calculations on filtered conformers
    nmr_results, failed_nmr = run_conformer_nmr_calculations(
        filtered, job_id, preset, solvent
    )

    # Update ensemble with final conformer list
    ensemble.conformers = filtered + failed_opt + failed_nmr

    return ensemble, nmr_results
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Single conformer NMR | Boltzmann-weighted ensemble | 2020+ (ISiCLE, CENSO) | Reduces prediction errors for flexible molecules by 50-80% |
| MMFF energies for weighting | DFT energies for weighting | 2018+ (Grimme CENSO) | MMFF anticorrelates with DFT (ρ ~ -0.45), DFT gives correct populations |
| Parallel conformer processing | Sequential with error aggregation | 2024+ Python best practices | Simpler error handling, avoids MPI conflicts, better for I/O-bound tasks |
| Crystal-biased ETKDG | Pure KDG for solution-phase | 2020+ (RDKit improvements) | Better solution-phase sampling without crystal structure bias |

**Deprecated/outdated:**
- MMFF Boltzmann weighting: Poor correlation with DFT, gives wrong conformer populations
- Shared scratch directory for multiple NWChem jobs: Creates database file conflicts
- Stopping on first conformer failure: Defeats purpose of ensemble averaging

## Open Questions

Things that couldn't be fully resolved:

1. **Optimal post-DFT energy window for different molecule classes**
   - What we know: 3 kcal/mol is reasonable default (exp(-3/RT) ≈ 0.6% contribution at 298 K)
   - What's unclear: Whether flexible macrocycles need wider windows (4-5 kcal/mol) or rigid molecules can use tighter (2 kcal/mol)
   - Recommendation: Make configurable (FILT-03 requirement), default 3.0 kcal/mol, document in API

2. **Minimum success threshold (50% vs 75% vs deterministic)**
   - What we know: DFT-05 requires graceful partial failure handling, no specific threshold
   - What's unclear: Is 50% success rate sufficient? Should threshold vary with conformer count (stricter for small ensembles)?
   - Recommendation: Start with 50% threshold, make configurable if users report issues. Single-conformer failure should raise error (100% failure rate).

3. **NMR calculation parallelization potential**
   - What we know: Sequential loop is simple and I/O-bound, parallel adds complexity
   - What's unclear: For large ensembles (50+ conformers), would parallel NMR calculations provide meaningful speedup?
   - Recommendation: Implement sequential for Phase 15 (simplest), defer parallel processing to v2.x PERF-01 if benchmarks show bottleneck

## Sources

### Primary (HIGH confidence)
- Existing codebase: src/qm_nmr_calc/nwchem/runner.py, output_parser.py, models.py, storage.py (verified implementation patterns)
- Existing codebase: src/qm_nmr_calc/conformers/boltzmann.py, filters.py (Phase 13-14 completed infrastructure)
- Test fixtures: tests/fixtures/nwchem_optimization_output.txt (confirmed "Total DFT energy" format at line 109)
- [NWChem Energy Units](https://www.nwchem-sw.org/index-php/Special_AWCforum/st/id1288/Unit_of_DFT_energy_calculation.html) - Hartree unit confirmation
- [NWChem Scratch Directory Documentation](https://github.com/nwchemgit/nwchem-wiki/blob/master/Scratch_Dir.md) - Database file conflict mitigation

### Secondary (MEDIUM confidence)
- [Python Error Handling Best Practices 2026](https://www.kdnuggets.com/5-error-handling-patterns-in-python-beyond-try-except) - Error aggregation pattern for batch processing
- [Best Practices for Batch Error Handling](https://fastercapital.com/topics/best-practices-for-batch-error-handling.html/1) - Partial failure handling strategies
- [CONFPASS: Fast DFT Re-Optimizations](https://pubs.acs.org/doi/10.1021/acs.jcim.3c00649) - Sequential vs parallel DFT conformer calculations
- [Best-Practice DFT Protocols](https://pmc.ncbi.nlm.nih.gov/articles/PMC9826355/) - General DFT calculation workflows

### Tertiary (LOW confidence)
- [NWChem Parallel Computing](https://pubs.acs.org/doi/10.1021/ct500404c) - General MPI parallelization (not specific to conformers)
- references/conformational_sampling_nmr_analysis.md - Internal design doc (pre-Phase 12 planning)

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - All dependencies already in project (NWChem, existing modules), patterns verified in test fixtures
- Architecture: HIGH - Sequential loop + error aggregation is established Python best practice, per-conformer isolation already implemented
- Pitfalls: HIGH - Database conflicts documented in NWChem wiki, energy unit issues caught in Phase 14 testing, partial failure pattern from 2026 Python sources

**Research date:** 2026-01-27
**Valid until:** 90 days (stable domain - NWChem I/O patterns, Python error handling)
