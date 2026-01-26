# Phase 10: Scaling Factors - Research

**Researched:** 2026-01-23
**Domain:** NMR chemical shift linear regression / statistical analysis
**Confidence:** HIGH

## Summary

This phase derives NWChem-specific scaling factors from DELTA50 benchmark data using linear regression. The standard methodology for NMR scaling factor derivation is well-established: fit calculated shielding values against experimental chemical shifts using ordinary least squares (OLS), optionally remove outliers via residual analysis, and validate on an external dataset.

The existing Python stack (scipy 1.17.0, statsmodels 0.14.6, matplotlib 3.10.8, pandas) provides all necessary tools. statsmodels OLS gives confidence intervals for slope/intercept directly, scipy.stats.bootstrap handles MAE/RMSD uncertainty, and matplotlib produces publication-quality plots. No additional dependencies are required.

The DELTA50 methodology (Grimblat et al. 2023) provides the template: derive factors from a training set, validate on external compounds. The validation set will be manually curated from nmrshiftdb2 (10-20 small molecules with high-quality experimental data).

**Primary recommendation:** Use statsmodels OLS for regression with 3-sigma residual-based outlier removal, scipy.stats.bootstrap for MAE/RMSD uncertainty, and store factors as Python constants with JSON export for external tools.

## Standard Stack

The established libraries/tools for this domain:

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| statsmodels | 0.14.6 | OLS regression with confidence intervals | Industry standard for statistical modeling in Python, provides conf_int() directly |
| scipy | 1.17.0 | Bootstrap resampling for uncertainty | scipy.stats.bootstrap is the canonical implementation |
| pandas | 2.3.3 | Data manipulation and aggregation | Already in project, handles multi-index grouping well |
| matplotlib | 3.10.8 | Publication-quality plots | Already in project, full control over plot styling |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| numpy | (via scipy) | Array operations, statistics | Mean, std, residual calculations |
| sklearn.metrics | (optional) | MAE/RMSE calculation | Alternative to manual calculation |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| statsmodels OLS | scipy.stats.linregress | linregress lacks conf_int(); would need manual calculation |
| statsmodels OLS | sklearn LinearRegression | sklearn lacks p-values and confidence intervals |
| Manual bootstrap | scipy.stats.bootstrap | scipy version is more robust, supports BCa method |
| seaborn regplot | matplotlib scatter + line | seaborn adds dependency; matplotlib sufficient for this use case |

**Installation:**
```bash
# No new dependencies needed - all already in pyproject.toml
# statsmodels, scipy, pandas, matplotlib are present
```

## Architecture Patterns

### Recommended Project Structure
```
src/qm_nmr_calc/benchmark/
    __init__.py
    data_loader.py          # Existing: loads DELTA50 molecules
    models.py               # Existing: BenchmarkResult, MoleculeData
    runner.py               # Existing: runs benchmark calculations
    analysis.py             # NEW: scaling factor derivation and analysis
    __main__.py             # Existing: CLI, extend with analyze command
```

### Pattern 1: Factor Organization by Lookup Key
**What:** Store scaling factors with composite key (functional, basis_set, nucleus, solvent)
**When to use:** All factor lookups in production code
**Example:**
```python
# Source: CONTEXT.md decision on factor organization
SCALING_FACTORS: dict[str, ScalingFactor] = {
    "B3LYP/6-311+G(2d,p)/1H/CHCl3": ScalingFactor(
        slope=-1.0234,
        intercept=31.8765,
        r_squared=0.9987,
        mae=0.12,
        rmsd=0.15,
        n_points=114,
        ci_slope=(−1.0289, −1.0179),
        ci_intercept=(31.8234, 31.9296),
    ),
    # ... additional factor sets
}

def get_factor_key(functional: str, basis_set: str, nucleus: str, solvent: str) -> str:
    """Generate lookup key for scaling factors."""
    return f"{functional}/{basis_set}/{nucleus}/{solvent}"
```

### Pattern 2: Data Pairing via Atom Assignments
**What:** Match calculated shielding values to experimental shifts using atom index assignments
**When to use:** Building regression datasets from shifts.json and experimental_shifts.json
**Example:**
```python
# Source: Existing experimental_shifts.json structure
def build_regression_data(
    shifts_json: dict,
    exp_data: MoleculeData,
    nucleus: str,  # "1H" or "13C"
) -> list[tuple[float, float]]:
    """Build (shielding, experimental_shift) pairs.

    Returns list of tuples for regression fitting.
    """
    pairs = []
    assignments = exp_data.h1_assignments if nucleus == "1H" else exp_data.c13_assignments

    if assignments is None:
        return pairs  # Skip molecules without assignments

    for atom_idx_str, exp_shift in assignments.items():
        atom_idx = int(atom_idx_str)
        # Find shielding for this atom index
        for i, idx in enumerate(shifts_json["shielding_data"]["index"]):
            if idx == atom_idx:
                shielding = shifts_json["shielding_data"]["shielding"][i]
                pairs.append((shielding, exp_shift))
                break

    return pairs
```

### Pattern 3: OLS with Residual-Based Outlier Removal
**What:** Fit OLS, identify outliers by residual > 3 std dev, refit without outliers
**When to use:** Standard scaling factor derivation
**Example:**
```python
# Source: statsmodels official documentation + CONTEXT.md decision
import statsmodels.api as sm
import numpy as np

def fit_with_outlier_removal(
    shielding: np.ndarray,
    shifts: np.ndarray,
    threshold: float = 3.0,
) -> tuple[sm.regression.linear_model.RegressionResultsWrapper, np.ndarray]:
    """Fit OLS regression with residual-based outlier removal.

    Args:
        shielding: Array of calculated shielding values (X)
        shifts: Array of experimental shift values (y)
        threshold: Number of std devs for outlier detection

    Returns:
        (fitted_model, outlier_mask) where outlier_mask is True for outliers
    """
    X = sm.add_constant(shielding)

    # Initial fit
    model = sm.OLS(shifts, X)
    results = model.fit()

    # Identify outliers
    residuals = results.resid
    std_resid = np.std(residuals)
    outlier_mask = np.abs(residuals) > threshold * std_resid

    if np.any(outlier_mask):
        # Refit without outliers
        X_clean = X[~outlier_mask]
        y_clean = shifts[~outlier_mask]
        model_clean = sm.OLS(y_clean, X_clean)
        results = model_clean.fit()

    return results, outlier_mask
```

### Pattern 4: Bootstrap Uncertainty for MAE/RMSD
**What:** Use scipy.stats.bootstrap to estimate confidence intervals for error metrics
**When to use:** Reporting uncertainty on MAE/RMSD statistics
**Example:**
```python
# Source: scipy.stats.bootstrap official documentation
from scipy.stats import bootstrap
import numpy as np

def bootstrap_mae_ci(
    residuals: np.ndarray,
    confidence_level: float = 0.95,
    n_resamples: int = 9999,
) -> tuple[float, float, float]:
    """Calculate MAE with bootstrap confidence interval.

    Returns:
        (mae, ci_low, ci_high)
    """
    def mae_stat(x, axis):
        return np.mean(np.abs(x), axis=axis)

    data = (residuals,)
    result = bootstrap(
        data,
        mae_stat,
        confidence_level=confidence_level,
        n_resamples=n_resamples,
        method='percentile',
    )

    mae = np.mean(np.abs(residuals))
    return mae, result.confidence_interval.low, result.confidence_interval.high
```

### Anti-Patterns to Avoid
- **Using sklearn for regression when CIs needed:** sklearn LinearRegression lacks confidence intervals; use statsmodels
- **TMS reference instead of empirical scaling:** Simple TMS subtraction has higher errors than linear scaling
- **Train/test split on DELTA50:** Dataset too small (50 compounds); use external validation instead
- **Ignoring solvent effects:** Factors differ significantly between CHCl3 and DMSO; never mix solvents in fitting
- **Combining nuclei in single fit:** 1H and 13C have different shielding ranges; always fit separately

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Confidence intervals for slope/intercept | Manual t-distribution calculation | `results.conf_int()` | statsmodels handles degrees of freedom correctly |
| Bootstrap resampling | Manual loop with random sampling | `scipy.stats.bootstrap` | Supports BCa method, handles edge cases |
| Residual diagnostics | Manual residual calculation | `results.get_influence()` | Provides DFBETAS, Cook's distance, etc. |
| RMSE calculation | `np.sqrt(np.mean(residuals**2))` | `statsmodels.tools.eval_measures.rmse` | Handles edge cases, consistent API |
| Regression plot with CI band | Manual polygon/fill_between | matplotlib scatter + line (simple enough) | No need for seaborn dependency |

**Key insight:** Statistical calculations have subtle edge cases (degrees of freedom, bias correction). Use statsmodels/scipy implementations rather than reimplementing formulas.

## Common Pitfalls

### Pitfall 1: Mismatched Atom Indices
**What goes wrong:** Calculated shielding paired with wrong experimental shift
**Why it happens:** Atom numbering in XYZ differs from assignment keys in JSON
**How to avoid:** Always use explicit atom index mapping from h1_assignments/c13_assignments
**Warning signs:** Implausible regression (R² < 0.99 for NMR data)

### Pitfall 2: Including Non-Assigned Atoms
**What goes wrong:** Atoms without experimental assignments included in regression
**Why it happens:** Iterating over all atoms instead of only assigned ones
**How to avoid:** Build pairs only from atoms with explicit experimental assignments
**Warning signs:** More data points than expected, poor fit quality

### Pitfall 3: Solvent Contamination
**What goes wrong:** Mixing CHCl3 and DMSO data in single regression
**Why it happens:** Not filtering by solvent before fitting
**How to avoid:** Explicit solvent parameter in all data loading functions
**Warning signs:** Higher than expected errors, bimodal residual distribution

### Pitfall 4: Wrong Regression Direction
**What goes wrong:** Fitting shielding = f(shift) instead of shift = f(shielding)
**Why it happens:** Confusion about which variable is X vs Y
**How to avoid:** Document clearly: X = shielding (calculated), Y = shift (experimental)
**Warning signs:** Slope has wrong sign (should be negative for this convention)

### Pitfall 5: Forgetting to Add Constant
**What goes wrong:** Regression through origin instead of with intercept
**Why it happens:** Not using sm.add_constant() before fitting
**How to avoid:** Always X = sm.add_constant(X) for OLS
**Warning signs:** Very high intercept error, poor R²

### Pitfall 6: Report Metrics on Training Data Only
**What goes wrong:** Overly optimistic accuracy claims
**Why it happens:** Not separating training (DELTA50) from validation (nmrshiftdb)
**How to avoid:** Always report both training and external validation metrics
**Warning signs:** Validation MAE significantly worse than training MAE

## Code Examples

Verified patterns from official sources:

### Complete OLS Regression with Statistics
```python
# Source: statsmodels official documentation
import statsmodels.api as sm
import numpy as np

# Prepare data
shielding = np.array([...])  # X values (calculated)
shifts = np.array([...])      # Y values (experimental)
X = sm.add_constant(shielding)

# Fit model
model = sm.OLS(shifts, X)
results = model.fit()

# Extract statistics
intercept, slope = results.params  # [const, shielding]
ci = results.conf_int(alpha=0.05)  # 95% CI for parameters
r_squared = results.rsquared
n_points = len(shifts)

# Calculate MAE and RMSD
predictions = results.fittedvalues
residuals = results.resid
mae = np.mean(np.abs(residuals))
rmsd = np.sqrt(np.mean(residuals**2))
```

### Pydantic Model for Scaling Factor
```python
# Source: Project convention from models.py
from pydantic import BaseModel

class ScalingFactor(BaseModel):
    """NMR scaling factor with metadata."""

    slope: float
    intercept: float
    r_squared: float
    mae: float
    rmsd: float
    n_points: int
    ci_slope: tuple[float, float]
    ci_intercept: tuple[float, float]
    mae_ci: tuple[float, float] | None = None
    rmsd_ci: tuple[float, float] | None = None
    outliers_removed: int = 0
```

### Matplotlib Regression Plot with Residuals
```python
# Source: matplotlib gallery patterns
import matplotlib.pyplot as plt
import numpy as np

def plot_regression_analysis(
    shielding: np.ndarray,
    shifts: np.ndarray,
    fitted: np.ndarray,
    residuals: np.ndarray,
    title: str,
    output_path: str,
):
    """Create regression scatter plot with residual subplot."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Scatter plot with regression line
    ax1.scatter(shielding, shifts, alpha=0.6, label='Data')
    sort_idx = np.argsort(shielding)
    ax1.plot(shielding[sort_idx], fitted[sort_idx], 'r-', label='Fit')
    ax1.set_xlabel('Calculated Shielding (ppm)')
    ax1.set_ylabel('Experimental Shift (ppm)')
    ax1.set_title(f'{title} - Regression')
    ax1.legend()

    # Residual plot
    ax2.scatter(fitted, residuals, alpha=0.6)
    ax2.axhline(y=0, color='r', linestyle='--')
    ax2.set_xlabel('Fitted Values (ppm)')
    ax2.set_ylabel('Residuals (ppm)')
    ax2.set_title(f'{title} - Residuals')

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
```

### Data Aggregation Across Compounds
```python
# Source: pandas patterns for multi-level grouping
import pandas as pd
from pathlib import Path
import orjson

def aggregate_all_data(
    results_dir: Path,
    exp_shifts: ExperimentalShifts,
    functional: str,
    solvent: str,
) -> pd.DataFrame:
    """Aggregate shielding-shift pairs across all compounds.

    Returns DataFrame with columns: compound, atom_idx, nucleus, shielding, exp_shift
    """
    records = []

    for mol_id, mol_data in exp_shifts.molecules.items():
        shifts_file = results_dir / mol_id / f"{functional}_{solvent}" / "shifts.json"
        if not shifts_file.exists():
            continue

        data = orjson.loads(shifts_file.read_bytes())

        # Process 1H assignments
        if mol_data.h1_assignments:
            for atom_idx_str, exp_shift in mol_data.h1_assignments.items():
                atom_idx = int(atom_idx_str)
                idx_in_list = data["shielding_data"]["index"].index(atom_idx)
                shielding = data["shielding_data"]["shielding"][idx_in_list]
                records.append({
                    "compound": mol_id,
                    "atom_idx": atom_idx,
                    "nucleus": "1H",
                    "shielding": shielding,
                    "exp_shift": exp_shift,
                })

        # Process 13C assignments
        if mol_data.c13_assignments:
            for atom_idx_str, exp_shift in mol_data.c13_assignments.items():
                atom_idx = int(atom_idx_str)
                idx_in_list = data["shielding_data"]["index"].index(atom_idx)
                shielding = data["shielding_data"]["shielding"][idx_in_list]
                records.append({
                    "compound": mol_id,
                    "atom_idx": atom_idx,
                    "nucleus": "13C",
                    "shielding": shielding,
                    "exp_shift": exp_shift,
                })

    return pd.DataFrame(records)
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Simple TMS reference | Empirical linear scaling | ~2005 onwards | 2-3x improvement in accuracy |
| Gas-phase calculations | Solvent models (PCM/SMD/COSMO) | ~2010 onwards | Better accuracy for polar compounds |
| Single scaling factor | Nucleus/solvent-specific factors | ~2014 (Pierens) | Accounts for systematic differences |
| CHESHIRE test set | DELTA50 curated set | 2023 | Avoids problematic conformers |

**Deprecated/outdated:**
- **Single TMS reference subtraction:** Replaced by empirical linear scaling which achieves ~0.1 ppm (1H) and ~2 ppm (13C) accuracy
- **Gas-phase scaling factors for solution NMR:** Solvent effects are significant; always use solvent-matched factors
- **Universal scaling factors:** Different functionals require different factors

## Open Questions

Things that couldn't be fully resolved:

1. **nmrshiftdb2 compound selection criteria**
   - What we know: Need 10-20 small molecules with high-quality experimental data for external validation
   - What's unclear: Specific selection criteria (molecule size, functional groups to include)
   - Recommendation: Select molecules similar in size/complexity to DELTA50, avoid known problematic cases (conformational flexibility, strong H-bonding)

2. **WP04 factor reliability for 13C**
   - What we know: WP04 is optimized for 1H, DELTA50 paper recommends it for 1H
   - What's unclear: Whether WP04 factors for 13C will be useful
   - Recommendation: Derive factors anyway; report with appropriate caveats if accuracy is lower

3. **DMSO experimental data availability**
   - What we know: DELTA50 experimental data is in CDCl3
   - What's unclear: Whether experimental DMSO data exists for validation
   - Recommendation: For DMSO, may need to use DELTA50 CHCl3 experimental data as proxy or note validation limitation

## Sources

### Primary (HIGH confidence)
- [statsmodels OLS documentation](https://www.statsmodels.org/dev/examples/notebooks/generated/ols.html) - OLS API, conf_int(), rsquared
- [scipy.stats.bootstrap documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.bootstrap.html) - Bootstrap API, BCa method
- DELTA50 paper PMC10051451 - Linear scaling methodology, validation approach

### Secondary (MEDIUM confidence)
- [CHESHIRE NMR scaling factors](http://cheshirenmr.info/ScalingFactors.htm) - Historical reference for scaling factor tables (site has expired certificate)
- [Pierens 2014](https://pubmed.ncbi.nlm.nih.gov/24854878/) - Solvent-specific scaling factor methodology
- [nmrshiftdb2](https://nmrshiftdb.nmr.uni-koeln.de/) - External validation data source

### Tertiary (LOW confidence)
- WebSearch results on robust regression (RANSAC, Huber) - Not needed if simple outlier removal works

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - All libraries already in project, well-documented
- Architecture: HIGH - Follows established NMR scaling factor methodology (DELTA50, CHESHIRE)
- Pitfalls: HIGH - Well-known issues in computational chemistry literature
- Validation approach: MEDIUM - nmrshiftdb selection criteria need manual curation

**Research date:** 2026-01-23
**Valid until:** 60 days (stable domain, libraries mature)

---

## Appendix: Key Formula Reference

**Linear scaling equation:**
```
shift = slope * shielding + intercept
```

Where:
- `shift` = predicted chemical shift (ppm, experimental convention)
- `shielding` = calculated isotropic shielding tensor (ppm)
- `slope` = negative (typically -0.95 to -1.05 for 13C, -1.02 to -1.06 for 1H)
- `intercept` = TMS-equivalent reference (typically 180-200 for 13C, 31-32 for 1H)

**Note on sign convention:** The DELTA50 paper uses `delta = intercept - slope * sigma` which is equivalent to our form with `slope` having opposite sign. Our convention (`shift = slope * shielding + intercept` with negative slope) matches the existing CHESHIRE factors in shifts.py.
