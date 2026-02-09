"""Scaling factor derivation from DELTA50 benchmark data.

Derives NMR scaling factors via linear regression of calculated shielding
values against experimental chemical shifts.

Usage:
    from qm_nmr_calc.benchmark.analysis import derive_all_factors

    factors = derive_all_factors()
    # Returns dict mapping "functional/basis_set/nucleus/solvent" to ScalingFactor
"""

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # Headless backend - MUST be before pyplot import
import matplotlib.pyplot as plt
import numpy as np
import orjson
import pandas as pd
import statsmodels.api as sm

from .data_loader import load_experimental_shifts
from .models import ScalingFactor
from .runner import BENCHMARK_PRESETS, get_results_dir

logger = logging.getLogger(__name__)


def get_factor_key(functional: str, basis_set: str, nucleus: str, solvent: str) -> str:
    """Generate lookup key for scaling factors.

    Args:
        functional: DFT functional (e.g., "B3LYP")
        basis_set: NMR basis set (e.g., "6-311+G(2d,p)")
        nucleus: "1H" or "13C"
        solvent: Solvent name (e.g., "CHCl3")

    Returns:
        Composite key string like "B3LYP/6-311+G(2d,p)/1H/CHCl3"
    """
    return f"{functional}/{basis_set}/{nucleus}/{solvent}"


def aggregate_regression_data(functional: str, solvent: str) -> pd.DataFrame:
    """Aggregate shielding-shift pairs across all DELTA50 compounds.

    Loads experimental assignments and pairs them with calculated shielding
    values from benchmark results.

    Args:
        functional: DFT functional (e.g., "B3LYP")
        solvent: Solvent name (e.g., "CHCl3")

    Returns:
        DataFrame with columns: compound, atom_idx, nucleus, shielding, exp_shift
        Empty DataFrame if no data found.
    """
    exp_shifts = load_experimental_shifts()
    results_dir = get_results_dir()
    records: list[dict] = []

    for mol_id, mol_data in exp_shifts.molecules.items():
        # Load shifts.json for this molecule/method combination
        shifts_file = results_dir / mol_id / f"{functional}_{solvent}" / "shifts.json"

        if not shifts_file.exists():
            logger.debug(f"No shifts.json for {mol_id}/{functional}_{solvent}, skipping")
            continue

        try:
            data = orjson.loads(shifts_file.read_bytes())
        except Exception as e:
            logger.warning(f"Failed to read {shifts_file}: {e}")
            continue

        shielding_data = data.get("shielding_data", {})
        indices = shielding_data.get("index", [])
        shieldings = shielding_data.get("shielding", [])

        # Build index -> shielding lookup
        idx_to_shielding = dict(zip(indices, shieldings))

        # Process 1H assignments
        if mol_data.h1_assignments:
            for atom_idx_str, exp_shift in mol_data.h1_assignments.items():
                atom_idx = int(atom_idx_str)
                if atom_idx in idx_to_shielding:
                    records.append({
                        "compound": mol_id,
                        "atom_idx": atom_idx,
                        "nucleus": "1H",
                        "shielding": idx_to_shielding[atom_idx],
                        "exp_shift": exp_shift,
                    })
                else:
                    logger.warning(
                        f"{mol_id}: H atom index {atom_idx} not in shielding data"
                    )

        # Process 13C assignments
        if mol_data.c13_assignments:
            for atom_idx_str, exp_shift in mol_data.c13_assignments.items():
                atom_idx = int(atom_idx_str)
                if atom_idx in idx_to_shielding:
                    records.append({
                        "compound": mol_id,
                        "atom_idx": atom_idx,
                        "nucleus": "13C",
                        "shielding": idx_to_shielding[atom_idx],
                        "exp_shift": exp_shift,
                    })
                else:
                    logger.warning(
                        f"{mol_id}: C atom index {atom_idx} not in shielding data"
                    )

    return pd.DataFrame(records)


def fit_scaling_factors(
    df: pd.DataFrame,
    nucleus: str,
    threshold: float = 3.0,
) -> ScalingFactor:
    """Fit OLS regression with residual-based outlier removal.

    Fits shift = slope * shielding + intercept using OLS regression.
    Outliers (|residual| > threshold * std) are removed and model is refit.

    Args:
        df: DataFrame with columns shielding and exp_shift
        nucleus: "1H" or "13C" to filter data
        threshold: Number of standard deviations for outlier detection

    Returns:
        ScalingFactor with regression statistics

    Raises:
        ValueError: If insufficient data points for regression
    """
    # Filter to specified nucleus
    data = df[df["nucleus"] == nucleus].copy()

    if len(data) < 3:
        raise ValueError(f"Insufficient data for {nucleus} regression: {len(data)} points")

    shielding = data["shielding"].values
    shifts = data["exp_shift"].values

    # Add constant for intercept
    X = sm.add_constant(shielding)

    # Initial fit
    model = sm.OLS(shifts, X)
    results = model.fit()

    # Identify outliers based on residuals
    residuals = results.resid
    std_resid = np.std(residuals)
    outlier_mask = np.abs(residuals) > threshold * std_resid
    n_outliers = int(np.sum(outlier_mask))

    # Refit without outliers if any found
    if n_outliers > 0:
        logger.info(f"Removing {n_outliers} outliers from {nucleus} regression")
        X_clean = X[~outlier_mask]
        y_clean = shifts[~outlier_mask]
        model_clean = sm.OLS(y_clean, X_clean)
        results = model_clean.fit()
        # Recalculate residuals for final metrics
        residuals = results.resid

    # Extract parameters: [intercept, slope]
    intercept, slope = results.params
    ci = results.conf_int(alpha=0.05)  # 95% CI, shape (2, 2)

    # Calculate error metrics
    mae = float(np.mean(np.abs(residuals)))
    rmsd = float(np.sqrt(np.mean(residuals**2)))

    return ScalingFactor(
        slope=float(slope),
        intercept=float(intercept),
        r_squared=float(results.rsquared),
        mae=mae,
        rmsd=rmsd,
        n_points=len(residuals),
        ci_slope=(float(ci[1, 0]), float(ci[1, 1])),  # Row 1 is slope
        ci_intercept=(float(ci[0, 0]), float(ci[0, 1])),  # Row 0 is intercept
        outliers_removed=n_outliers,
    )


def derive_all_factors(
    functionals: list[str] | None = None,
    solvents: list[str] | None = None,
) -> dict[str, ScalingFactor]:
    """Derive scaling factors for all functional/solvent/nucleus combinations.

    Args:
        functionals: List of functionals to process (default: ["B3LYP"])
        solvents: List of solvents to process (default: all 6 supported solvents)

    Returns:
        Dict mapping "functional/basis_set/nucleus/solvent" to ScalingFactor
    """
    if functionals is None:
        functionals = ["B3LYP"]  # WP04 not yet complete
    if solvents is None:
        solvents = ["CHCl3", "DMSO", "Methanol", "Water", "Acetone", "Benzene"]

    nuclei = ["1H", "13C"]
    factors: dict[str, ScalingFactor] = {}

    for functional in functionals:
        basis_set = BENCHMARK_PRESETS[functional]["nmr_basis_set"]

        for solvent in solvents:
            # Aggregate data for this functional/solvent combination
            df = aggregate_regression_data(functional, solvent)

            if df.empty:
                logger.warning(
                    f"No data for {functional}/{solvent}, skipping factor derivation"
                )
                continue

            for nucleus in nuclei:
                try:
                    factor = fit_scaling_factors(df, nucleus)
                    key = get_factor_key(functional, basis_set, nucleus, solvent)
                    factors[key] = factor
                    logger.info(
                        f"Derived {key}: slope={factor.slope:.4f}, "
                        f"intercept={factor.intercept:.2f}, R^2={factor.r_squared:.4f}"
                    )
                except ValueError as e:
                    logger.warning(f"Could not derive {functional}/{nucleus}/{solvent}: {e}")

    return factors


def plot_regression(
    df: pd.DataFrame,
    factor: ScalingFactor,
    nucleus: str,
    title: str,
    output_path: Path,
) -> None:
    """Create regression scatter plot with residual subplot.

    Left subplot: scatter of shielding vs experimental shift with regression line.
    Right subplot: residual plot showing fitted values vs residuals.

    Args:
        df: DataFrame with columns shielding and exp_shift (pre-filtered by nucleus)
        factor: ScalingFactor with regression parameters
        nucleus: "1H" or "13C" for filtering
        title: Plot title
        output_path: Path to save PNG file
    """
    # Filter to specified nucleus
    data = df[df["nucleus"] == nucleus].copy()
    if data.empty:
        logger.warning(f"No data for {nucleus} in plot_regression")
        return

    shielding = data["shielding"].values
    exp_shift = data["exp_shift"].values

    # Calculate fitted values and residuals
    fitted = factor.slope * shielding + factor.intercept
    residuals = exp_shift - fitted

    # Identify outliers for coloring (using stored outlier count as hint)
    # Re-compute outlier mask for visualization
    std_resid = np.std(residuals)
    outlier_mask = np.abs(residuals) > 3.0 * std_resid

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Left subplot: Scatter plot with regression line
    ax1.scatter(
        shielding[~outlier_mask],
        exp_shift[~outlier_mask],
        alpha=0.6,
        label="Data",
        color="C0",
    )
    if np.any(outlier_mask):
        ax1.scatter(
            shielding[outlier_mask],
            exp_shift[outlier_mask],
            alpha=0.8,
            label="Outliers",
            color="red",
            marker="x",
        )

    # Regression line
    sort_idx = np.argsort(shielding)
    ax1.plot(
        shielding[sort_idx], fitted[sort_idx], "r-", linewidth=2, label="Fit"
    )

    # Equation annotation
    eq_text = (
        f"shift = {factor.slope:.4f} * shielding + {factor.intercept:.2f}\n"
        f"R$^2$ = {factor.r_squared:.4f}"
    )
    ax1.annotate(
        eq_text,
        xy=(0.05, 0.95),
        xycoords="axes fraction",
        fontsize=10,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
    )

    ax1.set_xlabel("Calculated Shielding (ppm)")
    ax1.set_ylabel("Experimental Shift (ppm)")
    ax1.set_title(f"{title} - Regression")
    ax1.legend()

    # Right subplot: Residual plot
    ax2.scatter(fitted[~outlier_mask], residuals[~outlier_mask], alpha=0.6, color="C0")
    if np.any(outlier_mask):
        ax2.scatter(
            fitted[outlier_mask],
            residuals[outlier_mask],
            alpha=0.8,
            color="red",
            marker="x",
        )

    ax2.axhline(y=0, color="r", linestyle="--", linewidth=1)
    # Add +/- 3*std lines
    ax2.axhline(y=3 * std_resid, color="gray", linestyle="--", alpha=0.5)
    ax2.axhline(y=-3 * std_resid, color="gray", linestyle="--", alpha=0.5)

    ax2.set_xlabel("Fitted Values (ppm)")
    ax2.set_ylabel("Residuals (ppm)")
    ax2.set_title(f"{title} - Residuals")

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)  # CRITICAL: Prevent memory leaks


def plot_residual_histogram(
    residuals: np.ndarray,
    title: str,
    output_path: Path,
    mae: float | None = None,
    rmsd: float | None = None,
) -> None:
    """Create histogram of residual distribution.

    Args:
        residuals: Array of residual values
        title: Plot title
        output_path: Path to save PNG file
        mae: Mean absolute error (optional, computed if not provided)
        rmsd: Root mean square deviation (optional, computed if not provided)
    """
    if mae is None:
        mae = float(np.mean(np.abs(residuals)))
    if rmsd is None:
        rmsd = float(np.sqrt(np.mean(residuals**2)))

    fig, ax = plt.subplots(figsize=(8, 5))

    ax.hist(residuals, bins=30, alpha=0.7, edgecolor="black")
    ax.axvline(x=np.mean(residuals), color="red", linestyle="--", label="Mean")

    # Add MAE/RMSD text box
    stats_text = f"MAE = {mae:.3f} ppm\nRMSD = {rmsd:.3f} ppm"
    ax.annotate(
        stats_text,
        xy=(0.95, 0.95),
        xycoords="axes fraction",
        fontsize=10,
        verticalalignment="top",
        horizontalalignment="right",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
    )

    ax.set_xlabel("Residual (ppm)")
    ax.set_ylabel("Frequency")
    ax.set_title(title)
    ax.legend()

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)  # CRITICAL: Prevent memory leaks


def calculate_per_compound_stats(
    df: pd.DataFrame,
    factor: ScalingFactor,
    nucleus: str,
) -> pd.DataFrame:
    """Calculate per-compound prediction statistics.

    For each compound, calculates predicted shifts, residuals, and error metrics.

    Args:
        df: DataFrame with columns compound, shielding, exp_shift, nucleus
        factor: ScalingFactor with regression parameters
        nucleus: "1H" or "13C" to filter data

    Returns:
        DataFrame with columns: compound, mean_error, max_error, num_atoms
    """
    # Filter to specified nucleus
    data = df[df["nucleus"] == nucleus].copy()
    if data.empty:
        return pd.DataFrame(columns=["compound", "mean_error", "max_error", "num_atoms"])

    # Calculate predicted shift and residuals
    data["predicted"] = factor.slope * data["shielding"] + factor.intercept
    data["residual"] = data["exp_shift"] - data["predicted"]
    data["abs_error"] = data["residual"].abs()

    # Group by compound and aggregate
    stats = data.groupby("compound").agg(
        mean_error=("abs_error", "mean"),
        max_error=("abs_error", "max"),
        num_atoms=("atom_idx", "count"),
    ).reset_index()

    return stats.sort_values("mean_error", ascending=False)


def save_factors_json(
    factors: dict[str, ScalingFactor],
    output_path: Path,
) -> None:
    """Export scaling factors to JSON for external tools.

    Args:
        factors: Dict mapping factor keys to ScalingFactor objects
        output_path: Path to write JSON file
    """
    # Convert to serializable format
    export_data = {}
    for key, factor in factors.items():
        export_data[key] = {
            "slope": factor.slope,
            "intercept": factor.intercept,
            "r_squared": factor.r_squared,
            "mae": factor.mae,
            "rmsd": factor.rmsd,
            "n_points": factor.n_points,
            "ci_slope": list(factor.ci_slope),
            "ci_intercept": list(factor.ci_intercept),
            "outliers_removed": factor.outliers_removed,
        }

    output_path.write_bytes(
        orjson.dumps(export_data, option=orjson.OPT_INDENT_2)
    )
    logger.info(f"Exported {len(factors)} factors to {output_path}")


def generate_report(output_dir: Path) -> None:
    """Generate publication-quality scaling factor report.

    Creates SCALING-FACTORS.md with methodology, tables, and plot references.
    Generates regression and residual plots for each factor set.

    Args:
        output_dir: Directory to write report and plots
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    plots_dir = output_dir / "plots"
    plots_dir.mkdir(exist_ok=True)

    # Derive all factors
    factors = derive_all_factors()

    if not factors:
        logger.error("No scaling factors derived, cannot generate report")
        return

    # Save factors as JSON
    save_factors_json(factors, output_dir / "scaling_factors.json")

    # Group factors by functional/solvent for iteration
    factor_groups: dict[tuple[str, str], dict[str, ScalingFactor]] = {}
    for key, factor in factors.items():
        parts = key.split("/")
        functional, basis_set, nucleus, solvent = parts
        group_key = (functional, solvent)
        if group_key not in factor_groups:
            factor_groups[group_key] = {}
        factor_groups[group_key][nucleus] = factor

    # Generate plots and collect residuals for each factor
    all_residuals: dict[str, np.ndarray] = {}

    for (functional, solvent), nucleus_factors in factor_groups.items():
        df = aggregate_regression_data(functional, solvent)

        for nucleus, factor in nucleus_factors.items():
            # Calculate residuals for this factor
            data = df[df["nucleus"] == nucleus]
            shielding = data["shielding"].values
            exp_shift = data["exp_shift"].values
            predicted = factor.slope * shielding + factor.intercept
            residuals = exp_shift - predicted
            key = f"{functional}_{solvent}_{nucleus}"
            all_residuals[key] = residuals

            # Generate regression plot
            title = f"{nucleus} {functional} ({solvent})"
            reg_path = plots_dir / f"{functional}_{solvent}_{nucleus}_regression.png"
            plot_regression(df, factor, nucleus, title, reg_path)
            logger.info(f"Generated {reg_path}")

            # Generate residual histogram
            hist_title = f"{nucleus} {functional} ({solvent}) - Residual Distribution"
            hist_path = plots_dir / f"{functional}_{solvent}_{nucleus}_residuals.png"
            plot_residual_histogram(residuals, hist_title, hist_path, factor.mae, factor.rmsd)
            logger.info(f"Generated {hist_path}")

    # Build markdown report
    lines = [
        "# DELTA50 NMR Scaling Factors",
        "",
        "Derived from DELTA50 benchmark data using B3LYP/6-311+G(2d,p) calculations.",
        "",
        "## Methodology",
        "",
        "Scaling factors convert calculated NMR shielding values (sigma) to predicted",
        "chemical shifts (delta) using linear regression:",
        "",
        "```",
        "delta = slope * sigma + intercept",
        "```",
        "",
        "**Fitting procedure:**",
        "1. Aggregate shielding-shift pairs across all DELTA50 compounds",
        "2. Fit ordinary least squares (OLS) regression",
        "3. Identify outliers with residuals > 3 standard deviations",
        "4. Refit without outliers",
        "5. Report final statistics with 95% confidence intervals",
        "",
        "**Training data:** DELTA50 benchmark set (50 small organic molecules with",
        "experimentally assigned NMR spectra in CDCl3).",
        "",
        "## Scaling Factors",
        "",
        "| Nucleus | Solvent | Slope | Intercept | R^2 | MAE (ppm) | RMSD (ppm) | n |",
        "|---------|---------|-------|-----------|-----|-----------|------------|---|",
    ]

    # Add factor table rows
    for key, factor in sorted(factors.items()):
        parts = key.split("/")
        functional, basis_set, nucleus, solvent = parts
        slope_ci = f"{factor.slope:.4f} ({factor.ci_slope[0]:.4f}, {factor.ci_slope[1]:.4f})"
        int_ci = f"{factor.intercept:.2f} ({factor.ci_intercept[0]:.2f}, {factor.ci_intercept[1]:.2f})"
        lines.append(
            f"| {nucleus} | {solvent} | {slope_ci} | {int_ci} | "
            f"{factor.r_squared:.4f} | {factor.mae:.3f} | {factor.rmsd:.3f} | {factor.n_points} |"
        )

    lines.extend([
        "",
        "*Values in parentheses are 95% confidence intervals.*",
        "",
        "## Statistical Summary",
        "",
    ])

    # Per-compound statistics for each factor
    for (functional, solvent), nucleus_factors in factor_groups.items():
        df = aggregate_regression_data(functional, solvent)

        for nucleus, factor in nucleus_factors.items():
            stats = calculate_per_compound_stats(df, factor, nucleus)
            key = f"{functional}_{solvent}_{nucleus}"

            lines.extend([
                f"### {nucleus} {functional} ({solvent})",
                "",
                "| Compound | Mean Error (ppm) | Max Error (ppm) | Atoms |",
                "|----------|------------------|-----------------|-------|",
            ])

            for _, row in stats.iterrows():
                lines.append(
                    f"| {row['compound']} | {row['mean_error']:.3f} | "
                    f"{row['max_error']:.3f} | {int(row['num_atoms'])} |"
                )

            # Flag high-error compounds (>2x MAE)
            high_error = stats[stats["mean_error"] > 2 * factor.mae]
            if not high_error.empty:
                lines.extend([
                    "",
                    f"**Note:** Compounds with mean error > 2x MAE ({2*factor.mae:.3f} ppm):",
                    ", ".join(high_error["compound"].tolist()),
                ])

            lines.append("")

    # Add plot references
    lines.extend([
        "## Plots",
        "",
        "### Regression Analysis",
        "",
    ])

    for (functional, solvent), nucleus_factors in factor_groups.items():
        for nucleus in nucleus_factors:
            reg_file = f"{functional}_{solvent}_{nucleus}_regression.png"
            lines.append(f"![{nucleus} {functional} ({solvent}) Regression](plots/{reg_file})")
            lines.append("")

    lines.extend([
        "### Residual Distributions",
        "",
    ])

    for (functional, solvent), nucleus_factors in factor_groups.items():
        for nucleus in nucleus_factors:
            hist_file = f"{functional}_{solvent}_{nucleus}_residuals.png"
            lines.append(f"![{nucleus} {functional} ({solvent}) Residuals](plots/{hist_file})")
            lines.append("")

    # Usage section
    lines.extend([
        "## Usage",
        "",
        "Apply scaling factors to convert calculated shielding to predicted shift:",
        "",
        "```python",
        "from qm_nmr_calc.benchmark.analysis import derive_all_factors, get_factor_key",
        "",
        "# Get all factors",
        "factors = derive_all_factors()",
        "",
        "# Look up specific factor",
        "key = get_factor_key('B3LYP', '6-311+G(2d,p)', '1H', 'CHCl3')",
        "factor = factors[key]",
        "",
        "# Convert shielding to shift",
        "shielding = 30.5  # ppm (calculated)",
        "shift = factor.slope * shielding + factor.intercept",
        "print(f'Predicted shift: {shift:.2f} ppm')",
        "```",
        "",
        "Or load from JSON:",
        "",
        "```python",
        "import json",
        "",
        "with open('scaling_factors.json') as f:",
        "    factors = json.load(f)",
        "",
        "factor = factors['B3LYP/6-311+G(2d,p)/1H/CHCl3']",
        "shift = factor['slope'] * shielding + factor['intercept']",
        "```",
        "",
        "## Notes",
        "",
        "- **WP04 factors:** Not yet available (benchmark calculations incomplete)",
        "- **DMSO solvent:** Uses same experimental data as CHCl3 (from DELTA50 paper)",
        "  with COSMO solvent model applied during calculation",
        "- **Outlier removal:** Applied 3-sigma threshold; number removed varies by factor set",
        "",
        "---",
        "",
        "*Generated by qm_nmr_calc.benchmark.analysis*",
    ])

    # Write report
    report_path = output_dir / "SCALING-FACTORS.md"
    report_path.write_text("\n".join(lines))
    logger.info(f"Generated report: {report_path}")
