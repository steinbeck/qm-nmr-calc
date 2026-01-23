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
        solvents: List of solvents to process (default: ["CHCl3", "DMSO"])

    Returns:
        Dict mapping "functional/basis_set/nucleus/solvent" to ScalingFactor
    """
    if functionals is None:
        functionals = ["B3LYP"]  # WP04 not yet complete
    if solvents is None:
        solvents = ["CHCl3", "DMSO"]

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
