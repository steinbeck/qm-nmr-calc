"""Pydantic models for DELTA50 benchmark data."""
from pydantic import BaseModel


class ScalingFactor(BaseModel):
    """NMR scaling factor with regression metadata.

    Derived from linear regression of calculated shielding values
    against experimental chemical shifts: shift = slope * shielding + intercept
    """

    slope: float
    intercept: float
    r_squared: float
    mae: float  # Mean absolute error (ppm)
    rmsd: float  # Root mean square deviation (ppm)
    n_points: int  # Number of data points used in regression
    ci_slope: tuple[float, float]  # 95% confidence interval for slope
    ci_intercept: tuple[float, float]  # 95% confidence interval for intercept
    mae_ci: tuple[float, float] | None = None  # Bootstrap CI for MAE
    rmsd_ci: tuple[float, float] | None = None  # Bootstrap CI for RMSD
    outliers_removed: int = 0


class RegressionData(BaseModel):
    """Single data point for regression fitting.

    Pairs a calculated shielding value with its experimental shift.
    """

    compound: str  # e.g., "compound_01"
    atom_idx: int  # Atom index in molecule
    nucleus: str  # "1H" or "13C"
    shielding: float  # Calculated isotropic shielding (ppm)
    exp_shift: float  # Experimental chemical shift (ppm)


class MoleculeData(BaseModel):
    """DELTA50 molecule with experimental shifts."""

    id: str  # e.g., "molecule_01"
    name: str  # Human-readable name
    xyz_file: str  # Relative path from data/benchmark/delta50/molecules/
    h1_shifts: list[float]  # Experimental 1H shifts in ppm
    c13_shifts: list[float]  # Experimental 13C shifts in ppm
    # Optional atom assignments if available
    h1_assignments: dict[str, float] | None = None
    c13_assignments: dict[str, float] | None = None


class ExperimentalShifts(BaseModel):
    """Container for all DELTA50 experimental data."""

    source: dict[str, str]  # Paper citation info
    molecules: dict[str, MoleculeData]


class BenchmarkResult(BaseModel):
    """Single benchmark calculation result."""

    molecule_id: str
    functional: str
    basis_set: str
    solvent: str
    calculated_h1: list[float]  # Calculated 1H shifts
    calculated_c13: list[float]  # Calculated 13C shifts
    h1_mae: float | None = None  # MAE vs experimental (if available)
    c13_mae: float | None = None
    status: str = "complete"  # "complete", "failed", "skipped"
    error: str | None = None
