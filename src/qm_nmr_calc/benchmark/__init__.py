"""DELTA50 benchmark infrastructure for NMR shift validation."""

from .analysis import (
    aggregate_regression_data,
    calculate_per_compound_stats,
    derive_all_factors,
    fit_scaling_factors,
    generate_report,
    get_factor_key,
    plot_regression,
    plot_residual_histogram,
    save_factors_json,
)
from .data_loader import get_data_dir, load_delta50_molecules, load_experimental_shifts
from .models import BenchmarkResult, ExperimentalShifts, MoleculeData, RegressionData, ScalingFactor
from .runner import (
    FAILURE_THRESHOLD,
    aggregate_results,
    build_task_matrix,
    check_stop_requested,
    clear_stop_file,
    get_results_dir,
    run_benchmark,
    show_status,
    update_status,
)

__all__ = [
    # Models
    "MoleculeData",
    "ExperimentalShifts",
    "BenchmarkResult",
    "ScalingFactor",
    "RegressionData",
    # Data loading
    "load_delta50_molecules",
    "load_experimental_shifts",
    "get_data_dir",
    # Runner
    "run_benchmark",
    "build_task_matrix",
    "show_status",
    "aggregate_results",
    "get_results_dir",
    # Status and control
    "update_status",
    "check_stop_requested",
    "clear_stop_file",
    "FAILURE_THRESHOLD",
    # Analysis
    "derive_all_factors",
    "fit_scaling_factors",
    "aggregate_regression_data",
    "get_factor_key",
    "generate_report",
    "calculate_per_compound_stats",
    "plot_regression",
    "plot_residual_histogram",
    "save_factors_json",
]
