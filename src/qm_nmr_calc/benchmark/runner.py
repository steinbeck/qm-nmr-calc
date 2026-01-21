"""Benchmark runner for DELTA50 calculation matrix."""

import json
import logging
from pathlib import Path

from tqdm import tqdm

from qm_nmr_calc.nwchem import run_calculation
from qm_nmr_calc.shifts import shielding_to_shift

from .data_loader import get_data_dir, load_experimental_shifts
from .models import BenchmarkResult

logger = logging.getLogger(__name__)

# Benchmark configuration
FUNCTIONALS = ["B3LYP", "WP04"]
SOLVENTS = ["CHCl3", "DMSO"]

# Preset configurations for benchmark (not using existing presets.py - need WP04)
BENCHMARK_PRESETS = {
    "B3LYP": {
        "functional": "b3lyp",
        "basis_set": "6-31G*",
        "nmr_basis_set": "6-311+G(2d,p)",
        "max_iter": 150,
    },
    "WP04": {
        "functional": "wp04",  # Optimized for 1H shifts
        "basis_set": "6-31G*",
        "nmr_basis_set": "6-311++G(2d,p)",  # DELTA50 used this for WP04
        "max_iter": 150,
    },
}


def get_results_dir() -> Path:
    """Get path to benchmark results directory."""
    return get_data_dir().parent / "results"


def build_task_matrix(
    molecules: list[str] | None = None,
    functionals: list[str] | None = None,
    solvents: list[str] | None = None,
) -> list[dict]:
    """Build list of all benchmark tasks.

    Args:
        molecules: Specific molecule IDs to include (default: all 50)
        functionals: Functionals to test (default: B3LYP, WP04)
        solvents: Solvents to test (default: CHCl3, DMSO)

    Returns:
        List of task dicts with keys: molecule_id, functional, solvent, task_id
    """
    if functionals is None:
        functionals = FUNCTIONALS
    if solvents is None:
        solvents = SOLVENTS

    # Load molecule list if not specified
    if molecules is None:
        exp_shifts = load_experimental_shifts()
        molecules = sorted(exp_shifts.molecules.keys())

    tasks = []
    for mol_id in molecules:
        for func in functionals:
            for solv in solvents:
                task_id = f"{mol_id}_{func}_{solv}"
                tasks.append(
                    {
                        "task_id": task_id,
                        "molecule_id": mol_id,
                        "functional": func,
                        "solvent": solv,
                    }
                )

    return tasks


def is_task_complete(results_dir: Path, task: dict) -> bool:
    """Check if a task has already completed successfully.

    Checks for shifts.json marker file in the output directory.
    """
    output_dir = (
        results_dir / task["molecule_id"] / f"{task['functional']}_{task['solvent']}"
    )
    shifts_file = output_dir / "shifts.json"
    return shifts_file.exists()


def append_to_jsonl(filepath: Path, data: dict) -> None:
    """Append a result to JSONL progress file."""
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with filepath.open("a") as f:
        f.write(json.dumps(data) + "\n")


def run_single_calculation(
    mol_id: str,
    functional: str,
    solvent: str,
    results_dir: Path,
    processes: int = 4,
) -> BenchmarkResult:
    """Run a single benchmark calculation.

    Args:
        mol_id: Molecule ID (e.g., "compound_01")
        functional: Functional name (B3LYP or WP04)
        solvent: Solvent name (CHCl3 or DMSO)
        results_dir: Base results directory
        processes: Number of MPI processes

    Returns:
        BenchmarkResult with calculated shifts or error info
    """
    output_dir = results_dir / mol_id / f"{functional}_{solvent}"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load molecule geometry
    xyz_dir = get_data_dir() / "molecules"
    xyz_file = xyz_dir / f"{mol_id}.xyz"

    if not xyz_file.exists():
        return BenchmarkResult(
            molecule_id=mol_id,
            functional=functional,
            basis_set=BENCHMARK_PRESETS[functional]["nmr_basis_set"],
            solvent=solvent,
            calculated_h1=[],
            calculated_c13=[],
            status="failed",
            error=f"XYZ file not found: {xyz_file}",
        )

    try:
        # Get preset for this functional
        preset = BENCHMARK_PRESETS[functional]

        # Run NWChem calculation
        # Use skip_optimization=True with geometry_file to use pre-optimized coords
        result = run_calculation(
            smiles="",  # Not used when skip_optimization=True
            job_dir=output_dir,
            preset=preset,
            solvent=solvent.lower(),  # COSMO expects lowercase
            processes=processes,
            skip_optimization=True,
            geometry_file=xyz_file,
        )

        # Convert shielding to shifts
        # Use production scaling for now (will derive benchmark-specific later)
        shifts_data = shielding_to_shift(result["shielding_data"], preset="production")

        # Extract shift values
        h1_shifts = [s["shift"] for s in shifts_data["1H"]]
        c13_shifts = [s["shift"] for s in shifts_data["13C"]]

        # Save shifts to JSON
        shifts_output = {
            "molecule_id": mol_id,
            "functional": functional,
            "solvent": solvent,
            "h1_shifts": h1_shifts,
            "c13_shifts": c13_shifts,
            "shielding_data": result["shielding_data"],
        }
        shifts_file = output_dir / "shifts.json"
        with shifts_file.open("w") as f:
            json.dump(shifts_output, f, indent=2)

        return BenchmarkResult(
            molecule_id=mol_id,
            functional=functional,
            basis_set=preset["nmr_basis_set"],
            solvent=solvent,
            calculated_h1=h1_shifts,
            calculated_c13=c13_shifts,
            status="complete",
        )

    except Exception as e:
        logger.error(f"Calculation failed for {mol_id}/{functional}/{solvent}: {e}")
        return BenchmarkResult(
            molecule_id=mol_id,
            functional=functional,
            basis_set=BENCHMARK_PRESETS[functional]["nmr_basis_set"],
            solvent=solvent,
            calculated_h1=[],
            calculated_c13=[],
            status="failed",
            error=str(e),
        )


def run_benchmark(
    resume: bool = True,
    molecules: list[str] | None = None,
    functionals: list[str] | None = None,
    solvents: list[str] | None = None,
    processes: int = 4,
) -> list[BenchmarkResult]:
    """Execute DELTA50 benchmark calculation matrix.

    Args:
        resume: If True, skip already-completed calculations
        molecules: Specific molecule IDs to run (default: all 50)
        functionals: Functionals to test (default: B3LYP, WP04)
        solvents: Solvents to test (default: CHCl3, DMSO)
        processes: Number of MPI processes per calculation

    Returns:
        List of BenchmarkResult for all executed tasks
    """
    results_dir = get_results_dir()
    results_dir.mkdir(parents=True, exist_ok=True)
    progress_file = results_dir / "progress.jsonl"

    # Build task matrix
    tasks = build_task_matrix(molecules, functionals, solvents)

    # Filter if resuming
    if resume:
        tasks = [t for t in tasks if not is_task_complete(results_dir, t)]

    if not tasks:
        logger.info("No tasks to run (all complete or empty task list)")
        return []

    logger.info(f"Running {len(tasks)} benchmark calculations")

    results = []
    for task in tqdm(tasks, desc="DELTA50 Benchmark"):
        result = run_single_calculation(
            mol_id=task["molecule_id"],
            functional=task["functional"],
            solvent=task["solvent"],
            results_dir=results_dir,
            processes=processes,
        )
        results.append(result)

        # Log progress to JSONL
        append_to_jsonl(
            progress_file,
            {
                "task_id": task["task_id"],
                "status": result.status,
                "error": result.error,
            },
        )

    return results


def aggregate_results(output_file: Path | None = None) -> "pd.DataFrame":
    """Aggregate all benchmark results into summary DataFrame.

    Args:
        output_file: If provided, write CSV to this path

    Returns:
        DataFrame with all results
    """
    import pandas as pd

    results_dir = get_results_dir()
    records = []

    for mol_dir in sorted(results_dir.glob("compound_*")):
        for method_dir in mol_dir.glob("*_*"):
            shifts_file = method_dir / "shifts.json"
            if not shifts_file.exists():
                continue

            with shifts_file.open() as f:
                data = json.load(f)

            # Parse directory names
            molecule_id = mol_dir.name
            parts = method_dir.name.split("_", 1)
            if len(parts) == 2:
                functional, solvent = parts
            else:
                continue

            records.append(
                {
                    "molecule_id": molecule_id,
                    "functional": functional,
                    "solvent": solvent,
                    "num_h1": len(data.get("h1_shifts", [])),
                    "num_c13": len(data.get("c13_shifts", [])),
                }
            )

    df = pd.DataFrame(records)

    if output_file:
        output_file.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_file, index=False)
        logger.info(f"Summary written to {output_file} ({len(records)} results)")

    return df


def show_status() -> None:
    """Print benchmark progress status."""
    results_dir = get_results_dir()

    # Count completed tasks
    total_tasks = 50 * len(FUNCTIONALS) * len(SOLVENTS)  # 200

    completed = 0
    for mol_dir in results_dir.glob("compound_*"):
        for method_dir in mol_dir.glob("*_*"):
            if (method_dir / "shifts.json").exists():
                completed += 1

    print("DELTA50 Benchmark Status")
    print("========================")
    print(f"Completed: {completed}/{total_tasks} ({100*completed/total_tasks:.1f}%)")
    print(f"Remaining: {total_tasks - completed}")
    print(f"Results directory: {results_dir}")
