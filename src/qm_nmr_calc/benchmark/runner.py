"""Benchmark runner for DELTA50 calculation matrix."""

import logging
from collections import deque
from datetime import datetime, timezone
from pathlib import Path

import orjson
from tqdm import tqdm

from qm_nmr_calc.nwchem import run_calculation
from qm_nmr_calc.shifts import shielding_to_shift

from .data_loader import get_data_dir, load_experimental_shifts
from .models import BenchmarkResult

logger = logging.getLogger(__name__)

# Failure threshold: stop if >10% of molecules fail (5 unique molecules)
FAILURE_THRESHOLD = 5

# Benchmark configuration
FUNCTIONALS = ["B3LYP", "WP04"]
SOLVENTS = ["CHCl3", "DMSO", "Methanol", "Water", "Acetone", "Benzene"]

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


def check_stop_requested(results_dir: Path) -> bool:
    """Check if user requested graceful stop via marker file.

    User can create STOP file to request graceful shutdown between calculations.
    """
    stop_file = results_dir / "STOP"
    return stop_file.exists()


def clear_stop_file(results_dir: Path) -> None:
    """Remove STOP file after acknowledging stop request."""
    stop_file = results_dir / "STOP"
    if stop_file.exists():
        stop_file.unlink()


def update_status(
    results_dir: Path,
    state: str,
    total_tasks: int,
    completed: int,
    failed: int,
    current_task: dict | None,
    failures: list[dict],
    started_at: datetime,
    calc_times: deque[float],
    run_id: str | None = None,
) -> None:
    """Update global benchmark status file atomically.

    Args:
        results_dir: Base results directory
        state: Current state (running/stopped/complete/failed/paused)
        total_tasks: Total number of tasks to run
        completed: Number of completed tasks
        failed: Number of failed tasks
        current_task: Current task info (molecule_id, functional, solvent, started_at)
        failures: List of failure records (last 10 kept)
        started_at: When benchmark started
        calc_times: Rolling deque of recent calculation times in seconds
        run_id: Unique run identifier (defaults to started_at ISO string)
    """
    now = datetime.now(timezone.utc)
    status_file = results_dir / "status.json"

    # Calculate ETA from rolling average of last N calculation times
    avg_time: float | None = None
    eta_hours: float | None = None
    if calc_times:
        avg_time = sum(calc_times) / len(calc_times)
        remaining = total_tasks - completed - failed
        eta_seconds = remaining * avg_time
        eta_hours = eta_seconds / 3600 if eta_seconds > 0 else 0.0

    status = {
        "run_id": run_id or started_at.isoformat(),
        "started_at": started_at.isoformat(),
        "updated_at": now.isoformat(),
        "state": state,
        "total_tasks": total_tasks,
        "completed": completed,
        "failed": failed,
        "current_task": current_task,
        "failures": list(failures[-10:]),  # Keep last 10 failures
        "estimated_remaining_hours": round(eta_hours, 2) if eta_hours is not None else None,
        "avg_calc_time_seconds": round(avg_time, 1) if avg_time is not None else None,
    }

    # Atomic write: temp file then rename
    temp_file = status_file.with_suffix(".json.tmp")
    temp_file.write_bytes(orjson.dumps(status, option=orjson.OPT_INDENT_2))
    temp_file.rename(status_file)


def build_task_matrix(
    molecules: list[str] | None = None,
    functionals: list[str] | None = None,
    solvents: list[str] | None = None,
) -> list[dict]:
    """Build list of all benchmark tasks.

    Args:
        molecules: Specific molecule IDs to include (default: all 50)
        functionals: Functionals to test (default: B3LYP, WP04)
        solvents: Solvents to test (default: CHCl3, DMSO). Can include "vacuum"
            for gas-phase calculations.

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
    with filepath.open("ab") as f:
        f.write(orjson.dumps(data) + b"\n")


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
        solvent: Solvent name (CHCl3, DMSO, or vacuum)
        results_dir: Base results directory
        processes: Number of MPI processes

    Returns:
        BenchmarkResult with calculated shifts or error info.
        For solvents without scaling factors (e.g., vacuum), shifts will be empty
        but shielding data is still saved.
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

        # Convert shielding to shifts using DELTA50 regression factors
        # For solvents without scaling factors (e.g., vacuum), save shielding only
        try:
            shifts_data = shielding_to_shift(
                result["shielding_data"],
                functional=functional.upper(),
                basis_set=BENCHMARK_PRESETS[functional]["nmr_basis_set"],
                solvent=solvent,
            )
            h1_shifts = [s["shift"] for s in shifts_data["1H"]]
            c13_shifts = [s["shift"] for s in shifts_data["13C"]]
        except ValueError as e:
            # No scaling factor for this combination - save raw shielding only
            # This is expected for vacuum until factors are derived from this data
            logger.info(f"No scaling factor for {functional}/{solvent}: {e}")
            h1_shifts = []
            c13_shifts = []

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
        shifts_file.write_bytes(orjson.dumps(shifts_output, option=orjson.OPT_INDENT_2))

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
    headless: bool = False,
) -> tuple[list[BenchmarkResult], str]:
    """Execute DELTA50 benchmark calculation matrix.

    Args:
        resume: If True, skip already-completed calculations
        molecules: Specific molecule IDs to run (default: all 50)
        functionals: Functionals to test (default: B3LYP, WP04)
        solvents: Solvents to test (default: CHCl3, DMSO). Can include "vacuum"
            for gas-phase calculations without COSMO solvation.
        processes: Number of MPI processes per calculation
        headless: If True, disable tqdm and configure file logging

    Returns:
        Tuple of (results list, final state string)
        State is one of: "complete", "stopped", "paused", "failed"
    """
    results_dir = get_results_dir()
    results_dir.mkdir(parents=True, exist_ok=True)
    progress_file = results_dir / "progress.jsonl"

    # Configure logging for headless mode
    if headless:
        log_file = results_dir / "benchmark.log"
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(
            logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
        )
        logger.addHandler(file_handler)
        logger.info("Benchmark started in headless mode")

    # Build task matrix
    tasks = build_task_matrix(molecules, functionals, solvents)
    total_tasks = len(tasks)

    # Filter if resuming
    if resume:
        tasks = [t for t in tasks if not is_task_complete(results_dir, t)]

    if not tasks:
        logger.info("No tasks to run (all complete or empty task list)")
        return [], "complete"

    logger.info(f"Running {len(tasks)} benchmark calculations")

    # Initialize tracking
    started_at = datetime.now(timezone.utc)
    results: list[BenchmarkResult] = []
    failures: list[dict] = []
    failed_molecules: set[str] = set()  # Track unique failed molecules
    calc_times: deque[float] = deque(maxlen=10)  # Rolling average of last 10
    final_state = "complete"

    # Clear any existing STOP file from previous run
    clear_stop_file(results_dir)

    # Write initial status
    update_status(
        results_dir=results_dir,
        state="running",
        total_tasks=total_tasks,
        completed=len(results),
        failed=len(failures),
        current_task=None,
        failures=failures,
        started_at=started_at,
        calc_times=calc_times,
    )

    for task in tqdm(tasks, desc="DELTA50 Benchmark", disable=headless):
        # Check for stop request between calculations
        if check_stop_requested(results_dir):
            logger.info("Stop requested, finishing gracefully")
            final_state = "stopped"
            update_status(
                results_dir=results_dir,
                state="stopped",
                total_tasks=total_tasks,
                completed=len([r for r in results if r.status == "complete"]),
                failed=len(failures),
                current_task=None,
                failures=failures,
                started_at=started_at,
                calc_times=calc_times,
            )
            clear_stop_file(results_dir)
            return results, final_state

        # Record task start
        task_started = datetime.now(timezone.utc)
        current_task = {
            "molecule_id": task["molecule_id"],
            "functional": task["functional"],
            "solvent": task["solvent"],
            "started_at": task_started.isoformat(),
        }

        # Update status to show current task
        update_status(
            results_dir=results_dir,
            state="running",
            total_tasks=total_tasks,
            completed=len([r for r in results if r.status == "complete"]),
            failed=len(failures),
            current_task=current_task,
            failures=failures,
            started_at=started_at,
            calc_times=calc_times,
        )

        # Run calculation
        result = run_single_calculation(
            mol_id=task["molecule_id"],
            functional=task["functional"],
            solvent=task["solvent"],
            results_dir=results_dir,
            processes=processes,
        )
        results.append(result)

        # Record calculation time
        task_ended = datetime.now(timezone.utc)
        calc_time = (task_ended - task_started).total_seconds()
        calc_times.append(calc_time)

        # Track failures
        if result.status == "failed":
            failure_record = {
                "molecule_id": task["molecule_id"],
                "functional": task["functional"],
                "solvent": task["solvent"],
                "error": result.error,
                "timestamp": task_ended.isoformat(),
            }
            failures.append(failure_record)
            failed_molecules.add(task["molecule_id"])

            # Check failure threshold (unique molecules)
            if len(failed_molecules) > FAILURE_THRESHOLD:
                logger.warning(
                    f"Failure threshold exceeded ({len(failed_molecules)} unique molecules failed). Pausing."
                )
                final_state = "paused"
                update_status(
                    results_dir=results_dir,
                    state="paused",
                    total_tasks=total_tasks,
                    completed=len([r for r in results if r.status == "complete"]),
                    failed=len(failures),
                    current_task=None,
                    failures=failures,
                    started_at=started_at,
                    calc_times=calc_times,
                )
                return results, final_state

        # Log progress to JSONL
        append_to_jsonl(
            progress_file,
            {
                "task_id": task["task_id"],
                "status": result.status,
                "error": result.error,
                "duration_seconds": calc_time,
            },
        )

        # Update status after completion
        if headless:
            logger.info(
                f"Completed {task['task_id']}: {result.status} ({calc_time:.1f}s)"
            )

    # Final status update
    update_status(
        results_dir=results_dir,
        state="complete",
        total_tasks=total_tasks,
        completed=len([r for r in results if r.status == "complete"]),
        failed=len(failures),
        current_task=None,
        failures=failures,
        started_at=started_at,
        calc_times=calc_times,
    )

    return results, final_state


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

            data = orjson.loads(shifts_file.read_bytes())

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
    total_tasks = 50 * len(FUNCTIONALS) * len(SOLVENTS)

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
