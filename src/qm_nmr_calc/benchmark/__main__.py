"""CLI entry point for DELTA50 benchmark runner.

Usage:
    python -m qm_nmr_calc.benchmark run [--resume] [--headless] [--molecules M1 M2 ...]
    python -m qm_nmr_calc.benchmark status
    python -m qm_nmr_calc.benchmark stop
    python -m qm_nmr_calc.benchmark summary [--output FILE]

Headless execution (for long-running benchmarks):
    nohup python -m qm_nmr_calc.benchmark run --headless > /dev/null 2>&1 &
    python -m qm_nmr_calc.benchmark status   # Check progress
    python -m qm_nmr_calc.benchmark stop     # Request graceful stop
"""

import argparse
import logging
import sys
from pathlib import Path

import orjson

from .runner import aggregate_results, get_results_dir, run_benchmark, show_status


def setup_logging(verbose: bool = False) -> None:
    """Configure logging for CLI."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%H:%M:%S",
    )


def format_status(status: dict) -> str:
    """Format status.json for human-readable display.

    Args:
        status: Parsed status.json content

    Returns:
        Human-readable status string
    """
    lines = [
        "DELTA50 Benchmark Status",
        "=" * 40,
        f"State: {status['state'].upper()}",
    ]

    total = status.get("total_tasks", 0)
    completed = status.get("completed", 0)
    failed = status.get("failed", 0)
    pct = (100 * completed / total) if total > 0 else 0

    lines.append(f"Progress: {completed}/{total} ({pct:.1f}%)")
    lines.append(f"Failed: {failed}")

    # Current task (if running)
    current = status.get("current_task")
    if current:
        lines.append(
            f"Current: {current['molecule_id']} / {current['functional']} / {current['solvent']}"
        )

    # ETA
    eta_hours = status.get("estimated_remaining_hours")
    if eta_hours is not None:
        if eta_hours < 1:
            lines.append(f"ETA: {int(eta_hours * 60)} minutes")
        else:
            lines.append(f"ETA: {eta_hours:.1f} hours")

    # Timestamps
    started = status.get("started_at")
    if started:
        lines.append(f"Started: {started}")

    updated = status.get("updated_at")
    if updated:
        lines.append(f"Last update: {updated}")

    # Recent failures (last 3)
    failures = status.get("failures", [])
    if failures:
        lines.append("")
        lines.append(f"Recent failures ({len(failures)} total):")
        for failure in failures[-3:]:
            mol = failure.get("molecule_id", "unknown")
            func = failure.get("functional", "")
            solv = failure.get("solvent", "")
            err = failure.get("error", "")[:60]  # Truncate long errors
            lines.append(f"  - {mol}/{func}/{solv}: {err}")

    return "\n".join(lines)


def cmd_run(args: argparse.Namespace) -> int:
    """Execute benchmark calculations."""
    setup_logging(args.verbose)

    molecules = args.molecules if args.molecules else None
    functionals = args.functionals if args.functionals else None
    solvents = args.solvents if args.solvents else None

    if args.headless:
        results_dir = get_results_dir()
        print(f"Starting benchmark in headless mode...")
        print(f"Check progress with: python -m qm_nmr_calc.benchmark status")
        print(f"Request stop with: python -m qm_nmr_calc.benchmark stop")
        print(f"Log file: {results_dir / 'benchmark.log'}")

    results, final_state = run_benchmark(
        resume=args.resume,
        molecules=molecules,
        functionals=functionals,
        solvents=solvents,
        processes=args.processes,
        headless=args.headless,
    )

    # Print summary
    completed = sum(1 for r in results if r.status == "complete")
    failed = sum(1 for r in results if r.status == "failed")

    print(f"\nBenchmark run {final_state}:")
    print(f"  Completed: {completed}")
    print(f"  Failed: {failed}")

    if failed > 0:
        print("\nFailed calculations:")
        for r in results:
            if r.status == "failed":
                print(f"  - {r.molecule_id}/{r.functional}/{r.solvent}: {r.error}")

    # Return code based on state
    if final_state == "complete" and failed == 0:
        return 0
    elif final_state == "stopped":
        return 0  # Graceful stop is success
    elif final_state == "paused":
        return 2  # Paused due to threshold
    else:
        return 1  # Failures occurred


def cmd_status(args: argparse.Namespace) -> int:
    """Show benchmark progress."""
    results_dir = get_results_dir()
    status_file = results_dir / "status.json"

    # Try to read status.json (updated by running benchmark)
    if status_file.exists():
        try:
            status = orjson.loads(status_file.read_bytes())
            print(format_status(status))
            print(f"\nResults directory: {results_dir}")
            return 0
        except orjson.JSONDecodeError:
            # File might be mid-write, fall back to counting
            print("(status.json temporarily unavailable, showing file-based count)")

    # Fall back to legacy file-based counting
    show_status()
    return 0


def cmd_stop(args: argparse.Namespace) -> int:
    """Request graceful stop of benchmark run."""
    results_dir = get_results_dir()
    results_dir.mkdir(parents=True, exist_ok=True)
    stop_file = results_dir / "STOP"

    if stop_file.exists():
        print("Stop already requested.")
        print("Runner will stop after current calculation completes.")
    else:
        stop_file.touch()
        print("Stop requested.")
        print("Runner will stop after current calculation completes.")
        print("Check status with: python -m qm_nmr_calc.benchmark status")

    return 0


def cmd_summary(args: argparse.Namespace) -> int:
    """Generate summary CSV."""
    setup_logging(args.verbose)

    output = Path(args.output) if args.output else None
    if output is None:
        output = get_results_dir() / "summary.csv"

    df = aggregate_results(output)
    print(f"Summary generated: {len(df)} results")
    if not df.empty:
        print(df.to_string(index=False))
    return 0


def main() -> int:
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        prog="python -m qm_nmr_calc.benchmark",
        description="DELTA50 benchmark runner for NMR shift validation",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose logging",
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    # run subcommand
    run_parser = subparsers.add_parser(
        "run",
        help="Execute benchmark calculations",
        description=(
            "Execute DELTA50 benchmark calculations. "
            "For long-running benchmarks, use --headless with nohup: "
            "nohup python -m qm_nmr_calc.benchmark run --headless > /dev/null 2>&1 &"
        ),
    )
    run_parser.add_argument(
        "--resume",
        action="store_true",
        default=True,
        help="Skip completed calculations (default: True)",
    )
    run_parser.add_argument(
        "--no-resume",
        action="store_false",
        dest="resume",
        help="Re-run all calculations",
    )
    run_parser.add_argument(
        "--headless",
        action="store_true",
        help="Run in headless mode (disable tqdm, log to file). Use with nohup for background execution.",
    )
    run_parser.add_argument(
        "--molecules",
        nargs="+",
        help="Specific molecule IDs to run (default: all 50)",
    )
    run_parser.add_argument(
        "--functionals",
        nargs="+",
        choices=["B3LYP", "WP04"],
        help="Functionals to test (default: B3LYP WP04)",
    )
    run_parser.add_argument(
        "--solvents",
        nargs="+",
        choices=["CHCl3", "DMSO"],
        help="Solvents to test (default: CHCl3 DMSO)",
    )
    run_parser.add_argument(
        "--processes",
        type=int,
        default=4,
        help="Number of MPI processes per calculation (default: 4)",
    )
    run_parser.set_defaults(func=cmd_run)

    # status subcommand
    status_parser = subparsers.add_parser(
        "status",
        help="Show benchmark progress",
        description="Show current benchmark progress from status.json or by counting completed files.",
    )
    status_parser.set_defaults(func=cmd_status)

    # stop subcommand
    stop_parser = subparsers.add_parser(
        "stop",
        help="Request graceful stop of running benchmark",
        description=(
            "Create STOP marker file to request graceful shutdown. "
            "The runner will stop after completing the current calculation."
        ),
    )
    stop_parser.set_defaults(func=cmd_stop)

    # summary subcommand
    summary_parser = subparsers.add_parser(
        "summary",
        help="Generate summary CSV from results",
    )
    summary_parser.add_argument(
        "--output",
        "-o",
        type=str,
        help="Output CSV file path (default: data/benchmark/results/summary.csv)",
    )
    summary_parser.set_defaults(func=cmd_summary)

    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
