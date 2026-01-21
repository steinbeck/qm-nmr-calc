"""CLI entry point for DELTA50 benchmark runner.

Usage:
    python -m qm_nmr_calc.benchmark run [--resume] [--molecules M1 M2 ...]
    python -m qm_nmr_calc.benchmark status
    python -m qm_nmr_calc.benchmark summary [--output FILE]
"""

import argparse
import logging
import sys
from pathlib import Path

from .runner import aggregate_results, get_results_dir, run_benchmark, show_status


def setup_logging(verbose: bool = False) -> None:
    """Configure logging for CLI."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%H:%M:%S",
    )


def cmd_run(args: argparse.Namespace) -> int:
    """Execute benchmark calculations."""
    setup_logging(args.verbose)

    molecules = args.molecules if args.molecules else None
    functionals = args.functionals if args.functionals else None
    solvents = args.solvents if args.solvents else None

    results = run_benchmark(
        resume=args.resume,
        molecules=molecules,
        functionals=functionals,
        solvents=solvents,
        processes=args.processes,
    )

    # Print summary
    completed = sum(1 for r in results if r.status == "complete")
    failed = sum(1 for r in results if r.status == "failed")

    print(f"\nBenchmark run complete:")
    print(f"  Completed: {completed}")
    print(f"  Failed: {failed}")

    if failed > 0:
        print("\nFailed calculations:")
        for r in results:
            if r.status == "failed":
                print(f"  - {r.molecule_id}/{r.functional}/{r.solvent}: {r.error}")

    return 0 if failed == 0 else 1


def cmd_status(args: argparse.Namespace) -> int:
    """Show benchmark progress."""
    show_status()
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
    )
    status_parser.set_defaults(func=cmd_status)

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
