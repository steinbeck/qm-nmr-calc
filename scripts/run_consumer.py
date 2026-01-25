#!/usr/bin/env python
"""
Run the Huey consumer with startup validation.

Usage:
    uv run python scripts/run_consumer.py

Or directly:
    python scripts/run_consumer.py

The consumer will:
1. Validate environment (NWChem, directories)
2. Recover any interrupted jobs from previous runs
3. Start processing queued tasks
"""
import subprocess
import sys
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from qm_nmr_calc.startup import startup
from qm_nmr_calc.nwchem.runner import get_nwchem_version


def main():
    print("=" * 50)
    print("QM NMR Calculation Service - Consumer")
    print("=" * 50 + "\n")

    # Run startup checks
    startup()

    # Print version info
    nwchem_version = get_nwchem_version()
    print(f"NWChem version: {nwchem_version}")
    print()

    # Start Huey consumer
    # -w 1: Single worker (QM jobs are CPU-bound)
    # -k process: Process-based workers (not greenlet)
    print("Starting Huey consumer...")
    print("Press Ctrl+C to stop\n")

    subprocess.run([
        sys.executable, '-m', 'huey.bin.huey_consumer',
        'qm_nmr_calc.queue.huey',
        '-w', '1',
        '-k', 'process',
    ])


if __name__ == '__main__':
    main()
