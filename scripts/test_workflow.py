#!/usr/bin/env python
"""
Integration test for the complete calculation workflow.

This script:
1. Creates a job for a simple molecule (ethanol)
2. Enqueues the task
3. Waits for completion (with timeout)
4. Verifies the result

Usage:
    # First, start the consumer in another terminal:
    uv run python scripts/run_consumer.py

    # Then run this test:
    uv run python scripts/test_workflow.py

For a quick test without full DFT (force-field only):
    uv run python scripts/test_workflow.py --quick
"""
import argparse
import sys
import time
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from qm_nmr_calc import (
    create_job_directory,
    load_job_status,
    get_versions,
    run_optimization_task,
)
from qm_nmr_calc.storage import get_job_dir


def test_quick_workflow():
    """Quick test: just verify job creation and ISiCLE loading."""
    print("\n=== Quick Workflow Test ===\n")

    versions = get_versions()
    print(f"ISiCLE: {versions.isicle}, NWChem: {versions.nwchem}")

    # Test 1: Create job
    smiles = "CCO"  # Ethanol
    job = create_job_directory(smiles, versions.isicle, versions.nwchem)
    print(f"Created job: {job.job_id}")
    print(f"Status: {job.status}")

    # Test 2: Verify directory structure
    job_dir = get_job_dir(job.job_id)
    assert (job_dir / 'status.json').exists(), "status.json missing"
    assert (job_dir / 'output').is_dir(), "output/ missing"
    assert (job_dir / 'logs').is_dir(), "logs/ missing"
    print("Directory structure: OK")

    # Test 3: Load job back
    loaded = load_job_status(job.job_id)
    assert loaded is not None, "Failed to load job"
    assert loaded.input.smiles == smiles, "SMILES mismatch"
    print("Job persistence: OK")

    print("\n[PASS] Quick workflow test passed\n")
    return True


def test_full_workflow(timeout_seconds: int = 600):
    """
    Full test: submit job, run calculation, verify result.

    Requires consumer to be running in another process.
    """
    print("\n=== Full Workflow Test ===\n")
    print(f"Timeout: {timeout_seconds} seconds")
    print("NOTE: Consumer must be running in another terminal\n")

    versions = get_versions()

    # Test molecule: ethanol (small but real)
    smiles = "CCO"

    # Step 1: Create job
    job = create_job_directory(smiles, versions.isicle, versions.nwchem)
    print(f"Created job: {job.job_id}")
    print(f"SMILES: {smiles}")

    # Step 2: Enqueue task
    result = run_optimization_task(job.job_id)
    print(f"Task enqueued")

    # Step 3: Wait for completion
    print(f"\nWaiting for completion...")
    start_time = time.time()

    while True:
        status = load_job_status(job.job_id)

        if status.status == 'complete':
            elapsed = time.time() - start_time
            print(f"\n[COMPLETE] Job finished in {elapsed:.1f}s")
            break

        if status.status == 'failed':
            print(f"\n[FAILED] Job failed: {status.error_message}")
            if status.error_traceback:
                print(f"\nTraceback:\n{status.error_traceback}")
            return False

        elapsed = time.time() - start_time
        if elapsed > timeout_seconds:
            print(f"\n[TIMEOUT] Job did not complete within {timeout_seconds}s")
            print(f"Current status: {status.status}")
            return False

        # Print progress every 30 seconds
        if int(elapsed) % 30 == 0 and int(elapsed) > 0:
            print(f"  ... still {status.status} ({int(elapsed)}s)")

        time.sleep(5)

    # Step 4: Verify result
    job_dir = get_job_dir(job.job_id)
    output_file = job_dir / 'output' / 'optimized.xyz'

    if not output_file.exists():
        print(f"[FAIL] Output file missing: {output_file}")
        return False

    # Read and display output
    content = output_file.read_text()
    lines = content.strip().split('\n')
    n_atoms = int(lines[0]) if lines else 0
    print(f"\nOutput file: {output_file}")
    print(f"Atoms: {n_atoms}")
    print(f"First 5 lines:\n{chr(10).join(lines[:5])}")

    # Verify expected atoms for ethanol (C2H6O = 9 atoms)
    if n_atoms != 9:
        print(f"\n[WARN] Expected 9 atoms for ethanol, got {n_atoms}")

    print("\n[PASS] Full workflow test passed\n")
    return True


def main():
    parser = argparse.ArgumentParser(description='Test calculation workflow')
    parser.add_argument('--quick', action='store_true',
                        help='Quick test (no DFT calculation)')
    parser.add_argument('--timeout', type=int, default=600,
                        help='Timeout in seconds for full test (default: 600)')
    args = parser.parse_args()

    if args.quick:
        success = test_quick_workflow()
    else:
        success = test_full_workflow(args.timeout)

    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()
