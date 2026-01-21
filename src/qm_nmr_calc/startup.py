"""Startup validation and recovery logic."""
import sys
from datetime import datetime
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

from .nwchem import validate_nwchem
from .storage import DATA_DIR, list_jobs_by_status, update_job_status


def validate_environment() -> None:
    """
    Validate environment at startup. Exits if validation fails.

    Checks:
    1. NWChem is available and callable
    2. RDKit can initialize and generate 3D coordinates
    3. Data directory is writable
    """
    print("Validating environment...")

    # 1. Check NWChem
    validate_nwchem()
    print("  [OK] NWChem found")

    # 2. Check RDKit
    try:
        mol = Chem.MolFromSmiles("C")  # Simplest molecule
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        print("  [OK] RDKit working")
    except Exception as e:
        sys.exit(f"FATAL: RDKit initialization failed: {e}")

    # 3. Check data directory writable
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    test_file = DATA_DIR / '.write_test'
    try:
        test_file.write_text('test')
        test_file.unlink()
        print(f"  [OK] Data directory writable: {DATA_DIR}")
    except Exception as e:
        sys.exit(f"FATAL: Cannot write to {DATA_DIR}: {e}")

    print("Environment validation passed\n")


def recover_interrupted_jobs() -> int:
    """
    Mark any 'running' jobs as failed on startup.

    This handles the case where the consumer was killed (SIGKILL)
    and SIGNAL_INTERRUPTED didn't fire.

    Returns:
        Number of jobs recovered (marked as failed)
    """
    running_jobs = list_jobs_by_status('running')

    for job_id in running_jobs:
        try:
            update_job_status(
                job_id,
                status='failed',
                completed_at=datetime.utcnow(),
                error_message='Job interrupted - process restart (recovered on startup)'
            )
            print(f"Recovered interrupted job: {job_id}")
        except Exception as e:
            print(f"Warning: Failed to recover job {job_id}: {e}")

    if running_jobs:
        print(f"Recovered {len(running_jobs)} interrupted job(s)\n")

    return len(running_jobs)


def startup() -> None:
    """
    Full startup sequence: validate environment and recover interrupted jobs.
    Call this before starting the consumer.
    """
    validate_environment()
    recover_interrupted_jobs()
