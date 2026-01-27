"""NWChem execution and calculation orchestration.

This module provides the main entry point for running NWChem calculations,
replacing the ISiCLE wrapper with direct subprocess execution and our own
input generation/output parsing.
"""

import os
import subprocess
import sys
from pathlib import Path

from .geometry import smiles_to_xyz, load_geometry_file, mol_to_xyz_block
from .input_gen import generate_optimization_input, generate_shielding_input
from .output_parser import extract_optimized_geometry, parse_shielding_output


def validate_nwchem() -> None:
    """Validate NWChem is available and callable.

    Call at startup - exits if validation fails.

    Raises:
        SystemExit: If NWChem is not found or validation times out
    """
    try:
        result = subprocess.run(
            ["which", "nwchem"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        if result.returncode != 0:
            sys.exit("FATAL: nwchem not found in PATH. Ensure NWChem is installed.")
    except subprocess.TimeoutExpired:
        sys.exit("FATAL: 'which nwchem' timed out")
    except Exception as e:
        sys.exit(f"FATAL: Cannot verify nwchem: {e}")


def get_nwchem_version() -> str:
    """Get NWChem version string.

    Note: NWChem doesn't have a --version flag, so we hardcode the known version.

    Returns:
        Version string (e.g., "7.0.2")

    TODO: Parse dynamically from NWChem output header if needed.
    """
    return "7.0.2"


def run_nwchem(
    input_file: Path,
    output_file: Path,
    processes: int = 4,
) -> None:
    """Execute NWChem via subprocess with MPI.

    Uses shell execution like ISiCLE for compatibility.

    Args:
        input_file: Path to NWChem input file (.nw)
        output_file: Path to write NWChem output
        processes: Number of MPI processes

    Raises:
        RuntimeError: If NWChem returns non-zero exit code
    """
    # Match ISiCLE's approach: shell=True, --bind-to none, unset DISPLAY
    infile = str(input_file.resolve())
    outfile = str(output_file.resolve())
    logfile = str(output_file.with_suffix(".log").resolve())

    # Unset DISPLAY to prevent X11 errors in headless mode
    cmd = f"unset DISPLAY; mpirun --bind-to none -n {processes} nwchem {infile} > {outfile} 2> {logfile}"

    result = subprocess.call(cmd, shell=True)

    if result != 0:
        error_msg = f"NWChem failed with exit code {result}"
        # Read stderr from log file
        if Path(logfile).exists():
            stderr = Path(logfile).read_text()
            if stderr.strip():
                error_msg += f"\nstderr: {stderr}"
        # Include last 50 lines of stdout for debugging
        if Path(outfile).exists():
            stdout = Path(outfile).read_text()
            stdout_lines = stdout.strip().split("\n")
            if stdout_lines:
                error_msg += f"\nLast output lines:\n" + "\n".join(stdout_lines[-50:])
        raise RuntimeError(error_msg)


def run_calculation(
    smiles: str | None,
    job_dir: Path,
    preset: dict,
    solvent: str,
    processes: int = 4,
    skip_optimization: bool = False,
    geometry_file: Path | None = None,
    on_optimization_complete: callable = None,
    scratch_dir_override: Path | None = None,
) -> dict:
    """Run complete NMR calculation (geometry optimization + NMR shielding).

    This is the main entry point replacing isicle_wrapper functions.
    Runs a two-step DFT calculation:
    1. Geometry optimization (with COSMO solvation if solvent specified)
    2. NMR shielding calculation (with COSMO solvation if solvent specified)

    For vacuum/gas-phase calculations, pass solvent="vacuum" to run without COSMO.

    Args:
        smiles: SMILES string of molecule (None if geometry_file provided)
        job_dir: Job directory for outputs and scratch
        preset: Calculation preset dict with functional, basis_set, nmr_basis_set, max_iter
        solvent: Solvent name for COSMO (chcl3, dmso) or "vacuum" for gas-phase
        processes: Number of MPI processes
        skip_optimization: If True, skip geometry optimization and use geometry_file
        geometry_file: Path to pre-optimized geometry file (required if skip_optimization=True or smiles=None)
        on_optimization_complete: Optional callback after optimization completes
        scratch_dir_override: Optional custom scratch directory (for per-conformer isolation)

    Returns:
        Dict with:
            - 'geometry_file': Path to optimized.xyz
            - 'shielding_data': dict in shifts.py expected format
            - 'optimization_output': Path to optimize.out (or None if skipped)
            - 'shielding_output': Path to shielding.out

    Raises:
        ValueError: If both smiles and geometry_file are None, or if skip_optimization=True but geometry_file not provided
        RuntimeError: If NWChem calculation fails
    """
    # Validate inputs
    if smiles is None and geometry_file is None:
        raise ValueError("Either smiles or geometry_file must be provided")

    # Set up directories
    scratch_dir = scratch_dir_override if scratch_dir_override else job_dir / "scratch"
    output_dir = job_dir / "output"
    scratch_dir.mkdir(exist_ok=True)
    output_dir.mkdir(exist_ok=True)

    optimization_output = None

    if skip_optimization:
        # Load pre-optimized geometry
        if geometry_file is None:
            raise ValueError(
                "geometry_file is required when skip_optimization=True"
            )
        mol, xyz_block = load_geometry_file(geometry_file)
        optimized_xyz_path = output_dir / "optimized.xyz"
        # Copy geometry file to output dir if different location
        if geometry_file.resolve() != optimized_xyz_path.resolve():
            optimized_xyz_path.write_text(geometry_file.read_text())
    else:
        # Step 1: Generate initial 3D geometry
        if smiles is not None:
            # Generate from SMILES
            mol, xyz_block = smiles_to_xyz(smiles)
        else:
            # Load from geometry file
            if geometry_file is None:
                raise ValueError("geometry_file required when smiles is None and skip_optimization=False")
            mol, xyz_block = load_geometry_file(geometry_file)

        # Step 2: Generate optimization input with COSMO
        opt_input = generate_optimization_input(
            geometry_xyz=xyz_block,
            functional=preset["functional"],
            basis_set=preset["basis_set"],
            solvent=solvent,
            max_iter=preset.get("max_iter", 150),
        )

        # Write optimization input
        opt_input_file = scratch_dir / "optimize.nw"
        opt_input_file.write_text(opt_input)

        # Step 3: Run NWChem geometry optimization
        opt_output_file = scratch_dir / "optimize.out"
        run_nwchem(opt_input_file, opt_output_file, processes)
        optimization_output = opt_output_file

        # Step 4: Extract optimized geometry
        opt_output_text = opt_output_file.read_text()
        xyz_block = extract_optimized_geometry(opt_output_text)

        # Save optimized geometry
        optimized_xyz_path = output_dir / "optimized.xyz"
        # Write full XYZ format with atom count and title
        num_atoms = len(xyz_block.strip().split("\n"))
        full_xyz = f"{num_atoms}\nOptimized geometry from NWChem\n{xyz_block}"
        optimized_xyz_path.write_text(full_xyz)

    # Notify caller that optimization is complete (for step tracking)
    if on_optimization_complete:
        on_optimization_complete()

    # Step 5: Generate NMR shielding input with COSMO
    shielding_input = generate_shielding_input(
        geometry_xyz=xyz_block,
        functional=preset["functional"],
        basis_set=preset["nmr_basis_set"],
        solvent=solvent,
    )

    # Write shielding input
    shielding_input_file = scratch_dir / "shielding.nw"
    shielding_input_file.write_text(shielding_input)

    # Step 6: Run NWChem NMR calculation
    shielding_output_file = scratch_dir / "shielding.out"
    run_nwchem(shielding_input_file, shielding_output_file, processes)

    # Step 7: Parse shielding output
    shielding_output_text = shielding_output_file.read_text()
    shielding_data = parse_shielding_output(shielding_output_text)

    return {
        "geometry_file": optimized_xyz_path,
        "shielding_data": shielding_data,
        "optimization_output": optimization_output,
        "shielding_output": shielding_output_file,
    }


def run_conformer_dft_optimization(
    ensemble,
    job_id: str,
    preset: dict,
    solvent: str,
    processes: int = 4,
):
    """Run DFT geometry optimization on all conformers in an ensemble.

    Processes conformers sequentially. Individual conformer failures are caught
    and do not stop processing of remaining conformers. However, if more than
    50% of conformers fail, raises RuntimeError.

    Args:
        ensemble: ConformerEnsemble with conformers to optimize
        job_id: Job identifier for directory lookup
        preset: Calculation preset dict (functional, basis_set, nmr_basis_set, max_iter)
        solvent: Solvent name for COSMO
        processes: Number of MPI processes

    Returns:
        Tuple of (successful_conformers, failed_conformers) where each is a list of ConformerData

    Raises:
        RuntimeError: If less than 50% of conformers succeed
    """
    from qm_nmr_calc.models import ConformerData, ConformerEnsemble
    from qm_nmr_calc.storage import get_job_dir, get_conformer_scratch_dir
    from qm_nmr_calc.nwchem.output_parser import extract_dft_energy

    job_dir = get_job_dir(job_id)
    successful = []
    failed = []

    for conformer in ensemble.conformers:
        try:
            # Update status to optimizing
            conformer.status = "optimizing"

            # Get conformer-specific scratch directory
            scratch_dir = get_conformer_scratch_dir(job_id, conformer.conformer_id)

            # Resolve geometry file path
            geom_path = job_dir / conformer.geometry_file

            # Run DFT optimization (no SMILES, use geometry file)
            result = run_calculation(
                smiles=None,
                job_dir=job_dir,
                preset=preset,
                solvent=solvent,
                processes=processes,
                skip_optimization=False,
                geometry_file=geom_path,
                scratch_dir_override=scratch_dir,
            )

            # Extract DFT energy from optimization output
            opt_output_text = result["optimization_output"].read_text()
            dft_energy = extract_dft_energy(opt_output_text)

            # Update conformer with DFT results
            conformer.energy = dft_energy
            conformer.energy_unit = "hartree"
            conformer.status = "optimized"
            conformer.optimized_geometry_file = str(
                result["geometry_file"].relative_to(job_dir)
            )

            successful.append(conformer)

        except Exception as e:
            # Mark conformer as failed and continue
            import traceback
            conformer.status = "failed"
            error_msg = f"{type(e).__name__}: {str(e)}"
            conformer.error_message = error_msg[:200]
            failed.append(conformer)
            # DEBUG: Print traceback for test debugging
            # print(f"DEBUG: Conformer {conformer.conformer_id} failed: {error_msg}")
            # traceback.print_exc()

    # Check success rate
    success_rate = len(successful) / len(ensemble.conformers)
    if success_rate < 0.5:
        raise RuntimeError(
            f"DFT optimization failed for too many conformers: "
            f"{len(successful)}/{len(ensemble.conformers)} succeeded "
            f"({success_rate:.1%} success rate). "
            f"Need at least 50% success rate."
        )

    return (successful, failed)


def apply_post_dft_filter(
    optimized_conformers: list,
    window_kcal: float = 3.0,
) -> list:
    """Filter conformers by DFT energy window.

    Keeps only conformers within window_kcal of the lowest DFT energy.
    Uses existing filter_by_energy_window from conformers.filters module.

    Args:
        optimized_conformers: List of ConformerData with DFT energies in Hartree
        window_kcal: Energy window in kcal/mol above minimum (default: 3.0)

    Returns:
        List of ConformerData within energy window
    """
    from qm_nmr_calc.conformers.boltzmann import HARTREE_TO_KCAL
    from qm_nmr_calc.conformers.filters import filter_by_energy_window

    # Single conformer always passes
    if len(optimized_conformers) <= 1:
        return optimized_conformers

    # Extract conformer IDs and energies
    conf_ids = list(range(len(optimized_conformers)))
    energies_hartree = [c.energy for c in optimized_conformers]

    # Convert to kcal/mol for filtering
    energies_kcal = [e * HARTREE_TO_KCAL for e in energies_hartree]

    # Apply energy window filter
    kept_ids, _ = filter_by_energy_window(conf_ids, energies_kcal, window_kcal)

    # Return filtered conformers
    return [optimized_conformers[i] for i in kept_ids]
