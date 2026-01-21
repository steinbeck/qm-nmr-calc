"""NWChem execution and calculation orchestration.

This module provides the main entry point for running NWChem calculations,
replacing the ISiCLE wrapper with direct subprocess execution and our own
input generation/output parsing.
"""

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

    Args:
        input_file: Path to NWChem input file (.nw)
        output_file: Path to write NWChem output
        processes: Number of MPI processes

    Raises:
        RuntimeError: If NWChem returns non-zero exit code
    """
    cmd = ["mpirun", "-np", str(processes), "nwchem", str(input_file)]

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        cwd=input_file.parent,
    )

    # Write stdout to output file
    output_file.write_text(result.stdout)

    if result.returncode != 0:
        error_msg = f"NWChem failed with exit code {result.returncode}"
        if result.stderr:
            error_msg += f"\nstderr: {result.stderr}"
        # Include last 50 lines of stdout for debugging
        stdout_lines = result.stdout.strip().split("\n")
        if stdout_lines:
            error_msg += f"\nLast output lines:\n" + "\n".join(stdout_lines[-50:])
        raise RuntimeError(error_msg)


def run_calculation(
    smiles: str,
    job_dir: Path,
    preset: dict,
    solvent: str,
    processes: int = 4,
    skip_optimization: bool = False,
    geometry_file: Path | None = None,
) -> dict:
    """Run complete NMR calculation (geometry optimization + NMR shielding).

    This is the main entry point replacing isicle_wrapper functions.
    Runs a two-step DFT calculation:
    1. Geometry optimization with COSMO solvation
    2. NMR shielding calculation with COSMO solvation

    Args:
        smiles: SMILES string of molecule
        job_dir: Job directory for outputs and scratch
        preset: Calculation preset dict with functional, basis_set, nmr_basis_set, max_iter
        solvent: Solvent name for COSMO (chcl3 or dmso)
        processes: Number of MPI processes
        skip_optimization: If True, skip geometry optimization and use geometry_file
        geometry_file: Path to pre-optimized geometry file (required if skip_optimization=True)

    Returns:
        Dict with:
            - 'geometry_file': Path to optimized.xyz
            - 'shielding_data': dict in shifts.py expected format
            - 'optimization_output': Path to optimize.out (or None if skipped)
            - 'shielding_output': Path to shielding.out

    Raises:
        ValueError: If skip_optimization=True but geometry_file not provided
        RuntimeError: If NWChem calculation fails
    """
    # Set up directories
    scratch_dir = job_dir / "scratch"
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
        # Step 1: Generate initial 3D geometry from SMILES
        mol, xyz_block = smiles_to_xyz(smiles)

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
