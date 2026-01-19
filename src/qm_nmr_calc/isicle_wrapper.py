"""Thin wrapper around ISiCLE for geometry optimization and NMR calculations."""
import subprocess
import sys
from pathlib import Path
from typing import NamedTuple

import isicle

from .presets import CalculationPreset


class Versions(NamedTuple):
    """Software versions for reproducibility."""
    isicle: str
    nwchem: str


def validate_nwchem() -> None:
    """
    Validate NWChem is available and callable.
    Call at startup - exits if validation fails.
    """
    try:
        result = subprocess.run(
            ['which', 'nwchem'],
            capture_output=True,
            text=True,
            timeout=5
        )
        if result.returncode != 0:
            sys.exit("FATAL: nwchem not found in PATH. Ensure NWChem is installed.")
    except subprocess.TimeoutExpired:
        sys.exit("FATAL: 'which nwchem' timed out")
    except Exception as e:
        sys.exit(f"FATAL: Cannot verify nwchem: {e}")


def get_nwchem_version() -> str:
    """
    Get NWChem version string.
    Note: NWChem doesn't have a --version flag, so we hardcode the known version.
    In production, could parse from package manager or NWChem output header.
    """
    # TODO: Parse dynamically if needed
    return "7.0.2"


def get_versions() -> Versions:
    """Get software versions for reproducibility."""
    return Versions(
        isicle=isicle.__version__,
        nwchem=get_nwchem_version()
    )


def run_geometry_optimization(
    smiles: str,
    job_dir: Path,
    processes: int = 4,
    functional: str = 'b3lyp',
    basis_set: str = '6-31G*',
) -> Path:
    """
    Run geometry optimization via ISiCLE/NWChem.

    Args:
        smiles: SMILES string of molecule
        job_dir: Job directory for outputs and scratch
        processes: Number of MPI processes for NWChem
        functional: DFT functional (default: b3lyp)
        basis_set: Basis set (default: 6-31G*)

    Returns:
        Path to optimized geometry file (XYZ format)

    Raises:
        Exception: If calculation fails (propagate to Huey for error handling)
    """
    # 1. Load molecule from SMILES
    geom = isicle.load(smiles)

    # 2. Initial 3D embedding and force-field optimization
    geom = geom.initial_optimize(embed=True, forcefield='UFF', ff_iter=200)

    # 3. Set up scratch directory inside job directory
    scratch_dir = job_dir / 'scratch'
    scratch_dir.mkdir(exist_ok=True)

    # 4. DFT geometry optimization via NWChem
    wrapper = isicle.qm.dft(
        geom,
        backend='NWChem',
        tasks=['optimize'],
        functional=functional,
        basis_set=basis_set,
        scratch_dir=str(scratch_dir),
        processes=processes,
    )

    # 5. Parse results
    result = wrapper.parse()

    # 6. Save optimized geometry
    output_dir = job_dir / 'output'
    output_dir.mkdir(exist_ok=True)
    output_file = output_dir / 'optimized.xyz'
    isicle.save(str(output_file), result['geometry'])

    return output_file


def run_nmr_calculation(
    smiles: str,
    job_dir: Path,
    preset: CalculationPreset,
    solvent: str,
    processes: int = 4,
) -> dict:
    """
    Run geometry optimization + NMR shielding calculation via ISiCLE/NWChem.

    This is a two-step DFT workflow:
    1. Geometry optimization with the preset's basis_set
    2. NMR shielding calculation with the preset's nmr_basis_set

    Args:
        smiles: SMILES string of molecule
        job_dir: Job directory for outputs and scratch
        preset: Calculation preset dict with functional, basis_set, nmr_basis_set
        solvent: NWChem COSMO solvent name
        processes: Number of MPI processes

    Returns:
        dict with 'geometry_file', 'shielding_data', 'energy'

    Raises:
        RuntimeError: If calculation fails or doesn't produce shielding data
    """
    # 1. Load molecule from SMILES
    geom = isicle.load(smiles)

    # 2. Initial 3D embedding and force-field optimization
    geom = geom.initial_optimize(embed=True, forcefield='UFF', ff_iter=200)

    # 3. Set up scratch directory inside job directory
    scratch_dir = job_dir / 'scratch'
    scratch_dir.mkdir(exist_ok=True)

    # 4. DFT geometry optimization
    opt_wrapper = isicle.qm.dft(
        geom,
        backend='NWChem',
        tasks=['optimize'],
        functional=preset['functional'],
        basis_set=preset['basis_set'],
        cosmo=True,
        solvent=solvent,
        max_iter=preset['max_iter'],
        scratch_dir=str(scratch_dir),
        processes=processes,
    )
    opt_result = opt_wrapper.parse()
    optimized_geom = opt_result['geometry']

    # 5. Save optimized geometry to output directory
    output_dir = job_dir / 'output'
    output_dir.mkdir(exist_ok=True)
    output_file = output_dir / 'optimized.xyz'
    isicle.save(str(output_file), optimized_geom)

    # 6. DFT NMR shielding calculation on optimized geometry
    nmr_wrapper = isicle.qm.dft(
        optimized_geom,
        backend='NWChem',
        tasks=['shielding'],
        functional=preset['functional'],
        basis_set=preset['nmr_basis_set'],
        cosmo=True,
        solvent=solvent,
        scratch_dir=str(scratch_dir),
        processes=processes,
    )
    nmr_result = nmr_wrapper.parse()

    # 7. Extract and validate shielding data
    shielding_data = nmr_result.get('shielding')
    if shielding_data is None:
        raise RuntimeError("NMR calculation did not produce shielding data")

    return {
        'geometry_file': str(output_file),
        'shielding_data': shielding_data,
        'energy': nmr_result.get('energy'),
    }
