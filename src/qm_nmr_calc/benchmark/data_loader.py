"""Load DELTA50 benchmark data from repository files."""
from pathlib import Path
import json
from qm_nmr_calc.nwchem import load_geometry_file
from .models import MoleculeData, ExperimentalShifts


def get_data_dir() -> Path:
    """Get path to benchmark data directory."""
    # Relative to this file: src/qm_nmr_calc/benchmark -> data/benchmark
    return Path(__file__).parent.parent.parent.parent / "data" / "benchmark" / "delta50"


def load_experimental_shifts() -> ExperimentalShifts:
    """Load experimental shifts from JSON file.

    Returns:
        ExperimentalShifts object containing all experimental data

    Raises:
        FileNotFoundError: If experimental_shifts.json doesn't exist
        ValueError: If JSON structure is invalid
    """
    shifts_file = get_data_dir() / "experimental_shifts.json"

    if not shifts_file.exists():
        raise FileNotFoundError(
            f"Experimental shifts file not found: {shifts_file}\n"
            "See data/benchmark/delta50/README.md for setup instructions"
        )

    with shifts_file.open() as f:
        data = json.load(f)

    molecules = {}
    for mol_id, mol_data in data["molecules"].items():
        # Skip template and example entries
        if mol_id.startswith("_") or mol_id == "example":
            continue

        molecules[mol_id] = MoleculeData(
            id=mol_id,
            name=mol_data["name"],
            xyz_file=f"molecules/{mol_id}.xyz",
            h1_shifts=mol_data.get("h1_shifts", []),
            c13_shifts=mol_data.get("c13_shifts", []),
            h1_assignments=mol_data.get("h1_assignments"),
            c13_assignments=mol_data.get("c13_assignments"),
        )

    return ExperimentalShifts(source=data["source"], molecules=molecules)


def load_delta50_molecules() -> dict[str, tuple]:
    """Load all DELTA50 molecules from XYZ files.

    Returns:
        Dict mapping molecule_id to (RDKit Mol, XYZ block) tuple

    Raises:
        FileNotFoundError: If molecules directory doesn't exist
        RuntimeError: If no valid XYZ files are found
    """
    xyz_dir = get_data_dir() / "molecules"

    if not xyz_dir.exists():
        raise FileNotFoundError(
            f"Molecules directory not found: {xyz_dir}\n"
            "See data/benchmark/delta50/README.md for setup instructions"
        )

    molecules = {}
    xyz_files = sorted(xyz_dir.glob("compound_*.xyz"))

    if not xyz_files:
        raise RuntimeError(
            f"No compound_*.xyz files found in {xyz_dir}\n"
            "Expected: compound_01.xyz through compound_50.xyz\n"
            "See data/benchmark/delta50/README.md for setup instructions"
        )

    for xyz_file in xyz_files:
        molecule_id = xyz_file.stem
        mol, xyz_block = load_geometry_file(xyz_file, charge=0)
        molecules[molecule_id] = (mol, xyz_block)

    return molecules
