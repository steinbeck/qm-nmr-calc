"""Benchmark viewer router for DELTA50 data visualization."""

from pathlib import Path

from fastapi import APIRouter, HTTPException, Request
from fastapi.templating import Jinja2Templates

from ...benchmark.data_loader import load_delta50_molecules, load_experimental_shifts

# Template engine setup
templates = Jinja2Templates(
    directory=Path(__file__).resolve().parent.parent / "templates"
)

router = APIRouter(tags=["benchmark"])


@router.get("/benchmark/delta50")
async def benchmark_viewer(request: Request):
    """Render the DELTA50 benchmark data viewer page."""
    # Load experimental shifts to get molecule list
    shifts_data = load_experimental_shifts()

    # Build molecule list sorted by compound number
    molecule_list = []
    for mol_id, mol_data in shifts_data.molecules.items():
        # Extract compound number from id (e.g., "compound_01" -> 1)
        try:
            compound_number = int(mol_id.split("_")[1])
        except (IndexError, ValueError):
            compound_number = 0

        molecule_list.append({
            "id": mol_id,
            "name": mol_data.name,
            "compound_number": compound_number,
        })

    # Sort by compound number (1-50)
    molecule_list.sort(key=lambda x: x["compound_number"])

    return templates.TemplateResponse(
        request=request,
        name="benchmark_viewer.html",
        context={"molecule_list": molecule_list},
    )


@router.get("/benchmark/delta50/molecules/{compound_id}")
async def get_molecule_data(compound_id: str):
    """Get molecule structure and experimental shifts for 3D viewer.

    Args:
        compound_id: Molecule identifier (e.g., "compound_01")

    Returns:
        JSON with XYZ coordinates, shifts, and metadata
    """
    # Load experimental shifts
    shifts_data = load_experimental_shifts()

    # Validate compound_id exists
    if compound_id not in shifts_data.molecules:
        raise HTTPException(
            status_code=404,
            detail=f"Molecule '{compound_id}' not found in DELTA50 dataset",
        )

    mol_data = shifts_data.molecules[compound_id]

    # Load XYZ structure
    molecules = load_delta50_molecules()
    if compound_id not in molecules:
        raise HTTPException(
            status_code=404,
            detail=f"XYZ file for '{compound_id}' not found",
        )

    rdkit_mol, xyz_block = molecules[compound_id]

    # Extract compound number
    try:
        compound_number = int(compound_id.split("_")[1])
    except (IndexError, ValueError):
        compound_number = 0

    # Count atoms from XYZ block (H and C atoms)
    # Note: xyz_block from load_geometry_file has no header lines, just coordinates
    lines = xyz_block.strip().split("\n")
    num_h_atoms = 0
    num_c_atoms = 0
    total_atoms = 0
    for line in lines:
        parts = line.split()
        if parts:
            total_atoms += 1
            element = parts[0]
            if element == "H":
                num_h_atoms += 1
            elif element == "C":
                num_c_atoms += 1

    # Reconstruct full XYZ format with header for 3Dmol.js
    full_xyz = f"{total_atoms}\n{mol_data.name}\n{xyz_block}"

    # Build assignments dict with string keys (for JSON compatibility)
    h1_assignments = {}
    if mol_data.h1_assignments:
        h1_assignments = {str(k): v for k, v in mol_data.h1_assignments.items()}

    c13_assignments = {}
    if mol_data.c13_assignments:
        c13_assignments = {str(k): v for k, v in mol_data.c13_assignments.items()}

    return {
        "compound_id": compound_id,
        "name": mol_data.name,
        "compound_number": compound_number,
        "xyz": full_xyz,
        "h1_shifts": mol_data.h1_shifts,
        "c13_shifts": mol_data.c13_shifts,
        "h1_assignments": h1_assignments,
        "c13_assignments": c13_assignments,
        "num_h_atoms": num_h_atoms,
        "num_c_atoms": num_c_atoms,
    }
