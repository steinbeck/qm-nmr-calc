"""End-to-end conformer generation pipeline orchestration.

This module provides the high-level pipeline that connects conformer generation,
optimization, filtering, and file writing into a single function call.
"""

from rdkit.Chem import rdmolfiles

from ..models import ConformerData, ConformerEnsemble
from ..storage import create_conformer_directories, get_conformer_output_dir
from .filters import deduplicate_by_rmsd, filter_by_energy_window
from .generator import (
    calculate_num_conformers,
    generate_conformers_kdg,
    optimize_conformers_mmff,
)


def generate_conformer_ensemble(
    smiles: str,
    job_id: str,
    max_conformers: int | None = None,
    energy_window_kcal: float = 6.0,
    rmsd_threshold: float = 0.5,
    random_seed: int = 0xF00D,
) -> ConformerEnsemble:
    """Generate, optimize, filter, and write conformer ensemble from SMILES.

    This is the main pipeline function that orchestrates:
    1. Adaptive conformer count calculation
    2. RDKit KDG conformer generation
    3. MMFF optimization
    4. RMSD-based deduplication
    5. Energy window filtering
    6. Per-conformer directory creation
    7. XYZ file writing
    8. ConformerEnsemble model population

    Args:
        smiles: SMILES string of molecule
        job_id: Job identifier for directory creation
        max_conformers: Optional override for conformer count (None = adaptive)
        energy_window_kcal: Energy window in kcal/mol for filtering (default: 6.0)
        rmsd_threshold: RMSD threshold in Angstroms for deduplication (default: 0.5)
        random_seed: Random seed for reproducibility (default: 0xF00D)

    Returns:
        ConformerEnsemble with populated conformers, energies, and file paths

    Raises:
        ValueError: If SMILES is invalid or molecule lacks MMFF parameters
        RuntimeError: If conformer generation fails

    Example:
        >>> ensemble = generate_conformer_ensemble("CCO", "job123")
        >>> len(ensemble.conformers)
        5
        >>> ensemble.conformers[0].energy
        12.34
        >>> ensemble.conformers[0].geometry_file
        'output/conformers/conf_001/geometry.xyz'
    """
    # Step 1: Calculate adaptive conformer count
    num_confs = calculate_num_conformers(smiles, max_conformers)

    # Step 2: Generate conformers using KDG
    mol = generate_conformers_kdg(smiles, num_confs, random_seed)
    total_generated = mol.GetNumConformers()

    # Step 3: Optimize with MMFF
    optimization_results = optimize_conformers_mmff(mol)

    # Step 4: Extract conformer IDs and energies
    conf_ids = [conf.GetId() for conf in mol.GetConformers()]
    energies = [result[1] for result in optimization_results]  # (not_converged, energy)

    # Step 5: Deduplicate by RMSD (run dedup FIRST)
    kept_after_dedup = deduplicate_by_rmsd(mol, rmsd_threshold)

    # Step 6: Filter to only deduped conformers
    # Create parallel lists of only the deduped conformers
    deduped_conf_ids = []
    deduped_energies = []
    for conf_id, energy in zip(conf_ids, energies):
        if conf_id in kept_after_dedup:
            deduped_conf_ids.append(conf_id)
            deduped_energies.append(energy)

    # Step 7: Apply energy window filter on deduped subset
    final_conf_ids, final_energies = filter_by_energy_window(
        deduped_conf_ids, deduped_energies, energy_window_kcal
    )

    # Step 8: Create conformer string IDs (1-based, zero-padded)
    conformer_string_ids = [f"conf_{i+1:03d}" for i in range(len(final_conf_ids))]

    # Step 9: Create per-conformer storage directories
    create_conformer_directories(job_id, conformer_string_ids)

    # Step 10: Write XYZ files and build ConformerData list
    conformer_data_list = []
    for string_id, conf_id, energy in zip(
        conformer_string_ids, final_conf_ids, final_energies
    ):
        # Get output directory for this conformer
        output_dir = get_conformer_output_dir(job_id, string_id)
        xyz_path = output_dir / "geometry.xyz"

        # Write XYZ file using RDKit
        rdmolfiles.MolToXYZFile(mol, str(xyz_path), confId=conf_id)

        # Build ConformerData with relative path
        geometry_file_relative = f"output/conformers/{string_id}/geometry.xyz"
        conformer_data = ConformerData(
            conformer_id=string_id,
            energy=energy,
            energy_unit="kcal_mol",
            geometry_file=geometry_file_relative,
            status="pending",
        )
        conformer_data_list.append(conformer_data)

    # Step 11: Build and return ConformerEnsemble
    ensemble = ConformerEnsemble(
        method="rdkit_kdg",
        conformers=conformer_data_list,
        pre_dft_energy_window_kcal=energy_window_kcal,
        total_generated=total_generated,
        total_after_pre_filter=len(final_conf_ids),
    )

    return ensemble
