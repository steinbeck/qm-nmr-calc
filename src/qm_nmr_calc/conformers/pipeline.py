"""End-to-end conformer generation pipeline orchestration.

This module provides the high-level pipeline that connects conformer generation,
optimization, filtering, and file writing into a single function call.
"""

from rdkit.Chem import rdmolfiles

from ..models import ConformerData, ConformerEnsemble
from ..storage import create_conformer_directories, get_conformer_output_dir
from .clustering import cluster_and_select
from .filters import deduplicate_by_rmsd, filter_by_energy_window
from .generator import (
    calculate_num_conformers,
    generate_conformers_kdg,
    optimize_conformers_mmff,
)
from .xtb_ranking import detect_xtb_available, rank_conformers_by_xtb


def generate_conformer_ensemble(
    smiles: str,
    job_id: str,
    max_conformers: int | None = None,
    energy_window_kcal: float = 6.0,
    rmsd_threshold: float = 0.5,
    random_seed: int = 0xF00D,
    conformer_method: str = "rdkit_kdg",
    solvent: str = "chcl3",
    charge: int = 0,
    timeout_seconds: int = 7200,
    target_conformers_for_dft: int = 8,
    clustering_rmsd_threshold: float = 1.5,
) -> ConformerEnsemble:
    """Generate, optimize, filter, and write conformer ensemble from SMILES.

    This is the main pipeline function that orchestrates:
    1. Adaptive conformer count calculation (RDKit) OR CREST conformational search
    2. RDKit KDG conformer generation OR CREST ensemble generation
    3. MMFF optimization (RDKit only)
    4. RMSD-based deduplication (RDKit only - CREST handles internally)
    5. Energy window filtering
    6. Per-conformer directory creation
    7. XYZ file writing
    8. ConformerEnsemble model population

    Args:
        smiles: SMILES string of molecule
        job_id: Job identifier for directory creation
        max_conformers: Optional override for conformer count (None = adaptive, RDKit only)
        energy_window_kcal: Energy window in kcal/mol for filtering (default: 6.0)
        rmsd_threshold: RMSD threshold in Angstroms for deduplication (default: 0.5, RDKit only)
        random_seed: Random seed for reproducibility (default: 0xF00D, RDKit only)
        conformer_method: Method for conformer generation - "rdkit_kdg" or "crest" (default: "rdkit_kdg")
        solvent: Solvent for CREST ALPB solvation (default: "chcl3", CREST only)
        charge: Molecular charge (default: 0, CREST only)
        timeout_seconds: CREST timeout in seconds (default: 7200 = 2 hours, CREST only)
        target_conformers_for_dft: Target number of conformers to send to DFT (default: 8).
            After clustering, select up to this many diverse representatives.
        clustering_rmsd_threshold: RMSD threshold for diversity clustering (default: 1.5 Å).
            Larger than deduplication threshold to group similar conformations.

    Returns:
        ConformerEnsemble with populated conformers, energies, and file paths

    Raises:
        ValueError: If SMILES is invalid, molecule lacks MMFF parameters, CREST requested but not installed,
                    CREST requested for vacuum, or CREST requested for unsupported solvent
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
    # Dispatch to CREST or RDKit based on conformer_method
    if conformer_method == "crest":
        from .crest_generator import detect_crest_available, get_alpb_solvent, generate_crest_ensemble

        # Fail-fast: CREST requested but not installed
        if not detect_crest_available():
            raise ValueError(
                "CREST conformer method requested but CREST/xTB not installed. "
                "Install both crest and xtb binaries, or use 'rdkit_kdg' method."
            )

        # Fail-fast: CREST requires ALPB solvation
        alpb_solvent = get_alpb_solvent(solvent)
        if alpb_solvent is None:
            raise ValueError(
                f"CREST conformer method requires ALPB-compatible solvent. "
                f"Solvent '{solvent}' is not supported (vacuum or unsupported). "
                f"Use 'rdkit_kdg' method or change solvent to CHCl3 or DMSO."
            )

        return generate_crest_ensemble(
            smiles=smiles,
            job_id=job_id,
            solvent=alpb_solvent,
            charge=charge,
            energy_window_kcal=energy_window_kcal,
            timeout_seconds=timeout_seconds,
        )

    # RDKit path (default)
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

    # Step 8: Cluster conformers by RMSD for diversity selection
    # Only cluster if we have more conformers than the target
    if len(final_conf_ids) > target_conformers_for_dft:
        # Build energy lookup for ranking within clusters
        # Try xTB for better ranking correlation with DFT
        if detect_xtb_available():
            try:
                ranking_energies = rank_conformers_by_xtb(
                    mol, final_conf_ids,
                    charge=charge,
                    solvent=solvent,
                )
                ranking_method = "xtb_gfn2"
            except Exception:
                # Fall back to MMFF if xTB fails
                ranking_energies = dict(zip(final_conf_ids, final_energies))
                ranking_method = "mmff94s"
        else:
            ranking_energies = dict(zip(final_conf_ids, final_energies))
            ranking_method = "mmff94s"

        # Cluster and select representatives
        selected_conf_ids, num_clusters = cluster_and_select(
            mol, ranking_energies,
            rmsd_threshold=clustering_rmsd_threshold,
            max_conformers=target_conformers_for_dft,
        )

        # Filter to selected conformers
        selected_energies = [ranking_energies[cid] for cid in selected_conf_ids]

        # Track statistics
        pre_cluster_count = len(final_conf_ids)
        final_conf_ids = selected_conf_ids
        final_energies = selected_energies
    else:
        num_clusters = len(final_conf_ids)
        pre_cluster_count = len(final_conf_ids)
        ranking_method = "mmff94s"

    # Step 9: Create conformer string IDs (1-based, zero-padded)
    conformer_string_ids = [f"conf_{i+1:03d}" for i in range(len(final_conf_ids))]

    # Step 10: Create per-conformer storage directories
    create_conformer_directories(job_id, conformer_string_ids)

    # Step 11: Write XYZ files and build ConformerData list
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

    # Step 12: Build and return ConformerEnsemble
    ensemble = ConformerEnsemble(
        method="rdkit_kdg",
        conformers=conformer_data_list,
        pre_dft_energy_window_kcal=energy_window_kcal,
        total_generated=total_generated,
        total_after_pre_filter=pre_cluster_count,
    )

    return ensemble
