"""RMSD-based conformer clustering for pre-DFT selection.

Uses Butina algorithm to group structurally similar conformers,
then selects diverse representatives for expensive DFT calculations.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from rdkit.ML.Cluster import Butina


def cluster_conformers_by_rmsd(
    mol: Chem.Mol,
    rmsd_threshold: float = 1.5,
) -> list[list[int]]:
    """Cluster conformers using Butina algorithm with symmetry-aware RMSD.

    Groups structurally similar conformers into clusters. Each cluster
    represents a unique conformational basin.

    Args:
        mol: RDKit Mol with multiple conformers
        rmsd_threshold: RMSD threshold in Angstroms for clustering.
            Typical values: 0.5 (strict), 1.0 (balanced), 1.5 (aggressive)

    Returns:
        List of clusters, each cluster is a list of conformer IDs.
        Clusters are sorted by size (largest first).
        First conformer in each cluster is the centroid.
    """
    conf_ids = [conf.GetId() for conf in mol.GetConformers()]
    n_confs = len(conf_ids)

    if n_confs <= 1:
        return [[conf_ids[0]]] if conf_ids else []

    # Build distance matrix (lower triangle for Butina)
    # Uses symmetry-aware RMSD via GetBestRMS
    dists = []
    for i in range(1, n_confs):
        for j in range(i):
            rmsd = rdMolAlign.GetBestRMS(
                mol, mol,
                prbId=conf_ids[i],
                refId=conf_ids[j]
            )
            dists.append(rmsd)

    # Butina clustering returns tuple of tuples
    clusters = Butina.ClusterData(
        dists,
        n_confs,
        distThresh=rmsd_threshold,
        isDistData=True
    )

    # Convert indices back to conformer IDs
    return [[conf_ids[idx] for idx in cluster] for cluster in clusters]


def select_cluster_representatives(
    clusters: list[list[int]],
    energies: dict[int, float],
    max_representatives: int = 8,
) -> list[int]:
    """Select lowest-energy conformer from each cluster.

    Args:
        clusters: List of clusters from cluster_conformers_by_rmsd
        energies: Dict mapping conformer ID to energy (any unit, lower is better)
        max_representatives: Maximum number of representatives to return

    Returns:
        List of conformer IDs (one per cluster, lowest energy).
        Returns up to max_representatives conformers.
    """
    representatives = []

    for cluster in clusters:
        if len(representatives) >= max_representatives:
            break

        # Find lowest energy conformer in this cluster
        cluster_with_energy = [
            (cid, energies.get(cid, float('inf')))
            for cid in cluster
        ]
        cluster_with_energy.sort(key=lambda x: x[1])
        representatives.append(cluster_with_energy[0][0])

    return representatives


def cluster_and_select(
    mol: Chem.Mol,
    energies: dict[int, float],
    rmsd_threshold: float = 1.5,
    max_conformers: int = 8,
) -> tuple[list[int], int]:
    """Convenience function: cluster conformers and select representatives.

    Args:
        mol: RDKit Mol with multiple conformers
        energies: Dict mapping conformer ID to energy
        rmsd_threshold: RMSD threshold for clustering
        max_conformers: Target number of conformers to select

    Returns:
        Tuple of (selected_conf_ids, num_clusters)
    """
    clusters = cluster_conformers_by_rmsd(mol, rmsd_threshold)
    selected = select_cluster_representatives(clusters, energies, max_conformers)
    return selected, len(clusters)
