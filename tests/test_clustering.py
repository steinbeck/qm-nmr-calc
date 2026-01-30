"""Tests for conformer clustering functionality."""

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from qm_nmr_calc.conformers.clustering import (
    cluster_conformers_by_rmsd,
    select_cluster_representatives,
    cluster_and_select,
)


@pytest.fixture
def flexible_mol_with_conformers():
    """Create a flexible molecule (heptane) with multiple conformers."""
    mol = Chem.MolFromSmiles("CCCCCCC")
    mol = Chem.AddHs(mol)

    # Generate many conformers
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    AllChem.EmbedMultipleConfs(mol, numConfs=30, params=params)

    # Optimize with MMFF
    AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=200)

    return mol


@pytest.fixture
def rigid_mol_with_conformers():
    """Create a rigid molecule (benzene) with conformers."""
    mol = Chem.MolFromSmiles("c1ccccc1")
    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    AllChem.EmbedMultipleConfs(mol, numConfs=10, params=params)
    AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=200)

    return mol


class TestClusterConformersByRMSD:
    """Tests for RMSD clustering function."""

    def test_reduces_conformer_count(self, flexible_mol_with_conformers):
        """Clustering should reduce conformer count."""
        mol = flexible_mol_with_conformers
        n_initial = mol.GetNumConformers()

        clusters = cluster_conformers_by_rmsd(mol, rmsd_threshold=1.5)

        assert len(clusters) < n_initial
        assert len(clusters) >= 1

    def test_all_conformers_assigned(self, flexible_mol_with_conformers):
        """Every conformer should be in exactly one cluster."""
        mol = flexible_mol_with_conformers
        conf_ids = {conf.GetId() for conf in mol.GetConformers()}

        clusters = cluster_conformers_by_rmsd(mol, rmsd_threshold=1.5)

        clustered_ids = set()
        for cluster in clusters:
            for cid in cluster:
                assert cid not in clustered_ids, "Conformer in multiple clusters"
                clustered_ids.add(cid)

        assert clustered_ids == conf_ids

    def test_single_conformer(self):
        """Single conformer should return single cluster."""
        mol = Chem.MolFromSmiles("C")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)

        clusters = cluster_conformers_by_rmsd(mol, rmsd_threshold=1.5)

        assert len(clusters) == 1
        assert len(clusters[0]) == 1

    def test_threshold_affects_cluster_count(self, flexible_mol_with_conformers):
        """Tighter threshold should produce more clusters."""
        mol = flexible_mol_with_conformers

        clusters_tight = cluster_conformers_by_rmsd(mol, rmsd_threshold=0.5)
        clusters_loose = cluster_conformers_by_rmsd(mol, rmsd_threshold=2.0)

        assert len(clusters_tight) >= len(clusters_loose)

    def test_rigid_molecule_few_clusters(self, rigid_mol_with_conformers):
        """Rigid molecule should have few unique conformers."""
        mol = rigid_mol_with_conformers

        clusters = cluster_conformers_by_rmsd(mol, rmsd_threshold=0.5)

        # Benzene conformers should all be similar
        assert len(clusters) <= 3


class TestSelectClusterRepresentatives:
    """Tests for representative selection function."""

    def test_selects_lowest_energy(self):
        """Should select lowest energy from each cluster."""
        clusters = [[0, 1, 2], [3, 4], [5]]
        energies = {0: 10.0, 1: 8.0, 2: 12.0, 3: 5.0, 4: 6.0, 5: 15.0}

        reps = select_cluster_representatives(clusters, energies, max_representatives=3)

        assert reps == [1, 3, 5]  # Lowest from each cluster

    def test_respects_max_representatives(self):
        """Should not exceed max_representatives."""
        clusters = [[0], [1], [2], [3], [4]]
        energies = {i: float(i) for i in range(5)}

        reps = select_cluster_representatives(clusters, energies, max_representatives=3)

        assert len(reps) == 3

    def test_handles_missing_energies(self):
        """Should handle conformers without energy gracefully."""
        clusters = [[0, 1], [2]]
        energies = {0: 10.0}  # Missing energies for 1 and 2

        reps = select_cluster_representatives(clusters, energies, max_representatives=2)

        assert 0 in reps  # Should select known-energy conformer

    def test_empty_clusters(self):
        """Should handle empty cluster list."""
        reps = select_cluster_representatives([], {}, max_representatives=8)
        assert reps == []


class TestClusterAndSelect:
    """Tests for convenience function."""

    def test_full_pipeline(self, flexible_mol_with_conformers):
        """Test complete cluster and select workflow."""
        mol = flexible_mol_with_conformers

        # Create mock energies
        energies = {
            conf.GetId(): float(i)
            for i, conf in enumerate(mol.GetConformers())
        }

        selected, n_clusters = cluster_and_select(
            mol, energies,
            rmsd_threshold=1.5,
            max_conformers=8
        )

        assert len(selected) <= 8
        assert len(selected) <= n_clusters
        assert n_clusters >= 1

    def test_returns_cluster_count(self, flexible_mol_with_conformers):
        """Should return number of clusters found."""
        mol = flexible_mol_with_conformers
        energies = {conf.GetId(): 0.0 for conf in mol.GetConformers()}

        _, n_clusters = cluster_and_select(mol, energies)

        assert n_clusters > 0
        assert n_clusters <= mol.GetNumConformers()
