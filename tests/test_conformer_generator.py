"""
Unit tests for conformer generation module.

Tests KDG conformer generation, MMFF optimization, and adaptive conformer counting.
"""

import pytest
from rdkit import Chem

from qm_nmr_calc.conformers.generator import (
    calculate_num_conformers,
    generate_conformers_kdg,
    optimize_conformers_mmff,
)


class TestCalculateNumConformers:
    """Test adaptive conformer count calculation."""

    def test_rigid_molecule_returns_50(self):
        """Ethanol (0 rotatable bonds) should return 50."""
        result = calculate_num_conformers("CCO")
        assert result == 50

    def test_flexible_molecule_returns_200(self):
        """Dodecane (~9 rotatable bonds) should return 200."""
        result = calculate_num_conformers("CCCCCCCCCCCC")
        assert result == 200

    def test_max_conformers_override(self):
        """max_conformers parameter should override adaptive count."""
        result = calculate_num_conformers("CCO", max_conformers=10)
        assert result == 10

        result = calculate_num_conformers("CCCCCCCCCCCC", max_conformers=5)
        assert result == 5

    def test_invalid_smiles_raises_error(self):
        """Invalid SMILES should raise ValueError."""
        with pytest.raises(ValueError):
            calculate_num_conformers("INVALID_SMILES_123")


class TestGenerateConformersKDG:
    """Test KDG conformer generation."""

    def test_simple_molecule_generates_conformers(self):
        """Ethanol should generate requested conformers."""
        mol = generate_conformers_kdg("CCO", num_confs=5, random_seed=0xF00D)

        assert mol is not None
        assert mol.GetNumConformers() >= 1
        # May be fewer than 5 for small molecule, but should have at least 1
        assert mol.GetNumConformers() <= 5

    def test_generated_mol_has_hydrogens(self):
        """Generated molecule should include hydrogens."""
        mol = generate_conformers_kdg("CCO", num_confs=3)

        # Ethanol C2H6O should have 9 atoms total
        assert mol.GetNumAtoms() == 9

    def test_conformers_have_3d_coordinates(self):
        """Generated conformers should have 3D coordinates."""
        mol = generate_conformers_kdg("CCO", num_confs=3)

        conf = mol.GetConformer(0)
        pos = conf.GetAtomPosition(0)

        # Should have x, y, z coordinates
        assert hasattr(pos, 'x')
        assert hasattr(pos, 'y')
        assert hasattr(pos, 'z')

    def test_invalid_smiles_raises_error(self):
        """Invalid SMILES should raise ValueError."""
        with pytest.raises(ValueError):
            generate_conformers_kdg("INVALID_SMILES", num_confs=5)

    def test_reproducible_with_seed(self):
        """Same seed should produce same conformers."""
        mol1 = generate_conformers_kdg("CCCC", num_confs=5, random_seed=42)
        mol2 = generate_conformers_kdg("CCCC", num_confs=5, random_seed=42)

        assert mol1.GetNumConformers() == mol2.GetNumConformers()

        # Check first conformer coordinates match
        conf1 = mol1.GetConformer(0)
        conf2 = mol2.GetConformer(0)

        pos1 = conf1.GetAtomPosition(0)
        pos2 = conf2.GetAtomPosition(0)

        assert abs(pos1.x - pos2.x) < 1e-6
        assert abs(pos1.y - pos2.y) < 1e-6
        assert abs(pos1.z - pos2.z) < 1e-6


class TestOptimizeConformersMMFF:
    """Test MMFF conformer optimization."""

    def test_optimization_returns_results_list(self):
        """Optimization should return list of (not_converged, energy) tuples."""
        mol = generate_conformers_kdg("CCO", num_confs=3)
        results = optimize_conformers_mmff(mol, max_iters=200)

        assert isinstance(results, list)
        assert len(results) == mol.GetNumConformers()

        for not_converged, energy in results:
            assert isinstance(not_converged, int)
            assert isinstance(energy, float)

    def test_energies_are_reasonable(self):
        """MMFF energies should be in reasonable range for small molecules."""
        mol = generate_conformers_kdg("CCO", num_confs=3)
        results = optimize_conformers_mmff(mol)

        for not_converged, energy in results:
            # Ethanol MMFF energy should be in kcal/mol, roughly -10 to 50
            assert -100 < energy < 200

    def test_molecule_without_mmff_params_raises_error(self):
        """Molecule without MMFF params should raise ValueError."""
        # Create a molecule with unsupported atom types (e.g., boron with unusual valence)
        # Simple boron compounds often lack MMFF parameters
        mol = Chem.MolFromSmiles("[B]")
        mol = Chem.AddHs(mol)
        Chem.AllChem.EmbedMolecule(mol, randomSeed=42)

        with pytest.raises(ValueError, match="MMFF"):
            optimize_conformers_mmff(mol)

    def test_max_iters_parameter(self):
        """max_iters parameter should be respected."""
        mol = generate_conformers_kdg("CCO", num_confs=2)

        # Should run without error
        results = optimize_conformers_mmff(mol, max_iters=100)
        assert len(results) == 2


class TestIntegration:
    """Integration tests for full workflow."""

    def test_full_workflow_rigid_molecule(self):
        """Full workflow for rigid molecule."""
        smiles = "CCO"

        # Calculate conformer count
        num_confs = calculate_num_conformers(smiles)
        assert num_confs == 50

        # Generate conformers
        mol = generate_conformers_kdg(smiles, num_confs=num_confs)
        assert mol.GetNumConformers() >= 1

        # Optimize conformers
        results = optimize_conformers_mmff(mol)
        assert len(results) == mol.GetNumConformers()

    def test_full_workflow_flexible_molecule(self):
        """Full workflow for flexible molecule."""
        smiles = "CCCCCCCCCCCC"  # Dodecane, flexible (9 rotatable bonds)

        # Calculate conformer count
        num_confs = calculate_num_conformers(smiles)
        assert num_confs == 200

        # Generate conformers (use smaller count for test speed)
        mol = generate_conformers_kdg(smiles, num_confs=10)
        assert mol.GetNumConformers() >= 1

        # Optimize conformers
        results = optimize_conformers_mmff(mol)
        assert len(results) == mol.GetNumConformers()

        # Energies should vary for different conformers
        energies = [energy for _, energy in results]
        assert max(energies) - min(energies) > 0.1  # At least 0.1 kcal/mol variation
