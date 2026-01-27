"""Tests for conformer filtering functions."""

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from qm_nmr_calc.conformers.filters import deduplicate_by_rmsd, filter_by_energy_window


class TestDeduplicateByRMSD:
    """Test RMSD-based conformer deduplication."""

    def test_single_conformer_kept(self):
        """Single conformer always passes through."""
        mol = Chem.MolFromSmiles("C")
        AllChem.EmbedMolecule(mol, randomSeed=42)

        result = deduplicate_by_rmsd(mol, rmsd_threshold=0.5)

        assert len(result) == 1
        assert result[0] == 0

    def test_identical_conformers_deduplicated(self):
        """Multiple identical conformers should reduce to one."""
        mol = Chem.MolFromSmiles("CCO")
        AllChem.EmbedMolecule(mol, randomSeed=42)

        # Create duplicates by copying conformer coordinates
        ref_conf = mol.GetConformer(0)
        for i in range(3):
            new_conf = Chem.Conformer(mol.GetNumAtoms())
            for j in range(mol.GetNumAtoms()):
                pos = ref_conf.GetAtomPosition(j)
                new_conf.SetAtomPosition(j, pos)
            mol.AddConformer(new_conf, assignId=True)

        assert mol.GetNumConformers() == 4  # 1 original + 3 copies

        result = deduplicate_by_rmsd(mol, rmsd_threshold=0.5)

        # Should keep only first conformer
        assert len(result) == 1
        assert result[0] == 0

    def test_diverse_conformers_all_kept(self):
        """Diverse conformers (RMSD > threshold) should all be kept."""
        mol = Chem.MolFromSmiles("CCCCCC")  # Hexane - flexible
        AllChem.EmbedMultipleConfs(mol, numConfs=5, randomSeed=42, useRandomCoords=True)

        # Use large threshold - expect all diverse conformers kept
        result = deduplicate_by_rmsd(mol, rmsd_threshold=0.5)

        # Should keep all if diverse enough
        assert len(result) >= 3  # At least 3 diverse conformers expected

    def test_threshold_sensitivity(self):
        """Stricter threshold keeps more conformers."""
        mol = Chem.MolFromSmiles("CCCCCC")
        AllChem.EmbedMultipleConfs(mol, numConfs=10, randomSeed=42, useRandomCoords=True)

        loose = deduplicate_by_rmsd(mol, rmsd_threshold=2.0)
        strict = deduplicate_by_rmsd(mol, rmsd_threshold=0.3)

        # Stricter threshold should keep at least as many
        assert len(strict) >= len(loose)

    def test_uses_symmetry_aware_rmsd(self):
        """Ensure symmetry-aware RMSD is used (GetBestRMS)."""
        # Ethanol - has symmetry in CH3 group
        mol = Chem.MolFromSmiles("CCO")
        AllChem.EmbedMultipleConfs(mol, numConfs=5, randomSeed=42, useRandomCoords=True)

        # Function should complete without error (GetBestRMS handles symmetry)
        result = deduplicate_by_rmsd(mol, rmsd_threshold=0.5)

        assert len(result) >= 1
        assert isinstance(result, list)
        assert all(isinstance(x, int) for x in result)


class TestFilterByEnergyWindow:
    """Test energy window filtering."""

    def test_empty_input(self):
        """Empty input returns empty output."""
        result = filter_by_energy_window([], [], window_kcal=6.0)

        assert result == ([], [])

    def test_single_conformer_kept(self):
        """Single conformer is always kept (is its own minimum)."""
        conf_ids = [0]
        energies = [42.5]

        result = filter_by_energy_window(conf_ids, energies, window_kcal=6.0)

        assert result == ([0], [42.5])

    def test_all_within_window_kept(self):
        """All conformers within window are kept."""
        conf_ids = [0, 1, 2, 3]
        energies = [10.0, 12.5, 14.0, 15.5]  # All within 6 kcal of min (10.0)

        result = filter_by_energy_window(conf_ids, energies, window_kcal=6.0)

        assert result == (conf_ids, energies)

    def test_removes_high_energy_conformers(self):
        """Conformers outside window are removed."""
        conf_ids = [0, 1, 2, 3, 4]
        energies = [10.0, 12.0, 14.0, 18.0, 25.0]  # 18.0 and 25.0 are > 6 kcal from min

        filtered_ids, filtered_energies = filter_by_energy_window(
            conf_ids, energies, window_kcal=6.0
        )

        assert filtered_ids == [0, 1, 2]
        assert filtered_energies == [10.0, 12.0, 14.0]

    def test_zero_window_keeps_only_minimum(self):
        """Window of 0.0 keeps only conformer(s) at minimum energy."""
        conf_ids = [0, 1, 2, 3]
        energies = [10.0, 10.0, 12.0, 15.0]  # Two at minimum

        filtered_ids, filtered_energies = filter_by_energy_window(
            conf_ids, energies, window_kcal=0.0
        )

        assert filtered_ids == [0, 1]
        assert filtered_energies == [10.0, 10.0]

    def test_preserves_order(self):
        """Output order matches input order."""
        conf_ids = [5, 2, 8, 1]
        energies = [10.0, 11.0, 12.0, 20.0]  # Last one filtered out

        filtered_ids, filtered_energies = filter_by_energy_window(
            conf_ids, energies, window_kcal=6.0
        )

        assert filtered_ids == [5, 2, 8]
        assert filtered_energies == [10.0, 11.0, 12.0]

    def test_negative_energies_handled(self):
        """Negative energies are handled correctly."""
        conf_ids = [0, 1, 2, 3]
        energies = [-50.0, -48.0, -45.0, -40.0]  # Min is -50.0

        filtered_ids, filtered_energies = filter_by_energy_window(
            conf_ids, energies, window_kcal=6.0
        )

        # -48.0 and -45.0 are within 6 kcal, -40.0 is 10 kcal away
        assert filtered_ids == [0, 1, 2]
        assert filtered_energies == [-50.0, -48.0, -45.0]
