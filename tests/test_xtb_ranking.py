"""Tests for xTB-based conformer ranking.

These tests are designed to work whether xTB is installed or not.
Tests requiring xTB are skipped if the binary is not available.
"""

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from qm_nmr_calc.conformers.xtb_ranking import (
    detect_xtb_available,
    get_xtb_version,
    calculate_xtb_energy,
    rank_conformers_by_xtb,
    _map_solvent_to_alpb,
    HARTREE_TO_KCAL,
)


# Skip decorator for tests requiring xTB
requires_xtb = pytest.mark.skipif(
    not detect_xtb_available(),
    reason="xTB binary not available"
)


class TestXtbDetection:
    """Tests for xTB availability detection."""

    def test_detect_returns_bool(self):
        """Detection should return boolean."""
        result = detect_xtb_available()
        assert isinstance(result, bool)

    def test_version_returns_string_or_none(self):
        """Version should return string if available, None otherwise."""
        version = get_xtb_version()
        if detect_xtb_available():
            assert isinstance(version, str)
        else:
            assert version is None


class TestSolventMapping:
    """Tests for solvent name mapping."""

    def test_chloroform_variants(self):
        """All chloroform variants should map to chcl3."""
        assert _map_solvent_to_alpb("chcl3") == "chcl3"
        assert _map_solvent_to_alpb("chloroform") == "chcl3"
        assert _map_solvent_to_alpb("CDCl3") == "chcl3"

    def test_dmso_variants(self):
        """DMSO variants should map correctly."""
        assert _map_solvent_to_alpb("dmso") == "dmso"
        assert _map_solvent_to_alpb("DMSO-d6") == "dmso"

    def test_water_variants(self):
        """Water variants should map correctly."""
        assert _map_solvent_to_alpb("water") == "water"
        assert _map_solvent_to_alpb("H2O") == "water"
        assert _map_solvent_to_alpb("D2O") == "water"

    def test_unknown_solvent(self):
        """Unknown solvent should return None."""
        assert _map_solvent_to_alpb("unknown_solvent") is None

    def test_case_insensitive(self):
        """Mapping should be case-insensitive."""
        assert _map_solvent_to_alpb("CHCL3") == "chcl3"
        assert _map_solvent_to_alpb("Methanol") == "methanol"


@requires_xtb
class TestXtbEnergyCalculation:
    """Tests for xTB energy calculation (requires xTB installed)."""

    def test_simple_molecule(self):
        """Should calculate energy for simple molecule."""
        # Methane XYZ
        xyz = """5

C    0.000000    0.000000    0.000000
H    0.629118    0.629118    0.629118
H   -0.629118   -0.629118    0.629118
H   -0.629118    0.629118   -0.629118
H    0.629118   -0.629118   -0.629118
"""
        energy = calculate_xtb_energy(xyz)

        assert isinstance(energy, float)
        assert energy < 0  # Total energy should be negative

    def test_with_solvent(self):
        """Should run with solvation model."""
        xyz = """5

C    0.000000    0.000000    0.000000
H    0.629118    0.629118    0.629118
H   -0.629118   -0.629118    0.629118
H   -0.629118    0.629118   -0.629118
H    0.629118   -0.629118   -0.629118
"""
        energy_gas = calculate_xtb_energy(xyz, solvent=None)
        energy_solv = calculate_xtb_energy(xyz, solvent="chcl3")

        # Both should be valid energies
        assert isinstance(energy_gas, float)
        assert isinstance(energy_solv, float)
        # Solvation should change energy (slightly)
        assert energy_gas != energy_solv

    def test_charged_molecule(self):
        """Should handle charged molecules."""
        # Hydronium ion H3O+
        xyz = """4

O    0.000000    0.000000    0.000000
H    0.960000    0.000000    0.000000
H   -0.480000    0.831384    0.000000
H   -0.480000   -0.831384    0.000000
"""
        energy = calculate_xtb_energy(xyz, charge=1)
        assert isinstance(energy, float)


@requires_xtb
class TestRankConformersByXtb:
    """Tests for conformer ranking (requires xTB installed)."""

    @pytest.fixture
    def ethane_conformers(self):
        """Create ethane with conformers."""
        mol = Chem.MolFromSmiles("CC")
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        AllChem.EmbedMultipleConfs(mol, numConfs=5, params=params)
        return mol

    def test_returns_relative_energies(self, ethane_conformers):
        """Should return energies relative to minimum."""
        mol = ethane_conformers
        conf_ids = [conf.GetId() for conf in mol.GetConformers()]

        energies = rank_conformers_by_xtb(mol, conf_ids)

        # All energies should be non-negative (relative to min)
        assert all(e >= 0 for e in energies.values())
        # Minimum should be zero
        assert min(energies.values()) == pytest.approx(0.0, abs=1e-6)
        # Energies should be in kcal/mol (reasonable range)
        assert max(energies.values()) < 50  # No conformer should be 50 kcal/mol above min

    def test_returns_dict_for_all_conformers(self, ethane_conformers):
        """Should return energy for each conformer."""
        mol = ethane_conformers
        conf_ids = [conf.GetId() for conf in mol.GetConformers()]

        energies = rank_conformers_by_xtb(mol, conf_ids)

        assert len(energies) == len(conf_ids)
        assert all(cid in energies for cid in conf_ids)


class TestXtbNotAvailable:
    """Tests for behavior when xTB is not installed."""

    def test_rank_raises_when_unavailable(self):
        """Should raise clear error when xTB not available."""
        if detect_xtb_available():
            pytest.skip("xTB is available, cannot test unavailable behavior")

        mol = Chem.MolFromSmiles("C")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)

        with pytest.raises(RuntimeError, match="xTB not available"):
            rank_conformers_by_xtb(mol, [0])

    def test_calculate_raises_when_unavailable(self):
        """Should raise clear error when xTB not available."""
        if detect_xtb_available():
            pytest.skip("xTB is available, cannot test unavailable behavior")

        with pytest.raises(RuntimeError, match="xTB binary not found"):
            calculate_xtb_energy("1\n\nC 0 0 0")
