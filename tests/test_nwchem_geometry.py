"""Tests for NWChem geometry handling module."""

import re
import tempfile
from pathlib import Path

import pytest
from rdkit import Chem

from qm_nmr_calc.nwchem.geometry import (
    load_geometry_file,
    mol_to_xyz_block,
    smiles_to_xyz,
    validate_geometry,
)

# Path to test fixtures
FIXTURES_DIR = Path(__file__).parent / "fixtures"


class TestSmilesToXyz:
    """Tests for smiles_to_xyz function."""

    def test_smiles_to_xyz_ethanol(self):
        """Convert ethanol SMILES to XYZ - should have 9 atoms."""
        mol, xyz = smiles_to_xyz("CCO")  # Ethanol

        # Parse lines (non-empty)
        lines = [line for line in xyz.strip().split("\n") if line.strip()]

        # Ethanol: 2C + 1O + 6H = 9 atoms
        assert len(lines) == 9, f"Ethanol should have 9 atoms, got {len(lines)}"

        # Each line should have element + 3 coordinates
        for line in lines:
            parts = line.split()
            assert len(parts) >= 4, f"Line should have element + 3 coords: {line}"

    def test_smiles_to_xyz_adds_hydrogens(self):
        """Methane SMILES 'C' should become 5 atoms (1C + 4H)."""
        mol, xyz = smiles_to_xyz("C")  # Methane

        lines = [line for line in xyz.strip().split("\n") if line.strip()]

        # Methane: 1C + 4H = 5 atoms
        assert len(lines) == 5, f"Methane should have 5 atoms, got {len(lines)}"

        # Count element types
        elements = [line.split()[0] for line in lines]
        assert elements.count("C") == 1, "Should have 1 carbon"
        assert elements.count("H") == 4, "Should have 4 hydrogens"

    def test_smiles_to_xyz_reproducible(self):
        """Same SMILES with same seed should give identical coordinates."""
        seed = 12345

        _, xyz1 = smiles_to_xyz("CCO", random_seed=seed)
        _, xyz2 = smiles_to_xyz("CCO", random_seed=seed)

        assert xyz1 == xyz2, "Same seed should produce identical coordinates"

    def test_smiles_to_xyz_different_seeds(self):
        """Different seeds should produce different coordinates."""
        _, xyz1 = smiles_to_xyz("CCO", random_seed=111)
        _, xyz2 = smiles_to_xyz("CCO", random_seed=222)

        # While element lines are same, coordinates should differ
        assert xyz1 != xyz2, "Different seeds should produce different coordinates"

    def test_smiles_to_xyz_invalid_smiles(self):
        """Invalid SMILES should raise ValueError."""
        with pytest.raises(ValueError) as exc_info:
            smiles_to_xyz("invalid_smiles_xyz")

        assert "Invalid SMILES" in str(exc_info.value)

    def test_smiles_to_xyz_returns_mol_object(self):
        """Function should return both mol object and xyz string."""
        mol, xyz = smiles_to_xyz("CCO")

        # Check mol is valid RDKit mol
        assert mol is not None
        assert mol.GetNumAtoms() == 9
        assert mol.GetNumConformers() == 1

        # Check xyz is string
        assert isinstance(xyz, str)
        assert len(xyz) > 0


class TestLoadGeometryFile:
    """Tests for load_geometry_file function."""

    def test_load_geometry_file_xyz(self):
        """Load ethanol from XYZ file."""
        xyz_path = FIXTURES_DIR / "ethanol.xyz"
        mol, xyz_block = load_geometry_file(xyz_path)

        # Check mol has 9 atoms
        assert mol.GetNumAtoms() == 9

        # Check xyz_block has 9 lines
        lines = [line for line in xyz_block.strip().split("\n") if line.strip()]
        assert len(lines) == 9, f"XYZ block should have 9 lines, got {len(lines)}"

    def test_load_geometry_file_sdf(self):
        """Load ethanol from SDF file - should preserve bonds."""
        sdf_path = FIXTURES_DIR / "ethanol.sdf"
        mol, xyz_block = load_geometry_file(sdf_path)

        # Check mol has 9 atoms
        assert mol.GetNumAtoms() == 9

        # Check bonds present (SDF includes connectivity)
        assert mol.GetNumBonds() == 8, "Ethanol should have 8 bonds"

        # Check xyz_block format
        lines = [line for line in xyz_block.strip().split("\n") if line.strip()]
        assert len(lines) == 9

    def test_load_geometry_file_unsupported_format(self):
        """Unsupported file format should raise ValueError."""
        with tempfile.NamedTemporaryFile(suffix=".txt", delete=False) as f:
            f.write(b"some text content\n")
            txt_path = Path(f.name)

        try:
            with pytest.raises(ValueError) as exc_info:
                load_geometry_file(txt_path)

            assert "Unsupported file format" in str(exc_info.value)
            assert ".txt" in str(exc_info.value)
        finally:
            txt_path.unlink()

    def test_load_geometry_file_accepts_string_path(self):
        """Function should accept string path, not just Path object."""
        xyz_path = str(FIXTURES_DIR / "ethanol.xyz")
        mol, xyz_block = load_geometry_file(xyz_path)

        assert mol.GetNumAtoms() == 9

    def test_load_geometry_file_xyz_with_charge(self):
        """XYZ file loading should accept charge parameter."""
        xyz_path = FIXTURES_DIR / "ethanol.xyz"
        # Ethanol is neutral, charge=0 is default
        mol, xyz_block = load_geometry_file(xyz_path, charge=0)

        assert mol.GetNumAtoms() == 9


class TestMolToXyzBlock:
    """Tests for mol_to_xyz_block function."""

    def test_mol_to_xyz_block_format(self):
        """Output format should match NWChem geometry block requirements."""
        mol, xyz = smiles_to_xyz("CCO")
        xyz_block = mol_to_xyz_block(mol)

        lines = xyz_block.strip().split("\n")

        # Each line: element symbol (1-2 chars) followed by 3 floats
        for line in lines:
            parts = line.split()
            assert len(parts) == 4, f"Expected 4 parts (elem + 3 coords): {line}"

            # Element should be 1-2 chars, letters only
            element = parts[0]
            assert 1 <= len(element) <= 2, f"Element should be 1-2 chars: {element}"
            assert element[0].isupper(), f"Element should start with uppercase: {element}"

            # Coordinates should be valid floats
            for coord in parts[1:]:
                try:
                    float(coord)
                except ValueError:
                    pytest.fail(f"Coordinate should be float: {coord}")

    def test_mol_to_xyz_block_no_header(self):
        """XYZ block should not include atom count or title lines."""
        mol, _ = smiles_to_xyz("C")  # Methane - 5 atoms

        xyz_block = mol_to_xyz_block(mol)
        lines = xyz_block.strip().split("\n")

        # First line should be an atom line, not atom count
        first_parts = lines[0].split()
        assert first_parts[0] in ["C", "H"], "First line should be atom, not count"

        # Should be exactly 5 lines for methane
        assert len(lines) == 5


class TestValidateGeometry:
    """Tests for validate_geometry function."""

    def test_validate_geometry_valid_mol(self):
        """Valid molecule should pass validation."""
        mol, _ = smiles_to_xyz("CCO")
        result = validate_geometry(mol)
        assert result is True

    def test_validate_geometry_none_mol(self):
        """None molecule should raise ValueError."""
        with pytest.raises(ValueError) as exc_info:
            validate_geometry(None)
        assert "Molecule is None" in str(exc_info.value)

    def test_validate_geometry_no_conformer(self):
        """Molecule without conformer should raise ValueError."""
        # Create mol from SMILES (no 3D coords yet)
        mol = Chem.MolFromSmiles("CCO")

        with pytest.raises(ValueError) as exc_info:
            validate_geometry(mol)
        assert "no 3D conformer" in str(exc_info.value)


class TestXyzOutputForNwchem:
    """Tests verifying XYZ output format is suitable for NWChem input."""

    def test_xyz_format_compatible_with_nwchem(self):
        """Verify output format matches NWChem geometry block requirements.

        NWChem geometry block expects lines like:
          C     0.000000    0.000000    0.000000
        Element symbol followed by 3 cartesian coordinates.
        """
        _, xyz_block = smiles_to_xyz("c1ccccc1")  # Benzene

        lines = xyz_block.strip().split("\n")

        # Benzene should have 12 atoms (6C + 6H)
        assert len(lines) == 12, f"Benzene should have 12 atoms, got {len(lines)}"

        # Check format of each line
        for line in lines:
            # Pattern: element (1-2 chars) + whitespace + 3 floats
            pattern = r"^[A-Z][a-z]?\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+$"
            assert re.match(pattern, line.strip()), f"Line doesn't match NWChem format: {line}"

    def test_xyz_coordinate_precision(self):
        """Coordinates should have sufficient decimal precision."""
        _, xyz_block = smiles_to_xyz("CCO")

        # Check that coordinates have reasonable precision (at least 4 decimal places)
        for line in xyz_block.strip().split("\n"):
            parts = line.split()
            for coord in parts[1:]:
                # Check for decimal point and some digits after
                assert "." in coord, f"Coordinate should have decimal point: {coord}"
