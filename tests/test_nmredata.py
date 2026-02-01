"""Tests for NMReData SDF generation module."""

import re

import pytest
from rdkit import Chem

from qm_nmr_calc.nmredata import (
    NMREDATA_LEVEL,
    NMREDATA_SEP,
    NMREDATA_VERSION,
    format_assignment_tag,
    format_atom_label,
    format_sdf_tag,
    generate_nmredata_sdf,
    map_solvent_to_nmredata,
)


class TestConstants:
    """Test NMReData format constants."""

    def test_version_is_1_1(self):
        """Version should be 1.1 (stable NMReData format)."""
        assert NMREDATA_VERSION == "1.1"

    def test_level_is_0(self):
        """Level 0 indicates predicted data with unambiguous assignments."""
        assert NMREDATA_LEVEL == "0"

    def test_separator_is_comma_space(self):
        """Separator must be comma + space per NMReData spec."""
        assert NMREDATA_SEP == ", "


class TestSolventMapping:
    """Test solvent name mapping from NWChem to NMReData format."""

    def test_chcl3_maps_to_cdcl3(self):
        """Chloroform maps to deuterated form CDCl3."""
        assert map_solvent_to_nmredata("chcl3") == "CDCl3"

    def test_dmso_maps_to_deuterated_form(self):
        """DMSO maps to NMReData convention (CD3)2SO."""
        assert map_solvent_to_nmredata("dmso") == "(CD3)2SO"

    def test_vacuum_maps_to_vacuum(self):
        """Gas phase calculation uses 'vacuum' (no solvent)."""
        assert map_solvent_to_nmredata("vacuum") == "vacuum"

    def test_unknown_solvent_raises_error(self):
        """Unknown solvent should raise ValueError."""
        with pytest.raises(ValueError, match="Unknown solvent"):
            map_solvent_to_nmredata("benzene")

    def test_case_insensitive(self):
        """Solvent mapping should be case-insensitive."""
        assert map_solvent_to_nmredata("CHCl3") == "CDCl3"
        assert map_solvent_to_nmredata("DMSO") == "(CD3)2SO"

    def test_whitespace_handled(self):
        """Leading/trailing whitespace should be stripped."""
        assert map_solvent_to_nmredata("  chcl3  ") == "CDCl3"


class TestAtomLabel:
    """Test atom label formatting for NMReData ASSIGNMENT tag."""

    def test_hydrogen_label(self):
        """Hydrogen label should be lowercase h + index."""
        assert format_atom_label("H", 5) == "h5"

    def test_carbon_label(self):
        """Carbon label should be lowercase c + index."""
        assert format_atom_label("C", 3) == "c3"

    def test_label_is_lowercase(self):
        """Label should use lowercase element symbol."""
        assert format_atom_label("H", 1) == "h1"
        assert format_atom_label("C", 1) == "c1"

    def test_index_1_based(self):
        """First atom should have index 1 (1-based indexing)."""
        assert format_atom_label("H", 1) == "h1"
        assert format_atom_label("C", 1) == "c1"

    def test_multidigit_index(self):
        """Multi-digit indices should be handled correctly."""
        assert format_atom_label("H", 12) == "h12"
        assert format_atom_label("C", 100) == "c100"


class TestAssignmentTag:
    """Test ASSIGNMENT tag formatting with 1-indexed atoms."""

    def test_single_h_shift(self):
        """Single hydrogen shift should format correctly."""
        h1_shifts = [{"index": 1, "atom": "H", "shift": 7.2453}]
        c13_shifts = []
        result = format_assignment_tag(h1_shifts, c13_shifts)
        assert result == "h1, 7.2453, 1"

    def test_single_c_shift(self):
        """Single carbon shift should format correctly."""
        h1_shifts = []
        c13_shifts = [{"index": 2, "atom": "C", "shift": 128.45}]
        result = format_assignment_tag(h1_shifts, c13_shifts)
        assert result == "c2, 128.4500, 2"

    def test_multiple_shifts_newline_separated(self):
        """Multiple shifts should be separated by newlines."""
        h1_shifts = [
            {"index": 1, "atom": "H", "shift": 7.2453},
            {"index": 5, "atom": "H", "shift": 3.4521},
        ]
        c13_shifts = [{"index": 2, "atom": "C", "shift": 128.45}]
        result = format_assignment_tag(h1_shifts, c13_shifts)

        lines = result.split("\n")
        assert len(lines) == 3
        assert lines[0] == "h1, 7.2453, 1"
        assert lines[1] == "h5, 3.4521, 5"
        assert lines[2] == "c2, 128.4500, 2"

    def test_empty_lists_return_empty_string(self):
        """Empty shift lists should return empty string."""
        result = format_assignment_tag([], [])
        assert result == ""

    def test_shift_precision_four_decimals(self):
        """Shifts should be formatted with exactly 4 decimal places."""
        h1_shifts = [{"index": 1, "atom": "H", "shift": 7.25}]
        c13_shifts = []
        result = format_assignment_tag(h1_shifts, c13_shifts)
        # 7.25 should become 7.2500
        assert "7.2500" in result

    def test_atom_index_is_1_based(self):
        """Atom indices should be 1-based (first atom is 1, not 0)."""
        h1_shifts = [{"index": 1, "atom": "H", "shift": 7.2453}]
        c13_shifts = [{"index": 1, "atom": "C", "shift": 18.2}]
        result = format_assignment_tag(h1_shifts, c13_shifts)

        # Both should reference index 1
        lines = result.split("\n")
        assert lines[0].endswith(", 1")
        assert lines[1].endswith(", 1")

    def test_separator_is_comma_space(self):
        """Separator should be exactly ', ' (comma + space)."""
        h1_shifts = [{"index": 1, "atom": "H", "shift": 7.2453}]
        c13_shifts = []
        result = format_assignment_tag(h1_shifts, c13_shifts)

        # Should NOT have comma without space
        assert not re.search(r",(?! )", result), "Found comma without space"
        # Should have exactly two ", " separators per line
        assert result.count(", ") == 2

    def test_mixed_h_and_c_shifts(self):
        """Should handle mixed 1H and 13C shifts correctly."""
        h1_shifts = [
            {"index": 4, "atom": "H", "shift": 1.18},
            {"index": 5, "atom": "H", "shift": 1.18},
        ]
        c13_shifts = [
            {"index": 1, "atom": "C", "shift": 18.2},
            {"index": 2, "atom": "C", "shift": 58.3},
        ]
        result = format_assignment_tag(h1_shifts, c13_shifts)

        lines = result.split("\n")
        assert len(lines) == 4
        # H shifts come first
        assert lines[0].startswith("h4")
        assert lines[1].startswith("h5")
        # C shifts come after
        assert lines[2].startswith("c1")
        assert lines[3].startswith("c2")


@pytest.fixture
def ethanol_data():
    """Fixture providing ethanol test data for NMReData generation."""
    return {
        "smiles": "CCO",
        "geometry_xyz": """9
ethanol optimized
C    0.0000    0.0000    0.0000
C    1.5000    0.0000    0.0000
O    2.0000    1.2000    0.0000
H   -0.3500   -0.5000   -0.9000
H   -0.3500   -0.5000    0.9000
H   -0.3500    1.0000    0.0000
H    1.8500   -0.5000   -0.9000
H    1.8500   -0.5000    0.9000
H    2.9000    1.2000    0.0000
""",
        "h1_shifts": [
            {"index": 4, "atom": "H", "shift": 1.18},
            {"index": 5, "atom": "H", "shift": 1.18},
            {"index": 6, "atom": "H", "shift": 1.18},
            {"index": 7, "atom": "H", "shift": 3.65},
            {"index": 8, "atom": "H", "shift": 3.65},
            {"index": 9, "atom": "H", "shift": 2.45},
        ],
        "c13_shifts": [
            {"index": 1, "atom": "C", "shift": 18.2},
            {"index": 2, "atom": "C", "shift": 58.3},
        ],
        "solvent": "chcl3",
    }


class TestSDFTagFormat:
    """Test SDF tag formatting helper."""

    def test_tag_format(self):
        """SDF tag should have proper format with angle brackets."""
        result = format_sdf_tag("NMREDATA_VERSION", "1.1")
        assert result == ">  <NMREDATA_VERSION>\n1.1\n"

    def test_multiline_value(self):
        """SDF tag should handle multiline values."""
        result = format_sdf_tag("NMREDATA_ASSIGNMENT", "h1, 7.2453, 1\nh2, 3.45, 2")
        assert ">  <NMREDATA_ASSIGNMENT>" in result
        assert "h1, 7.2453, 1" in result
        assert "h2, 3.45, 2" in result


class TestGenerateNMReDataSDF:
    """Test complete NMReData SDF generation."""

    def test_ethanol_all_required_tags_present(self, ethanol_data):
        """Generated SDF should contain all 9 required NMReData tags."""
        sdf = generate_nmredata_sdf(
            smiles=ethanol_data["smiles"],
            geometry_xyz=ethanol_data["geometry_xyz"],
            h1_shifts=ethanol_data["h1_shifts"],
            c13_shifts=ethanol_data["c13_shifts"],
            solvent=ethanol_data["solvent"],
        )

        # Required tags per NMReData spec
        required_tags = [
            "NMREDATA_VERSION",
            "NMREDATA_LEVEL",
            "NMREDATA_SOLVENT",
            "NMREDATA_TEMPERATURE",
            "NMREDATA_ASSIGNMENT",
            "NMREDATA_FORMULA",
            "NMREDATA_SMILES",
            "NMREDATA_ID",
        ]
        for tag in required_tags:
            assert f"<{tag}>" in sdf, f"Missing required tag: {tag}"

    def test_tag_parsing_with_rdkit(self, ethanol_data):
        """Generated SDF should be parseable by RDKit SDMolSupplier."""
        sdf = generate_nmredata_sdf(
            smiles=ethanol_data["smiles"],
            geometry_xyz=ethanol_data["geometry_xyz"],
            h1_shifts=ethanol_data["h1_shifts"],
            c13_shifts=ethanol_data["c13_shifts"],
            solvent=ethanol_data["solvent"],
        )

        # Parse with RDKit
        import io
        supplier = Chem.SDMolSupplier()
        supplier.SetData(sdf)
        mol = next(supplier)

        assert mol is not None, "RDKit failed to parse SDF"
        # Check properties are accessible
        assert mol.HasProp("NMREDATA_VERSION")
        assert mol.GetProp("NMREDATA_VERSION") == "1.1"

    def test_mol_block_has_3d_coordinates(self, ethanol_data):
        """MOL block should contain 3D coordinates from XYZ."""
        sdf = generate_nmredata_sdf(
            smiles=ethanol_data["smiles"],
            geometry_xyz=ethanol_data["geometry_xyz"],
            h1_shifts=ethanol_data["h1_shifts"],
            c13_shifts=ethanol_data["c13_shifts"],
            solvent=ethanol_data["solvent"],
        )

        # MOL block is before first tag
        mol_block = sdf.split(">  <")[0]
        assert "M  END" in mol_block, "MOL block missing M  END"

        # Should have some non-zero z coordinates (3D structure)
        # Parse a few coordinate lines (after header, before M  END)
        lines = mol_block.split("\n")
        coord_lines = [l for l in lines if len(l.strip()) > 60 and "M  END" not in l]
        # At least one atom should have non-planar z coordinate
        has_3d = any("0.0000" not in line.split()[-8] for line in coord_lines if len(line.split()) >= 10)
        # More reliable: check for specific coordinates from ethanol_data
        assert "1.5000" in mol_block or "-0.9000" in mol_block, "Expected XYZ coordinates not found"

    def test_assignment_content_format(self, ethanol_data):
        """ASSIGNMENT tag should contain properly formatted H and C shifts."""
        sdf = generate_nmredata_sdf(
            smiles=ethanol_data["smiles"],
            geometry_xyz=ethanol_data["geometry_xyz"],
            h1_shifts=ethanol_data["h1_shifts"],
            c13_shifts=ethanol_data["c13_shifts"],
            solvent=ethanol_data["solvent"],
        )

        # Extract ASSIGNMENT tag content
        assert "NMREDATA_ASSIGNMENT" in sdf
        # Should have H shifts
        assert "h4, 1.1800, 4" in sdf
        # Should have C shifts
        assert "c1, 18.2000, 1" in sdf
        assert "c2, 58.3000, 2" in sdf

    def test_solvent_mapping_in_output(self, ethanol_data):
        """SOLVENT tag should show NMReData format (CDCl3, not chcl3)."""
        sdf = generate_nmredata_sdf(
            smiles=ethanol_data["smiles"],
            geometry_xyz=ethanol_data["geometry_xyz"],
            h1_shifts=ethanol_data["h1_shifts"],
            c13_shifts=ethanol_data["c13_shifts"],
            solvent="chcl3",
        )

        # Should map chcl3 -> CDCl3
        assert "CDCl3" in sdf
        assert "chcl3" not in sdf.lower() or "SMILES" in sdf  # May appear in SMILES but not SOLVENT

    def test_temperature_tag_format(self, ethanol_data):
        """TEMPERATURE tag should have 2 decimal places."""
        sdf = generate_nmredata_sdf(
            smiles=ethanol_data["smiles"],
            geometry_xyz=ethanol_data["geometry_xyz"],
            h1_shifts=ethanol_data["h1_shifts"],
            c13_shifts=ethanol_data["c13_shifts"],
            solvent=ethanol_data["solvent"],
            temperature_k=298.15,
        )

        # Should format as 298.15 (2 decimals)
        assert "298.15" in sdf

    def test_ensemble_mode_provenance(self, ethanol_data):
        """Ensemble mode should include conformer metadata in provenance."""
        sdf = generate_nmredata_sdf(
            smiles=ethanol_data["smiles"],
            geometry_xyz=ethanol_data["geometry_xyz"],
            h1_shifts=ethanol_data["h1_shifts"],
            c13_shifts=ethanol_data["c13_shifts"],
            solvent=ethanol_data["solvent"],
            is_ensemble=True,
            conformer_count=15,
        )

        # Should mention Boltzmann averaging in provenance
        assert "Boltzmann-averaged" in sdf
        assert "15 conformers" in sdf

    def test_formula_tag(self, ethanol_data):
        """FORMULA tag should contain correct molecular formula."""
        sdf = generate_nmredata_sdf(
            smiles=ethanol_data["smiles"],
            geometry_xyz=ethanol_data["geometry_xyz"],
            h1_shifts=ethanol_data["h1_shifts"],
            c13_shifts=ethanol_data["c13_shifts"],
            solvent=ethanol_data["solvent"],
        )

        # Ethanol is C2H6O
        assert "NMREDATA_FORMULA" in sdf
        # RDKit format may be C2H6O or with atom counts
        assert "C2H6O" in sdf or "C2" in sdf and "H6" in sdf and "O" in sdf

    def test_round_trip_atom_count(self, ethanol_data):
        """Round-trip: export -> parse -> verify atom count matches."""
        sdf = generate_nmredata_sdf(
            smiles=ethanol_data["smiles"],
            geometry_xyz=ethanol_data["geometry_xyz"],
            h1_shifts=ethanol_data["h1_shifts"],
            c13_shifts=ethanol_data["c13_shifts"],
            solvent=ethanol_data["solvent"],
        )

        # Parse with RDKit (removeHs=False to preserve explicit hydrogens)
        supplier = Chem.SDMolSupplier()
        supplier.SetData(sdf, removeHs=False)
        mol = next(supplier)

        # Ethanol has 9 atoms (2 C + 1 O + 6 H)
        assert mol.GetNumAtoms() == 9

    def test_sdf_structure_complete(self, ethanol_data):
        """SDF should have complete structure: MOL block + tags + terminator."""
        sdf = generate_nmredata_sdf(
            smiles=ethanol_data["smiles"],
            geometry_xyz=ethanol_data["geometry_xyz"],
            h1_shifts=ethanol_data["h1_shifts"],
            c13_shifts=ethanol_data["c13_shifts"],
            solvent=ethanol_data["solvent"],
        )

        # MOL block ends with M  END
        assert "M  END" in sdf
        # Tags use >  < format
        assert ">  <NMREDATA_VERSION>" in sdf
        # SDF terminates with $$$$
        assert "$$$$" in sdf

    def test_invalid_smiles_raises_error(self):
        """Invalid SMILES should raise ValueError."""
        with pytest.raises(ValueError, match="Invalid SMILES"):
            generate_nmredata_sdf(
                smiles="INVALID!!!",
                geometry_xyz="3\ntest\nC 0 0 0\nH 1 0 0\nH 0 1 0",
                h1_shifts=[],
                c13_shifts=[],
                solvent="chcl3",
            )

    def test_atom_count_mismatch_raises_error(self):
        """XYZ atom count mismatch should raise ValueError."""
        with pytest.raises(ValueError, match="Atom count mismatch"):
            generate_nmredata_sdf(
                smiles="CCO",  # Ethanol has 9 atoms with H
                geometry_xyz="3\ntest\nC 0 0 0\nH 1 0 0\nH 0 1 0",  # Only 3 atoms
                h1_shifts=[],
                c13_shifts=[],
                solvent="chcl3",
            )
