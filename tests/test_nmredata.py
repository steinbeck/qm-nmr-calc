"""Tests for NMReData SDF generation module."""

import re

import pytest

from qm_nmr_calc.nmredata import (
    NMREDATA_LEVEL,
    NMREDATA_SEP,
    NMREDATA_VERSION,
    format_assignment_tag,
    format_atom_label,
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
