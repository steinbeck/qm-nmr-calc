"""Tests for NMReData SDF generation module."""

import pytest

from qm_nmr_calc.nmredata import (
    NMREDATA_LEVEL,
    NMREDATA_SEP,
    NMREDATA_VERSION,
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
