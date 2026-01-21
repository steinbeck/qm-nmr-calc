"""Tests for NWChem output parsing module."""

from pathlib import Path

import pytest

from qm_nmr_calc.nwchem.output_parser import (
    extract_optimized_geometry,
    parse_shielding_output,
)
from qm_nmr_calc.shifts import shielding_to_shift


FIXTURES_DIR = Path(__file__).parent / "fixtures"


class TestExtractOptimizedGeometry:
    """Tests for extract_optimized_geometry function."""

    def test_parses_coordinates_from_fixture(self):
        """Extract geometry from mock NWChem optimization output."""
        output_text = (FIXTURES_DIR / "nwchem_optimization_output.txt").read_text()
        xyz = extract_optimized_geometry(output_text)

        # Should have 5 atom lines (methane: 1 C + 4 H)
        lines = xyz.strip().split("\n")
        assert len(lines) == 5

    def test_coordinate_format_is_valid(self):
        """Each line has element and 3 coordinates."""
        output_text = (FIXTURES_DIR / "nwchem_optimization_output.txt").read_text()
        xyz = extract_optimized_geometry(output_text)

        for line in xyz.strip().split("\n"):
            parts = line.split()
            assert len(parts) == 4, f"Expected 4 fields, got: {line}"

            element = parts[0]
            assert element.isalpha(), f"Element should be letters: {element}"
            assert len(element) <= 2, f"Element symbol too long: {element}"

            # Should be valid floats
            x, y, z = parts[1], parts[2], parts[3]
            float(x)  # Raises if invalid
            float(y)
            float(z)

    def test_extracts_correct_elements(self):
        """Elements match expected methane structure."""
        output_text = (FIXTURES_DIR / "nwchem_optimization_output.txt").read_text()
        xyz = extract_optimized_geometry(output_text)

        elements = [line.split()[0] for line in xyz.strip().split("\n")]
        assert elements.count("C") == 1
        assert elements.count("H") == 4

    def test_extracts_correct_coordinates(self):
        """Coordinates match expected values from fixture."""
        output_text = (FIXTURES_DIR / "nwchem_optimization_output.txt").read_text()
        xyz = extract_optimized_geometry(output_text)

        lines = xyz.strip().split("\n")

        # First atom is C at origin
        c_parts = lines[0].split()
        assert c_parts[0] == "C"
        assert float(c_parts[1]) == pytest.approx(0.0, abs=1e-6)
        assert float(c_parts[2]) == pytest.approx(0.0, abs=1e-6)
        assert float(c_parts[3]) == pytest.approx(0.0, abs=1e-6)

        # Second atom is H at (0.62918, 0.62918, 0.62918)
        h_parts = lines[1].split()
        assert h_parts[0] == "H"
        assert float(h_parts[1]) == pytest.approx(0.62918, abs=1e-6)

    def test_missing_section_raises_error(self):
        """RuntimeError raised when geometry section not found."""
        with pytest.raises(RuntimeError) as exc_info:
            extract_optimized_geometry("no geometry here at all")

        error_msg = str(exc_info.value).lower()
        assert "geometry" in error_msg or "coordinates" in error_msg

    def test_error_message_is_helpful(self):
        """Error message explains expected format."""
        with pytest.raises(RuntimeError) as exc_info:
            extract_optimized_geometry("some random output")

        error_msg = str(exc_info.value)
        assert "Output coordinates in angstroms" in error_msg


class TestParseShieldingOutput:
    """Tests for parse_shielding_output function."""

    def test_returns_expected_format(self):
        """Output dict has required keys: index, atom, shielding."""
        output_text = (FIXTURES_DIR / "nwchem_shielding_output.txt").read_text()
        result = parse_shielding_output(output_text)

        assert "index" in result
        assert "atom" in result
        assert "shielding" in result

    def test_lists_have_same_length(self):
        """All lists in result have same length."""
        output_text = (FIXTURES_DIR / "nwchem_shielding_output.txt").read_text()
        result = parse_shielding_output(output_text)

        n_atoms = len(result["index"])
        assert n_atoms > 0
        assert len(result["atom"]) == n_atoms
        assert len(result["shielding"]) == n_atoms

    def test_parses_correct_number_of_atoms(self):
        """Parses all 5 atoms in methane."""
        output_text = (FIXTURES_DIR / "nwchem_shielding_output.txt").read_text()
        result = parse_shielding_output(output_text)

        assert len(result["index"]) == 5

    def test_index_values_are_integers(self):
        """Atom indices are integers."""
        output_text = (FIXTURES_DIR / "nwchem_shielding_output.txt").read_text()
        result = parse_shielding_output(output_text)

        for idx in result["index"]:
            assert isinstance(idx, int)

    def test_atom_values_are_strings(self):
        """Element symbols are strings."""
        output_text = (FIXTURES_DIR / "nwchem_shielding_output.txt").read_text()
        result = parse_shielding_output(output_text)

        for atom in result["atom"]:
            assert isinstance(atom, str)
            assert atom.isalpha()

    def test_shielding_values_are_floats(self):
        """Shielding values are floats."""
        output_text = (FIXTURES_DIR / "nwchem_shielding_output.txt").read_text()
        result = parse_shielding_output(output_text)

        for shielding in result["shielding"]:
            assert isinstance(shielding, float)

    def test_extracts_correct_elements(self):
        """Extracts correct element types."""
        output_text = (FIXTURES_DIR / "nwchem_shielding_output.txt").read_text()
        result = parse_shielding_output(output_text)

        assert result["atom"].count("C") == 1
        assert result["atom"].count("H") == 4

    def test_extracts_correct_shielding_values(self):
        """Extracts correct isotropic shielding values."""
        output_text = (FIXTURES_DIR / "nwchem_shielding_output.txt").read_text()
        result = parse_shielding_output(output_text)

        # Carbon shielding should be around 195.9 ppm
        c_idx = result["atom"].index("C")
        c_shielding = result["shielding"][c_idx]
        assert c_shielding == pytest.approx(195.9423, abs=0.001)

        # Hydrogen shieldings should be around 31.2 ppm
        for i, atom in enumerate(result["atom"]):
            if atom == "H":
                assert result["shielding"][i] == pytest.approx(31.2, abs=0.1)

    def test_missing_section_raises_error(self):
        """RuntimeError raised when shielding section not found."""
        with pytest.raises(RuntimeError) as exc_info:
            parse_shielding_output("no shielding information here")

        error_msg = str(exc_info.value).lower()
        assert "shielding" in error_msg

    def test_error_message_is_helpful(self):
        """Error message explains expected format and task."""
        with pytest.raises(RuntimeError) as exc_info:
            parse_shielding_output("random text without shielding data")

        error_msg = str(exc_info.value)
        assert "GIAO" in error_msg or "shielding" in error_msg.lower()


class TestShieldingDataCompatibility:
    """Tests for compatibility with shifts.py module."""

    def test_compatible_with_shielding_to_shift(self):
        """Parsed shielding data works with shielding_to_shift function."""
        output_text = (FIXTURES_DIR / "nwchem_shielding_output.txt").read_text()
        shielding_data = parse_shielding_output(output_text)

        # Should not raise
        shifts = shielding_to_shift(shielding_data)

        # Should produce valid output
        assert "1H" in shifts
        assert "13C" in shifts

    def test_produces_h_shifts(self):
        """Produces 1H chemical shifts from methane."""
        output_text = (FIXTURES_DIR / "nwchem_shielding_output.txt").read_text()
        shielding_data = parse_shielding_output(output_text)
        shifts = shielding_to_shift(shielding_data)

        # Methane has 4 H atoms
        assert len(shifts["1H"]) == 4

        # Each shift dict has required keys
        for h_shift in shifts["1H"]:
            assert "index" in h_shift
            assert "atom" in h_shift
            assert "shielding" in h_shift
            assert "shift" in h_shift

    def test_produces_c_shifts(self):
        """Produces 13C chemical shifts from methane."""
        output_text = (FIXTURES_DIR / "nwchem_shielding_output.txt").read_text()
        shielding_data = parse_shielding_output(output_text)
        shifts = shielding_to_shift(shielding_data)

        # Methane has 1 C atom
        assert len(shifts["13C"]) == 1

        # The shift should be a reasonable value for methane carbon
        # Experimental methane 13C shift is around -2.3 ppm
        # Our calculated value depends on the scaling factors
        c_shift = shifts["13C"][0]["shift"]
        assert isinstance(c_shift, float)

    def test_shift_values_are_reasonable(self):
        """Chemical shift values are in reasonable range."""
        output_text = (FIXTURES_DIR / "nwchem_shielding_output.txt").read_text()
        shielding_data = parse_shielding_output(output_text)
        shifts = shielding_to_shift(shielding_data)

        # 1H shifts typically 0-15 ppm for organic molecules
        for h_shift in shifts["1H"]:
            assert -5 < h_shift["shift"] < 20

        # 13C shifts typically -20 to 220 ppm
        for c_shift in shifts["13C"]:
            assert -50 < c_shift["shift"] < 250
