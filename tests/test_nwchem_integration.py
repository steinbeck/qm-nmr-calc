"""Integration tests for NWChem execution.

These tests actually run NWChem to validate the full pipeline:
1. Generate input files
2. Execute NWChem
3. Parse output files
4. Verify results

Skip these tests if NWChem is not available (CI environments).

Run with: pytest tests/test_nwchem_integration.py -v
"""

import shutil
import subprocess
import tempfile
from pathlib import Path

import pytest

from qm_nmr_calc.nwchem import (
    extract_optimized_geometry,
    generate_optimization_input,
    generate_shielding_input,
    parse_shielding_output,
    smiles_to_xyz,
)
from qm_nmr_calc.shifts import shielding_to_shift


def nwchem_available() -> bool:
    """Check if NWChem is available in PATH."""
    return shutil.which("nwchem") is not None


# Skip all tests in this module if NWChem not available
pytestmark = pytest.mark.skipif(
    not nwchem_available(),
    reason="NWChem not available - install nwchem to run integration tests"
)


def run_nwchem(input_file: Path, output_file: Path, processes: int = 2) -> int:
    """Run NWChem and return exit code."""
    # Try mpirun first, fall back to direct nwchem
    if shutil.which("mpirun"):
        cmd = ["mpirun", "-np", str(processes), "nwchem", str(input_file)]
    else:
        cmd = ["nwchem", str(input_file)]

    with open(output_file, "w") as f:
        result = subprocess.run(
            cmd,
            stdout=f,
            stderr=subprocess.STDOUT,
            timeout=300,  # 5 minute timeout
        )
    return result.returncode


class TestNWChemGeometryOptimization:
    """Test geometry optimization with actual NWChem execution."""

    def test_methane_optimization(self, tmp_path):
        """Optimize methane geometry and verify output is parseable."""
        # Generate input
        mol, xyz_block = smiles_to_xyz("C")  # Methane
        input_text = generate_optimization_input(
            geometry_xyz=xyz_block,
            functional="b3lyp",
            basis_set="6-31G*",
            solvent="chcl3",
            max_iter=50,
        )

        # Write input file
        input_file = tmp_path / "optimize.nw"
        input_file.write_text(input_text)

        # Run NWChem
        output_file = tmp_path / "optimize.out"
        exit_code = run_nwchem(input_file, output_file)

        # Check execution succeeded
        assert exit_code == 0, f"NWChem failed with exit code {exit_code}. Output:\n{output_file.read_text()[-2000:]}"

        # Parse output
        output_text = output_file.read_text()
        optimized_xyz = extract_optimized_geometry(output_text)

        # Verify we got 5 atoms (1 C + 4 H)
        lines = [l for l in optimized_xyz.strip().split("\n") if l.strip()]
        assert len(lines) == 5, f"Expected 5 atoms, got {len(lines)}"

        # Verify elements
        elements = [line.split()[0] for line in lines]
        assert elements.count("C") == 1
        assert elements.count("H") == 4

    def test_ethanol_optimization(self, tmp_path):
        """Optimize ethanol geometry - slightly larger molecule."""
        # Generate input
        mol, xyz_block = smiles_to_xyz("CCO")  # Ethanol
        input_text = generate_optimization_input(
            geometry_xyz=xyz_block,
            functional="b3lyp",
            basis_set="6-31G*",
            solvent="dmso",
            max_iter=100,
        )

        # Write and run
        input_file = tmp_path / "optimize.nw"
        input_file.write_text(input_text)
        output_file = tmp_path / "optimize.out"
        exit_code = run_nwchem(input_file, output_file)

        assert exit_code == 0, f"NWChem failed. Output tail:\n{output_file.read_text()[-2000:]}"

        # Parse and verify
        output_text = output_file.read_text()
        optimized_xyz = extract_optimized_geometry(output_text)

        lines = [l for l in optimized_xyz.strip().split("\n") if l.strip()]
        assert len(lines) == 9, f"Ethanol should have 9 atoms (2C + 1O + 6H), got {len(lines)}"


class TestNWChemNMRShielding:
    """Test NMR shielding calculation with actual NWChem execution."""

    def test_methane_shielding(self, tmp_path):
        """Calculate NMR shielding for methane and verify parsing."""
        # First optimize geometry
        mol, xyz_block = smiles_to_xyz("C")
        opt_input = generate_optimization_input(
            geometry_xyz=xyz_block,
            functional="b3lyp",
            basis_set="6-31G*",
            solvent="chcl3",
            max_iter=50,
        )

        opt_input_file = tmp_path / "optimize.nw"
        opt_input_file.write_text(opt_input)
        opt_output_file = tmp_path / "optimize.out"

        exit_code = run_nwchem(opt_input_file, opt_output_file)
        assert exit_code == 0, "Optimization failed"

        # Get optimized geometry
        optimized_xyz = extract_optimized_geometry(opt_output_file.read_text())

        # Now run NMR shielding
        nmr_input = generate_shielding_input(
            geometry_xyz=optimized_xyz,
            functional="b3lyp",
            basis_set="6-311+G(2d,p)",
            solvent="chcl3",
        )

        nmr_input_file = tmp_path / "shielding.nw"
        nmr_input_file.write_text(nmr_input)
        nmr_output_file = tmp_path / "shielding.out"

        exit_code = run_nwchem(nmr_input_file, nmr_output_file)
        assert exit_code == 0, f"NMR calculation failed. Output tail:\n{nmr_output_file.read_text()[-2000:]}"

        # Parse shielding
        shielding_data = parse_shielding_output(nmr_output_file.read_text())

        # Verify structure
        assert len(shielding_data["index"]) == 5
        assert shielding_data["atom"].count("C") == 1
        assert shielding_data["atom"].count("H") == 4

        # Verify shielding values are reasonable
        for shielding in shielding_data["shielding"]:
            assert isinstance(shielding, float)
            assert -100 < shielding < 300  # Reasonable range for organic molecules

    def test_shielding_to_shift_conversion(self, tmp_path):
        """Full pipeline: SMILES -> optimization -> shielding -> chemical shifts."""
        # Optimize
        mol, xyz_block = smiles_to_xyz("C")
        opt_input = generate_optimization_input(xyz_block, "b3lyp", "6-31G*", "chcl3", 50)

        opt_file = tmp_path / "opt.nw"
        opt_file.write_text(opt_input)
        opt_out = tmp_path / "opt.out"

        assert run_nwchem(opt_file, opt_out) == 0
        optimized_xyz = extract_optimized_geometry(opt_out.read_text())

        # NMR
        nmr_input = generate_shielding_input(optimized_xyz, "b3lyp", "6-311+G(2d,p)", "chcl3")
        nmr_file = tmp_path / "nmr.nw"
        nmr_file.write_text(nmr_input)
        nmr_out = tmp_path / "nmr.out"

        assert run_nwchem(nmr_file, nmr_out) == 0
        shielding_data = parse_shielding_output(nmr_out.read_text())

        # Convert to chemical shifts
        shifts = shielding_to_shift(shielding_data, "B3LYP", "6-311+G(2d,p)", "CHCl3")

        # Verify we got shifts
        assert len(shifts["1H"]) == 4, "Should have 4 H shifts for methane"
        assert len(shifts["13C"]) == 1, "Should have 1 C shift for methane"

        # Methane H shifts should be around 0-1 ppm experimentally
        # Our calculated values may differ but should be in reasonable range
        for h in shifts["1H"]:
            assert -5 < h["shift"] < 10, f"H shift {h['shift']} out of expected range"

        # Methane C shift is around -2.3 ppm experimentally
        c_shift = shifts["13C"][0]["shift"]
        assert -30 < c_shift < 30, f"C shift {c_shift} out of expected range"

        print(f"\n=== Methane NMR Results ===")
        print(f"1H shifts: {[round(h['shift'], 2) for h in shifts['1H']]}")
        print(f"13C shift: {round(c_shift, 2)}")


class TestCOSMOSolvation:
    """Verify COSMO solvation is actually being applied."""

    def test_cosmo_in_output(self, tmp_path):
        """Verify COSMO appears in NWChem output."""
        mol, xyz_block = smiles_to_xyz("C")
        input_text = generate_optimization_input(xyz_block, "b3lyp", "6-31G*", "chcl3", 30)

        input_file = tmp_path / "cosmo_test.nw"
        input_file.write_text(input_text)
        output_file = tmp_path / "cosmo_test.out"

        exit_code = run_nwchem(input_file, output_file)
        assert exit_code == 0

        output_text = output_file.read_text()

        # Verify COSMO was recognized and applied
        assert "cosmo" in output_text.lower(), "COSMO not mentioned in output"
        # Look for dielectric constant in output
        assert "4.8" in output_text or "dielectric" in output_text.lower(), \
            "Dielectric constant not found in output"
