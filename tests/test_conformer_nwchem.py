"""Tests for conformer DFT optimization loop and post-DFT filtering.

All tests use mocks - no actual NWChem execution.
"""

import tempfile
import unittest.mock
from pathlib import Path

import pytest

from qm_nmr_calc.models import ConformerData, ConformerEnsemble


# Fixtures for test data creation


def make_ensemble(n: int = 3) -> ConformerEnsemble:
    """Create a ConformerEnsemble with n conformers for testing.

    Args:
        n: Number of conformers to create

    Returns:
        ConformerEnsemble with MMFF energies in kcal/mol
    """
    conformers = []
    for i in range(n):
        conf_id = f"conf_{i+1:03d}"
        conformers.append(
            ConformerData(
                conformer_id=conf_id,
                energy=10.0 + i * 0.5,  # MMFF energies in kcal/mol
                energy_unit="kcal_mol",
                geometry_file=f"output/conformers/{conf_id}/geometry.xyz",
                status="pending",
            )
        )
    return ConformerEnsemble(
        method="rdkit_kdg",
        conformers=conformers,
        temperature_k=298.15,
    )


def make_preset() -> dict:
    """Create a test calculation preset."""
    return {
        "functional": "b3lyp",
        "basis_set": "6-31G*",
        "nmr_basis_set": "6-311+G(2d,p)",
        "max_iter": 150,
    }


# Tests for run_calculation scratch_dir_override


def test_run_calculation_uses_scratch_dir_override():
    """Test that run_calculation uses scratch_dir_override parameter."""
    from qm_nmr_calc.nwchem.runner import run_calculation

    with tempfile.TemporaryDirectory() as tmpdir:
        job_dir = Path(tmpdir) / "job"
        job_dir.mkdir()
        scratch_override = Path(tmpdir) / "custom_scratch"
        scratch_override.mkdir()

        # Create a dummy geometry file
        geom_file = job_dir / "test.xyz"
        geom_file.write_text(
            "2\nTest molecule\nC  0.0  0.0  0.0\nH  1.0  0.0  0.0\n"
        )

        # Mock run_nwchem to capture where input files are written
        with unittest.mock.patch(
            "qm_nmr_calc.nwchem.runner.run_nwchem"
        ) as mock_run:
            # Set up mock to create expected output files
            def side_effect(input_file, output_file, processes):
                # Create output files in the scratch directory
                output_file.write_text(
                    "Total DFT energy =     -40.5186\n"
                    "Output coordinates in angstroms\n"
                    "No.       Tag         Charge          X              Y              Z\n"
                    "---- ---------------- ---------- -------------- -------------- --------------\n"
                    "   1 C                    6.0000     0.00000000     0.00000000     0.00000000\n"
                    "   2 H                    1.0000     1.09000000     0.00000000     0.00000000\n"
                    "\n"
                    "GIAO Chemical Shielding Tensors\n"
                    "Atom:    1  C\n"
                    "   isotropic  =     183.4567\n"
                    "Atom:    2  H\n"
                    "   isotropic  =     29.1234\n"
                )

            mock_run.side_effect = side_effect

            # Call with scratch_dir_override
            result = run_calculation(
                smiles="CC",
                job_dir=job_dir,
                preset=make_preset(),
                solvent="chcl3",
                processes=2,
                scratch_dir_override=scratch_override,
            )

            # Verify run_nwchem was called with input files in override directory
            assert mock_run.call_count >= 1
            first_call_input_file = mock_run.call_args_list[0][0][0]
            # Input file should be in scratch_override, not job_dir/scratch
            assert str(first_call_input_file).startswith(str(scratch_override))
            assert "scratch" not in str(first_call_input_file).split(str(scratch_override))[1].split("/")[1]


# Tests for run_conformer_dft_optimization


def test_conformer_dft_optimization_all_succeed():
    """Test successful DFT optimization of all conformers."""
    from qm_nmr_calc.nwchem.runner import run_conformer_dft_optimization

    ensemble = make_ensemble(n=3)

    with tempfile.TemporaryDirectory() as tmpdir:
        job_dir = Path(tmpdir) / "job"
        job_dir.mkdir()
        (job_dir / "output" / "conformers").mkdir(parents=True)

        # Create geometry files
        for conf in ensemble.conformers:
            geom_path = job_dir / conf.geometry_file
            geom_path.parent.mkdir(parents=True, exist_ok=True)
            geom_path.write_text("2\nTest\nC 0 0 0\nH 1 0 0\n")

        # Mock run_calculation to return success for all conformers
        with unittest.mock.patch(
            "qm_nmr_calc.nwchem.runner.run_calculation"
        ) as mock_run:
            # Create temp files for optimization output
            opt_output = Path(tmpdir) / "optimize.out"
            opt_output.write_text("Total DFT energy =     -40.5186\n")
            geom_output = Path(tmpdir) / "optimized.xyz"
            geom_output.write_text("2\nOptimized\nC 0 0 0\nH 1.1 0 0\n")

            mock_run.return_value = {
                "optimization_output": opt_output,
                "geometry_file": geom_output,
                "shielding_data": {},
                "shielding_output": Path(tmpdir) / "shielding.out",
            }

            # Mock storage functions
            with unittest.mock.patch(
                "qm_nmr_calc.nwchem.runner.get_job_dir", return_value=job_dir
            ):
                with unittest.mock.patch(
                    "qm_nmr_calc.nwchem.runner.get_conformer_scratch_dir",
                    side_effect=lambda job_id, conf_id: job_dir / "scratch" / conf_id,
                ):
                    # Run optimization
                    successful, failed = run_conformer_dft_optimization(
                        ensemble=ensemble,
                        job_id="test_job",
                        preset=make_preset(),
                        solvent="chcl3",
                        processes=2,
                    )

        # Verify results
        assert len(successful) == 3
        assert len(failed) == 0

        for conf in successful:
            assert conf.status == "optimized"
            assert conf.energy is not None
            assert isinstance(conf.energy, float)
            assert conf.energy_unit == "hartree"
            assert conf.optimized_geometry_file is not None


def test_conformer_dft_optimization_partial_failure():
    """Test DFT optimization with one conformer failing."""
    from qm_nmr_calc.nwchem.runner import run_conformer_dft_optimization

    ensemble = make_ensemble(n=4)

    with tempfile.TemporaryDirectory() as tmpdir:
        job_dir = Path(tmpdir) / "job"
        job_dir.mkdir()

        # Create geometry files
        for conf in ensemble.conformers:
            geom_path = job_dir / conf.geometry_file
            geom_path.parent.mkdir(parents=True, exist_ok=True)
            geom_path.write_text("2\nTest\nC 0 0 0\nH 1 0 0\n")

        # Mock run_calculation to fail for conformer index 1 (conf_002)
        call_count = [0]

        def mock_run_calculation(*args, **kwargs):
            idx = call_count[0]
            call_count[0] += 1

            if idx == 1:  # Second conformer (conf_002) fails
                raise RuntimeError("NWChem optimization failed for this conformer")

            # Success case
            opt_output = Path(tmpdir) / f"optimize_{idx}.out"
            opt_output.write_text("Total DFT energy =     -40.5186\n")
            geom_output = Path(tmpdir) / f"optimized_{idx}.xyz"
            geom_output.write_text("2\nOptimized\nC 0 0 0\nH 1.1 0 0\n")

            return {
                "optimization_output": opt_output,
                "geometry_file": geom_output,
                "shielding_data": {},
                "shielding_output": Path(tmpdir) / f"shielding_{idx}.out",
            }

        with unittest.mock.patch(
            "qm_nmr_calc.nwchem.runner.run_calculation",
            side_effect=mock_run_calculation,
        ):
            with unittest.mock.patch(
                "qm_nmr_calc.nwchem.runner.get_job_dir", return_value=job_dir
            ):
                with unittest.mock.patch(
                    "qm_nmr_calc.nwchem.runner.get_conformer_scratch_dir",
                    side_effect=lambda job_id, conf_id: job_dir / "scratch" / conf_id,
                ):
                    # Run optimization
                    successful, failed = run_conformer_dft_optimization(
                        ensemble=ensemble,
                        job_id="test_job",
                        preset=make_preset(),
                        solvent="chcl3",
                        processes=2,
                    )

        # Verify results
        assert len(successful) == 3
        assert len(failed) == 1

        # Check failed conformer
        assert failed[0].status == "failed"
        assert failed[0].error_message is not None
        assert "NWChem optimization failed" in failed[0].error_message

        # Check successful conformers
        for conf in successful:
            assert conf.status == "optimized"


def test_conformer_dft_optimization_too_many_failures():
    """Test that RuntimeError is raised when >50% of conformers fail."""
    from qm_nmr_calc.nwchem.runner import run_conformer_dft_optimization

    ensemble = make_ensemble(n=4)

    with tempfile.TemporaryDirectory() as tmpdir:
        job_dir = Path(tmpdir) / "job"
        job_dir.mkdir()

        # Create geometry files
        for conf in ensemble.conformers:
            geom_path = job_dir / conf.geometry_file
            geom_path.parent.mkdir(parents=True, exist_ok=True)
            geom_path.write_text("2\nTest\nC 0 0 0\nH 1 0 0\n")

        # Mock run_calculation to fail for 3 out of 4 conformers
        call_count = [0]

        def mock_run_calculation(*args, **kwargs):
            idx = call_count[0]
            call_count[0] += 1

            if idx != 0:  # Only first conformer succeeds
                raise RuntimeError("NWChem optimization failed")

            # Success case
            opt_output = Path(tmpdir) / "optimize.out"
            opt_output.write_text("Total DFT energy =     -40.5186\n")
            geom_output = Path(tmpdir) / "optimized.xyz"
            geom_output.write_text("2\nOptimized\nC 0 0 0\nH 1.1 0 0\n")

            return {
                "optimization_output": opt_output,
                "geometry_file": geom_output,
                "shielding_data": {},
                "shielding_output": Path(tmpdir) / "shielding.out",
            }

        with unittest.mock.patch(
            "qm_nmr_calc.nwchem.runner.run_calculation",
            side_effect=mock_run_calculation,
        ):
            with unittest.mock.patch(
                "qm_nmr_calc.nwchem.runner.get_job_dir", return_value=job_dir
            ):
                with unittest.mock.patch(
                    "qm_nmr_calc.nwchem.runner.get_conformer_scratch_dir",
                    side_effect=lambda job_id, conf_id: job_dir / "scratch" / conf_id,
                ):
                    # Run optimization - should raise RuntimeError
                    with pytest.raises(RuntimeError, match="success rate"):
                        run_conformer_dft_optimization(
                            ensemble=ensemble,
                            job_id="test_job",
                            preset=make_preset(),
                            solvent="chcl3",
                            processes=2,
                        )


def test_conformer_dft_optimization_status_tracking():
    """Test that conformer status transitions correctly through lifecycle."""
    from qm_nmr_calc.nwchem.runner import run_conformer_dft_optimization

    ensemble = make_ensemble(n=2)

    with tempfile.TemporaryDirectory() as tmpdir:
        job_dir = Path(tmpdir) / "job"
        job_dir.mkdir()

        # Create geometry files
        for conf in ensemble.conformers:
            geom_path = job_dir / conf.geometry_file
            geom_path.parent.mkdir(parents=True, exist_ok=True)
            geom_path.write_text("2\nTest\nC 0 0 0\nH 1 0 0\n")

        # Track status transitions
        status_transitions = []

        original_run_calculation = None

        def mock_run_calculation(*args, **kwargs):
            # Capture status of all conformers at this point
            for conf in ensemble.conformers:
                status_transitions.append((conf.conformer_id, conf.status))

            # Return success
            opt_output = Path(tmpdir) / "optimize.out"
            opt_output.write_text("Total DFT energy =     -40.5186\n")
            geom_output = Path(tmpdir) / "optimized.xyz"
            geom_output.write_text("2\nOptimized\nC 0 0 0\nH 1.1 0 0\n")

            return {
                "optimization_output": opt_output,
                "geometry_file": geom_output,
                "shielding_data": {},
                "shielding_output": Path(tmpdir) / "shielding.out",
            }

        with unittest.mock.patch(
            "qm_nmr_calc.nwchem.runner.run_calculation",
            side_effect=mock_run_calculation,
        ):
            with unittest.mock.patch(
                "qm_nmr_calc.nwchem.runner.get_job_dir", return_value=job_dir
            ):
                with unittest.mock.patch(
                    "qm_nmr_calc.nwchem.runner.get_conformer_scratch_dir",
                    side_effect=lambda job_id, conf_id: job_dir / "scratch" / conf_id,
                ):
                    # Run optimization
                    successful, failed = run_conformer_dft_optimization(
                        ensemble=ensemble,
                        job_id="test_job",
                        preset=make_preset(),
                        solvent="chcl3",
                        processes=2,
                    )

        # Verify final status
        assert all(conf.status == "optimized" for conf in successful)

        # Verify status transitions occurred
        # Should have transitions to "optimizing" for each conformer
        optimizing_transitions = [
            (cid, status)
            for cid, status in status_transitions
            if status == "optimizing"
        ]
        assert len(optimizing_transitions) >= 2


def test_conformer_dft_optimization_stores_optimized_geometry_path():
    """Test that optimized_geometry_file is set correctly."""
    from qm_nmr_calc.nwchem.runner import run_conformer_dft_optimization

    ensemble = make_ensemble(n=1)

    with tempfile.TemporaryDirectory() as tmpdir:
        job_dir = Path(tmpdir) / "job"
        job_dir.mkdir()

        # Create geometry file
        geom_path = job_dir / ensemble.conformers[0].geometry_file
        geom_path.parent.mkdir(parents=True, exist_ok=True)
        geom_path.write_text("2\nTest\nC 0 0 0\nH 1 0 0\n")

        # Mock run_calculation
        geom_output = job_dir / "output" / "optimized" / "conf_001_optimized.xyz"
        geom_output.parent.mkdir(parents=True, exist_ok=True)

        with unittest.mock.patch(
            "qm_nmr_calc.nwchem.runner.run_calculation"
        ) as mock_run:
            opt_output = Path(tmpdir) / "optimize.out"
            opt_output.write_text("Total DFT energy =     -40.5186\n")

            mock_run.return_value = {
                "optimization_output": opt_output,
                "geometry_file": geom_output,
                "shielding_data": {},
                "shielding_output": Path(tmpdir) / "shielding.out",
            }

            with unittest.mock.patch(
                "qm_nmr_calc.nwchem.runner.get_job_dir", return_value=job_dir
            ):
                with unittest.mock.patch(
                    "qm_nmr_calc.nwchem.runner.get_conformer_scratch_dir",
                    side_effect=lambda job_id, conf_id: job_dir / "scratch" / conf_id,
                ):
                    # Run optimization
                    successful, failed = run_conformer_dft_optimization(
                        ensemble=ensemble,
                        job_id="test_job",
                        preset=make_preset(),
                        solvent="chcl3",
                        processes=2,
                    )

        # Verify optimized_geometry_file is set and is relative to job_dir
        assert successful[0].optimized_geometry_file is not None
        assert successful[0].optimized_geometry_file.startswith("output/")


# Tests for apply_post_dft_filter


def test_apply_post_dft_filter_within_window():
    """Test that conformers within energy window are kept."""
    from qm_nmr_calc.nwchem.runner import apply_post_dft_filter

    # Create conformers with DFT energies in Hartree
    # Energy differences in milliHartree that are within 3 kcal/mol
    # 3 kcal/mol = 0.00478 Hartree
    conformers = [
        ConformerData(
            conformer_id="conf_001",
            energy=-40.5186,  # Minimum
            energy_unit="hartree",
            status="optimized",
        ),
        ConformerData(
            conformer_id="conf_002",
            energy=-40.5180,  # +0.0006 Ha = 0.38 kcal/mol
            energy_unit="hartree",
            status="optimized",
        ),
        ConformerData(
            conformer_id="conf_003",
            energy=-40.5170,  # +0.0016 Ha = 1.00 kcal/mol
            energy_unit="hartree",
            status="optimized",
        ),
    ]

    # Apply filter with default 3 kcal/mol window
    filtered = apply_post_dft_filter(conformers, window_kcal=3.0)

    # All should be kept
    assert len(filtered) == 3


def test_apply_post_dft_filter_removes_high_energy():
    """Test that conformers outside energy window are removed."""
    from qm_nmr_calc.nwchem.runner import apply_post_dft_filter

    # Create conformers where one is outside 3 kcal/mol window
    # 0.01 Hartree = 6.28 kcal/mol (outside 3 kcal/mol window)
    conformers = [
        ConformerData(
            conformer_id="conf_001",
            energy=-40.5186,  # Minimum
            energy_unit="hartree",
            status="optimized",
        ),
        ConformerData(
            conformer_id="conf_002",
            energy=-40.5180,  # +0.0006 Ha = 0.38 kcal/mol
            energy_unit="hartree",
            status="optimized",
        ),
        ConformerData(
            conformer_id="conf_003",
            energy=-40.5170,  # +0.0016 Ha = 1.00 kcal/mol
            energy_unit="hartree",
            status="optimized",
        ),
        ConformerData(
            conformer_id="conf_004",
            energy=-40.5086,  # +0.01 Ha = 6.28 kcal/mol (OUTSIDE)
            energy_unit="hartree",
            status="optimized",
        ),
    ]

    # Apply filter with default 3 kcal/mol window
    filtered = apply_post_dft_filter(conformers, window_kcal=3.0)

    # Only first 3 should be kept
    assert len(filtered) == 3
    assert all(c.conformer_id != "conf_004" for c in filtered)


def test_apply_post_dft_filter_custom_window():
    """Test post-DFT filter with custom energy window."""
    from qm_nmr_calc.nwchem.runner import apply_post_dft_filter

    # Create conformers where one is between 3 and 5 kcal/mol
    conformers = [
        ConformerData(
            conformer_id="conf_001",
            energy=-40.5186,  # Minimum
            energy_unit="hartree",
            status="optimized",
        ),
        ConformerData(
            conformer_id="conf_002",
            energy=-40.5106,  # +0.008 Ha = 5.02 kcal/mol
            energy_unit="hartree",
            status="optimized",
        ),
    ]

    # With 3 kcal/mol window, second is excluded
    filtered_3 = apply_post_dft_filter(conformers, window_kcal=3.0)
    assert len(filtered_3) == 1

    # With 6 kcal/mol window, second is included
    filtered_6 = apply_post_dft_filter(conformers, window_kcal=6.0)
    assert len(filtered_6) == 2


def test_apply_post_dft_filter_single_conformer():
    """Test that single conformer always passes filter."""
    from qm_nmr_calc.nwchem.runner import apply_post_dft_filter

    conformers = [
        ConformerData(
            conformer_id="conf_001",
            energy=-40.5186,
            energy_unit="hartree",
            status="optimized",
        )
    ]

    filtered = apply_post_dft_filter(conformers, window_kcal=3.0)
    assert len(filtered) == 1
