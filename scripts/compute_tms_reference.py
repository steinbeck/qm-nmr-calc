#!/usr/bin/env python3
"""Compute TMS reference shieldings for each preset.

Run this once to get TMS shielding values for calibrating chemical shifts.
"""

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import tempfile
import json
from qm_nmr_calc.isicle_wrapper import run_geometry_opt_dft, run_nmr_shielding
from qm_nmr_calc.presets import PRESETS, PresetName

TMS_SMILES = "C[Si](C)(C)C"
SOLVENT = "chcl3"  # CDCl3 - most common NMR solvent


def compute_tms_for_preset(preset_name: PresetName) -> dict:
    """Compute TMS shieldings for a given preset."""
    preset = PRESETS[preset_name]
    print(f"\n{'='*60}")
    print(f"Computing TMS for preset: {preset_name.value}")
    print(f"  Functional: {preset['functional']}")
    print(f"  Geometry basis: {preset['basis_set']}")
    print(f"  NMR basis: {preset['nmr_basis_set']}")
    print(f"  Solvent: {SOLVENT}")
    print(f"  Processes: {preset['processes']}")
    print(f"{'='*60}\n")

    with tempfile.TemporaryDirectory() as tmpdir:
        job_dir = Path(tmpdir)

        print("Step 1: Geometry optimization...")
        geometry_file, optimized_geom = run_geometry_opt_dft(
            smiles=TMS_SMILES,
            job_dir=job_dir,
            preset=preset,
            solvent=SOLVENT,
            processes=preset['processes'],
        )
        print(f"  Done. Geometry saved to {geometry_file}")

        print("\nStep 2: NMR shielding calculation...")
        nmr_result = run_nmr_shielding(
            optimized_geom=optimized_geom,
            job_dir=job_dir,
            preset=preset,
            solvent=SOLVENT,
            processes=preset['processes'],
        )
        print("  Done.")

        # Extract shielding values
        shielding_data = nmr_result['shielding_data']

        # Separate by atom type
        h_shieldings = []
        c_shieldings = []

        for idx, atom, shield in zip(
            shielding_data['index'],
            shielding_data['atom'],
            shielding_data['shielding'],
        ):
            if atom == 'H':
                h_shieldings.append(shield)
            elif atom == 'C':
                c_shieldings.append(shield)

        # TMS should have 12 equivalent H and 4 equivalent C
        # Average them (they should be very close)
        h_avg = sum(h_shieldings) / len(h_shieldings)
        c_avg = sum(c_shieldings) / len(c_shieldings)

        print(f"\nResults for {preset_name.value}:")
        print(f"  1H shieldings: {[round(s, 2) for s in h_shieldings]}")
        print(f"  1H average: {h_avg:.4f} ppm")
        print(f"  13C shieldings: {[round(s, 2) for s in c_shieldings]}")
        print(f"  13C average: {c_avg:.4f} ppm")

        return {
            "preset": preset_name.value,
            "functional": preset['functional'],
            "geometry_basis": preset['basis_set'],
            "nmr_basis": preset['nmr_basis_set'],
            "solvent": SOLVENT,
            "h1_shieldings": h_shieldings,
            "h1_reference": round(h_avg, 2),
            "c13_shieldings": c_shieldings,
            "c13_reference": round(c_avg, 2),
        }


def main():
    print("Computing TMS reference shieldings")
    print(f"TMS SMILES: {TMS_SMILES}")

    results = {}

    for preset_name in [PresetName.DRAFT, PresetName.PRODUCTION]:
        results[preset_name.value] = compute_tms_for_preset(preset_name)

    # Summary
    print("\n" + "="*60)
    print("SUMMARY - TMS Reference Shieldings")
    print("="*60)
    print("\nUse these values in shifts.py:\n")

    for name, data in results.items():
        print(f"{name.upper()} preset ({data['nmr_basis']}):")
        print(f"  1H TMS reference:  {data['h1_reference']} ppm")
        print(f"  13C TMS reference: {data['c13_reference']} ppm")
        print()

    # Save to JSON for reference
    output_file = Path(__file__).parent.parent / "data" / "tms_reference_values.json"
    output_file.parent.mkdir(exist_ok=True)
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"Full results saved to: {output_file}")


if __name__ == "__main__":
    main()
