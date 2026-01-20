"""Visualization functions for NMR spectrum plots and annotated structures."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # Non-interactive backend - MUST be before pyplot import
import matplotlib.pyplot as plt

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D


def generate_spectrum_plot(
    shifts: list[float],
    nucleus: str,
    output_dir: Path,
) -> tuple[Path, Path]:
    """Generate stick spectrum plot in SVG and PNG formats.

    Args:
        shifts: List of chemical shift values in ppm
        nucleus: "1H" or "13C" for labeling
        output_dir: Directory to write output files

    Returns:
        Tuple of (svg_path, png_path)
    """
    fig, ax = plt.subplots(figsize=(8, 4))

    # Stick plot: vertical lines at each shift, no markers, no baseline
    ax.stem(shifts, [1.0] * len(shifts), linefmt="k-", markerfmt=" ", basefmt=" ")

    # NMR convention: x-axis right to left (high ppm to low ppm)
    ax.invert_xaxis()

    # Minimal labeling
    ax.set_xlabel("ppm")
    ax.set_yticks([])  # No y-axis ticks (intensity is relative)

    # Add nucleus label in top-left corner
    ax.text(
        0.02,
        0.95,
        nucleus,
        transform=ax.transAxes,
        fontsize=12,
        verticalalignment="top",
    )

    # Auto-fit x-axis with padding
    if shifts:
        ppm_range = max(shifts) - min(shifts)
        padding = max(ppm_range * 0.1, 0.5)  # At least 0.5 ppm padding
        ax.set_xlim(max(shifts) + padding, min(shifts) - padding)

    # Save as both formats
    base_name = f"spectrum_{nucleus}"
    svg_path = output_dir / f"{base_name}.svg"
    png_path = output_dir / f"{base_name}.png"

    fig.savefig(svg_path, format="svg", bbox_inches="tight")
    fig.savefig(png_path, format="png", dpi=300, bbox_inches="tight")

    # CRITICAL: Close figure to prevent memory leaks
    plt.close(fig)

    return svg_path, png_path


def generate_annotated_structure(
    smiles: str,
    h1_shifts: list[dict],
    c13_shifts: list[dict],
    output_dir: Path,
) -> tuple[Path, Path]:
    """Generate structure with chemical shift annotations.

    Args:
        smiles: SMILES string for the molecule
        h1_shifts: List of {"index": int, "shift": float} for 1H (NWChem 1-based indices)
        c13_shifts: List of {"index": int, "shift": float} for 13C (NWChem 1-based indices)
        output_dir: Directory to write output files

    Returns:
        Tuple of (svg_path, png_path)
    """
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    # Set atomNote on each atom BEFORE calling PrepareMolForDrawing
    # CRITICAL: Convert NWChem 1-based indices to RDKit 0-based
    for shift_data in h1_shifts + c13_shifts:
        rdkit_idx = shift_data["index"] - 1  # Convert 1-based to 0-based
        atom = mol.GetAtomWithIdx(rdkit_idx)
        atom.SetProp("atomNote", f"{shift_data['shift']:.2f}")

    # Prepare molecule for drawing (after setting notes)
    rdMolDraw2D.PrepareMolForDrawing(mol)

    # Generate SVG
    d_svg = rdMolDraw2D.MolDraw2DSVG(800, 600)
    d_svg.drawOptions().annotationFontScale = 0.9
    d_svg.DrawMolecule(mol)
    d_svg.FinishDrawing()

    svg_path = output_dir / "structure_annotated.svg"
    svg_path.write_text(d_svg.GetDrawingText())

    # Generate PNG (300 DPI equivalent: 800x600 * 3 = 2400x1800)
    d_png = rdMolDraw2D.MolDraw2DCairo(2400, 1800)
    d_png.drawOptions().annotationFontScale = 0.9
    d_png.DrawMolecule(mol)
    d_png.FinishDrawing()

    png_path = output_dir / "structure_annotated.png"
    d_png.WriteDrawingText(str(png_path))

    return svg_path, png_path
