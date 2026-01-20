# Phase 5: Visualization - Research

**Researched:** 2026-01-20
**Domain:** Python scientific visualization (matplotlib, RDKit molecular drawing)
**Confidence:** HIGH

## Summary

Phase 5 requires generating two types of visualizations: (1) stick spectrum plots showing NMR chemical shifts as vertical lines, and (2) annotated molecular structure images with shift values labeled on atoms. Both outputs need to be provided in PNG (300 DPI) and SVG formats.

Research confirms the standard approach:
- **Spectrum plots**: Use matplotlib's `stem()` function with customization to hide markers/baseline, creating clean stick spectra. Matplotlib handles SVG/PNG export natively with `dpi` parameter.
- **Structure annotation**: RDKit's `atomNote` property combined with `MolDraw2DSVG` and `MolDraw2DCairo` provides native annotation capabilities without needing external compositing.

**Primary recommendation:** Use matplotlib for spectrum plots and RDKit's native annotation for structure images. Both libraries handle SVG/PNG export directly, avoiding complex format conversion pipelines.

## Standard Stack

The established libraries/tools for this domain:

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| matplotlib | 3.10.x | Spectrum stick plots | De facto Python plotting standard, excellent SVG/PNG export |
| rdkit | 2025.09.x | Structure annotation | Already in project, has built-in atom annotation support |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| pillow | 12.1.0 | Image manipulation (already installed) | Already a dependency, no extra install needed |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| matplotlib stem | plotly static export | Plotly adds complexity for simple stick plots |
| RDKit atomNote | PIL/Pillow compositing | atomNote is built-in, no custom code needed |
| cairosvg for PNG | matplotlib savefig dpi | Unnecessary complexity, matplotlib handles this |

**Installation:**
```bash
pip install matplotlib
# or
uv add matplotlib
```

Note: RDKit is already installed (2025.09.3) with Cairo support. Pillow is already installed (12.1.0).

## Architecture Patterns

### Recommended Project Structure
```
src/qm_nmr_calc/
├── visualization.py      # NEW: All visualization code
├── tasks.py              # Existing: Call visualization after NMR calculation
├── storage.py            # Existing: Add helper for visualization file paths
└── api/routers/
    └── jobs.py           # Existing: Add endpoints for visualization downloads
```

### Pattern 1: Visualization Module
**What:** Single module containing all visualization functions
**When to use:** Always - keeps visualization logic isolated and testable
**Example:**
```python
# src/qm_nmr_calc/visualization.py

from pathlib import Path
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for server use

def generate_spectrum_plot(
    shifts: list[float],
    nucleus: str,  # "1H" or "13C"
    output_dir: Path,
) -> tuple[Path, Path]:
    """Generate stick spectrum plot in SVG and PNG formats.

    Returns:
        Tuple of (svg_path, png_path)
    """
    fig, ax = plt.subplots(figsize=(8, 4))

    # Stick plot: vertical lines at each shift, no markers, no baseline
    ax.stem(shifts, [1.0] * len(shifts), linefmt='k-', markerfmt=' ', basefmt=' ')

    # NMR convention: x-axis right to left (high ppm to low ppm)
    ax.invert_xaxis()

    # Minimal labeling
    ax.set_xlabel('ppm')
    ax.set_yticks([])  # No y-axis ticks (intensity is relative)

    # Add nucleus label in corner
    ax.text(0.02, 0.95, nucleus, transform=ax.transAxes,
            fontsize=12, verticalalignment='top')

    # Auto-fit with padding
    if shifts:
        ppm_range = max(shifts) - min(shifts)
        padding = max(ppm_range * 0.1, 0.5)  # At least 0.5 ppm padding
        ax.set_xlim(max(shifts) + padding, min(shifts) - padding)

    # Save both formats
    base_name = f"spectrum_{nucleus.replace('^', '')}"
    svg_path = output_dir / f"{base_name}.svg"
    png_path = output_dir / f"{base_name}.png"

    fig.savefig(svg_path, format='svg', bbox_inches='tight')
    fig.savefig(png_path, format='png', dpi=300, bbox_inches='tight')
    plt.close(fig)

    return svg_path, png_path


def generate_annotated_structure(
    smiles: str,
    h1_shifts: list[dict],  # [{"index": 1, "shift": 7.26}, ...]
    c13_shifts: list[dict],
    output_dir: Path,
) -> tuple[Path, Path]:
    """Generate structure with chemical shift annotations.

    Returns:
        Tuple of (svg_path, png_path)
    """
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    # Build index->shift mapping (NWChem uses 1-based indices)
    shift_map = {}
    for s in h1_shifts:
        shift_map[s['index'] - 1] = f"{s['shift']:.2f}"  # Convert to 0-based
    for s in c13_shifts:
        shift_map[s['index'] - 1] = f"{s['shift']:.2f}"

    # Set atomNote on each atom
    for idx, atom in enumerate(mol.GetAtoms()):
        if idx in shift_map:
            atom.SetProp('atomNote', shift_map[idx])

    # Prepare for drawing
    rdMolDraw2D.PrepareMolForDrawing(mol)

    # Generate SVG
    d_svg = rdMolDraw2D.MolDraw2DSVG(800, 600)
    d_svg.drawOptions().annotationFontScale = 0.6
    d_svg.DrawMolecule(mol)
    d_svg.FinishDrawing()

    svg_path = output_dir / "structure_annotated.svg"
    svg_path.write_text(d_svg.GetDrawingText())

    # Generate PNG (300 DPI equivalent: 800x600 * 3 = 2400x1800)
    d_png = rdMolDraw2D.MolDraw2DCairo(2400, 1800)
    d_png.drawOptions().annotationFontScale = 0.6
    d_png.DrawMolecule(mol)
    d_png.FinishDrawing()

    png_path = output_dir / "structure_annotated.png"
    d_png.WriteDrawingText(str(png_path))

    return svg_path, png_path
```

### Pattern 2: Integration with Task Flow
**What:** Call visualization after NMR calculation completes
**When to use:** In `tasks.py` after results are computed
**Example:**
```python
# In tasks.py, after nmr_results is computed:

from .visualization import generate_spectrum_plot, generate_annotated_structure

# Generate visualizations
h1_svg, h1_png = generate_spectrum_plot(
    shifts=[s.shift for s in h1_shifts],
    nucleus="1H",
    output_dir=output_dir,
)
c13_svg, c13_png = generate_spectrum_plot(
    shifts=[s.shift for s in c13_shifts],
    nucleus="13C",
    output_dir=output_dir,
)
struct_svg, struct_png = generate_annotated_structure(
    smiles=smiles,
    h1_shifts=[{"index": s.index, "shift": s.shift} for s in h1_shifts],
    c13_shifts=[{"index": s.index, "shift": s.shift} for s in c13_shifts],
    output_dir=output_dir,
)
```

### Pattern 3: API File Serving
**What:** Serve visualization files via dedicated endpoints
**When to use:** Follows existing pattern in jobs.py for geometry/output files
**Example:**
```python
# In api/routers/jobs.py

@router.get("/{job_id}/spectrum/1h.png")
async def download_h1_spectrum_png(job_id: str):
    """Download 1H spectrum plot as PNG."""
    # ... validation ...
    return FileResponse(
        path=get_visualization_file(job_id, "spectrum_1H.png"),
        media_type="image/png",
        filename=f"{job_id}_1H_spectrum.png",
    )

@router.get("/{job_id}/spectrum/1h.svg")
async def download_h1_spectrum_svg(job_id: str):
    """Download 1H spectrum plot as SVG."""
    # ... validation ...
    return FileResponse(
        path=get_visualization_file(job_id, "spectrum_1H.svg"),
        media_type="image/svg+xml",
        filename=f"{job_id}_1H_spectrum.svg",
    )
```

### Anti-Patterns to Avoid
- **Generating images on-the-fly at request time:** Generate during job completion, not during API request - expensive computation should happen in worker
- **Using interactive matplotlib backend:** Always use `matplotlib.use('Agg')` on server - prevents display errors
- **Complex SVG-to-PNG conversion pipelines:** Both matplotlib and RDKit export PNG directly - no need for cairosvg

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Stick spectrum plot | Custom SVG generation | matplotlib `stem()` | Handles axis scaling, labels, export |
| Molecule annotation | PIL text overlay | RDKit `atomNote` | Auto-positions labels, handles bonds |
| PNG at 300 DPI | cairosvg conversion | matplotlib `dpi=300` | Direct export, no dependencies |
| SVG export | Manual SVG strings | Library `.savefig()` | Proper escaping, metadata |
| Inverted x-axis | Manual coordinate flip | `ax.invert_xaxis()` | Preserves all other axis behavior |

**Key insight:** Both matplotlib and RDKit are designed for scientific visualization and have solved these problems already. Custom solutions will miss edge cases (font scaling, coordinate systems, file format specifics).

## Common Pitfalls

### Pitfall 1: Matplotlib Backend Issues on Server
**What goes wrong:** `RuntimeError: cannot find a display to connect to` when running matplotlib on headless server
**Why it happens:** matplotlib defaults to interactive backend requiring display
**How to avoid:** Set backend before importing pyplot
```python
import matplotlib
matplotlib.use('Agg')  # MUST be before importing pyplot
import matplotlib.pyplot as plt
```
**Warning signs:** Works locally but fails in production/CI

### Pitfall 2: RDKit Atom Index Mismatch
**What goes wrong:** Shift values appear on wrong atoms
**Why it happens:** NWChem uses 1-based indices, RDKit uses 0-based
**How to avoid:** Convert indices when mapping shifts to atoms
```python
rdkit_idx = nwchem_idx - 1  # Convert 1-based to 0-based
```
**Warning signs:** Hydrogen shifts appearing on carbons or vice versa

### Pitfall 3: atomNote Not Appearing
**What goes wrong:** Structure renders but annotations are missing
**Why it happens:** Either `PrepareMolForDrawing` called after setting notes, or font scale too small
**How to avoid:** Set `atomNote` BEFORE calling `PrepareMolForDrawing`, use appropriate `annotationFontScale`
```python
mol.GetAtomWithIdx(0).SetProp('atomNote', '7.26')
rdMolDraw2D.PrepareMolForDrawing(mol)  # AFTER setting notes
```
**Warning signs:** SVG has no `class='note'` path elements

### Pitfall 4: PNG Resolution Confusion
**What goes wrong:** PNG looks fine in code but blurry when printed
**Why it happens:** Confusion between pixel dimensions and DPI
**How to avoid:** For 300 DPI at 8x4 inches, need 2400x1200 pixels
```python
# matplotlib handles this automatically:
fig.savefig('plot.png', dpi=300)  # Correct

# RDKit needs explicit pixel dimensions:
# For 300 DPI equivalent at 800x600 base:
d = MolDraw2DCairo(2400, 1800)  # 3x scale = 300 DPI equivalent
```
**Warning signs:** File size much smaller than expected, image looks pixelated at 100%

### Pitfall 5: Memory Leak from Unclosed Figures
**What goes wrong:** Worker process memory grows over time
**Why it happens:** matplotlib figures not closed after saving
**How to avoid:** Always call `plt.close(fig)` after saving
```python
fig, ax = plt.subplots()
# ... draw ...
fig.savefig('output.png')
plt.close(fig)  # CRITICAL: release memory
```
**Warning signs:** Worker process memory increases with each job

## Code Examples

Verified patterns from official sources:

### Matplotlib Stem Plot with Hidden Markers/Baseline
```python
# Source: https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.stem.html
import matplotlib.pyplot as plt

x = [7.26, 5.12, 3.89, 2.01]  # Chemical shifts
y = [1.0] * len(x)  # All same height for stick spectrum

# markerfmt=' ' hides markers, basefmt=' ' hides baseline
plt.stem(x, y, linefmt='k-', markerfmt=' ', basefmt=' ')
plt.gca().invert_xaxis()  # NMR convention: high ppm on left
plt.xlabel('ppm')
plt.savefig('spectrum.png', dpi=300, bbox_inches='tight')
plt.close()
```

### Matplotlib Figure Export
```python
# Source: https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.savefig.html
fig, ax = plt.subplots(figsize=(8, 4))
# ... plotting code ...

# PNG at 300 DPI (print quality)
fig.savefig('output.png', dpi=300, bbox_inches='tight', pad_inches=0.1)

# SVG (vector, scalable)
fig.savefig('output.svg', format='svg', bbox_inches='tight')

plt.close(fig)
```

### RDKit Atom Annotation
```python
# Source: https://www.rdkit.org/docs/GettingStartedInPython.html
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

mol = Chem.MolFromSmiles('CCO')

# Set annotation BEFORE PrepareMolForDrawing
mol.GetAtomWithIdx(0).SetProp('atomNote', '14.2')
mol.GetAtomWithIdx(1).SetProp('atomNote', '58.3')

rdMolDraw2D.PrepareMolForDrawing(mol)

# SVG output
d = rdMolDraw2D.MolDraw2DSVG(400, 300)
d.drawOptions().annotationFontScale = 0.6
d.DrawMolecule(mol)
d.FinishDrawing()
svg = d.GetDrawingText()
```

### RDKit PNG Export with Cairo
```python
# Source: https://www.rdkit.org/docs/source/rdkit.Chem.Draw.rdMolDraw2D.html
from rdkit.Chem.Draw import rdMolDraw2D

# For 300 DPI equivalent, multiply desired inches by 300
# 8" x 6" at 300 DPI = 2400 x 1800 pixels
d = rdMolDraw2D.MolDraw2DCairo(2400, 1800)
d.drawOptions().annotationFontScale = 0.6
d.DrawMolecule(mol)
d.FinishDrawing()
d.WriteDrawingText('output.png')
```

### Inverted X-Axis (NMR Convention)
```python
# Source: https://matplotlib.org/stable/gallery/subplots_axes_and_figures/invert_axes.html
ax.invert_xaxis()  # High ppm on left, low ppm on right
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Separate annotation library | RDKit native `atomNote` | RDKit 2020.03 | No external dependencies needed |
| matplotlib GTK backend | AGG backend (headless) | Long-standing | Server-compatible rendering |
| cairosvg for PNG | Direct library export | N/A | Simpler pipeline |

**Deprecated/outdated:**
- Using `MolToImage()` for annotated structures (limited customization) - use `MolDraw2D` classes instead
- Interactive backends on server (Tk, Qt) - use 'Agg' backend

## Open Questions

Things that couldn't be fully resolved:

1. **Annotation collision handling for crowded structures**
   - What we know: RDKit positions annotations automatically with basic collision avoidance
   - What's unclear: Behavior with very crowded structures (many H atoms)
   - Recommendation: Test with a complex molecule (glucose, caffeine) during implementation; may need `annotationFontScale` adjustment

2. **Color differentiation for 1H vs 13C labels**
   - What we know: `atomNoteColour` sets color for all annotations
   - What's unclear: Whether per-atom colors are possible
   - Recommendation: Start with single color; if distinction needed, could use two separate images or post-process

## Sources

### Primary (HIGH confidence)
- [matplotlib.pyplot.savefig](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.savefig.html) - DPI, format options
- [matplotlib.pyplot.stem](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.stem.html) - Stick plot parameters
- [matplotlib invert_axes](https://matplotlib.org/stable/gallery/subplots_axes_and_figures/invert_axes.html) - X-axis inversion
- [RDKit rdMolDraw2D](https://www.rdkit.org/docs/source/rdkit.Chem.Draw.rdMolDraw2D.html) - MolDraw2DSVG, MolDraw2DCairo
- [RDKit Drawing on drawings blog](https://greglandrum.github.io/rdkit-blog/posts/2025-03-07-drawing-on-drawings.html) - Annotation examples
- Local verification: RDKit 2025.09.3 atomNote functionality tested and confirmed working

### Secondary (MEDIUM confidence)
- [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html) - General drawing patterns
- [RDKit Getting Started](https://www.rdkit.org/docs/GettingStartedInPython.html) - atomNote property usage

### Tertiary (LOW confidence)
- None - all critical claims verified with official documentation

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - matplotlib and RDKit are well-documented standard tools
- Architecture: HIGH - follows existing codebase patterns (storage.py, tasks.py, jobs.py)
- Pitfalls: HIGH - verified through local testing (atomNote, matplotlib backend)

**Research date:** 2026-01-20
**Valid until:** 60 days (stable libraries, well-documented APIs)
