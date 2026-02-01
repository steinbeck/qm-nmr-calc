# Usage Guide

Complete usage documentation for the QM NMR Calculator, covering the web interface and calculation concepts.

**Target audience:** Researchers running NMR calculations and interpreting results.

**Prerequisites:** Complete the [Installation Guide](installation.md) before using this guide.

## Introduction

The QM NMR Calculator provides quantum mechanical prediction of 1H and 13C NMR chemical shifts using density functional theory (DFT). You can interact with the calculator through:

- **Web UI** - Visual interface at `http://localhost:8000` for molecule submission, job monitoring, and result visualization
- **REST API** - Programmatic access for automation and integration (documented separately)

This guide covers:
1. Web UI workflow (submit, status, results pages)
2. Calculation modes (single-conformer vs ensemble)
3. Solvent selection for COSMO solvation
4. Preset options (draft vs production quality)

## Web UI Workflow

The web interface consists of three main pages that guide you through the NMR calculation process.

### Submit Page (`/`)

The home page is where you submit molecules for NMR calculation.

**Molecule Input Options:**

You can provide molecular structure in two ways (mutually exclusive):

1. **SMILES String** - Enter the SMILES representation of your molecule
   - Example: `CCO` for ethanol, `c1ccccc1` for benzene
   - Live preview shows your molecule as you type (powered by SmilesDrawer)
   - The preview validates SMILES syntax before submission

2. **MOL/SDF File Upload** - Upload a pre-built structure file
   - Accepts `.mol` and `.sdf` formats
   - Single molecule per file
   - Useful when you have 3D coordinates from other software

**Required Settings:**

| Setting | Description |
|---------|-------------|
| **Solvent** | NMR solvent environment (see [Solvent Selection](#solvent-selection)) |
| **Calculation Mode** | Single conformer or ensemble (see [Calculation Modes](#calculation-modes)) |

**Optional Settings:**

| Setting | Description |
|---------|-------------|
| **Preset** | Draft (fast) or Production (accurate) - defaults to Production |
| **Conformer Method** | RDKit KDG (always available) or CREST/xTB (if installed) |
| **Max Conformers** | Override automatic conformer count (leave blank for auto) |
| **Molecule Name** | Label for your molecule (displayed in results) |
| **Email Notification** | Get notified when calculation completes |

**Submitting:**

Click "Calculate NMR" to start the calculation. You'll be redirected to the status page to monitor progress.

> **Note:** Invalid SMILES or unsupported file formats will show an error message. The molecule preview helps catch SMILES errors before submission.

### Status Page (`/status/{job_id}`)

After submission, you're redirected to the status page to monitor calculation progress.

**Job Details Panel:**

Shows your calculation settings:
- Job ID (unique identifier)
- Elapsed time since submission
- Input SMILES
- Solvent and preset selections

**3D Molecule Preview:**

An interactive 3D viewer shows the initial RDKit-generated geometry while the calculation runs. This is the starting structure before DFT optimization.

**Step Progress Tracker:**

Shows completion status for each calculation step:

For **single conformer** mode:
1. Geometry Optimization - DFT structure optimization
2. NMR Calculation - Shielding tensor computation
3. Post-processing - Spectrum generation and formatting

For **ensemble** mode:
1. Conformer Generation - Create multiple conformers
2. Geometry Optimization - DFT optimize each conformer
3. NMR Calculation - Compute shifts for each conformer
4. Averaging Shifts - Boltzmann-weighted averaging
5. Post-processing - Final spectrum generation

**Conformer Progress (Ensemble Only):**

For ensemble calculations, an additional panel shows:
- Current processing stage (optimizing, NMR running)
- Conformer count (X/N complete)
- Progress bar showing overall completion
- Estimated time remaining

**Auto-refresh Behavior:**

The page automatically polls for updates every 3 seconds and redirects to the results page upon completion. If the job fails, an error message is displayed with the option to submit a new job.

### Results Page (`/results/{job_id}`)

The results page provides comprehensive visualization and download options for completed calculations.

**Interactive 3D Viewer:**

The main feature is a rotatable 3D molecular structure with chemical shift labels:
- **Blue labels** - 1H (proton) chemical shifts in ppm
- **Orange labels** - 13C (carbon) chemical shifts in ppm
- Click and drag to rotate, scroll to zoom

For ensemble calculations, a dropdown selector lets you view different conformer geometries. The shift labels always show Boltzmann-averaged values regardless of which conformer geometry is displayed.

**Spectrum Images:**

- **1H NMR Spectrum** - Simulated proton spectrum (click to enlarge)
- **13C NMR Spectrum** - Simulated carbon spectrum (click to enlarge)

Both images can be downloaded as PNG files.

**Ensemble Metadata (Ensemble Only):**

For ensemble calculations, a summary panel shows:
- Conformers used in averaging
- Total conformers generated
- Conformer method (RDKit KDG or CREST/xTB)
- Temperature used for Boltzmann weighting (298 K)
- Energy range across conformers (kcal/mol)
- Top contributors table (highest-population conformers)

**Chemical Shift Tables:**

Separate tables for 1H and 13C shifts:
- Atom index (H1, H2... or C1, C2...)
- Chemical shift value in ppm

**Downloads Section:**

| Download | Description |
|----------|-------------|
| Geometry (XYZ) | Optimized 3D coordinates in XYZ format |
| Geometry (SDF) | Optimized structure with connectivity in SDF format |
| Raw Output (ZIP) | Complete NWChem output files for detailed analysis |
| 1H Spectrum (PNG) | Simulated 1H NMR spectrum image |
| 13C Spectrum (PNG) | Simulated 13C NMR spectrum image |
| Structure (PNG) | 2D structure image |

## Calculation Modes

The calculator supports two conformational approaches, selected on the submit page.

### Single Conformer Mode

Uses a single low-energy conformer for the entire calculation.

**How it works:**
1. RDKit generates an initial 3D structure
2. DFT optimizes this single geometry
3. NMR shieldings computed on the optimized structure
4. Shifts converted using TMS reference

**Best for:**
- Rigid molecules with few rotatable bonds
- Quick screening or validation of SMILES
- Large molecules where ensemble is too slow
- Initial exploration before production runs

**Typical timing:**
- Draft preset: 3-10 minutes
- Production preset: 15-45 minutes

**Limitations:**
- May miss conformer-dependent shift variations
- Less accurate for flexible molecules
- Uses arbitrary conformer if multiple minima exist

### Ensemble Mode (Recommended)

Generates multiple conformers and averages their NMR shifts using Boltzmann weighting.

**How it works:**
1. Generate conformer ensemble (RDKit or CREST)
2. Pre-select low-energy conformers using GFN2-xTB
3. DFT optimize each selected conformer
4. Compute NMR shifts for each optimized geometry
5. Calculate Boltzmann populations from DFT energies
6. Weight-average all shifts by population

**Best for:**
- Flexible molecules with rotatable bonds
- Molecules with multiple low-energy conformers
- Publication-quality predictions
- Accurate comparison with experimental spectra

**Conformer Methods:**

| Method | Speed | Thoroughness | Availability |
|--------|-------|--------------|--------------|
| RDKit KDG | Fast (seconds) | Good for most molecules | Always available |
| CREST/xTB | Slower (minutes) | More thorough search | Requires installation |

RDKit uses Knowledge-based Distance Geometry (KDG) for conformer generation. CREST performs metadynamics-based conformational search using the GFN2-xTB force field.

**When to choose each method:**
- **RDKit KDG**: Default choice, good balance of speed and coverage
- **CREST/xTB**: Highly flexible molecules, when RDKit misses important conformers

**Typical timing (ensemble):**
- Draft preset: 10-30 minutes (depends on conformer count)
- Production preset: 30-120 minutes

> **Note:** The conformer count is automatically determined based on molecular flexibility. More rotatable bonds = more conformers generated.

### Decision Guide

Use this flowchart to choose the right mode:

```
Is the molecule rigid (0-2 rotatable bonds)?
├── Yes → Single Conformer (faster, sufficient accuracy)
└── No → Is this a quick screening run?
    ├── Yes → Single Conformer (for speed)
    └── No → Ensemble Mode (for accuracy)
              └── Is CREST installed and available?
                  ├── Yes → Consider CREST for very flexible molecules
                  └── No → RDKit KDG is excellent for most cases
```

## Solvent Selection

Solvent environment affects both geometry optimization and NMR shielding calculations through the COSMO solvation model.

### Available Solvents

| Solvent | Display Name | Dielectric Constant | Common Use |
|---------|--------------|---------------------|------------|
| `chcl3` | CDCl3 | 4.81 | Most common NMR solvent |
| `dmso` | DMSO-d6 | 46.7 | Polar/ionic molecules, biomolecules |
| `vacuum` | Gas phase | 1.0 | Reference calculations, gas-phase comparisons |

### COSMO Solvation Model

The calculator uses the Conductor-like Screening Model (COSMO) to simulate solvent effects:

**What COSMO does:**
- Creates a cavity around the solute molecule
- Models the solvent as a dielectric continuum
- Calculates electrostatic screening from induced surface charges
- Affects both molecular geometry and electronic properties

**Impact on results:**
- Polar groups experience different shielding in polar vs non-polar solvents
- Geometry can change slightly between solvents
- Hydrogen bonding effects are approximated (not explicit)

### Choosing the Right Solvent

**Match your experimental conditions:**

The most important rule is to select the solvent that matches your experimental NMR spectrum. Solvent effects on chemical shifts can be significant (0.1-1.0 ppm for 1H).

**Chloroform (CDCl3):**
- Most common choice for organic molecules
- Good for non-polar to moderately polar compounds
- Use when your experimental spectrum was recorded in CDCl3

**DMSO (DMSO-d6):**
- Better for polar molecules with hydrogen bond donors/acceptors
- Common for pharmaceuticals, natural products
- Larger dielectric constant captures polar interactions

**Gas phase (vacuum):**
- Reference calculations without solvent effects
- Comparing with gas-phase experimental data
- Understanding intrinsic vs solvent-induced shifts

> **Note:** If you don't know which solvent to use, CDCl3 is a safe default for most organic molecules.

## Calculation Presets

Presets control the level of theory used for DFT calculations, trading accuracy for speed.

### Draft Preset

**Configuration:**
- Functional: B3LYP
- Basis set (geometry): 6-31G*
- Basis set (NMR): 6-31G*

**Characteristics:**
- Faster calculations (3-10 min for small molecules)
- Lower accuracy (typical MAE ~0.3-0.5 ppm for 1H)
- Good for initial screening

**Use when:**
- Validating SMILES before production run
- Quick conformational insights
- Testing calculation setup
- Large molecules where production is too slow

### Production Preset (Default)

**Configuration:**
- Functional: B3LYP
- Basis set (geometry): 6-31G*
- Basis set (NMR): 6-311+G(2d,p)

**Characteristics:**
- Higher accuracy (typical MAE ~0.15-0.25 ppm for 1H)
- Larger basis set for NMR captures electron distribution better
- Diffuse functions (+) improve description of electron density tails
- Polarization functions (2d,p) better for anisotropic shielding

**Use when:**
- Publication-quality predictions
- Structure verification against experimental spectra
- Comparing closely similar candidate structures
- Final results after draft screening

### Accuracy Comparison

| Metric | Draft | Production |
|--------|-------|------------|
| 1H MAE (ppm) | ~0.3-0.5 | ~0.15-0.25 |
| 13C MAE (ppm) | ~3-5 | ~2-3 |
| Calculation time | 1x | 3-5x |
| Memory usage | Lower | Higher |

*MAE = Mean Absolute Error compared to experimental values*

> **Note:** Both presets use the B3LYP functional with the GIAO (Gauge-Including Atomic Orbital) method for NMR shielding calculations. The difference is in basis set size for the NMR step.

## Tips for Best Results

### Before Submission

1. **Validate your SMILES** - Use the live preview to catch syntax errors
2. **Choose appropriate mode** - Ensemble for flexible molecules, single for rigid
3. **Match experimental solvent** - Use the same solvent as your NMR experiment

### Interpreting Results

1. **Check conformer populations** - If one conformer dominates (>90%), single conformer mode would give similar results
2. **Compare spectrum patterns** - Overall shape and relative positions matter more than absolute values
3. **Account for systematic errors** - DFT shifts may be uniformly offset from experiment

### Performance Tips

1. **Start with draft** - Verify your setup before running expensive production calculations
2. **Use single conformer for large molecules** - 50+ atoms benefit from faster single conformer mode
3. **Monitor job status** - Long-running jobs can be checked via the status page

## Related Documentation

- [Installation Guide](installation.md) - System setup and dependencies
- [Architecture](architecture.md) - System design and implementation details
- [Science](science.md) - NMR methodology, DP4+ probability, and references

