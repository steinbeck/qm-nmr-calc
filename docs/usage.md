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

## REST API Reference

The REST API provides programmatic access for automation and integration with other tools.

### Overview

| Property | Value |
|----------|-------|
| Base URL | `http://localhost:8000/api/v1` |
| Authentication | None required (local deployment) |
| Content-Type | `application/json` (requests and responses) |
| Workflow | Asynchronous: submit returns 202, poll for completion |

### Health Endpoints

**Liveness Probe:**

```bash
curl http://localhost:8000/health
```

Response:
```json
{"status": "alive"}
```

**Readiness Probe:**

```bash
curl http://localhost:8000/health/ready
```

Response:
```json
{
  "status": "ready",
  "checks": {
    "data_directory": "ok",
    "task_queue": "ok",
    "crest_available": true
  },
  "crest_available": true,
  "timestamp": "2026-02-01T12:00:00Z"
}
```

If any component is unavailable, returns HTTP 503 with `"status": "not ready"`.

### Job Submission

**Submit via SMILES (POST /api/v1/jobs):**

```bash
curl -X POST http://localhost:8000/api/v1/jobs \
  -H "Content-Type: application/json" \
  -d '{
    "smiles": "CCO",
    "solvent": "chcl3",
    "preset": "production",
    "name": "Ethanol",
    "conformer_mode": "ensemble",
    "conformer_method": "rdkit_kdg",
    "max_conformers": null,
    "notification_email": null
  }'
```

Request body fields:

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `smiles` | string | Yes | SMILES representation of molecule |
| `solvent` | string | Yes | NMR solvent (see [README](../README.md#supported-solvents) for all 13 codes, e.g., `chcl3`, `dmso`, `acetonitrile`) |
| `preset` | string | No | `draft` or `production` (default: `production`) |
| `name` | string | No | Optional molecule label (max 100 chars) |
| `conformer_mode` | string | No | `single` or `ensemble` (default: `single`) |
| `conformer_method` | string | No | `rdkit_kdg` or `crest` (ensemble only) |
| `max_conformers` | int | No | Override automatic conformer count |
| `notification_email` | string | No | Email for completion notification |

Response (HTTP 202 Accepted):
```json
{
  "job_id": "a1b2c3d4e5f6",
  "status": "queued",
  "created_at": "2026-02-01T12:00:00Z",
  "input_smiles": "CCO",
  "input_name": "Ethanol",
  "preset": "production",
  "solvent": "chcl3",
  "conformer_mode": "ensemble"
}
```

Headers:
- `Location: /api/v1/jobs/a1b2c3d4e5f6` - URL for status polling
- `Retry-After: 30` - Suggested polling interval

**Submit via File Upload (POST /api/v1/jobs/upload):**

```bash
curl -X POST http://localhost:8000/api/v1/jobs/upload \
  -F "file=@molecule.mol" \
  -F "solvent=chcl3" \
  -F "preset=production" \
  -F "conformer_mode=single"
```

Accepts `.mol` and `.sdf` files with a single molecule. Returns same response format as SMILES submission.

**List Available Solvents (GET /api/v1/jobs/solvents):**

```bash
curl http://localhost:8000/api/v1/jobs/solvents
```

Response:
```json
["acetone", "acetonitrile", "benzene", "chcl3", "dcm", "dmf", "dmso", "methanol", "pyridine", "thf", "toluene", "vacuum", "water"]
```

### Job Status

**Get Job Status (GET /api/v1/jobs/{job_id}):**

```bash
curl http://localhost:8000/api/v1/jobs/a1b2c3d4e5f6
```

Response varies by job status:

**Queued job:**
```json
{
  "job_id": "a1b2c3d4e5f6",
  "status": "queued",
  "created_at": "2026-02-01T12:00:00Z",
  "started_at": null,
  "current_step": null,
  "steps_completed": []
}
```

**Running job (single conformer):**
```json
{
  "job_id": "a1b2c3d4e5f6",
  "status": "running",
  "created_at": "2026-02-01T12:00:00Z",
  "started_at": "2026-02-01T12:00:05Z",
  "current_step": "geometry_optimization",
  "current_step_label": "Geometry Optimization",
  "step_started_at": "2026-02-01T12:00:05Z",
  "steps_completed": [],
  "conformer_mode": "single"
}
```

**Running job (ensemble with progress):**
```json
{
  "job_id": "a1b2c3d4e5f6",
  "status": "running",
  "current_step": "nmr_calculation",
  "current_step_label": "NMR Calculation",
  "conformer_mode": "ensemble",
  "conformer_method": "rdkit_kdg",
  "conformer_count": 5,
  "conformer_progress": [
    {"conformer_id": "conf_001", "status": "nmr_complete", "energy_kcal": 0.0, "population": 0.45},
    {"conformer_id": "conf_002", "status": "nmr_running", "energy_kcal": 0.82, "population": null},
    {"conformer_id": "conf_003", "status": "optimized", "energy_kcal": 1.15, "population": null}
  ]
}
```

**Complete job:**
```json
{
  "job_id": "a1b2c3d4e5f6",
  "status": "complete",
  "completed_at": "2026-02-01T12:15:00Z",
  "nmr_results": {
    "h1_shifts": [
      {"index": 1, "atom": "H", "shift": 3.65},
      {"index": 2, "atom": "H", "shift": 1.18}
    ],
    "c13_shifts": [
      {"index": 1, "atom": "C", "shift": 57.8},
      {"index": 2, "atom": "C", "shift": 18.2}
    ],
    "functional": "b3lyp",
    "basis_set": "6-311+G(2d,p)",
    "solvent": "chcl3",
    "scaling_factor_source": "DELTA50",
    "h1_expected_mae": "+/- 0.12 ppm",
    "c13_expected_mae": "+/- 1.95 ppm"
  }
}
```

**Failed job:**
```json
{
  "job_id": "a1b2c3d4e5f6",
  "status": "failed",
  "error_message": "NWChem geometry optimization failed: SCF did not converge"
}
```

**Polling pattern:**

```bash
# Poll every 30 seconds until complete
while true; do
  STATUS=$(curl -s http://localhost:8000/api/v1/jobs/a1b2c3d4e5f6 | jq -r '.status')
  echo "Status: $STATUS"
  if [ "$STATUS" = "complete" ] || [ "$STATUS" = "failed" ]; then
    break
  fi
  sleep 30
done
```

### Results Retrieval

All result endpoints require `status=complete`. Returns HTTP 409 if job is still running.

**NMR Shifts (GET /api/v1/jobs/{job_id}/results):**

```bash
curl http://localhost:8000/api/v1/jobs/a1b2c3d4e5f6/results
```

Response:
```json
{
  "h1_shifts": [
    {"index": 1, "atom": "H", "shift": 3.65},
    {"index": 2, "atom": "H", "shift": 1.18}
  ],
  "c13_shifts": [
    {"index": 1, "atom": "C", "shift": 57.8},
    {"index": 2, "atom": "C", "shift": 18.2}
  ],
  "functional": "b3lyp",
  "basis_set": "6-311+G(2d,p)",
  "solvent": "chcl3",
  "scaling_factor_source": "DELTA50",
  "h1_expected_mae": "+/- 0.12 ppm",
  "c13_expected_mae": "+/- 1.95 ppm",
  "ensemble_metadata": {
    "conformer_count": 5,
    "total_generated": 12,
    "method": "rdkit_kdg",
    "temperature_k": 298.15,
    "energy_range_kcal": 2.34,
    "top_populations": [
      {"id": "conf_001", "population": 0.45, "energy_kcal": 0.0},
      {"id": "conf_003", "population": 0.28, "energy_kcal": 0.52}
    ]
  }
}
```

**Optimized Geometry (XYZ format):**

```bash
curl http://localhost:8000/api/v1/jobs/a1b2c3d4e5f6/geometry \
  -o ethanol_optimized.xyz
```

**Optimized Geometry (SDF format with bonds):**

```bash
curl http://localhost:8000/api/v1/jobs/a1b2c3d4e5f6/geometry.sdf \
  -o ethanol_optimized.sdf
```

**3D Viewer Data (JSON with SDF + shift assignments):**

```bash
curl http://localhost:8000/api/v1/jobs/a1b2c3d4e5f6/geometry.json
```

Response includes XYZ, SDF, and atom-to-shift mappings for 3D visualization.

**Raw NWChem Output (ZIP archive):**

```bash
curl http://localhost:8000/api/v1/jobs/a1b2c3d4e5f6/output \
  -o nwchem_output.zip
```

**Spectrum Images:**

```bash
# 1H NMR spectrum
curl http://localhost:8000/api/v1/jobs/a1b2c3d4e5f6/spectrum/1h.png -o spectrum_1h.png
curl http://localhost:8000/api/v1/jobs/a1b2c3d4e5f6/spectrum/1h.svg -o spectrum_1h.svg

# 13C NMR spectrum
curl http://localhost:8000/api/v1/jobs/a1b2c3d4e5f6/spectrum/13c.png -o spectrum_13c.png
curl http://localhost:8000/api/v1/jobs/a1b2c3d4e5f6/spectrum/13c.svg -o spectrum_13c.svg
```

**Annotated Structure Images:**

```bash
curl http://localhost:8000/api/v1/jobs/a1b2c3d4e5f6/structure.png -o structure.png
curl http://localhost:8000/api/v1/jobs/a1b2c3d4e5f6/structure.svg -o structure.svg
```

### Error Handling

**HTTP Status Codes:**

| Code | Meaning | Example |
|------|---------|---------|
| 200 | Success | GET /jobs/{id}/results for complete job |
| 202 | Accepted | POST /jobs (job queued) |
| 404 | Not Found | Job ID doesn't exist |
| 409 | Conflict | Requesting results for incomplete job |
| 422 | Validation Error | Invalid SMILES or solvent |

**Error Response Format (RFC 7807 Problem Details):**

```json
{
  "type": "https://qm-nmr-calc.example/problems/invalid-smiles",
  "title": "Invalid SMILES String",
  "status": 422,
  "detail": "RDKit failed to parse SMILES: 'invalid-smiles'"
}
```

Common error types:
- `invalid-smiles` - SMILES string could not be parsed
- `invalid-solvent` - Unknown solvent name
- `invalid-file-type` - Uploaded file is not .mol or .sdf
- `job-not-found` - No job exists with that ID
- `job-not-complete` - Results requested before job finished

### Complete Workflow Example

This bash script submits a job, polls for completion, and downloads all results:

```bash
#!/bin/bash
# nmr_calculate.sh - Submit NMR calculation and download results
# Usage: ./nmr_calculate.sh "CCO" "ethanol"

SMILES="$1"
NAME="${2:-molecule}"
BASE_URL="http://localhost:8000"
TIMEOUT=3600  # 1 hour max

# Submit job
echo "Submitting $NAME ($SMILES)..."
RESPONSE=$(curl -s -X POST "$BASE_URL/api/v1/jobs" \
  -H "Content-Type: application/json" \
  -d "{\"smiles\": \"$SMILES\", \"solvent\": \"chcl3\", \"preset\": \"production\", \"name\": \"$NAME\"}")

JOB_ID=$(echo "$RESPONSE" | jq -r '.job_id')
if [ "$JOB_ID" = "null" ]; then
  echo "Error: $(echo "$RESPONSE" | jq -r '.detail')"
  exit 1
fi
echo "Job ID: $JOB_ID"

# Poll for completion
START_TIME=$(date +%s)
while true; do
  STATUS_RESPONSE=$(curl -s "$BASE_URL/api/v1/jobs/$JOB_ID")
  STATUS=$(echo "$STATUS_RESPONSE" | jq -r '.status')
  STEP=$(echo "$STATUS_RESPONSE" | jq -r '.current_step_label // "Waiting"')

  echo "Status: $STATUS - $STEP"

  if [ "$STATUS" = "complete" ]; then
    echo "Calculation complete!"
    break
  elif [ "$STATUS" = "failed" ]; then
    ERROR=$(echo "$STATUS_RESPONSE" | jq -r '.error_message')
    echo "Job failed: $ERROR"
    exit 1
  fi

  # Check timeout
  ELAPSED=$(($(date +%s) - START_TIME))
  if [ $ELAPSED -gt $TIMEOUT ]; then
    echo "Timeout after $TIMEOUT seconds"
    exit 1
  fi

  sleep 30
done

# Download results
mkdir -p "results_$NAME"
cd "results_$NAME"

echo "Downloading results..."
curl -s "$BASE_URL/api/v1/jobs/$JOB_ID/results" > nmr_shifts.json
curl -s "$BASE_URL/api/v1/jobs/$JOB_ID/geometry" -o optimized.xyz
curl -s "$BASE_URL/api/v1/jobs/$JOB_ID/geometry.sdf" -o optimized.sdf
curl -s "$BASE_URL/api/v1/jobs/$JOB_ID/spectrum/1h.png" -o spectrum_1h.png
curl -s "$BASE_URL/api/v1/jobs/$JOB_ID/spectrum/13c.png" -o spectrum_13c.png
curl -s "$BASE_URL/api/v1/jobs/$JOB_ID/output" -o nwchem_output.zip

# Display shifts
echo ""
echo "1H Chemical Shifts:"
jq -r '.h1_shifts[] | "  H\(.index): \(.shift) ppm"' nmr_shifts.json

echo ""
echo "13C Chemical Shifts:"
jq -r '.c13_shifts[] | "  C\(.index): \(.shift) ppm"' nmr_shifts.json

echo ""
echo "Results saved to results_$NAME/"
```

## Result Interpretation

Understanding how to interpret NMR calculation results is essential for making effective use of predicted chemical shifts.

### Chemical Shift Tables

The results page displays separate tables for 1H and 13C chemical shifts:

**1H (Proton) Shifts:**
- Typical range: 0-12 ppm for most organic molecules
- Downfield shifts (higher ppm): deshielded protons near electronegative atoms
- Upfield shifts (lower ppm): shielded protons in alkyl groups
- TMS reference: all shifts are relative to tetramethylsilane (0.0 ppm)

**13C (Carbon) Shifts:**
- Typical range: 0-220 ppm
- Carbonyl carbons: 160-220 ppm
- Aromatic/alkene carbons: 100-160 ppm
- Alkyl carbons: 0-60 ppm
- TMS reference: relative to tetramethylsilane (0.0 ppm)

**Atom Indexing:**
- Indices (H1, H2, C1, C2...) correspond to the canonical SMILES atom order
- For complex molecules, compare with the 3D viewer where shift labels are overlaid on atoms

### Expected Accuracy

The calculator provides Mean Absolute Error (MAE) estimates based on the DELTA50 benchmark dataset:

| Preset | 1H MAE | 13C MAE |
|--------|--------|---------|
| Production | ~0.12-0.15 ppm | ~1.9-2.5 ppm |
| Draft | ~0.3-0.5 ppm | ~4-6 ppm |

**Factors affecting accuracy:**
- **Molecule class:** MAE values derived from typical organic molecules; unusual functional groups may differ
- **Conformational flexibility:** Ensemble mode improves accuracy for flexible molecules
- **Solvent matching:** Use the same solvent as your experimental NMR
- **Scaling factors:** Applied automatically using DELTA50 regression parameters

**Interpreting MAE:**
- The displayed `+/- X.XX ppm` indicates typical deviation from experimental values
- Individual atom shifts may deviate more or less than the average
- Relative shift ordering is often more reliable than absolute values

### Spectrum Visualizations

The simulated spectra provide a visual representation of predicted NMR patterns:

**Spectrum Features:**
- X-axis: Chemical shift in ppm (reversed: high to low, left to right)
- Y-axis: Relative intensity (arbitrary units)
- Peak shape: Lorentzian line broadening centered at each predicted shift
- 1H spectrum: Typically displayed 0-12 ppm range
- 13C spectrum: Typically displayed 0-220 ppm range

**Using Spectra:**
- Compare overall pattern with experimental spectrum
- Peak positions should align within expected MAE
- Relative spacing between peaks is often more reliable than absolute positions
- Click spectrum images in web UI to enlarge

**Download Options:**
- PNG (300 DPI): Good for presentations and documents
- SVG: Vector format, ideal for publications and scaling

### 3D Structure Viewer

The interactive 3D viewer displays the optimized molecular geometry with shift labels:

**Viewer Controls:**
- **Rotate:** Click and drag
- **Zoom:** Scroll wheel
- **Reset view:** Double-click

**Shift Labels:**
- Blue labels: 1H chemical shifts in ppm
- Orange labels: 13C chemical shifts in ppm
- Labels positioned near corresponding atoms
- Ball-and-stick model with Jmol-style element colors

**For Ensemble Calculations:**
- Conformer dropdown: Select different conformer geometries to view
- Shift labels: Always show Boltzmann-averaged values (not per-conformer shifts)
- Default view: Lowest-energy conformer geometry
- Population percentages: Shown in conformer selector

### Ensemble Results

For ensemble calculations, shifts represent population-weighted averages across conformers:

**Boltzmann Averaging:**
- Each conformer contributes proportionally to its Boltzmann population
- Population = exp(-E/RT) / sum(exp(-E/RT)) where E is relative energy
- Temperature: 298.15 K (room temperature)
- Higher population conformers have larger influence on final shifts

**Ensemble Metadata Panel:**
- **Conformers used:** Number of conformers included in averaging
- **Total generated:** Conformers created before energy filtering
- **Energy range:** Spread of DFT energies in kcal/mol
- **Top contributors:** Highest-population conformers with their percentages
- **Conformer method:** RDKit KDG or CREST

**When to use ensemble:**
- Molecules with 3+ rotatable bonds
- When single-conformer result seems inconsistent with experiment
- For publication-quality predictions of flexible molecules

### Downloading Results

**Recommended file formats for different purposes:**

| Purpose | Format | Why |
|---------|--------|-----|
| Further calculations | XYZ | Universal, simple coordinates |
| Visualization tools | SDF | Contains bonds, opens in ChemDraw/Avogadro |
| Publications | SVG | Vector graphics, scales without pixelation |
| Presentations | PNG | Raster format, widely compatible |
| Reproducibility | ZIP | Full NWChem output for debugging |

**Structure files:**
- **XYZ:** Cartesian coordinates only (element + x, y, z)
- **SDF:** Includes bond connectivity from original SMILES

**For ensemble jobs:**
- Geometry files contain the lowest-energy conformer
- Full conformer ensemble available via geometry.json endpoint

### Troubleshooting Results

**Unexpected chemical shifts:**
1. Check solvent matches your experimental conditions
2. Try ensemble mode for flexible molecules
3. Verify SMILES correctly represents your structure
4. Consider that some functional groups have higher errors

**Missing atoms in 3D viewer:**
- Ensure SMILES includes all atoms (explicit hydrogens if needed)
- Check for disconnected fragments in SMILES

**Large MAE values:**
- Normal for unusual functional groups (phosphorus, sulfur, halogens)
- Consider that benchmark MAE was derived from specific molecule types
- Use relative shift ordering rather than absolute values

**Calculation failures:**
- Very large molecules may exceed memory limits
- Unusual bonding patterns may fail geometry optimization
- Check NWChem output (ZIP download) for specific error messages

For detailed methodology and scientific references, see the [Science documentation](science.md).

## Related Documentation

- [Installation Guide](installation.md) - System setup and dependencies
- [Architecture](architecture.md) - System design and implementation details
- [Science](science.md) - NMR methodology, DP4+ probability, and references

