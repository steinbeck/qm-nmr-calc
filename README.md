# QM NMR Calculator

Asynchronous web service for predicting NMR chemical shifts using quantum mechanical calculations. Submit a molecule (SMILES or structure file), and get back predicted 1H and 13C chemical shifts with spectrum visualizations.

## Features

- **Molecule Input**: Submit via SMILES string, MOL file, or SDF file
- **Async Processing**: Long-running QM calculations run in background
- **REST API**: Programmatic access to all functionality
- **Web UI**: Browser-based submission and results viewing
- **Visualizations**: Spectrum plots and annotated structure drawings
- **Email Notifications**: Optional notification when calculations complete
- **Calculation Presets**: Draft (fast) and Production (accurate) modes

## Architecture

```
┌─────────────────┐     ┌─────────────────┐     ┌─────────────────┐
│   Web UI /      │────>│   FastAPI       │────>│   Huey Queue    │
│   REST API      │<────│   Server        │<────│   (SQLite)      │
└─────────────────┘     └─────────────────┘     └────────┬────────┘
                                                         │
                                                         v
                                                ┌─────────────────┐
                                                │  Huey Consumer  │
                                                │  (1 worker)     │
                                                └────────┬────────┘
                                                         │
                                                         v
                                                ┌─────────────────┐
                                                │  RDKit/NWChem  │
                                                │  DFT Calcs      │
                                                └─────────────────┘
```

## Prerequisites

### System Requirements

- **OS**: Linux (tested on Debian/Ubuntu)
- **Python**: 3.11+
- **RAM**: 8GB minimum, 16GB+ recommended
- **Disk**: 10GB+ for NWChem scratch files
- **CPU**: Multi-core recommended (NWChem uses MPI)

### Required Software

1. **NWChem 7.0.2+** - Quantum chemistry package
2. **Open MPI** - For parallel NWChem execution
3. **uv** - Python package manager (recommended)

## Installation

### 1. Install System Dependencies

```bash
# Ubuntu/Debian
sudo apt-get update
sudo apt-get install -y \
    nwchem \
    openmpi-bin \
    libopenmpi-dev \
    python3.11 \
    python3.11-venv \
    git

# Verify NWChem is available
which nwchem
nwchem --version  # May not work, but binary should exist
```

### 2. Install uv (Python Package Manager)

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
source ~/.bashrc  # or restart shell
```

### 3. Clone and Setup QM NMR Calculator

```bash
cd ~/develop
git clone https://github.com/steinbeck/qm-nmr-calc.git
cd qm-nmr-calc

# Install dependencies
uv sync

# Verify installation
uv run python -c "import qm_nmr_calc; print('OK')"
```

### 4. Verify Environment

```bash
# Quick validation
uv run python -c "
from qm_nmr_calc.startup import validate_environment
validate_environment()
"
```

Expected output:
```
Validating environment...
  [OK] NWChem found
  [OK] RDKit working
  [OK] Data directory writable: data/jobs
Environment validation passed
```

## Running the Service

You need to run **two processes**: the API server and the Huey consumer.

### Terminal 1: Start the Huey Consumer (Job Processor)

```bash
cd ~/develop/qm-nmr-calc
uv run python scripts/run_consumer.py
```

This will:
1. Validate the environment (NWChem, RDKit, directories)
2. Recover any interrupted jobs from previous runs
3. Start processing queued calculations

### Terminal 2: Start the API Server

```bash
cd ~/develop/qm-nmr-calc
uv run python scripts/run_api.py
```

The API server runs on `http://localhost:8000` with:
- Web UI: `http://localhost:8000/`
- API docs: `http://localhost:8000/api/v1/openapi.json`
- Health check: `http://localhost:8000/health`

## Usage

### Web Interface

1. Open `http://localhost:8000/` in your browser
2. Enter a SMILES string or upload a MOL/SDF file
3. Select a solvent and calculation preset
4. Submit and wait for results

### REST API

**Submit a job:**
```bash
curl -X POST http://localhost:8000/api/v1/jobs \
  -H "Content-Type: application/json" \
  -d '{"smiles": "CCO", "solvent": "chcl3", "preset": "draft"}'
```

**Check job status:**
```bash
curl http://localhost:8000/api/v1/jobs/{job_id}
```

**Get results (when complete):**
```bash
# Chemical shifts JSON
curl http://localhost:8000/api/v1/jobs/{job_id}/results

# Optimized geometry (XYZ)
curl http://localhost:8000/api/v1/jobs/{job_id}/geometry

# Optimized geometry (SDF)
curl http://localhost:8000/api/v1/jobs/{job_id}/geometry.sdf

# Raw NWChem output (ZIP)
curl http://localhost:8000/api/v1/jobs/{job_id}/output -o output.zip

# Spectrum plots
curl http://localhost:8000/api/v1/jobs/{job_id}/spectrum/1H.svg -o spectrum_1H.svg
curl http://localhost:8000/api/v1/jobs/{job_id}/spectrum/13C.svg -o spectrum_13C.svg

# Annotated structure
curl http://localhost:8000/api/v1/jobs/{job_id}/structure.svg -o structure.svg
```

### Supported Solvents

| Code | Solvent |
|------|---------|
| `chcl3` | Chloroform (CDCl3) |
| `dmso` | DMSO (DMSO-d6) |

Solvent effects are applied via the COSMO solvation model to both geometry optimization and NMR shielding calculations.

### Calculation Presets

| Preset | Basis Set | NMR Basis | Use Case |
|--------|-----------|-----------|----------|
| `draft` | 6-31G* | 6-31G* | Quick checks, debugging |
| `production` | 6-31G* | 6-311+G(2d,p) | Publication quality (default) |

## Configuration

### Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `SMTP_HOST` | `localhost` | SMTP server for notifications |
| `SMTP_PORT` | `587` | SMTP port |
| `SMTP_USER` | (empty) | SMTP username (optional) |
| `SMTP_PASSWORD` | (empty) | SMTP password (optional) |
| `NOTIFICATION_FROM` | `noreply@qm-nmr-calc.example` | Sender email address |
| `BASE_URL` | `http://localhost:8000` | Base URL for links in emails |

Example:
```bash
export SMTP_HOST=smtp.gmail.com
export SMTP_PORT=587
export SMTP_USER=your-email@gmail.com
export SMTP_PASSWORD=your-app-password
export NOTIFICATION_FROM=your-email@gmail.com
export BASE_URL=https://your-domain.com
```

### Data Directory

Jobs are stored in `./data/jobs/{job_id}/`:
```
data/
├── huey.db           # Job queue database
├── jobs/
│   └── {job_id}/
│       ├── status.json      # Job metadata and results
│       ├── output/
│       │   ├── optimized.xyz
│       │   ├── spectrum_1H.svg
│       │   ├── spectrum_1H.png
│       │   ├── spectrum_13C.svg
│       │   ├── spectrum_13C.png
│       │   ├── structure_annotated.svg
│       │   └── structure_annotated.png
│       └── scratch/         # NWChem working files
│           ├── *.nw
│           └── *.out
└── tms_reference_values.json  # Pre-computed TMS shieldings
```

## Development

### Running Tests

```bash
# Quick tests (no NWChem required)
uv run pytest tests/ -v

# Full workflow test (requires running consumer)
uv run python scripts/test_workflow.py
```

### Project Structure

```
qm-nmr-calc/
├── src/qm_nmr_calc/
│   ├── api/
│   │   ├── app.py           # FastAPI application
│   │   ├── schemas.py       # Pydantic models
│   │   ├── routers/
│   │   │   ├── health.py    # Health check endpoints
│   │   │   ├── jobs.py      # Job API endpoints
│   │   │   └── web.py       # Web UI routes
│   │   ├── static/          # CSS, JS
│   │   └── templates/       # Jinja2 templates
│   ├── nwchem/              # NWChem integration module
│   │   ├── geometry.py      # SMILES to 3D, XYZ/SDF loading
│   │   ├── input_gen.py     # NWChem input file generation
│   │   ├── output_parser.py # Parse NWChem output files
│   │   └── runner.py        # Execute NWChem calculations
│   ├── models.py            # Data models
│   ├── notifications.py     # Email notifications
│   ├── presets.py           # Calculation presets
│   ├── queue.py             # Huey queue config
│   ├── shifts.py            # Chemical shift calculation
│   ├── solvents.py          # Solvent validation
│   ├── startup.py           # Environment validation
│   ├── storage.py           # Job file management
│   ├── tasks.py             # Huey task definitions
│   ├── validation.py        # Input validation
│   └── visualization.py     # Spectrum/structure plots
├── scripts/
│   ├── run_api.py           # Start API server
│   ├── run_consumer.py      # Start job consumer
│   └── test_workflow.py     # Integration test
├── tests/
├── data/                    # Job storage (gitignored)
├── .planning/               # Project planning docs
├── pyproject.toml
└── uv.lock
```

### Adjusting CPU Cores

Edit the preset in `src/qm_nmr_calc/presets.py`:

```python
PRESETS: dict[PresetName, CalculationPreset] = {
    PresetName.DRAFT: {
        ...
        "processes": 4,  # Change this to your core count
        ...
    },
    PresetName.PRODUCTION: {
        ...
        "processes": 4,  # Change this to your core count
        ...
    },
}
```

## Known Limitations

1. **Limited Solvents**: Only CHCl3 and DMSO are currently supported for COSMO solvation. Additional solvents can be added by specifying dielectric constants in `input_gen.py`.

2. **Chemical Shift Accuracy**: Predicted shifts may deviate from experimental values. Accuracy depends on basis set quality, conformational sampling, and reference compound calibration.

3. **Single Worker**: The Huey consumer runs with 1 worker because QM calculations are CPU-bound. Jobs are processed sequentially.

4. **Supported Nuclei**: Only 1H and 13C NMR shifts are calculated.

## Troubleshooting

### "FATAL: nwchem not found in PATH"

NWChem is not installed or not in PATH:
```bash
sudo apt-get install nwchem
which nwchem  # Should return path
```

### "FATAL: RDKit initialization failed"

RDKit is not installed correctly:
```bash
uv sync  # Reinstall dependencies
```

### Job stuck in "queued" status

The Huey consumer is not running:
```bash
uv run python scripts/run_consumer.py
```

### NWChem runs out of memory

Reduce the number of processes in presets.py or use the `draft` preset for large molecules.

## License

MIT

## Acknowledgments

This project's NWChem integration approach was informed by [ISiCLE](https://github.com/pnnl/isicle),
a Python package for high-throughput NMR chemical shift calculations developed at
Pacific Northwest National Laboratory.

- [NWChem](https://nwchemgit.github.io/) - Quantum chemistry package
- [RDKit](https://www.rdkit.org/) - Cheminformatics toolkit
