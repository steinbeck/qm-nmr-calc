# Stack Research: QM-NMR-Calc

**Project:** Async NMR calculation web service wrapping ISiCLE/NWChem
**Researched:** 2026-01-19
**Overall Confidence:** HIGH (verified via official documentation and PyPI)

## Executive Summary

This stack is designed for a single-VM deployment running long-running QM calculations. The key architectural decisions are:

1. **Huey with SQLite** for job queue — no Redis/Celery complexity, survives restarts, supports status tracking
2. **FastAPI** for async REST API — natural fit with Python ecosystem, excellent performance
3. **RDKit** for molecule visualization — industry standard, supports atom annotations for NMR shifts
4. **matplotlib** for spectrum plotting — simple, well-understood, sufficient for 1D NMR spectra
5. **ISiCLE + NWChem** as calculation engine — the established pipeline for automated NMR predictions

## Recommended Stack

### Core Framework

| Component | Choice | Version | Rationale |
|-----------|--------|---------|-----------|
| Web Framework | FastAPI | >=0.128.0 | Modern async Python, excellent performance (15k+ req/s), automatic OpenAPI docs, native Pydantic validation |
| ASGI Server | Uvicorn | >=0.34.0 | High-performance ASGI server, standard FastAPI deployment |
| Python | Python | >=3.10 | ISiCLE requires >=3.9, FastAPI drops 3.8, use 3.10+ for performance |

### Job Queue

| Component | Choice | Version | Rationale |
|-----------|--------|---------|-----------|
| Task Queue | Huey | >=2.5.0 | SQLite backend, no Redis required, supports status tracking, retries, scheduling |
| Storage | SqliteHuey | (built-in) | Zero external dependencies, WAL mode for concurrent reads, survives restarts |

**Why Huey over alternatives:**
- **vs Celery/RQ**: No Redis dependency, simpler deployment on single VM
- **vs FastAPI BackgroundTasks**: Needs job persistence, status tracking, retry support
- **vs litequeue/persist-queue**: Huey has mature worker management, better API

**SqliteHuey Configuration:**
```python
from huey import SqliteHuey

huey = SqliteHuey(
    'qm-nmr-calc',
    filename='/var/lib/qm-nmr-calc/huey.db',
    fsync=True,           # Durable writes for long calculations
    journal_mode='wal',   # Write-ahead logging for concurrent access
    timeout=30            # Higher timeout for busy periods
)
```

### Calculation Engine

| Component | Choice | Version | Rationale |
|-----------|--------|---------|-----------|
| NMR Pipeline | ISiCLE | >=2.0.0 | PNNL's automated NMR framework, handles conformer generation + DFT |
| QM Engine | NWChem | >=7.2.0 | Open-source DFT, ISiCLE's native backend, apt-installable |
| Molecule I/O | OpenBabel (pybel) | >=3.1.0 | Format conversion (SMILES, MOL, SDF), ISiCLE dependency |

**ISiCLE Dependencies (from requirements.txt):**
- numpy >=1.19.4
- pandas >=1.1.4
- snakemake >=6.3.0 (workflow engine)
- statsmodels >=0.11.1
- joblib (parallel processing)

### Visualization

| Component | Choice | Version | Rationale |
|-----------|--------|---------|-----------|
| Molecule Drawing | RDKit | >=2025.3.0 | Industry standard, supports `atomNote` for shift labels, PNG/SVG output |
| Spectrum Plotting | matplotlib | >=3.8.0 | Simple, well-known, sufficient for 1D NMR spectra |
| Peak Annotation | scipy | >=1.11.0 | `find_peaks()` for automatic peak detection |

**RDKit Atom Annotation Pattern:**
```python
from rdkit import Chem
from rdkit.Chem import Draw

mol = Chem.MolFromSmiles(smiles)
for atom in mol.GetAtoms():
    shift = shifts.get(atom.GetIdx())
    if shift:
        atom.SetProp('atomNote', f'{shift:.1f}')

img = Draw.MolToImage(mol, size=(400, 400))
```

### Email Notifications

| Component | Choice | Version | Rationale |
|-----------|--------|---------|-----------|
| Async Email | aiosmtplib | >=5.0.0 | Async SMTP client, works with FastAPI's async model |
| Email Building | email.message | (stdlib) | Standard library EmailMessage, no extra dependencies |

### Web UI

| Component | Choice | Version | Rationale |
|-----------|--------|---------|-----------|
| Templating | Jinja2 | >=3.1.0 | FastAPI native support, familiar syntax |
| Static Files | FastAPI StaticFiles | (built-in) | Built-in static file serving |
| CSS Framework | Pico CSS or Tailwind | latest | Lightweight, no-JS styling |

### Data Validation & Serialization

| Component | Choice | Version | Rationale |
|-----------|--------|---------|-----------|
| Validation | Pydantic | >=2.5.0 | FastAPI native, excellent for API schemas |
| JSON | orjson | >=3.9.0 | Fast JSON serialization for API responses |

### Storage

| Component | Choice | Rationale |
|-----------|--------|-----------|
| Job Metadata | JSON files | Filesystem-based per project constraints, one JSON per job |
| Calculation Files | Organized directories | `/jobs/{job_id}/` with inputs, outputs, logs |
| Job Queue | SQLite (via Huey) | Single file, concurrent-safe with WAL |

**Directory Structure:**
```
/var/lib/qm-nmr-calc/
  huey.db              # Job queue database
  jobs/
    {job_id}/
      input.json       # Job parameters
      molecule.sdf     # Input structure
      status.json      # Job status/progress
      output/
        nwchem.out     # Raw NWChem output
        shifts.json    # Parsed chemical shifts
        spectrum.png   # Generated spectrum plot
        structure.png  # Annotated structure
```

## Installation Notes

### NWChem Installation (Ubuntu/Debian)

**Simple installation:**
```bash
sudo apt-get update
sudo apt-get install nwchem
```

**Verify installation:**
```bash
nwchem --version
# Or run a test calculation
echo "start h2o; geometry; O 0 0 0; H 0 0 1; H 0 1 0; end; task scf" > test.nw
nwchem test.nw
```

**Important considerations:**
- Default apt package uses OpenMPI, fine for single-node
- For multi-core: `mpirun -np 4 nwchem input.nw`
- NWChem is memory-hungry: plan 4-8GB RAM minimum for small molecules

### ISiCLE Installation

**Recommended: Install from GitHub:**
```bash
# Create virtual environment
python -m venv venv
source venv/bin/activate

# Install ISiCLE
pip install git+https://github.com/pnnl/isicle.git

# Or if maintaining a fork:
pip install -e /path/to/isicle-fork
```

**Install additional dependencies:**
```bash
# OpenBabel for molecule I/O
conda install -c conda-forge openbabel
# Or via apt
sudo apt-get install openbabel python3-openbabel

# RDKit (commented out in ISiCLE but we need it)
pip install rdkit
```

**ISiCLE Configuration:**
ISiCLE uses YAML configuration files. Key settings:
```yaml
nwchem:
  scratch: /tmp/nwchem_scratch
  memory: "4 gb"

dft:
  functional: B3LYP
  basis: 6-31G*

nmr:
  nuclei: [H, C]
  reference: TMS  # Tetramethylsilane
```

### Python Environment

**Full requirements.txt:**
```
# Web framework
fastapi>=0.128.0
uvicorn[standard]>=0.34.0
python-multipart>=0.0.9  # File uploads

# Job queue
huey>=2.5.0

# Calculation engine
# ISiCLE installed from git, brings:
# - numpy, pandas, snakemake, statsmodels, joblib

# Visualization
rdkit>=2025.3.0
matplotlib>=3.8.0
scipy>=1.11.0

# Email
aiosmtplib>=5.0.0

# Utilities
pydantic>=2.5.0
orjson>=3.9.0
jinja2>=3.1.0
python-dotenv>=1.0.0
```

### Running the Stack

**Development:**
```bash
# Terminal 1: API server
uvicorn app.main:app --reload --host 0.0.0.0 --port 8000

# Terminal 2: Huey worker (single worker for QM calculations)
huey_consumer.py app.tasks.huey -w 1 -k process
```

**Production:**
```bash
# API server (multiple workers)
uvicorn app.main:app --host 0.0.0.0 --port 8000 --workers 4

# Huey worker (single worker - QM calculations are CPU-bound)
# Use process mode, not greenlet, for CPU-intensive work
huey_consumer.py app.tasks.huey -w 1 -k process
```

**Systemd service example:**
```ini
[Unit]
Description=QM-NMR-Calc Huey Worker
After=network.target

[Service]
Type=simple
User=www-data
WorkingDirectory=/opt/qm-nmr-calc
ExecStart=/opt/qm-nmr-calc/venv/bin/huey_consumer.py app.tasks.huey -w 1 -k process
Restart=always

[Install]
WantedBy=multi-user.target
```

## Alternatives Considered

### Web Framework

| Option | Why Not |
|--------|---------|
| Flask | Sync-first, would need WSGI adapter, less performant |
| Django | Overkill for API-focused service, brings ORM we don't need |
| Starlette | FastAPI adds Pydantic validation, OpenAPI docs for free |

### Job Queue

| Option | Why Not |
|--------|---------|
| Celery + Redis | Requires Redis, complex setup for single VM |
| RQ (Redis Queue) | Requires Redis |
| ARQ | Requires Redis, async-focused (our tasks are CPU-bound sync) |
| FastAPI BackgroundTasks | No persistence, no status tracking, no retry |
| APScheduler | Scheduling-focused, not task queue |
| Django-Q | Django dependency |
| Procrastinate | PostgreSQL dependency |

### Visualization

| Option | Why Not |
|--------|---------|
| py3Dmol | 3D interactive, overkill for simple 2D annotated structures |
| nmrglue | Great for experimental NMR data, but we're generating synthetic spectra |
| Bokeh/Plotly | Interactive plots, heavier than needed for static images |

### Molecule Handling

| Option | Why Not |
|--------|---------|
| CDK (via JPype) | Java dependency, RDKit is native Python |
| PyMOL | 3D visualization focus, licensing complexity |

## Dependencies Summary

### Core Dependencies

```
fastapi[standard]>=0.128.0   # Includes uvicorn, python-multipart
huey>=2.5.0                   # Task queue with SQLite
rdkit>=2025.3.0              # Molecule visualization
matplotlib>=3.8.0            # Spectrum plotting
scipy>=1.11.0                # Peak detection
aiosmtplib>=5.0.0            # Async email
pydantic>=2.5.0              # Data validation
jinja2>=3.1.0                # Templates
orjson>=3.9.0                # Fast JSON
python-dotenv>=1.0.0         # Environment config
```

### ISiCLE Dependencies (installed with ISiCLE)

```
numpy>=1.19.4
pandas>=1.1.4
snakemake>=6.3.0
statsmodels>=0.11.1
joblib
openbabel (system or conda)
```

### System Dependencies

```bash
# Ubuntu/Debian
sudo apt-get install nwchem openbabel python3-openbabel
```

## Configuration Management

**Recommended: Environment variables + .env file:**
```bash
# .env
QM_NMR_JOBS_DIR=/var/lib/qm-nmr-calc/jobs
QM_NMR_HUEY_DB=/var/lib/qm-nmr-calc/huey.db
QM_NMR_NWCHEM_SCRATCH=/tmp/nwchem
QM_NMR_NWCHEM_MEMORY=4gb

# Email (optional)
SMTP_HOST=smtp.example.com
SMTP_PORT=587
SMTP_USER=notifications@example.com
SMTP_PASSWORD=secret
NOTIFICATION_FROM=noreply@example.com
```

## Sources

### Official Documentation
- [FastAPI Documentation](https://fastapi.tiangolo.com/)
- [FastAPI PyPI](https://pypi.org/project/fastapi/) - Version 0.128.0
- [Huey Documentation](https://huey.readthedocs.io/)
- [Huey SQLite Storage](https://huey.readthedocs.io/en/1.11.0/sqlite.html)
- [RDKit Documentation](https://www.rdkit.org/docs/) - Version 2025.09.4
- [RDKit PyPI](https://pypi.org/project/rdkit/) - Version 2025.9.3
- [NWChem Download](https://nwchemgit.github.io/Download.html)
- [aiosmtplib Documentation](https://aiosmtplib.readthedocs.io/) - Version 5.0.0
- [OpenBabel Pybel](https://open-babel.readthedocs.io/en/latest/UseTheLibrary/Python_Pybel.html)

### ISiCLE
- [ISiCLE GitHub](https://github.com/pnnl/isicle)
- [ISiCLE Paper](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-018-0305-8) - Journal of Cheminformatics, 2018

### Background Research
- [FastAPI Background Tasks](https://fastapi.tiangolo.com/tutorial/background-tasks/)
- [Managing Background Tasks in FastAPI](https://leapcell.io/blog/managing-background-tasks-and-long-running-operations-in-fastapi)
- [Python Task Queue Comparison](https://judoscale.com/blog/choose-python-task-queue)
- [litequeue](https://github.com/litements/litequeue) - SQLite queue reference
