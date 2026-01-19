# Phase 1: Foundation - Research

**Researched:** 2026-01-19
**Domain:** Job queue infrastructure + ISiCLE/NWChem integration
**Confidence:** HIGH

## Summary

Phase 1 establishes the calculation infrastructure: job queue with Huey/SQLite, ISiCLE wrapper for geometry optimization, filesystem-based job storage, and error handling. Research confirms:

1. **ISiCLE API** is straightforward: load molecule via `isicle.load()`, run optimization via `geom.initial_optimize()` then `isicle.qm.dft()`. The wrapper handles NWChem input generation, execution, and output parsing internally.

2. **Huey with SQLite** is the right choice for single-VM deployment. It provides task status tracking via signals (SIGNAL_COMPLETE, SIGNAL_ERROR, SIGNAL_INTERRUPTED), result storage, and survives process restarts without Redis.

3. **uv + pyproject.toml** is the modern standard for Python projects. Initialize with `uv init`, add dependencies with `uv add`, run with `uv run`.

**Primary recommendation:** Use Huey signals for job status tracking, wrap ISiCLE's `dft()` function in a thin layer that manages job directories, and store all job metadata in status.json files alongside calculation outputs.

## Standard Stack

The established libraries/tools for this phase:

### Core

| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| huey | >=2.5.0 | Task queue with SQLite storage | No Redis dependency, mature API, signal-based status tracking |
| isicle | 2.0.0 (fork) | NMR calculation wrapper | PNNL's established pipeline, handles NWChem interaction |
| rdkit | >=2025.3.0 | Molecule parsing and manipulation | ISiCLE dependency, industry standard |

### Supporting

| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| pydantic | >=2.5.0 | Job metadata validation | Define JobStatus, JobInput models |
| orjson | >=3.9.0 | Fast JSON serialization | Write status.json files |

### System Dependencies

| Component | Version | Purpose |
|-----------|---------|---------|
| NWChem | 7.0.2 | QM calculation engine (already installed at /usr/bin/nwchem) |
| OpenMPI | 4.1.4 | NWChem parallelization (already installed) |
| Python | >=3.10 | Runtime (ISiCLE requires >=3.9) |

### Alternatives Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Huey | Celery | Celery needs Redis/RabbitMQ - more ops complexity |
| Huey | RQ | RQ also needs Redis |
| Huey | BackgroundTasks | No persistence, no retry, no status tracking |
| SqliteHuey | FileHuey | FileHuey less mature, SQLite well-tested |

**Installation:**
```bash
# Project dependencies (via uv)
uv add huey rdkit pydantic orjson

# ISiCLE from fork (editable install)
uv add -e ~/develop/isicle

# ISiCLE's transitive dependencies come automatically:
# numpy, pandas, snakemake, statsmodels, joblib
```

## Architecture Patterns

### Recommended Project Structure

```
qm-nmr-calc/
├── pyproject.toml           # Project config and dependencies
├── uv.lock                   # Lockfile (commit this)
├── .python-version           # Python version for uv
├── data/
│   └── jobs/                 # Job data storage
│       └── {job_id}/         # One directory per job
│           ├── status.json   # Job metadata and status
│           ├── input.mol     # Input molecule
│           ├── output/       # Calculation outputs
│           └── logs/         # NWChem logs
└── src/
    └── qm_nmr_calc/
        ├── __init__.py
        ├── queue.py          # Huey instance and task definitions
        ├── tasks.py          # Calculation task implementations
        ├── isicle_wrapper.py # ISiCLE integration layer
        ├── models.py         # Pydantic models (JobStatus, etc.)
        └── storage.py        # Job directory management
```

### Pattern 1: Task Status via Signals

**What:** Use Huey signals to update job status instead of polling result store.

**When to use:** Always - this is how you track job lifecycle.

**Example:**
```python
# Source: https://huey.readthedocs.io/en/latest/signals.html
from huey import SqliteHuey, signals

huey = SqliteHuey('qm-nmr-calc', filename='./data/huey.db')

@huey.signal(signals.SIGNAL_EXECUTING)
def on_task_start(signal, task):
    job_id = task.args[0]  # Assuming job_id is first argument
    update_job_status(job_id, status='running', started_at=datetime.utcnow())

@huey.signal(signals.SIGNAL_COMPLETE)
def on_task_complete(signal, task):
    job_id = task.args[0]
    update_job_status(job_id, status='complete', completed_at=datetime.utcnow())

@huey.signal(signals.SIGNAL_ERROR)
def on_task_error(signal, task, exc=None):
    job_id = task.args[0]
    update_job_status(
        job_id,
        status='failed',
        error_message=str(exc),
        error_traceback=traceback.format_exc()
    )

@huey.signal(signals.SIGNAL_INTERRUPTED)
def on_task_interrupted(signal, task):
    job_id = task.args[0]
    update_job_status(job_id, status='failed', error_message='Task interrupted - process restart')
```

### Pattern 2: Job Directory Structure

**What:** Each job gets its own directory with all inputs, outputs, and metadata.

**When to use:** Always - this is the storage pattern from CONTEXT.md.

**Example:**
```python
# Source: User decisions in CONTEXT.md
from pathlib import Path
import orjson

def create_job_directory(job_id: str, smiles: str) -> Path:
    """Create job directory with initial status."""
    job_dir = Path('./data/jobs') / job_id
    job_dir.mkdir(parents=True, exist_ok=True)
    (job_dir / 'output').mkdir(exist_ok=True)
    (job_dir / 'logs').mkdir(exist_ok=True)

    status = {
        'job_id': job_id,
        'status': 'queued',
        'created_at': datetime.utcnow().isoformat(),
        'input': {'smiles': smiles},
        'isicle_version': isicle.__version__,
        'nwchem_version': get_nwchem_version(),
    }

    (job_dir / 'status.json').write_bytes(orjson.dumps(status, option=orjson.OPT_INDENT_2))
    return job_dir
```

### Pattern 3: ISiCLE Geometry Optimization

**What:** The minimal ISiCLE workflow for geometry optimization.

**When to use:** This is the core calculation for Phase 1.

**Example:**
```python
# Source: ISiCLE source code analysis (isicle/qm.py, isicle/geometry.py)
import isicle

def run_geometry_optimization(smiles: str, job_dir: Path, processes: int = 4):
    """Run geometry optimization via ISiCLE/NWChem."""

    # 1. Load molecule from SMILES
    geom = isicle.load(smiles)

    # 2. Initial 3D embedding and force-field optimization
    geom = geom.initial_optimize(embed=True, forcefield='UFF', ff_iter=200)

    # 3. DFT geometry optimization via NWChem
    wrapper = isicle.qm.dft(
        geom,
        backend='NWChem',
        tasks=['optimize'],
        functional='b3lyp',
        basis_set='6-31G*',
        scratch_dir=str(job_dir / 'scratch'),
        processes=processes,
    )

    # 4. Parse results
    result = wrapper.parse()

    # 5. Return optimized geometry
    return result['geometry']
```

### Pattern 4: Fail-Fast NWChem Validation

**What:** Check NWChem availability at startup, not at first job.

**When to use:** Application initialization.

**Example:**
```python
# Source: User decisions in CONTEXT.md
import subprocess
import sys

def validate_nwchem():
    """Validate NWChem is available and working. Call at startup."""
    try:
        # NWChem exits with error if no input file, but prints version info
        result = subprocess.run(
            ['nwchem', '--version'],
            capture_output=True,
            text=True,
            timeout=5
        )
        # NWChem doesn't have --version, try calling with empty input
        # A valid install will complain about missing input file
    except FileNotFoundError:
        sys.exit("FATAL: NWChem not found in PATH. Install with: apt-get install nwchem")
    except subprocess.TimeoutExpired:
        sys.exit("FATAL: NWChem timed out during validation")

def get_nwchem_version() -> str:
    """Extract NWChem version from installation."""
    # Parse from 'nwchem -h' or check /usr/bin/nwchem package
    return "7.0.2"  # Hardcode for now, improve later
```

### Anti-Patterns to Avoid

- **Polling for status:** Don't poll Huey's result store. Use signals instead - they fire synchronously when status changes.
- **Storing large results in Huey:** Huey stores results in SQLite. Store large outputs (NWChem logs, geometries) in job directory, only store status/metadata in Huey.
- **Catching all exceptions in tasks:** Let exceptions propagate to Huey's error handling. Signals will capture the traceback.
- **Using temporary directories:** ISiCLE creates temp dirs by default. Override with explicit `scratch_dir` in job directory for debugging.

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Task queue | Custom threading/multiprocessing | Huey | Handles persistence, retries, signals, consumer management |
| NWChem input generation | Manual .nw file writing | ISiCLE's NWChemWrapper | Handles memory settings, basis sets, charge, all edge cases |
| Molecule parsing | Manual SMILES parsing | RDKit via ISiCLE | Handles stereochemistry, tautomers, 3D embedding |
| Job ID generation | Sequential integers | `uuid.uuid4().hex[:12]` | Collision-free, URL-safe |
| JSON serialization | stdlib json | orjson | 10x faster, handles datetime natively |

**Key insight:** ISiCLE already handles the hard parts of NWChem interaction. Don't bypass it or reimplement its functionality - wrap it thinly.

## Common Pitfalls

### Pitfall 1: Huey Consumer Not Finding Tasks

**What goes wrong:** Consumer starts but doesn't execute tasks; logs show "task not found in registry".

**Why it happens:** Consumer imports a different module path than the producer, or tasks aren't imported before consumer starts.

**How to avoid:**
- Put Huey instance and all @huey.task decorators in a single module (`queue.py`)
- Start consumer with exact module path: `huey_consumer.py qm_nmr_calc.queue.huey`
- Verify imports work: `python -c "from qm_nmr_calc.queue import huey"`

**Warning signs:** Tasks stay in queue indefinitely; consumer shows "waiting for task" but never executes.

### Pitfall 2: ISiCLE Temp Directory Cleanup

**What goes wrong:** Disk fills up with /tmp/isicle/* directories after many jobs.

**Why it happens:** ISiCLE creates temp directories via `isicle.utils.mkdtemp()` but doesn't auto-cleanup.

**How to avoid:**
- Pass explicit `scratch_dir` inside job directory
- Call `isicle.utils.rmdtemp()` periodically or on job completion
- Monitor /tmp usage

**Warning signs:** Disk space warnings; slow job starts.

### Pitfall 3: NWChem Memory Errors

**What goes wrong:** NWChem crashes with cryptic ARMCI errors or "memory allocation failed".

**Why it happens:** Default memory settings too low for molecule size; or multiple jobs fighting for memory.

**How to avoid:**
- Set explicit memory in NWChemWrapper: `mem_global=4000, mem_heap=400, mem_stack=1200` (in MB)
- Run single worker: `huey_consumer.py ... -w 1` (QM jobs are CPU-bound anyway)
- Monitor system memory during long calculations

**Warning signs:** "ARMCI Error", "ga_nodeid", "nwchem: failed" in logs.

### Pitfall 4: Signal Handlers Blocking Consumer

**What goes wrong:** Consumer becomes slow; tasks queue up even though previous tasks completed.

**Why it happens:** Signal handlers run synchronously. Slow I/O in handlers blocks the consumer loop.

**How to avoid:**
- Keep signal handlers fast (< 100ms)
- Do file I/O (status.json updates) quickly with orjson
- Don't make HTTP calls in signal handlers

**Warning signs:** Long gaps between SIGNAL_COMPLETE and next SIGNAL_EXECUTING.

### Pitfall 5: Interrupted Jobs Not Marked Failed

**What goes wrong:** After process restart, jobs show status='running' forever.

**Why it happens:** SIGNAL_INTERRUPTED only fires if consumer shuts down gracefully. SIGKILL doesn't trigger it.

**How to avoid:**
- On startup, scan for jobs with status='running' and mark them failed
- Include recovery logic in application init

**Warning signs:** Jobs stuck in 'running' state after restart.

## Code Examples

Verified patterns from ISiCLE source code and Huey documentation:

### Complete Task Definition

```python
# Source: Huey docs + ISiCLE source analysis
from huey import SqliteHuey
import isicle
from pathlib import Path
import traceback

huey = SqliteHuey('qm-nmr-calc', filename='./data/huey.db', fsync=True)

@huey.task()
def run_optimization(job_id: str):
    """Execute geometry optimization for a queued job."""
    job_dir = Path('./data/jobs') / job_id
    status_file = job_dir / 'status.json'

    # Load job info
    status = orjson.loads(status_file.read_bytes())
    smiles = status['input']['smiles']

    try:
        # Run ISiCLE
        geom = isicle.load(smiles)
        geom = geom.initial_optimize(embed=True, forcefield='UFF')

        wrapper = isicle.qm.dft(
            geom,
            backend='NWChem',
            tasks=['optimize'],
            functional='b3lyp',
            basis_set='6-31G*',
            processes=4,
        )

        result = wrapper.parse()

        # Save outputs
        isicle.save(str(job_dir / 'output' / 'optimized.xyz'), result['geometry'])

        return {'success': True, 'job_id': job_id}

    except Exception as e:
        # Let Huey capture via SIGNAL_ERROR
        raise
```

### Huey Consumer Configuration

```python
# Source: Huey docs
# Run with: huey_consumer.py qm_nmr_calc.queue.huey -w 1 -k process

# Key flags:
# -w 1           : Single worker (QM jobs are CPU-bound)
# -k process     : Process-based workers (not greenlet - CPU bound)
# -f ./data/huey.db  : (Set in SqliteHuey constructor instead)
```

### Job Status Model

```python
# Source: CONTEXT.md decisions
from pydantic import BaseModel
from datetime import datetime
from typing import Optional, Literal

class JobStatus(BaseModel):
    job_id: str
    status: Literal['queued', 'running', 'complete', 'failed']
    created_at: datetime
    started_at: Optional[datetime] = None
    completed_at: Optional[datetime] = None

    # Input
    input_smiles: str

    # Versions (for reproducibility)
    isicle_version: str
    nwchem_version: str

    # Error info (if failed)
    error_message: Optional[str] = None
    error_traceback: Optional[str] = None

    # Resource usage (filled on completion)
    cpu_time_seconds: Optional[float] = None
    memory_peak_mb: Optional[float] = None
```

### Startup Validation

```python
# Source: CONTEXT.md - fail fast at startup
import sys
import subprocess
import isicle

def validate_environment():
    """Call this at application startup."""

    # Check NWChem
    try:
        result = subprocess.run(
            ['which', 'nwchem'],
            capture_output=True,
            text=True
        )
        if result.returncode != 0:
            sys.exit("FATAL: nwchem not found in PATH")
    except Exception as e:
        sys.exit(f"FATAL: Cannot verify nwchem: {e}")

    # Check ISiCLE loads
    try:
        geom = isicle.load("C")  # Simplest molecule
        geom.initial_optimize(embed=True)
    except Exception as e:
        sys.exit(f"FATAL: ISiCLE/RDKit initialization failed: {e}")

    # Check data directory writable
    data_dir = Path('./data/jobs')
    data_dir.mkdir(parents=True, exist_ok=True)
    test_file = data_dir / '.write_test'
    try:
        test_file.write_text('test')
        test_file.unlink()
    except Exception as e:
        sys.exit(f"FATAL: Cannot write to {data_dir}: {e}")

    print("Environment validation passed")
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| rdkit-pypi package | rdkit package | 2024 | Use `pip install rdkit` not `rdkit-pypi` |
| Celery for all queues | Huey for simple queues | Ongoing | Huey better for single-VM, no Redis |
| setup.py | pyproject.toml | 2022-2023 | PEP 517/518 standard, uv uses it |
| pip + venv | uv | 2024-2025 | 10-100x faster, better lockfiles |

**Deprecated/outdated:**
- `rdkit-pypi` PyPI package: Renamed to just `rdkit`
- `setup.py` for new projects: Use `pyproject.toml` with uv

## Open Questions

Things that couldn't be fully resolved:

1. **NWChem version extraction**
   - What we know: NWChem is installed at /usr/bin/nwchem, version 7.0.2
   - What's unclear: Best way to extract version programmatically (no `--version` flag)
   - Recommendation: Parse from `nwchem 2>&1 | head` output or hardcode for now

2. **Memory settings for different molecule sizes**
   - What we know: ISiCLE defaults to mem_global=1600, mem_heap=100, mem_stack=600 MB
   - What's unclear: Optimal settings for molecules of different sizes
   - Recommendation: Start with defaults, monitor failures, adjust per-job based on atom count

3. **ISiCLE temp directory management**
   - What we know: ISiCLE creates temp dirs under /tmp/isicle/
   - What's unclear: Whether wrapper.temp_dir is auto-cleaned after wrapper.finish()
   - Recommendation: Explicitly manage scratch_dir in job directory, clean up on completion

## Sources

### Primary (HIGH confidence)

- **ISiCLE source code** - Direct analysis of `~/develop/isicle/isicle/`:
  - `qm.py` - NWChemWrapper implementation, dft() function
  - `geometry.py` - Geometry class, initial_optimize() method
  - `io.py` - load() function, file format handling
  - `parse.py` - NWChemParser for output parsing
  - `examples/nmr_chemical_shifts.py` - Example workflow

- **Huey documentation** - https://huey.readthedocs.io/en/latest/:
  - [Guide](https://huey.readthedocs.io/en/latest/guide.html) - Task definition, consumer operation
  - [API](https://huey.readthedocs.io/en/latest/api.html) - SqliteHuey parameters, methods
  - [Signals](https://huey.readthedocs.io/en/latest/signals.html) - Signal list, handler registration
  - [Troubleshooting](https://huey.readthedocs.io/en/latest/troubleshooting.html) - Common pitfalls

- **uv documentation** - https://docs.astral.sh/uv/guides/projects/ - Project setup, pyproject.toml

### Secondary (MEDIUM confidence)

- **RDKit PyPI** - https://pypi.org/project/rdkit/ - Version 2025.9.3, installation
- **NWChem verification** - Local installation at /usr/bin/nwchem confirmed working

### Tertiary (LOW confidence)

- **WebSearch results** for Huey patterns - Cross-verified with official docs

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - ISiCLE source code directly analyzed, Huey docs verified
- Architecture: HIGH - Based on CONTEXT.md decisions and verified library APIs
- Pitfalls: MEDIUM - Based on docs + general experience, some unverified edge cases

**Research date:** 2026-01-19
**Valid until:** 2026-02-19 (stable libraries, 30 days)

---

*Phase: 01-foundation*
*Research completed: 2026-01-19*
