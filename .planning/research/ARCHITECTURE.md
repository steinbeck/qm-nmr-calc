# Architecture Research: Async NMR Calculation Service

**Project:** qm-nmr-calc
**Researched:** 2026-01-19
**Confidence:** HIGH (verified with official documentation and established patterns)

## System Overview

The architecture follows the **async job submission pattern** - a well-established approach for long-running scientific computations. A FastAPI web service accepts molecule inputs, immediately queues jobs with Huey (SQLite backend), and returns job IDs. A separate Huey consumer process executes calculations via ISiCLE/NWChem, storing results to filesystem. Clients poll for status or receive email notifications on completion.

```
                                     +------------------+
                                     |   Email (SMTP)   |
                                     +--------^---------+
                                              |
+----------+      +-----------+      +--------+--------+      +-----------+
|  Client  | ---> |  FastAPI  | ---> |  Huey Consumer  | ---> |  ISiCLE/  |
|  (Web/   | <--- |  (API +   | <--- |  (Worker        | <--- |  NWChem   |
|   API)   |      |   Web UI) |      |   Process)      |      +-----------+
+----------+      +-----+-----+      +--------+--------+
                        |                     |
                        v                     v
                  +-----+-----+      +--------+--------+
                  |  SQLite   |      |   Filesystem    |
                  |  (Huey    |      |   (Job dirs,    |
                  |   Queue)  |      |   results)      |
                  +-----------+      +-----------------+
```

**Why this architecture:**
- **Separation of concerns:** API process stays responsive; heavy computation happens in dedicated worker
- **Resilience:** Jobs survive API restarts (persisted in SQLite)
- **Simplicity:** Single VM, no Redis/RabbitMQ, no Celery complexity
- **Scalability path:** Can add more workers later without architecture changes

## Components

### 1. FastAPI Application (API + Web UI)

**Purpose:** HTTP interface for job submission, status queries, and result retrieval.

**Responsibilities:**
- Accept molecule input (SMILES string or file upload)
- Validate input format
- Generate unique job ID (UUID)
- Create job directory with input files
- Enqueue calculation task to Huey
- Serve job status and results
- Render web UI for browser access

**Interfaces:**

| Endpoint | Method | Purpose |
|----------|--------|---------|
| `POST /api/jobs` | Submit | Accept molecule, return job_id |
| `GET /api/jobs/{id}` | Status | Return job status, metadata |
| `GET /api/jobs/{id}/results` | Results | Return calculated shifts (JSON) |
| `GET /api/jobs/{id}/files/{name}` | Download | Download specific result file |
| `GET /` | Web UI | Job submission form |
| `GET /jobs/{id}` | Web UI | Job status/results page |

**Key design decisions:**
- Job ID generated *before* enqueueing (allows immediate directory creation)
- Input files copied to job directory at submission time (ensures calculation has stable inputs)
- Status stored both in Huey result storage AND job directory (redundancy)

### 2. Huey Task Queue (SQLite Backend)

**Purpose:** Persistent job queue with task scheduling and result storage.

**Responsibilities:**
- Queue calculation tasks
- Track task state (pending, running, complete, error)
- Store task results (success/failure metadata)
- Support task priorities (optional)
- Persist across restarts

**Configuration:**

```python
from huey import SqliteHuey

huey = SqliteHuey(
    filename='data/huey.db',   # Persistent queue
    cache_mb=64,               # Page cache size
    fsync=True,                # Durable writes (important for jobs)
)
```

**Why Huey over alternatives:**
- Built-in SQLite support (no Redis required)
- Clean, simple API
- Task priorities supported
- Result storage included
- Active maintenance (huey 2.6.0+)
- Lightweight: ~2000 lines of code

**Why NOT:**
- `FastAPI.BackgroundTasks`: Tasks lost on restart, no persistence
- `Celery`: Requires Redis/RabbitMQ, complex for single-VM
- `RQ`: Requires Redis
- `ARQ`: Requires Redis

### 3. Huey Consumer (Worker Process)

**Purpose:** Execute queued calculation tasks.

**Responsibilities:**
- Pull tasks from queue
- Execute ISiCLE/NWChem calculations
- Write results to job directory
- Update job status file
- Trigger email notifications
- Handle errors and timeouts

**Execution model:**

```bash
# Single worker (default, recommended for NMR calculations)
huey_consumer.py app.tasks.huey -k process -w 1

# Multiple workers (if VM has resources)
huey_consumer.py app.tasks.huey -k process -w 2
```

**Why single worker initially:**
- NWChem calculations are CPU/memory intensive
- Single VM has limited resources
- Prevents resource contention
- Can increase `-w` count later if needed

### 4. ISiCLE Integration Layer

**Purpose:** Wrap ISiCLE library calls for the calculation pipeline.

**Responsibilities:**
- Convert input (SMILES/MOL) to ISiCLE-compatible format
- Configure calculation parameters (basis set, functional, solvent)
- Execute geometry optimization
- Execute NMR shielding calculation
- Parse NWChem output files
- Convert shielding values to chemical shifts
- Handle calculation failures gracefully

**Output files generated:**
- `input.smi` or `input.mol` - Original molecule input
- `*.xyz` - Generated 3D coordinates
- `*.nw` - NWChem input file
- `*.nwo` - NWChem output (raw calculation results)
- `*.log` - Calculation log/stderr
- `shifts.json` - Parsed chemical shifts (our format)

### 5. Filesystem Storage Layer

**Purpose:** Persistent storage for job data and results.

**Responsibilities:**
- Create job directories
- Store input files
- Store calculation outputs
- Store job metadata
- Enable direct file downloads

**Design rationale:**
- Filesystem is the simplest "database" for v1
- Easy inspection/debugging (just `ls` and `cat`)
- Natural fit for calculation output files
- Easy backup (rsync/tar)
- Supports multi-user extension (permissions)

### 6. Email Notification Service

**Purpose:** Notify users when calculations complete.

**Responsibilities:**
- Send completion email with result summary
- Send error email if calculation fails
- Include job URL for result access

**Implementation:**

```python
# Use aiosmtplib for async email within Huey task
import aiosmtplib
from email.message import EmailMessage

async def send_notification(job_id: str, email: str, status: str):
    message = EmailMessage()
    message["From"] = "nmr-calc@example.com"
    message["To"] = email
    message["Subject"] = f"NMR Calculation {status}: {job_id}"
    message.set_content(f"Your calculation is ready: https://example.com/jobs/{job_id}")

    await aiosmtplib.send(message, hostname=SMTP_HOST, port=SMTP_PORT)
```

**Why aiosmtplib:**
- Native async support
- No external dependencies (just SMTP server)
- Simple API
- Works with any SMTP provider (SendGrid, SES, local postfix)

## Data Flow

### Job Lifecycle: Submission to Completion

```
1. SUBMIT
   Client --> POST /api/jobs {smiles: "CCO", email: "user@example.com"}

2. VALIDATE
   API validates SMILES format (RDKit/OpenBabel)

3. CREATE JOB
   API generates job_id = uuid4()
   API creates /data/jobs/{job_id}/
   API writes input.smi with SMILES
   API writes meta.json with {status: "pending", email: "...", submitted: "..."}

4. ENQUEUE
   API calls huey task: run_nmr_calculation.schedule(args=(job_id,))
   Task ID stored in meta.json

5. RESPOND
   API returns {job_id: "abc-123", status: "pending"}

6. [Later] POLL STATUS
   Client --> GET /api/jobs/abc-123
   API reads meta.json, returns status

7. [Worker] EXECUTE
   Huey consumer picks up task
   Worker updates meta.json {status: "running", started: "..."}
   Worker invokes ISiCLE calculation
   ISiCLE writes output files to job directory
   Worker parses results, writes shifts.json
   Worker updates meta.json {status: "complete", finished: "..."}

8. [Worker] NOTIFY
   Worker sends email notification via aiosmtplib

9. RETRIEVE RESULTS
   Client --> GET /api/jobs/abc-123/results
   API returns shifts.json content
   Client --> GET /api/jobs/abc-123/files/output.nwo
   API returns raw NWChem output
```

### Status State Machine

```
pending --> running --> complete
              |            |
              v            v
           failed       (done)
              |
              v
           (done)
```

**Status values:**
- `pending`: Job queued, not yet started
- `running`: Worker executing calculation
- `complete`: Calculation finished successfully
- `failed`: Calculation failed (error in meta.json)

## File Structure

### Recommended Directory Layout

```
/data/
  huey.db                       # Huey SQLite queue database
  jobs/
    {job_id}/                   # One directory per job
      meta.json                 # Job metadata and status
      input.smi                 # Input SMILES (or input.mol)
      params.json               # Calculation parameters used

      # ISiCLE/NWChem working files
      mol.xyz                   # 3D coordinates
      mol.nw                    # NWChem input
      mol.nwo                   # NWChem output (raw)
      mol.log                   # Calculation stdout/stderr

      # Processed results
      shifts.json               # Parsed chemical shifts
      spectrum_1h.png           # Optional: spectrum plot
      spectrum_13c.png
      structure.svg             # Optional: annotated structure

      # Error handling
      error.txt                 # Present if status=failed

/app/
  main.py                       # FastAPI application
  tasks.py                      # Huey task definitions
  config.py                     # Configuration

  api/
    routes.py                   # API endpoints
    models.py                   # Pydantic models

  core/
    isicle_wrapper.py           # ISiCLE integration
    nmr_parser.py               # NWChem output parsing
    email_service.py            # Notification sender

  web/
    templates/                  # Jinja2 templates
    static/                     # CSS, JS
```

### Job Directory Rationale

**Why one directory per job:**
- Complete isolation between jobs
- Easy cleanup (rm -rf job_id)
- Easy backup (tar single job)
- Natural permission boundary for multi-user
- Simple debugging (all files together)

**Why flat (no date sharding):**
- Single VM, ~100s of jobs expected
- Filesystem handles 10K+ directories fine
- Simplicity over premature optimization
- Add sharding later if needed (yyyy/mm/job_id/)

### meta.json Schema

```json
{
  "job_id": "abc-123-def-456",
  "status": "complete",
  "email": "user@example.com",
  "molecule": {
    "smiles": "CCO",
    "name": "ethanol"
  },
  "parameters": {
    "preset": "standard",
    "basis": "6-311+G(2d,p)",
    "functional": "B3LYP",
    "solvent": "chloroform"
  },
  "timestamps": {
    "submitted": "2026-01-19T10:00:00Z",
    "started": "2026-01-19T10:00:05Z",
    "finished": "2026-01-19T10:15:30Z"
  },
  "error": null,
  "owner_id": null
}
```

### shifts.json Schema

```json
{
  "job_id": "abc-123-def-456",
  "molecule": {
    "smiles": "CCO",
    "formula": "C2H6O"
  },
  "method": {
    "basis": "6-311+G(2d,p)",
    "functional": "B3LYP",
    "solvent": "chloroform",
    "reference": "TMS"
  },
  "shifts": {
    "1H": [
      {"atom_index": 4, "atom_label": "H4", "shift_ppm": 1.23, "assignment": "CH3"},
      {"atom_index": 5, "atom_label": "H5", "shift_ppm": 1.23, "assignment": "CH3"},
      {"atom_index": 6, "atom_label": "H6", "shift_ppm": 1.23, "assignment": "CH3"},
      {"atom_index": 7, "atom_label": "H7", "shift_ppm": 3.65, "assignment": "CH2"},
      {"atom_index": 8, "atom_label": "H8", "shift_ppm": 3.65, "assignment": "CH2"},
      {"atom_index": 9, "atom_label": "H9", "shift_ppm": 2.45, "assignment": "OH"}
    ],
    "13C": [
      {"atom_index": 1, "atom_label": "C1", "shift_ppm": 18.2, "assignment": "CH3"},
      {"atom_index": 2, "atom_label": "C2", "shift_ppm": 58.1, "assignment": "CH2"}
    ]
  },
  "geometry": {
    "xyz_file": "mol.xyz",
    "optimized": true
  }
}
```

## Build Order

### Phase 1: Foundation (Build First)

Build these components first - everything else depends on them.

1. **Project structure and configuration**
   - Directory layout
   - Configuration loading (env vars, config file)
   - Data directory setup

2. **Job directory management**
   - Create job directories
   - Write/read meta.json
   - Status updates

3. **Huey setup with SQLite**
   - Configure SqliteHuey
   - Basic task definition (placeholder)
   - Consumer startup script

**Rationale:** These form the skeleton. API can't submit jobs without directories and queue.

### Phase 2: API Core (Build Second)

Build the HTTP interface once foundation exists.

4. **FastAPI application shell**
   - Basic app setup
   - Health check endpoint
   - CORS configuration

5. **Job submission endpoint**
   - POST /api/jobs
   - Input validation (SMILES)
   - Job directory creation
   - Task enqueueing

6. **Job status endpoint**
   - GET /api/jobs/{id}
   - Read meta.json
   - Return status

**Rationale:** API is the entry point for all interactions. Build it before calculation so you can test end-to-end with mock calculations.

### Phase 3: Calculation Pipeline (Build Third)

Build the actual computation once submission works.

7. **ISiCLE wrapper**
   - SMILES to 3D conversion
   - NWChem input generation
   - Calculation execution
   - Output file handling

8. **NMR parser**
   - Parse .nwo files
   - Extract shielding values
   - Convert to chemical shifts
   - Generate shifts.json

9. **Huey calculation task**
   - Wire up ISiCLE wrapper
   - Status updates during execution
   - Error handling
   - Result writing

**Rationale:** Core functionality. Must work before adding polish.

### Phase 4: Results and Notification (Build Fourth)

Build result delivery once calculations complete.

10. **Results endpoint**
    - GET /api/jobs/{id}/results
    - Return shifts.json

11. **File download endpoint**
    - GET /api/jobs/{id}/files/{name}
    - Stream raw files

12. **Email notification**
    - aiosmtplib integration
    - Completion notification
    - Error notification

**Rationale:** Users need to get results. Notification reduces polling.

### Phase 5: Web UI (Build Fifth)

Build browser interface once API is stable.

13. **Job submission page**
    - SMILES input form
    - File upload
    - Preset selection

14. **Job status page**
    - Status display
    - Auto-refresh/polling
    - Result display

15. **Results visualization**
    - Spectrum plot
    - Annotated structure
    - Download links

**Rationale:** Web UI is important but depends on stable API.

### Phase 6: Polish and Production (Build Last)

16. **Calculation presets**
    - Fast/draft mode
    - Standard mode
    - Publication quality
    - Solvation options

17. **Advanced parameters**
    - Custom basis sets
    - Custom functionals
    - Expert controls

18. **Production hardening**
    - Logging
    - Error recovery
    - Resource limits
    - Cleanup of old jobs

### Dependency Graph

```
[1. Config] --> [2. Job Dirs] --> [3. Huey] --> [4. FastAPI]
                                      |             |
                                      v             v
                                 [5. Submit] --> [6. Status]
                                      |
                                      v
                              [7. ISiCLE] --> [8. Parser] --> [9. Task]
                                                                  |
                                                                  v
                                                    [10. Results] --> [11. Files]
                                                                          |
                                                                          v
                                                                    [12. Email]
                                                                          |
                                                                          v
                                           [13. Submit UI] --> [14. Status UI] --> [15. Viz]
                                                                                       |
                                                                                       v
                                                               [16. Presets] --> [17. Advanced]
                                                                                       |
                                                                                       v
                                                                               [18. Production]
```

## Scalability Considerations

| Concern | Single VM (v1) | Multi-VM (future) |
|---------|----------------|-------------------|
| **Queue** | SqliteHuey (file-based) | Switch to RedisHuey |
| **Workers** | 1-2 workers, same VM | Multiple VMs, shared Redis |
| **Storage** | Local filesystem | Shared filesystem (NFS) or object storage |
| **API** | Single FastAPI process | Load balancer + multiple instances |
| **Database** | Filesystem meta.json | PostgreSQL for job metadata |

**v1 architecture supports growth:**
- Huey abstracts queue backend (swap SQLite for Redis later)
- Filesystem structure supports migration to shared storage
- API is stateless (can scale horizontally)
- Job ownership model ready for multi-user auth

## Sources

### Primary Sources (HIGH confidence)
- [Huey Documentation](https://huey.readthedocs.io/) - Task queue configuration and SQLite backend
- [Huey GitHub Repository](https://github.com/coleifer/huey) - SqliteHuey examples
- [FastAPI Background Tasks Documentation](https://fastapi.tiangolo.com/tutorial/background-tasks/) - When to use BackgroundTasks vs external queues
- [aiosmtplib Documentation](https://aiosmtplib.readthedocs.io/en/latest/usage.html) - Async email sending

### Supporting Sources (MEDIUM confidence)
- [Managing Background Tasks and Long-Running Operations in FastAPI](https://leapcell.io/blog/managing-background-tasks-and-long-running-operations-in-fastapi) - Job polling pattern
- [FastAPI Polling Strategy: Track Progress for Long-Running Tasks](https://openillumi.com/en/en-fastapi-long-task-progress-polling/) - Status endpoint design
- [Lightweight Django Task Queues in 2025: Beyond Celery](https://medium.com/@g.suryawanshi/lightweight-django-task-queues-in-2025-beyond-celery-74a95e0548ec) - Comparison of lightweight queue options
- [An automated framework for NMR chemical shift calculations](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-018-0305-8) - ISiCLE workflow and file structure

### Domain References
- [ISiCLE GitHub Repository](https://github.com/pnnl/isicle) - ISiCLE architecture and dependencies
- [NWChem Output Format](https://openbabel.org/docs/FileFormats/NWChem_output_format.html) - .nwo file format
- [DP5 NWChem Parser](https://github.com/Goodman-lab/DP5/blob/master/NWChem.py) - Example NMR shielding parser
