# Phase 2: Input and API - Research

**Researched:** 2026-01-19
**Domain:** FastAPI REST API + RDKit molecule validation
**Confidence:** HIGH

## Summary

Phase 2 builds the REST API layer for molecule submission and job status polling. Research confirms:

1. **FastAPI** is the correct choice for async REST APIs in Python. It provides automatic OpenAPI/Swagger documentation at `/docs` out of the box, native Pydantic v2 integration for request/response validation, and first-class support for async file uploads via `UploadFile`.

2. **RDKit molecule validation** uses `MolFromSmiles()` and `SDMolSupplier` which return `None` on invalid input rather than throwing exceptions. This enables strict pre-queue validation with specific error messages by capturing RDKit's diagnostic output.

3. **HTTP 202 Accepted** is the standard response for async job submission. The response should include a `Location` header pointing to the status endpoint and optionally a `Retry-After` header suggesting poll interval.

**Primary recommendation:** Use FastAPI with `/api/v1` prefix, resource-based endpoints (`/jobs`), RFC 7807-style error responses, and separate file upload endpoint for MOL/SDF files. Validate molecules with RDKit before queueing.

## Standard Stack

The established libraries/tools for this phase:

### Core

| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| fastapi | >=0.115.0 | REST API framework | Async, Pydantic v2, auto OpenAPI, industry standard |
| python-multipart | >=0.0.19 | Form/file upload parsing | Required by FastAPI for `UploadFile` |
| uvicorn | >=0.34.0 | ASGI server | FastAPI's recommended production server |
| rdkit | >=2025.9.3 | Molecule validation | Already installed, handles SMILES/MOL/SDF parsing |

### Supporting

| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| starlette | (via fastapi) | Request/response primitives | Status codes, background tasks |
| pydantic | >=2.5.0 | Request/response models | Already installed from Phase 1 |

### Alternatives Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| FastAPI | Flask | Flask lacks async, no built-in OpenAPI |
| FastAPI | Django REST | Too heavy for single-purpose API |
| uvicorn | gunicorn+uvicorn | Gunicorn adds process management, overkill for single VM |

**Installation:**
```bash
uv add fastapi python-multipart uvicorn
```

## Architecture Patterns

### Recommended Project Structure

```
src/
└── qm_nmr_calc/
    ├── __init__.py
    ├── api/                     # NEW: API layer
    │   ├── __init__.py
    │   ├── app.py               # FastAPI app instance
    │   ├── routers/
    │   │   ├── __init__.py
    │   │   ├── jobs.py          # /api/v1/jobs endpoints
    │   │   └── health.py        # /health endpoints
    │   ├── schemas.py           # API request/response models
    │   └── dependencies.py      # Shared dependencies
    ├── validation.py            # NEW: Molecule validation
    ├── models.py                # Existing: JobStatus, JobInput
    ├── storage.py               # Existing: Job directory management
    ├── queue.py                 # Existing: Huey instance
    ├── tasks.py                 # Existing: Calculation tasks
    ├── isicle_wrapper.py        # Existing: ISiCLE integration
    └── startup.py               # Existing: Startup validation
```

### Pattern 1: Resource-Based REST Endpoints

**What:** Endpoints organized around resources (jobs) not actions.

**When to use:** Always for REST APIs.

**Endpoints:**
```
POST   /api/v1/jobs           # Submit molecule (SMILES JSON)
POST   /api/v1/jobs/upload    # Submit molecule (MOL/SDF file)
GET    /api/v1/jobs/{job_id}  # Get job status
GET    /health                # Health check (no version prefix)
GET    /health/ready          # Readiness check
```

**Example:**
```python
# Source: FastAPI official docs - bigger applications
from fastapi import APIRouter, status
from ..schemas import JobSubmitRequest, JobStatusResponse

router = APIRouter(prefix="/jobs", tags=["jobs"])

@router.post(
    "",
    response_model=JobStatusResponse,
    status_code=status.HTTP_202_ACCEPTED,
)
async def submit_job(request: JobSubmitRequest):
    """Submit a molecule for NMR calculation."""
    # Validate, create job, queue task, return status
    pass
```

### Pattern 2: Separate File Upload Endpoint

**What:** Use different endpoint for file uploads vs JSON body.

**When to use:** When accepting both SMILES (JSON) and MOL/SDF (file upload).

**Why:** HTTP protocol limitation - cannot mix `multipart/form-data` (files) with `application/json` (body) in same request. FastAPI enforces this.

**Example:**
```python
# Source: FastAPI official docs - request files
from fastapi import APIRouter, UploadFile, Form, status
from typing import Annotated

router = APIRouter(prefix="/jobs", tags=["jobs"])

@router.post(
    "",
    response_model=JobStatusResponse,
    status_code=status.HTTP_202_ACCEPTED,
)
async def submit_smiles(request: JobSubmitRequest):
    """Submit molecule via SMILES string (JSON body)."""
    pass

@router.post(
    "/upload",
    response_model=JobStatusResponse,
    status_code=status.HTTP_202_ACCEPTED,
)
async def submit_file(
    file: UploadFile,
    name: Annotated[str | None, Form()] = None,
):
    """Submit molecule via MOL/SDF file upload."""
    # file.filename - original filename
    # file.content_type - MIME type
    # await file.read() - file contents
    pass
```

### Pattern 3: HTTP 202 with Location Header

**What:** Return job status immediately with headers pointing to status endpoint.

**When to use:** Always for async job submission.

**Example:**
```python
# Source: RFC 9110, Microsoft Azure async patterns
from fastapi import APIRouter, status
from fastapi.responses import JSONResponse

@router.post("", status_code=status.HTTP_202_ACCEPTED)
async def submit_job(request: JobSubmitRequest):
    job = create_job(request.smiles, request.name)
    queue_calculation(job.job_id)

    return JSONResponse(
        status_code=status.HTTP_202_ACCEPTED,
        content=job.model_dump(mode="json"),
        headers={
            "Location": f"/api/v1/jobs/{job.job_id}",
            "Retry-After": "30",  # Suggest 30 second poll interval
        },
    )
```

### Pattern 4: RFC 7807 Problem Details for Errors

**What:** Standardized error response format.

**When to use:** All API error responses.

**Example:**
```python
# Source: RFC 7807 / RFC 9457
from pydantic import BaseModel
from typing import Optional

class ProblemDetail(BaseModel):
    """RFC 7807 Problem Details response."""
    type: str = "about:blank"  # URI identifying error type
    title: str                 # Short human-readable summary
    status: int                # HTTP status code
    detail: Optional[str] = None  # Detailed explanation
    instance: Optional[str] = None  # URI of specific occurrence

# Example error response:
{
    "type": "https://qm-nmr-calc.example/problems/invalid-smiles",
    "title": "Invalid SMILES String",
    "status": 422,
    "detail": "Explicit valence for atom #1 O greater than permitted. Check carbon-oxygen bonds.",
    "instance": "/api/v1/jobs"
}
```

### Pattern 5: Molecule Validation with RDKit

**What:** Strict validation before queueing with specific error messages.

**When to use:** All molecule submissions.

**Example:**
```python
# Source: RDKit official docs - GettingStartedInPython
from rdkit import Chem
from rdkit import RDLogger
import io
import sys

def validate_smiles(smiles: str) -> tuple[Chem.Mol | None, str | None]:
    """
    Validate SMILES and return molecule or error message.

    Returns:
        (mol, None) on success
        (None, error_message) on failure
    """
    # Capture RDKit error messages
    stderr_capture = io.StringIO()
    RDLogger.logger().setLevel(RDLogger.ERROR)

    # Temporarily redirect stderr to capture RDKit messages
    old_stderr = sys.stderr
    sys.stderr = stderr_capture

    try:
        mol = Chem.MolFromSmiles(smiles)
    finally:
        sys.stderr = old_stderr

    if mol is None:
        error_msg = stderr_capture.getvalue().strip()
        if not error_msg:
            error_msg = "Invalid SMILES string"
        return None, error_msg

    return mol, None


def validate_mol_file(content: bytes, filename: str) -> tuple[Chem.Mol | None, str | None]:
    """
    Validate MOL/SDF file content.

    Returns:
        (mol, None) on success with single molecule
        (None, error_message) on failure or multiple molecules
    """
    content_str = content.decode('utf-8', errors='replace')

    # Check if SDF (contains $$$$ delimiter)
    if '$$$$' in content_str:
        # SDF file - check for multiple molecules
        suppl = Chem.SDMolSupplier()
        suppl.SetData(content_str)

        mols = [m for m in suppl if m is not None]

        if len(mols) == 0:
            return None, "No valid molecules found in SDF file"
        if len(mols) > 1:
            return None, f"SDF file contains {len(mols)} molecules. Only single-molecule files are accepted."

        return mols[0], None
    else:
        # MOL file
        mol = Chem.MolFromMolBlock(content_str)
        if mol is None:
            return None, "Invalid MOL file format"
        return mol, None
```

### Anti-Patterns to Avoid

- **Mixing file upload with JSON body:** HTTP doesn't support this. Use separate endpoints.
- **Validating after queueing:** Always validate molecules BEFORE adding to queue. Invalid jobs waste resources and confuse users.
- **Generic error messages:** "Invalid input" is unhelpful. Include specific RDKit error and position when possible.
- **Blocking long operations in endpoints:** Never run ISiCLE/NWChem in request handler. Use Huey queue.
- **Exposing internal errors:** Wrap exceptions, don't leak stack traces to API clients.

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| OpenAPI documentation | Custom Swagger setup | FastAPI built-in `/docs` | Zero config, auto-updates with code |
| Request validation | Manual type checking | Pydantic models | Type coercion, error messages, OpenAPI schema |
| File upload parsing | Manual multipart parsing | FastAPI UploadFile | Handles temp files, streaming, metadata |
| SMILES validation | Regex or custom parser | RDKit MolFromSmiles | Handles stereochemistry, aromaticity, valence |
| MOL/SDF parsing | Custom file parser | RDKit SDMolSupplier | Handles all format variations |
| Error response format | Custom JSON structure | RFC 7807 ProblemDetail | Standard, machine-readable, documented |

**Key insight:** FastAPI + Pydantic + RDKit handles all the hard parts. Focus on wiring them together correctly.

## Common Pitfalls

### Pitfall 1: Missing python-multipart Dependency

**What goes wrong:** `UploadFile` parameters cause cryptic import errors at runtime.

**Why it happens:** FastAPI requires `python-multipart` for form/file handling but doesn't auto-install it.

**How to avoid:**
- Always install with: `uv add fastapi python-multipart`
- Test file upload endpoint immediately after adding

**Warning signs:** "Could not import module 'multipart'" or similar import errors.

### Pitfall 2: RDKit Silent Failures

**What goes wrong:** Invalid molecules return `None` without explanation.

**Why it happens:** RDKit's `MolFromSmiles()` returns `None` for invalid input rather than raising exceptions. Error messages go to stderr, not the return value.

**How to avoid:**
- Always check `if mol is None`
- Capture stderr during RDKit calls to get error messages
- Never assume a SMILES string is valid without RDKit validation

**Warning signs:** Jobs created with invalid molecules, confusing "calculation failed" errors downstream.

### Pitfall 3: Async File Read Without Await

**What goes wrong:** File content is a coroutine object, not bytes.

**Why it happens:** `UploadFile.read()` is async in async path operations.

**How to avoid:**
```python
# Wrong
content = file.read()  # Returns coroutine!

# Correct
content = await file.read()  # Returns bytes
```

**Warning signs:** TypeError about coroutine objects.

### Pitfall 4: 404 vs 200 for Non-Existent Jobs

**What goes wrong:** Inconsistent API behavior confuses clients.

**Why it happens:** No clear decision on how to handle missing resources.

**Recommendation:** Return 404 Not Found for non-existent job IDs. This is standard REST behavior.

```python
@router.get("/{job_id}")
async def get_job_status(job_id: str):
    status = load_job_status(job_id)
    if status is None:
        raise HTTPException(
            status_code=404,
            detail=ProblemDetail(
                type="https://qm-nmr-calc.example/problems/job-not-found",
                title="Job Not Found",
                status=404,
                detail=f"No job exists with ID '{job_id}'",
            ).model_dump()
        )
    return status
```

### Pitfall 5: Forgetting Health Endpoint Authentication

**What goes wrong:** Health endpoints expose internal state or become DOS vectors.

**Why it happens:** Health checks need to be accessible but also protected.

**How to avoid:**
- `/health` (liveness) should be minimal and public
- `/health/ready` can check dependencies but keep it fast
- Don't expose sensitive details (database credentials, internal IPs)
- Consider rate limiting health endpoints

**Warning signs:** Slow health checks causing cascading failures; health endpoint revealing infrastructure details.

## Code Examples

Verified patterns from official documentation:

### Complete API Application Setup

```python
# Source: FastAPI official docs
# src/qm_nmr_calc/api/app.py

from fastapi import FastAPI
from .routers import jobs, health

app = FastAPI(
    title="QM NMR Calculator API",
    description="Asynchronous NMR quantum mechanical calculations",
    version="0.1.0",
    docs_url="/docs",
    redoc_url="/redoc",
    openapi_url="/api/v1/openapi.json",
)

# Health endpoints at root (no version prefix)
app.include_router(health.router)

# API v1 endpoints
app.include_router(jobs.router, prefix="/api/v1")
```

### Job Submission Request Schema

```python
# Source: CONTEXT.md decisions + Pydantic v2
# src/qm_nmr_calc/api/schemas.py

from pydantic import BaseModel, Field
from typing import Optional

class JobSubmitRequest(BaseModel):
    """Request body for SMILES job submission."""
    smiles: str = Field(
        ...,
        description="SMILES string of molecule to calculate",
        examples=["CCO", "c1ccccc1"],
    )
    name: Optional[str] = Field(
        None,
        description="Optional name/label for the molecule",
        examples=["Ethanol", "Benzene"],
        max_length=100,
    )

class JobStatusResponse(BaseModel):
    """Response for job status queries."""
    job_id: str = Field(..., description="12-character hex job ID")
    status: str = Field(..., description="Job status: queued, running, complete, failed")
    created_at: str = Field(..., description="ISO 8601 timestamp")
    started_at: Optional[str] = Field(None, description="When job started running")
    completed_at: Optional[str] = Field(None, description="When job completed")
    input_smiles: str = Field(..., description="Original SMILES input")
    input_name: Optional[str] = Field(None, description="User-provided molecule name")
    error_message: Optional[str] = Field(None, description="Error message if failed")
```

### Jobs Router Implementation

```python
# Source: FastAPI + Phase 1 integration
# src/qm_nmr_calc/api/routers/jobs.py

from fastapi import APIRouter, HTTPException, UploadFile, Form, status
from fastapi.responses import JSONResponse
from typing import Annotated

from ..schemas import JobSubmitRequest, JobStatusResponse, ProblemDetail
from ...validation import validate_smiles, validate_mol_file
from ...storage import create_job_directory, load_job_status
from ...isicle_wrapper import get_versions
from ...queue import huey
from ...tasks import run_optimization

router = APIRouter(prefix="/jobs", tags=["jobs"])

@router.post(
    "",
    response_model=JobStatusResponse,
    status_code=status.HTTP_202_ACCEPTED,
    responses={
        202: {"description": "Job accepted for processing"},
        422: {"model": ProblemDetail, "description": "Invalid molecule"},
    },
)
async def submit_smiles(request: JobSubmitRequest):
    """Submit molecule via SMILES string."""
    # Validate SMILES
    mol, error = validate_smiles(request.smiles)
    if mol is None:
        raise HTTPException(
            status_code=422,
            detail={
                "type": "https://qm-nmr-calc.example/problems/invalid-smiles",
                "title": "Invalid SMILES String",
                "status": 422,
                "detail": error,
            },
        )

    # Create job
    versions = get_versions()
    job = create_job_directory(
        smiles=request.smiles,
        name=request.name,
        isicle_version=versions.isicle,
        nwchem_version=versions.nwchem,
    )

    # Queue calculation
    run_optimization(job.job_id)

    return JSONResponse(
        status_code=status.HTTP_202_ACCEPTED,
        content=JobStatusResponse.model_validate(job).model_dump(mode="json"),
        headers={
            "Location": f"/api/v1/jobs/{job.job_id}",
            "Retry-After": "30",
        },
    )

@router.post(
    "/upload",
    response_model=JobStatusResponse,
    status_code=status.HTTP_202_ACCEPTED,
)
async def submit_file(
    file: UploadFile,
    name: Annotated[str | None, Form(description="Optional molecule name")] = None,
):
    """Submit molecule via MOL/SDF file upload."""
    # Validate file type
    allowed_types = {
        "chemical/x-mdl-molfile",
        "chemical/x-mdl-sdfile",
        "application/octet-stream",  # Common fallback
    }
    # Also check by extension
    filename = file.filename or ""
    if not (filename.endswith(".mol") or filename.endswith(".sdf")):
        if file.content_type not in allowed_types:
            raise HTTPException(
                status_code=422,
                detail={
                    "type": "https://qm-nmr-calc.example/problems/invalid-file-type",
                    "title": "Invalid File Type",
                    "status": 422,
                    "detail": f"Expected .mol or .sdf file, got '{filename}'",
                },
            )

    # Read and validate content
    content = await file.read()
    mol, error = validate_mol_file(content, filename)
    if mol is None:
        raise HTTPException(
            status_code=422,
            detail={
                "type": "https://qm-nmr-calc.example/problems/invalid-molecule-file",
                "title": "Invalid Molecule File",
                "status": 422,
                "detail": error,
            },
        )

    # Convert to SMILES for storage
    from rdkit import Chem
    smiles = Chem.MolToSmiles(mol)

    # Create job (same as SMILES submission)
    versions = get_versions()
    job = create_job_directory(
        smiles=smiles,
        name=name or filename,
        isicle_version=versions.isicle,
        nwchem_version=versions.nwchem,
    )

    run_optimization(job.job_id)

    return JSONResponse(
        status_code=status.HTTP_202_ACCEPTED,
        content=JobStatusResponse.model_validate(job).model_dump(mode="json"),
        headers={
            "Location": f"/api/v1/jobs/{job.job_id}",
            "Retry-After": "30",
        },
    )

@router.get(
    "/{job_id}",
    response_model=JobStatusResponse,
    responses={
        200: {"description": "Job status"},
        404: {"model": ProblemDetail, "description": "Job not found"},
    },
)
async def get_job_status(job_id: str):
    """Get status of a calculation job."""
    status_obj = load_job_status(job_id)
    if status_obj is None:
        raise HTTPException(
            status_code=404,
            detail={
                "type": "https://qm-nmr-calc.example/problems/job-not-found",
                "title": "Job Not Found",
                "status": 404,
                "detail": f"No job exists with ID '{job_id}'",
            },
        )
    return status_obj
```

### Health Router Implementation

```python
# Source: FastAPI health check best practices
# src/qm_nmr_calc/api/routers/health.py

from fastapi import APIRouter
from fastapi.responses import JSONResponse
from datetime import datetime

router = APIRouter(tags=["health"])

@router.get("/health")
async def liveness():
    """Liveness probe - is the service running?"""
    return {"status": "alive"}

@router.get("/health/ready")
async def readiness():
    """Readiness probe - can the service handle requests?"""
    # Check dependencies
    checks = {}

    # Check data directory writable
    from pathlib import Path
    data_dir = Path("./data/jobs")
    try:
        data_dir.mkdir(parents=True, exist_ok=True)
        test_file = data_dir / ".health_check"
        test_file.write_text("ok")
        test_file.unlink()
        checks["data_directory"] = "ok"
    except Exception as e:
        checks["data_directory"] = f"error: {e}"
        return JSONResponse(
            status_code=503,
            content={"status": "not ready", "checks": checks},
        )

    # Check Huey queue accessible
    try:
        from ...queue import huey
        # Just verify import works
        checks["task_queue"] = "ok"
    except Exception as e:
        checks["task_queue"] = f"error: {e}"
        return JSONResponse(
            status_code=503,
            content={"status": "not ready", "checks": checks},
        )

    return {
        "status": "ready",
        "checks": checks,
        "timestamp": datetime.utcnow().isoformat() + "Z",
    }
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Flask + flask-restful | FastAPI | 2019-2024 | Async, type hints, auto-docs |
| Pydantic v1 | Pydantic v2 | 2023 | 5-50x faster, better errors |
| Manual OpenAPI | FastAPI auto-generation | 2019+ | Zero config documentation |
| Custom error formats | RFC 7807/9457 | 2016/2023 | Standardized machine-readable errors |
| Sync file upload | Async UploadFile | 2020+ | Memory efficient, streaming |

**Deprecated/outdated:**
- Flask for new API projects: Use FastAPI for better DX and performance
- Pydantic v1 syntax (`@validator`): Use v2 syntax (`@field_validator`)
- RFC 7807: Superseded by RFC 9457 (compatible, minor improvements)

## Open Questions

Things that couldn't be fully resolved:

1. **RDKit stderr capture reliability**
   - What we know: RDKit writes errors to stderr, not return values
   - What's unclear: Whether all error messages are captured across all failure modes
   - Recommendation: Test with known-bad molecules; add fallback generic message

2. **File size limits for uploads**
   - What we know: FastAPI/Starlette use `SpooledTemporaryFile` (memory up to limit, then disk)
   - What's unclear: Optimal limit for MOL/SDF files
   - Recommendation: Start with 1MB limit, adjust based on real-world usage

3. **Multi-molecule SDF detection edge cases**
   - What we know: SDF uses `$$$$` as molecule delimiter
   - What's unclear: Whether all SDF variations use this consistently
   - Recommendation: Count molecules via SDMolSupplier iteration, not string parsing

## Sources

### Primary (HIGH confidence)

- **FastAPI Official Documentation** - https://fastapi.tiangolo.com/
  - [Request Files](https://fastapi.tiangolo.com/tutorial/request-files/) - UploadFile API
  - [Bigger Applications](https://fastapi.tiangolo.com/tutorial/bigger-applications/) - Router patterns
  - [Metadata and Docs URLs](https://fastapi.tiangolo.com/tutorial/metadata/) - OpenAPI configuration

- **RDKit Official Documentation** - https://www.rdkit.org/docs/
  - [Getting Started](https://www.rdkit.org/docs/GettingStartedInPython.html) - MolFromSmiles, SDMolSupplier
  - [rdmolfiles module](https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html) - File parsing API

- **RFC 9457** (obsoletes RFC 7807) - https://www.rfc-editor.org/rfc/rfc9457.html - Problem Details for HTTP APIs

### Secondary (MEDIUM confidence)

- **Microsoft Azure Async Request-Reply Pattern** - https://learn.microsoft.com/en-us/azure/architecture/patterns/async-request-reply - Location header, polling patterns

- **fastapi-rfc7807** GitHub - https://github.com/vapor-ware/fastapi-rfc7807 - RFC 7807 implementation reference

### Tertiary (LOW confidence)

- WebSearch results for FastAPI versioning patterns - Cross-verified with official docs
- WebSearch results for health check best practices - Based on Kubernetes patterns

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - FastAPI is dominant Python API framework, well-documented
- Architecture: HIGH - Based on official FastAPI patterns and RFC standards
- Molecule validation: HIGH - RDKit behavior verified in official docs
- Pitfalls: MEDIUM - Based on common issues, some edge cases unverified

**Research date:** 2026-01-19
**Valid until:** 2026-02-19 (stable libraries, 30 days)

---

*Phase: 02-input-and-api*
*Research completed: 2026-01-19*
