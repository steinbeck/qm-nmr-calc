# Phase 4: Results Delivery - Research

**Researched:** 2026-01-19
**Domain:** FastAPI file downloads + async email notifications
**Confidence:** HIGH

## Summary

Phase 4 delivers calculation results through multiple channels: JSON API for structured NMR data, file downloads for geometry and raw output, and email notifications on job completion. Research confirms:

1. **FastAPI FileResponse** is the standard way to serve file downloads. It streams files efficiently without loading them entirely into memory, automatically sets Content-Length/Last-Modified/ETag headers, and supports custom filenames for the Content-Disposition header.

2. **HTTP status codes for results** should follow REST conventions: 200 OK with data when complete, 404 Not Found for non-existent jobs, and 409 Conflict when results are requested before job completion (the resource exists but is not in the right state for the requested operation).

3. **aiosmtplib** provides async email sending that integrates well with FastAPI. Python's stdlib `email.message.EmailMessage` with `add_alternative()` creates proper multipart messages with both plain text and HTML versions.

**Primary recommendation:** Add dedicated endpoints for each result type (`/jobs/{job_id}/results`, `/jobs/{job_id}/geometry`, `/jobs/{job_id}/output`), use FileResponse for downloads with appropriate MIME types, implement optional email field in submission request with notification on Huey task completion.

## Standard Stack

The established libraries/tools for this phase:

### Core

| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| fastapi.responses.FileResponse | (built-in) | Serve file downloads | Async streaming, auto headers, Content-Disposition |
| aiosmtplib | >=3.0.0 | Async SMTP client | Non-blocking email in FastAPI async context |
| email.message.EmailMessage | (stdlib) | Build email messages | Standard library, supports multipart/alternative |

### Supporting

| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| email.mime.multipart | (stdlib) | Multipart messages | HTML + plain text emails |
| pathlib | (stdlib) | Path handling | File existence checks, secure path construction |

### Alternatives Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| aiosmtplib | smtplib | smtplib is synchronous, blocks event loop |
| FileResponse | StreamingResponse | StreamingResponse is more complex, FileResponse handles files better |
| email.message | third-party email lib | stdlib sufficient, no extra dependencies |

**Installation:**
```bash
uv add aiosmtplib
```

## Architecture Patterns

### Recommended Project Structure

```
src/
└── qm_nmr_calc/
    ├── api/
    │   ├── routers/
    │   │   ├── jobs.py           # Add results/download endpoints
    │   │   └── health.py
    │   ├── schemas.py            # Existing + no new schemas needed
    │   └── app.py
    ├── notifications.py          # NEW: Email notification logic
    ├── storage.py                # Add file path helpers
    ├── tasks.py                  # Add email notification on completion
    └── models.py                 # Add email field to JobInput
```

### Pattern 1: Dedicated Result Endpoints

**What:** Separate endpoints for each result type rather than query parameters.

**When to use:** When result types have different content types and behaviors.

**Endpoints:**
```
GET /api/v1/jobs/{job_id}                 # Existing: full status with NMR results
GET /api/v1/jobs/{job_id}/results         # NEW: NMR results only (JSON)
GET /api/v1/jobs/{job_id}/geometry        # NEW: Download optimized XYZ
GET /api/v1/jobs/{job_id}/geometry.sdf    # NEW: Download as SDF format
GET /api/v1/jobs/{job_id}/output          # NEW: Download raw NWChem output (zip)
```

**Example:**
```python
# Source: FastAPI custom responses + project patterns
from fastapi import APIRouter, HTTPException, status
from fastapi.responses import FileResponse
from pathlib import Path

router = APIRouter(prefix="/jobs", tags=["jobs"])

@router.get(
    "/{job_id}/geometry",
    response_class=FileResponse,
    responses={
        200: {"description": "Optimized molecular geometry (XYZ format)"},
        404: {"model": ProblemDetail, "description": "Job not found"},
        409: {"model": ProblemDetail, "description": "Job not complete"},
    },
)
async def download_geometry(job_id: str):
    """Download optimized molecular geometry as XYZ file."""
    job_status = load_job_status(job_id)
    if job_status is None:
        raise HTTPException(status_code=404, detail={...})

    if job_status.status != "complete":
        raise HTTPException(
            status_code=409,
            detail={
                "type": "https://qm-nmr-calc.example/problems/job-not-complete",
                "title": "Job Not Complete",
                "status": 409,
                "detail": f"Job '{job_id}' is {job_status.status}. Results not available.",
            },
        )

    geometry_path = Path(job_status.optimized_geometry_file)
    if not geometry_path.exists():
        raise HTTPException(status_code=404, detail={...})

    return FileResponse(
        path=geometry_path,
        filename=f"{job_id}_optimized.xyz",
        media_type="chemical/x-xyz",
    )
```

### Pattern 2: FileResponse with Appropriate MIME Types

**What:** Use correct chemical MIME types for molecule files.

**When to use:** All file downloads.

**MIME types:**
```python
MIME_TYPES = {
    ".xyz": "chemical/x-xyz",
    ".sdf": "chemical/x-mdl-sdfile",
    ".mol": "chemical/x-mdl-molfile",
    ".out": "text/plain",  # NWChem output
    ".nw": "text/plain",   # NWChem input
    ".zip": "application/zip",
}
```

**Example:**
```python
# Source: FastAPI FileResponse docs
from fastapi.responses import FileResponse

return FileResponse(
    path="/path/to/file.xyz",
    filename="molecule_optimized.xyz",  # Suggested download name
    media_type="chemical/x-xyz",
    headers={"Content-Disposition": "attachment"},  # Force download
)
```

### Pattern 3: HTTP 409 Conflict for Not Ready

**What:** Return 409 when resource exists but operation cannot be performed due to current state.

**When to use:** When requesting results from incomplete jobs.

**Why 409:** The job resource exists (not 404), and the request is valid (not 400), but the current state prevents the operation. RFC 9110: "indicates a request conflict with the current state of the target resource."

**Example:**
```python
# Source: RFC 9110 Section 15.5.10
if job_status.status in ("queued", "running"):
    raise HTTPException(
        status_code=409,
        detail={
            "type": "https://qm-nmr-calc.example/problems/job-not-complete",
            "title": "Job Not Complete",
            "status": 409,
            "detail": f"Job is currently '{job_status.status}'. Please poll until status is 'complete'.",
            "instance": f"/api/v1/jobs/{job_id}",
        },
    )
```

### Pattern 4: Async Email on Task Completion

**What:** Send email notification from Huey task after job completes.

**When to use:** When user provided notification email at submission.

**Integration point:** Huey signal handlers (already in queue.py).

**Example:**
```python
# Source: aiosmtplib docs + Python email.message
import asyncio
from email.message import EmailMessage
import aiosmtplib

async def send_completion_email(
    to_email: str,
    job_id: str,
    status: str,
    base_url: str = "http://localhost:8000",
):
    """Send job completion notification email."""
    msg = EmailMessage()
    msg["From"] = "noreply@qm-nmr-calc.example"
    msg["To"] = to_email
    msg["Subject"] = f"NMR Calculation {status.title()}: {job_id}"

    results_url = f"{base_url}/api/v1/jobs/{job_id}"

    # Plain text version
    msg.set_content(f"""\
Your NMR calculation has {status}.

Job ID: {job_id}
Status: {status}

View results: {results_url}
""")

    # HTML version
    msg.add_alternative(f"""\
<html>
<body>
<p>Your NMR calculation has <strong>{status}</strong>.</p>
<p><strong>Job ID:</strong> {job_id}<br>
<strong>Status:</strong> {status}</p>
<p><a href="{results_url}">View Results</a></p>
</body>
</html>
""", subtype='html')

    await aiosmtplib.send(
        msg,
        hostname=SMTP_HOST,
        port=SMTP_PORT,
        username=SMTP_USER,
        password=SMTP_PASSWORD,
        start_tls=True,
    )


# In Huey signal handler (sync context), use asyncio.run()
def on_task_complete(signal, task, *args):
    job_status = load_job_status(task.args[0])
    if job_status.input.notification_email:
        asyncio.run(send_completion_email(
            job_status.input.notification_email,
            job_status.job_id,
            job_status.status,
        ))
```

### Pattern 5: Optional Email Field in Submission

**What:** Add optional email field to job submission for opt-in notifications.

**When to use:** When implementing NOTF-01.

**Example:**
```python
# Source: Pydantic EmailStr type
from pydantic import BaseModel, EmailStr, Field
from typing import Optional

class JobSubmitRequest(BaseModel):
    smiles: str = Field(...)
    name: Optional[str] = Field(None)
    preset: Literal["draft", "production"] = Field(default="production")
    solvent: str = Field(...)
    notification_email: Optional[EmailStr] = Field(
        None,
        description="Email address for completion notification (opt-in)",
        examples=["user@example.com"],
    )
```

### Anti-Patterns to Avoid

- **Serving files from user-controlled paths:** Always construct paths server-side using job_id, never accept file paths from client.
- **Blocking email in request handlers:** Use Huey signal or BackgroundTasks, not inline in endpoint.
- **Missing file existence checks:** Always verify file exists before FileResponse; return 404 if missing.
- **Exposing scratch directory contents:** Only serve specific output files, not arbitrary scratch files.
- **Returning 404 for incomplete jobs:** Use 409 Conflict to distinguish "job doesn't exist" from "job not ready".

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| File download headers | Manual Content-Disposition | FileResponse filename param | Handles encoding, browser quirks |
| Async email | Threading + smtplib | aiosmtplib | Non-blocking, event loop friendly |
| HTML + text email | Manual MIME assembly | EmailMessage.add_alternative() | Correct multipart structure |
| Email validation | Regex patterns | Pydantic EmailStr | Comprehensive validation, standard |
| Path traversal prevention | Manual ../ checking | Construct paths from job_id only | No user-controlled paths |

**Key insight:** FileResponse and EmailMessage handle the complexity. Focus on business logic (when to send, what to include), not transport mechanics.

## Common Pitfalls

### Pitfall 1: Path Traversal in File Downloads

**What goes wrong:** Attacker requests `../../../etc/passwd` via job_id.

**Why it happens:** Directly using user input in file paths.

**How to avoid:**
```python
# WRONG - vulnerable to path traversal
file_path = f"./data/jobs/{job_id}/{filename}"

# RIGHT - construct path from known components only
job_dir = get_job_dir(job_id)  # Uses DATA_DIR / job_id
if not job_dir.exists():
    raise HTTPException(status_code=404, ...)
geometry_file = job_dir / "output" / "optimized.xyz"
```

**Warning signs:** File paths contain user-provided values beyond validated job_id.

### Pitfall 2: Blocking Event Loop with Email

**What goes wrong:** Email sending blocks all other requests.

**Why it happens:** Using synchronous smtplib in async endpoint.

**How to avoid:**
- Use aiosmtplib for async email
- Send email from Huey task (separate process), not from API request
- If must send from API, use BackgroundTasks

**Warning signs:** API latency spikes when emails are sent.

### Pitfall 3: Missing Files for Completed Jobs

**What goes wrong:** Job shows complete but geometry file missing.

**Why it happens:** Calculation completed but file save failed; or file was cleaned up.

**How to avoid:**
- Check file exists before returning FileResponse
- Return 404 with specific error message
- Log missing files for debugging

```python
if not geometry_path.exists():
    raise HTTPException(
        status_code=404,
        detail={
            "type": "https://qm-nmr-calc.example/problems/file-not-found",
            "title": "Output File Not Found",
            "status": 404,
            "detail": f"Geometry file for job '{job_id}' is missing. This may indicate a calculation error.",
        },
    )
```

**Warning signs:** Users report 500 errors on download endpoints.

### Pitfall 4: Email Validation at Wrong Time

**What goes wrong:** Invalid email stored, then email send fails later.

**Why it happens:** Email validated only at send time, not submission time.

**How to avoid:**
- Use Pydantic EmailStr in request schema
- Validation happens at submission, before job is queued
- Invalid email returns 422 immediately

**Warning signs:** Jobs complete but notifications fail silently.

### Pitfall 5: Hardcoded Base URL in Emails

**What goes wrong:** Email links point to localhost or wrong domain.

**Why it happens:** Base URL hardcoded during development.

**How to avoid:**
- Configure BASE_URL via environment variable
- Use Request.url_for() if sending from API context
- Test email links in staging environment

```python
# config.py
from pydantic_settings import BaseSettings

class Settings(BaseSettings):
    base_url: str = "http://localhost:8000"
    # In production: BASE_URL=https://qm-nmr-calc.example
```

**Warning signs:** Emails contain localhost URLs in production.

## Code Examples

Verified patterns from official documentation:

### NMR Results Endpoint

```python
# Source: Existing project patterns
from fastapi import APIRouter, HTTPException, status
from ..schemas import NMRResultsResponse, ProblemDetail
from ...storage import load_job_status

router = APIRouter(prefix="/jobs", tags=["jobs"])

@router.get(
    "/{job_id}/results",
    response_model=NMRResultsResponse,
    responses={
        200: {"description": "NMR calculation results"},
        404: {"model": ProblemDetail, "description": "Job not found"},
        409: {"model": ProblemDetail, "description": "Job not complete"},
    },
)
async def get_nmr_results(job_id: str):
    """Get NMR chemical shifts for a completed job.

    Returns 1H and 13C chemical shifts with atom indices,
    plus calculation metadata (functional, basis set, solvent).
    """
    job_status = load_job_status(job_id)
    if job_status is None:
        raise HTTPException(
            status_code=404,
            detail={
                "type": "https://qm-nmr-calc.example/problems/job-not-found",
                "title": "Job Not Found",
                "status": 404,
                "detail": f"No job exists with ID '{job_id}'",
            },
        )

    if job_status.status != "complete":
        raise HTTPException(
            status_code=409,
            detail={
                "type": "https://qm-nmr-calc.example/problems/job-not-complete",
                "title": "Job Not Complete",
                "status": 409,
                "detail": f"Job is '{job_status.status}'. Poll /api/v1/jobs/{job_id} until complete.",
            },
        )

    if job_status.nmr_results is None:
        raise HTTPException(
            status_code=404,
            detail={
                "type": "https://qm-nmr-calc.example/problems/no-results",
                "title": "Results Not Available",
                "status": 404,
                "detail": f"Job '{job_id}' completed but no NMR results were stored.",
            },
        )

    return {
        "h1_shifts": [
            {"index": s.index, "atom": s.atom, "shift": s.shift}
            for s in job_status.nmr_results.h1_shifts
        ],
        "c13_shifts": [
            {"index": s.index, "atom": s.atom, "shift": s.shift}
            for s in job_status.nmr_results.c13_shifts
        ],
        "functional": job_status.nmr_results.functional,
        "basis_set": job_status.nmr_results.basis_set,
        "solvent": job_status.nmr_results.solvent,
    }
```

### Geometry Download Endpoint

```python
# Source: FastAPI FileResponse docs
from fastapi import APIRouter, HTTPException
from fastapi.responses import FileResponse
from pathlib import Path
from ...storage import load_job_status, get_job_dir

@router.get(
    "/{job_id}/geometry",
    response_class=FileResponse,
    responses={
        200: {
            "description": "Optimized molecular geometry",
            "content": {"chemical/x-xyz": {}},
        },
        404: {"model": ProblemDetail},
        409: {"model": ProblemDetail},
    },
)
async def download_geometry(job_id: str):
    """Download optimized molecular geometry as XYZ file.

    The XYZ file contains the DFT-optimized 3D coordinates
    that were used for the NMR shielding calculation.
    """
    job_status = load_job_status(job_id)
    if job_status is None:
        raise HTTPException(status_code=404, detail={...})

    if job_status.status != "complete":
        raise HTTPException(status_code=409, detail={...})

    # Construct path securely (no user input in path)
    job_dir = get_job_dir(job_id)
    geometry_file = job_dir / "output" / "optimized.xyz"

    if not geometry_file.exists():
        raise HTTPException(
            status_code=404,
            detail={
                "type": "https://qm-nmr-calc.example/problems/file-not-found",
                "title": "Geometry File Not Found",
                "status": 404,
                "detail": f"Optimized geometry file not found for job '{job_id}'.",
            },
        )

    return FileResponse(
        path=geometry_file,
        filename=f"{job_id}_optimized.xyz",
        media_type="chemical/x-xyz",
    )
```

### Email Notification Module

```python
# Source: aiosmtplib docs + Python email.message docs
# src/qm_nmr_calc/notifications.py

import asyncio
from email.message import EmailMessage
from typing import Optional
import aiosmtplib
import os

# Configuration from environment
SMTP_HOST = os.getenv("SMTP_HOST", "localhost")
SMTP_PORT = int(os.getenv("SMTP_PORT", "587"))
SMTP_USER = os.getenv("SMTP_USER", "")
SMTP_PASSWORD = os.getenv("SMTP_PASSWORD", "")
NOTIFICATION_FROM = os.getenv("NOTIFICATION_FROM", "noreply@qm-nmr-calc.example")
BASE_URL = os.getenv("BASE_URL", "http://localhost:8000")


async def send_job_notification(
    to_email: str,
    job_id: str,
    status: str,
    error_message: Optional[str] = None,
) -> bool:
    """Send job completion/failure notification email.

    Args:
        to_email: Recipient email address
        job_id: Job identifier
        status: Job status (complete, failed)
        error_message: Error details if failed

    Returns:
        True if sent successfully, False otherwise
    """
    msg = EmailMessage()
    msg["From"] = NOTIFICATION_FROM
    msg["To"] = to_email
    msg["Subject"] = f"NMR Calculation {status.title()}: {job_id}"

    results_url = f"{BASE_URL}/api/v1/jobs/{job_id}"

    # Plain text version
    if status == "complete":
        plain_body = f"""\
Your NMR calculation has completed successfully.

Job ID: {job_id}
Status: Complete

View your results: {results_url}

Available downloads:
- Chemical shifts (JSON): {results_url}/results
- Optimized geometry (XYZ): {results_url}/geometry
- Raw NWChem output: {results_url}/output
"""
    else:
        plain_body = f"""\
Your NMR calculation has failed.

Job ID: {job_id}
Status: Failed
Error: {error_message or "Unknown error"}

View details: {results_url}
"""

    msg.set_content(plain_body)

    # HTML version
    if status == "complete":
        html_body = f"""\
<html>
<body style="font-family: Arial, sans-serif; line-height: 1.6;">
<h2>NMR Calculation Complete</h2>
<p>Your NMR calculation has completed successfully.</p>
<table style="border-collapse: collapse;">
<tr><td><strong>Job ID:</strong></td><td>{job_id}</td></tr>
<tr><td><strong>Status:</strong></td><td style="color: green;">Complete</td></tr>
</table>
<h3>Downloads</h3>
<ul>
<li><a href="{results_url}/results">Chemical Shifts (JSON)</a></li>
<li><a href="{results_url}/geometry">Optimized Geometry (XYZ)</a></li>
<li><a href="{results_url}/output">Raw NWChem Output</a></li>
</ul>
<p><a href="{results_url}" style="background-color: #4CAF50; color: white; padding: 10px 20px; text-decoration: none; border-radius: 5px;">View Results</a></p>
</body>
</html>
"""
    else:
        html_body = f"""\
<html>
<body style="font-family: Arial, sans-serif; line-height: 1.6;">
<h2>NMR Calculation Failed</h2>
<p>Your NMR calculation encountered an error.</p>
<table style="border-collapse: collapse;">
<tr><td><strong>Job ID:</strong></td><td>{job_id}</td></tr>
<tr><td><strong>Status:</strong></td><td style="color: red;">Failed</td></tr>
<tr><td><strong>Error:</strong></td><td>{error_message or "Unknown error"}</td></tr>
</table>
<p><a href="{results_url}">View Details</a></p>
</body>
</html>
"""

    msg.add_alternative(html_body, subtype='html')

    try:
        await aiosmtplib.send(
            msg,
            hostname=SMTP_HOST,
            port=SMTP_PORT,
            username=SMTP_USER if SMTP_USER else None,
            password=SMTP_PASSWORD if SMTP_PASSWORD else None,
            start_tls=bool(SMTP_USER),  # Use TLS if authenticating
        )
        return True
    except Exception as e:
        # Log error but don't fail the job
        print(f"Failed to send notification email: {e}")
        return False


def send_job_notification_sync(
    to_email: str,
    job_id: str,
    status: str,
    error_message: Optional[str] = None,
) -> bool:
    """Synchronous wrapper for use in Huey tasks."""
    return asyncio.run(send_job_notification(
        to_email, job_id, status, error_message
    ))
```

### Integration in Huey Signal Handler

```python
# In queue.py - extend existing signal handler
from huey.signals import SIGNAL_COMPLETE, SIGNAL_ERROR
from .storage import load_job_status
from .notifications import send_job_notification_sync

@huey.signal(SIGNAL_COMPLETE)
def on_task_complete(signal, task, *args, **kwargs):
    """Handle task completion - update status and send notification."""
    if task.name == "run_nmr_task":
        job_id = task.args[0]
        job_status = load_job_status(job_id)

        # Send notification if email was provided
        if job_status and job_status.input.notification_email:
            send_job_notification_sync(
                to_email=job_status.input.notification_email,
                job_id=job_id,
                status="complete",
            )

@huey.signal(SIGNAL_ERROR)
def on_task_error(signal, task, exc, *args, **kwargs):
    """Handle task error - send failure notification."""
    if task.name == "run_nmr_task":
        job_id = task.args[0]
        job_status = load_job_status(job_id)

        if job_status and job_status.input.notification_email:
            send_job_notification_sync(
                to_email=job_status.input.notification_email,
                job_id=job_id,
                status="failed",
                error_message=str(exc),
            )
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Sync smtplib | aiosmtplib | 2020+ | Non-blocking email in async apps |
| Manual MIME | EmailMessage.add_alternative() | Python 3.6+ | Cleaner multipart message construction |
| 503 for "not ready" | 409 Conflict | Best practice | Clearer semantics (resource state vs service unavailable) |
| Query params for format | Dedicated endpoints | REST best practice | Cleaner URLs, better caching |

**Deprecated/outdated:**
- Using smtplib directly in async code: Blocks event loop, use aiosmtplib
- email.MIMEText/MIMEMultipart directly: EmailMessage is cleaner for most cases
- Returning 202 for result polling: 202 is for initial submission, use 200/409 for polling

## Open Questions

Things that couldn't be fully resolved:

1. **Raw NWChem output location**
   - What we know: ISiCLE writes to scratch_dir inside job_dir
   - What's unclear: Exact file names and what's worth keeping after calculation
   - Recommendation: Check actual scratch contents after running a job; likely want `.out` files

2. **SDF conversion from XYZ**
   - What we know: RDKit can convert formats; XYZ doesn't have connectivity
   - What's unclear: Whether to regenerate SDF from original SMILES or convert from XYZ
   - Recommendation: Generate SDF from original mol + optimized coordinates if needed

3. **Email delivery reliability**
   - What we know: aiosmtplib sends emails asynchronously
   - What's unclear: Retry policy, bounce handling, delivery tracking
   - Recommendation: Log failures, don't retry automatically for v1, consider email service later

## Sources

### Primary (HIGH confidence)

- [FastAPI Custom Response](https://fastapi.tiangolo.com/advanced/custom-response/) - FileResponse documentation
- [aiosmtplib Documentation](https://aiosmtplib.readthedocs.io/en/latest/usage.html) - Async SMTP client usage
- [Python email.examples](https://docs.python.org/3/library/email.examples.html) - EmailMessage with add_alternative
- [RFC 9110 Section 15.5.10](https://www.rfc-editor.org/rfc/rfc9110.html#section-15.5.10) - HTTP 409 Conflict

### Secondary (MEDIUM confidence)

- [Async Request-Reply Pattern](https://learn.microsoft.com/en-us/azure/architecture/patterns/async-request-reply) - 202/303 polling pattern
- [REST API Long-Running Tasks](https://restfulapi.net/rest-api-design-for-long-running-tasks/) - Status code conventions
- [XYZ File Format](https://en.wikipedia.org/wiki/XYZ_file_format) - MIME type chemical/x-xyz
- [Chemical Table File](https://en.wikipedia.org/wiki/Chemical_table_file) - SDF MIME type chemical/x-mdl-sdfile

### Tertiary (LOW confidence)

- WebSearch results for NWChem output files - Needs verification with actual calculation
- WebSearch results for FastAPI + email patterns - Verified against official docs

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - FastAPI FileResponse and aiosmtplib well-documented
- Architecture: HIGH - Following established REST patterns and existing project structure
- File downloads: HIGH - FileResponse behavior verified in official docs
- Email: MEDIUM - aiosmtplib API verified, integration patterns need testing
- Pitfalls: MEDIUM - Based on common issues, security patterns well-established

**Research date:** 2026-01-19
**Valid until:** 2026-02-19 (stable libraries, 30 days)

---

*Phase: 04-results-delivery*
*Research completed: 2026-01-19*
