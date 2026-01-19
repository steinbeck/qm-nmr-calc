"""Pydantic models for job status and input validation."""

from datetime import datetime
from typing import Literal, Optional

from pydantic import BaseModel, ConfigDict


class JobInput(BaseModel):
    """Input parameters for a calculation job."""

    model_config = ConfigDict(strict=True)

    smiles: str
    name: Optional[str] = None  # User-provided molecule name/label


class JobStatus(BaseModel):
    """Complete status of a calculation job."""

    model_config = ConfigDict(
        # Use 'iso8601' serialization for datetime fields
        ser_json_timedelta="iso8601",
    )

    job_id: str
    status: Literal["queued", "running", "complete", "failed"]
    created_at: datetime
    started_at: Optional[datetime] = None
    completed_at: Optional[datetime] = None

    # Input
    input: JobInput

    # Versions (for reproducibility)
    isicle_version: str
    nwchem_version: str

    # Error info (if failed)
    error_message: Optional[str] = None
    error_traceback: Optional[str] = None

    # Resource usage (filled on completion)
    cpu_time_seconds: Optional[float] = None
    memory_peak_mb: Optional[float] = None
