"""Pydantic models for job status and input validation."""

from datetime import datetime
from typing import Literal, Optional

from pydantic import BaseModel, ConfigDict


class AtomShift(BaseModel):
    """NMR chemical shift for a single atom."""

    model_config = ConfigDict(strict=True)

    index: int  # NWChem 1-based atom index
    atom: str  # Element symbol: 'H' or 'C'
    shielding: float  # Raw isotropic shielding in ppm
    shift: float  # Chemical shift relative to TMS in ppm


class NMRResults(BaseModel):
    """NMR calculation results."""

    model_config = ConfigDict(strict=True)

    h1_shifts: list[AtomShift]  # 1H chemical shifts, sorted by shift descending
    c13_shifts: list[AtomShift]  # 13C chemical shifts, sorted by shift descending
    functional: str  # DFT functional used
    basis_set: str  # Basis set used for NMR
    solvent: str  # COSMO solvent used


class StepTiming(BaseModel):
    """Timing information for a completed calculation step."""

    # Note: strict=False to allow datetime string coercion from JSON
    model_config = ConfigDict(strict=False)

    step: str  # Step identifier
    label: str  # Human-readable label
    started_at: datetime
    completed_at: datetime
    duration_seconds: float


class JobInput(BaseModel):
    """Input parameters for a calculation job."""

    model_config = ConfigDict(strict=True)

    smiles: str
    name: Optional[str] = None  # User-provided molecule name/label
    preset: Literal["draft", "production"] = "production"
    solvent: str  # Required field - no default, user must specify
    notification_email: Optional[str] = None  # Opt-in email for completion notification


class JobStatus(BaseModel):
    """Complete status of a calculation job."""

    model_config = ConfigDict(
        # Use 'iso8601' serialization for datetime fields
        ser_json_timedelta="iso8601",
        extra="ignore",  # Allow old files with isicle_version to load
    )

    job_id: str
    status: Literal["queued", "running", "complete", "failed"]
    created_at: datetime
    started_at: Optional[datetime] = None
    completed_at: Optional[datetime] = None

    # Input
    input: JobInput

    # Version (for reproducibility)
    nwchem_version: str

    # Step progress tracking
    current_step: Optional[str] = None  # Current step identifier
    current_step_label: Optional[str] = None  # Human-readable label
    step_started_at: Optional[datetime] = None  # When current step started
    steps_completed: list[StepTiming] = []  # Completed steps with timings

    # Error info (if failed)
    error_message: Optional[str] = None
    error_traceback: Optional[str] = None

    # Resource usage (filled on completion)
    cpu_time_seconds: Optional[float] = None
    memory_peak_mb: Optional[float] = None

    # NMR results (populated on completion)
    nmr_results: Optional[NMRResults] = None
    optimized_geometry_file: Optional[str] = None  # Path to XYZ file
