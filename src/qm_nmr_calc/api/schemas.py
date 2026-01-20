"""Pydantic models for API requests and responses."""

from typing import Literal, Optional

from pydantic import BaseModel, EmailStr, Field


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
    preset: Literal["draft", "production"] = Field(
        default="production",
        description="Calculation preset: draft (fast) or production (accurate)",
    )
    solvent: str = Field(
        ...,
        description="NMR solvent for COSMO solvation model (e.g., chcl3, dmso)",
        examples=["chcl3", "dmso", "acetone"],
    )
    notification_email: Optional[EmailStr] = Field(
        None,
        description="Email address for completion notification (opt-in)",
        examples=["user@example.com"],
    )


class AtomShiftResponse(BaseModel):
    """API response model for a single atom chemical shift."""

    index: int = Field(..., description="Atom index (1-based)")
    atom: str = Field(..., description="Element symbol (H or C)")
    shift: float = Field(..., description="Chemical shift in ppm (TMS reference)")


class NMRResultsResponse(BaseModel):
    """API response model for NMR calculation results."""

    h1_shifts: list[AtomShiftResponse] = Field(
        ..., description="1H chemical shifts sorted by ppm descending"
    )
    c13_shifts: list[AtomShiftResponse] = Field(
        ..., description="13C chemical shifts sorted by ppm descending"
    )
    functional: str = Field(..., description="DFT functional used")
    basis_set: str = Field(..., description="Basis set used for NMR calculation")
    solvent: str = Field(..., description="Solvent used for COSMO solvation")


class StepTimingResponse(BaseModel):
    """API response model for a completed calculation step."""

    step: str = Field(..., description="Step identifier")
    label: str = Field(..., description="Human-readable step label")
    started_at: str = Field(..., description="ISO 8601 timestamp when step started")
    completed_at: str = Field(..., description="ISO 8601 timestamp when step completed")
    duration_seconds: float = Field(..., description="Step duration in seconds")


class JobStatusResponse(BaseModel):
    """Response for job status queries."""

    job_id: str = Field(..., description="12-character hex job ID")
    status: str = Field(
        ..., description="Job status: queued, running, complete, failed"
    )
    created_at: str = Field(..., description="ISO 8601 timestamp")
    started_at: Optional[str] = Field(None, description="When job started running")
    completed_at: Optional[str] = Field(None, description="When job completed")
    input_smiles: str = Field(..., description="Original SMILES input")
    input_name: Optional[str] = Field(None, description="User-provided molecule name")
    preset: str = Field(..., description="Calculation preset used")
    solvent: str = Field(..., description="Solvent used for calculation")
    # Step progress tracking
    current_step: Optional[str] = Field(None, description="Current step identifier")
    current_step_label: Optional[str] = Field(None, description="Current step label")
    step_started_at: Optional[str] = Field(None, description="When current step started")
    steps_completed: list[StepTimingResponse] = Field(
        default_factory=list, description="Completed steps with timings"
    )
    error_message: Optional[str] = Field(None, description="Error message if failed")
    nmr_results: Optional[NMRResultsResponse] = Field(
        None, description="NMR results (available when status is complete)"
    )


class ProblemDetail(BaseModel):
    """RFC 7807 Problem Details response for API errors."""

    type: str = Field(
        default="about:blank",
        description="URI identifying the error type",
    )
    title: str = Field(
        ...,
        description="Short human-readable summary of the problem",
    )
    status: int = Field(
        ...,
        description="HTTP status code",
    )
    detail: Optional[str] = Field(
        None,
        description="Detailed explanation of the problem",
    )
    instance: Optional[str] = Field(
        None,
        description="URI of the specific occurrence of the problem",
    )
