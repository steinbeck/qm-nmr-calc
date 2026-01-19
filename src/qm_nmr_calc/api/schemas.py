"""Pydantic models for API requests and responses."""

from typing import Optional

from pydantic import BaseModel, Field


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
    status: str = Field(
        ..., description="Job status: queued, running, complete, failed"
    )
    created_at: str = Field(..., description="ISO 8601 timestamp")
    started_at: Optional[str] = Field(None, description="When job started running")
    completed_at: Optional[str] = Field(None, description="When job completed")
    input_smiles: str = Field(..., description="Original SMILES input")
    input_name: Optional[str] = Field(None, description="User-provided molecule name")
    error_message: Optional[str] = Field(None, description="Error message if failed")


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
