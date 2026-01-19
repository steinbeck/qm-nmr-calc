"""Job submission and status endpoints."""

from typing import Annotated, Optional

from fastapi import APIRouter, Form, HTTPException, UploadFile, status
from fastapi.responses import JSONResponse
from rdkit import Chem

from ..schemas import JobStatusResponse, JobSubmitRequest, ProblemDetail
from ...isicle_wrapper import get_versions
from ...models import JobStatus
from ...solvents import validate_solvent, get_supported_solvents
from ...storage import create_job_directory, load_job_status
from ...tasks import run_nmr_task
from ...validation import validate_mol_file, validate_smiles

router = APIRouter(prefix="/jobs", tags=["jobs"])


def job_status_to_response(job_status: JobStatus) -> dict:
    """Convert JobStatus model to JobStatusResponse dict.

    Flattens input.smiles -> input_smiles, input.name -> input_name.
    Converts datetime fields to ISO strings.
    Converts NMR results if present.
    """
    # Convert NMR results if present
    nmr_results = None
    if job_status.nmr_results is not None:
        nmr_results = {
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

    return {
        "job_id": job_status.job_id,
        "status": job_status.status,
        "created_at": job_status.created_at.isoformat() + "Z",
        "started_at": (
            job_status.started_at.isoformat() + "Z"
            if job_status.started_at
            else None
        ),
        "completed_at": (
            job_status.completed_at.isoformat() + "Z"
            if job_status.completed_at
            else None
        ),
        "input_smiles": job_status.input.smiles,
        "input_name": job_status.input.name,
        "preset": job_status.input.preset,
        "solvent": job_status.input.solvent,
        "error_message": job_status.error_message,
        "nmr_results": nmr_results,
    }


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
    """Submit molecule via SMILES string for NMR calculation.

    Returns immediately with job ID. Use GET /api/v1/jobs/{job_id} to poll status.
    """
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

    # Validate solvent
    normalized_solvent = validate_solvent(request.solvent)
    if normalized_solvent is None:
        supported = ", ".join(get_supported_solvents())
        raise HTTPException(
            status_code=422,
            detail={
                "type": "https://qm-nmr-calc.example/problems/invalid-solvent",
                "title": "Invalid Solvent",
                "status": 422,
                "detail": f"Unknown solvent '{request.solvent}'. Supported solvents: {supported}",
            },
        )

    # Get software versions for reproducibility
    versions = get_versions()

    # Create job directory and initial status
    job_status = create_job_directory(
        smiles=request.smiles,
        solvent=normalized_solvent,
        isicle_version=versions.isicle,
        nwchem_version=versions.nwchem,
        name=request.name,
        preset=request.preset,
    )

    # Queue the NMR calculation task
    run_nmr_task(job_status.job_id)

    return JSONResponse(
        status_code=status.HTTP_202_ACCEPTED,
        content=job_status_to_response(job_status),
        headers={
            "Location": f"/api/v1/jobs/{job_status.job_id}",
            "Retry-After": "30",
        },
    )


@router.post(
    "/upload",
    response_model=JobStatusResponse,
    status_code=status.HTTP_202_ACCEPTED,
    responses={
        202: {"description": "Job accepted for processing"},
        422: {"model": ProblemDetail, "description": "Invalid molecule file"},
    },
)
async def submit_file(
    file: UploadFile,
    solvent: Annotated[
        str, Form(description="NMR solvent for COSMO solvation model")
    ],
    name: Annotated[
        Optional[str], Form(description="Optional molecule name")
    ] = None,
    preset: Annotated[
        str, Form(description="Calculation preset: draft or production")
    ] = "production",
):
    """Submit molecule via MOL/SDF file upload for NMR calculation.

    Accepts .mol or .sdf files with a single molecule.
    Returns immediately with job ID. Use GET /api/v1/jobs/{job_id} to poll status.
    """
    # Validate file type by extension
    filename = file.filename or ""
    if not (filename.endswith(".mol") or filename.endswith(".sdf")):
        # Also check by content-type
        allowed_types = {
            "chemical/x-mdl-molfile",
            "chemical/x-mdl-sdfile",
            "application/octet-stream",
        }
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

    # Validate solvent
    normalized_solvent = validate_solvent(solvent)
    if normalized_solvent is None:
        supported = ", ".join(get_supported_solvents())
        raise HTTPException(
            status_code=422,
            detail={
                "type": "https://qm-nmr-calc.example/problems/invalid-solvent",
                "title": "Invalid Solvent",
                "status": 422,
                "detail": f"Unknown solvent '{solvent}'. Supported solvents: {supported}",
            },
        )

    # Convert to SMILES for storage
    smiles = Chem.MolToSmiles(mol)

    # Get software versions
    versions = get_versions()

    # Create job (use filename as name if not provided)
    job_status = create_job_directory(
        smiles=smiles,
        solvent=normalized_solvent,
        isicle_version=versions.isicle,
        nwchem_version=versions.nwchem,
        name=name or filename,
        preset=preset,
    )

    # Queue NMR calculation
    run_nmr_task(job_status.job_id)

    return JSONResponse(
        status_code=status.HTTP_202_ACCEPTED,
        content=job_status_to_response(job_status),
        headers={
            "Location": f"/api/v1/jobs/{job_status.job_id}",
            "Retry-After": "30",
        },
    )


@router.get(
    "/solvents",
    response_model=list[str],
    responses={
        200: {"description": "List of supported solvents"},
    },
)
async def list_solvents():
    """List supported NMR solvents for COSMO solvation model.

    Returns list of NWChem COSMO solvent names that can be used
    in job submission requests.
    """
    return get_supported_solvents()


@router.get(
    "/{job_id}",
    response_model=JobStatusResponse,
    responses={
        200: {"description": "Job status"},
        404: {"model": ProblemDetail, "description": "Job not found"},
    },
)
async def get_job_status(job_id: str):
    """Get status of a calculation job.

    Returns full job status including timestamps and any error message.
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
    return job_status_to_response(job_status)
