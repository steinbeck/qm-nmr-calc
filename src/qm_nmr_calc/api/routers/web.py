"""Web UI router for browser-based interface."""

from pathlib import Path
from typing import Annotated, Optional

from fastapi import APIRouter, Form, Request, UploadFile
from fastapi.responses import HTMLResponse, RedirectResponse
from fastapi.templating import Jinja2Templates
from rdkit import Chem

from ...nwchem import get_nwchem_version
from ...presets import PRESETS
from ...solvents import (
    SUPPORTED_SOLVENTS,
    get_solvent_display_name,
    get_supported_solvents,
    validate_solvent,
)
from ...storage import create_job_directory, get_job_dir, load_job_status
from ...tasks import run_nmr_task, run_ensemble_nmr_task, _generate_initial_xyz
from ...validation import validate_mol_file, validate_smiles

# Template engine setup
templates = Jinja2Templates(
    directory=Path(__file__).resolve().parent.parent / "templates"
)

router = APIRouter(tags=["web"])


def _get_form_context() -> dict:
    """Build common form context with solvents, presets, and CREST availability."""
    from ...conformers import detect_crest_available

    # Build solvent options with descriptions
    solvents = [
        {"value": key, "label": desc}
        for key, desc in sorted(SUPPORTED_SOLVENTS.items(), key=lambda x: x[1])
    ]

    # Build preset options with descriptions
    presets = [
        {"value": preset.value, "label": config["description"]}
        for preset, config in PRESETS.items()
    ]

    # Check CREST availability for form display
    crest_available = detect_crest_available()

    return {"solvents": solvents, "presets": presets, "crest_available": crest_available}


@router.get("/", response_class=HTMLResponse)
async def home(request: Request) -> HTMLResponse:
    """Render the submission form page."""
    context = _get_form_context()
    context["form_data"] = {}  # Empty form data for fresh submission
    return templates.TemplateResponse(
        request=request,
        name="submit.html",
        context=context,
    )


@router.post("/submit", response_class=HTMLResponse)
async def submit_job(
    request: Request,
    smiles: Annotated[Optional[str], Form()] = None,
    file: Optional[UploadFile] = None,
    solvent: Annotated[str, Form()] = "",
    preset: Annotated[str, Form()] = "production",
    name: Annotated[Optional[str], Form()] = None,
    notification_email: Annotated[Optional[str], Form()] = None,
    conformer_mode: Annotated[str, Form()] = "ensemble",  # Default to ensemble per CONTEXT.md
    conformer_method: Annotated[Optional[str], Form()] = None,  # None = use default (rdkit_kdg)
    max_conformers: Annotated[Optional[str], Form()] = None,  # String from form, convert to int
):
    """Process job submission from web form."""
    # Preserve form data for re-render on error
    form_data = {
        "smiles": smiles or "",
        "solvent": solvent,
        "preset": preset,
        "name": name or "",
        "notification_email": notification_email or "",
        "conformer_mode": conformer_mode,
        "conformer_method": conformer_method or "",
        "max_conformers": max_conformers or "",
    }

    def render_error(error: str) -> HTMLResponse:
        """Helper to render form with error message."""
        context = _get_form_context()
        context["error"] = error
        context["form_data"] = form_data
        return templates.TemplateResponse(
            request=request,
            name="submit.html",
            context=context,
        )

    # Validate that at least one input is provided
    smiles_provided = smiles and smiles.strip()
    file_provided = file and file.filename

    if not smiles_provided and not file_provided:
        return render_error("Please provide either a SMILES string or upload a MOL/SDF file.")

    # Process input and get molecule
    final_smiles = None
    final_name = name

    if smiles_provided:
        # Validate SMILES
        mol, error = validate_smiles(smiles.strip())
        if mol is None:
            return render_error(f"Invalid SMILES: {error}")
        final_smiles = smiles.strip()

    elif file_provided:
        # Read and validate file
        filename = file.filename or "uploaded.mol"
        if not (filename.endswith(".mol") or filename.endswith(".sdf")):
            return render_error("File must be a .mol or .sdf file.")

        content = await file.read()
        mol, error = validate_mol_file(content, filename)
        if mol is None:
            return render_error(f"Invalid molecule file: {error}")

        # Convert to SMILES for storage
        final_smiles = Chem.MolToSmiles(mol)

        # Use filename as name if not provided
        if not final_name:
            final_name = filename

    # Validate solvent
    if not solvent:
        return render_error("Please select a solvent.")

    normalized_solvent = validate_solvent(solvent)
    if normalized_solvent is None:
        supported = ", ".join(get_supported_solvents())
        return render_error(f"Unknown solvent '{solvent}'. Supported: {supported}")

    # Get software versions
    nwchem_version = get_nwchem_version()

    # Convert max_conformers from string to int if provided
    max_conformers_int = None
    if max_conformers and max_conformers.strip():
        try:
            max_conformers_int = int(max_conformers.strip())
            if max_conformers_int < 1 or max_conformers_int > 500:
                return render_error("Max conformers must be between 1 and 500.")
        except ValueError:
            return render_error("Max conformers must be a number.")

    # If mode is single, ignore conformer_method
    effective_method = None if conformer_mode == "single" else (conformer_method or None)

    # Create job directory and initial status
    job_status = create_job_directory(
        smiles=final_smiles,
        solvent=normalized_solvent,
        nwchem_version=nwchem_version,
        name=final_name,
        preset=preset,
        notification_email=notification_email,
        conformer_mode=conformer_mode,
        conformer_method=effective_method,  # Convert empty string to None
        max_conformers=max_conformers_int,
    )

    # Generate initial 3D geometry for immediate visualization
    job_dir = get_job_dir(job_status.job_id)
    initial_xyz_path = job_dir / "output" / "initial.xyz"
    _generate_initial_xyz(final_smiles, initial_xyz_path)

    # Queue the appropriate NMR calculation task based on conformer mode
    if job_status.input.conformer_mode == "ensemble":
        run_ensemble_nmr_task(job_status.job_id)
    else:
        run_nmr_task(job_status.job_id)

    # Redirect to status page with 303 See Other
    return RedirectResponse(
        url=f"/status/{job_status.job_id}",
        status_code=303,
    )


@router.get("/status/{job_id}", response_class=HTMLResponse)
async def job_status_page(request: Request, job_id: str) -> HTMLResponse:
    """Render the job status page."""
    job_status = load_job_status(job_id)

    if job_status is None:
        # Redirect to home with error
        context = _get_form_context()
        context["error"] = f"Job '{job_id}' not found."
        context["form_data"] = {}
        return templates.TemplateResponse(
            request=request,
            name="submit.html",
            context=context,
            status_code=404,
        )

    # Convert job status to template-friendly dict
    job_data = {
        "job_id": job_status.job_id,
        "status": job_status.status,
        "created_at": job_status.created_at.isoformat() + "Z",
        "input_smiles": job_status.input.smiles,
        "input_name": job_status.input.name,
        "solvent": get_solvent_display_name(job_status.input.solvent),
        "preset": job_status.input.preset,
        "error_message": job_status.error_message,
        "conformer_mode": job_status.input.conformer_mode,
    }

    return templates.TemplateResponse(
        request=request,
        name="status.html",
        context={"job": job_data},
    )


@router.get("/results/{job_id}", response_class=HTMLResponse)
async def results_page(request: Request, job_id: str) -> HTMLResponse:
    """Render the results page for a completed job."""
    job_status = load_job_status(job_id)

    # Job not found
    if job_status is None:
        context = _get_form_context()
        context["error"] = f"Job '{job_id}' not found."
        context["form_data"] = {}
        return templates.TemplateResponse(
            request=request,
            name="submit.html",
            context=context,
            status_code=404,
        )

    # Job not complete - redirect to status page
    if job_status.status != "complete":
        return RedirectResponse(url=f"/status/{job_id}", status_code=303)

    # Job complete but no results
    if job_status.nmr_results is None:
        context = _get_form_context()
        context["error"] = "Job completed but no results available."
        context["form_data"] = {}
        return templates.TemplateResponse(
            request=request,
            name="submit.html",
            context=context,
            status_code=500,
        )

    # Build context for results template
    job_context = {
        "job_id": job_status.job_id,
        "input_smiles": job_status.input.smiles,
        "input_name": job_status.input.name,
        "preset": job_status.input.preset,
        "solvent": get_solvent_display_name(job_status.input.solvent),
        "status": job_status.status,
        "conformer_mode": job_status.input.conformer_mode,
    }

    results_context = {
        "functional": job_status.nmr_results.functional,
        "basis_set": job_status.nmr_results.basis_set,
        "solvent": get_solvent_display_name(job_status.nmr_results.solvent),
    }

    return templates.TemplateResponse(
        request=request,
        name="results.html",
        context={
            "job": job_context,
            "results": results_context,
        },
    )
