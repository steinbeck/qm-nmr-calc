"""Web UI router for browser-based interface."""

from pathlib import Path

from fastapi import APIRouter, Request
from fastapi.responses import HTMLResponse, RedirectResponse
from fastapi.templating import Jinja2Templates

from ...presets import PRESETS
from ...solvents import SUPPORTED_SOLVENTS
from ...storage import load_job_status

# Template engine setup
templates = Jinja2Templates(
    directory=Path(__file__).resolve().parent.parent / "templates"
)

router = APIRouter(tags=["web"])


@router.get("/", response_class=HTMLResponse)
async def home(request: Request) -> HTMLResponse:
    """Render the submission form page."""
    # Build solvent options with descriptions
    solvents = [
        {"value": key, "label": desc}
        for key, desc in sorted(SUPPORTED_SOLVENTS.items(), key=lambda x: x[1])
    ]

    # Build preset options with descriptions
    presets = [
        {"value": preset.name, "label": config["description"]}
        for preset, config in PRESETS.items()
    ]

    return templates.TemplateResponse(
        request=request,
        name="submit.html",
        context={
            "solvents": solvents,
            "presets": presets,
        },
    )


@router.get("/results/{job_id}", response_class=HTMLResponse)
async def results(request: Request, job_id: str) -> HTMLResponse:
    """Render the results page for a completed job."""
    job_status = load_job_status(job_id)

    # Job not found
    if job_status is None:
        return templates.TemplateResponse(
            request=request,
            name="submit.html",
            context={
                "solvents": [
                    {"value": key, "label": desc}
                    for key, desc in sorted(SUPPORTED_SOLVENTS.items(), key=lambda x: x[1])
                ],
                "presets": [
                    {"value": preset.name, "label": config["description"]}
                    for preset, config in PRESETS.items()
                ],
                "error": f"Job '{job_id}' not found.",
            },
            status_code=404,
        )

    # Job not complete - redirect to status page
    if job_status.status != "complete":
        return RedirectResponse(url=f"/status/{job_id}", status_code=303)

    # Job complete but no results
    if job_status.nmr_results is None:
        return templates.TemplateResponse(
            request=request,
            name="submit.html",
            context={
                "solvents": [
                    {"value": key, "label": desc}
                    for key, desc in sorted(SUPPORTED_SOLVENTS.items(), key=lambda x: x[1])
                ],
                "presets": [
                    {"value": preset.name, "label": config["description"]}
                    for preset, config in PRESETS.items()
                ],
                "error": "Job completed but no results available.",
            },
            status_code=500,
        )

    # Build context for results template
    job_context = {
        "job_id": job_status.job_id,
        "input_smiles": job_status.input.smiles,
        "input_name": job_status.input.name,
        "preset": job_status.input.preset,
        "solvent": job_status.input.solvent,
        "status": job_status.status,
    }

    results_context = {
        "functional": job_status.nmr_results.functional,
        "basis_set": job_status.nmr_results.basis_set,
        "solvent": job_status.nmr_results.solvent,
    }

    return templates.TemplateResponse(
        request=request,
        name="results.html",
        context={
            "job": job_context,
            "results": results_context,
        },
    )
