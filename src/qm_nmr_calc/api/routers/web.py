"""Web UI router for browser-based interface."""

from pathlib import Path

from fastapi import APIRouter, Request
from fastapi.responses import HTMLResponse
from fastapi.templating import Jinja2Templates

from ...presets import PRESETS
from ...solvents import SUPPORTED_SOLVENTS

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
