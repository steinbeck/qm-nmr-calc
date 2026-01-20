"""FastAPI application for QM NMR Calculator API."""

from pathlib import Path

from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles

from .routers import health, jobs, web

# Base directory for static files and templates
BASE_DIR = Path(__file__).resolve().parent

app = FastAPI(
    title="QM NMR Calculator API",
    description="Asynchronous NMR quantum mechanical calculations via ISiCLE/NWChem",
    version="0.1.0",
    docs_url="/docs",
    redoc_url="/redoc",
    openapi_url="/api/v1/openapi.json",
)

# Mount static files before router includes to avoid route conflicts
app.mount("/static", StaticFiles(directory=BASE_DIR / "static"), name="static")

# Health endpoints at root (no version prefix)
app.include_router(health.router)

# Web UI routes at root (no version prefix)
app.include_router(web.router)

# API v1 endpoints
app.include_router(jobs.router, prefix="/api/v1")
