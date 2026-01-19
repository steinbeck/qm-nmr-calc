"""FastAPI application for QM NMR Calculator API."""

from fastapi import FastAPI

from .routers import health, jobs

app = FastAPI(
    title="QM NMR Calculator API",
    description="Asynchronous NMR quantum mechanical calculations via ISiCLE/NWChem",
    version="0.1.0",
    docs_url="/docs",
    redoc_url="/redoc",
    openapi_url="/api/v1/openapi.json",
)

# Health endpoints at root (no version prefix)
app.include_router(health.router)

# API v1 endpoints
app.include_router(jobs.router, prefix="/api/v1")
