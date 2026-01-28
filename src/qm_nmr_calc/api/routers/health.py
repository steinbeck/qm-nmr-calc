"""Health check endpoints for liveness and readiness probes."""

from datetime import datetime
from pathlib import Path

from fastapi import APIRouter
from fastapi.responses import JSONResponse

router = APIRouter(tags=["health"])


@router.get("/health")
async def liveness():
    """Liveness probe - is the service running?

    Returns minimal response for load balancer probes.
    No dependency checks - just confirms the service is up.
    """
    return {"status": "alive"}


@router.get("/health/ready")
async def readiness():
    """Readiness probe - can the service handle requests?

    Checks:
    - Data directory is writable
    - Huey task queue is importable
    - CREST/xTB binary availability
    """
    checks = {}

    # Check data directory writable
    from ...storage import DATA_DIR

    try:
        DATA_DIR.mkdir(parents=True, exist_ok=True)
        test_file = DATA_DIR / ".health_check"
        test_file.write_text("ok")
        test_file.unlink()
        checks["data_directory"] = "ok"
    except Exception as e:
        checks["data_directory"] = f"error: {e}"
        return JSONResponse(
            status_code=503,
            content={"status": "not ready", "checks": checks},
        )

    # Check Huey queue accessible
    try:
        from ...queue import huey  # noqa: F401

        checks["task_queue"] = "ok"
    except Exception as e:
        checks["task_queue"] = f"error: {e}"
        return JSONResponse(
            status_code=503,
            content={"status": "not ready", "checks": checks},
        )

    # Check CREST availability (non-blocking, informational only)
    from ...conformers.crest_generator import detect_crest_available

    checks["crest_available"] = detect_crest_available()

    return {
        "status": "ready",
        "checks": checks,
        "crest_available": checks["crest_available"],
        "timestamp": datetime.utcnow().isoformat() + "Z",
    }
