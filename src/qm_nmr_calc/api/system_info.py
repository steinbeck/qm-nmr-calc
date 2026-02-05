"""System information for template context."""

import os
import socket

import httpx

from ..presets import _get_available_memory_gb, get_default_processes

# Memory per MPI process (NWChem recommendation from presets.py)
_MEMORY_PER_PROCESS_GB = 2


def _get_machine_type() -> str:
    """Get GCP machine type or hostname for local development.

    Queries GCP metadata server with 0.5s timeout.
    Falls back to hostname if not on GCP.
    """
    try:
        response = httpx.get(
            "http://metadata.google.internal/computeMetadata/v1/instance/machine-type",
            headers={"Metadata-Flavor": "Google"},
            timeout=0.5,
        )
        if response.status_code == 200:
            # Response is like "projects/123456/zones/us-central1-a/machineTypes/e2-standard-4"
            # Extract just the machine type
            return response.text.split("/")[-1]
    except (httpx.TimeoutException, httpx.ConnectError, httpx.RequestError):
        pass

    return socket.gethostname()


def _get_total_cores() -> int:
    """Get total CPU cores on the machine."""
    return os.cpu_count() or 4


# Cache system info at module load (same pattern as presets.py)
_nwchem_processes = get_default_processes()
_SYSTEM_INFO = {
    # Machine specs
    "machine_type": _get_machine_type(),
    "total_cores": _get_total_cores(),
    "total_memory_gb": round(_get_available_memory_gb()),
    # NWChem allocation
    "nwchem_processes": _nwchem_processes,
    "nwchem_memory_gb": _nwchem_processes * _MEMORY_PER_PROCESS_GB,
}


def get_cpu_cores() -> int:
    """Get total CPU cores on the machine."""
    return _SYSTEM_INFO["total_cores"]


def get_memory_gb() -> int:
    """Get total available memory in GB."""
    return _SYSTEM_INFO["total_memory_gb"]


def get_machine_type() -> str:
    """Get GCP machine type or hostname."""
    return _SYSTEM_INFO["machine_type"]


def get_nwchem_processes() -> int:
    """Get number of MPI processes NWChem will use."""
    return _SYSTEM_INFO["nwchem_processes"]


def get_nwchem_memory_gb() -> int:
    """Get memory NWChem will use (2GB per process)."""
    return _SYSTEM_INFO["nwchem_memory_gb"]


def get_system_info() -> dict:
    """Get combined system info dict for template context."""
    return _SYSTEM_INFO.copy()
