"""Calculation preset configurations for NMR calculations."""

import os
from enum import Enum
from typing import TypedDict


def _get_available_memory_gb() -> float:
    """Get available memory in GB (container limit or system memory)."""
    # Try cgroup v2 first (modern Docker/Kubernetes)
    try:
        with open("/sys/fs/cgroup/memory.max", "r") as f:
            val = f.read().strip()
            if val != "max":
                return int(val) / (1024**3)
    except (FileNotFoundError, ValueError, PermissionError):
        pass

    # Try cgroup v1 (older systems)
    try:
        with open("/sys/fs/cgroup/memory/memory.limit_in_bytes", "r") as f:
            val = int(f.read().strip())
            # Very high value means no limit set
            if val < 2**62:
                return val / (1024**3)
    except (FileNotFoundError, ValueError, PermissionError):
        pass

    # Fall back to system memory
    try:
        import shutil
        total, _, _ = shutil.disk_usage("/")  # dummy call to check we can import
        with open("/proc/meminfo", "r") as f:
            for line in f:
                if line.startswith("MemTotal:"):
                    # MemTotal is in kB
                    return int(line.split()[1]) / (1024**2)
    except (FileNotFoundError, ValueError, PermissionError):
        pass

    return 8.0  # Conservative default


def get_default_processes() -> int:
    """Get default number of MPI processes from environment or auto-detect.

    Uses NWCHEM_NPROC env var if set, otherwise auto-detects based on:
    1. Available CPUs: uses 80% (to leave headroom for system)
    2. Available memory: requires 2GB per MPI process (NWChem recommendation)
    3. Hard cap at 40 (diminishing returns beyond that for typical molecules)

    The minimum of CPU-based and memory-based limits is used.
    """
    env_val = os.environ.get("NWCHEM_NPROC")
    if env_val:
        try:
            return int(env_val)
        except ValueError:
            pass

    try:
        # CPU-based limit: 80% of available cores
        cpu_count = os.cpu_count() or 4
        cpu_limit = int(cpu_count * 0.8)

        # Memory-based limit: 2GB per MPI process
        memory_gb = _get_available_memory_gb()
        memory_limit = int(memory_gb / 2)

        # Use the more restrictive limit, cap at 40, minimum 1
        target = min(cpu_limit, memory_limit)
        return max(1, min(target, 40))
    except Exception:
        return 4


class PresetName(str, Enum):
    """Available calculation presets."""

    DRAFT = "draft"
    PRODUCTION = "production"


class CalculationPreset(TypedDict):
    """Configuration for a calculation preset.

    Attributes:
        name: Human-readable preset name
        description: Description of the preset's purpose
        functional: DFT functional (e.g., 'b3lyp')
        basis_set: Basis set for geometry optimization
        nmr_basis_set: Basis set for NMR shielding calculation
        processes: Number of parallel processes
        max_iter: Maximum geometry optimization iterations
    """

    name: str
    description: str
    functional: str
    basis_set: str
    nmr_basis_set: str
    processes: int
    max_iter: int


# Compute default processes once at module load
_DEFAULT_PROCESSES = get_default_processes()

PRESETS: dict[PresetName, CalculationPreset] = {
    PresetName.DRAFT: {
        "name": "draft",
        "description": "Fast calculations for quick checks",
        "functional": "b3lyp",
        "basis_set": "6-31G*",
        "nmr_basis_set": "6-31G*",
        "processes": _DEFAULT_PROCESSES,
        "max_iter": 100,
    },
    PresetName.PRODUCTION: {
        "name": "production",
        "description": "Balanced accuracy and compute time (default)",
        "functional": "b3lyp",
        "basis_set": "6-31G*",
        "nmr_basis_set": "6-311+G(2d,p)",
        "processes": _DEFAULT_PROCESSES,
        "max_iter": 300,
    },
}

DEFAULT_PRESET = PresetName.PRODUCTION
