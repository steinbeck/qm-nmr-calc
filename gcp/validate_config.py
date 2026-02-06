"""GCP deployment config validation using Pydantic v2 and TOML.

This module validates GCP deployment configuration from a TOML file,
ensuring project IDs, CPU/RAM ratios, and resource limits are valid
before any GCP operations run.

Usage as module:
    from gcp.validate_config import GCPConfig, load_config, format_exports

Usage as script:
    python3 gcp/validate_config.py [--config path/to/config.toml]
"""

from __future__ import annotations

import argparse
import re
import sys
import tomllib
from pathlib import Path

from pydantic import BaseModel, Field, ValidationError, field_validator


class GCPConfig(BaseModel):
    """Validated GCP deployment configuration."""

    project_id: str
    resource_prefix: str = "qm-nmr-calc"
    cpu_cores: int = Field(ge=4, le=224)
    ram_gb: int = Field(ge=1, le=896)
    disk_size_gb: int = Field(default=100, ge=10, le=65536)

    @field_validator("project_id")
    @classmethod
    def validate_project_id(cls, v: str) -> str:
        """Validate GCP project ID format (6-30 chars, lowercase/digits/hyphens)."""
        if len(v) < 6:
            raise ValueError(
                f"Project ID must be at least 6 characters, got {len(v)}"
            )
        if len(v) > 30:
            raise ValueError(
                f"Project ID must be at most 30 characters, got {len(v)}"
            )
        if v != v.lower():
            raise ValueError(
                "Project ID must be lowercase (no uppercase letters allowed)"
            )
        if not re.match(r"^[a-z][a-z0-9-]*[a-z0-9]$", v):
            if v.startswith("-"):
                raise ValueError("Project ID must not start with a hyphen")
            if v.endswith("-"):
                raise ValueError("Project ID must not end with a hyphen")
            raise ValueError(
                "Project ID must contain only lowercase letters, digits, and hyphens"
            )
        return v

    @field_validator("ram_gb")
    @classmethod
    def validate_ram_ratio(cls, v: int, info) -> int:
        """Validate RAM/CPU ratio is between 0.5 and 8.0 GB per core."""
        cpu_cores = info.data.get("cpu_cores")
        if cpu_cores is not None:
            ratio = v / cpu_cores
            if ratio < 0.5:
                raise ValueError(
                    f"RAM/CPU ratio too low: {v}GB / {cpu_cores} cores = "
                    f"{ratio:.1f} GB/core (minimum 0.5 GB/core)"
                )
            if ratio > 8.0:
                raise ValueError(
                    f"RAM/CPU ratio too high: {v}GB / {cpu_cores} cores = "
                    f"{ratio:.1f} GB/core (maximum 8.0 GB/core)"
                )
        return v


def load_config(path: Path) -> GCPConfig:
    """Load and validate a GCP config from a TOML file.

    Args:
        path: Path to the TOML config file.

    Returns:
        Validated GCPConfig instance.

    Raises:
        FileNotFoundError: If the config file doesn't exist.
        ValueError: If the TOML is invalid or missing [gcp] section.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")

    try:
        with open(path, "rb") as f:
            data = tomllib.load(f)
    except tomllib.TOMLDecodeError as e:
        raise ValueError(f"Invalid TOML syntax: {e}") from e

    if "gcp" not in data:
        raise ValueError(
            "Config file missing [gcp] section. "
            "Add a [gcp] section header before your settings."
        )

    return GCPConfig(**data["gcp"])


def format_exports(config: GCPConfig) -> str:
    """Format config as bash-compatible export statements.

    Uses single quotes around string values to prevent shell injection.

    Args:
        config: Validated GCPConfig instance.

    Returns:
        Multi-line string of bash export statements.
    """
    lines = [
        f"export GCP_PROJECT_ID='{config.project_id}'",
        f"export CPU_CORES={config.cpu_cores}",
        f"export RAM_GB={config.ram_gb}",
        f"export DISK_SIZE_GB={config.disk_size_gb}",
        f"export RESOURCE_PREFIX='{config.resource_prefix}'",
    ]
    return "\n".join(lines)


def main() -> None:
    """CLI entrypoint: validate config and print exports."""
    parser = argparse.ArgumentParser(
        description="Validate GCP deployment configuration"
    )
    script_dir = Path(__file__).parent
    default_config = script_dir / "config.toml"
    parser.add_argument(
        "--config",
        type=Path,
        default=default_config,
        help=f"Path to TOML config file (default: {default_config})",
    )
    args = parser.parse_args()

    try:
        config = load_config(args.config)
    except FileNotFoundError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        print(
            "Create from template: cp gcp/config.toml.example gcp/config.toml",
            file=sys.stderr,
        )
        sys.exit(1)
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)
    except ValidationError as e:
        for error in e.errors():
            field = ".".join(str(loc) for loc in error["loc"])
            msg = error["msg"]
            print(f"ERROR: {field}: {msg}", file=sys.stderr)
        sys.exit(1)

    print(format_exports(config))


if __name__ == "__main__":
    main()
