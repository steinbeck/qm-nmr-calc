"""GCP machine type selection, validation, and resource calculation.

This module maps user CPU/RAM requirements to GCP machine types, validates
availability across zones with fallback, calculates Docker memory limits,
and generates startup scripts for VM deployment.

Core workflow:
    1. find_available_zone() calls get_ranked_regions() to get sorted regions
    2. Iterates through regions trying select_machine_type() until one succeeds
    3. calculate_docker_resources() computes Docker limits (VM RAM - 8GB OS overhead)
    4. generate_startup_script() creates the VM startup script

Usage as module:
    from gcp.select_machine import find_available_zone, calculate_docker_resources

Usage as script:
    python3 gcp/select_machine.py --cpu-cores 8 --ram-gb 32
    python3 gcp/select_machine.py --generate-startup-script --cpu-cores 8 --ram-gb 32
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys

from gcp.query_pricing import get_ranked_regions

# General-purpose and compute-optimized families safe for NMR workloads.
# Excludes exotic families: ultramem, hypermem, GPU-attached (a2/a3/g2), m1/m2/m3.
ALLOWED_MACHINE_FAMILIES = {
    "c2", "c3", "c3d", "c4",
    "e2",
    "n1", "n2", "n2d", "n4",
    "t2d", "t2a",
}


# ---------------------------------------------------------------------------
# Internal gcloud wrapper (single mock target for all tests)
# ---------------------------------------------------------------------------


def _run_gcloud(args: list[str]) -> str:
    """Run gcloud command and return stdout.

    Args:
        args: Command-line arguments to pass to gcloud.

    Returns:
        Command stdout stripped of whitespace.

    Raises:
        RuntimeError: If gcloud exits with non-zero status.
    """
    cmd = ["gcloud"] + args
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)

    if result.returncode != 0:
        raise RuntimeError(
            f"gcloud command failed with exit code {result.returncode}: "
            f"{result.stderr.strip()}"
        )

    return result.stdout.strip()


# ---------------------------------------------------------------------------
# Machine type selection
# ---------------------------------------------------------------------------


def select_machine_type(cpu_cores: int, ram_gb: int, zone: str) -> str:
    """Select the smallest GCP machine type matching CPU/RAM requirements.

    Args:
        cpu_cores: Minimum CPU cores required.
        ram_gb: Minimum RAM in GB required.
        zone: GCP zone to search (e.g., "us-central1-a").

    Returns:
        Machine type name (e.g., "n2-standard-8").

    Raises:
        ValueError: If no machine types match requirements.
        RuntimeError: If gcloud command fails.
    """
    ram_mb = ram_gb * 1024

    args = [
        "compute",
        "machine-types",
        "list",
        f"--zones={zone}",
        f"--filter=guestCpus>={cpu_cores} AND memoryMb>={ram_mb}",
        "--format=json",
        "--quiet",
    ]

    output = _run_gcloud(args)
    machines = json.loads(output)

    # Filter to general-purpose/compute-optimized families only
    machines = [
        m for m in machines
        if m["name"].split("-")[0] in ALLOWED_MACHINE_FAMILIES
    ]

    if not machines:
        raise ValueError(
            f"No machine types found matching {cpu_cores} CPU cores "
            f"and {ram_gb}GB RAM in zone {zone}"
        )

    # Sort by estimated cost: CPUs + RAM_GB approximates spot pricing
    # gcloud --sort-by sorts lexicographically, not numerically
    # This avoids ultramem (high RAM, low CPU) and highcpu (high CPU, low RAM)
    machines.sort(key=lambda m: m["guestCpus"] + m["memoryMb"] / 1024)

    return machines[0]["name"]


def validate_machine_in_zone(machine_type: str, zone: str) -> bool:
    """Check if a machine type exists in a specific zone.

    Args:
        machine_type: Machine type name to validate.
        zone: GCP zone to check.

    Returns:
        True if machine type exists in zone, False otherwise.
    """
    args = [
        "compute",
        "machine-types",
        "list",
        f"--zones={zone}",
        f"--filter=name={machine_type}",
        "--format=value(name)",
        "--quiet",
    ]

    try:
        output = _run_gcloud(args)
        return output.strip() == machine_type
    except RuntimeError:
        return False


def find_available_zone(cpu_cores: int, ram_gb: int) -> dict:
    """Find the cheapest zone with a machine type matching requirements.

    Calls get_ranked_regions(cpu_cores, ram_gb) internally to obtain the
    sorted region list from query_pricing.py (cache/API/fallback).

    Iterates through regions in price order, trying to select a machine type
    in each zone. Returns the first successful match.

    Args:
        cpu_cores: Minimum CPU cores required.
        ram_gb: Minimum RAM in GB required.

    Returns:
        Dict with keys: zone, region, machine_type.

    Raises:
        RuntimeError: If all regions exhausted without finding a match.
    """
    ranked_regions = get_ranked_regions(cpu_cores=cpu_cores, ram_gb=ram_gb)

    for region_entry in ranked_regions:
        zone = region_entry["zone"]
        region = region_entry["region"]

        try:
            machine_type = select_machine_type(cpu_cores, ram_gb, zone)
            return {
                "zone": zone,
                "region": region,
                "machine_type": machine_type,
            }
        except (ValueError, RuntimeError) as e:
            print(
                f"WARNING: Failed to find machine in {zone}: {e}",
                file=sys.stderr,
            )
            continue

    raise RuntimeError(
        f"All regions exhausted: no zone has a machine type matching "
        f"{cpu_cores} CPU cores and {ram_gb}GB RAM"
    )


# ---------------------------------------------------------------------------
# Resource calculation
# ---------------------------------------------------------------------------


def calculate_docker_resources(machine_type: str, zone: str) -> dict:
    """Calculate Docker resource limits for a given machine type.

    Subtracts 8GB OS overhead from total VM RAM for worker container limit.

    Args:
        machine_type: GCP machine type name.
        zone: GCP zone where machine type exists.

    Returns:
        Dict with keys: worker_memory_limit, nwchem_nproc, total_ram_gb, total_cpus.

    Raises:
        ValueError: If RAM after OS overhead is less than 4GB.
        RuntimeError: If gcloud command fails.
    """
    args = [
        "compute",
        "machine-types",
        "describe",
        machine_type,
        f"--zone={zone}",
        "--format=json",
        "--quiet",
    ]

    output = _run_gcloud(args)
    machine_info = json.loads(output)

    total_cpus = machine_info["guestCpus"]
    memory_mb = machine_info["memoryMb"]
    total_ram_gb = memory_mb // 1024

    # Subtract 8GB OS overhead
    worker_memory_gb = total_ram_gb - 8

    if worker_memory_gb < 4:
        raise ValueError(
            f"Insufficient RAM after OS overhead: {total_ram_gb}GB total - 8GB overhead "
            f"= {worker_memory_gb}GB available (minimum 4GB required)"
        )

    return {
        "worker_memory_limit": f"{worker_memory_gb}g",
        "nwchem_nproc": total_cpus,
        "total_ram_gb": total_ram_gb,
        "total_cpus": total_cpus,
    }


# ---------------------------------------------------------------------------
# Startup script generation
# ---------------------------------------------------------------------------


def generate_startup_script(
    worker_memory_limit: str, resource_prefix: str, disk_size_gb: int
) -> str:
    r"""Generate GCP VM startup script for qm-nmr-calc deployment.

    The script:
    - Installs Docker if not present
    - Mounts persistent disk
    - Pulls images from GHCR
    - Downloads docker-compose.yml from GitHub
    - Creates HTTP-only docker-compose.gcp.yml override (port 80, no Caddy/HTTPS)
    - Writes .env with WORKER_MEMORY_LIMIT and NWCHEM_NPROC=$(nproc)
    - Starts services with docker compose
    - Disables Caddy for HTTP-only deployment (API serves port 80 directly)

    Args:
        worker_memory_limit: Docker memory limit (e.g., "24g").
        resource_prefix: Resource name prefix (e.g., "qm-nmr-calc").
        disk_size_gb: Persistent disk size in GB.

    Returns:
        Complete startup script as a string.
    """
    return f"""#!/bin/bash
# GCP VM startup script for {resource_prefix} deployment
#
# This script runs automatically when the VM starts (including after preemption).
# It installs Docker, mounts the persistent disk, pulls images from GHCR, and starts containers.
#
# All output is logged to /var/log/startup-script.log for debugging.

set -euo pipefail

# Log all output to startup script log
exec > >(tee -a /var/log/startup-script.log)
exec 2>&1

# Colors for output
RED='\\033[0;31m'
GREEN='\\033[0;32m'
YELLOW='\\033[1;33m'
NC='\\033[0m' # No Color

echo_info() {{ echo -e "${{GREEN}}[INFO]${{NC}} $(date '+%Y-%m-%d %H:%M:%S') $1"; }}
echo_warn() {{ echo -e "${{YELLOW}}[WARN]${{NC}} $(date '+%Y-%m-%d %H:%M:%S') $1"; }}
echo_error() {{ echo -e "${{RED}}[ERROR]${{NC}} $(date '+%Y-%m-%d %H:%M:%S') $1"; }}

echo_info "Starting VM setup..."

# ============================================================================
# Install Docker and Docker Compose
# ============================================================================
echo_info "Installing Docker and Docker Compose..."
apt-get update
apt-get install -y ca-certificates curl git

# Install Docker using official convenience script
echo_info "Installing Docker via official script..."
curl -fsSL https://get.docker.com | sh

# Enable and start Docker service
echo_info "Enabling Docker service..."
systemctl enable docker
systemctl start docker

# Wait for Docker to be ready
echo_info "Waiting for Docker to be ready..."
for i in {{1..10}}; do
    if docker info >/dev/null 2>&1; then
        echo_info "Docker is ready"
        break
    fi
    if [[ $i -eq 10 ]]; then
        echo_error "Docker failed to start after 10 attempts"
        exit 1
    fi
    echo_warn "Waiting for Docker... ($i/10)"
    sleep 2
done

# ============================================================================
# Mount persistent disk
# ============================================================================
DEVICE_NAME="/dev/disk/by-id/google-data-disk"
MOUNT_POINT="/mnt/disks/data"

echo_info "Mounting persistent disk..."
mkdir -p "$MOUNT_POINT"

# Format only if first boot (no existing filesystem)
if ! blkid "$DEVICE_NAME" >/dev/null 2>&1; then
    echo_info "Formatting disk (first boot)..."
    mkfs.ext4 -F "$DEVICE_NAME"
else
    echo_info "Disk already formatted, skipping format step"
fi

# Mount disk
if ! mountpoint -q "$MOUNT_POINT"; then
    echo_info "Mounting disk to $MOUNT_POINT..."
    mount -o discard,defaults "$DEVICE_NAME" "$MOUNT_POINT"
else
    echo_info "Disk already mounted at $MOUNT_POINT"
fi

# Add to fstab for automatic remount on restart (if not already present)
if ! grep -q "$DEVICE_NAME" /etc/fstab; then
    echo_info "Adding disk to /etc/fstab for automatic mounting..."
    echo "$DEVICE_NAME $MOUNT_POINT ext4 discard,defaults,nofail 0 2" >> /etc/fstab
else
    echo_info "Disk already in /etc/fstab"
fi

# ============================================================================
# Setup application directory
# ============================================================================
APP_DIR="/opt/{resource_prefix}"
echo_info "Creating application directory: $APP_DIR"
mkdir -p "$APP_DIR"
cd "$APP_DIR"

# ============================================================================
# Setup data directory with correct permissions
# ============================================================================
DATA_DIR="$MOUNT_POINT/jobs"
echo_info "Setting up data directory: $DATA_DIR"
mkdir -p "$DATA_DIR"
# Ensure UID 999 (appuser in container) can write
chown -R 999:999 "$MOUNT_POINT" 2>/dev/null || true
chmod -R 755 "$MOUNT_POINT" 2>/dev/null || true

# ============================================================================
# Detect CPU count and create .env file with resource limits
# ============================================================================
echo_info "Detecting CPU count for resource configuration..."
CPU_COUNT=$(nproc)
echo_info "Detected $CPU_COUNT CPU cores"

# Create .env file with computed memory limit and runtime CPU detection
cat > .env <<EOF
WORKER_MEMORY_LIMIT={worker_memory_limit}
NWCHEM_NPROC=$(nproc)
EOF
echo_info "Created .env with WORKER_MEMORY_LIMIT={worker_memory_limit}, NWCHEM_NPROC=$CPU_COUNT"

# ============================================================================
# Download docker-compose.yml from GitHub
# ============================================================================
echo_info "Downloading docker-compose.yml from GitHub..."
REPO_URL="https://raw.githubusercontent.com/steinbeck/qm-nmr-calc/master"

if ! curl -fsSL "${{REPO_URL}}/docker-compose.yml" -o docker-compose.yml; then
    echo_error "Failed to download docker-compose.yml"
    exit 1
fi

echo_info "Downloaded docker-compose.yml"

# ============================================================================
# Create HTTP-only docker-compose.gcp.yml override (no Caddy, port 80)
# ============================================================================
echo_info "Creating docker-compose.gcp.yml with GCP-specific overrides..."

cat > docker-compose.gcp.yml <<'EOF'
# Docker Compose GCP overrides for {resource_prefix}
# HTTP-only deployment with persistent disk storage (no Caddy/HTTPS)
#
# Usage: docker compose -f docker-compose.yml -f docker-compose.gcp.yml up -d

services:
  init:
    volumes:
      - /mnt/disks/data:/app/data

  api:
    volumes:
      - /mnt/disks/data:/app/data
    environment:
      - ENVIRONMENT=production
      - LOG_LEVEL=info
    ports:
      - "80:8000"

  worker:
    volumes:
      - /mnt/disks/data:/app/data
    environment:
      - ENVIRONMENT=production
      - LOG_LEVEL=info
      - NWCHEM_NPROC=${{NWCHEM_NPROC}}
    deploy:
      resources:
        limits:
          memory: ${{WORKER_MEMORY_LIMIT}}

  # Disable Caddy — move to inactive profile so it never starts
  # (HTTP-only deployment, API serves directly on port 80)
  caddy:
    profiles:
      - caddy-enabled
EOF

echo_info "Created docker-compose.gcp.yml (HTTP-only, no Caddy)"

# ============================================================================
# Pull images from GHCR
# ============================================================================
echo_info "Pulling Docker images from GHCR..."
if ! docker compose pull; then
    echo_error "Failed to pull Docker images"
    exit 1
fi

echo_info "Images pulled successfully"

# ============================================================================
# Start services
# ============================================================================
echo_info "Starting Docker Compose services..."
if ! docker compose -f docker-compose.yml -f docker-compose.gcp.yml up -d; then
    echo_error "Failed to start Docker Compose services"
    exit 1
fi

echo_info "Services started successfully"

# ============================================================================
# Show running containers
# ============================================================================
echo_info "Running containers:"
docker compose -f docker-compose.yml -f docker-compose.gcp.yml ps

# ============================================================================
# Wait for services to be healthy
# ============================================================================
echo_info "Waiting for services to become healthy..."

MAX_WAIT=120
WAIT_INTERVAL=5
WAITED=0

while [[ $WAITED -lt $MAX_WAIT ]]; do
    # Check if API is healthy
    API_HEALTH=$(docker compose -f docker-compose.yml -f docker-compose.gcp.yml ps api --format "{{{{.Health}}}}" 2>/dev/null || echo "unknown")

    if [[ "$API_HEALTH" == "healthy" ]]; then
        echo_info "API service is healthy!"
        break
    fi

    echo_info "Waiting for API to become healthy... ($WAITED/${{MAX_WAIT}}s)"
    sleep $WAIT_INTERVAL
    WAITED=$((WAITED + WAIT_INTERVAL))
done

if [[ $WAITED -ge $MAX_WAIT ]]; then
    echo_warn "API did not become healthy within ${{MAX_WAIT}}s - check logs with: docker compose logs api"
fi

# ============================================================================
# Final status
# ============================================================================
echo ""
echo "=============================================="
echo_info "VM Setup Complete!"
echo "=============================================="
echo ""
echo "HTTP endpoint: http://<VM_IP>"
echo "Data dir:      $MOUNT_POINT"
echo ""
echo "Container status:"
docker compose -f docker-compose.yml -f docker-compose.gcp.yml ps --format "table {{{{.Service}}}}\\t{{{{.Status}}}}\\t{{{{.Health}}}}"
echo ""
echo_info "Startup script finished at $(date '+%Y-%m-%d %H:%M:%S')"
"""


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main() -> None:
    """CLI entry point for machine selection and startup script generation."""
    parser = argparse.ArgumentParser(
        description="Select GCP machine type and calculate resources"
    )
    parser.add_argument(
        "--cpu-cores",
        type=int,
        default=8,
        help="Minimum CPU cores (default: 8)",
    )
    parser.add_argument(
        "--ram-gb",
        type=int,
        default=32,
        help="Minimum RAM in GB (default: 32)",
    )
    parser.add_argument(
        "--zone",
        type=str,
        help="Specific zone to use (optional, auto-selects cheapest if omitted)",
    )
    parser.add_argument(
        "--generate-startup-script",
        action="store_true",
        help="Generate startup script instead of JSON output",
    )
    parser.add_argument(
        "--resource-prefix",
        type=str,
        default="qm-nmr-calc",
        help="Resource name prefix (default: qm-nmr-calc)",
    )
    parser.add_argument(
        "--disk-size-gb",
        type=int,
        default=100,
        help="Persistent disk size in GB (default: 100)",
    )

    args = parser.parse_args()

    # Find available zone and machine type
    if args.zone:
        machine_type = select_machine_type(args.cpu_cores, args.ram_gb, args.zone)
        zone_info = {"zone": args.zone, "machine_type": machine_type}
    else:
        zone_info = find_available_zone(args.cpu_cores, args.ram_gb)

    # Calculate Docker resources
    resources = calculate_docker_resources(zone_info["machine_type"], zone_info["zone"])

    if args.generate_startup_script:
        # Output startup script
        script = generate_startup_script(
            resources["worker_memory_limit"],
            args.resource_prefix,
            args.disk_size_gb,
        )
        print(script)
    else:
        # Output JSON
        output = {
            "zone": zone_info["zone"],
            "region": zone_info.get("region", zone_info["zone"].rsplit("-", 1)[0]),
            "machine_type": zone_info["machine_type"],
            "resources": resources,
        }
        print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
