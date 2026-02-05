#!/bin/bash
# GCP VM startup script for qm-nmr-calc deployment
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
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo_info() { echo -e "${GREEN}[INFO]${NC} $(date '+%Y-%m-%d %H:%M:%S') $1"; }
echo_warn() { echo -e "${YELLOW}[WARN]${NC} $(date '+%Y-%m-%d %H:%M:%S') $1"; }
echo_error() { echo -e "${RED}[ERROR]${NC} $(date '+%Y-%m-%d %H:%M:%S') $1"; }

echo_info "Starting VM setup..."

# ============================================================================
# Install Docker and Docker Compose
# ============================================================================
echo_info "Installing Docker and Docker Compose..."
apt-get update
apt-get install -y docker.io docker-compose-v2 curl git

# Enable and start Docker service
echo_info "Enabling Docker service..."
systemctl enable docker
systemctl start docker

# Wait for Docker to be ready
echo_info "Waiting for Docker to be ready..."
for i in {1..10}; do
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
APP_DIR="/opt/qm-nmr-calc"
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
# Fetch DOMAIN from instance metadata
# ============================================================================
echo_info "Fetching DOMAIN from instance metadata..."
DOMAIN=$(curl -s "http://metadata.google.internal/computeMetadata/v1/instance/attributes/DOMAIN" -H "Metadata-Flavor: Google" || true)
if [[ -z "$DOMAIN" ]]; then
    echo_warn "DOMAIN not set in metadata, using 'localhost'"
    DOMAIN="localhost"
fi
echo_info "Domain: $DOMAIN"

# Create .env file with DOMAIN
echo "DOMAIN=${DOMAIN}" > .env

# ============================================================================
# Download docker-compose.yml and Caddyfile from GitHub
# ============================================================================
echo_info "Downloading docker-compose.yml and Caddyfile from GitHub..."
REPO_URL="https://raw.githubusercontent.com/Steinbeck-Lab/qm-nmr-calc/master"

if ! curl -fsSL "${REPO_URL}/docker-compose.yml" -o docker-compose.yml; then
    echo_error "Failed to download docker-compose.yml"
    exit 1
fi

if ! curl -fsSL "${REPO_URL}/Caddyfile" -o Caddyfile; then
    echo_error "Failed to download Caddyfile"
    exit 1
fi

echo_info "Downloaded docker-compose.yml and Caddyfile"

# ============================================================================
# Create docker-compose.gcp.yml override
# ============================================================================
echo_info "Creating docker-compose.gcp.yml with GCP-specific overrides..."

cat > docker-compose.gcp.yml <<'EOF'
# Docker Compose GCP overrides for qm-nmr-calc
# Uses persistent disk for data storage and sets production environment

services:
  api:
    volumes:
      # Override to use mounted persistent disk instead of named volume
      - /mnt/disks/data:/app/data
    environment:
      - ENVIRONMENT=production
      - LOG_LEVEL=info

  worker:
    volumes:
      # Override to use mounted persistent disk instead of named volume
      - /mnt/disks/data:/app/data
    environment:
      - ENVIRONMENT=production
      - LOG_LEVEL=info
      # Let auto-detection handle CPU count based on VM size
      - NWCHEM_NPROC=${NWCHEM_NPROC:-}

  # Caddy DOMAIN comes from .env file (already set in base compose)
  # Named volumes (caddy_data, caddy_config) remain as Docker volumes for certificate storage
EOF

echo_info "Created docker-compose.gcp.yml"

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
    API_HEALTH=$(docker compose -f docker-compose.yml -f docker-compose.gcp.yml ps api --format "{{.Health}}" 2>/dev/null || echo "unknown")

    if [[ "$API_HEALTH" == "healthy" ]]; then
        echo_info "API service is healthy!"
        break
    fi

    echo_info "Waiting for API to become healthy... ($WAITED/${MAX_WAIT}s)"
    sleep $WAIT_INTERVAL
    WAITED=$((WAITED + WAIT_INTERVAL))
done

if [[ $WAITED -ge $MAX_WAIT ]]; then
    echo_warn "API did not become healthy within ${MAX_WAIT}s - check logs with: docker compose logs api"
fi

# ============================================================================
# Final status
# ============================================================================
echo ""
echo "=============================================="
echo_info "VM Setup Complete!"
echo "=============================================="
echo ""
echo "Domain:     https://${DOMAIN}"
echo "Data dir:   $MOUNT_POINT"
echo ""
echo "Container status:"
docker compose -f docker-compose.yml -f docker-compose.gcp.yml ps --format "table {{.Service}}\t{{.Status}}\t{{.Health}}"
echo ""
echo_info "Startup script finished at $(date '+%Y-%m-%d %H:%M:%S')"
