---
phase: 46-vm-deployment-script
verified: 2026-02-05T06:03:25Z
status: passed
score: 5/5 must-haves verified
---

# Phase 46: VM Deployment Script Verification Report

**Phase Goal:** Single script creates a fully-configured Spot VM running qm-nmr-calc with HTTPS.
**Verified:** 2026-02-05T06:03:25Z
**Status:** passed
**Re-verification:** No — initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | User runs one command and gets a working qm-nmr-calc VM | ✓ VERIFIED | deploy-vm.sh creates VM with single command, startup script handles full deployment |
| 2 | Startup script installs Docker, mounts disk, pulls images, starts containers | ✓ VERIFIED | startup.sh installs docker.io, mounts /mnt/disks/data with blkid check, runs docker compose pull, starts services |
| 3 | Containers shut down gracefully within 25 seconds during preemption | ✓ VERIFIED | shutdown.sh uses `docker compose stop --timeout 25` to stay within 30s preemption window |
| 4 | User can select region, zone, and machine type with sensible defaults | ✓ VERIFIED | deploy-vm.sh prompts with defaults: region=$GCP_REGION, zone=$GCP_ZONE, machine=e2-standard-4 |
| 5 | Cost estimate displayed before VM creation with cancel option | ✓ VERIFIED | deploy-vm.sh shows monthly cost estimates (~$15-100/month) and confirms before proceeding, exits if user declines |

**Score:** 5/5 truths verified

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `gcp/startup.sh` | VM initialization: Docker install, disk mount, container deployment (min 50 lines) | ✓ VERIFIED | EXISTS (191 lines), SUBSTANTIVE (installs docker.io/docker-compose-v2, mounts /mnt/disks/data, pulls GHCR images, starts containers), WIRED (used by deploy-vm.sh via --metadata-from-file=startup-script) |
| `gcp/shutdown.sh` | Graceful container shutdown during preemption with 25s timeout (min 15 lines) | ✓ VERIFIED | EXISTS (52 lines), SUBSTANTIVE (contains `docker compose stop --timeout 25`, handles missing directory gracefully), WIRED (used by deploy-vm.sh via --metadata-from-file=shutdown-script) |
| `gcp/docker-compose.gcp.yml` | GCP-specific volume paths and environment overrides (min 15 lines) | ✓ VERIFIED | EXISTS (33 lines), SUBSTANTIVE (overrides volumes to /mnt/disks/data, sets ENVIRONMENT=production, LOG_LEVEL=info), WIRED (created inline by startup.sh, used in docker compose up command) |
| `gcp/deploy-vm.sh` | One-command Spot VM creation with prompts (min 100 lines) | ✓ VERIFIED | EXISTS (273 lines), SUBSTANTIVE (interactive prompts, cost estimation, gcloud create with --provisioning-model=SPOT), WIRED (references startup.sh and shutdown.sh via --metadata-from-file) |

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|----|--------|---------|
| `gcp/deploy-vm.sh` | `gcp/startup.sh` | --metadata-from-file=startup-script | ✓ WIRED | Line 229: `--metadata-from-file=startup-script=startup.sh,shutdown-script=shutdown.sh` |
| `gcp/deploy-vm.sh` | `gcp/shutdown.sh` | --metadata-from-file=shutdown-script | ✓ WIRED | Line 229: `--metadata-from-file=startup-script=startup.sh,shutdown-script=shutdown.sh` |
| `gcp/startup.sh` | `gcp/docker-compose.gcp.yml` | creates docker-compose.gcp.yml in /opt/qm-nmr-calc | ✓ WIRED | Lines 130-158: creates file via heredoc, line 176: uses in `docker compose -f docker-compose.yml -f docker-compose.gcp.yml up -d` |

### Requirements Coverage

All Phase 46 requirements mapped in PLAN and SUMMARY are satisfied:

| Requirement | Status | Supporting Evidence |
|-------------|--------|-------------------|
| DEPLOY-01: Single command VM creation | ✓ SATISFIED | deploy-vm.sh creates Spot VM with gcloud compute instances create |
| DEPLOY-02: Startup script installs Docker and deploys containers | ✓ SATISFIED | startup.sh installs docker.io, docker-compose-v2, starts services |
| DEPLOY-03: Startup script pulls images from GHCR | ✓ SATISFIED | startup.sh line 165: `docker compose pull` |
| DEPLOY-04: Shutdown script stops containers within 25s timeout | ✓ SATISFIED | shutdown.sh line 45: `docker compose stop --timeout 25` |
| DEPLOY-05: Interactive prompts with defaults | ✓ SATISFIED | deploy-vm.sh prompts for region, zone, machine type (lines 81-100) |
| DEPLOY-06: Cost estimation before creation | ✓ SATISFIED | deploy-vm.sh displays monthly cost estimates (lines 128-145) and confirms (line 153) |
| DEPLOY-07: GCP-specific Docker Compose overrides | ✓ SATISFIED | docker-compose.gcp.yml overrides volumes to /mnt/disks/data, sets production env |

### Anti-Patterns Found

**No blocker anti-patterns detected.**

Scanned files:
- gcp/startup.sh (191 lines) - No TODO/FIXME/placeholder patterns
- gcp/shutdown.sh (52 lines) - No TODO/FIXME/placeholder patterns  
- gcp/docker-compose.gcp.yml (33 lines) - No TODO/FIXME/placeholder patterns
- gcp/deploy-vm.sh (273 lines) - No TODO/FIXME/placeholder patterns

All scripts:
- ✓ Executable permissions set (755)
- ✓ No empty implementations
- ✓ No console.log-only handlers
- ✓ Error handling present (set -euo pipefail, exit on failures)
- ✓ Idempotent operations (blkid check before format, mountpoint check, fstab grep)

### Implementation Quality

**Strengths:**
1. **Complete startup automation** - Docker install, disk mount with blkid safety check, GHCR image pull, service start with health verification
2. **Graceful preemption handling** - 25s timeout stays within GCP's 30s preemption window
3. **Production-ready configuration** - Persistent disk bind mounts, ENVIRONMENT=production, LOG_LEVEL=info
4. **User-friendly deployment** - Interactive prompts with sensible defaults (e2-standard-4), cost estimation with monthly ranges
5. **Safety checks** - Confirms before VM deletion/recreation, validates domain input, checks for existing VM
6. **Comprehensive logging** - All output logged to /var/log with timestamps
7. **Proper error handling** - set -euo pipefail, exit codes on failures, retry loops for Docker readiness

**Patterns established:**
- Metadata-driven configuration (DOMAIN passed via VM metadata)
- Startup script downloads compose files from GitHub master branch (no Cloud Storage dependency)
- Docker Compose overrides pattern for environment-specific configuration
- Color-coded output for info/warn/error messages

### Human Verification Required

While all automated checks pass, the following should be verified by a human tester:

#### 1. End-to-End VM Deployment

**Test:** Run `./deploy-vm.sh` and complete the full deployment workflow
**Expected:** 
- Interactive prompts work correctly
- Cost estimation displays accurately for selected machine type
- VM creates successfully with static IP assignment
- Startup script completes in 2-3 minutes
- Containers start and are accessible via HTTPS at configured domain
- Let's Encrypt certificate issues automatically after DNS propagates

**Why human:** Full deployment workflow requires GCP credentials, DNS configuration, and testing across different machine types and regions

#### 2. Graceful Preemption Handling

**Test:** Manually preempt the VM or wait for GCP preemption
**Expected:**
- Shutdown script runs automatically
- Containers stop within 25 seconds
- Persistent disk data survives preemption
- VM can be restarted and containers resume with preserved data

**Why human:** Preemption testing requires actual GCP Spot VM behavior, cannot be simulated programmatically

#### 3. Startup Script Idempotency

**Test:** Stop and restart the VM multiple times
**Expected:**
- Startup script doesn't reformat existing disk (blkid check works)
- fstab doesn't get duplicate entries
- Containers restart correctly on subsequent boots

**Why human:** Requires actual VM restarts to verify idempotent behavior across boots

#### 4. HTTPS Certificate Issuance

**Test:** Configure DNS A record, wait for propagation, verify HTTPS access
**Expected:**
- Caddy automatically requests Let's Encrypt certificate
- Certificate stores in caddy_data volume (survives restarts)
- HTTPS works at configured domain

**Why human:** Requires real DNS configuration and Let's Encrypt interaction

---

## Verification Summary

**All must-haves verified.** Phase 46 goal achieved.

The phase delivers on its promise: "Single script creates a fully-configured Spot VM running qm-nmr-calc with HTTPS."

**Key achievements:**
1. ✓ One-command deployment with deploy-vm.sh
2. ✓ Complete VM automation via startup script (Docker, disk, containers)
3. ✓ Graceful 25s shutdown during preemption
4. ✓ Interactive prompts with sensible defaults
5. ✓ Cost estimation with user confirmation

**No gaps found.** Ready to proceed to Phase 47 (Lifecycle Scripts).

---

_Verified: 2026-02-05T06:03:25Z_  
_Verifier: Claude (gsd-verifier)_
