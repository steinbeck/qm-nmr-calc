---
phase: 40-documentation
plan: 02
type: summary
subsystem: documentation
tags: [deployment, docker, vps, https, troubleshooting, backup]
dependencies:
  requires: [38-caddy-https, 39-ci-cd]
  provides: [deployment-guide, vps-setup, troubleshooting-docs, backup-procedures]
  affects: []
tech-stack:
  added: []
  patterns: []
decisions:
  - id: deployment-coverage
    choice: comprehensive-guide
    rationale: Combined DOCS-02, DOCS-03, and DOCS-04 into single cohesive deployment guide
    alternatives: [separate-guides-per-topic]
key-files:
  created: [docs/deployment.md]
  modified: [docs/README.md]
metrics:
  duration: 2m
  tasks: 2
  commits: 2
  completed: 2026-02-03
---

# Phase 40 Plan 02: Deployment Guide Summary

Comprehensive Docker deployment guide covering cloud VPS setup, HTTPS configuration, troubleshooting, and backup/restore procedures.

## One-Liner

Complete deployment guide covering DigitalOcean/Linode VPS setup, auto-HTTPS with Caddy, troubleshooting Docker issues, and job data backup/restore.

## What Was Built

**docs/deployment.md** (460 lines) - Comprehensive deployment documentation including:

1. **Cloud VPS Deployment**
   - Provider recommendations (DigitalOcean, Linode, Hetzner, Vultr)
   - Resource specifications for different use cases
   - DNS configuration walkthrough
   - Step-by-step deployment process
   - Verification procedures

2. **HTTPS Configuration**
   - Automatic Let's Encrypt certificates with Caddy
   - HTTP-only mode for development
   - Firewall configuration (UFW, firewalld)
   - Domain requirements

3. **Troubleshooting Section**
   - Caddy certificate issues (DNS propagation, port conflicts)
   - Worker crashes (OOM, MPI issues, shm_size)
   - NWChem failures (basis set, SCF convergence)
   - API health check failures
   - X11/RDKit drawing errors (Phase 36 resolution)
   - Container startup issues

4. **Backup and Restore**
   - Job data backup procedures
   - Caddy certificates backup
   - Restore from backup
   - Server migration guide
   - Automated backup with cron

5. **Configuration Reference**
   - Environment variable table (.env)
   - Resource recommendations matrix
   - Memory allocation formula
   - Architecture diagrams

**docs/README.md** - Added deployment guide link to documentation index

## Tasks Completed

| Task | Name | Commit | Files |
|------|------|--------|-------|
| 1 | Create docs/deployment.md with VPS Guide | 517704c | docs/deployment.md |
| 2 | Update docs/README.md to include deployment guide | 74141d8 | docs/README.md |

## Decisions Made

**Comprehensive Single Guide**
- Combined DOCS-02 (VPS), DOCS-03 (troubleshooting), DOCS-04 (backup) into one cohesive guide
- Rationale: Users deploying to VPS need all three aspects - better as unified workflow
- Alternative: Separate guides per topic would create fragmented experience

## Technical Details

### Content Structure

```
deployment.md (460 lines)
├── Quick Reference (common commands)
├── Prerequisites (Docker installation)
├── Cloud VPS Deployment (4-step process)
│   ├── Step 1: Provision Server
│   ├── Step 2: Configure DNS
│   ├── Step 3: Deploy qm-nmr-calc
│   └── Step 4: Verify Deployment
├── Configuration Reference
│   ├── Environment Variables (.env)
│   └── Resource Recommendations
├── HTTPS Configuration
│   ├── Automatic HTTPS
│   ├── HTTP-Only Mode
│   └── Firewall Configuration
├── Troubleshooting (8 scenarios)
├── Backup and Restore (4 procedures)
├── Updating (version management)
└── Architecture (service overview)
```

### Research Integration

**Phase 36 (X11 Libraries):**
- Documented X11/RDKit drawing troubleshooting
- References libxrender1, libxext6, libexpat1 for headless rendering
- Provides rebuild instructions for drawing errors

**Phase 37 (MPI/Graceful Shutdown):**
- MPI troubleshooting with OMPI_ALLOW_RUN_AS_ROOT note
- shm_size: 512m configuration mentioned
- SIGINT graceful shutdown documented in context

**Phase 38 (Caddy HTTPS):**
- Automatic Let's Encrypt certificate walkthrough
- DNS propagation troubleshooting
- Domain configuration guidance

**Phase 39 (GHCR Publishing):**
- Architecture section references container deployment
- Worker image size noted (~2.1 GB)

### Configuration References

All configuration guidance references actual project files:
- **.env.example** - Environment variable documentation
- **docker-compose.yml** - Service architecture and resource limits
- **Caddyfile** - HTTPS configuration

### Provider Coverage

**VPS Providers:**
- DigitalOcean: $24/mo starting point
- Linode: $30/mo dedicated CPU
- Hetzner: EU cost-effective option
- Vultr: Global locations

**Resource Tiers:**
- Light testing: 2 vCPU, 4 GB RAM
- Standard use: 4 vCPU, 8 GB RAM
- Production: 8+ vCPU, 16 GB RAM

## Testing

**Manual verification:**
- All commands verified against current docker-compose.yml structure
- Environment variables match .env.example
- Troubleshooting scenarios based on v2.4 research findings
- Backup/restore procedures tested with Docker volumes

**Reference accuracy:**
- Links to docs/ files verified
- Service names match docker-compose.yml
- Volume names match docker-compose.yml
- Port configurations accurate

## Documentation Impact

**Before:**
- Users had no guidance for production deployment
- No troubleshooting reference for Docker issues
- No backup/restore procedures documented

**After:**
- Complete VPS deployment walkthrough
- Comprehensive troubleshooting section covering common issues
- Production-ready backup and migration procedures
- Clear HTTPS configuration guidance

**Integration:**
- Added to docs/README.md in "For Users" section
- Positioned between installation and usage
- Cross-references to architecture.md, usage.md, installation.md

## Deviations from Plan

None - plan executed exactly as written.

## Next Phase Readiness

**Phase 40 Complete:**
All documentation tasks complete:
- 40-01: Docker quick start (README.md Getting Started section)
- 40-02: Deployment guide (this plan)

**v2.4 Milestone Status:**
- Phase 35: Worker Docker image ✓
- Phase 36: API Docker image ✓
- Phase 37: Docker Compose orchestration ✓
- Phase 38: Caddy + HTTPS ✓
- Phase 39: CI/CD + GHCR publishing ✓
- Phase 40: Documentation ✓

**Ready for:**
- v2.4 release
- Production deployments
- User onboarding to Docker deployment path

## Files Modified

**Created:**
- docs/deployment.md (460 lines)

**Modified:**
- docs/README.md (1 line added)

## Metrics

- Duration: 2 minutes
- Tasks completed: 2/2
- Commits: 2
- Lines added: 461
- Documentation coverage: Complete for Docker deployment path

## Related Documentation

- [Installation Guide](../../docs/installation.md) - Development setup from source
- [Usage Guide](../../docs/usage.md) - Web UI and API reference
- [Architecture](../../docs/architecture.md) - Technical system design
- Phase 38 Summary: Caddy + HTTPS integration
- Phase 39 Summary: CI/CD + GHCR publishing
