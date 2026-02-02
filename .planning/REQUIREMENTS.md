# Requirements: qm-nmr-calc v2.4 Docker Deployment

**Defined:** 2026-02-02
**Core Value:** Reliable async NMR predictions with full control over calculation parameters

## v2.4 Requirements

Requirements for Docker deployment milestone. Each maps to roadmap phases.

### Container Stack

- [ ] **DOCK-01**: User can deploy app with `docker compose up -d` command
- [ ] **DOCK-02**: Worker container includes NWChem, CREST, and xTB pre-installed
- [ ] **DOCK-03**: App container runs FastAPI server with all Python dependencies
- [ ] **DOCK-04**: Job data persists across container restarts via named volume
- [ ] **DOCK-05**: Huey queue database persists across container restarts
- [ ] **DOCK-06**: All services have health check configuration
- [ ] **DOCK-07**: All services restart automatically on failure (unless-stopped policy)
- [ ] **DOCK-08**: User can configure deployment via `.env` file
- [ ] **DOCK-09**: `.env.example` documents all configuration options with comments

### Reverse Proxy & HTTPS

- [ ] **HTTPS-01**: Caddy reverse proxy serves app on ports 80/443
- [ ] **HTTPS-02**: HTTPS certificates obtained automatically via Let's Encrypt
- [ ] **HTTPS-03**: HTTP requests redirect to HTTPS automatically
- [ ] **HTTPS-04**: User can configure domain via `DOMAIN` environment variable
- [ ] **HTTPS-05**: Deployment works on localhost without domain (HTTP mode)

### Image Publishing

- [ ] **GHCR-01**: Pre-built images available on GitHub Container Registry
- [ ] **GHCR-02**: GitHub Actions workflow builds and publishes on release
- [ ] **GHCR-03**: Images tagged with version (e.g., `v2.4.0`) and `latest`
- [ ] **GHCR-04**: Multi-architecture builds support amd64 and arm64

### Operations

- [ ] **OPS-01**: Worker handles SIGTERM gracefully (completes current job)
- [ ] **OPS-02**: Worker container has memory limits configured
- [ ] **OPS-03**: User can configure NWChem parallel processes via environment variable
- [ ] **OPS-04**: User can view all service logs with `docker compose logs`

### Documentation

- [ ] **DOCS-01**: Quick start section in README (5-minute deployment)
- [ ] **DOCS-02**: Deployment guide for cloud VPS setup
- [ ] **DOCS-03**: Troubleshooting section for common issues
- [ ] **DOCS-04**: Backup and restore instructions for job data

## Future Requirements

Deferred to post-v2.4. Tracked but not in current roadmap.

### Kubernetes

- **K8S-01**: Kubernetes manifest for deployment
- **K8S-02**: Helm chart for customizable deployment

### Monitoring

- **MON-01**: Prometheus metrics endpoint
- **MON-02**: Grafana dashboard template

### Security

- **SEC-01**: Docker secrets integration for sensitive config
- **SEC-02**: Rate limiting at reverse proxy level

## Out of Scope

Explicitly excluded. Documented to prevent scope creep.

| Feature | Reason |
|---------|--------|
| Kubernetes manifests | Overkill for target users (single VM academic deployments) |
| Redis/message broker | Huey SQLite is adequate for single-node |
| Built-in monitoring stack | Grafana/Prometheus adds complexity; use external tools if needed |
| Horizontal autoscaling | Single VM deployment; manual worker count adjustment |
| Docker Swarm mode | Single node doesn't need orchestration |
| Multi-user authentication | Separate milestone if needed |
| Database container | Filesystem-based storage works well |
| Singularity/Apptainer | HPC container format; defer to post-v2.4 |

## Traceability

Which phases cover which requirements. Updated during roadmap creation.

| Requirement | Phase | Status |
|-------------|-------|--------|
| DOCK-01 | Phase 37 | Pending |
| DOCK-02 | Phase 35 | Pending |
| DOCK-03 | Phase 36 | Pending |
| DOCK-04 | Phase 37 | Pending |
| DOCK-05 | Phase 37 | Pending |
| DOCK-06 | Phase 37 | Pending |
| DOCK-07 | Phase 37 | Pending |
| DOCK-08 | Phase 37 | Pending |
| DOCK-09 | Phase 37 | Pending |
| HTTPS-01 | Phase 38 | Pending |
| HTTPS-02 | Phase 38 | Pending |
| HTTPS-03 | Phase 38 | Pending |
| HTTPS-04 | Phase 38 | Pending |
| HTTPS-05 | Phase 38 | Pending |
| GHCR-01 | Phase 39 | Pending |
| GHCR-02 | Phase 39 | Pending |
| GHCR-03 | Phase 39 | Pending |
| GHCR-04 | Phase 39 | Pending |
| OPS-01 | Phase 37 | Pending |
| OPS-02 | Phase 37 | Pending |
| OPS-03 | Phase 37 | Pending |
| OPS-04 | Phase 37 | Pending |
| DOCS-01 | Phase 40 | Pending |
| DOCS-02 | Phase 40 | Pending |
| DOCS-03 | Phase 40 | Pending |
| DOCS-04 | Phase 40 | Pending |

**Coverage:**
- v2.4 requirements: 26 total
- Mapped to phases: 26
- Unmapped: 0

---
*Requirements defined: 2026-02-02*
*Last updated: 2026-02-02 after roadmap creation*
