# Requirements: qm-nmr-calc v2.6

**Defined:** 2026-02-04
**Core Value:** Reliable async NMR predictions with full control over calculation parameters

## v2.6 Requirements

Requirements for GCP Spot VM deployment. Enables cost-effective cloud deployment.

### Infrastructure

- [x] **INFRA-01**: Script creates static external IP for stable DNS
- [x] **INFRA-02**: Script creates firewall rules (HTTP 80, HTTPS 443, SSH 22)
- [x] **INFRA-03**: Script creates persistent disk for job data and certificates
- [x] **INFRA-04**: Persistent disk survives VM deletion/recreation

### Deployment

- [x] **DEPLOY-01**: One-command VM creation with Spot configuration
- [x] **DEPLOY-02**: Startup script installs Docker and deploys containers
- [x] **DEPLOY-03**: Startup script pulls images from GHCR
- [x] **DEPLOY-04**: Graceful container shutdown during preemption (25s timeout)
- [x] **DEPLOY-05**: Interactive prompts for region, zone, machine type with defaults
- [x] **DEPLOY-06**: Cost estimation displayed before VM creation
- [x] **DEPLOY-07**: docker-compose.gcp.yml override for GCP-specific settings

### Lifecycle

- [ ] **LIFE-01**: Stop command halts VM (preserves data, stops billing)
- [ ] **LIFE-02**: Start command resumes stopped VM
- [ ] **LIFE-03**: Delete command removes VM (preserves persistent disk)
- [ ] **LIFE-04**: Status command shows VM state and IP address
- [ ] **LIFE-05**: SSH command provides shell access to VM
- [ ] **LIFE-06**: Logs command streams container logs
- [ ] **LIFE-07**: Configuration persistence (remembers VM name, zone between commands)

### Documentation

- [ ] **DOCS-01**: README section on GCP deployment option
- [ ] **DOCS-02**: Prerequisites documented (GCP account, gcloud CLI, domain)
- [ ] **DOCS-03**: Cost estimates documented (spot vs on-demand)
- [ ] **DOCS-04**: Preemption limitations documented (job loss on interrupt)
- [ ] **DOCS-05**: DNS configuration guide for common providers

## Future Requirements

Deferred to later milestones.

### Preemption Resilience

- **PREEMPT-01**: NWChem restart file persistence to survive preemption
- **PREEMPT-02**: Automatic job retry after preemption
- **PREEMPT-03**: Managed Instance Group for auto-recreation

### Multi-Cloud

- **CLOUD-01**: AWS Spot Instance support
- **CLOUD-02**: Azure Spot VM support

## Out of Scope

Explicitly excluded. Documented to prevent scope creep.

| Feature | Reason |
|---------|--------|
| Terraform/IaC | Overkill for single VM, manual lifecycle |
| Auto-scaling | Manual lifecycle, user controls when to run |
| Load balancer | Single VM, Caddy handles HTTPS |
| Cloud Monitoring integration | docker logs sufficient for this scope |
| Custom VM images | Startup script is simpler to maintain |
| Job checkpoint/restart | Complex, accept job loss for v2.6 |

## Traceability

Which phases cover which requirements. Updated during roadmap creation.

| Requirement | Phase | Status |
|-------------|-------|--------|
| INFRA-01 | Phase 45 | Complete |
| INFRA-02 | Phase 45 | Complete |
| INFRA-03 | Phase 45 | Complete |
| INFRA-04 | Phase 45 | Complete |
| DEPLOY-01 | Phase 46 | Complete |
| DEPLOY-02 | Phase 46 | Complete |
| DEPLOY-03 | Phase 46 | Complete |
| DEPLOY-04 | Phase 46 | Complete |
| DEPLOY-05 | Phase 46 | Complete |
| DEPLOY-06 | Phase 46 | Complete |
| DEPLOY-07 | Phase 46 | Complete |
| LIFE-01 | Phase 47 | Pending |
| LIFE-02 | Phase 47 | Pending |
| LIFE-03 | Phase 47 | Pending |
| LIFE-04 | Phase 47 | Pending |
| LIFE-05 | Phase 47 | Pending |
| LIFE-06 | Phase 47 | Pending |
| LIFE-07 | Phase 47 | Pending |
| DOCS-01 | Phase 48 | Pending |
| DOCS-02 | Phase 48 | Pending |
| DOCS-03 | Phase 48 | Pending |
| DOCS-04 | Phase 48 | Pending |
| DOCS-05 | Phase 48 | Pending |

**Coverage:**
- v2.6 requirements: 23 total
- Mapped to phases: 23
- Unmapped: 0

---
*Requirements defined: 2026-02-04*
*Last updated: 2026-02-04 after Phase 46 complete*
