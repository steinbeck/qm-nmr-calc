# Requirements: qm-nmr-calc v2.6

**Defined:** 2026-02-04
**Core Value:** Reliable async NMR predictions with full control over calculation parameters

## v2.6 Requirements

Requirements for GCP Spot VM deployment. Enables cost-effective cloud deployment.

### Infrastructure

- [ ] **INFRA-01**: Script creates static external IP for stable DNS
- [ ] **INFRA-02**: Script creates firewall rules (HTTP 80, HTTPS 443, SSH 22)
- [ ] **INFRA-03**: Script creates persistent disk for job data and certificates
- [ ] **INFRA-04**: Persistent disk survives VM deletion/recreation

### Deployment

- [ ] **DEPLOY-01**: One-command VM creation with Spot configuration
- [ ] **DEPLOY-02**: Startup script installs Docker and deploys containers
- [ ] **DEPLOY-03**: Startup script pulls images from GHCR
- [ ] **DEPLOY-04**: Graceful container shutdown during preemption (25s timeout)
- [ ] **DEPLOY-05**: Interactive prompts for region, zone, machine type with defaults
- [ ] **DEPLOY-06**: Cost estimation displayed before VM creation
- [ ] **DEPLOY-07**: docker-compose.gcp.yml override for GCP-specific settings

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
| INFRA-01 | TBD | Pending |
| INFRA-02 | TBD | Pending |
| INFRA-03 | TBD | Pending |
| INFRA-04 | TBD | Pending |
| DEPLOY-01 | TBD | Pending |
| DEPLOY-02 | TBD | Pending |
| DEPLOY-03 | TBD | Pending |
| DEPLOY-04 | TBD | Pending |
| DEPLOY-05 | TBD | Pending |
| DEPLOY-06 | TBD | Pending |
| DEPLOY-07 | TBD | Pending |
| LIFE-01 | TBD | Pending |
| LIFE-02 | TBD | Pending |
| LIFE-03 | TBD | Pending |
| LIFE-04 | TBD | Pending |
| LIFE-05 | TBD | Pending |
| LIFE-06 | TBD | Pending |
| LIFE-07 | TBD | Pending |
| DOCS-01 | TBD | Pending |
| DOCS-02 | TBD | Pending |
| DOCS-03 | TBD | Pending |
| DOCS-04 | TBD | Pending |
| DOCS-05 | TBD | Pending |

**Coverage:**
- v2.6 requirements: 23 total
- Mapped to phases: 0
- Unmapped: 23 (pending roadmap creation)

---
*Requirements defined: 2026-02-04*
*Last updated: 2026-02-04 after initial definition*
