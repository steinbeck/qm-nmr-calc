# Roadmap: qm-nmr-calc

## Milestones

- [x] **v1.0 MVP** - Phases 1-6 (shipped 2026-01-20)
- [x] **v1.1 Accurate Chemical Shifts** - Phases 7-11.2 (shipped 2026-01-25)
- [x] **v2.0 Conformational Sampling** - Phases 12-17 (shipped 2026-01-28)
- [x] **v2.0.1 Conformer Pre-selection** - Phase 24 (shipped 2026-01-30)
- [x] **v2.1 UI Redesign** - Phases 18-23 (shipped 2026-01-31)
- [x] **v2.2 Documentation** - Phases 25-31 (shipped 2026-02-01)
- [x] **v2.3 NMReData Export** - Phases 32-34 (shipped 2026-02-01)
- [x] **v2.4 Docker Deployment** - Phases 35-40 (shipped 2026-02-03)
- [x] **v2.5 ARM64 Docker Support** - Phases 41-44 (shipped 2026-02-04)
- [x] **v2.6 Google Cloud Spot Deployment** - Phases 45-48.1 (shipped 2026-02-05)
- [ ] **v2.7 Automated GCP Deployment** - Phases 49-53 (in progress)

## Overview

**Current milestone:** v2.7 Automated GCP Deployment

v2.7 replaces v2.6's interactive deployment scripts with a fully automated, config-driven system that finds the cheapest spot instance across all GCP regions and deploys end-to-end without manual intervention. HTTP-only deployment pattern for fire-up-and-burn usage. Also fixes conformer progress tracking display bug.

## Phases

<details>
<summary>v1.0 MVP (Phases 1-6) - SHIPPED 2026-01-20</summary>

### Phase 1: Project Setup and Foundation
**Goal**: Runnable service with basic architecture
**Plans**: 3 plans
**Status**: Complete

### Phase 2: NWChem Integration
**Goal**: Execute NWChem calculations via Huey
**Plans**: 3 plans
**Status**: Complete

### Phase 3: NMR Calculation Pipeline
**Goal**: End-to-end NMR predictions
**Plans**: 2 plans
**Status**: Complete

### Phase 4: Web Interface
**Goal**: Usable web UI for job submission
**Plans**: 2 plans
**Status**: Complete

### Phase 5: Result Visualization
**Goal**: Visual outputs (spectrum, structure)
**Plans**: 2 plans
**Status**: Complete

### Phase 6: Polish and Documentation
**Goal**: Production-ready v1.0
**Plans**: 2 plans
**Status**: Complete

</details>

<details>
<summary>v1.1 Accurate Chemical Shifts (Phases 7-11.2) - SHIPPED 2026-01-25</summary>

### Phase 7: Remove ISiCLE Dependency
**Goal**: Direct NWChem integration
**Plans**: 2 plans
**Status**: Complete

### Phase 8: COSMO Solvation
**Goal**: Accurate solvent modeling
**Plans**: 2 plans
**Status**: Complete

### Phase 9: DELTA50 Benchmark
**Goal**: NWChem-specific scaling factors
**Plans**: 3 plans
**Status**: Complete

### Phase 10: Enhanced Visualization
**Goal**: Interactive 3D viewer with annotations
**Plans**: 2 plans
**Status**: Complete

### Phase 11: Geometry Input Options
**Goal**: Accept pre-optimized structures
**Plans**: 2 plans
**Status**: Complete

### Phase 11.1: Solvent Expansion Research
**Goal**: Research additional solvents
**Plans**: 1 plan
**Status**: Complete

### Phase 11.2: Vacuum Benchmark
**Goal**: Gas-phase scaling factors
**Plans**: 3 plans
**Status**: Complete

</details>

<details>
<summary>v2.0 Conformational Sampling (Phases 12-17) - SHIPPED 2026-01-28</summary>

### Phase 12: Conformer Data Model and Storage
**Goal**: Foundation for multi-conformer calculations
**Plans**: 3 plans
**Status**: Complete

### Phase 13: RDKit KDG Conformer Generation
**Goal**: Generate conformer ensembles using RDKit distance geometry
**Plans**: 3 plans
**Status**: Complete

### Phase 14: Boltzmann Averaging Implementation
**Goal**: Numerically stable Boltzmann weighting and ensemble averaging
**Plans**: 2 plans
**Status**: Complete

### Phase 15: Multi-Conformer NWChem Integration
**Goal**: Sequential DFT optimization and NMR calculation loop
**Plans**: 3 plans
**Status**: Complete

### Phase 16: CREST Integration (Optional High-Accuracy Mode)
**Goal**: Production-quality conformer generation using CREST metadynamics
**Plans**: 3 plans
**Status**: Complete

### Phase 17: API Integration and Progress Tracking
**Goal**: User-facing ensemble mode with conformer selection and progress
**Plans**: 5 plans
**Status**: Complete

</details>

<details>
<summary>v2.0.1 Conformer Pre-selection Hotfix - SHIPPED 2026-01-30</summary>

### Phase 24: Conformer Pre-selection with xTB
**Goal**: Efficient conformer filtering using RMSD clustering and xTB ranking
**Plans**: 3 plans
**Status**: Complete

</details>

<details>
<summary>v2.1 UI Redesign (Phases 18-23) - SHIPPED 2026-01-31</summary>

### Phase 18: CSS Foundation and Design System
**Goal**: Reusable CSS architecture with design tokens
**Plans**: 4 plans
**Status**: Complete

### Phase 19: Results Page Redesign
**Goal**: Bento grid layout with prominent visualizations
**Plans**: 3 plans
**Status**: Complete

### Phase 20: Submit Page Redesign
**Goal**: Clean form layout with logical grouping
**Plans**: 2 plans
**Status**: Complete

### Phase 21: Status Page Redesign
**Goal**: Job progress visualization with step tracking
**Plans**: 2 plans
**Status**: Complete

### Phase 22: Responsive and Layout Polish
**Goal**: Mobile-first breakpoints with optimized glass effects
**Plans**: 2 plans
**Status**: Complete

### Phase 23: Accessibility and Testing
**Goal**: WCAG compliance and cross-browser validation
**Plans**: 2 plans
**Status**: Complete

</details>

<details>
<summary>v2.2 Documentation (Phases 25-31) - SHIPPED 2026-02-01</summary>

### Phase 25: README and Documentation Structure
**Goal**: Overhaul README with architecture overview
**Plans**: 1 plan
**Status**: Complete

### Phase 26: Installation Guide
**Goal**: Comprehensive installation documentation
**Plans**: 1 plan
**Status**: Complete

### Phase 27: Usage Guide
**Goal**: User-facing documentation for web UI and API
**Plans**: 2 plans
**Status**: Complete

### Phase 28: Technical Architecture
**Goal**: Developer-facing system architecture docs
**Plans**: 2 plans
**Status**: Complete

### Phase 29: Library Documentation
**Goal**: Third-party library integrations documentation
**Plans**: 2 plans
**Status**: Complete

### Phase 30: DP4+ Science Documentation
**Goal**: Scientific writeup with full derivations
**Plans**: 2 plans
**Status**: Complete

### Phase 31: Documentation Polish
**Goal**: Consistency and cross-linking review
**Plans**: 1 plan
**Status**: Complete

</details>

<details>
<summary>v2.3 NMReData Export (Phases 32-34) - SHIPPED 2026-02-01</summary>

### Phase 32: Core NMReData Generation Module
**Goal**: Generate NMReData-compliant SDF files
**Plans**: 1 plan
**Status**: Complete

### Phase 33: API and UI Integration
**Goal**: NMReData download via REST API and web UI
**Plans**: 1 plan
**Status**: Complete

### Phase 34: Testing and Validation
**Goal**: Format compliance and round-trip testing
**Plans**: 1 plan
**Status**: Complete

</details>

<details>
<summary>v2.4 Docker Deployment (Phases 35-40) - SHIPPED 2026-02-03</summary>

### Phase 35: Worker Container
**Goal**: NWChem, CREST, and xTB run correctly in a Docker container with all environment configuration.
**Plans**: 1 plan
**Status**: Complete

### Phase 36: API Container
**Goal**: FastAPI application runs in a minimal, secure container with health checks.
**Plans**: 1 plan
**Status**: Complete

### Phase 37: Docker Compose Integration
**Goal**: Complete deployment with single `docker compose up -d` command, persistent data, and operational controls.
**Plans**: 2 plans
**Status**: Complete

### Phase 38: Caddy + HTTPS
**Goal**: Production-ready HTTPS with automatic certificate management via Let's Encrypt.
**Plans**: 1 plan
**Status**: Complete

### Phase 39: CI/CD + GHCR Publishing
**Goal**: Pre-built images available on GitHub Container Registry, built automatically on release.
**Plans**: 1 plan
**Status**: Complete

### Phase 40: Documentation
**Goal**: Users can deploy qm-nmr-calc in 5 minutes with clear guidance for production setup and troubleshooting.
**Plans**: 2 plans
**Status**: Complete

</details>

<details>
<summary>v2.5 ARM64 Docker Support (Phases 41-44) - SHIPPED 2026-02-04</summary>

### Phase 41: ARM64 Dockerfile Creation
**Goal**: ARM64 worker container exists with all computational chemistry packages installed via conda-forge.
**Plans**: 1 plan
**Status**: Complete

### Phase 42: Local Validation
**Goal**: ARM64 worker container passes all computational chemistry tests on Apple Silicon.
**Plans**: 1 plan
**Status**: Complete

### Phase 43: CI/CD Integration
**Goal**: GitHub Actions builds and publishes multi-arch images with native ARM64 runner.
**Plans**: 1 plan
**Status**: Complete

### Phase 44: Documentation and Release
**Goal**: Users know ARM64 is supported and can deploy without architecture-specific instructions.
**Plans**: 1 plan
**Status**: Complete

</details>

<details>
<summary>v2.6 Google Cloud Spot Deployment (Phases 45-48.1) - SHIPPED 2026-02-05</summary>

### Phase 45: GCP Infrastructure Setup
**Goal**: GCP project configured with networking and storage prerequisites for Spot VM deployment.
**Depends on**: Nothing (first phase of v2.6)
**Requirements**: INFRA-01, INFRA-02, INFRA-03, INFRA-04
**Success Criteria** (what must be TRUE):
  1. Static external IP is reserved and can be used for DNS configuration
  2. Firewall rules allow HTTP (80), HTTPS (443), and SSH (22) traffic
  3. Persistent disk exists for job data and Let's Encrypt certificates
  4. Persistent disk survives VM deletion (verified by delete/recreate cycle)
**Plans**: 1 plan
**Status**: Complete

Plans:
- [x] 45-01-PLAN.md - Infrastructure setup and teardown scripts

### Phase 46: VM Deployment Script
**Goal**: Single script creates a fully-configured Spot VM running qm-nmr-calc with HTTPS.
**Depends on**: Phase 45
**Requirements**: DEPLOY-01, DEPLOY-02, DEPLOY-03, DEPLOY-04, DEPLOY-05, DEPLOY-06, DEPLOY-07
**Success Criteria** (what must be TRUE):
  1. User runs one command and gets a working qm-nmr-calc instance
  2. Startup script installs Docker, pulls GHCR images, and starts containers
  3. Containers shut down gracefully within 25 seconds during preemption
  4. User can select region, zone, and machine type with sensible defaults
  5. Cost estimate displayed before VM creation (user can cancel)
**Plans**: 1 plan
**Status**: Complete

Plans:
- [x] 46-01-PLAN.md - VM deployment with startup/shutdown scripts and Docker Compose override

### Phase 47: Lifecycle Management Scripts
**Goal**: Users can stop, start, check status, and access their GCP VM without memorizing gcloud commands.
**Depends on**: Phase 46
**Requirements**: LIFE-01, LIFE-02, LIFE-03, LIFE-04, LIFE-05, LIFE-06, LIFE-07
**Success Criteria** (what must be TRUE):
  1. Stop command halts VM (stops compute billing, preserves data)
  2. Start command resumes VM (services auto-start via startup script)
  3. Delete command removes VM but preserves persistent disk
  4. Status command shows VM state, IP address, and running containers
  5. SSH and logs commands provide easy access for debugging
**Plans**: 1 plan
**Status**: Complete

Plans:
- [x] 47-01-PLAN.md - Lifecycle management scripts (stop, start, delete, status, ssh, logs)

### Phase 48: Documentation and Testing
**Goal**: Users can deploy to GCP with clear guidance on prerequisites, costs, and limitations.
**Depends on**: Phase 47
**Requirements**: DOCS-01, DOCS-02, DOCS-03, DOCS-04, DOCS-05
**Success Criteria** (what must be TRUE):
  1. README includes GCP deployment as an option alongside Docker Compose
  2. Prerequisites documented (GCP account, gcloud CLI, domain for HTTPS)
  3. Cost estimates documented (spot vs on-demand, static IP charges)
  4. Preemption behavior documented (jobs in progress will be lost)
  5. DNS configuration guide covers common providers (Cloudflare, Namecheap)
**Plans**: 1 plan
**Status**: Complete

Plans:
- [x] 48-01-PLAN.md - GCP deployment documentation

### Phase 48.1: Machine Info Display (INSERTED)
**Goal**: Web UI displays machine information (type, CPU cores, memory) for debugging and user awareness.
**Depends on**: Phase 48
**Requirements**: INFO-01
**Success Criteria** (what must be TRUE):
  1. Results page shows machine type (if on GCP) or hostname
  2. Results page shows number of CPU cores NWChem is configured to use
  3. Results page shows memory allocation for NWChem
**Plans**: 0 plans
**Status**: Not started

Plans:
- [ ] TBD (run /gsd:plan-phase 48.1 to break down)

</details>

## v2.7 Automated GCP Deployment (In Progress)

**Milestone Goal:** Fully non-interactive GCP deployment that reads config from a file, auto-discovers the cheapest spot instance across all regions, and deploys end-to-end without manual intervention. HTTP-only. Replaces v2.6 interactive scripts. Also fixes conformer progress tracking display bug.

### Phase 49: Config Foundation and Pricing Query
**Goal**: Non-interactive deployment foundation with reliable pricing data and validated configuration.
**Depends on**: Nothing (first phase of v2.7)
**Requirements**: CFG-01, CFG-02, CFG-03, CFG-04, CFG-05, PRC-01, PRC-03, PRC-04
**Success Criteria** (what must be TRUE):
  1. User creates TOML config file specifying CPU cores, RAM, GCP project, and disk size
  2. Config validation catches errors before any GCP operations (missing project ID, impossible CPU/RAM combinations)
  3. System queries CloudPrice.net API for spot pricing across all GCP regions
  4. Pricing data cached with 24-hour TTL to avoid redundant queries
  5. Hardcoded regional fallback rankings used when pricing API unavailable
  6. All gcloud commands use --quiet and --format=json for non-interactive execution
**Plans**: 2 plans
**Status**: Complete (2026-02-06)

Plans:
- [x] 49-01-PLAN.md -- TOML config validation with Pydantic (TDD)
- [x] 49-02-PLAN.md -- Spot pricing query with API, caching, and fallback (TDD)

### Phase 50: Machine Selection and Resource Calculation
**Goal**: Correct machine type mapping and dynamic Docker resource limit calculation.
**Depends on**: Phase 49
**Requirements**: PRC-02, MCH-01, MCH-02, MCH-03, MCH-04
**Success Criteria** (what must be TRUE):
  1. System maps CPU/RAM requirements to appropriate GCP machine types via gcloud filtering
  2. Machine type availability validated in target zone before VM creation
  3. System falls back to next-cheapest region if primary zone lacks capacity
  4. Docker memory limit calculated dynamically from selected machine type (VM_RAM - 8GB for OS)
  5. NWCHEM_NPROC calculated from actual CPU count on VM host (not inside container)
  6. Startup script template generated with computed WORKER_MEMORY_LIMIT and NWCHEM_NPROC
**Plans**: 2 plans
**Status**: Not started

Plans:
- [ ] 50-01-PLAN.md -- Machine type selection and resource calculation Python module (TDD)
- [ ] 50-02-PLAN.md -- Bash library wrapper for shell integration

### Phase 51: Deployment Orchestration
**Goal**: End-to-end automated deployment with progressive feedback and error handling.
**Depends on**: Phase 50
**Requirements**: PRC-05, DEP-01, DEP-02, DEP-03, DEP-04, DEP-05, DEP-07, DEP-08, DEP-09, RPL-01, RPL-02, RPL-03
**Success Criteria** (what must be TRUE):
  1. Single command deploys end-to-end with zero interactive prompts
  2. Infrastructure operations are idempotent (create if missing, reuse if exists)
  3. Deployment progress displayed with timestamped feedback
  4. Cost estimate displayed before VM creation
  5. Dry-run mode (--dry-run) shows planned actions without executing
  6. Failed deployment cleans up orphaned resources automatically
  7. deploy-auto.sh orchestrator replaces v2.6 interactive deploy-vm.sh
**Plans**: 0 plans
**Status**: Not started

Plans:
- [ ] TBD (run /gsd:plan-phase 51 to break down)

### Phase 52: HTTP-Only Container Deployment
**Goal**: Container deployment with HTTP-only configuration and correct resource limits.
**Depends on**: Phase 51
**Requirements**: DEP-06, LCY-01, LCY-02, RPL-04
**Success Criteria** (what must be TRUE):
  1. docker-compose.gcp.yml exposes HTTP on port 80 (no Caddy HTTPS configuration)
  2. Dynamic .env file generated with computed WORKER_MEMORY_LIMIT and NWCHEM_NPROC
  3. Container startup validated (health checks pass)
  4. Existing lifecycle scripts (start, stop, delete, status, ssh, logs) continue to work
  5. Teardown script removes all created resources cleanly
  6. HTTPS/domain/Caddy TLS configuration removed from GCP deployment
**Plans**: 0 plans
**Status**: Not started

Plans:
- [ ] TBD (run /gsd:plan-phase 52 to break down)

### Phase 53: Conformer Progress Bug Fix
**Goal**: Conformer progress tracking displays correctly during processing.
**Depends on**: Nothing (independent UI fix)
**Requirements**: BUG-01, BUG-02
**Success Criteria** (what must be TRUE):
  1. Status bar shows accurate conformer count during processing (e.g., "1/2" not "0/2")
  2. Conformer statuses update correctly as conformers complete
  3. Progress bar reflects actual conformer completion state
**Plans**: 0 plans
**Status**: Not started

Plans:
- [ ] TBD (run /gsd:plan-phase 53 to break down)

## Progress

**Execution Order:**
Phases execute in numeric order: 49 -> 50 -> 51 -> 52 -> 53
(Phase 53 can run in parallel with 49-52 as it's independent UI fix)

| Phase | Milestone | Plans Complete | Status | Completed |
|-------|-----------|----------------|--------|-----------|
| 1. Project Setup | v1.0 | 3/3 | Complete | 2026-01-20 |
| 2. NWChem Integration | v1.0 | 3/3 | Complete | 2026-01-20 |
| 3. NMR Pipeline | v1.0 | 2/2 | Complete | 2026-01-20 |
| 4. Web Interface | v1.0 | 2/2 | Complete | 2026-01-20 |
| 5. Visualization | v1.0 | 2/2 | Complete | 2026-01-20 |
| 6. Polish | v1.0 | 2/2 | Complete | 2026-01-20 |
| 7. Direct NWChem | v1.1 | 2/2 | Complete | 2026-01-25 |
| 8. COSMO Solvation | v1.1 | 2/2 | Complete | 2026-01-25 |
| 9. DELTA50 Benchmark | v1.1 | 3/3 | Complete | 2026-01-25 |
| 10. 3D Visualization | v1.1 | 2/2 | Complete | 2026-01-25 |
| 11. Geometry Input | v1.1 | 2/2 | Complete | 2026-01-25 |
| 11.1 Solvent Research | v1.1 | 1/1 | Complete | 2026-01-25 |
| 11.2 Vacuum Benchmark | v1.1 | 3/3 | Complete | 2026-01-25 |
| 12. Data Model | v2.0 | 3/3 | Complete | 2026-01-26 |
| 13. RDKit Generation | v2.0 | 3/3 | Complete | 2026-01-27 |
| 14. Boltzmann Averaging | v2.0 | 2/2 | Complete | 2026-01-27 |
| 15. NWChem Integration | v2.0 | 3/3 | Complete | 2026-01-27 |
| 16. CREST Integration | v2.0 | 3/3 | Complete | 2026-01-28 |
| 17. API Integration | v2.0 | 5/5 | Complete | 2026-01-28 |
| 18. CSS Foundation | v2.1 | 4/4 | Complete | 2026-01-29 |
| 19. Results Redesign | v2.1 | 3/3 | Complete | 2026-01-29 |
| 20. Submit Redesign | v2.1 | 2/2 | Complete | 2026-01-29 |
| 21. Status Redesign | v2.1 | 2/2 | Complete | 2026-01-31 |
| 22. Responsive Polish | v2.1 | 2/2 | Complete | 2026-01-31 |
| 23. Accessibility | v2.1 | 2/2 | Complete | 2026-01-31 |
| 24. Conformer Pre-selection | v2.0.1 | 3/3 | Complete | 2026-01-30 |
| 25. README & Docs Structure | v2.2 | 1/1 | Complete | 2026-01-31 |
| 26. Installation Guide | v2.2 | 1/1 | Complete | 2026-01-31 |
| 27. Usage Guide | v2.2 | 2/2 | Complete | 2026-02-01 |
| 28. Technical Architecture | v2.2 | 2/2 | Complete | 2026-02-01 |
| 29. Library Documentation | v2.2 | 2/2 | Complete | 2026-02-01 |
| 30. DP4+ Science Documentation | v2.2 | 2/2 | Complete | 2026-02-01 |
| 31. Documentation Polish | v2.2 | 1/1 | Complete | 2026-02-01 |
| 32. NMReData Module | v2.3 | 1/1 | Complete | 2026-02-01 |
| 33. API/UI Integration | v2.3 | 1/1 | Complete | 2026-02-01 |
| 34. Testing & Validation | v2.3 | 1/1 | Complete | 2026-02-01 |
| 35. Worker Container | v2.4 | 1/1 | Complete | 2026-02-02 |
| 36. API Container | v2.4 | 1/1 | Complete | 2026-02-02 |
| 37. Docker Compose Integration | v2.4 | 2/2 | Complete | 2026-02-03 |
| 38. Caddy + HTTPS | v2.4 | 1/1 | Complete | 2026-02-03 |
| 39. CI/CD + GHCR Publishing | v2.4 | 1/1 | Complete | 2026-02-03 |
| 40. Documentation | v2.4 | 2/2 | Complete | 2026-02-03 |
| 41. ARM64 Dockerfile | v2.5 | 1/1 | Complete | 2026-02-04 |
| 42. Local Validation | v2.5 | 1/1 | Complete | 2026-02-04 |
| 43. CI/CD Integration | v2.5 | 1/1 | Complete | 2026-02-04 |
| 44. Documentation | v2.5 | 1/1 | Complete | 2026-02-04 |
| 45. GCP Infrastructure | v2.6 | 1/1 | Complete | 2026-02-04 |
| 46. VM Deployment | v2.6 | 1/1 | Complete | 2026-02-04 |
| 47. Lifecycle Scripts | v2.6 | 1/1 | Complete | 2026-02-05 |
| 48. Documentation | v2.6 | 1/1 | Complete | 2026-02-05 |
| 48.1 Machine Info | v2.6 | 0/TBD | Not started | - |
| **49. Config & Pricing** | **v2.7** | **2/2** | **Complete** | **2026-02-06** |
| **50. Machine Selection** | **v2.7** | **0/2** | **Not started** | - |
| **51. Orchestration** | **v2.7** | **0/TBD** | **Not started** | - |
| **52. HTTP Container** | **v2.7** | **0/TBD** | **Not started** | - |
| **53. Conformer Bug Fix** | **v2.7** | **0/TBD** | **Not started** | - |

## Coverage (v2.6)

**v2.6 Requirements: 23 total**

| Requirement | Phase | Description |
|-------------|-------|-------------|
| INFRA-01 | 45 | Static external IP for stable DNS |
| INFRA-02 | 45 | Firewall rules (HTTP, HTTPS, SSH) |
| INFRA-03 | 45 | Persistent disk for job data and certificates |
| INFRA-04 | 45 | Persistent disk survives VM deletion |
| DEPLOY-01 | 46 | One-command VM creation with Spot configuration |
| DEPLOY-02 | 46 | Startup script installs Docker and deploys containers |
| DEPLOY-03 | 46 | Startup script pulls images from GHCR |
| DEPLOY-04 | 46 | Graceful container shutdown during preemption |
| DEPLOY-05 | 46 | Interactive prompts for region/zone/machine type |
| DEPLOY-06 | 46 | Cost estimation before VM creation |
| DEPLOY-07 | 46 | docker-compose.gcp.yml override for GCP settings |
| LIFE-01 | 47 | Stop command halts VM |
| LIFE-02 | 47 | Start command resumes VM |
| LIFE-03 | 47 | Delete command removes VM (preserves disk) |
| LIFE-04 | 47 | Status command shows VM state and IP |
| LIFE-05 | 47 | SSH command provides shell access |
| LIFE-06 | 47 | Logs command streams container logs |
| LIFE-07 | 47 | Configuration persistence between commands |
| DOCS-01 | 48 | README section on GCP deployment |
| DOCS-02 | 48 | Prerequisites documented |
| DOCS-03 | 48 | Cost estimates documented |
| DOCS-04 | 48 | Preemption limitations documented |
| DOCS-05 | 48 | DNS configuration guide |

**Mapped: 23/23 (100%)**

## Coverage (v2.7)

**v2.7 Requirements: 31 total**

| Requirement | Phase | Description |
|-------------|-------|-------------|
| CFG-01 | 49 | CPU/RAM in TOML config |
| CFG-02 | 49 | GCP project and prefix in config |
| CFG-03 | 49 | Persistent disk size in config |
| CFG-04 | 49 | Config validation with Pydantic |
| CFG-05 | 49 | Clear error messages for invalid config |
| PRC-01 | 49 | Query spot pricing across all regions |
| PRC-02 | 50 | Select cheapest region/zone for spec |
| PRC-03 | 49 | Pricing data cached with TTL |
| PRC-04 | 49 | Hardcoded fallback when API unavailable |
| PRC-05 | 51 | Cost estimate before deployment |
| MCH-01 | 50 | Map CPU/RAM to machine types |
| MCH-02 | 50 | Validate machine type availability |
| MCH-03 | 50 | Fallback to next-cheapest region |
| MCH-04 | 50 | Dynamic Docker resource limits |
| DEP-01 | 51 | Single command, zero prompts |
| DEP-02 | 51 | Non-interactive gcloud commands |
| DEP-03 | 51 | Idempotent infrastructure operations |
| DEP-04 | 51 | VM created as Spot with correct config |
| DEP-05 | 51 | Dynamic startup script generation |
| DEP-06 | 52 | HTTP-only on port 80 |
| DEP-07 | 51 | Timestamped progress feedback |
| DEP-08 | 51 | Auto-cleanup on failed deployment |
| DEP-09 | 51 | Dry-run mode |
| LCY-01 | 52 | Existing lifecycle scripts work |
| LCY-02 | 52 | Teardown script removes all resources |
| RPL-01 | 51 | Replaces v2.6 interactive scripts |
| RPL-02 | 51 | No interactive prompts |
| RPL-03 | 51 | DNS check removed |
| RPL-04 | 52 | HTTPS/Caddy removed |
| BUG-01 | 53 | Conformer progress updates correctly |
| BUG-02 | 53 | Status bar shows accurate count |

**Mapped: 31/31 (100%)**

---
*Last updated: 2026-02-06 - Phase 50 planned (2 plans in 2 waves)*
