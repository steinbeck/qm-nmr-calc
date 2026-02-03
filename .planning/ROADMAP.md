# Roadmap: qm-nmr-calc

## Milestones

- [x] **v1.0 MVP** - Phases 1-6 (shipped 2026-01-20)
- [x] **v1.1 Accurate Chemical Shifts** - Phases 7-11.2 (shipped 2026-01-25)
- [x] **v2.0 Conformational Sampling** - Phases 12-17 (shipped 2026-01-28)
- [x] **v2.0.1 Conformer Pre-selection** - Phase 24 (shipped 2026-01-30)
- [x] **v2.1 UI Redesign** - Phases 18-23 (shipped 2026-01-31)
- [x] **v2.2 Documentation** - Phases 25-31 (shipped 2026-02-01)
- [x] **v2.3 NMReData Export** - Phases 32-34 (shipped 2026-02-01)
- [ ] **v2.4 Docker Deployment** - Phases 35-40 (in progress)

## Overview

**Current milestone:** v2.4 Docker Deployment

Docker deployment transforms qm-nmr-calc from a manual installation into a `docker compose up` experience. The roadmap progresses from the most complex, blocking component (worker container with NWChem/CREST/xTB) through orchestration and HTTPS, culminating in CI/CD publishing and documentation. Each phase delivers a testable increment toward production-ready containerization.

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

<details open>
<summary>v2.4 Docker Deployment (Phases 35-40) - IN PROGRESS</summary>

**Milestone Goal:** Make qm-nmr-calc deployable with `docker compose up` -- pre-built images, auto-HTTPS, production-ready defaults.

### Phase 35: Worker Container
**Goal**: NWChem, CREST, and xTB run correctly in a Docker container with all environment configuration.
**Depends on**: Nothing (first phase of v2.4)
**Requirements**: DOCK-02
**Success Criteria** (what must be TRUE):
  1. User can build worker image with `docker build -f Dockerfile.worker`
  2. NWChem DFT calculation completes successfully inside container
  3. CREST conformer search completes successfully inside container
  4. xTB energy calculation completes successfully inside container
  5. Huey consumer runs and processes queued jobs in container
**Plans**: 1 plan
**Status**: Complete

Plans:
- [x] 35-01-PLAN.md -- Create and validate worker container with NWChem, CREST, xTB, and Huey

### Phase 36: API Container
**Goal**: FastAPI application runs in a minimal, secure container with health checks.
**Depends on**: Nothing (parallel with Phase 35)
**Requirements**: DOCK-03
**Success Criteria** (what must be TRUE):
  1. User can build API image with `docker build -f Dockerfile.api`
  2. FastAPI server responds to HTTP requests at port 8000
  3. Health check endpoint returns 200 OK
  4. Container runs as non-root user
**Plans**: 1 plan
**Status**: Complete

Plans:
- [x] 36-01-PLAN.md -- Create Dockerfile.api with multi-stage build and validation script

### Phase 37: Docker Compose Integration
**Goal**: Complete deployment with single `docker compose up -d` command, persistent data, and operational controls.
**Depends on**: Phase 35, Phase 36
**Requirements**: DOCK-01, DOCK-04, DOCK-05, DOCK-06, DOCK-07, DOCK-08, DOCK-09, OPS-01, OPS-02, OPS-03, OPS-04
**Success Criteria** (what must be TRUE):
  1. User can start entire stack with `docker compose up -d`
  2. Job data persists after `docker compose down && docker compose up -d`
  3. Huey queue state persists across container restarts
  4. All services restart automatically after simulated failure
  5. User can configure deployment via `.env` file with documented options
  6. Worker completes current job on SIGTERM before stopping
  7. User can view logs with `docker compose logs`
**Plans**: 2 plans
**Status**: Complete

Plans:
- [x] 37-01-PLAN.md -- Create docker-compose.yml and .env.example configuration
- [x] 37-02-PLAN.md -- Validation script and integration testing

### Phase 38: Caddy + HTTPS
**Goal**: Production-ready HTTPS with automatic certificate management via Let's Encrypt.
**Depends on**: Phase 37
**Requirements**: HTTPS-01, HTTPS-02, HTTPS-03, HTTPS-04, HTTPS-05
**Success Criteria** (what must be TRUE):
  1. Caddy reverse proxy serves application on ports 80 and 443
  2. HTTPS certificate obtained automatically when domain is configured
  3. HTTP requests redirect to HTTPS automatically
  4. User can set domain via `DOMAIN` environment variable
  5. Deployment works on localhost without domain (HTTP-only mode)
**Plans**: 1 plan
**Status**: Complete

Plans:
- [x] 38-01-PLAN.md -- Add Caddy reverse proxy with conditional HTTPS

### Phase 39: CI/CD + GHCR Publishing
**Goal**: Pre-built images available on GitHub Container Registry, built automatically on release.
**Depends on**: Phase 38
**Requirements**: GHCR-01, GHCR-02, GHCR-03, GHCR-04
**Success Criteria** (what must be TRUE):
  1. Pre-built images exist on ghcr.io/[owner]/qm-nmr-calc
  2. GitHub Actions workflow triggers on release tag
  3. Images tagged with version (e.g., v2.4.0) and `latest`
  4. Images support both amd64 and arm64 architectures
**Plans**: 1 plan
**Status**: Complete

Plans:
- [x] 39-01-PLAN.md -- Create GitHub Actions workflow for GHCR publishing with OCI labels

### Phase 40: Documentation
**Goal**: Users can deploy qm-nmr-calc in 5 minutes with clear guidance for production setup and troubleshooting.
**Depends on**: Phase 39
**Requirements**: DOCS-01, DOCS-02, DOCS-03, DOCS-04
**Success Criteria** (what must be TRUE):
  1. README contains quick start section with 5-minute deployment instructions
  2. Deployment guide covers cloud VPS setup (DigitalOcean, Linode, etc.)
  3. Troubleshooting section addresses common issues from research pitfalls
  4. Backup and restore instructions document volume management
**Plans**: TBD

Plans:
- [ ] 40-01: TBD

</details>

## Progress

**Execution Order:**
Phases execute in numeric order: 35 -> 36 -> 37 -> 38 -> 39 -> 40

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
| **35. Worker Container** | **v2.4** | **1/1** | **Complete** | 2026-02-02 |
| **36. API Container** | **v2.4** | **1/1** | **Complete** | 2026-02-02 |
| **37. Docker Compose Integration** | **v2.4** | **2/2** | **Complete** | 2026-02-03 |
| **38. Caddy + HTTPS** | **v2.4** | **1/1** | **Complete** | 2026-02-03 |
| **39. CI/CD + GHCR Publishing** | **v2.4** | **1/1** | **Complete** | 2026-02-03 |
| **40. Documentation** | **v2.4** | **0/TBD** | **Not started** | - |

## Coverage (v2.4)

**v2.4 Requirements: 26 total**

| Requirement | Phase | Description |
|-------------|-------|-------------|
| DOCK-01 | 37 | Deploy with `docker compose up -d` |
| DOCK-02 | 35 | Worker includes NWChem, CREST, xTB |
| DOCK-03 | 36 | App container runs FastAPI |
| DOCK-04 | 37 | Job data persists via volume |
| DOCK-05 | 37 | Huey queue persists via volume |
| DOCK-06 | 37 | Health check configuration |
| DOCK-07 | 37 | Restart policies |
| DOCK-08 | 37 | `.env` file configuration |
| DOCK-09 | 37 | `.env.example` documentation |
| HTTPS-01 | 38 | Caddy on ports 80/443 |
| HTTPS-02 | 38 | Auto-HTTPS via Let's Encrypt |
| HTTPS-03 | 38 | HTTP-to-HTTPS redirect |
| HTTPS-04 | 38 | DOMAIN env var |
| HTTPS-05 | 38 | Localhost HTTP mode |
| GHCR-01 | 39 | Images on GHCR |
| GHCR-02 | 39 | GitHub Actions build |
| GHCR-03 | 39 | Version + latest tags |
| GHCR-04 | 39 | Multi-arch (amd64 + arm64) |
| OPS-01 | 37 | Graceful SIGTERM handling |
| OPS-02 | 37 | Worker memory limits |
| OPS-03 | 37 | NWCHEM_NPROC env var |
| OPS-04 | 37 | `docker compose logs` |
| DOCS-01 | 40 | Quick start (5-minute) |
| DOCS-02 | 40 | Cloud VPS guide |
| DOCS-03 | 40 | Troubleshooting |
| DOCS-04 | 40 | Backup/restore |

**Mapped: 26/26 (100%)**

---
*Last updated: 2026-02-03 - Phase 39 complete*
