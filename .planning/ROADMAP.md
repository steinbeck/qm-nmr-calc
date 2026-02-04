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
- [ ] **v2.5 ARM64 Docker Support** - Phases 41-44 (in progress)

## Overview

**Current milestone:** v2.5 ARM64 Docker Support

v2.5 enables native ARM64 execution for the worker container, unlocking local development on Apple Silicon Macs and deployment to ARM-based cloud instances (AWS Graviton). The roadmap progresses from creating the ARM64-specific Dockerfile using conda-forge packages, through local validation on Apple Silicon, to CI/CD integration with native ARM64 runners, culminating in documentation and release. Each phase builds on the previous, with the critical dependency being a working Dockerfile before any validation or automation can proceed.

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

<details open>
<summary>v2.5 ARM64 Docker Support (Phases 41-44) - IN PROGRESS</summary>

**Milestone Goal:** Native ARM64 support for the worker container, enabling local development on Apple Silicon Macs and deployment to ARM-based cloud instances without emulation.

### Phase 41: ARM64 Dockerfile Creation
**Goal**: ARM64 worker container exists with all computational chemistry packages installed via conda-forge.
**Depends on**: Nothing (first phase of v2.5)
**Requirements**: CONT-05
**Success Criteria** (what must be TRUE):
  1. Dockerfile.worker.arm64 exists and builds successfully on ARM64 host
  2. Container uses micromamba base image with conda-forge packages (not x86 binaries)
  3. NWChem, xTB, and CREST are installed from conda-forge linux-aarch64 channel
  4. Environment variables configured for NWChem basis sets and OpenBLAS threading
**Plans**: 1 plan
**Status**: Complete

Plans:
- [x] 41-01-PLAN.md -- Create Dockerfile.worker.arm64, env-worker-arm64.yaml, and validation script

### Phase 42: Local Validation
**Goal**: ARM64 worker container passes all computational chemistry tests on Apple Silicon.
**Depends on**: Phase 41
**Requirements**: CONT-01, CONT-02, CONT-03, CONT-04, CONT-06
**Success Criteria** (what must be TRUE):
  1. NWChem DFT geometry optimization completes without SIGILL or crashes
  2. NWChem NMR shielding calculation produces valid chemical shifts
  3. xTB energy calculation completes successfully
  4. CREST conformer search finds multiple conformers
  5. Results match x86_64 output within tolerance (0.5 ppm 1H, 2.0 ppm 13C)
  6. Full NMR prediction pipeline produces results matching x86 within tolerance
**Plans**: 1 plan
**Status**: Not started

Plans:
- [ ] 42-01-PLAN.md -- Create validation script and test fixtures, build and validate ARM64 container on Apple Silicon

### Phase 43: CI/CD Integration
**Goal**: GitHub Actions builds and publishes multi-arch images with native ARM64 runner.
**Depends on**: Phase 42
**Requirements**: CICD-01, CICD-02, CICD-03, CICD-04
**Success Criteria** (what must be TRUE):
  1. GitHub Actions workflow builds ARM64 worker image on ubuntu-24.04-arm runner
  2. Multi-arch manifest created combining amd64 and arm64 images
  3. Single image tag (latest, vX.Y.Z) pulls correct architecture automatically
  4. ARM64 build completes in reasonable time (under 15 minutes, not QEMU-slow)
**Plans**: TBD
**Status**: Not started

Plans:
- [ ] 43-01-PLAN.md -- Update publish-images.yml with ARM64 build job and manifest merge

### Phase 44: Documentation and Release
**Goal**: Users know ARM64 is supported and can deploy without architecture-specific instructions.
**Depends on**: Phase 43
**Requirements**: DOCS-01, DOCS-02, DOCS-03
**Success Criteria** (what must be TRUE):
  1. README mentions ARM64/Apple Silicon support in Docker deployment section
  2. Known issues section documents any ARM64-specific caveats discovered during validation
  3. docker-compose.yml works unchanged on ARM64 (auto-pulls correct architecture)
**Plans**: TBD
**Status**: Not started

Plans:
- [ ] 44-01-PLAN.md -- Update README and deployment docs with ARM64 support

</details>

## Progress

**Execution Order:**
Phases execute in numeric order: 41 -> 42 -> 43 -> 44

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
| **41. ARM64 Dockerfile** | **v2.5** | **1/1** | **Complete** | 2026-02-04 |
| **42. Local Validation** | **v2.5** | **0/1** | **Not started** | - |
| **43. CI/CD Integration** | **v2.5** | **0/1** | **Not started** | - |
| **44. Documentation** | **v2.5** | **0/1** | **Not started** | - |

## Coverage (v2.5)

**v2.5 Requirements: 13 total**

| Requirement | Phase | Description |
|-------------|-------|-------------|
| CONT-01 | 42 | ARM64 worker runs NWChem DFT geometry optimization |
| CONT-02 | 42 | ARM64 worker runs NWChem NMR shielding calculation |
| CONT-03 | 42 | ARM64 worker runs xTB energy calculations |
| CONT-04 | 42 | ARM64 worker runs CREST conformer search |
| CONT-05 | 41 | ARM64 worker uses conda-forge packages |
| CONT-06 | 42 | ARM64 worker produces numerically equivalent results |
| CICD-01 | 43 | GitHub Actions builds ARM64 worker with native runner |
| CICD-02 | 43 | Multi-arch manifest combining amd64 and arm64 |
| CICD-03 | 43 | Single image tag works on both architectures |
| CICD-04 | 43 | ARM64 build uses native runner (not QEMU) |
| DOCS-01 | 44 | README documents ARM64/Apple Silicon support |
| DOCS-02 | 44 | Known issues section covers ARM64-specific caveats |
| DOCS-03 | 44 | docker-compose.yml works unchanged on ARM64 |

**Mapped: 13/13 (100%)**

---
*Last updated: 2026-02-04 - Phase 42 plan created*
