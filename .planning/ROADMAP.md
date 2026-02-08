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
- [x] **v2.7 Automated GCP Deployment** - Phases 49-53 (shipped 2026-02-06)
- [ ] **v2.8 Expanded Solvent Support** - Phases 54-58 (in progress)

## Overview

**Current milestone:** v2.8 Expanded Solvent Support

Add 4 new NMR solvents (methanol, water, acetone, benzene) with DELTA50-derived B3LYP scaling factors. This extends solvent coverage from 3 to 7, following the same benchmark-derive-integrate pipeline established in v1.1 but now applied to 4 solvents in parallel. The compute-intensive benchmark phase (200 NWChem calculations) dominates the timeline.

- [x] **Phase 54: Benchmark Infrastructure** - Extend CLI and input_gen to accept 4 new solvents (2026-02-07)
- [ ] **Phase 55: DELTA50 Benchmark Calculations** - Run 200 NWChem calculations (50 molecules x 4 solvents)
- [ ] **Phase 56: Scaling Factor Derivation** - Derive OLS factors and validate quality gates
- [ ] **Phase 57: Solvent Integration** - Wire all 4 solvents into solvents.py, shifts.py, and UI/API
- [ ] **Phase 58: Documentation** - Update SCALING-FACTORS.md and README for 7-solvent support

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
**Plans**: 1 plan
**Status**: Complete

### Phase 46: VM Deployment Script
**Goal**: Single script creates a fully-configured Spot VM running qm-nmr-calc with HTTPS.
**Plans**: 1 plan
**Status**: Complete

### Phase 47: Lifecycle Management Scripts
**Goal**: Users can stop, start, check status, and access their GCP VM without memorizing gcloud commands.
**Plans**: 1 plan
**Status**: Complete

### Phase 48: Documentation and Testing
**Goal**: Users can deploy to GCP with clear guidance on prerequisites, costs, and limitations.
**Plans**: 1 plan
**Status**: Complete

### Phase 48.1: Machine Info Display (INSERTED)
**Goal**: Web UI displays machine information (type, CPU cores, memory) for debugging and user awareness.
**Plans**: 0 plans
**Status**: Not started

</details>

<details>
<summary>v2.7 Automated GCP Deployment (Phases 49-53) - SHIPPED 2026-02-06</summary>

### Phase 49: Config Foundation and Pricing Query
**Goal**: Non-interactive deployment foundation with reliable pricing data and validated configuration.
**Plans**: 2 plans
**Status**: Complete

### Phase 50: Machine Selection and Resource Calculation
**Goal**: Correct machine type mapping and dynamic Docker resource limit calculation.
**Plans**: 2 plans
**Status**: Complete

### Phase 51: Deployment Orchestration
**Goal**: End-to-end automated deployment with progressive feedback and error handling.
**Plans**: 2 plans
**Status**: Complete

### Phase 52: HTTP-Only Container Deployment
**Goal**: Migrate lifecycle and teardown scripts to v2.7 TOML config with runtime zone detection.
**Plans**: 2 plans
**Status**: Complete

### Phase 53: Conformer Progress Bug Fix
**Goal**: Conformer progress tracking displays correctly during processing.
**Plans**: 1 plan
**Status**: Complete

</details>

## v2.8 Expanded Solvent Support (Phases 54-58)

**Milestone Goal:** Add 4 new NMR solvents (methanol, water, acetone, benzene) with DELTA50-derived B3LYP scaling factors, extending supported solvents from 3 to 7.

### Phase 54: Benchmark Infrastructure
**Goal**: Benchmark tooling accepts all 4 new solvents so calculations can be dispatched
**Depends on**: Nothing (first phase of v2.8)
**Requirements**: BENCH-01, INTG-05
**Success Criteria** (what must be TRUE):
  1. Running `python -m qm_nmr_calc.benchmark run --solvents Methanol --functionals B3LYP --headless` is accepted by the CLI (no "invalid choice" error)
  2. Running the benchmark CLI with Water, Acetone, and Benzene as solvents is also accepted
  3. NWChem COSMO input files generated for each new solvent contain correct solvent parameters (verified by inspecting a generated .nw file)
  4. Benzene is present in input_gen.py SUPPORTED_SOLVENTS (it needs to be added; methanol, water, acetone already exist there)
**Plans**: 1 plan

Plans:
- [x] 54-01-PLAN.md -- Extend benchmark CLI and NWChem input generation for 4 new solvents

### Phase 55: DELTA50 Benchmark Calculations
**Goal**: All 200 NWChem benchmark calculations complete successfully (50 molecules x 4 solvents)
**Depends on**: Phase 54
**Requirements**: BENCH-02, BENCH-03, BENCH-04, BENCH-05
**Success Criteria** (what must be TRUE):
  1. Methanol benchmark data directory contains 50 completed NWChem output files with NMR shielding tensors
  2. Water benchmark data directory contains 50 completed NWChem output files with NMR shielding tensors
  3. Acetone benchmark data directory contains 50 completed NWChem output files with NMR shielding tensors
  4. Benzene benchmark data directory contains 50 completed NWChem output files with NMR shielding tensors
  5. No benchmark molecule fails with a COSMO convergence error (all 200 calculations parse cleanly)
**Plans**: 2 plans

Plans:
- [ ] 55-01-PLAN.md -- Run DELTA50 benchmarks for Methanol and Water solvents (100 calculations)
- [ ] 55-02-PLAN.md -- Run DELTA50 benchmarks for Acetone and Benzene solvents (100 calculations)

### Phase 56: Scaling Factor Derivation
**Goal**: OLS-derived scaling factors for all 4 new solvents pass quality gates and are stored in package data
**Depends on**: Phase 55
**Requirements**: BENCH-06, VALID-01, VALID-02, VALID-03
**Success Criteria** (what must be TRUE):
  1. `python -m qm_nmr_calc.benchmark analyze` produces scaling factors for all 4 new solvents (8 new factor sets: 4 solvents x 2 nuclei)
  2. All 8 new factor sets have R-squared > 0.99
  3. 1H MAE is below 0.2 ppm for each new solvent
  4. 13C MAE is below 3.0 ppm for each new solvent
  5. scaling_factors.json contains 14 total factor sets (6 existing + 8 new) with slope, intercept, R-squared, and MAE for each
**Plans**: TBD

Plans:
- [ ] 56-01-PLAN.md - Run analysis, validate quality gates, update scaling_factors.json

### Phase 57: Solvent Integration
**Goal**: Users can select any of the 4 new solvents in the web UI and API and get accurate NMR predictions
**Depends on**: Phase 56
**Requirements**: INTG-01, INTG-02, INTG-03, INTG-04, INTG-06
**Success Criteria** (what must be TRUE):
  1. Submitting a molecule via API with `solvent=methanol` returns predicted shifts using methanol-specific scaling factors
  2. Submitting a molecule via API with `solvent=water`, `solvent=acetone`, or `solvent=benzene` each returns predicted shifts using the correct solvent-specific scaling factors
  3. Web UI solvent dropdown lists all 7 solvents (CHCl3, DMSO, Vacuum, Methanol, Water, Acetone, Benzene)
  4. Existing CHCl3, DMSO, and vacuum calculations produce identical results to before (no regression)
**Plans**: TBD

Plans:
- [ ] 57-01-PLAN.md - Add 4 solvents to solvents.py, shifts.py, and verify end-to-end

### Phase 58: Documentation
**Goal**: Documentation reflects 7-solvent support with accuracy statistics for each solvent
**Depends on**: Phase 57
**Requirements**: DOCS-01, DOCS-02
**Success Criteria** (what must be TRUE):
  1. SCALING-FACTORS.md lists all 14 factor sets with slope, intercept, R-squared, MAE, and outlier count for each
  2. README states 7 supported solvents and lists them by name
  3. Factor statistics in documentation match values in scaling_factors.json exactly
**Plans**: TBD

Plans:
- [ ] 58-01-PLAN.md - Update SCALING-FACTORS.md and README for 7-solvent support

## Progress

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
| 49. Config & Pricing | v2.7 | 2/2 | Complete | 2026-02-06 |
| 50. Machine Selection | v2.7 | 2/2 | Complete | 2026-02-06 |
| 51. Orchestration | v2.7 | 2/2 | Complete | 2026-02-06 |
| 52. HTTP Container | v2.7 | 2/2 | Complete | 2026-02-06 |
| 53. Conformer Bug Fix | v2.7 | 1/1 | Complete | 2026-02-06 |
| **54. Benchmark Infrastructure** | **v2.8** | **1/1** | **Complete** | **2026-02-07** |
| **55. DELTA50 Calculations** | **v2.8** | **0/2** | **Not started** | **-** |
| **56. Scaling Factor Derivation** | **v2.8** | **0/TBD** | **Not started** | **-** |
| **57. Solvent Integration** | **v2.8** | **0/TBD** | **Not started** | **-** |
| **58. Documentation** | **v2.8** | **0/TBD** | **Not started** | **-** |

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

## Coverage (v2.8)

**v2.8 Requirements: 17 total**

| Requirement | Phase | Description |
|-------------|-------|-------------|
| BENCH-01 | 54 | Benchmark CLI accepts 4 new solvents |
| INTG-05 | 54 | NWChem COSMO correct for all 4 new solvents |
| BENCH-02 | 55 | DELTA50 benchmark for methanol (50 molecules) |
| BENCH-03 | 55 | DELTA50 benchmark for water (50 molecules) |
| BENCH-04 | 55 | DELTA50 benchmark for acetone (50 molecules) |
| BENCH-05 | 55 | DELTA50 benchmark for benzene (50 molecules) |
| BENCH-06 | 56 | OLS scaling factors for 8 new factor sets |
| VALID-01 | 56 | All 8 factor sets R-squared > 0.99 |
| VALID-02 | 56 | 1H MAE below 0.2 ppm per solvent |
| VALID-03 | 56 | 13C MAE below 3.0 ppm per solvent |
| INTG-01 | 57 | Methanol selectable in UI and API |
| INTG-02 | 57 | Water selectable in UI and API |
| INTG-03 | 57 | Acetone selectable in UI and API |
| INTG-04 | 57 | Benzene selectable in UI and API |
| INTG-06 | 57 | Scaling factors loaded correctly for new solvents |
| DOCS-01 | 58 | SCALING-FACTORS.md updated |
| DOCS-02 | 58 | README updated for 7 solvents |

**Mapped: 17/17 (100%)**

---
*Last updated: 2026-02-07 - Phase 54 complete*
