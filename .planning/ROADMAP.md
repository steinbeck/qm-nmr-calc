# Roadmap: qm-nmr-calc

## Milestones

- ✅ **v1.0 Core NMR Service** - Phases 1-6 (shipped 2026-01-20)
- ✅ **v1.1 Accurate Chemical Shifts** - Phases 7-11.2 (shipped 2026-01-25)
- ✅ **v2.0 Conformational Sampling** - Phases 12-17 (shipped 2026-01-28)
- ✅ **v2.0.1 Conformer Pre-selection** - Phase 24 (shipped 2026-01-30)
- ✅ **v2.1 UI Redesign** - Phases 18-23 (shipped 2026-01-31)
- ✅ **v2.2 Documentation** - Phases 25-31 (shipped 2026-02-01)
- ✅ **v2.3 NMReData Export** - Phases 32-34 (shipped 2026-02-01)
- ✅ **v2.4 Docker Deployment** - Phases 35-40 (shipped 2026-02-03)
- ✅ **v2.5 ARM64 Docker Support** - Phases 41-44 (shipped 2026-02-04)
- ✅ **v2.6 GCP Spot Deployment** - Phases 45-49 (shipped 2026-02-05)
- ✅ **v2.7 Automated GCP Deployment** - Phases 50-53 (shipped 2026-02-06)
- ✅ **v2.8 Expanded Solvent Support** - Phases 54-58 (shipped 2026-02-09)
- 🚧 **v2.9 Extended Solvent Coverage** - Phases 59-63 (in progress)

## Phases

<details>
<summary>✅ v1.0 Core NMR Service (Phases 1-6) - SHIPPED 2026-01-20</summary>

### Phase 1: Project Scaffold
**Goal**: Working FastAPI application with basic endpoints
**Plans**: 3 plans

Plans:
- [x] 01-01: Initial setup and FastAPI skeleton
- [x] 01-02: Job submission endpoint
- [x] 01-03: Job status and results endpoints

### Phase 2: Job Queue
**Goal**: Reliable background task processing with Huey
**Plans**: 2 plans

Plans:
- [x] 02-01: Huey integration and task definitions
- [x] 02-02: SQLite persistence and error handling

### Phase 3: NMR Calculation Pipeline
**Goal**: ISiCLE integration for NWChem calculations
**Plans**: 4 plans

Plans:
- [x] 03-01: ISiCLE wrapper and geometry optimization
- [x] 03-02: NMR shielding calculation
- [x] 03-03: Result parsing and storage
- [x] 03-04: Draft/production presets

### Phase 4: Results Visualization
**Goal**: Visual outputs for NMR predictions
**Plans**: 3 plans

Plans:
- [x] 04-01: Spectrum plot generation
- [x] 04-02: Annotated structure drawing
- [x] 04-03: 3D molecule viewer

### Phase 5: Web UI
**Goal**: User-friendly interface for job submission and results
**Plans**: 2 plans

Plans:
- [x] 05-01: Submit page and SMILES input
- [x] 05-02: Status and results pages

### Phase 6: Email Notifications
**Goal**: Notify users when calculations complete
**Plans**: 2 plans

Plans:
- [x] 06-01: SMTP integration and email templates
- [x] 06-02: Opt-in notification system

</details>

<details>
<summary>✅ v1.1 Accurate Chemical Shifts (Phases 7-11.2) - SHIPPED 2026-01-25</summary>

### Phase 7: NWChem Direct Integration
**Goal**: Replace ISiCLE with direct NWChem I/O
**Plans**: 4 plans

Plans:
- [x] 07-01: Input file generation for geometry optimization
- [x] 07-02: Input file generation for NMR shielding
- [x] 07-03: Output parsing for energies and shifts
- [x] 07-04: COSMO solvation support

### Phase 8: Benchmark Infrastructure
**Goal**: DELTA50 benchmark tooling for scaling factor derivation
**Plans**: 3 plans

Plans:
- [x] 08-01: Benchmark CLI and molecule loader
- [x] 08-02: Calculation orchestration
- [x] 08-03: Result aggregation and analysis

### Phase 9: DELTA50 Calculations
**Goal**: Run 150 NWChem benchmark calculations
**Plans**: 3 plans

Plans:
- [x] 09-01: CHCl3 calculations (50 molecules)
- [x] 09-02: DMSO calculations (50 molecules)
- [x] 09-03: Vacuum calculations (50 molecules)

### Phase 10: Scaling Factor Derivation
**Goal**: OLS regression factors for all solvent/nucleus combinations
**Plans**: 2 plans

Plans:
- [x] 10-01: Regression implementation and validation
- [x] 10-02: Factor storage and loading

### Phase 11: Chemical Shift Integration
**Goal**: Apply scaling factors in production pipeline
**Plans**: 4 plans

Plans:
- [x] 11-01: Factor loading and application
- [x] 11-02: API metadata exposure
- [x] 11-03: Vacuum (gas-phase) support
- [x] 11-04: Benchmark viewer page

### Phase 11.1: COSMO Bug Fix (INSERTED)
**Goal**: Fix CPHF convergence failure in NMR+COSMO
**Plans**: 2 plans

Plans:
- [x] 11.1-01: Add `direct` keyword to NWChem NMR input
- [x] 11.1-02: Rerun failed COSMO calculations

### Phase 11.2: 3D Visualization Enhancement (INSERTED)
**Goal**: Interactive molecule viewer on results page
**Plans**: 3 plans

Plans:
- [x] 11.2-01: 3Dmol.js integration
- [x] 11.2-02: Shift label overlays
- [x] 11.2-03: View controls and styling

</details>

<details>
<summary>✅ v2.0 Conformational Sampling (Phases 12-17) - SHIPPED 2026-01-28</summary>

### Phase 12: RDKit Conformer Generation
**Goal**: Distance geometry conformer generation pipeline
**Plans**: 3 plans

Plans:
- [x] 12-01: KDG conformer generation
- [x] 12-02: Energy window filtering
- [x] 12-03: SDF export for NWChem

### Phase 13: CREST Integration
**Goal**: Production-quality conformer searching with CREST/xTB
**Plans**: 3 plans

Plans:
- [x] 13-01: CREST wrapper and detection
- [x] 13-02: ALPB solvation model
- [x] 13-03: RDKit fallback when unavailable

### Phase 14: Multi-Conformer NWChem
**Goal**: DFT calculations on conformer ensembles
**Plans**: 4 plans

Plans:
- [x] 14-01: Parallel geometry optimization
- [x] 14-02: Parallel NMR calculations
- [x] 14-03: Post-DFT energy filtering
- [x] 14-04: Result aggregation

### Phase 15: Boltzmann Weighting
**Goal**: Energy-weighted average chemical shifts
**Plans**: 3 plans

Plans:
- [x] 15-01: Boltzmann population calculation
- [x] 15-02: Weighted shift averaging
- [x] 15-03: Numerical stability (exp-normalize)

### Phase 16: Ensemble API
**Goal**: Expose conformer metadata in responses
**Plans**: 2 plans

Plans:
- [x] 16-01: Conformer count and populations
- [x] 16-02: Energy range and provenance

### Phase 17: Ensemble UI
**Goal**: Progress tracking for ensemble jobs
**Plans**: 2 plans

Plans:
- [x] 17-01: Conformer progress visualization
- [x] 17-02: Ensemble mode selector

</details>

<details>
<summary>✅ v2.0.1 Conformer Pre-selection (Phase 24) - SHIPPED 2026-01-30</summary>

### Phase 24: RMSD Clustering and xTB Ranking
**Goal**: Reduce DFT workload with conformer pre-selection
**Plans**: 3 plans

Plans:
- [x] 24-01: RMSD Butina clustering
- [x] 24-02: xTB energy ranking with ALPB
- [x] 24-03: MMFF fallback when xTB unavailable

</details>

<details>
<summary>✅ v2.1 UI Redesign (Phases 18-23) - SHIPPED 2026-01-31</summary>

### Phase 18: Custom CSS Architecture
**Goal**: Replace Pico CSS with custom design system
**Plans**: 4 plans

Plans:
- [x] 18-01: CSS reset and cascade layers
- [x] 18-02: Design tokens and variables
- [x] 18-03: Responsive grid system
- [x] 18-04: Component library

### Phase 19: Glassmorphism Effects
**Goal**: Backdrop blur and transparency styling
**Plans**: 3 plans

Plans:
- [x] 19-01: Card glassmorphism with backdrop-filter
- [x] 19-02: Safari fallback and accessibility
- [x] 19-03: Hover states and transitions

### Phase 20: Bento Grid Layout
**Goal**: Modern asymmetric card arrangements
**Plans**: 2 plans

Plans:
- [x] 20-01: Grid template system
- [x] 20-02: Breakpoints and mobile optimization

### Phase 21: Results Page Redesign
**Goal**: Hero 3D viewer with organized data cards
**Plans**: 3 plans

Plans:
- [x] 21-01: Full-width 3D viewer
- [x] 21-02: Spectrum and data cards
- [x] 21-03: Download section

### Phase 22: Submit Page Redesign
**Goal**: Clean two-column form with molecule preview
**Plans**: 2 plans

Plans:
- [x] 22-01: Form layout and inputs
- [x] 22-02: SmilesDrawer preview card

### Phase 23: Status Page Redesign
**Goal**: Step tracker with conformer progress
**Plans**: 3 plans

Plans:
- [x] 23-01: Step tracker component
- [x] 23-02: Conformer progress visualization
- [x] 23-03: Status card and metadata

</details>

<details>
<summary>✅ v2.2 Documentation (Phases 25-31) - SHIPPED 2026-02-01</summary>

### Phase 25: README Overhaul
**Goal**: Comprehensive project overview with quick start
**Plans**: 2 plans

Plans:
- [x] 25-01: Architecture diagram and quick start
- [x] 25-02: Feature list and development guide

### Phase 26: Installation Guide
**Goal**: Multi-distro installation instructions
**Plans**: 2 plans

Plans:
- [x] 26-01: NWChem installation per distro
- [x] 26-02: CREST/xTB setup

### Phase 27: Usage Guide
**Goal**: Web UI and API usage examples
**Plans**: 2 plans

Plans:
- [x] 27-01: Web UI workflow
- [x] 27-02: 26 curl examples for REST API

### Phase 28: Architecture Documentation
**Goal**: Technical architecture with diagrams
**Plans**: 2 plans

Plans:
- [x] 28-01: Stack and data flow diagrams
- [x] 28-02: Job lifecycle and error handling

### Phase 29: Library Documentation
**Goal**: Document all external dependencies
**Plans**: 1 plan

Plans:
- [x] 29-01: RDKit, NWChem, Huey, 3Dmol.js, CREST/xTB

### Phase 30: Science Writeup
**Goal**: DP4+ methodology with full derivations
**Plans**: 2 plans

Plans:
- [x] 30-01: DP4+ theory and equations
- [x] 30-02: Literature citations (9 DOIs)

### Phase 31: Cross-Reference Audit
**Goal**: Standardized links across all docs
**Plans**: 1 plan

Plans:
- [x] 31-01: Verify and standardize all doc links

</details>

<details>
<summary>✅ v2.3 NMReData Export (Phases 32-34) - SHIPPED 2026-02-01</summary>

### Phase 32: NMReData SDF Generator
**Goal**: Generate standard-compliant NMReData files
**Plans**: 1 plan

Plans:
- [x] 32-01: NMReData tag generation and SDF export

### Phase 33: API Endpoint
**Goal**: Download endpoint for NMReData files
**Plans**: 1 plan

Plans:
- [x] 33-01: GET /api/v1/jobs/{job_id}/nmredata.sdf

### Phase 34: Web UI Integration
**Goal**: Download button on results page
**Plans**: 1 plan

Plans:
- [x] 34-01: NMReData download button

</details>

<details>
<summary>✅ v2.4 Docker Deployment (Phases 35-40) - SHIPPED 2026-02-03</summary>

### Phase 35: Worker Container
**Goal**: Docker image with NWChem, CREST, xTB
**Plans**: 2 plans

Plans:
- [x] 35-01: Multi-stage Dockerfile with conda
- [x] 35-02: Worker entrypoint and health checks

### Phase 36: API Container
**Goal**: FastAPI application container
**Plans**: 1 plan

Plans:
- [x] 36-01: Multi-stage build with slim Python

### Phase 37: Docker Compose
**Goal**: Orchestrate API, worker, and reverse proxy
**Plans**: 2 plans

Plans:
- [x] 37-01: Compose file with persistent volumes
- [x] 37-02: Caddy reverse proxy with auto-HTTPS

### Phase 38: GHCR Publishing
**Goal**: Pre-built images on GitHub Container Registry
**Plans**: 1 plan

Plans:
- [x] 38-01: GitHub Actions CI/CD workflow

### Phase 39: Deployment Documentation
**Goal**: Comprehensive VPS deployment guide
**Plans**: 1 plan

Plans:
- [x] 39-01: 460-line deployment guide

### Phase 40: Quick Start Integration
**Goal**: Docker as primary getting started path
**Plans**: 1 plan

Plans:
- [x] 40-01: Update README with Docker quick start

</details>

<details>
<summary>✅ v2.5 ARM64 Docker Support (Phases 41-44) - SHIPPED 2026-02-04</summary>

### Phase 41: ARM64 Worker Container
**Goal**: ARM64-compatible worker via conda-forge
**Plans**: 1 plan

Plans:
- [x] 41-01: ARM64 Dockerfile with conda NWChem

### Phase 42: Multi-Arch Images
**Goal**: Publish amd64 + arm64 images to GHCR
**Plans**: 1 plan

Plans:
- [x] 42-01: Update GitHub Actions for multi-arch

### Phase 43: Architecture Detection
**Goal**: Platform-specific docker-compose
**Plans**: 1 plan

Plans:
- [x] 43-01: docker-compose.arm64.yml override

### Phase 44: ARM64 Documentation
**Goal**: Document ARM64 deployment process
**Plans**: 1 plan

Plans:
- [x] 44-01: Update DEPLOYMENT.md with ARM64 section

</details>

<details>
<summary>✅ v2.6 GCP Spot VM Deployment (Phases 45-49) - SHIPPED 2026-02-05</summary>

### Phase 45: Spot VM Provisioning
**Goal**: Create GCP Spot VM with lifecycle scripts
**Plans**: 1 plan

Plans:
- [x] 45-01: gcloud VM creation with startup script

### Phase 46: Persistent Disk
**Goal**: Persistent storage across VM lifecycle
**Plans**: 1 plan

Plans:
- [x] 46-01: Disk attach/detach in startup/shutdown

### Phase 47: Deployment Script
**Goal**: End-to-end GCP deployment automation
**Plans**: 1 plan

Plans:
- [x] 47-01: deploy-gcp.sh orchestration script

### Phase 48: GCP Documentation
**Goal**: Document GCP spot deployment process
**Plans**: 1 plan

Plans:
- [x] 48-01: GCP_DEPLOYMENT.md guide

### Phase 49: Cleanup Automation
**Goal**: Teardown script for resource cleanup
**Plans**: 1 plan

Plans:
- [x] 49-01: destroy-gcp.sh script

</details>

<details>
<summary>✅ v2.7 Automated GCP Deployment (Phases 50-53) - SHIPPED 2026-02-06</summary>

### Phase 50: TOML Configuration
**Goal**: Config-driven deployment with zero prompts
**Plans**: 2 plans

Plans:
- [x] 50-01: TOML schema and Pydantic validation
- [x] 50-02: Config parser integration

### Phase 51: Spot Pricing Discovery
**Goal**: Auto-discover cheapest GCP region
**Plans**: 2 plans

Plans:
- [x] 51-01: CloudPrice.net API integration
- [x] 51-02: Fallback to gcloud pricing

### Phase 52: Dynamic Machine Selection
**Goal**: Auto-select machine type from CPU/RAM requirements
**Plans**: 2 plans

Plans:
- [x] 52-01: Machine type selection logic
- [x] 52-02: Docker resource limits

### Phase 53: Modular Bash Library
**Goal**: Composable deployment library architecture
**Plans**: 3 plans

Plans:
- [x] 53-01: config.sh, pricing.sh, machine.sh
- [x] 53-02: infrastructure.sh orchestration
- [x] 53-03: Dry-run mode

</details>

<details>
<summary>✅ v2.8 Expanded Solvent Support (Phases 54-58) - SHIPPED 2026-02-09</summary>

### Phase 54: Benchmark Infrastructure
**Goal**: Extend CLI for methanol, water, acetone, benzene
**Plans**: 1 plan

Plans:
- [x] 54-01: Update CLI to accept 4 new solvents

### Phase 55: DELTA50 Benchmark Calculations
**Goal**: Run 200 NWChem calculations for 4 new solvents
**Plans**: 2 plans

Plans:
- [x] 55-01: Methanol + Water (100 calculations)
- [x] 55-02: Acetone + Benzene (100 calculations)

### Phase 56: Scaling Factor Derivation
**Goal**: Derive 8 new OLS factor sets with validation
**Plans**: 1 plan

Plans:
- [x] 56-01: Regression and quality gate validation

### Phase 57: Solvent Integration
**Goal**: Wire 4 new solvents through all gatekeeper modules
**Plans**: 1 plan

Plans:
- [x] 57-01: solvents.py, shifts.py, nmredata.py integration

### Phase 58: Documentation
**Goal**: Update factor and solvent documentation
**Plans**: 1 plan

Plans:
- [x] 58-01: SCALING-FACTORS.md and README updates

</details>

### 🚧 v2.9 Extended Solvent Coverage (In Progress)

**Milestone Goal:** Add 6 more NMR solvents (pyridine, THF, toluene, DCM, acetonitrile, DMF) with DELTA50-derived scaling factors, extending from 7 to 13 solvents.

#### Phase 59: Benchmark Infrastructure
**Goal**: Extend CLI and add NWChem COSMO name mapping for 6 new solvents
**Depends on**: Phase 58
**Requirements**: BENCH-01
**Success Criteria** (what must be TRUE):
  1. Benchmark CLI accepts pyridine, thf, toluene, dcm, acetonitrile, and dmf as valid solvents
  2. Input generation correctly maps user-friendly names to NWChem COSMO names (acetonitrile → acetntrl)
  3. Test calculations run successfully with all 6 new solvents
**Plans**: 1 plan

Plans:
- [ ] 59-01: Extend CLI validation and add COSMO name mapping layer

#### Phase 60: DELTA50 Pyridine + THF
**Goal**: Complete 100 benchmark calculations for pyridine and THF
**Depends on**: Phase 59
**Requirements**: BENCH-02, BENCH-03
**Success Criteria** (what must be TRUE):
  1. All 50 pyridine DELTA50 molecules calculate successfully with COSMO solvation
  2. All 50 THF DELTA50 molecules calculate successfully with COSMO solvation
  3. Calculated 1H and 13C shifts extracted and stored for both solvents
**Plans**: 1 plan

Plans:
- [ ] 60-01: Run pyridine (50 molecules) + THF (50 molecules) benchmark calculations

#### Phase 61: DELTA50 Toluene + DCM
**Goal**: Complete 100 benchmark calculations for toluene and DCM
**Depends on**: Phase 60
**Requirements**: BENCH-04, BENCH-05
**Success Criteria** (what must be TRUE):
  1. All 50 toluene DELTA50 molecules calculate successfully with COSMO solvation
  2. All 50 DCM DELTA50 molecules calculate successfully with COSMO solvation
  3. Calculated 1H and 13C shifts extracted and stored for both solvents
**Plans**: 1 plan

Plans:
- [ ] 61-01: Run toluene (50 molecules) + DCM (50 molecules) benchmark calculations

#### Phase 62: DELTA50 Acetonitrile + DMF
**Goal**: Complete 100 benchmark calculations for acetonitrile and DMF
**Depends on**: Phase 61
**Requirements**: BENCH-06, BENCH-07
**Success Criteria** (what must be TRUE):
  1. All 50 acetonitrile DELTA50 molecules calculate successfully with COSMO solvation (using acetntrl COSMO name)
  2. All 50 DMF DELTA50 molecules calculate successfully with COSMO solvation
  3. Calculated 1H and 13C shifts extracted and stored for both solvents
**Plans**: 1 plan

Plans:
- [ ] 62-01: Run acetonitrile (50 molecules) + DMF (50 molecules) benchmark calculations

#### Phase 63: Scaling Factor Derivation
**Goal**: Derive and validate 12 new OLS scaling factor sets for all 6 solvents
**Depends on**: Phase 62
**Requirements**: BENCH-08, VALID-01, VALID-02, VALID-03
**Success Criteria** (what must be TRUE):
  1. OLS regression produces 12 scaling factor sets (6 solvents × 2 nuclei) with slope and intercept
  2. All 12 factor sets pass R² > 0.99 quality gate
  3. All 6 solvents pass 1H MAE < 0.2 ppm quality gate
  4. All 6 solvents pass 13C MAE < 3.0 ppm quality gate
  5. Scaling factors serialized to JSON for production use
**Plans**: 1 plan

Plans:
- [ ] 63-01: Derive OLS factors and validate quality gates

#### Phase 64: Solvent Integration
**Goal**: Wire all 6 new solvents through production pipeline (API, UI, shifts, NMReData)
**Depends on**: Phase 63
**Requirements**: INTG-01, INTG-02, INTG-03, INTG-04, INTG-05, INTG-06, INTG-07, INTG-08
**Success Criteria** (what must be TRUE):
  1. User can select pyridine, THF, toluene, DCM, acetonitrile, and DMF in web UI solvent dropdown
  2. API accepts all 6 new solvents in POST /api/v1/jobs/submit
  3. Input generation applies correct NWChem COSMO names for all 6 solvents
  4. Scaling factors load and apply correctly for all 6 new solvents
  5. NMReData export uses correct deuterated solvent names for all 6 solvents
**Plans**: 1 plan

Plans:
- [ ] 64-01: Update solvents.py, shifts.py, nmredata.py, and UI/API for 6 new solvents

#### Phase 65: Documentation
**Goal**: Update documentation to reflect 13-solvent support
**Depends on**: Phase 64
**Requirements**: DOCS-01, DOCS-02
**Success Criteria** (what must be TRUE):
  1. SCALING-FACTORS.md documents all 26 factor sets (13 solvents × 2 nuclei) with R² and MAE
  2. README lists all 13 supported solvents
  3. Documentation cross-references remain valid
**Plans**: 1 plan

Plans:
- [ ] 65-01: Update SCALING-FACTORS.md and README

## Progress

**Execution Order:**
Phases execute in numeric order: 59 → 60 → 61 → 62 → 63 → 64 → 65

| Phase | Milestone | Plans Complete | Status | Completed |
|-------|-----------|----------------|--------|-----------|
| 1. Project Scaffold | v1.0 | 3/3 | Complete | 2026-01-19 |
| 2. Job Queue | v1.0 | 2/2 | Complete | 2026-01-19 |
| 3. NMR Calculation Pipeline | v1.0 | 4/4 | Complete | 2026-01-19 |
| 4. Results Visualization | v1.0 | 3/3 | Complete | 2026-01-20 |
| 5. Web UI | v1.0 | 2/2 | Complete | 2026-01-20 |
| 6. Email Notifications | v1.0 | 2/2 | Complete | 2026-01-20 |
| 7. NWChem Direct Integration | v1.1 | 4/4 | Complete | 2026-01-21 |
| 8. Benchmark Infrastructure | v1.1 | 3/3 | Complete | 2026-01-22 |
| 9. DELTA50 Calculations | v1.1 | 3/3 | Complete | 2026-01-23 |
| 10. Scaling Factor Derivation | v1.1 | 2/2 | Complete | 2026-01-24 |
| 11. Chemical Shift Integration | v1.1 | 4/4 | Complete | 2026-01-25 |
| 11.1. COSMO Bug Fix | v1.1 | 2/2 | Complete | 2026-01-25 |
| 11.2. 3D Visualization Enhancement | v1.1 | 3/3 | Complete | 2026-01-25 |
| 12. RDKit Conformer Generation | v2.0 | 3/3 | Complete | 2026-01-26 |
| 13. CREST Integration | v2.0 | 3/3 | Complete | 2026-01-26 |
| 14. Multi-Conformer NWChem | v2.0 | 4/4 | Complete | 2026-01-27 |
| 15. Boltzmann Weighting | v2.0 | 3/3 | Complete | 2026-01-27 |
| 16. Ensemble API | v2.0 | 2/2 | Complete | 2026-01-28 |
| 17. Ensemble UI | v2.0 | 2/2 | Complete | 2026-01-28 |
| 18. Custom CSS Architecture | v2.1 | 4/4 | Complete | 2026-01-29 |
| 19. Glassmorphism Effects | v2.1 | 3/3 | Complete | 2026-01-30 |
| 20. Bento Grid Layout | v2.1 | 2/2 | Complete | 2026-01-30 |
| 21. Results Page Redesign | v2.1 | 3/3 | Complete | 2026-01-31 |
| 22. Submit Page Redesign | v2.1 | 2/2 | Complete | 2026-01-31 |
| 23. Status Page Redesign | v2.1 | 3/3 | Complete | 2026-01-31 |
| 24. RMSD Clustering and xTB Ranking | v2.0.1 | 3/3 | Complete | 2026-01-30 |
| 25. README Overhaul | v2.2 | 2/2 | Complete | 2026-01-31 |
| 26. Installation Guide | v2.2 | 2/2 | Complete | 2026-01-31 |
| 27. Usage Guide | v2.2 | 2/2 | Complete | 2026-02-01 |
| 28. Architecture Documentation | v2.2 | 2/2 | Complete | 2026-02-01 |
| 29. Library Documentation | v2.2 | 1/1 | Complete | 2026-02-01 |
| 30. Science Writeup | v2.2 | 2/2 | Complete | 2026-02-01 |
| 31. Cross-Reference Audit | v2.2 | 1/1 | Complete | 2026-02-01 |
| 32. NMReData SDF Generator | v2.3 | 1/1 | Complete | 2026-02-01 |
| 33. API Endpoint | v2.3 | 1/1 | Complete | 2026-02-01 |
| 34. Web UI Integration | v2.3 | 1/1 | Complete | 2026-02-01 |
| 35. Worker Container | v2.4 | 2/2 | Complete | 2026-02-02 |
| 36. API Container | v2.4 | 1/1 | Complete | 2026-02-02 |
| 37. Docker Compose | v2.4 | 2/2 | Complete | 2026-02-02 |
| 38. GHCR Publishing | v2.4 | 1/1 | Complete | 2026-02-03 |
| 39. Deployment Documentation | v2.4 | 1/1 | Complete | 2026-02-03 |
| 40. Quick Start Integration | v2.4 | 1/1 | Complete | 2026-02-03 |
| 41. ARM64 Worker Container | v2.5 | 1/1 | Complete | 2026-02-04 |
| 42. Multi-Arch Images | v2.5 | 1/1 | Complete | 2026-02-04 |
| 43. Architecture Detection | v2.5 | 1/1 | Complete | 2026-02-04 |
| 44. ARM64 Documentation | v2.5 | 1/1 | Complete | 2026-02-04 |
| 45. Spot VM Provisioning | v2.6 | 1/1 | Complete | 2026-02-05 |
| 46. Persistent Disk | v2.6 | 1/1 | Complete | 2026-02-05 |
| 47. Deployment Script | v2.6 | 1/1 | Complete | 2026-02-05 |
| 48. GCP Documentation | v2.6 | 1/1 | Complete | 2026-02-05 |
| 49. Cleanup Automation | v2.6 | 1/1 | Complete | 2026-02-05 |
| 50. TOML Configuration | v2.7 | 2/2 | Complete | 2026-02-06 |
| 51. Spot Pricing Discovery | v2.7 | 2/2 | Complete | 2026-02-06 |
| 52. Dynamic Machine Selection | v2.7 | 2/2 | Complete | 2026-02-06 |
| 53. Modular Bash Library | v2.7 | 3/3 | Complete | 2026-02-06 |
| 54. Benchmark Infrastructure | v2.8 | 1/1 | Complete | 2026-02-07 |
| 55. DELTA50 Benchmark Calculations | v2.8 | 2/2 | Complete | 2026-02-08 |
| 56. Scaling Factor Derivation | v2.8 | 1/1 | Complete | 2026-02-08 |
| 57. Solvent Integration | v2.8 | 1/1 | Complete | 2026-02-09 |
| 58. Documentation | v2.8 | 1/1 | Complete | 2026-02-09 |
| 59. Benchmark Infrastructure | v2.9 | 0/1 | Not started | - |
| 60. DELTA50 Pyridine + THF | v2.9 | 0/1 | Not started | - |
| 61. DELTA50 Toluene + DCM | v2.9 | 0/1 | Not started | - |
| 62. DELTA50 Acetonitrile + DMF | v2.9 | 0/1 | Not started | - |
| 63. Scaling Factor Derivation | v2.9 | 0/1 | Not started | - |
| 64. Solvent Integration | v2.9 | 0/1 | Not started | - |
| 65. Documentation | v2.9 | 0/1 | Not started | - |

## Coverage

**v2.9 Requirements:**
- Total: 22 requirements
- Mapped: 22 requirements (100%)
- Unmapped: 0 requirements

**Requirement-to-Phase Mapping:**
- BENCH-01 → Phase 59
- BENCH-02 → Phase 60
- BENCH-03 → Phase 60
- BENCH-04 → Phase 61
- BENCH-05 → Phase 61
- BENCH-06 → Phase 62
- BENCH-07 → Phase 62
- BENCH-08 → Phase 63
- VALID-01 → Phase 63
- VALID-02 → Phase 63
- VALID-03 → Phase 63
- INTG-01 → Phase 64
- INTG-02 → Phase 64
- INTG-03 → Phase 64
- INTG-04 → Phase 64
- INTG-05 → Phase 64
- INTG-06 → Phase 64
- INTG-07 → Phase 64
- INTG-08 → Phase 64
- DOCS-01 → Phase 65
- DOCS-02 → Phase 65
