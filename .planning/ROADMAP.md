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
- ✅ **v2.9 Extended Solvent Coverage** - Phases 59-65 (shipped 2026-02-11)
- 🚧 **v3.0 Publication & Dataset Release** - Phases 66-70 (in progress)

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

<details>
<summary>✅ v2.9 Extended Solvent Coverage (Phases 59-65) - SHIPPED 2026-02-11</summary>

### Phase 59: Benchmark Infrastructure
**Goal**: Extend CLI and add NWChem COSMO name mapping for 6 new solvents
**Plans**: 1 plan

Plans:
- [x] 59-01: Add 6 solvents to input_gen + CLI with COSMO name mapping for acetonitrile

### Phase 60: DELTA50 Pyridine + THF
**Goal**: Complete 100 benchmark calculations for pyridine and THF
**Plans**: 1 plan

Plans:
- [x] 60-01: Run Pyridine (50 molecules) + THF (50 molecules) benchmark calculations

### Phase 61: DELTA50 Toluene + DCM
**Goal**: Complete 100 benchmark calculations for toluene and DCM
**Plans**: 1 plan

Plans:
- [x] 61-01: Verify pre-existing Toluene + DCM results and generate BENCHMARK-RESULTS-TD.md

### Phase 62: DELTA50 Acetonitrile + DMF
**Goal**: Complete 100 benchmark calculations for acetonitrile and DMF
**Plans**: 1 plan

Plans:
- [x] 62-01: Verify pre-existing Acetonitrile + DMF results and generate BENCHMARK-RESULTS-AD.md

### Phase 63: Scaling Factor Derivation
**Goal**: Derive and validate 12 new OLS scaling factor sets for all 6 solvents
**Plans**: 1 plan

Plans:
- [x] 63-01: Derive OLS factors and validate quality gates

### Phase 64: Solvent Integration
**Goal**: Wire all 6 new solvents through production pipeline (API, UI, shifts, NMReData)
**Plans**: 1 plan

Plans:
- [x] 64-01: Update solvents.py, shifts.py, nmredata.py, and UI/API for 6 new solvents

### Phase 65: Documentation
**Goal**: Update documentation to reflect 13-solvent support
**Plans**: 1 plan

Plans:
- [x] 65-01: Add vacuum to SCALING-FACTORS.md (26 total), update README to 13 solvents, fix stale cross-references in docs/

</details>

---

### Phase 66: Dataset Archive Preparation
**Goal**: DELTA50 dataset ready for repository upload with FAIR-compliant metadata
**Dependencies**: None (starts milestone)
**Requirements**: DATA-01, DATA-02, DATA-03, DATA-04, DATA-05, DATA-06, DATA-07, DATA-08, DATA-09

**Success Criteria**:
1. User can navigate publications/dataset/ directory and find all 650 NWChem calculations organized by molecule/solvent hierarchy with .nw inputs and .out outputs
2. User can open processed data exports (CSV for scaling factors, JSON for shielding tensors) and verify molecule metadata completeness
3. User can read README.md and understand dataset scope, citation format, usage instructions without external documentation
4. User can examine metadata.json and verify all DataCite 4.6 mandatory fields plus chemistry extensions (SMILES/InChI/InChIKey for 50 molecules) are complete
5. User can verify computational provenance checklist documents all non-default parameters (NWChem version with patch, COSMO grid settings, compiler info, all settings)

### Phase 67: Repository Upload + DOI
**Goal**: Dataset uploaded to Radar4Chem with reserved DOI before manuscript writing
**Dependencies**: Phase 66 (requires complete archive)
**Requirements**: REPO-01, REPO-02, REPO-03, REPO-04

**Success Criteria**:
1. User can access dataset on Radar4Chem (or Zenodo/Figshare fallback) and see immediate publication status (no embargo, open science)
2. User can copy reserved DOI from repository page and verify it resolves to live dataset
3. User can examine repository metadata form and confirm all DataCite mandatory fields plus chemistry-specific extensions are populated
4. User can cite dataset DOI in manuscript drafts from day one

### Phase 68: Data Descriptor Manuscript
**Goal**: Nature Scientific Data manuscript submission-ready describing DELTA50 dataset
**Dependencies**: Phase 67 (needs dataset DOI for citation)
**Requirements**: DESC-01, DESC-02, DESC-03, DESC-04, DESC-05, DESC-06, DESC-07, DESC-08, DESC-09, DESC-10

**Success Criteria**:
1. User can read Background & Summary section and immediately understand dataset motivation, scope, and reuse value in 2-3 paragraphs
2. User can follow Methods section step-by-step to reproduce calculations (DFT functional, basis set, COSMO params with grid settings, SCF criteria documented)
3. User can examine Technical Validation section and find quantitative metrics for all 13 solvents (RMSE, R², MAE per solvent/nucleus, outlier analysis, literature comparison)
4. User can review 5-6 publication-quality figures and understand molecule diversity, workflow, cross-solvent accuracy patterns without reading main text
5. User can open submission/ directory and find flattened structure ready for journal upload portal

### Phase 69: Application Note Manuscript
**Goal**: Journal of Cheminformatics manuscript submission-ready describing qm-nmr-calc tool
**Dependencies**: Phase 67 (needs dataset DOI for citation)
**Requirements**: APP-01, APP-02, APP-03, APP-04, APP-05, APP-06, APP-07, APP-08, APP-09

**Success Criteria**:
1. User can read Background section and understand problem context (why DELTA50 scaling factors are needed for NMR structure elucidation)
2. User can follow Operation section and successfully submit test job via web UI or API using provided examples
3. User can examine Use Cases section and find 3-5 practical examples showing real-world applications (structure elucidation, diastereomer assignment, solvent selection)
4. User can review benchmarking section and see qm-nmr-calc compared against 4-6 state-of-the-art methods with quantitative metrics (accuracy, computational cost)
5. User can open submission/ directory and find flattened structure ready for journal upload

### Phase 70: AI Writing Validation
**Goal**: Both manuscripts pass audit for AI writing patterns and demonstrate authorial voice
**Dependencies**: Phase 68, Phase 69 (requires completed manuscripts)
**Requirements**: VALID-01, VALID-02, VALID-03, VALID-04, VALID-05

**Success Criteria**:
1. User can search both manuscripts for formulaic transitions (Moreover/Additionally/Furthermore) and find varied, natural language instead
2. User can verify all citations resolve to real papers with correct DOIs (no phantom references, all metadata accurate)
3. User can read both manuscripts and hear distinct authorial voice with domain-specific insights and critical analysis (not generic descriptions)
4. User can examine paragraph openings across each paper and find natural variation in sentence structure (no rigid patterns)
5. User can compare technical terminology usage across both papers and find consistent, correct application throughout

---

## Progress

| Milestone | Phases | Plans | Status |
|-----------|--------|-------|--------|
| v1.0 Core NMR Service | 6 | 16 | Shipped 2026-01-20 |
| v1.1 Accurate Chemical Shifts | 8 | 21 | Shipped 2026-01-25 |
| v2.0 Conformational Sampling | 6 | 18 | Shipped 2026-01-28 |
| v2.0.1 Conformer Pre-selection | 1 | 3 | Shipped 2026-01-30 |
| v2.1 UI Redesign | 6 | 17 | Shipped 2026-01-31 |
| v2.2 Documentation | 7 | 10 | Shipped 2026-02-01 |
| v2.3 NMReData Export | 3 | 3 | Shipped 2026-02-01 |
| v2.4 Docker Deployment | 6 | 8 | Shipped 2026-02-03 |
| v2.5 ARM64 Docker Support | 4 | 4 | Shipped 2026-02-04 |
| v2.6 GCP Spot Deployment | 5 | 4 | Shipped 2026-02-05 |
| v2.7 Automated GCP Deployment | 5 | 9 | Shipped 2026-02-06 |
| v2.8 Expanded Solvent Support | 5 | 6 | Shipped 2026-02-09 |
| v2.9 Extended Solvent Coverage | 7 | 7 | Shipped 2026-02-11 |
| **v3.0 Publication & Dataset Release** | **5** | **0** | **In Progress** |

---
*Last updated: 2026-02-11*
