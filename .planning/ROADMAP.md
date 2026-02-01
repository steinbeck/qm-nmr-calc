# Roadmap: qm-nmr-calc

## Milestones

- ✅ **v1.0 MVP** - Phases 1-6 (shipped 2026-01-20)
- ✅ **v1.1 Accurate Chemical Shifts** - Phases 7-11.2 (shipped 2026-01-25)
- ✅ **v2.0 Conformational Sampling** - Phases 12-17 (shipped 2026-01-28)
- ✅ **v2.0.1 Conformer Pre-selection** - Phase 24 (shipped 2026-01-30)
- ✅ **v2.1 UI Redesign** - Phases 18-23 (shipped 2026-01-31)
- 🔄 **v2.2 Documentation** - Phases 25-31

## Overview

**Current milestone:** v2.2 Documentation - Comprehensive documentation including README overhaul, installation/usage guides, technical architecture, library documentation, and thorough DP4+/NMR science writeup with full derivation.

**Audience:** Both academic researchers and developers/contributors.

See `.planning/REQUIREMENTS.md` for detailed requirements.

## Phases

<details>
<summary>✅ v1.0 MVP (Phases 1-6) - SHIPPED 2026-01-20</summary>

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
<summary>✅ v1.1 Accurate Chemical Shifts (Phases 7-11.2) - SHIPPED 2026-01-25</summary>

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
<summary>✅ v2.0 Conformational Sampling (Phases 12-17) - SHIPPED 2026-01-28</summary>

### Phase 12: Conformer Data Model and Storage
**Goal**: Foundation for multi-conformer calculations -- data structures, file organization, backward compatibility
**Depends on**: Phase 11.2
**Requirements**: DFT-03, API-04 (partial)
**Success Criteria** (what must be TRUE):
  1. System tracks conformer IDs, energies (with units), and geometries in structured format
  2. Each conformer uses isolated scratch directory to prevent NWChem file conflicts
  3. Canonical atom ordering established and verified across conformer lifecycle
  4. Existing v1.x single-conformer jobs load and run without modification
  5. Job directory structure supports per-conformer outputs (`output/conformers/`, `output/optimized/`)
**Plans**: 3 plans
**Status**: Complete

### Phase 13: RDKit KDG Conformer Generation
**Goal**: Generate conformer ensembles using RDKit distance geometry (no external dependencies)
**Depends on**: Phase 12
**Requirements**: CONF-03, CONF-06, FILT-01
**Success Criteria** (what must be TRUE):
  1. System generates conformers using KDG method (pure distance geometry, no crystal bias)
  2. MMFF optimization and RMSD deduplication produces diverse low-energy set
  3. Pre-DFT energy window filter (default 6 kcal/mol) removes high-energy conformers
  4. Adaptive conformer count scales with rotatable bond count
  5. App works fully without CREST/xTB installed (RDKit-only mode validated)
**Plans**: 3 plans
**Status**: Complete

### Phase 14: Boltzmann Averaging Implementation
**Goal**: Numerically stable Boltzmann weighting and ensemble averaging
**Depends on**: Phase 13
**Requirements**: BOLTZ-01, BOLTZ-02, BOLTZ-03, BOLTZ-04
**Success Criteria** (what must be TRUE):
  1. System calculates Boltzmann weights from DFT energies using exp-normalize trick (no overflow/underflow)
  2. System computes population-weighted average chemical shifts across conformer ensemble
  3. Temperature parameter (default 298.15 K) correctly applied in Boltzmann calculation
  4. Unit tests validate numerics with wide energy ranges (0-20 kcal/mol) and known test cases
  5. Single-conformer edge case (ensemble of 1) returns correct degeneracy (weight = 1.0)
**Plans**: 2 plans
**Status**: Complete

### Phase 15: Multi-Conformer NWChem Integration
**Goal**: Sequential DFT optimization and NMR calculation loop for conformer ensemble
**Depends on**: Phase 14
**Requirements**: DFT-01, DFT-02, DFT-04, DFT-05, FILT-02, FILT-03
**Success Criteria** (what must be TRUE):
  1. System runs DFT geometry optimization on each conformer surviving pre-DFT filter
  2. System extracts DFT energies from optimization step for Boltzmann weighting
  3. Post-DFT energy window filter (default 3 kcal/mol) applied before NMR calculations
  4. System runs NMR shielding calculation on each conformer surviving post-DFT filter
  5. Partial conformer failures handled gracefully (job continues with successful conformers if >50% succeed)
  6. User can configure energy window parameters (pre-DFT and post-DFT thresholds)
**Plans**: 3 plans
**Status**: Complete

### Phase 16: CREST Integration (Optional High-Accuracy Mode)
**Goal**: Production-quality conformer generation using CREST metadynamics (when available)
**Depends on**: Phase 15
**Requirements**: CONF-04, CONF-05
**Success Criteria** (what must be TRUE):
  1. System generates conformers using CREST/xTB when binaries available on PATH
  2. System auto-detects CREST/xTB availability and reports status in API health endpoint
  3. CREST includes ALPB solvation model matching job solvent parameter
  4. CREST timeout (default 7200s) with clear error message prevents hanging on macrocycles (fail-fast, no auto-fallback)
  5. Environment variables (OMP_STACKSIZE, GFORTRAN_UNBUFFERED_ALL) configured for stability
**Plans**: 3 plans
**Status**: Complete

### Phase 17: API Integration and Progress Tracking
**Goal**: User-facing ensemble mode with conformer selection, metadata, and progress updates
**Depends on**: Phase 16
**Requirements**: CONF-01, CONF-02, API-01, API-02, API-03, PROG-01, PROG-02
**Success Criteria** (what must be TRUE):
  1. User can choose between single-conformer mode (v1.x behavior) and ensemble mode when submitting job
  2. User can choose RDKit KDG or CREST as conformer generation method in ensemble mode
  3. API returns Boltzmann-weighted average shifts (not per-conformer detail)
  4. API response includes ensemble metadata (conformer count, energy range, top 3 populations, method used, temperature)
  5. 3D viewer shows lowest-energy conformer geometry with shift labels
  6. Job status includes ensemble-specific states (generating_conformers, optimizing_conformers X/N, calculating_nmr X/N, averaging_shifts)
  7. Web UI displays conformer progress for ensemble jobs
**Plans**: 5 plans
**Status**: Complete

</details>

<details>
<summary>✅ v2.1 UI Redesign (Phases 18-23) - SHIPPED 2026-01-31</summary>

### Phase 18: CSS Foundation and Design System
**Goal**: Reusable CSS architecture with design tokens, glass components, and bento grid utilities
**Plans**: 4 plans
**Status**: Complete

### Phase 19: Results Page Redesign
**Goal**: Bento grid layout for results page with prominent visualizations
**Plans**: 3 plans
**Status**: Complete

### Phase 20: Submit Page Redesign
**Goal**: Clean form layout with logical grouping and solid backgrounds
**Plans**: 2 plans
**Status**: Complete

### Phase 21: Status Page Redesign
**Goal**: Job progress visualization with step tracking
**Plans**: 2 plans
**Status**: Complete

### Phase 22: Responsive and Layout Polish
**Goal**: Mobile-first breakpoints with performance-optimized glass effects
**Plans**: 2 plans
**Status**: Complete

### Phase 23: Accessibility and Testing
**Goal**: WCAG compliance, cross-browser validation, and performance optimization
**Plans**: 2 plans
**Status**: Complete

*Full details: .planning/milestones/v2.1-ROADMAP.md*

</details>

<details>
<summary>✅ v2.0.1 Conformer Pre-selection Hotfix - SHIPPED 2026-01-30</summary>

**Hotfix Goal:** Reduce conformers sent to expensive DFT optimization from ~40 to ~8 using RMSD clustering and xTB energy ranking.

### Phase 24: Conformer Pre-selection with xTB
**Goal**: Efficient conformer filtering using RMSD clustering and xTB semi-empirical ranking
**Depends on**: Phase 17 (v2.0 complete)
**Requirements**: Performance optimization for ensemble calculations
**Success Criteria** (what must be TRUE):
  1. RMSD-based Butina clustering reduces redundant conformers (40 → ~8-12 clusters)
  2. xTB (GFN2) single-point energy ranking available when xTB binary detected
  3. System selects ~8 diverse, low-energy conformers for DFT optimization
  4. Fallback to MMFF ranking when xTB not available (RDKit-only mode works)
  5. Flexible molecules (hexanol, decane) complete in <3 hours instead of 10+ hours
**Plans**: 3 plans

Plans:
- [x] 24-01-PLAN.md — RMSD clustering with Butina algorithm (RDKit)
- [x] 24-02-PLAN.md — xTB integration for energy ranking
- [x] 24-03-PLAN.md — Pipeline integration and testing
**Status**: Complete

</details>

## Progress

**Execution Order:**
Phases execute in numeric order: 1 -> 6 (v1.0) -> 7 -> 11.2 (v1.1) -> 12 -> 17 (v2.0) -> 18 -> 23 (v2.1)

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

## v2.2 Documentation (Phases 25-31)

### Phase 25: README and Documentation Structure
**Goal**: Overhaul README.md with clear introduction, architecture overview, and links to detailed docs
**Requirements**: README-01, README-02, README-03, README-04, STRUCT-01, STRUCT-02
**Success Criteria** (what must be TRUE):
  1. README.md provides compelling introduction explaining what the app does and its value
  2. High-level architecture diagram shows system components and data flow
  3. docs/ directory created with logical organization
  4. README links to all detailed documentation pages
  5. Quick start section enables users to get running immediately
**Plans**: 1 plan

Plans:
- [x] 25-01-PLAN.md — README overhaul and docs/ structure
**Status**: Complete

### Phase 26: Installation Guide
**Goal**: Comprehensive installation documentation covering all dependencies and configurations
**Requirements**: INSTALL-01, INSTALL-02, INSTALL-03, INSTALL-04, INSTALL-05
**Success Criteria** (what must be TRUE):
  1. System dependencies (NWChem, MPI, Python) fully documented with version requirements
  2. uv package manager installation explained
  3. Optional CREST/xTB setup documented for ensemble mode
  4. Environment validation steps included
  5. Common installation issues and troubleshooting covered
**Plans**: 1 plan

Plans:
- [x] 26-01-PLAN.md — Complete installation documentation
**Status**: Complete

### Phase 27: Usage Guide
**Goal**: User-facing documentation for web UI and REST API workflows
**Requirements**: USAGE-01, USAGE-02, USAGE-03, USAGE-04, USAGE-05, USAGE-06
**Success Criteria** (what must be TRUE):
  1. Web UI workflow documented with descriptions of each page
  2. REST API fully documented with curl examples for all endpoints
  3. Single-conformer vs ensemble mode selection explained
  4. Solvent selection and implications documented
  5. Preset differences (draft vs production) explained
  6. Result interpretation guide (shifts, spectra, 3D viewer)
**Plans**: 2 plans

Plans:
- [x] 27-01-PLAN.md — Web UI workflow and calculation concepts
- [x] 27-02-PLAN.md — REST API reference and result interpretation
**Status**: Complete

### Phase 28: Technical Architecture
**Goal**: Developer-facing documentation of system architecture and data flows
**Requirements**: TECH-01, TECH-02, TECH-03, TECH-04, TECH-05, TECH-06
**Success Criteria** (what must be TRUE):
  1. Full stack explained (FastAPI, Huey, NWChem, RDKit, 3Dmol.js, CSS)
  2. Data flow diagram from job submission to results
  3. Job lifecycle states and transitions documented
  4. File storage structure and job directory layout explained
  5. Conformer ensemble pipeline stages documented
  6. CSS architecture (layers, components, tokens) explained
**Plans**: 2 plans

Plans:
- [x] 28-01-PLAN.md — Core architecture (stack, data flow, job lifecycle, file storage)
- [x] 28-02-PLAN.md — Conformer pipeline and CSS architecture
**Status**: Complete

### Phase 29: Library Documentation
**Goal**: Documentation of third-party library integrations and usage patterns
**Requirements**: LIB-01, LIB-02, LIB-03, LIB-04, LIB-05, LIB-06
**Success Criteria** (what must be TRUE):
  1. RDKit usage documented (SMILES, conformers, visualization)
  2. NWChem integration documented (input generation, output parsing)
  3. Huey task queue usage documented
  4. 3Dmol.js integration documented (viewer, shift labels)
  5. CREST/xTB integration documented (optional ensemble mode)
  6. SmilesDrawer usage documented
**Plans**: 2 plans

Plans:
- [ ] 29-01-PLAN.md — Python libraries (RDKit, NWChem, Huey)
- [ ] 29-02-PLAN.md — JavaScript libraries and CREST/xTB (3Dmol.js, SmilesDrawer, CREST)
**Status**: Not Started

### Phase 30: DP4+ Science Documentation
**Goal**: Thorough scientific writeup of NMR prediction methodology with full derivations
**Requirements**: DP4-01, DP4-02, DP4-03, DP4-04, DP4-05, DP4-06, DP4-07, DP4-08, DP4-09
**Success Criteria** (what must be TRUE):
  1. NMR chemical shift prediction fundamentals explained
  2. DFT theory basis (B3LYP, basis sets, GIAO) documented
  3. COSMO solvation model explained
  4. Linear scaling methodology fully derived
  5. DELTA50 benchmark and scaling factor derivation documented
  6. Boltzmann weighting mathematically derived
  7. Conformational sampling importance explained with references
  8. Literature citations included (ISiCLE, DELTA50, CREST, Grimblat)
  9. Expected accuracy (MAE values) and limitations documented
**Plans**: TBD
**Status**: Not Started

### Phase 31: Documentation Polish and Cross-References
**Goal**: Final review ensuring consistency and proper cross-linking
**Requirements**: STRUCT-03, STRUCT-04
**Success Criteria** (what must be TRUE):
  1. All docs cross-reference related pages appropriately
  2. Consistent formatting and heading structure across all docs
  3. All links verified working
  4. Table of contents accurate
**Plans**: TBD
**Status**: Not Started

---

*Last updated: 2026-02-01 - Phase 29 planned*
