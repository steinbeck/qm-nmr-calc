# Roadmap: qm-nmr-calc

## Milestones

- ✅ **v1.0 MVP** - Phases 1-6 (shipped 2026-01-20)
- ✅ **v1.1 Accurate Chemical Shifts** - Phases 7-11.2 (shipped 2026-01-25)
- 🚧 **v2.0 Conformational Sampling** - Phases 12-17 (in progress)

## Overview

v2.0 adds Boltzmann-weighted ensemble NMR averaging to enable accurate predictions for flexible molecules. Flexible molecules adopt multiple low-energy conformations in solution; experimental NMR spectra are population-weighted averages across the conformer ensemble. This milestone addresses the single-conformer limitation identified in v1.x.

The roadmap follows a risk-mitigation strategy: establish data foundations first (Phase 12), validate end-to-end workflow with RDKit-only path (Phases 13-15, 17), then add CREST complexity as optional enhancement (Phase 16). Critical architectural decisions implemented early: per-conformer scratch directories prevent NWChem conflicts, canonical atom ordering ensures consistency, numerically stable Boltzmann implementation handles wide energy ranges.

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

### v2.0 Conformational Sampling (In Progress)

**Milestone Goal:** Boltzmann-weighted ensemble averaging for flexible molecules -- generate conformers, compute NMR on each, weight by DFT energies for population-averaged shifts.

#### Phase 12: Conformer Data Model and Storage
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

Plans:
- [x] 12-01-PLAN.md -- Conformer data models and per-conformer storage structure
- [x] 12-02-PLAN.md -- Canonical atom ordering (TDD)
- [x] 12-03-PLAN.md -- Backward compatibility verification and API schema updates

#### Phase 13: RDKit KDG Conformer Generation
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

Plans:
- [x] 13-01-PLAN.md -- KDG conformer generation and MMFF optimization (TDD)
- [x] 13-02-PLAN.md -- RMSD deduplication and energy window filtering (TDD)
- [x] 13-03-PLAN.md -- Pipeline orchestration, XYZ writing, and integration tests

#### Phase 14: Boltzmann Averaging Implementation
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

Plans:
- [x] 14-01-PLAN.md -- Boltzmann weight calculation with numerical stability (TDD)
- [x] 14-02-PLAN.md -- Ensemble NMR shift averaging and orchestration (TDD)

#### Phase 15: Multi-Conformer NWChem Integration
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

Plans:
- [x] 15-01-PLAN.md -- DFT energy extraction from NWChem output (TDD)
- [x] 15-02-PLAN.md -- Per-conformer DFT optimization loop with partial failure handling (TDD)
- [x] 15-03-PLAN.md -- Per-conformer NMR loop and full ensemble orchestrator (TDD)

#### Phase 16: CREST Integration (Optional High-Accuracy Mode)
**Goal**: Production-quality conformer generation using CREST metadynamics (when available)
**Depends on**: Phase 15
**Requirements**: CONF-04, CONF-05
**Success Criteria** (what must be TRUE):
  1. System generates conformers using CREST/xTB when binaries available on PATH
  2. System auto-detects CREST/xTB availability and reports status in API health endpoint
  3. CREST includes ALPB solvation model matching job solvent parameter
  4. CREST timeout (default 3600s) with automatic fallback to RDKit prevents hanging on macrocycles
  5. Environment variables (OMP_STACKSIZE, GFORTRAN_UNBUFFERED_ALL) configured for stability
**Plans**: TBD

Plans:
- [ ] 16-01: TBD

#### Phase 17: API Integration and Progress Tracking
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
**Plans**: TBD

Plans:
- [ ] 17-01: TBD

## Progress

**Execution Order:**
Phases execute in numeric order: 12 -> 13 -> 14 -> 15 -> 16 -> 17

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
| 16. CREST Integration | v2.0 | 0/? | Not started | - |
| 17. API Integration | v2.0 | 0/? | Not started | - |

---

*Last updated: 2026-01-27 after Phase 15 execution complete*
