# Roadmap: qm-nmr-calc

## Overview

This roadmap delivers a working async NMR prediction service across two milestones. Milestone v1.0 (phases 1-6) established the core service with ISiCLE/NWChem integration, REST API, job queue, and web UI. Milestone v1.1 (phases 7-11) replaces ISiCLE dependency, derives NWChem-specific scaling factors from DELTA50 benchmark data, and delivers accurate solvated NMR predictions with publication-quality results.

## Phases

**Phase Numbering:**
- Integer phases (1, 2, 3): Planned milestone work
- Decimal phases (2.1, 2.2): Urgent insertions (marked with INSERTED)

Decimal phases appear between their surrounding integers in numeric order.

### Milestone v1.0: Core NMR Service (Complete)

- [x] **Phase 1: Foundation** - Project structure, job queue, ISiCLE wrapper
- [x] **Phase 2: Input and API** - Molecule submission, validation, status polling
- [x] **Phase 3: NMR Calculations** - Chemical shift calculations and presets
- [x] **Phase 4: Results Delivery** - JSON results, file downloads, email notifications
- [x] **Phase 5: Visualization** - Spectrum plots and annotated structures
- [x] **Phase 6: Web UI** - Browser interface for submission and results

### Milestone v1.1: Accurate Chemical Shifts (Active)

- [ ] **Phase 7: NWChem Integration** - Direct NWChem I/O handling and COSMO solvation
- [ ] **Phase 8: DELTA50 Setup** - Benchmark dataset infrastructure
- [ ] **Phase 9: Benchmark Calculations** - Execute DELTA50 calculation matrix
- [ ] **Phase 10: Scaling Factors** - Derive and validate NWChem-specific scaling factors
- [ ] **Phase 11: Production Integration** - Apply scaling factors and new functional to production

## Phase Details

### Phase 1: Foundation
**Goal**: Establish working calculation infrastructure with job queue and ISiCLE integration
**Depends on**: Nothing (first phase)
**Requirements**: CALC-01, CALC-04, CALC-05
**Success Criteria** (what must be TRUE):
  1. Job can be queued and executed in background (Huey consumer picks up and runs it)
  2. ISiCLE/NWChem runs geometry optimization on a test molecule
  3. Failed calculations produce clear error status and message (not silent failure)
  4. Job state persists across process restarts (Huey/SQLite)
**Plans**: 3 plans

Plans:
- [x] 01-01-PLAN.md — Project setup with uv, Pydantic models, job storage layer
- [x] 01-02-PLAN.md — ISiCLE wrapper and Huey task queue with signal handlers
- [x] 01-03-PLAN.md — Startup validation, recovery logic, and integration test

### Phase 2: Input and API
**Goal**: Users can submit molecules and check job status via REST API
**Depends on**: Phase 1
**Requirements**: INP-01, INP-02, INP-03, INP-04, RES-04, UI-01
**Success Criteria** (what must be TRUE):
  1. User can POST a SMILES string and receive a job ID immediately
  2. User can POST a MOL/SDF file upload and receive a job ID immediately
  3. Invalid molecules return clear validation error (before queueing)
  4. User can GET job status (queued/running/complete/failed)
  5. OpenAPI/Swagger documentation is served at /docs
**Plans**: 3 plans

Plans:
- [x] 02-01-PLAN.md — FastAPI dependencies, molecule validation, API schemas
- [x] 02-02-PLAN.md — Job submission and health check routers
- [x] 02-03-PLAN.md — App assembly, startup script, integration tests

### Phase 3: NMR Calculations
**Goal**: System produces accurate NMR chemical shifts with configurable quality levels
**Depends on**: Phase 2
**Requirements**: CALC-02, CALC-03, CALC-06
**Success Criteria** (what must be TRUE):
  1. Completed job includes 1H NMR chemical shifts for all hydrogen atoms
  2. Completed job includes 13C NMR chemical shifts for all carbon atoms
  3. User can select calculation preset (draft/production) at submission time
  4. Different presets produce different calculation parameters (basis set, functional)
**Plans**: 3 plans

Plans:
- [x] 03-01-PLAN.md — Presets module, shielding-to-shift conversion, solvent validation
- [x] 03-02-PLAN.md — Extended models with NMR results, API schemas for preset/solvent
- [x] 03-03-PLAN.md — NMR calculation pipeline (geometry opt + shielding), API updates

### Phase 4: Results Delivery
**Goal**: Users can retrieve all calculation results in multiple formats
**Depends on**: Phase 3
**Requirements**: RES-01, RES-02, RES-03, NOTF-01
**Success Criteria** (what must be TRUE):
  1. User can GET chemical shifts as JSON with atom assignments (atom index to shift mapping)
  2. User can download optimized molecular geometry (XYZ or SDF format)
  3. User can download raw NWChem output files
  4. User can opt-in to email notification and receives email when job completes
**Plans**: 2 plans

Plans:
- [x] 04-01-PLAN.md — Results endpoint and file download endpoints (XYZ, SDF, raw output ZIP)
- [x] 04-02-PLAN.md — Email notification system (opt-in field, notifications module, signal handlers)

### Phase 5: Visualization
**Goal**: System generates visual representations of NMR results
**Depends on**: Phase 4
**Requirements**: RES-05, RES-06
**Success Criteria** (what must be TRUE):
  1. Completed job includes downloadable spectrum plot image (PNG showing peaks)
  2. Completed job includes annotated structure image (molecule with shift values on atoms)
  3. Spectrum plot shows peaks at correct chemical shift positions
**Plans**: 2 plans

Plans:
- [x] 05-01-PLAN.md — Visualization module (matplotlib spectrum plots, RDKit structure annotation)
- [x] 05-02-PLAN.md — Integration with tasks and API download endpoints

### Phase 6: Web UI
**Goal**: Users can interact with the service through a browser without API calls
**Depends on**: Phase 5
**Requirements**: UI-02, UI-03
**Success Criteria** (what must be TRUE):
  1. User can submit a molecule via web form (SMILES input or file upload)
  2. User can view job status page with auto-refresh while waiting
  3. User can view results page with spectrum plot, structure image, and download links
  4. Web UI is clean and presentable (not raw unstyled HTML)
**Plans**: 3 plans

Plans:
- [x] 06-01-PLAN.md — Template infrastructure, base layout, static files, web router setup
- [x] 06-02-PLAN.md — Submission form and status page with polling
- [x] 06-03-PLAN.md — Results page with image grid, downloads, and lightbox modal

### Phase 7: NWChem Integration
**Goal**: System directly handles NWChem I/O and enables accurate solvated calculations
**Depends on**: Phase 6
**Requirements**: NW-01, NW-02, NW-03, NW-04, NW-05, NW-06
**Success Criteria** (what must be TRUE):
  1. System generates valid NWChem input files for geometry optimization without ISiCLE runtime
  2. System generates valid NWChem input files for NMR shielding with COSMO solvation
  3. System parses NWChem output files to extract shielding tensors for H and C atoms
  4. User can submit pre-optimized geometry and system skips geometry optimization step
  5. ISiCLE attribution is visible in code comments and documentation where code adapted
**Plans**: 3-4 plans

### Phase 8: DELTA50 Setup
**Goal**: DELTA50 benchmark infrastructure ready to execute calculation matrix
**Depends on**: Phase 7
**Requirements**: BENCH-01, BENCH-02, BENCH-03
**Success Criteria** (what must be TRUE):
  1. All 50 DELTA50 molecule structures are available in SDF or XYZ format
  2. Experimental 1H and 13C shift values are stored with molecule ID associations
  3. Benchmark runner can execute calculations for all combinations: 50 molecules x 2 functionals x 2 solvents
  4. Benchmark runner outputs raw calculation results in organized directory structure
**Plans**: 2-3 plans

### Phase 9: Benchmark Calculations
**Goal**: Complete DELTA50 calculation matrix for scaling factor derivation
**Depends on**: Phase 8
**Requirements**: None (execution phase using Phase 8 infrastructure)
**Success Criteria** (what must be TRUE):
  1. All 50 DELTA50 molecules calculated with B3LYP in CHCl3 (50 calculations)
  2. All 50 DELTA50 molecules calculated with B3LYP in DMSO (50 calculations)
  3. All 50 DELTA50 molecules calculated with WP04 in CHCl3 (50 calculations for 1H)
  4. All 50 DELTA50 molecules calculated with WP04 in DMSO (50 calculations for 1H)
  5. Failed calculations documented with error analysis (if any)
**Plans**: 2 plans

### Phase 10: Scaling Factors
**Goal**: Derive and validate NWChem-specific scaling factors from benchmark data
**Depends on**: Phase 9
**Requirements**: SCALE-01, SCALE-02, SCALE-03, SCALE-04
**Success Criteria** (what must be TRUE):
  1. Scaling factors (slope, intercept) derived via linear regression for B3LYP/CHCl3
  2. Scaling factors derived for B3LYP/DMSO, WP04/CHCl3, WP04/DMSO
  3. Each scaling factor set validated with MAE and RMSD statistics vs experimental shifts
  4. Scaling factors stored in code-accessible format (JSON, Python constants, or config file)
  5. Validation report shows accuracy improvement over current CHESHIRE factors
**Plans**: 2-3 plans

### Phase 11: Production Integration
**Goal**: Production calculations use NWChem-derived factors and enable WP04 functional
**Depends on**: Phase 10
**Requirements**: PROD-01, PROD-02, PROD-03, PROD-04
**Success Criteria** (what must be TRUE):
  1. Production calculations apply COSMO solvation when solvent is specified (bug fixed)
  2. Production calculations use NWChem-derived scaling factors instead of CHESHIRE
  3. User can select WP04 functional at submission time (new preset or parameter)
  4. ISiCLE is no longer a runtime dependency (removed from production imports and dependencies)
  5. Production calculation with WP04/CHCl3 produces more accurate 1H shifts than B3LYP/gas
**Plans**: 2-3 plans

## Progress

**Execution Order:**
Phases execute in numeric order: 1 -> 2 -> 3 -> 4 -> 5 -> 6 -> 7 -> 8 -> 9 -> 10 -> 11

### Milestone v1.0: Core NMR Service

| Phase | Plans Complete | Status | Completed |
|-------|----------------|--------|-----------|
| 1. Foundation | 3/3 | Complete | 2026-01-19 |
| 2. Input and API | 3/3 | Complete | 2026-01-19 |
| 3. NMR Calculations | 3/3 | Complete | 2026-01-19 |
| 4. Results Delivery | 2/2 | Complete | 2026-01-19 |
| 5. Visualization | 2/2 | Complete | 2026-01-20 |
| 6. Web UI | 3/3 | Complete | 2026-01-20 |

### Milestone v1.1: Accurate Chemical Shifts

| Phase | Plans Complete | Status | Completed |
|-------|----------------|--------|-----------|
| 7. NWChem Integration | 0/3-4 | Pending | — |
| 8. DELTA50 Setup | 0/2-3 | Pending | — |
| 9. Benchmark Calculations | 0/2 | Pending | — |
| 10. Scaling Factors | 0/2-3 | Pending | — |
| 11. Production Integration | 0/2-3 | Pending | — |
