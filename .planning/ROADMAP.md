# Roadmap: qm-nmr-calc

## Overview

This roadmap delivers a working async NMR prediction service in 6 phases. We start with the calculation foundation (ISiCLE/NWChem integration and job queue), then build the API layer for input and status, add NMR calculations and presets, deliver results with notifications, generate visualizations, and finally wrap it in a web UI. Each phase delivers a coherent capability that builds on the previous.

## Phases

**Phase Numbering:**
- Integer phases (1, 2, 3): Planned milestone work
- Decimal phases (2.1, 2.2): Urgent insertions (marked with INSERTED)

Decimal phases appear between their surrounding integers in numeric order.

- [x] **Phase 1: Foundation** - Project structure, job queue, ISiCLE wrapper
- [x] **Phase 2: Input and API** - Molecule submission, validation, status polling
- [x] **Phase 3: NMR Calculations** - Chemical shift calculations and presets
- [ ] **Phase 4: Results Delivery** - JSON results, file downloads, email notifications
- [ ] **Phase 5: Visualization** - Spectrum plots and annotated structures
- [ ] **Phase 6: Web UI** - Browser interface for submission and results

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
- [ ] 04-01-PLAN.md — Results endpoint and file download endpoints (XYZ, SDF, raw output ZIP)
- [ ] 04-02-PLAN.md — Email notification system (opt-in field, notifications module, signal handlers)

### Phase 5: Visualization
**Goal**: System generates visual representations of NMR results
**Depends on**: Phase 4
**Requirements**: RES-05, RES-06
**Success Criteria** (what must be TRUE):
  1. Completed job includes downloadable spectrum plot image (PNG showing peaks)
  2. Completed job includes annotated structure image (molecule with shift values on atoms)
  3. Spectrum plot shows peaks at correct chemical shift positions
**Plans**: TBD

Plans:
- [ ] 05-01: TBD

### Phase 6: Web UI
**Goal**: Users can interact with the service through a browser without API calls
**Depends on**: Phase 5
**Requirements**: UI-02, UI-03
**Success Criteria** (what must be TRUE):
  1. User can submit a molecule via web form (SMILES input or file upload)
  2. User can view job status page with auto-refresh while waiting
  3. User can view results page with spectrum plot, structure image, and download links
  4. Web UI is clean and presentable (not raw unstyled HTML)
**Plans**: TBD

Plans:
- [ ] 06-01: TBD

## Progress

**Execution Order:**
Phases execute in numeric order: 1 -> 2 -> 3 -> 4 -> 5 -> 6

| Phase | Plans Complete | Status | Completed |
|-------|----------------|--------|-----------|
| 1. Foundation | 3/3 | Complete | 2026-01-19 |
| 2. Input and API | 3/3 | Complete | 2026-01-19 |
| 3. NMR Calculations | 3/3 | Complete | 2026-01-19 |
| 4. Results Delivery | 0/2 | Not started | - |
| 5. Visualization | 0/? | Not started | - |
| 6. Web UI | 0/? | Not started | - |
