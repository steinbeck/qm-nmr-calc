# Requirements: qm-nmr-calc

**Defined:** 2026-01-19
**Core Value:** Reliable async NMR predictions with full control over calculation parameters — submit a molecule, get back accurate ¹H/¹³C shifts without babysitting long-running calculations.

## v1 Requirements

Requirements for initial release. Each maps to roadmap phases.

### Input

- [ ] **INP-01**: User can submit molecule via SMILES string
- [ ] **INP-02**: User can upload molecule as MOL/SDF file
- [ ] **INP-03**: System validates molecule structure before queueing
- [ ] **INP-04**: System returns job ID immediately on submission

### Calculation

- [ ] **CALC-01**: System runs geometry optimization via ISiCLE/NWChem
- [ ] **CALC-02**: System calculates ¹H NMR chemical shifts
- [ ] **CALC-03**: System calculates ¹³C NMR chemical shifts
- [ ] **CALC-04**: Jobs queue for background processing
- [ ] **CALC-05**: System handles calculation failures gracefully (status, error message)
- [ ] **CALC-06**: User can select calculation preset (fast/standard/publication)

### Results

- [ ] **RES-01**: User can retrieve chemical shifts as JSON with atom assignments
- [ ] **RES-02**: User can retrieve optimized molecular geometry
- [ ] **RES-03**: User can download raw NWChem output files
- [ ] **RES-04**: User can poll for job status
- [ ] **RES-05**: System generates visual spectrum plot
- [ ] **RES-06**: System generates annotated structure drawing with shifts on atoms

### Notifications

- [ ] **NOTF-01**: User can opt-in to email notification on job completion

### Interface

- [ ] **UI-01**: REST API with OpenAPI/Swagger documentation
- [ ] **UI-02**: Web UI for submitting jobs
- [ ] **UI-03**: Web UI for viewing job status and results

## v2 Requirements

Deferred to future release. Tracked but not in current roadmap.

### Calculation

- **CALC-ADV**: Advanced parameter control (custom basis sets, functionals, solvation models)
- **CALC-CONF**: Conformer averaging for flexible molecules

### Input

- **INP-PREVIEW**: Molecule preview before submission

### Results

- **RES-SHARE**: Shareable result links

## Out of Scope

Explicitly excluded. Documented to prevent scope creep.

| Feature | Reason |
|---------|--------|
| Multi-user authentication | Single user for v1, architecture ready for later |
| Other nuclei (¹⁵N, ³¹P, ¹⁹F) | ¹H/¹³C only for v1 |
| Batch submission UI | API could support later |
| Real-time WebSocket updates | Polling sufficient for v1 |
| Spectrum processing/editing | NMRium's domain, we display only |
| Mobile interface | Web UI is desktop-focused |

## Traceability

Which phases cover which requirements. Updated during roadmap creation.

| Requirement | Phase | Status |
|-------------|-------|--------|
| INP-01 | Phase 2 | Complete |
| INP-02 | Phase 2 | Complete |
| INP-03 | Phase 2 | Complete |
| INP-04 | Phase 2 | Complete |
| CALC-01 | Phase 1 | Complete |
| CALC-02 | Phase 3 | Complete |
| CALC-03 | Phase 3 | Complete |
| CALC-04 | Phase 1 | Complete |
| CALC-05 | Phase 1 | Complete |
| CALC-06 | Phase 3 | Complete |
| RES-01 | Phase 4 | Complete |
| RES-02 | Phase 4 | Complete |
| RES-03 | Phase 4 | Complete |
| RES-04 | Phase 2 | Complete |
| RES-05 | Phase 5 | Pending |
| RES-06 | Phase 5 | Pending |
| NOTF-01 | Phase 4 | Complete |
| UI-01 | Phase 2 | Complete |
| UI-02 | Phase 6 | Pending |
| UI-03 | Phase 6 | Pending |

**Coverage:**
- v1 requirements: 20 total
- Mapped to phases: 20
- Unmapped: 0

---
*Requirements defined: 2026-01-19*
*Last updated: 2026-01-19 after roadmap creation*
