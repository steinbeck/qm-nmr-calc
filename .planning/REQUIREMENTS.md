# Requirements: qm-nmr-calc

**Defined:** 2026-01-19 (v1.0), 2026-01-21 (v1.1)
**Core Value:** Reliable async NMR predictions with full control over calculation parameters — submit a molecule, get back accurate ¹H/¹³C shifts without babysitting long-running calculations.

## v1.0 Requirements ✓ (Complete)

All v1.0 requirements delivered 2026-01-20.

### Input ✓

- [x] **INP-01**: User can submit molecule via SMILES string
- [x] **INP-02**: User can upload molecule as MOL/SDF file
- [x] **INP-03**: System validates molecule structure before queueing
- [x] **INP-04**: System returns job ID immediately on submission

### Calculation ✓

- [x] **CALC-01**: System runs geometry optimization via ISiCLE/NWChem
- [x] **CALC-02**: System calculates ¹H NMR chemical shifts
- [x] **CALC-03**: System calculates ¹³C NMR chemical shifts
- [x] **CALC-04**: Jobs queue for background processing
- [x] **CALC-05**: System handles calculation failures gracefully (status, error message)
- [x] **CALC-06**: User can select calculation preset (fast/standard/publication)

### Results ✓

- [x] **RES-01**: User can retrieve chemical shifts as JSON with atom assignments
- [x] **RES-02**: User can retrieve optimized molecular geometry
- [x] **RES-03**: User can download raw NWChem output files
- [x] **RES-04**: User can poll for job status
- [x] **RES-05**: System generates visual spectrum plot
- [x] **RES-06**: System generates annotated structure drawing with shifts on atoms

### Notifications ✓

- [x] **NOTF-01**: User can opt-in to email notification on job completion

### Interface ✓

- [x] **UI-01**: REST API with OpenAPI/Swagger documentation
- [x] **UI-02**: Web UI for submitting jobs
- [x] **UI-03**: Web UI for viewing job status and results

## v1.1 Requirements (Active)

Requirements for accurate chemical shift predictions with NWChem-derived scaling factors.

### NWChem Integration

- [ ] **NW-01**: System generates NWChem input files for geometry optimization (adapted from ISiCLE)
- [ ] **NW-02**: System generates NWChem input files for NMR shielding calculation (adapted from ISiCLE)
- [ ] **NW-03**: System parses NWChem output to extract shielding tensors (adapted from ISiCLE)
- [ ] **NW-04**: System supports COSMO solvation model with CHCl3 and DMSO
- [ ] **NW-05**: System accepts pre-optimized XYZ/SDF geometries (skip geometry optimization)
- [ ] **NW-06**: ISiCLE attribution included in code/docs where code adapted

### Benchmark Infrastructure

- [ ] **BENCH-01**: DELTA50 dataset structures available in standard format
- [ ] **BENCH-02**: DELTA50 experimental shifts stored with molecule associations
- [ ] **BENCH-03**: Benchmark runner can execute calculations for dataset × method × solvent

### Scaling Factors

- [ ] **SCALE-01**: System derives scaling factors via linear regression from benchmark data
- [ ] **SCALE-02**: Scaling factors validated with MAE/RMSD statistics
- [ ] **SCALE-03**: Scaling factors stored for B3LYP in CHCl3 and DMSO
- [ ] **SCALE-04**: Scaling factors stored for WP04 (¹H) in CHCl3 and DMSO

### Production Integration

- [ ] **PROD-01**: Production calculations use COSMO solvation (fix current bug)
- [ ] **PROD-02**: Production calculations use NWChem-derived scaling factors
- [ ] **PROD-03**: User can select WP04 functional for improved ¹H accuracy
- [ ] **PROD-04**: ISiCLE dependency removed from production code

## v2 Requirements

Deferred to future release. Tracked but not in current roadmap.

### Calculation

- **CALC-ADV**: Advanced parameter control (custom basis sets, functionals, solvation models)
- **CALC-CONF**: Conformer averaging for flexible molecules

### Input

- **INP-PREVIEW**: Molecule preview before submission

### Results

- **RES-SHARE**: Shareable result links

### Nuclei

- **NUC-15N**: ¹⁵N NMR chemical shift predictions
- **NUC-19F**: ¹⁹F NMR chemical shift predictions
- **NUC-31P**: ³¹P NMR chemical shift predictions

## Out of Scope

Explicitly excluded. Documented to prevent scope creep.

| Feature | Reason |
|---------|--------|
| Multi-user authentication | Single user for now, architecture ready for later |
| Other nuclei (¹⁵N, ³¹P, ¹⁹F) | ¹H/¹³C only for v1.x |
| Batch submission UI | API could support later |
| Real-time WebSocket updates | Polling sufficient |
| Public benchmark API | Benchmark tooling is internal only |
| Mobile interface | Web UI is desktop-focused |

## Traceability

Which phases cover which requirements.

| Requirement | Phase | Status |
|-------------|-------|--------|
| NW-01 | Phase 7 | Pending |
| NW-02 | Phase 7 | Pending |
| NW-03 | Phase 7 | Pending |
| NW-04 | Phase 7 | Pending |
| NW-05 | Phase 7 | Pending |
| NW-06 | Phase 7 | Pending |
| BENCH-01 | Phase 8 | Pending |
| BENCH-02 | Phase 8 | Pending |
| BENCH-03 | Phase 8 | Pending |
| SCALE-01 | Phase 10 | Pending |
| SCALE-02 | Phase 10 | Pending |
| SCALE-03 | Phase 10 | Pending |
| SCALE-04 | Phase 10 | Pending |
| PROD-01 | Phase 11 | Pending |
| PROD-02 | Phase 11 | Pending |
| PROD-03 | Phase 11 | Pending |
| PROD-04 | Phase 11 | Pending |

**Coverage:**
- v1.1 requirements: 17 total
- Mapped to phases: 17
- Unmapped: 0 ✓

---
*Requirements defined: 2026-01-19 (v1.0), 2026-01-21 (v1.1)*
*Last updated: 2026-01-21 after v1.1 milestone definition*
