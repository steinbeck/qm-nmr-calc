# Requirements: qm-nmr-calc v2.9

**Defined:** 2026-02-10
**Core Value:** Reliable async NMR predictions with full control over calculation parameters -- submit a molecule, get back accurate 1H/13C shifts without babysitting long-running calculations.

## v2.9 Requirements

Requirements for extended solvent coverage. Same benchmark-derive-integrate pipeline as v2.8, applied to 6 additional solvents. Each maps to roadmap phases.

### Benchmark Infrastructure

- [x] **BENCH-01**: Benchmark CLI accepts pyridine, thf, toluene, dcm, acetonitrile, dmf as valid solvents
- [x] **BENCH-02**: DELTA50 benchmark runs for all 50 molecules in pyridine solvent
- [x] **BENCH-03**: DELTA50 benchmark runs for all 50 molecules in THF solvent
- [x] **BENCH-04**: DELTA50 benchmark runs for all 50 molecules in toluene solvent
- [x] **BENCH-05**: DELTA50 benchmark runs for all 50 molecules in DCM solvent
- [x] **BENCH-06**: DELTA50 benchmark runs for all 50 molecules in acetonitrile solvent
- [x] **BENCH-07**: DELTA50 benchmark runs for all 50 molecules in DMF solvent
- [x] **BENCH-08**: OLS scaling factors derived for 1H and 13C in each new solvent (12 new factor sets)

### Validation

- [x] **VALID-01**: All 12 new scaling factor sets have R² > 0.99
- [x] **VALID-02**: 1H MAE is below 0.2 ppm for each new solvent
- [x] **VALID-03**: 13C MAE is below 3.0 ppm for each new solvent

### Code Integration

- [ ] **INTG-01**: User can select pyridine as solvent in web UI and API
- [ ] **INTG-02**: User can select THF as solvent in web UI and API
- [ ] **INTG-03**: User can select toluene as solvent in web UI and API
- [ ] **INTG-04**: User can select DCM as solvent in web UI and API
- [ ] **INTG-05**: User can select acetonitrile as solvent in web UI and API
- [ ] **INTG-06**: User can select DMF as solvent in web UI and API
- [ ] **INTG-07**: NWChem COSMO generates correct solvation for all 6 new solvents (including acetntrl mapping)
- [ ] **INTG-08**: Scaling factors loaded and applied correctly for all new solvents

### Documentation

- [ ] **DOCS-01**: SCALING-FACTORS.md updated with all 26 factor sets (13 solvents x 2 nuclei)
- [ ] **DOCS-02**: README updated to list 13 supported solvents

## Future Requirements

Deferred to later milestones.

### Publication

- **PUB-01**: DELTA50 benchmark dataset paper describing methodology and results
- **PUB-02**: Public dataset release with raw NWChem outputs

## Out of Scope

| Feature | Reason |
|---------|--------|
| WP04 functional factors for new solvents | NWChem doesn't support WP04 without custom compilation |
| Experimental shifts re-measurement per solvent | Standard practice uses CDCl3 reference set for all solvents |
| Solvent-specific conformer generation | COSMO handles solvation at DFT level, conformers generated in vacuum/MMFF |

## Traceability

| Requirement | Phase | Status |
|-------------|-------|--------|
| BENCH-01 | Phase 59 | Complete |
| BENCH-02 | Phase 60 | Complete |
| BENCH-03 | Phase 60 | Complete |
| BENCH-04 | Phase 61 | Complete |
| BENCH-05 | Phase 61 | Complete |
| BENCH-06 | Phase 62 | Complete |
| BENCH-07 | Phase 62 | Complete |
| BENCH-08 | Phase 63 | Complete |
| VALID-01 | Phase 63 | Complete |
| VALID-02 | Phase 63 | Complete |
| VALID-03 | Phase 63 | Complete |
| INTG-01 | Phase 64 | Pending |
| INTG-02 | Phase 64 | Pending |
| INTG-03 | Phase 64 | Pending |
| INTG-04 | Phase 64 | Pending |
| INTG-05 | Phase 64 | Pending |
| INTG-06 | Phase 64 | Pending |
| INTG-07 | Phase 64 | Pending |
| INTG-08 | Phase 64 | Pending |
| DOCS-01 | Phase 65 | Pending |
| DOCS-02 | Phase 65 | Pending |

**Coverage:**
- v2.9 requirements: 22 total
- Mapped to phases: 22
- Unmapped: 0

---
*Requirements defined: 2026-02-10*
*Last updated: 2026-02-10 after roadmap creation*
