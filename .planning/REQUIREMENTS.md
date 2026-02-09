# Requirements: qm-nmr-calc v2.8

**Defined:** 2026-02-07
**Core Value:** Reliable async NMR predictions with full control over calculation parameters -- submit a molecule, get back accurate 1H/13C shifts without babysitting long-running calculations.

## v2.8 Requirements

Requirements for expanded solvent support. Each maps to roadmap phases.

### Benchmark Infrastructure

- [x] **BENCH-01**: Benchmark CLI accepts methanol, water, acetone, benzene as valid solvents
- [x] **BENCH-02**: DELTA50 benchmark runs for all 50 molecules in methanol solvent
- [x] **BENCH-03**: DELTA50 benchmark runs for all 50 molecules in water solvent
- [x] **BENCH-04**: DELTA50 benchmark runs for all 50 molecules in acetone solvent
- [x] **BENCH-05**: DELTA50 benchmark runs for all 50 molecules in benzene solvent
- [x] **BENCH-06**: OLS scaling factors derived for 1H and 13C in each new solvent (8 new factor sets)

### Code Integration

- [x] **INTG-01**: User can select methanol as solvent in web UI and API
- [x] **INTG-02**: User can select water as solvent in web UI and API
- [x] **INTG-03**: User can select acetone as solvent in web UI and API
- [x] **INTG-04**: User can select benzene as solvent in web UI and API
- [x] **INTG-05**: NWChem COSMO generates correct solvation for all 4 new solvents
- [x] **INTG-06**: Scaling factors loaded and applied correctly for all new solvents

### Validation

- [x] **VALID-01**: All 8 new scaling factor sets have R² > 0.99
- [x] **VALID-02**: 1H MAE is below 0.2 ppm for each new solvent
- [x] **VALID-03**: 13C MAE is below 3.0 ppm for each new solvent

### Documentation

- [ ] **DOCS-01**: SCALING-FACTORS.md updated with new solvent factor statistics
- [ ] **DOCS-02**: README updated to list 7 supported solvents

## Future Requirements

Deferred to later milestones.

### Additional Solvents

- **SOLV-01**: Pyridine-d5 COSMO solvation and scaling factors
- **SOLV-02**: THF-d8 COSMO solvation and scaling factors
- **SOLV-03**: Toluene-d8 COSMO solvation and scaling factors

## Out of Scope

| Feature | Reason |
|---------|--------|
| WP04 functional factors for new solvents | NWChem doesn't support WP04 without custom compilation |
| Experimental shifts re-measurement per solvent | Standard practice uses CDCl3 reference set for all solvents |
| Solvent-specific conformer generation | COSMO handles solvation at DFT level, conformers generated in vacuum/MMFF |

## Traceability

| Requirement | Phase | Status |
|-------------|-------|--------|
| BENCH-01 | 54 | Complete |
| BENCH-02 | 55 | Complete |
| BENCH-03 | 55 | Complete |
| BENCH-04 | 55 | Complete |
| BENCH-05 | 55 | Complete |
| BENCH-06 | 56 | Complete |
| INTG-01 | 57 | Complete |
| INTG-02 | 57 | Complete |
| INTG-03 | 57 | Complete |
| INTG-04 | 57 | Complete |
| INTG-05 | 54 | Complete |
| INTG-06 | 57 | Complete |
| VALID-01 | 56 | Complete |
| VALID-02 | 56 | Complete |
| VALID-03 | 56 | Complete |
| DOCS-01 | 58 | Pending |
| DOCS-02 | 58 | Pending |

**Coverage:**
- v2.8 requirements: 17 total
- Mapped to phases: 17
- Unmapped: 0

---
*Requirements defined: 2026-02-07*
*Last updated: 2026-02-09 — INTG-01, INTG-02, INTG-03, INTG-04, INTG-06 complete (Phase 57)*
