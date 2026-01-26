# Requirements: qm-nmr-calc v2.0

**Defined:** 2026-01-26
**Core Value:** Reliable async NMR predictions with full control over calculation parameters -- submit a molecule, get back accurate 1H/13C shifts without babysitting long-running calculations.

## v2.0 Requirements

Requirements for v2.0 Conformational Sampling. Each maps to roadmap phases.

### Conformer Generation

- [ ] **CONF-01**: User can choose between single-conformer mode (v1.x behavior) and ensemble mode when submitting a job
- [ ] **CONF-02**: User can choose RDKit KDG or CREST as conformer generation method in ensemble mode
- [ ] **CONF-03**: System generates conformers using RDKit KDG (pure distance geometry, no crystal bias)
- [ ] **CONF-04**: System generates conformers using CREST/xTB when CREST binary is available on PATH
- [ ] **CONF-05**: System auto-detects CREST/xTB availability and reports it in API health/status
- [ ] **CONF-06**: App works fully without CREST/xTB installed (RDKit-only mode)

### Energy Filtering

- [ ] **FILT-01**: System applies pre-DFT energy window filter (default 6 kcal/mol) to conformer set
- [ ] **FILT-02**: System applies post-DFT energy window filter (default 3 kcal/mol) after geometry optimization
- [ ] **FILT-03**: User can configure energy window parameters (pre-DFT and post-DFT thresholds)

### DFT Calculations

- [ ] **DFT-01**: System runs DFT geometry optimization on each conformer surviving pre-DFT filter
- [ ] **DFT-02**: System runs NMR shielding calculation on each conformer surviving post-DFT filter
- [ ] **DFT-03**: Each conformer uses isolated scratch directory to prevent NWChem file conflicts
- [ ] **DFT-04**: System extracts DFT energies from optimization step for Boltzmann weighting
- [ ] **DFT-05**: System handles partial conformer failures gracefully (continues with successful conformers)

### Boltzmann Averaging

- [ ] **BOLTZ-01**: System calculates Boltzmann weights from DFT optimization energies
- [ ] **BOLTZ-02**: System computes population-weighted average chemical shifts across conformer ensemble
- [ ] **BOLTZ-03**: Boltzmann implementation is numerically stable (handles wide energy ranges without overflow/underflow)
- [ ] **BOLTZ-04**: User can set temperature parameter for Boltzmann weighting (default 298.15 K)

### API and Output

- [ ] **API-01**: API returns Boltzmann-weighted average shifts (not per-conformer detail)
- [ ] **API-02**: API response includes ensemble metadata (conformer count, energy range, top 3 populations, method used, temperature)
- [ ] **API-03**: 3D viewer shows lowest-energy conformer geometry with shift labels
- [ ] **API-04**: Existing v1.x single-conformer API behavior preserved when conformer_mode=single (backward compatible)

### Progress Tracking

- [ ] **PROG-01**: Job status includes ensemble-specific states (generating_conformers, optimizing_conformers X/N, calculating_nmr X/N, averaging_shifts)
- [ ] **PROG-02**: Web UI displays conformer progress for ensemble jobs

## v2.x Requirements

Deferred to future release. Tracked but not in current roadmap.

### Performance

- **PERF-01**: Parallel conformer DFT processing within single Huey task
- **PERF-02**: Conformer ensemble caching by (SMILES, method, energy_window)
- **PERF-03**: Adaptive post-DFT filtering by cumulative Boltzmann population (95% threshold)

### Enhanced Output

- **OUT-01**: Per-conformer detail export endpoint (/jobs/{job_id}/conformers)
- **OUT-02**: Rotatable bond count display in web UI
- **OUT-03**: Multi-temperature Boltzmann weighting

### Accuracy

- **ACC-01**: GFN2-xTB quasi-RRHO entropy corrections for Boltzmann weighting
- **ACC-02**: Solvent-specific CREST conformer generation (ALPB matching NMR solvent)

## Out of Scope

Explicitly excluded. Documented to prevent scope creep.

| Feature | Reason |
|---------|--------|
| DP4+ probability scoring | Requires experimental spectrum input, scope change for v3.0 |
| Real-time WebSocket progress | Marginal UX improvement over polling, adds infrastructure complexity |
| Automatic rigidity detection | User explicit choice sufficient; heuristics misclassify too often |
| Per-conformer shifts in default API | Overwhelming data, rarely actionable for typical use case |
| MMFF energy Boltzmann weighting | Anticorrelates with DFT, would give wrong conformer populations |
| Automatic method selection by molecule size | Size is poor proxy for flexibility |
| WP04 functional | NWChem doesn't support WP04 without custom compilation |

## Traceability

Which phases cover which requirements. Updated during roadmap creation.

| Requirement | Phase | Status |
|-------------|-------|--------|
| CONF-01 | Phase 17 | Pending |
| CONF-02 | Phase 17 | Pending |
| CONF-03 | Phase 13 | Pending |
| CONF-04 | Phase 16 | Pending |
| CONF-05 | Phase 16 | Pending |
| CONF-06 | Phase 13 | Pending |
| FILT-01 | Phase 13 | Pending |
| FILT-02 | Phase 15 | Pending |
| FILT-03 | Phase 15 | Pending |
| DFT-01 | Phase 15 | Pending |
| DFT-02 | Phase 15 | Pending |
| DFT-03 | Phase 12 | Pending |
| DFT-04 | Phase 15 | Pending |
| DFT-05 | Phase 15 | Pending |
| BOLTZ-01 | Phase 14 | Pending |
| BOLTZ-02 | Phase 14 | Pending |
| BOLTZ-03 | Phase 14 | Pending |
| BOLTZ-04 | Phase 14 | Pending |
| API-01 | Phase 17 | Pending |
| API-02 | Phase 17 | Pending |
| API-03 | Phase 17 | Pending |
| API-04 | Phase 12 | Pending |
| PROG-01 | Phase 17 | Pending |
| PROG-02 | Phase 17 | Pending |

**Coverage:**
- v2.0 requirements: 24 total
- Mapped to phases: 24 (100%)
- Unmapped: 0

**Phase distribution:**
- Phase 12 (Data Model): 2 requirements
- Phase 13 (RDKit Generation): 3 requirements
- Phase 14 (Boltzmann): 4 requirements
- Phase 15 (NWChem Integration): 6 requirements
- Phase 16 (CREST): 2 requirements
- Phase 17 (API): 7 requirements

---
*Requirements defined: 2026-01-26*
*Last updated: 2026-01-26 after roadmap creation*
