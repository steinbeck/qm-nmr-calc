# Requirements: qm-nmr-calc v2.3 NMReData Export

**Defined:** 2026-02-01
**Core Value:** Enable machine-readable export of NMR prediction results in the NMReData standard format.

## v2.3 Requirements

Requirements for NMReData export feature. Each maps to roadmap phases.

### NMReData Format

- [x] **NMRD-01**: Generated file contains valid MOL block with optimized 3D coordinates
- [x] **NMRD-02**: File includes `NMREDATA_VERSION` tag (1.1 or 2.0)
- [x] **NMRD-03**: File includes `NMREDATA_LEVEL` tag (0 for predicted data)
- [x] **NMRD-04**: File includes `NMREDATA_SOLVENT` tag with proper mapping (CHCl3→CDCl3, DMSO→(CD3)2SO, vacuum→vacuum)
- [x] **NMRD-05**: File includes `NMREDATA_TEMPERATURE` tag (298.15 K or ensemble temperature)
- [x] **NMRD-06**: File includes `NMREDATA_ASSIGNMENT` tag with 1H and 13C chemical shifts mapped to 1-indexed atom numbers
- [x] **NMRD-07**: File includes `NMREDATA_FORMULA` tag with molecular formula
- [x] **NMRD-08**: File includes `NMREDATA_SMILES` tag with input SMILES
- [x] **NMRD-09**: File includes provenance metadata (calculation method, basis set, scaling factors applied)

### REST API

- [ ] **API-01**: User can download NMReData file via `GET /api/v1/jobs/{job_id}/nmredata.sdf`
- [ ] **API-02**: Response includes proper HTTP headers (media type `chemical/x-mdl-sdfile`, filename `{job_id}_nmredata.sdf`)
- [ ] **API-03**: Endpoint returns 404 if job not found, 409 if job not complete

### Web UI

- [ ] **UI-01**: Results page includes download button for NMReData file alongside existing PNG/XYZ downloads

### Testing

- [ ] **TEST-01**: Unit tests validate NMReData tag formatting (separators, atom numbering)
- [ ] **TEST-02**: Integration tests validate complete endpoint flow
- [ ] **TEST-03**: Exported file can be parsed by RDKit SDMolSupplier and atom assignments verified

## Future Requirements

Deferred to later milestones.

### Enhanced NMReData

- **NMRD-10**: File includes `NMREDATA_INCHI` tag for structure matching
- **NMRD-11**: File includes `NMREDATA_1D_1H` and `NMREDATA_1D_13C` pseudo-spectrum tags for visualization tool compatibility
- **NMRD-12**: Per-conformer export option for ensemble calculations

### J-coupling

- **NMRD-13**: File includes `NMREDATA_J` tag with calculated coupling constants (requires additional computation)

## Out of Scope

Explicitly excluded from v2.3. Documented to prevent scope creep.

| Feature | Reason |
|---------|--------|
| J-coupling constants | Not computed by qm-nmr-calc; would require additional QM calculation |
| 2D correlation data | Experimental technique only; not applicable to predictions |
| Per-conformer export | NMReData format has limited multi-conformer support; use Boltzmann-averaged shifts |
| External validation tools | Java parser requires JRE; RDKit round-trip sufficient for v2.3 |

## Traceability

Which phases cover which requirements. Updated during roadmap creation.

| Requirement | Phase | Status |
|-------------|-------|--------|
| NMRD-01 | Phase 32 | Complete |
| NMRD-02 | Phase 32 | Complete |
| NMRD-03 | Phase 32 | Complete |
| NMRD-04 | Phase 32 | Complete |
| NMRD-05 | Phase 32 | Complete |
| NMRD-06 | Phase 32 | Complete |
| NMRD-07 | Phase 32 | Complete |
| NMRD-08 | Phase 32 | Complete |
| NMRD-09 | Phase 32 | Complete |
| API-01 | Phase 33 | Pending |
| API-02 | Phase 33 | Pending |
| API-03 | Phase 33 | Pending |
| UI-01 | Phase 33 | Pending |
| TEST-01 | Phase 34 | Pending |
| TEST-02 | Phase 34 | Pending |
| TEST-03 | Phase 34 | Pending |

**Coverage:**
- v2.3 requirements: 16 total
- Mapped to phases: 16
- Unmapped: 0 ✓

---
*Requirements defined: 2026-02-01*
*Last updated: 2026-02-01 after Phase 32 completion*
