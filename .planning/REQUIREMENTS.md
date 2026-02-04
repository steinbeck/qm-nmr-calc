# Requirements: qm-nmr-calc v2.5

**Defined:** 2026-02-04
**Core Value:** Reliable async NMR predictions with full control over calculation parameters

## v2.5 Requirements

Requirements for ARM64 Docker support. Each maps to roadmap phases.

### Container Build

- [x] **CONT-01**: ARM64 worker container runs NWChem DFT geometry optimization
- [x] **CONT-02**: ARM64 worker container runs NWChem NMR shielding calculation
- [x] **CONT-03**: ARM64 worker container runs xTB energy calculations
- [x] **CONT-04**: ARM64 worker container runs CREST conformer search
- [x] **CONT-05**: ARM64 worker uses conda-forge packages (not x86 binaries)
- [x] **CONT-06**: ARM64 worker produces numerically equivalent results to x86_64

### CI/CD

- [ ] **CICD-01**: GitHub Actions builds ARM64 worker image with native runner
- [ ] **CICD-02**: Multi-arch manifest created combining amd64 and arm64
- [ ] **CICD-03**: Single image tag (`latest`, `vX.Y.Z`) works on both architectures
- [ ] **CICD-04**: ARM64 build uses native runner (not QEMU emulation)

### Documentation

- [ ] **DOCS-01**: README documents ARM64/Apple Silicon support
- [ ] **DOCS-02**: Known issues section covers ARM64-specific caveats
- [ ] **DOCS-03**: docker-compose.yml works unchanged on ARM64

## Future Requirements

Deferred to later milestones.

### Performance Optimization

- **PERF-01**: OpenBLAS TARGET optimization for specific ARM hardware (Graviton, Apple M-series)
- **PERF-02**: Performance benchmarking ARM64 vs x86_64

## Out of Scope

Explicitly excluded. Documented to prevent scope creep.

| Feature | Reason |
|---------|--------|
| Native macOS binaries | Docker is deployment target, not native installation |
| Windows ARM support | Linux containers only for now |
| Architecture-specific compose files | Single compose file must work everywhere |
| QEMU fallback documentation | Defeats purpose of native ARM64 support |

## Traceability

Which phases cover which requirements. Updated during roadmap creation.

| Requirement | Phase | Status |
|-------------|-------|--------|
| CONT-01 | Phase 42 | Complete |
| CONT-02 | Phase 42 | Complete |
| CONT-03 | Phase 42 | Complete |
| CONT-04 | Phase 42 | Complete |
| CONT-05 | Phase 41 | Complete |
| CONT-06 | Phase 42 | Complete |
| CICD-01 | Phase 43 | Pending |
| CICD-02 | Phase 43 | Pending |
| CICD-03 | Phase 43 | Pending |
| CICD-04 | Phase 43 | Pending |
| DOCS-01 | Phase 44 | Pending |
| DOCS-02 | Phase 44 | Pending |
| DOCS-03 | Phase 44 | Pending |

**Coverage:**
- v2.5 requirements: 13 total
- Mapped to phases: 13
- Unmapped: 0

---
*Requirements defined: 2026-02-04*
*Last updated: 2026-02-04 after roadmap creation*
