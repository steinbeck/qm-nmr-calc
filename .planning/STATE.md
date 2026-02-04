# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-03)

**Core value:** Reliable async NMR predictions with full control over calculation parameters -- submit a molecule, get back accurate 1H/13C shifts without babysitting long-running calculations.
**Current focus:** v2.5 ARM64 Docker Support - Phase 41 (ARM64 Dockerfile Creation)

## Current Position

Milestone: v2.5 ARM64 Docker Support
Phase: 41 of 44 (ARM64 Dockerfile Creation)
Plan: Ready to plan
Status: Ready to plan
Last activity: 2026-02-04 -- v2.5 roadmap created

Progress: [####################] 96 plans complete (v1.0-v2.4)
         [                    ] 0% of v2.5 (0/4 plans)

## Performance Metrics

**Velocity:**
- Total plans completed: 96 (v1.0: 16, v1.1: 21, v2.0: 18, v2.0.1: 3, v2.1: 17, v2.2: 10, v2.3: 3, v2.4: 8)
- Average duration: ~7 min
- Total execution time: ~663 min (~11 hours)

**By Milestone:**

| Milestone | Phases | Plans | Duration | Status |
|-----------|--------|-------|----------|--------|
| v1.0 Core NMR Service | 6 | 16 | 2 days | Shipped 2026-01-20 |
| v1.1 Accurate Chemical Shifts | 8 | 21 | 5 days | Shipped 2026-01-25 |
| v2.0 Conformational Sampling | 6 | 18 | ~2 days | Shipped 2026-01-28 |
| v2.0.1 Conformer Pre-selection | 1 | 3 | 18 min | Shipped 2026-01-30 |
| v2.1 UI Redesign | 6 | 17 | ~3 days | Shipped 2026-01-31 |
| v2.2 Documentation | 7 | 10 | 2 days | Shipped 2026-02-01 |
| v2.3 NMReData Export | 3 | 3 | 1 day | Shipped 2026-02-01 |
| v2.4 Docker Deployment | 6 | 8 | ~2 hours | Shipped 2026-02-03 |
| v2.5 ARM64 Docker Support | 4 | 4 | - | In progress |

## Accumulated Context

### Decisions

All prior decisions logged in PROJECT.md Key Decisions table.
v2.4 decisions archived to MILESTONES.md.

**v2.5 Decisions (pending):**
- Separate Dockerfile.worker.arm64 (not unified Dockerfile with build args)
- conda-forge packages for NWChem, xTB, CREST on ARM64 (no pre-compiled binaries available)
- micromamba base image for faster package resolution
- GitHub Actions ubuntu-24.04-arm runner for native ARM64 builds

### Roadmap Evolution

- v1.0: 6 phases (1-6), shipped 2026-01-20
- v1.1: 8 phases (7-11.2, including 3 inserted), shipped 2026-01-25
- v2.0: 6 phases (12-17), shipped 2026-01-28
- v2.0.1: 1 phase (24), shipped 2026-01-30
- v2.1: 6 phases (18-23), shipped 2026-01-31
- v2.2: 7 phases (25-31), shipped 2026-02-01
- v2.3: 3 phases (32-34), shipped 2026-02-01
- v2.4: 6 phases (35-40), shipped 2026-02-03
- v2.5: 4 phases (41-44), in progress

### Pending Todos

- Dark mode (deferred from v2.1)
- Enhanced interactivity (deferred from v2.1)
- User accounts and calculation history
- NMReData enhancements: INCHI tag, 1D pseudo-spectrum tags, J-coupling

### Blockers/Concerns

**Active:**
None

**Research Flags (from v2.5 research):**
- Phase 42: NWChem basis set paths may differ in conda-forge package (verify NWCHEM_BASIS_LIBRARY)
- Phase 42: OpenMPI configuration needs OMPI_ALLOW_RUN_AS_ROOT
- Phase 42: OpenBLAS threading needs explicit OPENBLAS_NUM_THREADS=4

## Session Continuity

Last session: 2026-02-04
Stopped at: v2.5 roadmap created, ready for plan-phase 41
Resume file: None
Next: `/gsd:plan-phase 41` to create ARM64 Dockerfile plan
Tests: All tests passing (356 tests)
Codebase: ~6,400 LOC Python, ~2,450 LOC tests, ~950 LOC templates, ~2,400 LOC CSS, ~4,560 LOC docs
Docker: Worker image 2.1GB (amd64), API image ~733MB (multi-arch), Caddy reverse proxy
