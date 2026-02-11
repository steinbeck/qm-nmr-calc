# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-10)

**Core value:** Reliable async NMR predictions with full control over calculation parameters -- submit a molecule, get back accurate 1H/13C shifts without babysitting long-running calculations.
**Current focus:** v2.9 Extended Solvent Coverage

## Current Position

Milestone: v2.9 Extended Solvent Coverage
Phase: Phase 61 - DELTA50 Toluene & DCM (Complete)
Plan: 61-01 complete (1/1)
Status: Phase 61 complete, ready for Phase 62
Last activity: 2026-02-11 — Completed 61-01-PLAN.md (100 benchmark verifications, 100% success)

Progress: [###########################] 123 plans complete across 12 milestones (v1.0-v2.8, v2.9-partial), 4 plans pending (v2.9)

## Performance Metrics

**Velocity:**
- Total plans completed: 123 (v1.0: 16, v1.1: 21, v2.0: 18, v2.0.1: 3, v2.1: 17, v2.2: 10, v2.3: 3, v2.4: 8, v2.5: 4, v2.6: 4, v2.7: 9, v2.8: 6, v2.9: 3)
- Average duration: ~6.5 min (excluding benchmark compute time)
- Total execution time: ~817 min (~13.6 hours) + ~27.5 hours benchmark compute

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
| v2.5 ARM64 Docker Support | 4 | 4 | ~1 day | Shipped 2026-02-04 |
| v2.6 GCP Spot Deployment | 5 | 4 | ~1 day | Shipped 2026-02-05 |
| v2.7 Automated GCP Deployment | 5 | 9 | ~32 min | Shipped 2026-02-06 |
| v2.8 Expanded Solvent Support | 5 | 6 | ~17h compute + 53min | Shipped 2026-02-09 |
| v2.9 Extended Solvent Coverage | 7 | 7 | In progress | Started 2026-02-10 (3/7 complete) |

## Accumulated Context

### Decisions

All decisions logged in PROJECT.md Key Decisions table.
v2.8 decisions archived to .planning/milestones/v2.8-ROADMAP.md.

Recent decisions (v2.9):

| ID | Phase | Decision | Impact |
|----|-------|----------|--------|
| COSMO-NAME-MAPPING | 59-01 | Map user-friendly solvent names to NWChem COSMO names | Acetonitrile accepted in CLI/code, mapped to "acetntrl" for NWChem COSMO |
| OPT-IN-NEW-SOLVENTS | 59-01 | New solvents opt-in only via --solvents flag | Users must explicitly request pyridine/thf/toluene/dcm/acetonitrile/dmf until scaling factors exist |
| PILOT-BEFORE-FULL | 60-01 | Run pilot (10 calc) before full benchmark (90 calc) | Validates COSMO convergence before committing ~10 hours compute |
| CHECKPOINT-VERIFICATION | 60-01 | Human verification at pilot and full completion checkpoints | Ensures calculations actually completed before marking phase done |

### Pending Todos

- Dark mode (deferred from v2.1)
- Enhanced interactivity (deferred from v2.1)
- User accounts and calculation history
- NMReData enhancements: INCHI tag, 1D pseudo-spectrum tags, J-coupling
- Phase 48.1 implementation (machine info display - v2.6)
- Additional solvents: pyridine-d5, THF-d8, toluene-d8

### Blockers/Concerns

**Active:**
None

## Session Continuity

Last session: 2026-02-11 06:51 UTC
Stopped at: Completed 61-01-PLAN.md
Resume file: None
Next: /gsd:plan-phase 62
Tests: 441 tests (382 passing, 2 skipped)
Codebase: ~7,320 LOC Python, ~3,115 LOC tests, ~950 LOC templates, ~2,400 LOC CSS, ~4,800 LOC docs, ~2,170 LOC GCP scripts
Benchmark data: 300 calculations (50 benzene × 4 solvents + 50 pyridine + 50 thf + 50 toluene + 50 dcm) with shielding tensors
