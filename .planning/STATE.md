# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-10)

**Core value:** Reliable async NMR predictions with full control over calculation parameters -- submit a molecule, get back accurate 1H/13C shifts without babysitting long-running calculations.
**Current focus:** v2.9 Extended Solvent Coverage

## Current Position

Milestone: v2.9 Extended Solvent Coverage
Phase: Not started (defining requirements)
Plan: —
Status: Defining requirements
Last activity: 2026-02-10 — Milestone v2.9 started

Progress: [########################] 120 plans complete across 12 milestones (v1.0-v2.8)

## Performance Metrics

**Velocity:**
- Total plans completed: 120 (v1.0: 16, v1.1: 21, v2.0: 18, v2.0.1: 3, v2.1: 17, v2.2: 10, v2.3: 3, v2.4: 8, v2.5: 4, v2.6: 4, v2.7: 9, v2.8: 6)
- Average duration: ~6.5 min (excluding benchmark compute time)
- Total execution time: ~784 min (~13.1 hours) + ~17 hours benchmark compute

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

## Accumulated Context

### Decisions

All decisions logged in PROJECT.md Key Decisions table.
v2.8 decisions archived to .planning/milestones/v2.8-ROADMAP.md.

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

Last session: 2026-02-10
Stopped at: v2.8 milestone archived
Resume file: None
Next: /gsd:new-milestone to define next focus area
Tests: 434 tests (375 passing, 2 skipped)
Codebase: ~7,300 LOC Python, ~3,050 LOC tests, ~950 LOC templates, ~2,400 LOC CSS, ~4,800 LOC docs, ~2,170 LOC GCP scripts
