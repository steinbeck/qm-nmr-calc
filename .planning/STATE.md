# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-07)

**Core value:** Reliable async NMR predictions with full control over calculation parameters -- submit a molecule, get back accurate 1H/13C shifts without babysitting long-running calculations.
**Current focus:** v2.8 Expanded Solvent Support

## Current Position

Milestone: v2.8 Expanded Solvent Support
Phase: Not started (defining requirements)
Status: Defining requirements
Last activity: 2026-02-07 — Milestone v2.8 started

Progress: [#####################] 113 plans complete across 11 milestones (v1.0-v2.7)

## Performance Metrics

**Velocity:**
- Total plans completed: 113 (v1.0: 16, v1.1: 21, v2.0: 18, v2.0.1: 3, v2.1: 17, v2.2: 10, v2.3: 3, v2.4: 8, v2.5: 4, v2.6: 4, v2.7: 9)
- Average duration: ~6.3 min
- Total execution time: ~719 min (~12.0 hours)

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

## Accumulated Context

### Decisions

All prior decisions logged in PROJECT.md Key Decisions table.
v2.7 decisions archived to .planning/milestones/v2.7-ROADMAP.md.

### Roadmap Evolution

- v1.0: 6 phases (1-6), shipped 2026-01-20
- v1.1: 8 phases (7-11.2, including 3 inserted), shipped 2026-01-25
- v2.0: 6 phases (12-17), shipped 2026-01-28
- v2.0.1: 1 phase (24), shipped 2026-01-30
- v2.1: 6 phases (18-23), shipped 2026-01-31
- v2.2: 7 phases (25-31), shipped 2026-02-01
- v2.3: 3 phases (32-34), shipped 2026-02-01
- v2.4: 6 phases (35-40), shipped 2026-02-03
- v2.5: 4 phases (41-44), shipped 2026-02-04
- v2.6: 5 phases (45-48.1), shipped 2026-02-05
- v2.7: 5 phases (49-53), shipped 2026-02-06

### Pending Todos

- Dark mode (deferred from v2.1)
- Enhanced interactivity (deferred from v2.1)
- User accounts and calculation history
- NMReData enhancements: INCHI tag, 1D pseudo-spectrum tags, J-coupling
- Phase 48.1 implementation (machine info display - v2.6)

### Blockers/Concerns

**Active:**
None

## Session Continuity

Last session: 2026-02-06
Stopped at: v2.7 milestone completed and archived
Resume file: None
Next: No active milestone. Run /gsd:new-milestone to start next version.
Tests: 415 tests (377 pre-existing + 19 config + 19 pricing + 19 machine)
Codebase: ~7,300 LOC Python, ~3,050 LOC tests, ~950 LOC templates, ~2,400 LOC CSS, ~4,800 LOC docs, ~2,170 LOC GCP scripts
Docker: Worker image 2.1GB (amd64), API image ~733MB (multi-arch), ARM64 worker 2.1GB (arm64), multi-arch manifests on GHCR
GCP: v2.7 shipped - TOML config, dynamic pricing/machine selection, HTTP-only deployment
