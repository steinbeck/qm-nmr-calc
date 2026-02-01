# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-01)

**Core value:** Reliable async NMR predictions with full control over calculation parameters -- submit a molecule, get back accurate 1H/13C shifts without babysitting long-running calculations.
**Current focus:** v2.3 NMReData Export - Enable machine-readable export of NMR prediction results in the NMReData standard format

## Current Position

Milestone: v2.3 NMReData Export
Phase: Phase 32 - Core NMReData Generation Module
Plan: 01 of 01 (complete)
Status: Phase 32 complete - ready for Phase 33
Last activity: 2026-02-01 — Completed 32-01-PLAN.md

Progress: █████░░░░░░░░░░░░░░░ 33% (Phase 32/34 complete, Phase 33/34 next)

## Performance Metrics

**Velocity:**
- Total plans completed: 86 (v1.0: 16, v1.1: 21, v2.0: 18, v2.0.1: 3, v2.1: 17, v2.2: 10, v2.3: 1)
- Average duration: ~7 min
- Total execution time: ~590 min (~9.8 hours)

**By Milestone:**

| Milestone | Phases | Plans | Duration | Status |
|-----------|--------|-------|----------|--------|
| v1.0 Core NMR Service | 6 | 16 | 2 days | Shipped 2026-01-20 |
| v1.1 Accurate Chemical Shifts | 8 | 21 | 5 days | Shipped 2026-01-25 |
| v2.0 Conformational Sampling | 6 | 18 | ~2 days | Shipped 2026-01-28 |
| v2.0.1 Conformer Pre-selection | 1 | 3 | 18 min | Shipped 2026-01-30 |
| v2.1 UI Redesign | 6 | 17 | ~3 days | Shipped 2026-01-31 |
| v2.2 Documentation | 7 | 10 | 2 days | Shipped 2026-02-01 |
| v2.3 NMReData Export | 3 | 1 (so far) | In progress | Active |

## Accumulated Context

### Decisions

All prior decisions logged in PROJECT.md Key Decisions table.
v2.2 decisions archived to milestones/v2.2-ROADMAP.md.

**Phase 32 decisions:**
- Lowercase atom labels (h5, c3) for NMReData format
- Preserve 1-based atom indices from NWChem (no conversion)
- Solvent mapping: chcl3→CDCl3, dmso→(CD3)2SO
- 4 decimal precision for shifts per NMReData spec

### Roadmap Evolution

- v1.0: 6 phases (1-6), shipped 2026-01-20
- v1.1: 8 phases (7-11.2, including 3 inserted), shipped 2026-01-25
- v2.0: 6 phases (12-17), shipped 2026-01-28
- v2.0.1: 1 phase (24), shipped 2026-01-30
- v2.1: 6 phases (18-23), shipped 2026-01-31
- v2.2: 7 phases (25-31), shipped 2026-02-01
- v2.3: 3 phases (32-34), in progress

### v2.3 Phase Structure

**Phase 32: Core NMReData Generation Module**
- Goal: Generate NMReData-compliant SDF files with predicted shifts
- Requirements: NMRD-01 to NMRD-09 (9 requirements)
- Key deliverables: nmredata.py module, tag formatting, atom numbering conversion
- Critical: Avoid off-by-one errors (RDKit 0-indexed → SDF 1-indexed)

**Phase 33: API and UI Integration**
- Goal: Enable download via REST endpoint and web UI button
- Requirements: API-01 to API-03, UI-01 (4 requirements)
- Key deliverables: GET /api/v1/jobs/{job_id}/nmredata.sdf endpoint, download button
- Follows existing download pattern from /geometry.sdf

**Phase 34: Testing and Validation**
- Goal: Comprehensive testing and round-trip validation
- Requirements: TEST-01 to TEST-03 (3 requirements)
- Key deliverables: Unit tests, integration tests, RDKit round-trip validation
- Validates format compliance and atom assignment correctness

### Pending Todos

- Deploy to production and test with real workloads
- Consider future enhancements: user accounts, calculation history, batch processing
- Dark mode (deferred from v2.1)
- Enhanced interactivity (deferred from v2.1)

### Blockers/Concerns

**Active:**
None

## Session Continuity

Last session: 2026-02-01
Stopped at: Completed Phase 32 Plan 01 (NMReData core module)
Resume file: None
Next: `/gsd:plan-phase 33` to begin Phase 33 (API and UI Integration)
Tests: All tests passing (285 prior + 36 nmredata = 321 tests)
Codebase: ~6,270 LOC Python, ~2,230 LOC tests, ~940 LOC templates, ~2,400 LOC CSS, ~4,100 LOC docs
