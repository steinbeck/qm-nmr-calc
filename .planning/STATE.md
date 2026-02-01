# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-01)

**Core value:** Reliable async NMR predictions with full control over calculation parameters -- submit a molecule, get back accurate 1H/13C shifts without babysitting long-running calculations.
**Current focus:** v2.3 NMReData Export - Enable machine-readable export of NMR prediction results in the NMReData standard format

## Current Position

Milestone: v2.3 NMReData Export
Phase: Phase 32 - Core NMReData Generation Module
Plan: —
Status: Roadmap created, ready to begin Phase 32
Last activity: 2026-02-01 — Milestone v2.3 roadmap created

Progress: ░░░░░░░░░░░░░░░░░░░░ 0% (Phase 32/34)

## Performance Metrics

**Velocity:**
- Total plans completed: 85 (v1.0: 16, v1.1: 21, v2.0: 18, v2.0.1: 3, v2.1: 17, v2.2: 10)
- Average duration: ~7 min
- Total execution time: ~580 min (~9.7 hours)

**By Milestone:**

| Milestone | Phases | Plans | Duration | Status |
|-----------|--------|-------|----------|--------|
| v1.0 Core NMR Service | 6 | 16 | 2 days | Shipped 2026-01-20 |
| v1.1 Accurate Chemical Shifts | 8 | 21 | 5 days | Shipped 2026-01-25 |
| v2.0 Conformational Sampling | 6 | 18 | ~2 days | Shipped 2026-01-28 |
| v2.0.1 Conformer Pre-selection | 1 | 3 | 18 min | Shipped 2026-01-30 |
| v2.1 UI Redesign | 6 | 17 | ~3 days | Shipped 2026-01-31 |
| v2.2 Documentation | 7 | 10 | 2 days | Shipped 2026-02-01 |
| v2.3 NMReData Export | 3 | TBD | In progress | Active |

## Accumulated Context

### Decisions

All prior decisions logged in PROJECT.md Key Decisions table.
v2.2 decisions archived to milestones/v2.2-ROADMAP.md.

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
Stopped at: v2.3 roadmap created
Resume file: None
Next: `/gsd:plan-phase 32` to begin Phase 32 implementation
Tests: All tests passing (257 unit + 28 conformer/xTB = 285 tests)
Codebase: ~6,000 LOC Python, ~1,800 LOC tests, ~940 LOC templates, ~2,400 LOC CSS, ~4,100 LOC docs
