# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-01-31)

**Core value:** Reliable async NMR predictions with full control over calculation parameters -- submit a molecule, get back accurate 1H/13C shifts without babysitting long-running calculations.
**Current focus:** v2.2 Documentation milestone

## Current Position

Milestone: v2.2 Documentation
Phase: 30 of 7 (DP4+ Science Documentation)
Plan: 1 of 2 in current phase
Status: In progress
Last activity: 2026-02-01 -- Completed 30-01-PLAN.md (NMR fundamentals, DFT theory, COSMO)

Progress: ████████████████░░░░ 78% (6.5/7 phases complete)

## Performance Metrics

**Velocity:**
- Total plans completed: 83 (v1.0: 16, v1.1: 21, v2.0: 18, v2.0.1: 3, v2.1: 17, v2.2: 8)
- Average duration: ~7 min
- Total execution time: ~573 min (~9.5 hours)

**By Milestone:**

| Milestone | Phases | Plans | Duration | Status |
|-----------|--------|-------|----------|--------|
| v1.0 Core NMR Service | 6 | 16 | 2 days | Shipped 2026-01-20 |
| v1.1 Accurate Chemical Shifts | 8 | 21 | 5 days | Shipped 2026-01-25 |
| v2.0 Conformational Sampling | 6 | 18 | ~2 days | Shipped 2026-01-28 |
| v2.0.1 Conformer Pre-selection | 1 | 3 | 18 min | Complete 2026-01-30 |
| v2.1 UI Redesign | 6 | 17 | ~8 min avg | Complete 2026-01-31 |
| v2.2 Documentation | 7 | 8+ | - | In Progress |

## Accumulated Context

### Decisions

All prior decisions logged in PROJECT.md Key Decisions table.

**v2.2 Documentation decisions:**

- Audience: Both academic researchers and developers/contributors
- Science depth: Full derivation of DP4+ methodology with literature references
- Version: v2.2 (next minor after v2.1)
- Research approach: Leverage planning docs + web research for DP4+ content
- Documentation structure: docs/ directory with README linking to detailed pages
- Installation guide: Multi-distro instructions, conda + manual CREST/xTB paths
- Validation approach: Health endpoint as primary environment check
- Troubleshooting: Problem/Symptom/Fix format for common issues
- Usage guide: Web UI first, REST API in separate plan (27-02)
- API docs: 26 curl examples covering all endpoints with complete workflow script
- Architecture docs: 7 Mermaid diagrams (stack, data flow, job states, conformer states, pipeline stages)
- Conformer pipeline: 9 stages documented with state machine and RDKit/CREST comparison
- CSS architecture: Cascade layers, design tokens, BEM naming, file organization
- Library docs: Code-first documentation with source file references, integration points tables
- Library docs: Index conversion warning (3Dmol.js 0-based vs NWChem 1-based)
- Science docs: MathJax equations, DOI citations, code-to-theory links

### Roadmap Evolution

- v1.0: 6 phases (1-6), shipped 2026-01-20
- v1.1: 8 phases (7-11.2, including 3 inserted), shipped 2026-01-25
- v2.0: 6 phases (12-17), shipped 2026-01-28
- v2.0.1: Hotfix for conformer pre-selection (complete 2026-01-30)
- v2.1: 6 phases (18-23), shipped 2026-01-31
- v2.2: 7 phases (25-31), documentation milestone

### Pending Todos

- Deploy to production and test with real workloads
- Consider future enhancements: user accounts, calculation history, batch processing
- Dark mode (deferred from v2.1)
- Enhanced interactivity (deferred from v2.1)

### Recently Fixed

- UX: Conformer progress bar shows 0% until first conformer completes (fixed in bf27668)

### Blockers/Concerns

**Active:**
None

## Session Continuity

Last session: 2026-02-01
Stopped at: Completed 30-01-PLAN.md (NMR fundamentals, DFT theory, COSMO)
Resume file: None
Next: Phase 30 Plan 02 (Linear scaling, Boltzmann, conformers, accuracy)
Tests: All tests passing (257 unit + 28 conformer/xTB = 285 tests)
Codebase: ~6,000 LOC Python, ~1,800 LOC tests, ~940 LOC templates, ~2,400 LOC CSS, ~3,200 LOC docs
