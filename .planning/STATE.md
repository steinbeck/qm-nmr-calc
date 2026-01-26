# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-01-26)

**Core value:** Reliable async NMR predictions with full control over calculation parameters -- submit a molecule, get back accurate 1H/13C shifts without babysitting long-running calculations.
**Current focus:** v2.0 Conformational Sampling -- Phase 12 ready to plan

## Current Position

Phase: 12 of 17 (Conformer Data Model and Storage)
Plan: 02 of 03 (Phase 12)
Status: In progress
Last activity: 2026-01-26 -- Completed 12-02-PLAN.md (Canonical Atom Ordering)

Progress: █░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░ 33% (Phase 12: 1/3 plans)

## Performance Metrics

**Velocity:**
- Total plans completed: 37
- Average duration: 8.0 min
- Total execution time: 302 min (5.0 hours)

**By Milestone:**

| Milestone | Phases | Plans | Duration |
|-----------|--------|-------|----------|
| v1.0 Core NMR Service | 6 | 16 | 2 days |
| v1.1 Accurate Chemical Shifts | 8 | 21 | 5 days |
| v2.0 Conformational Sampling | 6 | TBD | In progress |

**Recent Trend:**
- Last 5 plans: 3-12 min
- Trend: Fast (12-02: TDD module, 3 min)

## Accumulated Context

### Decisions

All decisions logged in PROJECT.md Key Decisions table.
Recent v2.0 decisions affecting current work:

- KDG over ETKDG: Avoids crystal structure bias for solution-phase NMR
- CREST/xTB optional: App works without them, enables high-accuracy mode when detected
- Boltzmann weight by DFT energies: Most accurate readily available energy level for weighting
- RDKit CanonicalRankAtoms for atom ordering: Built-in deterministic canonicalization (12-02)
- Mapping as metadata not transformation: Don't physically reorder atoms (breaks RDKit ops) (12-02)

### Roadmap Evolution

- v1.0: 6 phases (1-6), shipped 2026-01-20
- v1.1: 8 phases (7-11.2, including 3 inserted), shipped 2026-01-25
- v2.0: 6 phases (12-17), started 2026-01-26
  - Phase structure: Data Model → RDKit → Boltzmann → NWChem → CREST → API
  - Risk mitigation: RDKit-only path (12-15, 17) complete before CREST complexity (16)

### Pending Todos

None.

### Blockers/Concerns

**v2.0 Architecture:**
- Numerical stability: Boltzmann weighting with wide energy ranges needs exp-normalize trick (Phase 14)
- Scratch directory isolation: Per-conformer directories critical to prevent NWChem database corruption (Phase 12)
- CREST timeouts: Macrocycles can hang, need timeout with RDKit fallback (Phase 16)

**Resolved (v2.0):**
- Atom ordering consistency: Canonical indexing established via CanonicalRankAtoms (12-02)

**Resolved (v1.x):**
- RDKit stderr capture limitation (known, fallback messages used)
- Single-conformer limitation (being addressed in v2.0)

## Session Continuity

Last session: 2026-01-26
Stopped at: Completed 12-02-PLAN.md (Canonical Atom Ordering via TDD)
Resume file: None
Next: Continue Phase 12 (2 of 3 plans remaining)
