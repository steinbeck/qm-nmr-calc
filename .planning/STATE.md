# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-01-26)

**Core value:** Reliable async NMR predictions with full control over calculation parameters -- submit a molecule, get back accurate 1H/13C shifts without babysitting long-running calculations.
**Current focus:** v2.0 Conformational Sampling -- Phase 12 complete, Phase 13 next

## Current Position

Phase: 13 of 17 (RDKit Conformer Generation)
Plan: 02 of 03 (Phase 13)
Status: In progress
Last activity: 2026-01-27 -- Completed 13-02-PLAN.md (Conformer filters)

Progress: ███░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░ 33% (Phase 13: 1/3 plans)

## Performance Metrics

**Velocity:**
- Total plans completed: 40
- Average duration: 7.7 min
- Total execution time: 316 min (5.3 hours)

**By Milestone:**

| Milestone | Phases | Plans | Duration |
|-----------|--------|-------|----------|
| v1.0 Core NMR Service | 6 | 16 | 2 days |
| v1.1 Accurate Chemical Shifts | 8 | 21 | 5 days |
| v2.0 Conformational Sampling | 6 | TBD | In progress |

**Recent Trend:**
- Last 5 plans: 3-6 min
- Trend: Fast (12-02: TDD module, 3 min; 12-03: Compat + API, 5.6 min; 13-02: Filters TDD, 4 min)

## Accumulated Context

### Decisions

All decisions logged in PROJECT.md Key Decisions table.
Recent v2.0 decisions affecting current work:

- KDG over ETKDG: Avoids crystal structure bias for solution-phase NMR
- CREST/xTB optional: App works without them, enables high-accuracy mode when detected
- Boltzmann weight by DFT energies: Most accurate readily available energy level for weighting
- RDKit CanonicalRankAtoms for atom ordering: Built-in deterministic canonicalization (12-02)
- Mapping as metadata not transformation: Don't physically reorder atoms (breaks RDKit ops) (12-02)
- Per-conformer scratch directory isolation: Each conformer gets unique scratch dir to prevent NWChem database conflicts (12-01)
- Explicit energy unit tracking: ConformerData stores energy_unit alongside energy to prevent conversion bugs (12-01)
- str type for API conformer_mode: API schemas use str (not Literal) for flexibility, validation at service layer (12-03)
- Pairwise GetBestRMS over matrix approach: Clearer code, acceptable performance for typical ensemble sizes (13-02)
- Decoupled energy filtering: Operates on extracted data, reusable for MMFF and DFT (13-02)

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
- CREST timeouts: Macrocycles can hang, need timeout with RDKit fallback (Phase 16)

**Resolved (v2.0):**
- Atom ordering consistency: Canonical indexing established via CanonicalRankAtoms (12-02)
- Scratch directory isolation: Per-conformer directories prevent NWChem database corruption (12-01)
- Backward compatibility: v1.x JSON loads verified with 8 fixture tests (12-03)

**Resolved (v1.x):**
- RDKit stderr capture limitation (known, fallback messages used)
- Single-conformer limitation (being addressed in v2.0)

## Session Continuity

Last session: 2026-01-27
Stopped at: Completed 13-02-PLAN.md (Conformer filters)
Resume file: None
Next: Continue Phase 13 (13-03: Pipeline integration)
