# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-01-26)

**Core value:** Reliable async NMR predictions with full control over calculation parameters -- submit a molecule, get back accurate 1H/13C shifts without babysitting long-running calculations.
**Current focus:** v2.0 Conformational Sampling -- Phase 14 in progress

## Current Position

Phase: 14 of 17 (Boltzmann Averaging)
Plan: 01 of 01 (Phase 14)
Status: Phase complete
Last activity: 2026-01-27 -- Completed 14-01-PLAN.md (Boltzmann weight calculation)

Progress: ████████████████████████████████████████████ 100% (Phase 14: 1/1 plans complete)

## Performance Metrics

**Velocity:**
- Total plans completed: 44
- Average duration: 8.0 min
- Total execution time: 355 min (5.9 hours)

**By Milestone:**

| Milestone | Phases | Plans | Duration |
|-----------|--------|-------|----------|
| v1.0 Core NMR Service | 6 | 16 | 2 days |
| v1.1 Accurate Chemical Shifts | 8 | 21 | 5 days |
| v2.0 Conformational Sampling | 6 | TBD | In progress |

**Recent Trend:**
- Last 5 plans: 4-29 min
- Trend: Fast TDD (13-02: 4 min, 13-03: 29 min, 14-01: 6 min TDD with 20 tests)

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
- KDG params for solution-phase: pruneRmsThresh=-1.0 (no pre-opt pruning), numThreads=0 (all threads), random coords fallback (13-01)
- 8 rotatable bonds threshold: Separates rigid (50 confs) from flexible (200 confs) based on conformational space needs (13-01)
- MMFF native units: Return energies in kcal/mol without conversion (13-01)
- Dedup before energy filter: Run RMSD deduplication first, then energy window filter on deduped subset (13-03)
- Sequential 1-based conformer IDs: conf_001, conf_002, ... for user-facing consistency (13-03)
- Relative geometry_file paths: Stored relative to job dir for portability (13-03)
- exp-normalize trick for Boltzmann: Subtract min energy before exp to prevent overflow/underflow (14-01)
- Pure Python math.exp: No numpy dependency for Boltzmann weights, sufficient performance for typical conformer counts (14-01)

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
- CREST timeouts: Macrocycles can hang, need timeout with RDKit fallback (Phase 16)

**Resolved (v2.0):**
- Numerical stability: Boltzmann weighting implemented with exp-normalize trick, tested with extreme ranges (14-01)
- Atom ordering consistency: Canonical indexing established via CanonicalRankAtoms (12-02)
- Scratch directory isolation: Per-conformer directories prevent NWChem database corruption (12-01)
- Backward compatibility: v1.x JSON loads verified with 8 fixture tests (12-03)

**Resolved (v1.x):**
- RDKit stderr capture limitation (known, fallback messages used)
- Single-conformer limitation (being addressed in v2.0)

## Session Continuity

Last session: 2026-01-27
Stopped at: Completed 14-01-PLAN.md (Boltzmann weight calculation)
Resume file: None
Next: Phase 15 (Conformer NMR Integration) - Integrate Boltzmann weighting into NMR averaging
