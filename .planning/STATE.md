# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-01-29)

**Core value:** Reliable async NMR predictions with full control over calculation parameters -- submit a molecule, get back accurate 1H/13C shifts without babysitting long-running calculations.
**Current focus:** v2.1 UI Redesign -- Modern bento grid layout with glassmorphism effects

## Current Position

Milestone: v2.1 UI Redesign
Phase: 18 - CSS Foundation and Design System
Plan: Not started (roadmap created)
Status: Ready to begin Phase 18
Last activity: 2026-01-29 -- Created v2.1 roadmap (6 phases, 35 requirements)

Progress: ▓░░░░░ Phase 18 of 23 (22%)

## Performance Metrics

**Velocity:**
- Total plans completed: 55 (all from v1.0, v1.1, v2.0)
- Average duration: ~9 min
- Total execution time: ~500 min (~8.3 hours)

**By Milestone:**

| Milestone | Phases | Plans | Duration | Status |
|-----------|--------|-------|----------|--------|
| v1.0 Core NMR Service | 6 | 16 | 2 days | Shipped 2026-01-20 |
| v1.1 Accurate Chemical Shifts | 8 | 21 | 5 days | Shipped 2026-01-25 |
| v2.0 Conformational Sampling | 6 | 18 | ~2 days | Shipped 2026-01-28 |
| v2.1 UI Redesign | 6 | 0 | - | In progress |

**Recent Trend:**
- Last 5 plans: 9-20 min (v2.0 Phase 17)
- Trend: Consistent fast execution

## Accumulated Context

### Decisions

All decisions logged in PROJECT.md Key Decisions table.

**v2.1 UI Redesign decisions:**

- Pure CSS approach: No build tools, frameworks, or preprocessors needed
- Replace Pico CSS entirely: Custom CSS simpler than overriding framework for bento grids
- CSS Grid with grid-template-areas: Semantic bento layouts with named regions
- Native backdrop-filter: 92% browser support, 96% with fallbacks (webkit prefix)
- CSS Cascade Layers: Architectural organization without specificity wars
- BEM naming convention: Component boundaries prevent naming conflicts
- Multi-file CSS: HTTP/2 makes multiple files performant without concatenation
- Design tokens as CSS custom properties: Replaces preprocessor variables
- Glass opacity 85-95%: Higher opacity for WCAG contrast compliance
- Mobile performance budget: 2-3 glass elements max, 6px blur vs 10px desktop

**v2.0 Conformational Sampling decisions (inherited context):**

- KDG over ETKDG: Avoids crystal structure bias for solution-phase NMR
- CREST/xTB optional: App works without them, enables high-accuracy mode when detected
- Boltzmann weight by DFT energies: Most accurate readily available energy level for weighting
- RDKit CanonicalRankAtoms for atom ordering: Built-in deterministic canonicalization
- Per-conformer scratch directory isolation: Each conformer gets unique scratch dir to prevent NWChem database conflicts
- Explicit energy unit tracking: ConformerData stores energy_unit alongside energy to prevent conversion bugs
- exp-normalize trick for Boltzmann: Subtract min energy before exp to prevent overflow/underflow
- Sequential conformer processing: Process conformers one at a time for simplicity
- 50% success threshold: Fail if more than half of conformers fail DFT optimization

### Roadmap Evolution

- v1.0: 6 phases (1-6), shipped 2026-01-20
- v1.1: 8 phases (7-11.2, including 3 inserted), shipped 2026-01-25
- v2.0: 6 phases (12-17), shipped 2026-01-28
  - Phase structure: Data Model -> RDKit -> Boltzmann -> NWChem -> CREST -> API
  - Risk mitigation: RDKit-only path (12-15, 17) complete before CREST complexity (16)
- v2.1: 6 phases (18-23), in progress
  - Phase structure: CSS Foundation -> Results -> Submit -> Status -> Responsive -> Accessibility
  - Foundation-first: Design system before page implementations
  - Complexity order: Results (most complex) -> Submit -> Status (simplest)

### Pending Todos

None. Ready to start Phase 18.

### Blockers/Concerns

**None.**

**Resolved (v2.0):**
- Ensemble filtering before averaging: run_ensemble_nmr_task now filters to nmr_complete conformers before calling average_ensemble_nmr
- CREST timeouts: Implemented subprocess timeout with user-friendly fallback message
- Numerical stability: Boltzmann weighting implemented with exp-normalize trick, tested with extreme ranges
- Atom ordering consistency: Canonical indexing established via CanonicalRankAtoms
- Scratch directory isolation: Per-conformer directories prevent NWChem database corruption
- Backward compatibility: v1.x JSON loads verified with 8 fixture tests
- DFT energy regex mismatch: Fixed to accept both `=` and `:` in "Total DFT energy" line
- Functional case sensitivity: Added `.upper()` normalization in `get_scaling_factor`
- Shared optimized geometry path: Fixed by copying to per-conformer output dir after each optimization

**Resolved (v1.x):**
- RDKit stderr capture limitation (known, fallback messages used)
- Single-conformer limitation (addressed in v2.0)

## Session Continuity

Last session: 2026-01-29
Stopped at: v2.1 roadmap created, ready to plan Phase 18
Resume file: None
Next: `/gsd:plan-phase 18` to begin CSS Foundation and Design System
Tests: All tests passing (226 unit + 3 integration = 229 tests)
Codebase: 5,432 LOC Python, 1,417 LOC tests, 892 LOC templates
