# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-01-29)

**Core value:** Reliable async NMR predictions with full control over calculation parameters -- submit a molecule, get back accurate 1H/13C shifts without babysitting long-running calculations.
**Current focus:** v2.1 UI Redesign -- Modern bento grid layout with glassmorphism effects

## Current Position

Milestone: v2.1 UI Redesign
Phase: 21 of 23 (Status Page Redesign) - COMPLETE
Plan: 2 of 2 in current phase
Status: Phase complete, verified
Last activity: 2026-01-31 -- Completed Phase 21 (Status Page Redesign)

Progress: ██████████████░░░░░░ 70% (4/6 v2.1 phases complete)

## Performance Metrics

**Velocity:**
- Total plans completed: 60 (v1.0: 16, v1.1: 21, v2.0: 18, v2.0.1: 3, v2.1: 2)
- Average duration: ~9 min
- Total execution time: ~538 min (~9 hours)

**By Milestone:**

| Milestone | Phases | Plans | Duration | Status |
|-----------|--------|-------|----------|--------|
| v1.0 Core NMR Service | 6 | 16 | 2 days | Shipped 2026-01-20 |
| v1.1 Accurate Chemical Shifts | 8 | 21 | 5 days | Shipped 2026-01-25 |
| v2.0 Conformational Sampling | 6 | 18 | ~2 days | Shipped 2026-01-28 |
| v2.0.1 Conformer Pre-selection | 1 | 3 | 18 min | Complete 2026-01-30 |
| v2.1 UI Redesign | 6 | 11/? | - | In Progress |

## Accumulated Context

### Decisions

All decisions logged in PROJECT.md Key Decisions table.

**v2.0.1 Conformer Pre-selection (complete 2026-01-30):**

- Target: ~8 conformers for DFT optimization (down from 40-50)
- Butina algorithm for RMSD clustering: Handles symmetry via GetBestRMS
- Default threshold 1.5 Angstrom: Aggressive reduction, balanced diversity
- Lowest-energy conformer per cluster: Best representative of conformational basin
- xTB as optional dependency: App works without it, uses MMFF fallback
- GFN2-xTB method: Balance of speed and accuracy for conformer ranking
- ALPB implicit solvation: Faster than full solvent models, adequate for ranking
- Subprocess with timeout: Prevent hanging on problematic structures
- Relative energies in kcal/mol: Consistent with MMFF output format
- Conditional clustering: Only cluster if conformer count exceeds target (avoid over-reduction)
- Graceful xTB fallback: Try xTB ranking first, catch exceptions and fall back to MMFF
- Metadata preservation: Track pre-clustering count in total_after_pre_filter field

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
- CSS load order: layers.css -> reset -> tokens -> base -> layout -> components -> utilities -> legacy
- Legacy styles in @layer components: Controlled migration path for existing custom.css
- Solid white background for viewer-card: WebGL/3Dmol.js has performance issues with glass/blur effects
- 44px min-height for download buttons: WCAG 2.5.5 touch target recommendation
- tabular-nums for shift tables: Ensures decimal points align in numeric columns
- Page-specific CSS in pages/ subdirectory: Separates page-level from generic components
- Jinja2 conditional rendering for ensemble elements: Cleaner than JS style.display toggles
- page_css block in base.html: Child templates inject page-specific CSS after components
- SmilesDrawer CDN via unpkg: Simple deployment without npm build complexity
- 400ms debounce delay: Balances responsiveness with preview performance
- devicePixelRatio canvas handling: Crisp rendering on retina displays

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
- v2.0.1: Hotfix for conformer pre-selection efficiency (complete 2026-01-30, 3 plans, 18 min)
- v2.1: 6 phases (18-23), paused at Phase 21
  - Phase structure: CSS Foundation -> Results -> Submit -> Status -> Responsive -> Accessibility
  - Foundation-first: Design system before page implementations
  - Complexity order: Results (most complex) -> Submit -> Status (simplest)

### Pending Todos

- Plan and execute Phase 22 (Responsive and Layout Polish)
- Plan and execute Phase 23 (Accessibility and Testing)
- Consider v2.0.1/v2.1 deployment and testing with production workloads

### Blockers/Concerns

**Active:**
None

**Resolved (v2.0.1):**
- Conformer pre-selection integration: Complete (plan 24-03)
  - Pipeline successfully integrates clustering and xTB ranking
  - Reduces DFT workload from 40-50 to ~8 conformers
  - Real-world test: hexanol reduced from 28 to 3 conformers
  - xTB ranking with graceful MMFF fallback
  - All tests passing (11 clustering + 14 xTB + 3 integration = 28 new tests)

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

Last session: 2026-01-31
Stopped at: Completed Phase 21 (Status Page Redesign)
Resume file: None
Next: Plan Phase 22 (Responsive and Layout Polish)
Tests: All tests passing (251 unit + 6 integration + 28 clustering/xTB = 285 tests)
Codebase: 5,797 LOC Python, 1,775 LOC tests, 941 LOC templates + 514+ lines CSS + SmilesDrawer integration
