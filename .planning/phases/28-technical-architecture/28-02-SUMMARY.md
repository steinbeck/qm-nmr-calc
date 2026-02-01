---
phase: 28-technical-architecture
plan: 02
subsystem: documentation
tags: [mermaid, conformers, css, bem, cascade-layers, tokens]

# Dependency graph
requires:
  - phase: 28-01
    provides: Core architecture documentation (tech stack, data flow, job lifecycle, file storage)
provides:
  - Conformer ensemble pipeline documentation with 9-stage diagram
  - Conformer state machine diagram
  - RDKit vs CREST comparison
  - CSS cascade layers architecture
  - Design tokens reference
  - BEM naming convention guide
  - File organization diagram
affects: [contributor-onboarding, css-development, conformer-modifications]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Mermaid flowchart for pipeline visualization"
    - "Mermaid stateDiagram-v2 for conformer states"
    - "CSS layer documentation with @layer syntax examples"

key-files:
  created:
    - ".planning/phases/28-technical-architecture/28-02-SUMMARY.md"
  modified:
    - "docs/architecture.md"

key-decisions:
  - "Documented actual file names from codebase (generator.py, not rdkit_generator.py)"
  - "Included xTB pre-ranking as optional enhancement"
  - "Showed complete token categories from actual tokens.css"

patterns-established:
  - "Architecture docs use Mermaid for all diagrams"
  - "Pipeline stages shown as flowchart with method annotations"
  - "CSS docs include example component code"

# Metrics
duration: 2min
completed: 2026-02-01
---

# Phase 28 Plan 02: Pipeline and CSS Architecture Summary

**9-stage conformer ensemble pipeline with state machine plus CSS cascade layers, design tokens, and BEM conventions**

## Performance

- **Duration:** 2 min
- **Started:** 2026-02-01T08:11:48Z
- **Completed:** 2026-02-01T08:14:02Z
- **Tasks:** 3
- **Files modified:** 1

## Accomplishments
- Documented complete 9-stage conformer pipeline from SMILES to Boltzmann averaging
- Added conformer state machine diagram showing pending -> optimizing -> nmr_complete flow
- Documented RDKit ETKDGv3 vs CREST comparison with use case guidance
- Documented CSS cascade layer system with priority order
- Included design token categories with actual values from tokens.css
- Created BEM naming guide with component examples
- Documented CSS file organization with components/ and pages/ directories

## Task Commits

Each task was committed atomically:

1. **Task 1: Conformer Ensemble Pipeline Documentation** - Combined with Task 2
2. **Task 2: CSS Architecture Documentation** - Combined with Task 1
3. **Task 3: Commit Pipeline and CSS Architecture Documentation** - `ea0dc8e`

**Plan metadata:** (pending)

## Files Created/Modified
- `docs/architecture.md` - Added Conformer Ensemble Pipeline (~120 lines) and CSS Architecture (~180 lines) sections

## Decisions Made
- Used actual file names from codebase inspection (generator.py not rdkit_generator.py)
- Documented xTB pre-ranking as optional enhancement in pipeline
- Included complete token categories matching actual tokens.css content
- Documented both components/ and pages/ directories per actual file structure

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
None

## User Setup Required

None - documentation only, no external service configuration required.

## Next Phase Readiness
- Phase 28 technical architecture documentation complete
- All 6 sections documented: Technology Stack, Data Flow, Job Lifecycle, File Storage, Conformer Pipeline, CSS Architecture
- docs/architecture.md at 717 lines (exceeds 350 minimum)
- Ready for Phase 29 (Contributing Guide)

---
*Phase: 28-technical-architecture*
*Completed: 2026-02-01*
