---
phase: 05-visualization
plan: 01
subsystem: visualization
tags: [matplotlib, rdkit, svg, png, nmr-spectrum]

# Dependency graph
requires:
  - phase: 04-results-delivery
    provides: NMRResults model with h1_shifts and c13_shifts
provides:
  - generate_spectrum_plot() for 1H/13C stick spectra (SVG/PNG)
  - generate_annotated_structure() for molecule images with shift labels
affects: [05-02 (integration), 06-web-ui]

# Tech tracking
tech-stack:
  added: [matplotlib>=3.10.0]
  patterns: [Agg backend for server rendering, atomNote for RDKit annotations]

key-files:
  created:
    - src/qm_nmr_calc/visualization.py
  modified:
    - pyproject.toml

key-decisions:
  - "Agg backend set before pyplot import for headless rendering"
  - "NWChem 1-based to RDKit 0-based index conversion"
  - "2400x1800 PNG for 300 DPI equivalent at 8x6 inches"
  - "atomNote set before PrepareMolForDrawing()"

patterns-established:
  - "plt.close(fig) after savefig to prevent memory leaks"
  - "stem() with empty markerfmt/basefmt for clean stick plots"

# Metrics
duration: 2min
completed: 2026-01-20
---

# Phase 5 Plan 1: Visualization Module Summary

**Matplotlib stem plots for NMR spectra and RDKit atomNote for annotated molecular structures, outputting SVG and PNG formats**

## Performance

- **Duration:** 2 min
- **Started:** 2026-01-20T10:27:26Z
- **Completed:** 2026-01-20T10:29:17Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments
- Added matplotlib 3.10.8 dependency for spectrum plotting
- Created generate_spectrum_plot() with inverted x-axis (NMR convention)
- Created generate_annotated_structure() using RDKit atomNote property
- Both functions output SVG and PNG (300 DPI) formats

## Task Commits

Each task was committed atomically:

1. **Task 1: Add matplotlib dependency** - `b224bc8` (chore)
2. **Task 2: Create visualization module** - `224562a` (feat)

## Files Created/Modified
- `pyproject.toml` - Added matplotlib>=3.10.0 to dependencies
- `src/qm_nmr_calc/visualization.py` - New visualization module with two main functions

## Decisions Made
- Use Agg backend (must be set before pyplot import) for server-side rendering
- Convert NWChem 1-based indices to RDKit 0-based when mapping shifts to atoms
- Set atomNote BEFORE calling PrepareMolForDrawing() per RDKit documentation
- Use 2400x1800 pixels for PNG to achieve 300 DPI equivalent at 8x6 inches
- Use stem() with empty markerfmt=' ' and basefmt=' ' for clean stick plots

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- Visualization functions ready for integration with tasks.py
- API endpoints can be added to serve generated images
- Functions tested with ethanol (CCO) structure successfully

---
*Phase: 05-visualization*
*Completed: 2026-01-20*
