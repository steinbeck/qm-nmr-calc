---
phase: 20-submit-redesign
plan: 02
subsystem: ui
tags: [html, jinja2, javascript, smiles-drawer, form-validation, molecule-preview]

# Dependency graph
requires:
  - phase: 20-01
    provides: submit-page.css layout and component styles
provides:
  - redesigned submit.html template with two-column layout
  - SmilesDrawer real-time molecule preview
  - improved form accessibility with aria attributes
  - debounced SMILES validation
affects: [21-status-page, future form enhancements]

# Tech tracking
tech-stack:
  added: [SmilesDrawer 1.2.0 via CDN]
  patterns: [debounced input validation, CDN script loading, retina canvas handling]

key-files:
  modified:
    - src/qm_nmr_calc/api/templates/submit.html

key-decisions:
  - "SmilesDrawer CDN via unpkg for simple deployment (no npm build)"
  - "400ms debounce delay balances responsiveness with performance"
  - "Canvas uses devicePixelRatio for retina display support"
  - "Input mutual exclusion preserved (SMILES disables file, file disables SMILES)"

patterns-established:
  - "CDN scripts in block scripts at end of template"
  - "Debounce pattern for real-time validation feedback"
  - "Status text with semantic classes for states (empty, valid, invalid)"

# Metrics
duration: ~5min
completed: 2026-01-29
---

# Phase 20 Plan 02: Submit Template and SmilesDrawer Integration Summary

**Two-column form layout with real-time SmilesDrawer molecule preview, required field indicators, and improved accessibility**

## Performance

- **Duration:** ~5 min (including checkpoint verification)
- **Started:** 2026-01-29T15:26:00Z
- **Completed:** 2026-01-29T15:31:00Z
- **Tasks:** 3 (2 auto + 1 checkpoint)
- **Files modified:** 1

## Accomplishments
- Restructured submit.html with two-column grid layout using submit-page.css
- Integrated SmilesDrawer 1.2.0 for real-time 2D molecule rendering
- Added required field indicators with accessibility attributes (aria-required)
- Implemented debounced SMILES validation with visual feedback
- Preserved all existing form functionality (input mutual exclusion, conformer mode/method interaction)

## Task Commits

Each task was committed atomically:

1. **Task 1: Update submit.html with two-column layout** - `f05bc44` (feat)
2. **Task 2: Add SmilesDrawer integration** - `f05bc44` (feat)
3. **Task 3: Visual verification checkpoint** - User approved

**Note:** Tasks 1 and 2 were committed together as a single atomic change.

## Files Created/Modified
- `src/qm_nmr_calc/api/templates/submit.html` - Redesigned with two-column layout, form groups, preview card, SmilesDrawer integration, and improved accessibility

## Decisions Made
- Used SmilesDrawer CDN (unpkg) instead of npm package for simpler deployment without build tools
- 400ms debounce delay chosen to balance typing responsiveness with preview performance
- devicePixelRatio handling ensures crisp molecule rendering on retina displays
- Status text uses semantic CSS classes for styling: --empty (gray), --valid (green), --invalid (red)
- File upload shows "no preview available" message since SmilesDrawer only handles SMILES strings

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- Submit page redesign complete (both CSS and template)
- Phase 20 complete - ready for Phase 21 (Status Page Redesign)
- SmilesDrawer pattern can be reused for any future SMILES input fields

---
*Phase: 20-submit-redesign*
*Completed: 2026-01-29*
