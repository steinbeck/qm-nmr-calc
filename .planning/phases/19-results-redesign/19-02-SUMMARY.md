---
phase: 19-results-redesign
plan: 02
subsystem: ui
tags: [jinja2, bento-grid, glass-card, viewer-card, templates]

# Dependency graph
requires:
  - phase: 19-01
    provides: results-page.css with viewer-card, shift-table, download-btn components
  - phase: 18
    provides: CSS architecture (layers, tokens, layout, glass-card component)
provides:
  - Bento grid results page template using Phase 18/19-01 CSS classes
  - Jinja2 conditional rendering for ensemble-only elements
  - Updated JavaScript with populateShiftTables() function
  - Page-specific CSS loading via page_css block
affects: [19-03, 20-submit-page, 21-status-page, 22-responsive, 23-accessibility]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Page-specific CSS via {% block page_css %} in child templates"
    - "Jinja2 conditionals for conformer_mode-based DOM rendering"
    - "Null checks in JavaScript for conditionally rendered DOM elements"

key-files:
  created: []
  modified:
    - src/qm_nmr_calc/api/templates/base.html
    - src/qm_nmr_calc/api/templates/results.html

key-decisions:
  - "Jinja2 conditional rendering instead of JavaScript style.display toggles for ensemble elements"
  - "Shift tables always visible (not in <details> accordion) for better accessibility"
  - "populateShiftTables() called for all jobs, not just ensemble"

patterns-established:
  - "page_css block: Child templates inject page-specific CSS after components, before utilities"
  - "setIfExists() helper: Safe DOM updates for conditionally rendered elements"
  - "bento-grid__item--span-N: Grid column span utility classes"

# Metrics
duration: 2min
completed: 2026-01-29
---

# Phase 19 Plan 02: Results Template Summary

**Bento grid layout for results.html with hero 3D viewer, glass cards for spectra/shifts/downloads, and Jinja2 conditional rendering for ensemble metadata**

## Performance

- **Duration:** 2 min
- **Started:** 2026-01-29T14:15:50Z
- **Completed:** 2026-01-29T14:17:43Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments
- Results page uses modern bento grid layout with 6-column desktop grid
- 3D viewer in prominent span-4 hero position spanning 2 rows
- Spectrum images, shift tables, and calculation details in span-2 glass cards
- Ensemble metadata card only renders for ensemble jobs (Jinja2 conditional, not JS toggle)
- Downloads in span-6 full-width card with styled download buttons
- JavaScript updated with null-safe DOM operations for conditional elements

## Task Commits

Each task was committed atomically:

1. **Task 1: Update base.html to load results-page.css conditionally** - `d9bab40` (feat)
2. **Task 2: Rewrite results.html with bento grid layout** - `f0141d7` (feat)

## Files Created/Modified
- `src/qm_nmr_calc/api/templates/base.html` - Added {% block page_css %} for page-specific CSS injection
- `src/qm_nmr_calc/api/templates/results.html` - Complete rewrite with bento grid layout

## Decisions Made
- **Jinja2 over JavaScript for conditional rendering:** Ensemble-only elements use `{% if job.conformer_mode == 'ensemble' %}` instead of `style="display: none"` with JS toggles. This is cleaner (no FOUC), more accessible (no hidden content), and reduces JavaScript complexity.
- **Shift tables always visible:** Removed `<details>` accordion wrapper from shift tables. Tables now display in dedicated cards for better discoverability.
- **Universal shift table population:** `populateShiftTables()` called for all jobs (not just ensemble), using existing `nmr_results` data from API.

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None - template rewrite and JavaScript updates worked as specified.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- Template uses all CSS classes from 19-01 (viewer-card, shift-table, download-btn)
- Ready for 19-03 manual visual verification checkpoint
- Same bento grid pattern can be applied to submit page (Phase 20) and status page (Phase 21)

---
*Phase: 19-results-redesign*
*Completed: 2026-01-29*
