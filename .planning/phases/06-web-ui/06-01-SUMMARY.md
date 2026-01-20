---
phase: 06-web-ui
plan: 01
subsystem: ui
tags: [jinja2, pico-css, fastapi, templates, static-files]

# Dependency graph
requires:
  - phase: 02-input-and-api
    provides: FastAPI app and routers structure
provides:
  - Jinja2 template infrastructure
  - Base HTML template with Pico CSS blue theme
  - Static file serving at /static
  - Web router for browser routes
affects: [06-02, 06-03, 06-04]

# Tech tracking
tech-stack:
  added: [pico-css-2, jinja2-templates]
  patterns: [template-inheritance, static-file-serving]

key-files:
  created:
    - src/qm_nmr_calc/api/templates/base.html
    - src/qm_nmr_calc/api/templates/submit.html
    - src/qm_nmr_calc/api/static/css/custom.css
    - src/qm_nmr_calc/api/routers/web.py
  modified:
    - src/qm_nmr_calc/api/app.py
    - src/qm_nmr_calc/api/routers/__init__.py

key-decisions:
  - "Pico CSS blue theme via CDN for minimal/clean scientific aesthetic"
  - "Static files mounted at /static before routers to avoid route conflicts"
  - "Web router at root (no /api prefix) for browser-friendly URLs"

patterns-established:
  - "Template inheritance: all pages extend base.html"
  - "url_for('static', path=...) for static asset references"
  - "Template context includes solvents and presets from domain modules"

# Metrics
duration: 2min
completed: 2026-01-20
---

# Phase 6 Plan 1: Template Infrastructure Summary

**Jinja2 template system with Pico CSS blue theme and static file serving for web UI foundation**

## Performance

- **Duration:** 2 min
- **Started:** 2026-01-20T12:09:20Z
- **Completed:** 2026-01-20T12:11:43Z
- **Tasks:** 3
- **Files modified:** 6

## Accomplishments

- Base HTML template with Pico CSS blue theme, header/footer, and Jinja2 blocks
- Custom CSS with article styling, modal support, results grid, and status indicators
- Web router serving home page at / with solvents and presets context
- Static file serving integrated into FastAPI app

## Task Commits

Each task was committed atomically:

1. **Task 1: Create template and static directory structure** - `c8b89f3` (feat)
2. **Task 2: Create web router with home route** - `d24aa6d` (feat)
3. **Task 3: Integrate templates and web router into app** - `efc6198` (feat)

## Files Created/Modified

- `src/qm_nmr_calc/api/templates/base.html` - Base layout with Pico CSS, header nav, footer
- `src/qm_nmr_calc/api/templates/submit.html` - Placeholder page extending base (replaced in 06-02)
- `src/qm_nmr_calc/api/static/css/custom.css` - Custom styles for articles, images, modals, results grid
- `src/qm_nmr_calc/api/routers/web.py` - Web UI router with / home route
- `src/qm_nmr_calc/api/routers/__init__.py` - Updated exports to include web module
- `src/qm_nmr_calc/api/app.py` - Added StaticFiles mount and web router inclusion

## Decisions Made

- **Pico CSS via CDN:** Chose CDN delivery over local copy for simplicity and caching benefits
- **Static mount order:** Mounted static files before router includes to avoid route conflicts with catch-all patterns
- **Template path resolution:** Used `Path(__file__).resolve().parent.parent / "templates"` for reliable relative path resolution

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- Template infrastructure ready for submission form implementation in 06-02
- Custom CSS includes styles for status indicators, results grid, and modals that will be used in later plans
- Solvents and presets already passed to template context, ready for form dropdowns

---
*Phase: 06-web-ui*
*Completed: 2026-01-20*
