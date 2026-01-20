# Phase 6: Web UI - Context

**Gathered:** 2026-01-20
**Status:** Ready for planning

<domain>
## Phase Boundary

Browser interface for submitting molecules and viewing results without API calls. Three pages: submission form, status page with auto-refresh, and results page with visualizations and downloads. The web UI wraps the existing API — no new backend functionality.

</domain>

<decisions>
## Implementation Decisions

### Visual style and layout
- Minimalist/clean aesthetic — lots of whitespace, subtle colors, professional scientific feel
- Blue/white color scheme — blue accents on white/light gray (classic scientific software look)
- Simple text header with app name (e.g., "QM NMR Calculator") — no logo needed
- CSS framework: Claude's discretion (whatever achieves clean look efficiently)

### Submission form
- Both SMILES input and file upload visible on same form (not tabs)
- Standard dropdown menus for preset and solvent selection
- Email notification field visible but clearly marked as optional
- On successful submission: redirect immediately to status page

### Status page behavior
- Auto-refresh via polling while job is running (no manual refresh button needed)
- Display job status plus elapsed time (how long job has been running)
- On job completion: auto-redirect to results page
- On job failure: stay on status page, show error message clearly

### Results presentation
- All images displayed side by side: annotated structure + 1H spectrum + 13C spectrum
- Click any image to enlarge in modal/lightbox
- Download section grouped at bottom of page (not per-image buttons)
- Show calculation metadata section: preset, solvent, functional, basis set used

### Claude's Discretion
- Specific CSS framework choice (Tailwind, Bootstrap, plain CSS)
- Exact layout grid/responsive behavior
- Polling interval for auto-refresh
- Modal/lightbox implementation approach
- Form validation UX details

</decisions>

<specifics>
## Specific Ideas

- No specific product references — standard clean scientific web app aesthetic

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 06-web-ui*
*Context gathered: 2026-01-20*
