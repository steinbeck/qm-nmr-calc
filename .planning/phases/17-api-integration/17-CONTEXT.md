# Phase 17: API Integration and Progress Tracking - Context

**Gathered:** 2026-01-28
**Status:** Ready for planning

<domain>
## Phase Boundary

User-facing ensemble mode for v2.0 conformational sampling. Users can choose between single-conformer (v1.x behavior) and ensemble mode, select RDKit KDG or CREST as the conformer generation method, and track progress through long-running ensemble calculations. API returns Boltzmann-weighted average shifts with full per-conformer detail. Web UI displays progress and allows conformer selection in 3D viewer.

</domain>

<decisions>
## Implementation Decisions

### Job Submission Options
- Default mode: **Ensemble mode** (not single-conformer) — users opt OUT to single-conformer if they want speed
- Default conformer method: **RDKit KDG** — always available, fast, good enough for most molecules
- Mode and method dropdowns: **Always visible** on submission form (not hidden in advanced section)
- Energy window thresholds (pre-DFT 6 kcal/mol, post-DFT 3 kcal/mol): **Show in advanced section** — power users can tweak

### Progress Tracking
- Granularity: **Detailed with ETA** — stage + conformer counts + estimated time remaining
- API status response: **Per-conformer status** — array with each conformer's ID and status (pending/running/complete/failed)
- Update mechanism: **Polling** — UI polls every few seconds (existing pattern)
- Poll interval: **5 seconds** — current default, responsive but not excessive

### Result Presentation
- API response detail: **Full per-conformer data** — metadata + per-conformer shifts for advanced users (larger response)
- 3D viewer geometry: **User-selectable** — dropdown to switch between conformers in the viewer
- Conformer selector info: **ID + energy + population** — e.g., "conf_001: 0.0 kcal/mol (45.2%)"
- Shift comparison: **Side-by-side table** — full table with columns for each conformer's shifts

### Error Handling
- CREST unavailable: **Warn and fallback** — proceed with RDKit, include warning in response that CREST wasn't used
- Partial conformer failures: **Show results with warning** — return averaged shifts + "N conformers failed, results based on M successful"
- CREST timeout message: **User-friendly guidance** — "Conformer search took too long. Try RDKit method or a smaller molecule."
- Failed job retry: **Yes, pre-filled form** — link to submission form with molecule pre-filled, suggest alternative settings

### Claude's Discretion
- Exact ETA calculation algorithm
- Progress bar visual style
- API response structure for per-conformer data
- Form layout and spacing
- Temperature parameter exposure (default 298.15 K works for most cases)

</decisions>

<specifics>
## Specific Ideas

No specific product references — open to standard approaches based on existing UI patterns.

</specifics>

<deferred>
## Deferred Ideas

- **UI Redesign with Bento Grid + Glassmorphism** — User wants a visual overhaul using Bento Grid layout with Glassmorphism effects. This is a distinct design/presentation phase, separate from the functional API integration work. Recommend as Phase 18 after v2.0 functionality is complete.

</deferred>

---

*Phase: 17-api-integration*
*Context gathered: 2026-01-28*
