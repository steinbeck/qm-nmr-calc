# Phase 5: Visualization - Context

**Gathered:** 2026-01-20
**Status:** Ready for planning

<domain>
## Phase Boundary

Generate visual representations of NMR calculation results. Two outputs: stick spectrum plots showing chemical shifts as vertical lines, and annotated molecular structure images with shift values labeled on atoms. Both available for download as static image files.

</domain>

<decisions>
## Implementation Decisions

### Spectrum plot style
- Stick spectrum (vertical lines at each chemical shift) for both 1H and 13C
- Black peaks on white background (classic publication style)
- No peak labels on spectrum — values shown only on annotated structure
- Auto-fit x-axis range to data (with padding), not fixed ranges
- X-axis goes right to left (high ppm to low ppm) — standard NMR convention
- Small nucleus label ("1H" or "13C") in corner of each plot

### Axis labeling
- Minimal: just "ppm" on x-axis
- No y-axis label (intensity is relative)
- No calculation metadata on images (functional, basis set, solvent) — keep clean

### Structure annotation
- Single image showing both 1H and 13C shift values
- Labels positioned adjacent to atoms, with leader lines if needed for clarity
- 2 decimal places for shift values (e.g., 7.26 ppm)
- Label every atom individually, even chemically equivalent atoms

### Image output format
- Provide both PNG and SVG formats for all images
- PNG at 300 dpi (print/publication quality)
- Three separate files: 1H spectrum, 13C spectrum, annotated structure

### Claude's Discretion
- Exact peak height/intensity representation
- Leader line styling and collision avoidance algorithm
- Font choices and sizing
- Padding amount for auto-fit axis ranges
- Color differentiation for 1H vs 13C labels on structure (if needed for readability)

</decisions>

<specifics>
## Specific Ideas

- Interactive hover-to-reveal shift values — deferred to Phase 6 (Web UI), not part of static images
- Future sophisticated 1H spectrum with proper Lorentzian line shapes and coupling patterns — separate future phase

</specifics>

<deferred>
## Deferred Ideas

- **Sophisticated 1H spectrum visualization** — Lorentzian peak shapes with proper linewidths and coupling patterns (J-coupling multiplets). Add as a "nice-to-have" phase after the v1.0 milestone is complete.
- **Interactive spectrum** — Hover over peaks to see shift values. Belongs in Phase 6 (Web UI) with JavaScript rendering.

</deferred>

---

*Phase: 05-visualization*
*Context gathered: 2026-01-20*
