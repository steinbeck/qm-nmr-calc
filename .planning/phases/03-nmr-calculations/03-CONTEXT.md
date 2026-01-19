# Phase 3: NMR Calculations - Context

**Gathered:** 2026-01-19
**Status:** Ready for planning

<domain>
## Phase Boundary

Produce accurate ¹H and ¹³C NMR chemical shifts with configurable quality presets. Users select a preset at submission time, and completed jobs include atom-assigned shifts. Visualization of results and result file downloads are separate phases.

</domain>

<decisions>
## Implementation Decisions

### Preset definitions
- Two presets: draft and production
- **Draft**: Prioritizes speed over accuracy — minimal basis set, fast settings for quick checks
- **Production**: Good accuracy with reasonable compute time — balanced settings, not extreme
- Default preset: **production** (users get reliable results by default, opt-in to draft for speed)
- Research needed: Investigate what ISiCLE provides natively for presets/parameters and how much can be customized

### Calculation pipeline
- Geometry optimization: **Always run before NMR shielding** — NMR accuracy depends on good 3D structure
- Conformer handling: ISiCLE handles conformer generation and Boltzmann averaging natively; also expose a **single-conformer option** for simpler/faster calculations
- Solvation: **Implicit solvent model required** — user must specify solvent per job (no default)
- Research needed: Verify ISiCLE's conformer and solvation capabilities, identify available solvent models

### Claude's Discretion
- Specific basis sets and functionals for each preset
- NWChem input file generation details
- Result parsing and shift extraction logic
- Reference compound handling for ppm conversion

</decisions>

<specifics>
## Specific Ideas

- ISiCLE already handles much of the calculation pipeline (conformers, Boltzmann averaging) — leverage rather than reinvent
- Need to verify ISiCLE's parameter customization capabilities during research phase

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 03-nmr-calculations*
*Context gathered: 2026-01-19*
