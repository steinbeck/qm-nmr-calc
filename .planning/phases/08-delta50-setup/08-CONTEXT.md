# Phase 8: DELTA50 Setup - Context

**Gathered:** 2026-01-21
**Status:** Ready for planning

<domain>
## Phase Boundary

Establish DELTA50 benchmark infrastructure: 50 molecules with experimental shifts, and a runner capable of executing the full calculation matrix (50 molecules × 2 functionals × 2 solvents). This phase sets up the data and tooling; actual execution happens in Phase 9.

</domain>

<decisions>
## Implementation Decisions

### Data source
- Download from original DELTA50 publication's Supporting Information
- Research will identify the exact paper and data source
- Deal with what's actually available — don't overcomplicate for hypotheticals

### Storage format
- Benchmark data committed to repository (not external download at runtime)
- Output must support chemistry-friendly formats: SDF, NMReData
- Internal storage format: Claude's discretion

### Benchmark runner
- Failure handling: Skip failed molecules and continue (log failures for review)
- Execution mode: Sequential — one calculation at a time
- Resume support: Yes — check for existing results and skip completed calculations
- Invocation: Both CLI (`python -m qm_nmr_calc.benchmark run`) and importable function

### Output organization
- Directory structure: By molecule (`molecule_01/B3LYP_CHCl3/`, `molecule_01/WP04_DMSO/`, etc.)
- Raw files: Keep everything — .nw input, .out output, and parsed JSON
- Summary: Yes — aggregated CSV/JSON with all molecules, methods, and calculated shifts
- Location: `data/benchmark/` (data/benchmark/delta50/, data/benchmark/results/)

### Claude's Discretion
- Internal data structure for molecule/shift storage
- Exact CLI argument design
- Progress reporting verbosity
- Summary file format details (CSV vs JSON vs both)

</decisions>

<specifics>
## Specific Ideas

No specific requirements — open to standard approaches for benchmark infrastructure.

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 08-delta50-setup*
*Context gathered: 2026-01-21*
