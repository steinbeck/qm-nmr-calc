# Phase 7: NWChem Integration - Context

**Gathered:** 2026-01-21
**Status:** Ready for planning

<domain>
## Phase Boundary

Replace ISiCLE wrapper with direct NWChem handling — input file generation, output parsing, and working COSMO solvation. Users can submit pre-optimized geometries to skip geometry optimization. ISiCLE runtime dependency removed.

</domain>

<decisions>
## Implementation Decisions

### ISiCLE Code Adaptation
- Extract as needed — start writing our version, refer to ISiCLE when stuck on specifics
- Minimal extraction — only what we need for our methods (input generation, output parsing)
- Remove ISiCLE dependency completely — clean break, no dev dependency either
- Attribution in README only — single acknowledgment, no in-code comments

### Module Organization
- Replace isicle_wrapper.py immediately — no coexistence period, no feature flag
- SMILES-to-3D conversion lives in the NWChem module — full pipeline in one place
- Redesign API freely — best interface for new implementation, update callers as needed

### Pre-optimized Geometry Handling
- Explicit flag required — API parameter like `skip_optimization=true`, not inferred from format
- Supported formats: XYZ and SDF
- NMR shielding always runs — can only skip geometry optimization, not NMR calculation
- Basic validation — check file parses correctly and has valid atom coordinates

### COSMO Solvation
- Support CHCl3 and DMSO only — matches v1.1 milestone scope (scaling factors derived for these)
- Solvent is required parameter — no default, forces user to explicitly choose
- Reject unsupported solvents with 400 error — list valid options in error message
- Apply COSMO to both steps — geometry optimization AND NMR shielding (more accurate)

### Claude's Discretion
- Code organization (single module vs split files)
- Exact function signatures and internal API design
- NWChem input file formatting details
- Error handling and logging approach

</decisions>

<specifics>
## Specific Ideas

- Study ISiCLE's `isicle/qm/` module for NWChem input/output patterns when needed
- Current isicle_wrapper.py has `cosmo=False` hardcoded — this is the bug we're fixing

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 07-nwchem-integration*
*Context gathered: 2026-01-21*
