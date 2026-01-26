# Phase 11: Production Integration - Context

**Gathered:** 2026-01-23
**Status:** Ready for planning

<domain>
## Phase Boundary

Wire the NWChem-derived B3LYP scaling factors from Phase 10 into the production calculation pipeline and remove ISiCLE as a runtime dependency. Users get accurate solvated NMR predictions using factors derived from DELTA50 benchmark data.

**Scope adjustment:** WP04 functional is deferred since NWChem doesn't support it natively. This phase focuses on B3LYP scaling factors only.

</domain>

<decisions>
## Implementation Decisions

### Factor Application
- Replace TMS reference entirely with regression factors (shift = slope * shielding + intercept)
- Load scaling factors from JSON at runtime (data/benchmark/delta50/scaling_factors.json)
- Reject calculations if no scaling factor exists for the solvent/functional combination (no fallback to TMS)
- Include factor metadata in API response (factor source, R², MAE)

### ISiCLE Removal
- Delete all ISiCLE imports, wrapper code, and dependency completely
- Remove isicle_version field from models and API responses
- Keep attribution in LICENSE/README only (not scattered through code)

### API/UI Changes
- B3LYP-only functional (no selector needed - single option)
- Show expected MAE in results metadata (e.g., "1H: ±0.12 ppm typical")
- Methodology documentation accessible via help menu/button (not prominent banner)

### Claude's Discretion
- Exact structure of factor loading code
- How to handle edge cases in factor lookup
- Documentation content and placement

</decisions>

<specifics>
## Specific Ideas

- Factor metadata should include R² and MAE from the scaling factor derivation
- Accuracy shown as "±X.XX ppm typical" format
- Documentation should explain the DELTA50 benchmark methodology

</specifics>

<deferred>
## Deferred Ideas

- WP04 functional support — requires custom implementation or different QM package
- Other solvents beyond CHCl3 and DMSO — need benchmark data first
- Per-atom uncertainty estimates — could be added based on residual analysis

</deferred>

---

*Phase: 11-production-integration*
*Context gathered: 2026-01-23*
