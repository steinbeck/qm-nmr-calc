# Phase 10: Scaling Factors - Context

**Gathered:** 2026-01-23
**Status:** Ready for planning

<domain>
## Phase Boundary

Derive NWChem-specific scaling factors from DELTA50 benchmark data and validate against external nmrshiftdb compounds. Factors convert calculated shielding values to predicted chemical shifts via linear regression. This phase produces validated scaling factors ready for production integration (Phase 11).

</domain>

<decisions>
## Implementation Decisions

### Regression approach
- Simple linear regression: shift = slope × shielding + intercept
- Remove outliers using residual-based method (fit, identify >3 std dev, refit)
- Fit on all DELTA50 data (no train/test split)
- External validation against nmrshiftdb compounds (10-20 manually curated)

### Factor organization
- Per nucleus/solvent granularity: separate factors for 1H-CHCl3, 1H-DMSO, 13C-CHCl3, 13C-DMSO
- Functional-aware lookup: key is (functional, basis_set, nucleus, solvent)
- Include metadata with each factor set: slope, intercept, R², MAE, RMSD, n_points
- Design supports future WP04 factors alongside B3LYP

### Validation metrics
- Report both MAE and RMSD
- Primary validation: external accuracy on nmrshiftdb compounds (small set, 10-20)
- Validation workflow: derive factors from DELTA50 → apply to nmrshiftdb → report prediction accuracy
- Include 95% confidence intervals for slope/intercept
- Bootstrap for MAE/RMSD uncertainty estimates

### Output format
- Analysis code in qm_nmr_calc.benchmark module (analysis.py)
- Human-readable SCALING-FACTORS.md report with tables and plots
- Plots: regression scatter plots, residual distributions (PNG)
- Claude's discretion on factor storage (Python constants vs JSON)

### Publication data
- Capture per-compound statistics (individual MAE/deviations for supplementary tables)
- Correlation matrices across molecule classes/functional groups
- Full statistical rigor: confidence intervals, p-values where appropriate
- Data structured for main text tables and supplementary material

### Claude's Discretion
- Factor storage format (Python dict vs JSON file)
- Specific outlier threshold (starting point: 3 std dev)
- Plot styling and layout
- nmrshiftdb compound selection criteria

</decisions>

<specifics>
## Specific Ideas

- Planning to write a publication about this work — analysis should capture comprehensive, publishable data
- Follow DELTA50 methodology: derive from benchmark, validate on external set
- Current system uses CHESHIRE TMS reference factors — new factors replace this approach

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 10-scaling-factors*
*Context gathered: 2026-01-23*
