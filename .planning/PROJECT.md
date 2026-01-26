# qm-nmr-calc

## What This Is

An asynchronous web service for running NMR quantum mechanical calculations on organic molecules. Users submit molecular structures (SMILES or file uploads), calculations run in the background using NWChem with COSMO solvation, and users retrieve predicted chemical shifts with interactive 3D visualization and optimized geometries when ready.

## Core Value

Reliable async NMR predictions with full control over calculation parameters -- submit a molecule, get back accurate 1H/13C shifts without babysitting long-running calculations.

## Current State

**Shipped:** v1.1 Accurate Chemical Shifts (2026-01-25)

**Codebase:** 5,432 LOC Python, 1,417 LOC tests, 892 LOC templates
**Tech stack:** FastAPI, Huey (SQLite), NWChem, RDKit, 3Dmol.js, Pico CSS
**Test suite:** 95 tests (92 unit + 3 NWChem integration)

**What works today:**
- Submit molecules (SMILES/MOL/SDF) via REST API or web UI
- NWChem DFT calculations with B3LYP/6-311+G(2d,p)
- COSMO solvation for CHCl3, DMSO, or vacuum (gas phase)
- DELTA50-derived scaling factors (1H MAE: 0.12 ppm, 13C MAE: 2.0 ppm)
- Interactive 3D molecule viewer with shift annotations (3Dmol.js)
- Spectrum plots, annotated structure drawings, file downloads
- Calculation presets (draft/production)
- Job status polling + email notifications

**Known limitation:** Single-conformer predictions only (being addressed in v2.0). Flexible molecules may have inaccurate predictions because the experimental NMR spectrum is a population-weighted average across all thermally accessible conformers.

## Current Milestone: v2.0 Conformational Sampling

**Goal:** Boltzmann-weighted ensemble averaging for flexible molecules -- generate multiple conformers, compute NMR on each, weight by DFT energies for population-averaged shifts.

**Target features:**
- RDKit KDG conformer generation (built-in, no new deps)
- CREST/xTB conformer generation (optional, auto-detected)
- Two-stage energy filtering (pre-DFT wide window, post-DFT tight window)
- Full DFT geometry optimization on each conformer
- Boltzmann-weighted averaging of NMR shifts
- User choice: single-conformer (v1.x behavior) or ensemble mode
- User choice: RDKit or CREST conformer method

## Requirements

### Validated

- Submit molecules via SMILES string or structure file (SDF/MOL) -- v1.0
- Queue calculations for background processing -- v1.0
- Check job status and retrieve results via API -- v1.0
- Return predicted 1H and 13C NMR chemical shifts with atom assignments -- v1.0
- Return optimized molecular geometry -- v1.0
- Return raw NWChem output files -- v1.0
- Visual spectrum plot generation -- v1.0
- Annotated structure drawing showing shifts on atoms -- v1.0
- Calculation presets (draft/production) -- v1.0
- Email notification when calculation completes -- v1.0
- Clean, usable web UI for submitting jobs and viewing results -- v1.0
- REST API for programmatic access -- v1.0
- Custom NWChem input/output handling (no ISiCLE runtime dependency) -- v1.1
- Working COSMO solvation for CHCl3, DMSO, and vacuum -- v1.1
- NWChem-derived scaling factors from DELTA50 benchmark -- v1.1
- Accept pre-optimized XYZ/SDF geometries (skip geometry opt) -- v1.1
- Interactive 3D molecule visualization with shift labels -- v1.1
- Scaling factor metadata in API responses (source, expected MAE) -- v1.1

### Active

- [ ] Conformer generation via RDKit KDG (pure distance geometry, no crystal bias)
- [ ] Conformer generation via CREST/xTB (optional system deps, auto-detected)
- [ ] Energy window filtering (initial pre-DFT and post-DFT)
- [ ] DFT geometry optimization on each conformer
- [ ] NMR shielding calculation on each conformer
- [ ] Boltzmann weighting by DFT energies
- [ ] Weighted-average chemical shift output
- [ ] User choice between single-conformer and ensemble modes
- [ ] User choice between RDKit and CREST conformer methods
- [ ] API and web UI updates for conformational sampling parameters

### Out of Scope

- Multi-user authentication -- single user for now, architecture supports adding later
- Other nuclei (15N, 31P, 19F) -- 1H/13C only for now
- Batch submission UI -- can be added later, API could support
- Mobile interface -- web UI is desktop-focused
- Real-time calculation streaming -- poll for completion
- Public benchmark API -- benchmark tooling is internal only
- WP04 functional -- NWChem doesn't support WP04 without custom compilation

## Context

**NWChem**: Open-source quantum chemistry package that performs the actual DFT calculations. Installed system dependency. Direct I/O integration (input generation, output parsing) replaces former ISiCLE dependency.

**Calculation time**: NMR QM calculations take minutes to hours depending on molecule size and theory level. Async architecture is essential -- API calls return immediately with job IDs. Conformational sampling multiplies this by the number of conformers.

**DELTA50 benchmark**: 50-molecule dataset used to derive NWChem-specific scaling factors via OLS regression. Factors stored in JSON, loaded lazily via importlib.resources.

**Conformational sampling reference**: `references/conformational_sampling_nmr_analysis.md` -- comprehensive analysis of RDKit vs CREST approaches, Boltzmann weighting implementation, energy window guidelines, and NWChem input templates for ensemble calculations.

**CREST/xTB**: Optional system dependencies for production-quality conformer searching. CREST uses metadynamics with GFN2-xTB for thorough conformational sampling. Better conformer ranking than force fields (Spearman rho ~0.39-0.47 vs DFT). ALPB solvation model available.

**Project structure**:
- `qm-nmr-calc` (this repo): API server, job queue, web UI, NWChem integration, benchmark tooling
- NWChem: Installed system dependency
- CREST/xTB: Optional system dependencies (for conformer generation)

## Constraints

- **Deployment**: Single cloud VM -- architecture should work without distributed infrastructure
- **Storage**: Filesystem-based -- no database requirement, job metadata in JSON files
- **Python ecosystem**: Backend stays in Python for NWChem/RDKit integration
- **Single conformer**: Being addressed in v2.0 -- conformational sampling with Boltzmann weighting
- **CREST/xTB optional**: App must work without CREST/xTB installed (RDKit-only fallback)

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| FastAPI for backend | Modern async Python, good for APIs | Good -- clean API with OpenAPI docs |
| Filesystem storage | Simple for single VM, avoids database ops overhead | Good -- JSON files work well |
| Huey with SQLite | Single VM doesn't need Celery/Redis complexity | Good -- crash-safe, simple |
| Direct NWChem I/O | ISiCLE development stalled, more control needed | Good -- simpler, COSMO bug fixed |
| DELTA50 regression factors | Literature-standard benchmark dataset | Good -- 1H MAE 0.12 ppm, 13C MAE 2.0 ppm |
| OLS with 3-sigma outlier removal | Standard approach for scaling factor derivation | Good -- R^2 > 0.99 for all factor sets |
| 3Dmol.js from CDN | Minimal build complexity for 3D visualization | Good -- zero build step, works well |
| Vacuum as solvent value | Cleaner than separate gas-phase parameter | Good -- consistent API surface |

| KDG over ETKDG for conformers | Avoids crystal structure bias for solution-phase NMR | -- Pending |
| CREST/xTB as optional deps | App works without them, enables when detected | -- Pending |
| Boltzmann weight by DFT energies | Most accurate readily available energy level | -- Pending |

---
*Last updated: 2026-01-26 after v2.0 milestone initialization*
