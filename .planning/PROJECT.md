# qm-nmr-calc

## What This Is

An asynchronous web service for running NMR quantum mechanical calculations on organic molecules. Users submit molecular structures (SMILES or file uploads), calculations run in the background using NWChem with COSMO solvation, and users retrieve predicted chemical shifts with interactive 3D visualization and optimized geometries when ready.

## Core Value

Reliable async NMR predictions with full control over calculation parameters -- submit a molecule, get back accurate 1H/13C shifts without babysitting long-running calculations.

## Current Milestone: v2.4 Docker Deployment

**Goal:** Make qm-nmr-calc deployable with `docker compose up` — pre-built images, auto-HTTPS, production-ready defaults.

**Target features:**
- Docker Compose stack (app + worker + Caddy reverse proxy)
- Pre-built images published to GitHub Container Registry (GHCR)
- NWChem, CREST, xTB pre-installed in worker image
- Auto-HTTPS via Caddy with Let's Encrypt
- Health checks and restart policies for auto-recovery
- Volume persistence for job data
- Environment configuration via `.env` file
- Deployment documentation

## Current State

**Active:** v2.4 Docker Deployment
**Last shipped:** v2.3 NMReData Export (2026-02-01)

**Codebase:** ~6,000 LOC Python, ~1,800 LOC tests, ~940 LOC templates, ~2,400 LOC CSS
**Tech stack:** FastAPI, Huey (SQLite), NWChem, RDKit, 3Dmol.js, Custom CSS
**Test suite:** 285 tests (257 unit + 28 conformer/xTB)

**What works today:**
- Submit molecules (SMILES/MOL/SDF) via REST API or modern web UI
- Single-conformer or ensemble (Boltzmann-weighted) NMR predictions
- RDKit KDG conformer generation, optional CREST/xTB for production quality
- RMSD clustering + xTB ranking for efficient conformer pre-selection
- NWChem DFT calculations with B3LYP/6-311+G(2d,p)
- COSMO solvation for CHCl3, DMSO, or vacuum (gas phase)
- DELTA50-derived scaling factors (1H MAE: 0.12 ppm, 13C MAE: 2.0 ppm)
- Modern glassmorphism UI with bento grid layouts
- Interactive 3D molecule viewer with shift annotations
- Spectrum plots, annotated structure drawings, file downloads
- Step tracker with conformer progress visualization
- WCAG-compliant accessibility (keyboard nav, reduced motion, contrast)
- Mobile-optimized responsive design

## Next Milestone

**TBD** — run `/gsd:new-milestone` after v2.4 ships

Future considerations tracked in `.planning/BACKLOG.md`:
- Dark mode (color scheme, system preference detection)
- Enhanced interactivity (card expansion, drag-and-drop)
- User accounts and calculation history

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
- Conformer generation via RDKit KDG (pure distance geometry) -- v2.0
- Conformer generation via CREST/xTB (optional, auto-detected) -- v2.0
- Two-stage energy window filtering (pre-DFT and post-DFT) -- v2.0
- DFT geometry optimization on each conformer -- v2.0
- NMR shielding calculation on each conformer -- v2.0
- Boltzmann weighting by DFT energies -- v2.0
- Weighted-average chemical shift output -- v2.0
- User choice between single-conformer and ensemble modes -- v2.0
- Ensemble metadata in API responses (conformer count, populations) -- v2.0
- Conformer progress tracking in web UI -- v2.0

- Bento grid layout system for all pages -- v2.1
- Glassmorphism card styling (backdrop blur, transparency, borders) -- v2.1
- Results page redesign with prominent 3D viewer and spectra -- v2.1
- Submit page redesign with clean form layout -- v2.1
- Status page redesign with progress visualization -- v2.1
- Custom CSS framework replacing Pico CSS -- v2.1
- Responsive breakpoints for tablet and mobile -- v2.1
- Hover states and smooth transitions -- v2.1
- WCAG accessibility (contrast, focus, reduced motion) -- v2.1
- RMSD clustering for conformer pre-selection -- v2.0.1
- xTB energy ranking for conformer filtering -- v2.0.1
- Comprehensive documentation (README, installation, usage, architecture, libraries, science) -- v2.2
- DP4+ methodology writeup with full derivations and literature citations -- v2.2
- NMReData SDF export with chemical shifts and atom assignments -- v2.3
- NMReData download via REST API and web UI -- v2.3

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
- **Conformer sampling addressed**: v2.0 added Boltzmann-weighted ensemble predictions
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

| KDG over ETKDG for conformers | Avoids crystal structure bias for solution-phase NMR | Good -- E2E validated on 2-penten-1-ol |
| CREST/xTB as optional deps | App works without them, enables when detected | Good -- RDKit-only pipeline works end-to-end |
| Boltzmann weight by DFT energies | Most accurate readily available energy level | Good -- chemically reasonable shifts from E2E test |
| Pure CSS (no framework) | Modern CSS has variables, grid, layers | Good -- 2,400 LOC, no build step |
| CSS Cascade Layers | Organize styles without specificity wars | Good -- clean architecture |
| Glassmorphism with high opacity | 85-95% for WCAG contrast compliance | Good -- accessible and attractive |
| RMSD clustering for conformers | Reduce redundant DFT calculations | Good -- 40→8 conformers, 10x faster |
| xTB for conformer ranking | Better than MMFF, faster than DFT | Good -- optional with MMFF fallback |

---
*Last updated: 2026-02-02 after v2.4 Docker Deployment milestone started*
