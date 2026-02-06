# qm-nmr-calc

## What This Is

An asynchronous web service for running NMR quantum mechanical calculations on organic molecules. Users submit molecular structures (SMILES or file uploads), calculations run in the background using NWChem with COSMO solvation, and users retrieve predicted chemical shifts with interactive 3D visualization and optimized geometries when ready.

## Core Value

Reliable async NMR predictions with full control over calculation parameters -- submit a molecule, get back accurate 1H/13C shifts without babysitting long-running calculations.

## Current State

**Last shipped:** v2.7 Automated GCP Deployment (2026-02-06)

**Codebase:** ~7,300 LOC Python, ~3,050 LOC tests, ~950 LOC templates, ~2,400 LOC CSS, ~4,800 LOC docs, ~2,170 LOC GCP scripts
**Tech stack:** FastAPI, Huey (SQLite), NWChem, RDKit, 3Dmol.js, Custom CSS, Docker, Caddy
**Test suite:** 415 tests
**Docker:** Worker 2.1GB, API ~733MB, Caddy reverse proxy, GHCR publishing
**GCP:** Automated TOML-config deployment, dynamic pricing/machine selection, lifecycle management

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
- Automated GCP spot deployment from TOML config (zero prompts)
- Dynamic spot pricing discovery across all GCP regions
- Auto machine type selection matching CPU/RAM requirements
- HTTP-only fire-up-and-burn cloud deployment pattern

## Current Milestone

No active milestone. v2.7 shipped 2026-02-06.

## Future Considerations

- Dark mode (color scheme, system preference detection)
- Enhanced interactivity (card expansion, drag-and-drop)
- User accounts and calculation history
- Kubernetes/Helm charts for HPC environments
- Prometheus metrics and monitoring

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
- Deploy with `docker compose up -d` command -- v2.4
- Worker container with NWChem, CREST, xTB pre-installed -- v2.4
- API container with FastAPI and health checks -- v2.4
- Persistent volumes for job data and queue -- v2.4
- Auto-HTTPS via Caddy and Let's Encrypt -- v2.4
- Pre-built images on GHCR with CI/CD publishing -- v2.4
- Comprehensive deployment documentation -- v2.4
- ARM64 Docker worker container via conda-forge -- v2.5
- Multi-arch images on GHCR (amd64 + arm64) -- v2.5
- GCP Spot VM deployment with lifecycle scripts -- v2.6
- Persistent disk across VM lifecycle -- v2.6
- TOML config-driven GCP deployment (zero prompts) -- v2.7
- Auto-discover cheapest spot region via CloudPrice.net -- v2.7
- Dynamic machine type selection and Docker resource limits -- v2.7
- Dry-run mode for deployment preview -- v2.7
- HTTP-only deployment (no HTTPS/domain needed) -- v2.7
- Conformer progress tracking display fix -- v2.7

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
| Worker amd64-only | CREST/xTB lack arm64 binaries | Good -- API supports both |
| Caddy for reverse proxy | Auto-HTTPS, zero config | Good -- Let's Encrypt works seamlessly |
| GITHUB_TOKEN for GHCR | No PAT management needed | Good -- automatic authentication |
| Docker as primary install | Most users want easy deployment | Good -- 5-minute quick start |

| Non-interactive GCP deploy | v2.6 scripts required manual intervention every time | Good -- single command, zero prompts |
| Config-file-driven deployment | Replaces interactive prompts with declarative config | Good -- TOML with Pydantic validation |
| Auto-discover cheapest spot region | Queries GCP pricing API instead of hardcoded defaults | Good -- CloudPrice.net + fallback |
| HTTP-only GCP deployment | Fire-up-and-burn pattern, no domain needed | Good -- simpler, faster setup |
| Modular bash library architecture | Composable config/pricing/machine/infra libraries | Good -- clean, reusable |

---
*Last updated: 2026-02-06 after v2.7 milestone shipped*
