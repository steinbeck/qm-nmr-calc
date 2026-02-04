# Project Milestones: qm-nmr-calc

## v2.4 Docker Deployment (Shipped: 2026-02-03)

**Delivered:** Production-ready Docker deployment with `docker compose up`, pre-built GHCR images, auto-HTTPS via Caddy, and comprehensive deployment documentation.

**Phases completed:** 35-40 (6 phases, 8 plans total)

**Key accomplishments:**

- Worker container with NWChem, CREST, xTB pre-installed (2.1GB image, amd64)
- API container with FastAPI multi-stage build (~733MB image, amd64+arm64)
- Docker Compose orchestration with persistent volumes, health checks, graceful shutdown
- Caddy reverse proxy with automatic Let's Encrypt HTTPS certificates
- GitHub Actions CI/CD workflow publishing to GHCR on release tags
- Comprehensive deployment guide (460 lines) covering VPS setup, troubleshooting, backup/restore
- Docker Quick Start as primary README getting started path

**Stats:**

- ~1,065 lines of Docker/deployment artifacts
- 6 phases, 8 plans
- ~2 days (2026-02-02 to 2026-02-03)

**Git range:** `docs(35)` → `docs(40)`

**What's next:** TBD — run `/gsd:new-milestone` to define next focus area.

---

## v2.3 NMReData Export (Shipped: 2026-02-01)

**Delivered:** Machine-readable export of NMR prediction results in NMReData standard format for interoperability with NMR analysis tools.

**Phases completed:** 32-34 (3 phases, 3 plans total)

**Key accomplishments:**

- NMReData-compliant SDF generation with all 8 required tags (VERSION, LEVEL, SOLVENT, TEMPERATURE, ASSIGNMENT, FORMULA, SMILES, ID)
- REST API download endpoint `GET /api/v1/jobs/{job_id}/nmredata.sdf` with proper HTTP headers (chemical/x-mdl-sdfile)
- Web UI download button on results page following glassmorphism design system
- Automatic solvent mapping (chcl3→CDCl3, dmso→(CD3)2SO, vacuum→vacuum)
- Ensemble mode support with Boltzmann-averaged shifts and provenance metadata
- Comprehensive test coverage: 44 tests (36 unit + 8 integration)

**Stats:**

- ~1,023 lines added (270 LOC module, 641 LOC tests, 112 LOC endpoint/UI)
- 3 phases, 3 plans
- 1 day (2026-02-01)

**Git range:** `ac9409a` → `094fda8`

**What's next:** Deploy to production. Future considerations: NMREDATA_J (coupling constants), NMREDATA_INCHI, per-conformer export.

---

## v2.2 Documentation (Shipped: 2026-02-01)

**Delivered:** Comprehensive documentation for academic researchers and developers, including DP4+ science writeup with full derivations.

**Phases completed:** 25-31 (7 phases, 10 plans total)

**Key accomplishments:**

- README overhaul with Mermaid architecture diagram and quick start section
- Comprehensive installation guide with multi-distro instructions and CREST/xTB setup
- Usage guide with web UI workflow and 26 curl API examples
- Technical architecture docs with 7 Mermaid diagrams (stack, data flow, job lifecycle)
- Library documentation for RDKit, NWChem, Huey, 3Dmol.js, SmilesDrawer, CREST/xTB
- DP4+ science writeup with full derivations, 9 DOI citations, and accuracy metrics
- Standardized cross-references and verified links across all documentation

**Stats:**

- ~4,100 lines of documentation
- 7 phases, 10 plans
- 2 days from start to ship (2026-01-31 to 2026-02-01)

**Git range:** `feat(25-01)` -> `feat(31-01)`

**What's next:** Deploy to production. Future considerations: dark mode, enhanced interactivity, user accounts.

---

## v2.1 UI Redesign (Shipped: 2026-01-31)

**Delivered:** Modern bento grid layout with glassmorphism effects, WCAG-compliant accessibility, and mobile-optimized performance.

**Phases completed:** 18-23 (6 phases, 15 plans total)

**Key accomplishments:**

- Custom CSS architecture replacing Pico CSS (cascade layers, 62 design tokens, BEM components)
- Glassmorphism effects with backdrop-filter blur and Safari/accessibility fallbacks
- Bento grid layout system with asymmetric card arrangements
- Results page redesign with hero 3D viewer and organized data cards
- Submit page redesign with two-column form and SmilesDrawer molecule preview
- Status page step tracker with conformer progress visualization
- WCAG accessibility compliance (4.5:1 contrast, focus indicators, reduced motion support)
- Mobile-optimized performance (reduced blur effects, responsive breakpoints, touch-safe interactions)

**Stats:**

- 164 files modified (+17,302 lines)
- 2,443 lines of CSS
- 6 phases, 15 plans
- 3 days from start to ship (2026-01-29 to 2026-01-31)

**Git range:** `feat(18-01)` -> `feat(23-01)`

**What's next:** Deploy to production. Future considerations: dark mode, card interactivity, user accounts.

---

## v2.0 Conformational Sampling (Shipped: 2026-01-28)

**Delivered:** Boltzmann-weighted ensemble NMR predictions with RDKit and optional CREST conformer generation.

**Phases completed:** 12-17 (6 phases, 19 plans total)

**Key accomplishments:**

- RDKit KDG conformer generation (pure distance geometry, no crystal bias)
- CREST/xTB integration for production-quality conformer searching (when available)
- Boltzmann weighting by DFT energies with exp-normalize numerical stability
- Multi-conformer NWChem pipeline with post-DFT energy filtering
- Ensemble metadata in API (conformer count, populations, energy range)
- Web UI progress tracking for ensemble jobs

**Stats:**

- 6 phases, 19 plans
- 3 days from start to ship (2026-01-26 to 2026-01-28)

**Git range:** `feat(12-01)` -> `feat(17-05)`

**What's next:** v2.1 UI Redesign

---

## v2.0.1 Conformer Pre-selection (Shipped: 2026-01-30)

**Delivered:** Performance optimization reducing DFT workload from ~40 to ~8 conformers using RMSD clustering and xTB ranking.

**Phases completed:** 24 (1 phase, 3 plans total)

**Key accomplishments:**

- RMSD-based Butina clustering for conformer deduplication
- xTB (GFN2) semi-empirical energy ranking with ALPB solvation
- Graceful MMFF fallback when xTB not available
- Hexanol ensemble: 28 → 3 conformers, completing in <1 hour vs 10+ hours

**Stats:**

- 1 phase, 3 plans
- 18 minutes total execution

**Git range:** `feat(24-01)` -> `feat(24-03)`

**What's next:** v2.1 UI Redesign

---

## v1.1 Accurate Chemical Shifts (Shipped: 2026-01-25)

**Delivered:** Direct NWChem integration with DELTA50-derived scaling factors achieving literature-comparable NMR prediction accuracy (1H: 0.12 ppm, 13C: 2.0 ppm MAE).

**Phases completed:** 7-11.2 (8 phases, 21 plans total)

**Key accomplishments:**

- Direct NWChem I/O replacing ISiCLE runtime dependency (input generation, output parsing, COSMO solvation)
- DELTA50 benchmark infrastructure with 150 NMR calculations (100 solvated + 50 vacuum)
- OLS regression-derived scaling factors for 6 solvent/nucleus combinations (CHCl3/DMSO/vacuum x 1H/13C)
- Interactive 3Dmol.js visualization on benchmark viewer and job results pages
- Vacuum (gas-phase) calculation support with comparable accuracy to solvated
- API transparency with scaling factor metadata (source, expected MAE)

**Stats:**

- 185 files created/modified (+22,609 / -1,561 lines)
- 5,432 lines Python, 1,417 lines tests, 892 lines templates
- 8 phases, 21 plans, 116 commits
- 5 days from start to ship (2026-01-21 to 2026-01-25)

**Git range:** `701f613` -> `f7e1b13`

**What's next:** v1.2 Conformational Sampling -- Boltzmann-weighted ensemble averaging for flexible molecules

---

## v1.0 Core NMR Service (Shipped: 2026-01-20)

**Delivered:** Working async NMR prediction service with REST API, web UI, calculation presets, and email notifications.

**Phases completed:** 1-6 (6 phases, 16 plans total)

**Key accomplishments:**

- Async job queue with SQLite-backed Huey for crash-safe background calculations
- REST API with OpenAPI docs for molecule submission (SMILES/MOL/SDF) and result retrieval
- NMR calculation pipeline with draft/production presets
- Spectrum plots and annotated structure drawings
- Clean web UI with Pico CSS for submission, status polling, and results viewing
- Email notification system (opt-in)

**Stats:**

- 6 phases, 16 plans
- 2 days from start to ship (2026-01-19 to 2026-01-20)

**What's next:** v1.1 Accurate Chemical Shifts

---
