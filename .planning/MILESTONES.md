# Project Milestones: qm-nmr-calc

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
