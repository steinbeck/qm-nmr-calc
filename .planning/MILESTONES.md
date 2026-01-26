# Project Milestones: qm-nmr-calc

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
