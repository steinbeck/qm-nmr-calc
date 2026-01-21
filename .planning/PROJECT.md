# qm-nmr-calc

## What This Is

An asynchronous web service for running NMR quantum mechanical calculations on organic molecules. Users submit molecular structures (SMILES or file uploads), calculations run in the background using ISiCLE/NWChem, and users retrieve predicted chemical shifts and optimized geometries when ready.

## Core Value

Reliable async NMR predictions with full control over calculation parameters — submit a molecule, get back accurate ¹H/¹³C shifts without babysitting long-running calculations.

## Current Milestone: v1.1 Accurate Chemical Shifts

**Goal:** Replace ISiCLE dependency and derive NWChem-specific scaling factors from DELTA50 benchmark, enabling accurate solvated NMR predictions.

**Target features:**
- Custom NWChem wrapper (no ISiCLE runtime dependency)
- Working COSMO solvation (currently broken - solvent param ignored)
- NWChem-derived scaling factors for CHCl3 and DMSO solvents
- B3LYP (general) + WP04 (optimized ¹H) functionals
- Publication-quality benchmark data from DELTA50 dataset

## Previous Milestone: v1.0 Core NMR Service ✓

**Delivered:** Working async NMR prediction service with REST API and web UI.
- Submit molecules (SMILES/MOL) via API or web UI
- ISiCLE/NWChem calculation pipeline
- Results: raw files, JSON with atom-assigned shifts, spectrum plot, annotated structure drawing
- Calculation presets (draft/production)
- Job status polling + email notifications
- Clean web interface with Pico CSS

## Requirements

### Validated (v1.0)

- ✓ Submit molecules via SMILES string or structure file (SDF/MOL)
- ✓ Queue calculations for background processing
- ✓ Check job status and retrieve results via API
- ✓ Return predicted ¹H and ¹³C NMR chemical shifts with atom assignments
- ✓ Return optimized molecular geometry
- ✓ Return raw NWChem output files
- ✓ Visual spectrum plot generation
- ✓ Annotated structure drawing showing shifts on atoms
- ✓ Calculation presets (draft/production)
- ✓ Email notification when calculation completes
- ✓ Clean, usable web UI for submitting jobs and viewing results
- ✓ REST API for programmatic access

### Active (v1.1)

- [ ] Custom NWChem input/output handling (no ISiCLE runtime dependency)
- [ ] Working COSMO solvation for CHCl3 and DMSO
- [ ] NWChem-derived scaling factors from DELTA50 benchmark
- [ ] WP04 functional option for improved ¹H accuracy
- [ ] Accept pre-optimized XYZ/SDF geometries (skip geometry opt)

### Out of Scope

- Multi-user authentication — single user for now, architecture supports adding later
- Other nuclei (¹⁵N, ³¹P, ¹⁹F) — ¹H/¹³C only for now
- Batch submission UI — can be added later, API could support
- Mobile interface — web UI is desktop-focused
- Real-time calculation streaming — poll for completion
- Public benchmark API — benchmark tooling is internal only

## Context

**ISiCLE**: Pacific Northwest National Laboratory's in silico chemical library engine (https://github.com/pnnl/isicle/). Handles geometry optimization and NMR shielding calculations via NWChem. Development has slowed; we'll maintain a minimal fork for bug fixes while keeping close to upstream.

**NWChem**: Open-source quantum chemistry package that performs the actual DFT calculations. Needs to be installed on the target VM.

**Calculation time**: NMR QM calculations can take minutes to hours depending on molecule size and theory level. Async architecture is essential — API calls return immediately with job IDs.

**Project structure**:
- `qm-nmr-calc` (this repo): API server, job queue, web UI, orchestration scripts
- ISiCLE fork: Calculation engine dependency
- NWChem: Installed system dependency

## Constraints

- **Deployment**: Single cloud VM — architecture should work without distributed infrastructure
- **Storage**: Filesystem-based — no database requirement for v1, but job metadata needs organization
- **Python ecosystem**: ISiCLE is Python, backend should stay in Python for clean integration
- **ISiCLE compatibility**: Fork stays minimal to ease upstream tracking

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| FastAPI for backend | Modern async Python, good for APIs, natural fit with ISiCLE | — Pending |
| Filesystem storage | Simple for single VM, avoids database ops overhead | — Pending |
| Lightweight job queue | Single VM doesn't need Celery/Redis complexity | — Pending |
| ISiCLE as dependency | Minimal fork, wrap rather than modify | — Pending |

---
*Last updated: 2026-01-21 after v1.1 milestone definition*
