# qm-nmr-calc

## What This Is

An asynchronous web service for running NMR quantum mechanical calculations on organic molecules. Users submit molecular structures (SMILES or file uploads), calculations run in the background using ISiCLE/NWChem, and users retrieve predicted chemical shifts and optimized geometries when ready.

## Core Value

Reliable async NMR predictions with full control over calculation parameters — submit a molecule, get back accurate ¹H/¹³C shifts without babysitting long-running calculations.

## Current Milestone: v1.0 Core NMR Service

**Goal:** Build a working async NMR prediction service with REST API and clean web UI.

**Target features:**
- Submit molecules (SMILES/MOL) via API or web UI
- ISiCLE/NWChem calculation pipeline
- Results: raw files, JSON with atom-assigned shifts, spectrum plot, annotated structure drawing
- Calculation presets + advanced parameter control
- Job status polling + email notifications
- Clean, presentable web interface

## Requirements

### Validated

(None yet — ship to validate)

### Active

- [ ] Submit molecules via SMILES string or structure file (SDF/MOL)
- [ ] Queue calculations for background processing
- [ ] Check job status and retrieve results via API
- [ ] Return predicted ¹H and ¹³C NMR chemical shifts with atom assignments
- [ ] Return optimized molecular geometry
- [ ] Return raw NWChem output files
- [ ] Visual spectrum plot generation
- [ ] Annotated structure drawing showing shifts on atoms
- [ ] Calculation presets (fast/draft, publication quality, solvation options)
- [ ] Advanced parameter control (basis sets, functionals, solvation models)
- [ ] Email notification when calculation completes
- [ ] Clean, usable web UI for submitting jobs and viewing results
- [ ] REST API for programmatic access
- [ ] Job ownership model (ready for multi-user extension)

### Out of Scope

- Multi-user authentication — single user for v1, architecture supports adding later
- Other nuclei (¹⁵N, ³¹P, ¹⁹F) — ¹H/¹³C only for v1
- Batch submission UI — can be added later, API could support
- Mobile interface — web UI is desktop-focused
- Real-time calculation streaming — poll for completion

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
*Last updated: 2026-01-19 after v1.0 milestone definition*
