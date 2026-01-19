# Project Research Summary

**Project:** qm-nmr-calc
**Domain:** Async NMR QM Calculation Web Service
**Researched:** 2026-01-19
**Confidence:** MEDIUM-HIGH

## Executive Summary

This project fills a clear gap: no publicly available web service offers QM-level NMR predictions with a clean REST API. The recommended approach is a FastAPI backend with Huey/SQLite task queue wrapping ISiCLE/NWChem for calculations. This architecture keeps dependencies minimal (no Redis, no database) while supporting the async job pattern essential for calculations that run minutes to hours.

The stack is well-established: FastAPI for the API layer, Huey with SQLite for job persistence without Redis complexity, RDKit for molecule handling and visualization, and ISiCLE/NWChem as the calculation engine. ISiCLE provides validated accuracy (MAE 0.33 ppm for 1H, 3.93 ppm for 13C with B3LYP/cc-pVDZ) but has API stability issues that require an abstraction layer.

The primary risks are NWChem resource management (memory exhaustion, scratch disk overflow) and ISiCLE API instability. These must be addressed in Phase 1 through proper memory/disk configuration and defensive wrapper code. The single-VM deployment constraint is manageable with strict concurrency control (1 worker initially) and proper job state persistence.

## Key Findings

### Recommended Stack

The stack prioritizes simplicity for single-VM deployment while maintaining clear upgrade paths. Huey with SQLite is the key architectural choice — it provides job persistence, status tracking, and retry support without requiring Redis or external services.

**Core technologies:**
- **FastAPI (>=0.128.0):** Modern async Python with automatic OpenAPI docs and Pydantic validation
- **Huey (>=2.5.0) with SqliteHuey:** Persistent job queue without Redis, survives restarts, supports priorities
- **ISiCLE (>=2.0.0) + NWChem (>=7.2.0):** Validated NMR calculation pipeline from PNNL
- **RDKit (>=2025.3.0):** Industry-standard molecule handling, supports atom annotation for shift labels
- **matplotlib (>=3.8.0):** Spectrum plotting, simple and sufficient for 1D NMR
- **aiosmtplib (>=5.0.0):** Async email notifications

### Expected Features

**Must have (table stakes):**
- SMILES string and MOL/SDF file input with validation
- Three calculation presets (Fast/Standard/Publication quality)
- Job submission returning immediate job ID
- Status polling endpoint with queued/running/complete/failed states
- JSON results with atom-assigned chemical shifts
- Optimized geometry download (XYZ/SDF)
- Raw NWChem output files download

**Should have (competitive):**
- Simulated spectrum plot (static PNG with Lorentzian peaks)
- Annotated structure image with shifts on atoms
- Email notification on completion
- Job history/list endpoint
- Molecule preview before submission

**Defer (v2+):**
- Interactive spectrum viewer (consider NMRium embed)
- Structure drawing widget (JSME or Ketcher)
- Webhook callbacks
- Batch submission API
- Additional nuclei (15N, 19F, 31P)
- User accounts and authentication

### Architecture Approach

The architecture follows the async job submission pattern with clear separation: FastAPI handles HTTP requests and immediately queues jobs to Huey, a separate Huey consumer process executes ISiCLE/NWChem calculations, and results are written to filesystem in job-specific directories. This separation ensures the API stays responsive while heavy computation happens in dedicated workers.

**Major components:**
1. **FastAPI Application:** HTTP interface for job submission, status queries, result retrieval, web UI
2. **Huey Task Queue (SQLite):** Persistent job queue with state tracking, survives restarts
3. **Huey Consumer (Worker):** Executes calculations, writes results, sends notifications
4. **ISiCLE Integration Layer:** Abstraction around ISiCLE calls for defensive error handling
5. **Filesystem Storage:** Job directories with input files, calculation outputs, metadata JSON
6. **Email Service:** Optional completion notifications via aiosmtplib

### Critical Pitfalls

1. **NWChem Memory Exhaustion** — Always set explicit memory (`memory total 4 gb` per core), use `semidirect` mode, implement pre-flight memory estimation based on molecule size

2. **NWChem Scratch Disk Overflow** — Use per-job scratch directories with auto-cleanup, enable `direct` SCF, monitor filesystem at 80% alert threshold

3. **FastAPI BackgroundTasks Unsuitable** — Use Huey queue instead; BackgroundTasks loses jobs on restart, has no persistence/retry. This is a foundational architecture decision for Phase 1

4. **SMILES Parsing Variations** — Validate with RDKit immediately on submission, canonicalize before processing, return clear error messages with specific parsing failures

5. **ISiCLE API Instability** — Pin version strictly, build abstraction layer with try/catch around all ISiCLE operations, comprehensive test suite before deployment

## Implications for Roadmap

Based on research, suggested phase structure:

### Phase 1: Foundation and ISiCLE Integration
**Rationale:** Everything depends on job management and calculation pipeline working correctly. NWChem configuration (memory, scratch) must be baked in from the start to avoid debugging failures later.
**Delivers:** Project structure, job directory management, Huey queue setup, ISiCLE wrapper with proper error handling, NWChem memory/scratch configuration
**Addresses:** Core infrastructure for all subsequent phases
**Avoids:** Memory exhaustion, scratch overflow, BackgroundTasks trap, zombie processes

### Phase 2: API Core
**Rationale:** Once calculation pipeline works, expose it via HTTP. Building API second allows testing end-to-end with real calculations.
**Delivers:** Job submission endpoint (POST /api/jobs), status endpoint (GET /api/jobs/{id}), input validation with RDKit, SMILES + file upload support
**Uses:** FastAPI, Pydantic, RDKit for validation
**Implements:** API layer, validation layer

### Phase 3: Results and Delivery
**Rationale:** Calculations complete; users need to retrieve results. Email notification reduces polling load.
**Delivers:** Results endpoint, file download endpoints, JSON shifts format, email notification service
**Uses:** aiosmtplib, filesystem storage patterns
**Implements:** Results delivery component

### Phase 4: Visualization
**Rationale:** Visual output differentiates from raw API. Static images are achievable without frontend complexity.
**Delivers:** Spectrum plot generation (matplotlib), annotated structure images (RDKit atomNote), downloadable PNG/SVG
**Uses:** matplotlib, RDKit drawing, scipy peak detection
**Implements:** Visualization component

### Phase 5: Web UI
**Rationale:** Browser interface after API is stable. Forms submit to existing endpoints.
**Delivers:** Job submission form, status page with auto-refresh, results display page, download links
**Uses:** Jinja2, FastAPI static files, Pico CSS or Tailwind
**Implements:** Web UI layer

### Phase 6: Polish and Presets
**Rationale:** Refinement after core flow works. Calculation presets map user-friendly options to DFT parameters.
**Delivers:** Fast/Standard/Publication presets, solvent selection, advanced parameter exposure, job history endpoint
**Implements:** Configuration layer, production hardening

### Phase Ordering Rationale

- **Foundation first:** ISiCLE/NWChem integration is the highest-risk, most complex piece. Discovering integration issues early avoids rework. Memory and scratch configuration must be correct from the start.
- **API before UI:** API is simpler to test and debug. UI depends on stable API contracts.
- **Visualization separate:** Independent of API correctness; can be parallelized or deferred if timeline is tight.
- **Presets last:** Polish layer that maps to existing functionality, low risk.

### Research Flags

Phases likely needing deeper research during planning:
- **Phase 1:** ISiCLE integration details, NWChem output parsing format, conformer handling — ISiCLE-specific documentation is sparse
- **Phase 4:** RDKit atom annotation specifics, spectrum Lorentzian broadening parameters

Phases with standard patterns (skip research-phase):
- **Phase 2:** FastAPI REST patterns well-documented
- **Phase 3:** File serving, JSON responses are standard
- **Phase 5:** Jinja2 templating is straightforward

## Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| Stack | HIGH | Verified via PyPI, official docs; all components actively maintained |
| Features | MEDIUM | Based on competitive analysis; feature expectations inferred from similar tools |
| Architecture | HIGH | Async job pattern well-established; Huey/SQLite verified in documentation |
| Pitfalls | MEDIUM | ISiCLE-specific pitfalls have sparse documentation; NWChem pitfalls well-documented |

**Overall confidence:** MEDIUM-HIGH

### Gaps to Address

- **ISiCLE conformer handling:** CREST failures reported but unclear how ISiCLE handles them; need integration testing with various molecule sizes during Phase 1
- **NWChem output parsing:** Multiple parser implementations exist (DP5, OpenBabel) but exact format for NMR shieldings needs verification during implementation
- **Calculation time estimation:** No clear formula found; will need empirical calibration based on atom count and basis set
- **ISiCLE API surface:** Documentation incomplete; need exploratory testing to understand actual API behavior vs documented behavior

## Sources

### Primary (HIGH confidence)
- [FastAPI Documentation](https://fastapi.tiangolo.com/) — API patterns, background tasks limitations
- [Huey Documentation](https://huey.readthedocs.io/) — SQLite backend configuration, worker setup
- [NWChem FAQ](https://nwchemgit.github.io/FAQ.html) — Memory configuration, AUTOZ issues
- [NWChem Memory Wiki](https://github.com/nwchemgit/nwchem/wiki/Memory) — Memory allocation model
- [RDKit Documentation](https://www.rdkit.org/docs/) — Molecule parsing, atom annotation

### Secondary (MEDIUM confidence)
- [ISiCLE Paper (J. Cheminform. 2018)](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-018-0305-8) — Accuracy benchmarks, workflow overview
- [ISiCLE GitHub](https://github.com/pnnl/isicle) — API patterns, known issues
- [Azure Async Request-Reply Pattern](https://learn.microsoft.com/en-us/azure/architecture/patterns/async-request-reply) — Job polling design
- [NMRium Features](https://www.nmrium.com/features) — UI/UX patterns for NMR tools

### Tertiary (LOW confidence)
- [CREST GitHub Issue #203](https://github.com/crest-lab/crest/issues/203) — Conformer SIGSEGV, workarounds unclear
- [ISiCLE GitHub Issues](https://github.com/pnnl/isicle/issues) — API bugs, needs direct testing to verify current state

---
*Research completed: 2026-01-19*
*Ready for roadmap: yes*
