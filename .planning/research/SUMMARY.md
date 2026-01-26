# Research Summary: Conformational Sampling for NMR (v2.0)

**Project:** qm-nmr-calc v2.0 Conformational Sampling
**Domain:** Computational Chemistry / NMR Prediction
**Researched:** 2026-01-26
**Confidence:** HIGH

## Executive Summary

This milestone adds Boltzmann-weighted ensemble NMR averaging to qm-nmr-calc's existing single-conformer prediction service. Flexible molecules adopt multiple low-energy conformations in solution; accurate NMR prediction requires averaging shifts across the thermally populated conformer ensemble. Research shows this is standard practice for computational NMR in both academic (ISiCLE, CENSO) and commercial (Spartan) tools.

The recommended approach uses RDKit KDG conformer generation as an always-available baseline, with optional CREST/xTB for high-accuracy mode when binaries are detected. No new Python dependencies are needed—scipy and RDKit are already in the stack. The architecture fits cleanly into the existing Huey task queue: one task per ensemble, sequential conformer processing (matching the single-worker environment), per-conformer scratch directories to avoid NWChem file conflicts. The critical implementation decisions are already validated: DFT energies for Boltzmann weighting (not MMFF/xTB), KDG for solution-phase sampling (not ETKDG's crystal bias), and weighted-average-only output (no per-conformer detail, cleaner API).

The main risks are numerical (Boltzmann overflow with wide energy ranges), operational (CREST hanging on macrocycles), and data integrity (atom ordering consistency across conformers). All have proven mitigation strategies from the literature: exp-normalize trick, timeouts with fallback, and canonical atom tracking. The most important architectural decision is scratch directory isolation per conformer—violating this causes NWChem database corruption when multiple optimizations run. Implementation should prioritize the sequential single-task path in Phase 1, deferring parallel optimization and CREST integration to later phases.

## Key Findings

### Recommended Stack

**Status: Zero new Python dependencies required.** RDKit (>=2025.09.3) and scipy (>=1.17.0) are already in the project stack and provide everything needed for conformer generation and Boltzmann statistics. CREST/xTB are optional system binaries that enable high-accuracy mode when detected on PATH.

**Core technologies:**
- **RDKit KDG:** Solution-phase conformer generation without crystal bias (ETKDG unsuitable for solution NMR) — already available, use KDG method explicitly
- **scipy/numpy:** Boltzmann weight calculation using custom exp-normalize implementation (scipy.stats.boltzmann is wrong distribution) — already available
- **CREST + xTB binaries:** Optional metadynamics-based conformer generation with better DFT correlation (ρ ~ 0.4 vs MMFF ρ ~ -0.1) — auto-detect via `shutil.which('crest')`
- **NWChem:** DFT optimization and NMR for each conformer — existing integration works as-is, just call N times

**Critical version notes:**
- RDKit: Use KDG not ETKDG (ETKDG has CSD crystal structure bias unsuitable for solution-phase NMR)
- CREST: 3.0+ integrates tblite, reducing xTB dependency for some features
- xTB: Required by CREST for QCG and other features despite tblite integration

**What NOT to use:**
- MMFF energies for Boltzmann weighting (anticorrelate with DFT, ρ ~ -0.1 to -0.45) — use for pre-filtering only
- xTB energies for final weighting (better than MMFF but not DFT-accurate) — use for pre-filtering only
- scipy.stats.boltzmann (wrong distribution, PMF not conformer weighting) — implement custom numpy version
- CENSO workflow (excellent but requires ORCA/Turbomole) — incompatible with NWChem backend

### Expected Features

**Must have (table stakes):**
- Conformer generation method choice (RDKit vs CREST) — users need speed/accuracy tradeoff
- Boltzmann-weighted average shifts — core feature, non-negotiable
- Energy window filtering (pre-DFT and post-DFT) — reduces compute without losing accuracy
- Temperature parameter — experimental condition matching essential
- Single-conformer mode (backward compatibility) — don't force ensemble on rigid molecules
- Per-atom shift assignments — existing v1.x feature, ensemble returns weighted average per atom
- Lowest-energy geometry visualization — standard output for structure viewing

**Should have (competitive advantages):**
- Automatic CREST detection with RDKit fallback — zero-config high-accuracy mode, better UX than manual backend selection
- Weighted-average-only API (no per-conformer detail dump) — simpler, cleaner API focused on experimental comparison (Spartan/CENSO show detail, we don't)
- Population metadata (conformer count, energy range, top 3 populations) — transparency without overwhelming detail
- Two-stage energy filtering (wide pre-DFT, tight post-DFT) — saves DFT compute on MMFF/xTB artifacts
- Solvent-specific DELTA50 scaling factors — existing v1.1 advantage carries forward to ensemble

**Defer (v2+):**
- Per-conformer detail export (geometries, individual shifts) — niche advanced use case, adds API complexity
- DP4+ probability scoring — requires experimental spectrum input, scope change for v3.0
- Real-time progress streaming (WebSocket) — marginal UX improvement over polling, adds infrastructure complexity
- Multi-temperature ensemble — variable-temperature NMR experiments (niche)
- Automatic rigidity heuristics — user explicit choice sufficient (rotatable bond count is crude, nConf20 adds complexity)
- Conformer caching — optimization for parameter sweeps, defer to v2.2

**Anti-features (explicitly avoid):**
- Automatic rigidity detection forcing ensemble/single mode — misclassification wastes compute or gives poor results, user knows their molecule better
- Per-conformer shift detail in default API — overwhelming data (20 conf × 30 atoms = 600 values), rarely actionable
- MMFF energy Boltzmann weighting — anticorrelates with DFT, discards important conformers
- Automatic method selection by molecule size — size poor proxy for flexibility (18-atom cyclohexane needs ensemble, 30-atom rigid steroid doesn't)

### Architecture Approach

**One Huey task per ensemble, sequential conformer processing.** Current single-worker environment (one VM, one Huey process) makes parallel task spawning pointless. The entire conformer set is processed in a single task with sequential DFT optimization and NMR calculation loops. This matches the existing v1.x pattern (single task per job) and simplifies orchestration, progress tracking, and error handling.

**Per-conformer scratch directories prevent NWChem conflicts.** Critical architectural requirement: each conformer's NWChem calculation must use isolated scratch directory (`scratch/conf_{id:04d}/`) with unique start name to avoid database file collisions. Current v1.x code reuses `scratch/` which works for single-conformer but causes "database corrupted" errors in ensemble mode.

**Major components:**
1. **conformers.py** (new) — Conformer generation with `generate_rdkit_conformers()` and `generate_crest_conformers()`, auto-detection via `has_crest()`, adaptive conformer count by rotatable bond count
2. **averaging.py** (new) — Boltzmann weighting with `boltzmann_weights()` using exp-normalize trick, `average_shieldings()` for weighted ensemble average
3. **ensemble.py** (new) — Orchestration helpers: energy filtering (pre/post DFT), conformer data models, progress tracking updates
4. **tasks.py** (extend) — Add `run_ensemble_nmr_task()` alongside existing `run_nmr_task()`, route based on `conformer_mode` field
5. **models.py** (extend) — Add `ConformerData`, `EnsembleData`, extend `JobInput` and `JobStatus` with ensemble fields (all optional, backward compatible)
6. **nwchem/runner.py** (minor change) — Accept `conformer_id` parameter for unique scratch directory and start name generation

**Data flow:**
Conformer generation (RDKit/CREST) → Pre-DFT energy filter (6 kcal/mol) → DFT optimization loop (parse energies) → Post-DFT filter (3 kcal/mol or 95% cumulative) → NMR shielding loop → Boltzmann average → Apply DELTA50 scaling → Visualize lowest-energy conformer

**Job directory structure changes:**
- Add `output/conformers/` for initial geometries
- Add `output/optimized/` for DFT-optimized geometries
- Add `output/conformer_ensemble.json` for full ensemble metadata
- Change `scratch/` to `scratch/conf_{id:04d}/` subdirectories per conformer
- `output/optimized.xyz` remains lowest-energy conformer for visualization (backward compatible)

### Critical Pitfalls

1. **Boltzmann numerical overflow/underflow** — Energy differences >10 kcal/mol cause `exp(-ΔE/RT)` underflow to 0.0, breaking normalization. **Solution:** Exp-normalize trick (shift energies by minimum before exponentiation, ensuring max exponent is 0).

2. **Wrong energy units for Boltzmann weighting** — CREST outputs Hartrees, NWChem outputs Hartrees, Boltzmann expects kcal/mol. Forgetting conversion gives nonsensical weights (off by factor of 627). **Solution:** Convert to kcal/mol at parsing stage (HARTREE_TO_KCAL = 627.5095), enforce unit tracking in conformer data structures.

3. **CREST hanging/timeout on macrocycles** — CREST enters infinite loop on >12-membered rings or highly flexible molecules, no progress output, eventually fills disk. **Solution:** Implement timeout (default 3600s) with automatic fallback to RDKit, set `GFORTRAN_UNBUFFERED_ALL=1` for macOS compatibility.

4. **Atom ordering inconsistency across conformers** — Atom indices change between CREST generation, NWChem optimization, and NMR parsing. Averaging shifts for "atom 5" mixes carbon in one conformer with hydrogen in another. **Solution:** Establish canonical atom order from SMILES, verify atom sequence consistency across conformers, use canonical indexing for averaging.

5. **NWChem scratch directory conflicts** — Running multiple conformers in parallel causes database file collisions (`molecule.db` overwritten), producing "database corrupted" errors or wrong results. **Solution:** Per-conformer scratch directories (`scratch/conf_{id:04d}/`) with unique NWChem start names (`start conf_0` instead of `start molecule`).

6. **Using MMFF/xTB energies for Boltzmann weighting** — Force field energies anticorrelate with DFT (ρ ~ -0.1 to -0.45), causing wrong conformer populations and averaged shifts. **Solution:** Always use DFT energies from NWChem optimization for Boltzmann weighting; use MMFF/xTB only for pre-DFT filtering.

## Implications for Roadmap

Based on research, suggested phase structure prioritizes risk mitigation and incremental validation:

### Phase 1: Conformer Data Model and Storage
**Rationale:** Establish foundation before generation/calculation logic. Canonical atom ordering must be designed correctly from the start—retrofitting is expensive and error-prone (Pitfall 4).

**Delivers:** `ConformerData`, `EnsembleData` models with unit tracking, per-conformer job directory structure, status tracking schema

**Addresses:** Atom ordering consistency (Pitfall 4), energy unit tracking (Pitfall 2), backward compatibility (existing v1.x jobs load with defaults)

**Avoids:** Retrofitting data structures after implementation discovers ordering inconsistencies

**Research needs:** Standard phase (established patterns), skip `/gsd:research-phase`

---

### Phase 2: RDKit KDG Conformer Generation
**Rationale:** RDKit-only path is simpler (no external binaries), validates end-to-end workflow before adding CREST complexity. KDG method critical for solution-phase NMR (not ETKDG).

**Delivers:** `generate_rdkit_conformers()` with KDG method, MMFF optimization, RMSD deduplication, adaptive conformer count by rotatable bonds

**Addresses:** Solution-phase sampling (KDG not ETKDG per research), conformer deduplication (Technical Debt Pattern 2), insufficient sampling (Pattern 1)

**Uses:** RDKit (existing), MMFF for pre-filtering energies only

**Avoids:** Crystal structure bias from ETKDG (Pattern 3)

**Research needs:** Standard phase (RDKit well-documented), skip `/gsd:research-phase`

---

### Phase 3: Boltzmann Averaging Implementation
**Rationale:** Core algorithm must be numerically stable before integrating with expensive DFT calculations. Unit tests with wide energy ranges (0-20 kcal/mol) catch overflow issues early.

**Delivers:** `boltzmann_weights()` with exp-normalize trick, `average_shieldings()` for weighted ensemble, unit tests with known test cases

**Addresses:** Numerical overflow (Pitfall 1), energy unit consistency (Pitfall 2), single-conformer degeneracy (Missing Feature 3)

**Uses:** numpy (existing), custom implementation (not scipy.stats.boltzmann per research)

**Avoids:** Underflow to 0.0 breaking normalization, wrong scipy distribution usage

**Research needs:** Standard phase (statistical mechanics textbook material), skip `/gsd:research-phase`

---

### Phase 4: Multi-Conformer NWChem Integration
**Rationale:** Most complex phase—parallel DFT jobs with scratch isolation, energy parsing, partial failure handling. Sequential processing validated before attempting parallelization.

**Delivers:** `run_ensemble_nmr_task()`, per-conformer scratch directories, DFT energy extraction, two-stage energy filtering (pre/post DFT)

**Addresses:** Scratch directory conflicts (Pitfall 5), wrong energy source (Pitfall 6), partial failures (Missing Feature 1), disk space explosion (Trap 2)

**Uses:** NWChem (existing integration), unique scratch dirs per conformer, DFT energies only for Boltzmann weighting

**Implements:** Sequential conformer optimization loop (matches single-worker environment), per-conformer progress updates

**Avoids:** Database corruption from shared scratch, MMFF energies for final weighting

**Research needs:** Phase-specific research likely needed for NWChem output parsing edge cases (different functionals, convergence failures)

---

### Phase 5: CREST Integration (Optional High-Accuracy Mode)
**Rationale:** Deferred until RDKit path fully validated. CREST adds operational complexity (hanging, timeouts, memory issues) that shouldn't block core functionality.

**Delivers:** `generate_crest_conformers()`, CREST/xTB auto-detection, timeout with RDKit fallback, ALPB solvation model support

**Addresses:** CREST hanging (Pitfall 3), xTB memory explosion (Gotcha 1), ALPB vs GBSA confusion (Gotcha 2)

**Uses:** CREST + xTB binaries (optional), subprocess timeout, environment variable setup (OMP_STACKSIZE, GFORTRAN_UNBUFFERED_ALL)

**Implements:** Graceful fallback to RDKit on timeout/failure, molecule size-based method recommendation (<200 atoms)

**Avoids:** Mandatory CREST dependency blocking app functionality, infinite loops on macrocycles

**Research needs:** Phase-specific research needed for CREST command-line options, output parsing edge cases, solvent model compatibility

---

### Phase 6: API Integration and Progress Tracking
**Rationale:** User-facing layer after backend validated. Progress tracking essential for long-running ensemble jobs (hours).

**Delivers:** Ensemble fields in API schemas, job submission routing, fine-grained progress updates, time estimation

**Addresses:** Progress reporting (Missing Feature 2), timeout handling (Trap 3), conformer mode selection

**Uses:** FastAPI (existing), status.json updates after each conformer, estimated completion time calculation

**Implements:** `conformer_mode` field (single/ensemble), `conformer_method` field (rdkit/crest), energy window parameters, population metadata in response

**Avoids:** Users seeing "Running" for hours with no updates, duplicate job submissions due to perceived timeouts

**Research needs:** Standard phase (REST API patterns), skip `/gsd:research-phase`

---

### Phase Ordering Rationale

**Dependencies drive order:**
- Data models (Phase 1) before generation logic (Phase 2) — canonical atom ordering must be designed correctly upfront
- Boltzmann averaging (Phase 3) before DFT integration (Phase 4) — validate numerics with cheap test data before expensive calculations
- RDKit path (Phases 2-4) before CREST (Phase 5) — simpler path proves architecture before adding operational complexity
- Backend working (Phases 1-5) before API (Phase 6) — validate core functionality before exposing to users

**Risk mitigation:**
- Most critical pitfalls (numerical overflow, atom ordering, scratch conflicts) addressed in early phases
- CREST operational issues (hanging, memory) isolated in Phase 5, can be deferred if blocked
- Expensive DFT integration (Phase 4) happens after cheap validation (Phases 2-3)

**Incremental validation:**
- Each phase delivers testable functionality
- RDKit-only path (Phases 1-4, 6) is complete feature without CREST
- CREST (Phase 5) is pure enhancement, not blocking

### Research Flags

**Phases needing deeper research during planning:**
- **Phase 4 (Multi-Conformer NWChem):** NWChem output parsing edge cases (different functionals, convergence failures, error messages), scratch directory cleanup strategies, disk space monitoring
- **Phase 5 (CREST Integration):** CREST command-line options for edge cases (charged molecules, radicals), output format variations, timeout recovery strategies, xTB environment configuration for different OS/hardware

**Phases with standard patterns (skip research-phase):**
- **Phase 1:** Data modeling (established Pydantic patterns)
- **Phase 2:** RDKit conformer generation (well-documented, existing examples)
- **Phase 3:** Boltzmann statistics (textbook implementation)
- **Phase 6:** REST API patterns (existing FastAPI stack)

## Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| Stack | HIGH | All components verified in official docs and existing project dependencies; no new Python packages needed |
| Features | HIGH | Validated against ISiCLE, CENSO, Spartan implementations; clear table stakes vs differentiators |
| Architecture | HIGH | Single-task sequential pattern matches existing v1.x architecture and single-worker environment |
| Pitfalls | HIGH | All pitfalls sourced from GitHub issues, documented failures, and computational chemistry literature |

**Overall confidence:** HIGH

### Gaps to Address

**Numerical stability validation:** Unit tests needed for Boltzmann weighting with extreme energy ranges (0-20 kcal/mol, 0-40 kcal/mol) to verify exp-normalize trick handles edge cases. Validate against known test cases (e.g., 2 conformers at 0 and 1 kcal/mol → weights 0.843 and 0.157 at 298K).

**CREST timeout tuning:** Default 3600s timeout may need adjustment based on actual molecule complexity in production. Monitor CREST runtime during beta testing to set appropriate timeout (1 hour for typical drugs, may need 2-3 hours for complex macrocycles).

**Atom ordering verification:** Implement consistency checks during development, but real-world validation needs diverse test set (different SMILES formats, explicit hydrogens vs implicit, charged species). Add assertion that fails loudly if atom sequence differs between conformers.

**Disk space monitoring:** Current code doesn't track disk usage. For 50-conformer ensembles with 500 MB scratch per conformer, need monitoring to prevent job failure mid-run. Add pre-flight check (require 10 GB free) and cleanup strategy (delete scratch files after parsing results).

**Partial failure UX:** Research indicates partial failures common (5-10% of conformers fail DFT optimization). Need to define threshold: continue if >90% succeed? Warn user if <80%? Fail if <50%? Requires UX decision during Phase 4 implementation.

**CREST solvent model coverage:** ALPB supports limited solvent list. If user requests ensemble with unsupported solvent, should job: (a) fail, (b) run CREST in vacuum and warn, or (c) fall back to RDKit? Defer decision to Phase 5 implementation, likely option (c) for best UX.

## Sources

### Primary (HIGH confidence)

**Official Documentation:**
- RDKit rdDistGeom module: KDG method documentation and parameters
- RDKit AllChem module: EmbedMultipleConfs, MMFFOptimizeMoleculeConfs APIs
- CREST Documentation: Command-line keywords, installation, file formats
- CREST File Formats: Multi-XYZ structure, energy comment line format
- xTB Setup and Installation: Binary installation, environment variables
- xTB GBSA Documentation: ALPB vs GBSA solvation models, parameterization
- NWChem Wiki Scratch_Dir: Scratch directory guidelines, start name conflicts

**Repository Releases:**
- CREST GitHub Releases: Version 3.0.2 stable, continuous pre-release
- xTB GitHub Releases: Version 6.7.1 stable
- RDKit Releases: 2025.09.1 allene/cumulene support discussion

**Scientific Literature:**
- Riniker & Landrum (2015) J. Chem. Inf. Model.: KDG method development, crystal bias in ETKDG
- Pracht et al. (2024) J. Chem. Phys.: CREST program description and validation
- Bannwarth et al. (2019) J. Chem. Theory Comput.: GFN2-xTB method parameterization
- Yesiltepe et al. (2018) J. Cheminform.: ISiCLE NWChem-based NMR with conformer averaging

**Existing Codebase:**
- `/home/chris/develop/qm-nmr-calc/src/qm_nmr_calc/tasks.py`: Current Huey task structure
- `/home/chris/develop/qm-nmr-calc/src/qm_nmr_calc/storage.py`: Job directory management patterns
- `/home/chris/develop/qm-nmr-calc/src/qm_nmr_calc/nwchem/runner.py`: NWChem execution and output parsing
- `/home/chris/develop/qm-nmr-calc/pyproject.toml`: Existing dependencies (RDKit >=2025.09.3, scipy >=1.17.0)

### Secondary (MEDIUM confidence)

**Community Resources:**
- RDKit Discussion #8226: Best practices for conformer generation (KDG vs ETKDG trade-offs)
- BLOPIG Blog: Advances in conformer generation (ETKDG crystal bias analysis)
- Corin Wagen Blog: Evaluating Boltzmann weighting error (exp-normalize trick, numerical stability)
- Tim Vieira Blog: Exp-normalize trick derivation (log-sum-exp for numerical stability)

**GitHub Issues (Documented Failures):**
- CREST Issue #242, #105, #313, #382: Cregen elimination, optimization failures, macOS hanging
- xTB Issue #700, #439, Discussion #811: Memory explosion on large molecules, stack overflow
- RDKit Issue #1504, #5250, #8001, Discussion #5368: Stereochemistry loss, generation failures, RMSD pruning

**Methodological Reviews:**
- Nature Reviews Chemistry: NMR ensemble determination by data deconvolution
- J. Chem. Inf. Model.: nConf20 descriptor for conformational flexibility quantification
- J. Chem. Phys.: GEOM energy-annotated conformer dataset (large-scale validation)
- PMC: Conformational search importance for NMR (force field vs DFT correlation analysis)

### Tertiary (LOW confidence, needs validation)

**Performance estimates:** DFT timing (5-15 min optimization, 10-30 min NMR for 30-atom molecule) based on typical NWChem performance, not measured for this specific stack. Validate during Phase 4 implementation with benchmark molecules.

**Conformer count recommendations:** Guidelines (5-10 for rigid, 50-100 for flexible) from literature, not validated for RDKit KDG with this specific DFT level. May need tuning based on actual conformer diversity observed during beta testing.

**CREST timeout values:** Default 3600s timeout suggested based on community reports ("30-60 min typical"), not measured for project's target molecule classes. Monitor during Phase 5 and adjust based on actual data.

---

*Research completed: 2026-01-26*
*Ready for roadmap: YES*
*Recommended starting phases: Phase 1 (Data Model), Phase 2 (RDKit), Phase 3 (Boltzmann)*
