# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-01-26)

**Core value:** Reliable async NMR predictions with full control over calculation parameters -- submit a molecule, get back accurate 1H/13C shifts without babysitting long-running calculations.
**Current focus:** v2.0 Conformational Sampling -- Phase 15 in progress

## Current Position

Phase: 15 of 17 (Multi-Conformer NWChem Integration)
Plan: 03 of 03 (Phase 15)
Status: Phase complete
Last activity: 2026-01-27 -- Completed 15-03-PLAN.md (NMR loop and full ensemble orchestrator)

Progress: █████████████████████████████████████████████ 100% (47/47 plans complete)

## Performance Metrics

**Velocity:**
- Total plans completed: 47
- Average duration: 8.1 min
- Total execution time: 390 min (6.5 hours)

**By Milestone:**

| Milestone | Phases | Plans | Duration |
|-----------|--------|-------|----------|
| v1.0 Core NMR Service | 6 | 16 | 2 days |
| v1.1 Accurate Chemical Shifts | 8 | 21 | 5 days |
| v2.0 Conformational Sampling | 6 | TBD | In progress |

**Recent Trend:**
- Last 5 plans: 6-10 min
- Trend: Consistent fast TDD (14-01: 6 min, 14-02: 7 min, 15-01: 9 min, 15-02: 10 min, 15-03: 9 min)

## Accumulated Context

### Decisions

All decisions logged in PROJECT.md Key Decisions table.
Recent v2.0 decisions affecting current work:

- KDG over ETKDG: Avoids crystal structure bias for solution-phase NMR
- CREST/xTB optional: App works without them, enables high-accuracy mode when detected
- Boltzmann weight by DFT energies: Most accurate readily available energy level for weighting
- RDKit CanonicalRankAtoms for atom ordering: Built-in deterministic canonicalization (12-02)
- Mapping as metadata not transformation: Don't physically reorder atoms (breaks RDKit ops) (12-02)
- Per-conformer scratch directory isolation: Each conformer gets unique scratch dir to prevent NWChem database conflicts (12-01)
- Explicit energy unit tracking: ConformerData stores energy_unit alongside energy to prevent conversion bugs (12-01)
- str type for API conformer_mode: API schemas use str (not Literal) for flexibility, validation at service layer (12-03)
- KDG params for solution-phase: pruneRmsThresh=-1.0 (no pre-opt pruning), numThreads=0 (all threads), random coords fallback (13-01)
- 8 rotatable bonds threshold: Separates rigid (50 confs) from flexible (200 confs) based on conformational space needs (13-01)
- MMFF native units: Return energies in kcal/mol without conversion (13-01)
- Dedup before energy filter: Run RMSD deduplication first, then energy window filter on deduped subset (13-03)
- Sequential 1-based conformer IDs: conf_001, conf_002, ... for user-facing consistency (13-03)
- Relative geometry_file paths: Stored relative to job dir for portability (13-03)
- exp-normalize trick for Boltzmann: Subtract min energy before exp to prevent overflow/underflow (14-01)
- Pure Python math.exp: No numpy dependency for Boltzmann weights, sufficient performance for typical conformer counts (14-01)
- Weighted averaging by atom index: Cross-conformer matching via index field (NWChem 1-based), not implicit ordering (14-02)
- Descending shift sort: NMR convention (highest shift first) for all returned results (14-02)
- In-place weight population: average_ensemble_nmr mutates ConformerData.weight for tracking (14-02)
- re.findall last-occurrence for optimization: Takes final energy from multi-cycle NWChem optimization (15-01)
- Hartree native units for DFT energy: No conversion from NWChem output, matches Boltzmann functions (15-01)
- Sequential conformer processing: Process conformers one at a time for simplicity, parallelize later if needed (15-02)
- 50% success threshold: Fail if more than half of conformers fail DFT optimization (15-02)
- smiles=None support: Allow geometry-file-only input for conformer optimization path (15-02)
- No NMR minimum success threshold: Any successful NMR results usable since DFT already caught systematic failures (15-03)
- nmr_basis_set substitution: NMR step uses nmr_basis_set instead of optimization basis_set for higher accuracy (15-03)
- In-place ensemble mutation: Conformer status fields updated by reference during pipeline (15-03)

### Roadmap Evolution

- v1.0: 6 phases (1-6), shipped 2026-01-20
- v1.1: 8 phases (7-11.2, including 3 inserted), shipped 2026-01-25
- v2.0: 6 phases (12-17), started 2026-01-26
  - Phase structure: Data Model → RDKit → Boltzmann → NWChem → CREST → API
  - Risk mitigation: RDKit-only path (12-15, 17) complete before CREST complexity (16)

### Pending Todos

None.

### Blockers/Concerns

**v2.0 Architecture:**
- CREST timeouts: Macrocycles can hang, need timeout with RDKit fallback (Phase 16)

**Resolved (v2.0):**
- Numerical stability: Boltzmann weighting implemented with exp-normalize trick, tested with extreme ranges (14-01)
- Atom ordering consistency: Canonical indexing established via CanonicalRankAtoms (12-02)
- Scratch directory isolation: Per-conformer directories prevent NWChem database corruption (12-01)
- Backward compatibility: v1.x JSON loads verified with 8 fixture tests (12-03)

**Resolved (v1.x):**
- RDKit stderr capture limitation (known, fallback messages used)
- Single-conformer limitation (being addressed in v2.0)

## Session Continuity

Last session: 2026-01-27
Stopped at: Completed 15-03-PLAN.md (NMR loop and full ensemble orchestrator) -- Phase 15 complete
Resume file: None
Next: Phase 16 (CREST Integration) - High-accuracy conformer generation with timeout handling
