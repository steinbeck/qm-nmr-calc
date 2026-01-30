---
phase: 24-conformer-preselection
plan: 01
subsystem: conformers
tags: [rdkit, clustering, butina, rmsd, conformer-selection]

# Dependency graph
requires:
  - phase: 12-conformer-generation
    provides: RDKit conformer generation and MMFF optimization
provides:
  - RMSD-based clustering to reduce conformer count before DFT
  - Butina algorithm implementation with symmetry-aware RMSD
  - Cluster representative selection by energy
affects: [25-conformer-integration, dft-optimization, ensemble-workflow]

# Tech tracking
tech-stack:
  added: [rdkit.ML.Cluster.Butina, rdkit.Chem.rdMolAlign]
  patterns: [clustering-for-preselection, symmetry-aware-rmsd]

key-files:
  created:
    - src/qm_nmr_calc/conformers/clustering.py
    - tests/test_clustering.py
  modified: []

key-decisions:
  - "Butina algorithm for RMSD clustering (handles symmetry via GetBestRMS)"
  - "Default threshold 1.5 Angstrom (aggressive reduction, balanced diversity)"
  - "Lowest-energy conformer per cluster as representative"
  - "Configurable max_conformers parameter (default 8)"

patterns-established:
  - "Symmetry-aware RMSD: Use GetBestRMS instead of basic RMSD for correct symmetry handling"
  - "Cluster-then-select pattern: Cluster first by structure, then select by energy"

# Metrics
duration: 3 min
completed: 2026-01-30
---

# Phase 24 Plan 01: RMSD Clustering Summary

**Butina RMSD clustering with symmetry-aware distance calculation reduces 40+ conformers to ~8 diverse representatives for DFT optimization**

## Performance

- **Duration:** 3 min
- **Started:** 2026-01-30T08:43:55Z
- **Completed:** 2026-01-30T08:46:43Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments

- Created clustering.py module with Butina algorithm implementation
- Implemented symmetry-aware RMSD via GetBestRMS for correct molecular symmetry handling
- Cluster representative selection picks lowest-energy conformer from each cluster
- Comprehensive test suite (11 tests) covering edge cases and full pipeline
- Verified reduction: hexane 30→2 conformers, heptane 30→8-12 clusters

## Task Commits

Each task was committed atomically:

1. **Task 1: Create clustering.py module** - `1107bd7` (feat)
2. **Task 2: Create unit tests for clustering** - `dffb026` (test)

**Plan metadata:** (pending final commit)

## Files Created/Modified

- `src/qm_nmr_calc/conformers/clustering.py` - RMSD clustering functions using Butina algorithm
- `tests/test_clustering.py` - Comprehensive tests with flexible/rigid molecule fixtures

## Decisions Made

**Butina algorithm for clustering:**
- Standard algorithm for molecular clustering in cheminformatics
- Handles arbitrary distance matrices, works well with RMSD
- Returns clusters sorted by size (largest first)

**Symmetry-aware RMSD (GetBestRMS):**
- Critical for molecules with symmetry (e.g., benzene, substituted rings)
- GetBestRMS considers all symmetry-equivalent atom mappings
- Basic RMSD would overestimate distances for symmetric molecules

**Default RMSD threshold 1.5 Angstrom:**
- Research (CONFORMER_PRESELECTION.md) shows 1.5Å reduces 40 conformers to 8-12 clusters
- Aggressive reduction suitable for expensive DFT calculations
- Configurable parameter allows tuning per molecule

**Lowest-energy representative per cluster:**
- Each cluster represents a conformational basin
- Selecting lowest-energy ensures best representative of that basin
- Alternative (centroid) would require additional RMSD calculations

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

**Ready for integration (Phase 25):**
- Clustering module complete and tested
- Needs integration into conformer generation workflow
- Should be called after MMFF optimization, before DFT

**Implementation path for Phase 25:**
1. Add clustering step to generator.py or create new preselection.py module
2. Extract MMFF energies from optimization results
3. Call cluster_and_select with mol and energies
4. Filter conformers to selected representatives before DFT
5. Update ensemble workflow to use clustered conformers

**Performance expectations:**
- Clustering is fast (RMSD calculation scales O(n²), n=30-50 is negligible)
- Should reduce DFT workload from 40-50 conformers to 8-12
- Maintains conformational diversity while eliminating redundancy

---
*Phase: 24-conformer-preselection*
*Completed: 2026-01-30*
