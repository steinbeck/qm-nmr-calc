---
phase: 12-conformer-data-model
plan: "01"
subsystem: data-models
status: complete
tags: [pydantic, data-models, storage, conformers, v2.0-foundation]
dependencies:
  requires: []
  provides:
    - conformer-data-models
    - per-conformer-storage
    - backward-compatible-jobstatus
  affects:
    - 13-rdkit-conformer-generation
    - 14-boltzmann-averaging
    - 15-nwchem-conformer-loop
    - 16-crest-integration
    - 17-api-conformer-support
tech-stack:
  added: []
  patterns:
    - per-resource-directory-isolation
    - optional-fields-for-backward-compatibility
key-files:
  created:
    - src/qm_nmr_calc/models.py (ConformerData, ConformerEnsemble, EnergyUnit)
    - tests/test_conformer_models.py
  modified:
    - src/qm_nmr_calc/models.py (JobStatus, JobInput extended)
    - src/qm_nmr_calc/storage.py (per-conformer directory helpers)
decisions:
  - id: CONF-001
    title: Per-conformer scratch directory isolation
    choice: Each conformer gets unique scratch/conformers/{conf_id}/ directory
    rationale: Prevents NWChem database file conflicts when running parallel conformer optimizations
    alternatives: ["Shared scratch with unique filenames", "Sequential processing only"]
    impact: Eliminates Pitfall 5 from research (database corruption)
  - id: CONF-002
    title: Backward compatibility via Optional defaults
    choice: All new JobStatus/JobInput fields use Optional with defaults
    rationale: v1.x status.json files must load without modification
    alternatives: ["Version migrations", "Separate v2 schema"]
    impact: Zero breaking changes, smooth v1→v2 transition
  - id: CONF-003
    title: Energy unit explicit tracking
    choice: ConformerData stores energy_unit alongside energy value
    rationale: Different stages use different units (RDKit kcal/mol, NWChem Hartree)
    alternatives: ["Always convert on read", "Implicit unit assumption"]
    impact: Prevents unit conversion bugs, explicit is better than implicit
metrics:
  tests-added: 20
  tests-passing: 124
  duration: "4.3 min"
  completed: "2026-01-26"
---

# Phase 12 Plan 01: Conformer Data Model and Storage Summary

**One-liner:** Pydantic models for conformer ensembles (ConformerData, ConformerEnsemble) with per-conformer scratch directory isolation and full v1.x backward compatibility.

## What Was Built

### Data Models (src/qm_nmr_calc/models.py)

**EnergyUnit type:**
- Literal["hartree", "kcal_mol", "kj_mol"]
- Explicit tracking prevents unit conversion errors

**ConformerData model:**
- `conformer_id`: Unique identifier (e.g., "conf_001")
- `energy`: Optional[float] = None (populated after DFT)
- `energy_unit`: Optional[EnergyUnit] = None
- `weight`: Optional[float] = None (Boltzmann weight)
- `geometry_file`: Optional[str] = None (relative path)
- `optimized_geometry_file`: Optional[str] = None
- `rmsd_from_ref`: Optional[float] = None (from lowest-energy conformer)
- `status`: Literal[pending/optimizing/optimized/nmr_running/nmr_complete/failed]
- `error_message`: Optional[str] = None
- Tracks lifecycle from generation through NMR calculation

**ConformerEnsemble model:**
- `method`: Literal["rdkit_kdg", "crest"]
- `conformers`: list[ConformerData]
- `temperature_k`: float = 298.15 (for Boltzmann weighting)
- `pre_dft_energy_window_kcal`: float = 6.0 (RDKit MMFF filter)
- `post_dft_energy_window_kcal`: float = 3.0 (DFT energy filter)
- `total_generated`: int = 0 (before any filtering)
- `total_after_pre_filter`: int = 0 (after MMFF filter)
- `total_after_post_filter`: int = 0 (after DFT filter)
- Captures full conformer generation and filtering pipeline

**JobStatus extensions:**
- `conformer_mode`: Literal["single", "ensemble"] = "single"
- `conformer_ensemble`: Optional[ConformerEnsemble] = None
- Both fields optional with defaults → v1.x status.json loads unchanged

**JobInput extensions:**
- `conformer_mode`: Literal["single", "ensemble"] = "single"
- `conformer_method`: Optional[Literal["rdkit_kdg", "crest"]] = None
- `max_conformers`: Optional[int] = None (None = adaptive default)
- Enables API to specify conformer generation parameters

### Storage Helpers (src/qm_nmr_calc/storage.py)

**create_conformer_directories(job_id, conformer_ids):**
- Creates isolated scratch directories: `scratch/conformers/{conf_id}/`
- Creates per-conformer output: `output/conformers/{conf_id}/`
- Creates optimized geometry directory: `output/optimized/`
- Returns dict mapping conformer_id → scratch_dir Path
- **Critical for NWChem isolation:** Prevents database file conflicts

**get_conformer_scratch_dir(job_id, conformer_id):**
- Returns `scratch/conformers/{conformer_id}/` path
- Does NOT create (use create_conformer_directories first)

**get_conformer_output_dir(job_id, conformer_id):**
- Returns `output/conformers/{conformer_id}/` path

**get_optimized_conformers_dir(job_id):**
- Returns `output/optimized/` path for final geometries

### Test Coverage (tests/test_conformer_models.py)

**20 comprehensive tests covering:**
- ConformerData creation, defaults, status transitions, serialization
- ConformerEnsemble creation, defaults, custom parameters, serialization
- JobStatus backward compatibility (v1.x JSON loads with defaults)
- JobInput conformer field defaults
- Per-conformer directory creation and isolation
- Storage path helpers

**All tests pass.** No regressions in existing 104 tests.

## Technical Details

### Directory Structure (Multi-Conformer Job)

```
data/jobs/{job_id}/
├── status.json                   # JobStatus with conformer_ensemble
├── output/
│   ├── conformers/
│   │   ├── conf_001/            # Per-conformer outputs
│   │   ├── conf_002/
│   │   └── conf_003/
│   └── optimized/               # Final optimized geometries
│       ├── conf_001.xyz
│       ├── conf_002.xyz
│       └── conf_003.xyz
├── scratch/
│   └── conformers/
│       ├── conf_001/            # Isolated NWChem scratch
│       ├── conf_002/
│       └── conf_003/
└── logs/
```

**Why per-conformer scratch?** NWChem creates database files (.db, .lock) in scratch. Multiple NWChem processes sharing scratch causes corruption. Each conformer gets isolated scratch to prevent conflicts (Pitfall 5 from research).

### Backward Compatibility Strategy

**v1.x status.json (single conformer):**
```json
{
  "job_id": "abc123",
  "status": "complete",
  "input": {"smiles": "CCO", "solvent": "chcl3"},
  "nwchem_version": "7.0.2",
  "nmr_results": {...}
}
```

**Loads as:**
- `conformer_mode`: "single" (default)
- `conformer_ensemble`: None (default)
- All other fields unchanged

**v2.0 status.json (conformer ensemble):**
```json
{
  "job_id": "xyz789",
  "status": "running",
  "conformer_mode": "ensemble",
  "conformer_ensemble": {
    "method": "rdkit_kdg",
    "conformers": [
      {"conformer_id": "conf_001", "status": "optimized", "energy": -154.32},
      {"conformer_id": "conf_002", "status": "optimizing", "energy": null}
    ],
    "temperature_k": 298.15,
    "total_generated": 20,
    "total_after_pre_filter": 8
  }
}
```

**No migration needed.** Pydantic Optional defaults handle missing fields automatically.

## Decisions Made

### Per-Conformer Scratch Directory Isolation
**Decision:** Each conformer gets unique `scratch/conformers/{conf_id}/` directory.

**Rationale:** Research revealed NWChem creates database files (.db, .lock) in scratch. Multiple NWChem processes sharing scratch causes corruption. Isolated directories prevent conflicts.

**Impact:** Eliminates database corruption risk (Pitfall 5). Enables safe parallel conformer processing in future phases.

**Alternatives considered:**
1. Shared scratch with unique filenames → Complex naming, still risk of collision
2. Sequential processing only → Slow, wastes CPU cores

### Backward Compatibility via Optional Defaults
**Decision:** All new JobStatus/JobInput fields use `Optional` with defaults.

**Rationale:** Existing v1.x jobs have status.json files without conformer fields. Must load without error or migration script.

**Impact:**
- v1.x jobs continue working unchanged
- v2.0 can coexist with v1.x
- Smooth transition without data migration

**Alternatives considered:**
1. Version migrations → Complex, error-prone, requires downtime
2. Separate v2 schema → Breaks existing API clients

### Energy Unit Explicit Tracking
**Decision:** ConformerData stores `energy_unit` alongside `energy` value.

**Rationale:** Different stages use different units:
- RDKit MMFF: kcal/mol
- NWChem DFT: Hartree
- Boltzmann calculation: needs consistent units

**Impact:** Prevents unit conversion bugs. Explicit is better than implicit.

**Alternatives considered:**
1. Always convert on read → Requires knowing source, error-prone
2. Implicit unit assumption → Fragile, breaks when adding CREST (uses kJ/mol)

## Next Phase Readiness

### What This Enables

**Phase 13 (RDKit Conformer Generation):**
- Can populate ConformerEnsemble with generated conformers
- Can use create_conformer_directories() to set up storage
- Can track generation method and filter counts

**Phase 14 (Boltzmann Averaging):**
- Can read conformer energies and energy_unit
- Can calculate and store Boltzmann weights
- Can use temperature_k from ensemble

**Phase 15 (NWChem Conformer Loop):**
- Can iterate over conformers list
- Can use get_conformer_scratch_dir() for isolated NWChem runs
- Can update per-conformer status (pending → optimizing → optimized)

**Phase 16 (CREST Integration):**
- ConformerEnsemble supports method="crest"
- Can compare RDKit vs CREST conformer counts
- Can track different generation methods

**Phase 17 (API Conformer Support):**
- JobInput has conformer_mode and conformer_method fields
- API can accept conformer generation parameters
- JobStatus can return ensemble data in responses

### Blockers/Dependencies

**None.** This is the foundation phase. All downstream phases depend on these models.

### Known Risks

**Atom ordering consistency:** ConformerData stores geometry files, but doesn't yet track canonical atom order. Phase 13 must establish canonical ordering when generating conformers. If atom order changes between geometry optimization and NMR calculation, chemical shifts will map to wrong atoms.

**Mitigation:** Phase 13 must use canonical SMILES and save atom mapping alongside geometries.

## Deviations from Plan

None - plan executed exactly as written.

## Performance

**Duration:** 4.3 minutes (258 seconds)

**Breakdown:**
- Task 1 (models.py): ~1 min
- Task 2 (storage.py): ~1 min
- Task 3 (tests): ~2 min
- Verification: ~0.3 min

**Tests added:** 20 (all passing)
**Total tests:** 124 passing, 4 NWChem integration failures (pre-existing)

## Quality Checks

- All new models use `strict=True` (strict type checking)
- All new fields have explicit types (no Any)
- All Optional fields have explicit defaults
- All functions have docstrings with Args/Returns
- All directory operations use pathlib.Path (not strings)
- All tests follow existing patterns (pytest classes, descriptive names)
- No existing tests broken
- No existing functionality changed (additive only)

## Commits

| Task | Commit | Message |
|------|--------|---------|
| 1 | 6741e46 | feat(12-01): add conformer data models to models.py |
| 2 | e5d6bba | feat(12-01): add per-conformer directory helpers to storage.py |
| 3 | aaa4958 | test(12-01): add unit tests for conformer models and storage |

**Total changes:**
- 2 files modified (models.py, storage.py)
- 1 file created (test_conformer_models.py)
- 125 lines added (models.py)
- 83 lines added (storage.py)
- 342 lines added (tests)
- 550 total lines added

---

**Status:** Complete. Foundation for v2.0 conformational sampling established. All downstream phases can now build on ConformerData/ConformerEnsemble models and per-conformer storage structure.
