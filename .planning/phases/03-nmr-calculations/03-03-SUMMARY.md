---
phase: 03-nmr-calculations
plan: 03
subsystem: calculation-pipeline
tags: [nmr, dft, isicle, nwchem, huey, api]

dependency-graph:
  requires: ["03-01", "03-02"]
  provides: ["run_nmr_calculation", "run_nmr_task", "solvent-validation", "/solvents-endpoint"]
  affects: ["04-output", "05-testing"]

tech-stack:
  added: []
  patterns: ["two-step-dft-workflow", "shielding-to-shift-conversion", "solvent-validation"]

key-files:
  created: []
  modified:
    - src/qm_nmr_calc/isicle_wrapper.py
    - src/qm_nmr_calc/tasks.py
    - src/qm_nmr_calc/api/routers/jobs.py

decisions: []

metrics:
  duration: 5 min
  completed: 2026-01-19
---

# Phase 3 Plan 3: NMR Calculation Pipeline Summary

Two-step DFT workflow (geometry optimization + NMR shielding), Huey task orchestration, and API solvent validation with /solvents endpoint.

## What Was Built

### 1. run_nmr_calculation Function (isicle_wrapper.py)

Added `run_nmr_calculation()` function implementing the two-step DFT workflow:

```python
def run_nmr_calculation(
    smiles: str,
    job_dir: Path,
    preset: CalculationPreset,
    solvent: str,
    processes: int = 4,
) -> dict:
    """
    Run geometry optimization + NMR shielding calculation via ISiCLE/NWChem.
    Returns: dict with 'geometry_file', 'shielding_data', 'energy'
    """
```

Workflow steps:
1. Load molecule from SMILES via `isicle.load()`
2. Initial 3D embedding with UFF force field
3. DFT geometry optimization with preset's `basis_set`
4. Save optimized geometry to output/optimized.xyz
5. DFT NMR shielding calculation with preset's `nmr_basis_set`
6. Validate shielding data present (raises RuntimeError if not)

### 2. run_nmr_task Huey Task (tasks.py)

Added `run_nmr_task()` Huey task that orchestrates the complete NMR calculation:

```python
@huey.task()
def run_nmr_task(job_id: str) -> dict:
    """Execute NMR calculation (geometry optimization + shielding) for a queued job."""
```

Task workflow:
1. Load job status and validate it's queued
2. Get preset configuration from job input
3. Run `run_nmr_calculation()` for two-step DFT
4. Convert shielding to shifts via `shielding_to_shift()`
5. Build `NMRResults` with `AtomShift` objects
6. Save `nmr_results.json` to output directory
7. Update job status with NMR results

### 3. API Solvent Validation (api/routers/jobs.py)

Updated both submission endpoints to validate solvents:

- `POST /api/v1/jobs` - SMILES submission with solvent validation
- `POST /api/v1/jobs/upload` - File upload with solvent validation

Invalid solvent returns 422 with RFC 7807 ProblemDetail:
```json
{
  "type": "https://qm-nmr-calc.example/problems/invalid-solvent",
  "title": "Invalid Solvent",
  "status": 422,
  "detail": "Unknown solvent 'xyz'. Supported solvents: acetntrl, acetone, ..."
}
```

### 4. Solvents List Endpoint

Added `GET /api/v1/jobs/solvents` endpoint:

```python
@router.get("/solvents", response_model=list[str])
async def list_solvents():
    """List supported NMR solvents for COSMO solvation model."""
    return get_supported_solvents()  # Returns 11 common NMR solvents
```

## Technical Decisions

| Decision | Rationale |
|----------|-----------|
| Keep run_geometry_optimization | Existing function retained for potential standalone use |
| Validate shielding data presence | RuntimeError if NMR calculation doesn't produce data |
| Normalize solvent before storing | Ensures consistent lowercase storage |
| run_nmr_task updates job status | Signal handlers handle running/complete/failed transitions |

## Verification Results

All verification tests passed:

```
1. Invalid solvent returns 422: PASS
2. Valid submission returns 202 with preset=production, solvent=chcl3: PASS
3. /solvents returns list with 11 solvents including chcl3: PASS

API tests passed
```

## Deviations from Plan

None - plan executed exactly as written.

## Commits

| Commit | Description |
|--------|-------------|
| a91cf58 | feat(03-03): add run_nmr_calculation for two-step DFT workflow |
| b81b05e | feat(03-03): add run_nmr_task for complete NMR calculation |
| 4a976da | feat(03-03): add solvent validation and /solvents endpoint |

## Next Phase Readiness

### Blockers
None identified.

### Prerequisites Met
- [x] run_nmr_calculation performs two-step DFT
- [x] run_nmr_task orchestrates calculation and stores results
- [x] API accepts preset and solvent parameters
- [x] Invalid solvent returns 422
- [x] /solvents endpoint returns supported list

### Ready For
- Phase 04: Output formatting and download endpoints
- Phase 05: Testing (integration tests for NMR workflow)

---

*Phase: 03-nmr-calculations*
*Plan: 03*
*Completed: 2026-01-19*
