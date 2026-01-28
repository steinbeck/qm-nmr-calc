---
phase: 17-api-integration
plan: 01
subsystem: tasks
tags: [huey, ensemble, api, dispatch, progress]
depends_on:
  requires:
    - 15-03: run_ensemble_dft_and_nmr function
    - 16-03: generate_conformer_ensemble with CREST dispatch
    - 14-02: average_ensemble_nmr function
  provides:
    - run_ensemble_nmr_task Huey task
    - Mode-aware task dispatch in API endpoints
    - Progress callback support in runner functions
  affects:
    - 17-02: Schema extensions already done (out of order execution)
    - 17-03: Web UI will use new endpoints

tech-stack:
  added: []
  patterns:
    - Mode-aware task dispatch
    - Progress callback pattern for long-running operations
    - CREST fallback with warning mechanism

key-files:
  created:
    - tests/test_ensemble_task.py
  modified:
    - src/qm_nmr_calc/tasks.py
    - src/qm_nmr_calc/storage.py
    - src/qm_nmr_calc/nwchem/runner.py
    - src/qm_nmr_calc/models.py
    - src/qm_nmr_calc/api/routers/jobs.py
    - src/qm_nmr_calc/api/routers/web.py

decisions:
  - key: progress-callback-pattern
    choice: "Callback with (step, current, total) tuple"
    reason: "Flexible for both DFT optimization and NMR calculation steps"
  - key: crest-fallback-warning
    choice: "Store warning in job status, continue with RDKit"
    reason: "User sees CREST unavailable but job still completes"
  - key: web-default-ensemble
    choice: "Web UI defaults to conformer_mode=ensemble"
    reason: "CONTEXT.md decision - users opt OUT to single-conformer"

metrics:
  duration: "17 min"
  completed: "2026-01-28"
---

# Phase 17 Plan 01: Ensemble Task and API Dispatch Summary

Ensemble Huey task with mode-aware dispatch, progress callbacks, and CREST fallback handling.

## What Was Built

### 1. run_ensemble_nmr_task (tasks.py)
New Huey task that orchestrates ensemble NMR calculations:
- Generates conformer ensemble via pipeline (RDKit KDG or CREST)
- Runs DFT optimization + NMR on all conformers via run_ensemble_dft_and_nmr
- Filters to nmr_complete conformers before Boltzmann averaging (critical: prevents length mismatch)
- Generates spectrum plots and annotated structure from averaged shifts
- Handles CREST fallback with warning when CREST unavailable

### 2. Progress Callback Support (runner.py)
Added `progress_callback` parameter to:
- `run_ensemble_dft_and_nmr()`
- `run_conformer_dft_optimization()`
- `run_conformer_nmr_calculations()`

Callback signature: `(step: str, current: int, total: int)`
Called after each conformer completes (success or failure) for real-time progress updates.

### 3. Storage Conformer Parameters (storage.py)
Extended `create_job_directory()` to accept:
- `conformer_mode: str = "single"` (v1.x default, backward compatible)
- `conformer_method: Optional[str] = None` (rdkit_kdg or crest)
- `max_conformers: Optional[int] = None` (None = adaptive)

### 4. API Endpoint Mode-Aware Dispatch (jobs.py, web.py)
Both API and web endpoints now:
- Pass conformer params to `create_job_directory()`
- Dispatch to `run_ensemble_nmr_task` when `conformer_mode="ensemble"`
- Dispatch to `run_nmr_task` when `conformer_mode="single"`

Web UI defaults to ensemble mode per CONTEXT.md decision.

### 5. Model Updates (models.py)
Added `conformer_method_warning: Optional[str]` to JobStatus for CREST fallback tracking.

## Key Implementation Details

### Conformer Filtering Before Averaging
```python
# CRITICAL: Filter to only nmr_complete conformers before averaging
nmr_complete_conformers = [c for c in ensemble.conformers if c.status == "nmr_complete"]
filtered_ensemble = ensemble.model_copy(update={"conformers": nmr_complete_conformers})
avg_nmr = average_ensemble_nmr(filtered_ensemble, nmr_results)
```

This prevents the `len(conformers) != len(nmr_results)` error documented in STATE.md blockers.

### CREST Fallback Pattern
```python
try:
    ensemble = generate_conformer_ensemble(conformer_method="crest", ...)
except ValueError as e:
    if "CREST" in str(e) and conformer_method == "crest":
        warning_msg = str(e) + " Falling back to RDKit KDG method."
        update_job_status(job_id, conformer_method_warning=warning_msg)
        ensemble = generate_conformer_ensemble(conformer_method="rdkit_kdg", ...)
```

## Commits

| Commit | Description |
|--------|-------------|
| 12df9ba | feat(17-01): add run_ensemble_nmr_task and progress tracking |
| 00d1fcd | feat(17-01): wire API endpoints for mode-aware task dispatch |
| 29866fc | test(17-01): add ensemble task and mode-aware dispatch tests |

## Test Coverage

13 new tests in `tests/test_ensemble_task.py`:
- Import and signature tests for task and runner functions
- Storage parameter acceptance tests
- Ensemble task filters nmr_complete conformers before averaging
- API dispatch tests (ensemble and single modes)
- CREST fallback warning test

All 13 tests passing.

## Deviations from Plan

None - plan executed exactly as written.

## Verification Results

All verification criteria from plan met:
- [x] `run_ensemble_nmr_task` exists and can be imported
- [x] `run_ensemble_nmr_task` generates conformers using pipeline
- [x] `run_ensemble_nmr_task` runs DFT+NMR using `run_ensemble_dft_and_nmr`
- [x] `run_ensemble_nmr_task` filters to nmr_complete before Boltzmann averaging
- [x] `run_ensemble_nmr_task` handles CREST fallback with warning
- [x] `run_ensemble_nmr_task` reports progress via callback (X/N counts)
- [x] `run_ensemble_dft_and_nmr` accepts `progress_callback` parameter
- [x] API `submit_smiles` dispatches to correct task based on `conformer_mode`
- [x] Web `submit_job` dispatches to correct task based on `conformer_mode`
- [x] `create_job_directory` accepts and stores conformer parameters
- [x] Visualizations generated using `generate_spectrum_plot` and `generate_annotated_structure`
- [x] New tests pass (13/13)
- [x] Existing tests still pass

## Next Phase Readiness

Plan 17-01 complete. Foundation for v2.0 API integration established:
- Ensemble task ready for end-to-end testing
- Progress callback infrastructure in place for UI polling
- CREST fallback mechanism operational

Note: Plan 17-02 (schema extensions) was already completed in a previous session.
Plan 17-03 (web UI) can now use the mode-aware dispatch.
