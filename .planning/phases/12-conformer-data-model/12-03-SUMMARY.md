---
phase: 12-conformer-data-model
plan: "03"
subsystem: api-schemas
status: complete
tags: [backward-compat, api-schemas, pydantic, conformer-mode, testing]
dependencies:
  requires:
    - 12-01-conformer-data-models
    - 12-02-canonical-atom-ordering
  provides:
    - backward-compatibility-verified
    - api-conformer-mode-schema
  affects:
    - 17-api-conformer-support
tech-stack:
  added: []
  patterns:
    - fixture-based-backward-compat-testing
    - optional-api-fields-with-defaults
key-files:
  created:
    - tests/test_backward_compat.py
  modified:
    - src/qm_nmr_calc/api/schemas.py
    - src/qm_nmr_calc/api/routers/jobs.py
decisions:
  - id: API-COMPAT-001
    title: Use str type for conformer_mode in API schemas (not Literal)
    choice: conformer_mode is str with default "single" in API schemas
    rationale: API flexibility -- internal models use Literal for validation but API schemas use str to avoid coupling API contract to internal enum
    alternatives: ["Literal type in API schema", "Enum type"]
    impact: Future conformer modes can be added without API schema changes
metrics:
  tests-added: 8
  tests-passing: 132
  duration: "5.6 min"
  completed: "2026-01-26"
---

# Phase 12 Plan 03: Backward Compatibility and API Schema Update Summary

**One-liner:** 8-test backward compatibility suite proving v1.x JSON loads correctly, plus conformer_mode added to API request/response schemas for v2.0 readiness.

## What Was Built

### Backward Compatibility Test Suite (tests/test_backward_compat.py)

**8 tests with realistic v1.x JSON fixtures:**

1. `test_v1_status_json_loads` -- Complete v1.x status.json loads via model_validate
2. `test_v1_status_defaults` -- Loaded v1 job has conformer_mode="single", conformer_ensemble=None
3. `test_v1_input_defaults` -- Loaded v1 input has correct conformer defaults
4. `test_v1_status_preserves_all_fields` -- All existing fields (job_id, NMR results, steps, resources) preserved
5. `test_v1_status_roundtrip` -- Load, dump, reload cycle preserves all values
6. `test_v1_queued_job_loads` -- Minimal queued job (no results) loads
7. `test_v1_failed_job_loads` -- Failed job with error message loads
8. `test_v1_job_input_loads` -- Minimal JobInput loads with conformer defaults

**Fixtures cover three v1.x job states:**
- Complete job with NMR results, steps, and geometry
- Queued job (minimal fields)
- Failed job with error details

### API Schema Updates (src/qm_nmr_calc/api/schemas.py)

**JobSubmitRequest -- new fields:**
- `conformer_mode: str = "single"` -- Sampling mode (single or ensemble)
- `conformer_method: Optional[str] = None` -- Generation method (rdkit_kdg or crest)
- `max_conformers: Optional[int] = None` -- Max conformers (None = adaptive default)
- All fields have Field descriptions for OpenAPI documentation

**JobStatusResponse -- new field:**
- `conformer_mode: str = "single"` -- Reflects job's conformer mode in responses

### Response Conversion Update (src/qm_nmr_calc/api/routers/jobs.py)

**job_status_to_response():**
- Added `conformer_mode: job_status.conformer_mode` to response dict
- Single-line addition, minimal change surface

## Technical Details

### Backward Compatibility Strategy

The test suite uses three inline dict fixtures that exactly match real v1.x status.json structure (no conformer fields present). The tests verify that Pydantic's Optional defaults handle the missing fields automatically:

- Missing `conformer_mode` defaults to `"single"`
- Missing `conformer_ensemble` defaults to `None`
- Missing `conformer_method` defaults to `None`
- Missing `max_conformers` defaults to `None`

### API Schema Design Choice

API schemas use `str` type for `conformer_mode` rather than `Literal["single", "ensemble"]`:
- Internal models (models.py) use Literal for strict validation
- API schemas use str for flexibility -- avoids coupling API contract to internal type
- Validation of allowed values happens at the service layer, not the API schema layer
- This means future conformer modes can be added without API schema version bump

### Test Suite Status

| Category | Count | Status |
|----------|-------|--------|
| Backward compat (new) | 8 | All passing |
| Conformer models (12-01) | 20 | All passing |
| Atom ordering (12-02) | 13 | All passing |
| API tests | 11 | All passing |
| Other existing tests | 80 | All passing |
| NWChem integration | 4 | Pre-existing failures (env) |
| **Total** | **136** | **132 passing** |

## Decisions Made

### str Type for API conformer_mode

**Decision:** Use `str` with default `"single"` in API schemas rather than `Literal`.

**Rationale:** API schemas define the external contract. Using `str` allows adding new conformer modes (e.g., "targeted") without changing the API schema. Internal model validation (via `Literal` in models.py) catches invalid values before they reach storage.

**Impact:** Looser API coupling. Validation still enforced at service layer.

## Deviations from Plan

None -- plan executed exactly as written.

## Next Phase Readiness

### What This Completes

**Phase 12 is now fully complete:**
- Plan 01: ConformerData, ConformerEnsemble models, storage helpers
- Plan 02: Canonical atom ordering (atom_ordering.py)
- Plan 03: Backward compatibility verified, API schemas updated

**Requirements verified:**
- DFT-03 (partial): Conformer data structures ready for per-conformer DFT
- API-04 (partial): conformer_mode in API schemas, v1.x backward compat proven

### What Phase 13 Needs

Phase 13 (RDKit Conformer Generation) can now:
- Create ConformerEnsemble and ConformerData instances
- Use create_conformer_directories() for storage setup
- Accept conformer_mode from API via JobSubmitRequest
- Return conformer_mode in API responses via JobStatusResponse
- Use canonical_atom_order() from atom_ordering.py for consistent indexing

## Commits

| Task | Commit | Message |
|------|--------|---------|
| 1 | 1dba9a3 | test(12-03): add backward compatibility tests for v1.x JSON fixtures |
| 2 | 53d4049 | feat(12-03): add conformer_mode to API schemas and response |

**Total changes:**
- 1 file created (test_backward_compat.py, 204 lines)
- 2 files modified (schemas.py +17 lines, jobs.py +1 line)
- 222 total lines added

---

**Status:** Complete. Phase 12 (Conformer Data Model and Storage) is fully finished. All data models, atom ordering, backward compatibility, and API schemas are ready for Phase 13 (RDKit Conformer Generation).
