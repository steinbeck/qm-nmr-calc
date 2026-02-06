---
phase: 49-config-foundation-and-pricing-query
plan: 01
subsystem: gcp-deployment
tags: [toml, pydantic, config-validation, tdd, bash]
dependency-graph:
  requires: []
  provides: [gcp-config-model, config-toml-format, bash-config-loader, config-validation-tests]
  affects: [49-02, 50, 51, 52]
tech-stack:
  added: []
  patterns: [pydantic-v2-validation, toml-config, bash-python-bridge, tdd-red-green]
key-files:
  created:
    - gcp/__init__.py
    - gcp/validate_config.py
    - gcp/config.toml.example
    - gcp/lib/config.sh
    - tests/test_gcp_config.py
  modified:
    - .gitignore
decisions:
  - id: D-49-01-01
    description: "Used ram_gb ge=1 (not ge=8) in Pydantic model because RAM/CPU ratio validator already enforces minimum effective RAM"
  - id: D-49-01-02
    description: "Added !gcp/lib/ negation to .gitignore since top-level lib/ pattern blocks gcp/lib/ directory"
metrics:
  duration: ~8 min (actual work time)
  completed: 2026-02-06
---

# Phase 49 Plan 01: TOML Config Validation Summary

**One-liner:** Pydantic v2 config model with TOML parsing, bash export bridge, and 19 TDD tests for GCP deployment parameters.

## What Was Done

### Task 1: Write failing tests (RED phase)
Created `tests/test_gcp_config.py` with 19 test cases covering:
- Valid config parsing (full and minimal with defaults)
- Project ID validation: missing, too short, too long, uppercase, leading/trailing hyphens
- RAM/CPU ratio validation: too low (0.25), too high (16), boundary values (0.5, 8.0)
- Field bounds: cpu_cores minimum (4), disk_size_gb minimum (10)
- TOML loading: valid file, missing [gcp] section, invalid syntax, file not found
- Export formatting: all 5 bash export statements present

All 19 tests failed with `ModuleNotFoundError` (RED confirmed).

### Task 2: Implement and make tests pass (GREEN phase)

**gcp/validate_config.py** - Python module (importable + CLI script):
- `GCPConfig(BaseModel)` with field validators for project_id format and RAM/CPU ratio
- `load_config(path)` - TOML parsing with section and syntax error handling
- `format_exports(config)` - bash-compatible export statements with shell injection prevention
- CLI mode with `--config` argument, user-friendly error messages

**gcp/config.toml.example** - Commented template with all fields documented

**gcp/lib/config.sh** - Bash library bridging Python validation to shell scripts via `eval $(python3 validate_config.py)`

**gcp/__init__.py** - Empty file enabling package imports for testing

**.gitignore** - Added `gcp/config.toml`, `gcp/.cache/`, and `!gcp/lib/` negation

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] .gitignore lib/ pattern blocks gcp/lib/**
- **Found during:** Task 2, Step 3
- **Issue:** Top-level `lib/` in `.gitignore` (Python packaging pattern) matches `gcp/lib/` recursively
- **Fix:** Added `!gcp/lib/` negation at end of `.gitignore`
- **Files modified:** `.gitignore`
- **Commit:** `6fad586`

## Verification Results

| Check | Result |
|-------|--------|
| 19 new tests pass | PASS |
| Existing tests pass (249+ checked) | PASS (1 pre-existing failure in test_conformer_nwchem unrelated) |
| CLI rejects invalid config | PASS |
| CLI produces bash exports | PASS |
| config.toml in .gitignore | PASS |
| gcp/lib/config.sh trackable | PASS |

## Commits

| Commit | Type | Description |
|--------|------|-------------|
| `2110d82` | test | Failing tests for GCP config validation (RED) |
| `6fad586` | feat | Implement TOML config validation with Pydantic v2 (GREEN) |

## Next Phase Readiness

Plan 49-02 (spot pricing query) can proceed. It will import from `gcp.validate_config` and build the pricing module alongside it. The `gcp/lib/` directory and `gcp/__init__.py` are ready for additional modules.
