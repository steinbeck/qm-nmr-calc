---
phase: 08-delta50-setup
plan: 02
subsystem: benchmark
tags: [benchmark, cli, runner, tdd, nwchem]

dependency-graph:
  requires: ["08-01"]
  provides: ["benchmark-runner", "benchmark-cli"]
  affects: ["08-03", "09-*", "10-*"]

tech-stack:
  added: []
  patterns: ["CLI subcommands", "resume-on-marker-file", "JSONL progress tracking"]

key-files:
  created:
    - src/qm_nmr_calc/benchmark/runner.py
    - src/qm_nmr_calc/benchmark/__main__.py
    - tests/test_benchmark.py
  modified:
    - src/qm_nmr_calc/benchmark/__init__.py

decisions:
  - key: WP04 functional preset
    choice: "Separate BENCHMARK_PRESETS dict with WP04 config"
    reason: "WP04 optimized for 1H shifts, not in existing presets.py"
  - key: Resume detection
    choice: "Check for shifts.json marker file"
    reason: "Simple existence check, atomic write ensures completeness"
  - key: compound_XX naming
    choice: "Use compound_XX not molecule_XX"
    reason: "Match actual file naming from DELTA50 dataset"

metrics:
  duration: "3 min"
  completed: "2026-01-21"
---

# Phase 8 Plan 2: Benchmark Runner and CLI Summary

**One-liner:** CLI benchmark runner with resume support, WP04 functional, and 200-task matrix (50 molecules x 2 functionals x 2 solvents)

## What Was Built

### 1. Benchmark Runner (`runner.py`)
- `run_benchmark()`: Execute DELTA50 calculation matrix with tqdm progress
- `build_task_matrix()`: Generate 200 tasks (50 x B3LYP/WP04 x CHCl3/DMSO)
- `is_task_complete()`: Resume support via shifts.json marker file
- `run_single_calculation()`: Single molecule/method/solvent calculation
- `aggregate_results()`: Collect results into pandas DataFrame
- `show_status()`: Display progress (0/200 initially)
- `BENCHMARK_PRESETS`: B3LYP and WP04 configurations

### 2. CLI Entry Point (`__main__.py`)
- `python -m qm_nmr_calc.benchmark run`: Execute with filters and resume
- `python -m qm_nmr_calc.benchmark status`: Show progress
- `python -m qm_nmr_calc.benchmark summary`: Generate CSV summary

### 3. Unit Tests (`test_benchmark.py`)
- 16 tests covering data loader, models, and runner
- Tests for task matrix building, completion detection, preset validation

## Key Technical Details

**WP04 Functional Preset:**
```python
"WP04": {
    "functional": "wp04",  # Optimized for 1H shifts
    "basis_set": "6-31G*",
    "nmr_basis_set": "6-311++G(2d,p)",
    "max_iter": 150,
}
```

**Task Matrix Structure:**
- 50 molecules: compound_01 through compound_50
- 2 functionals: B3LYP, WP04
- 2 solvents: CHCl3, DMSO
- Total: 200 calculations

**Resume Logic:**
```python
def is_task_complete(results_dir: Path, task: dict) -> bool:
    output_dir = results_dir / task["molecule_id"] / f"{task['functional']}_{task['solvent']}"
    return (output_dir / "shifts.json").exists()
```

## Files Changed

| File | Change | Purpose |
|------|--------|---------|
| `src/qm_nmr_calc/benchmark/runner.py` | Created | Core runner with resume support |
| `src/qm_nmr_calc/benchmark/__main__.py` | Created | CLI entry point with subcommands |
| `src/qm_nmr_calc/benchmark/__init__.py` | Modified | Export runner functions |
| `tests/test_benchmark.py` | Created | 16 unit tests |

## Decisions Made

1. **Separate BENCHMARK_PRESETS**: WP04 functional not in main presets.py, kept separate for benchmark use
2. **shifts.json marker**: Resume checks for this specific file, atomic write ensures completeness
3. **compound_XX naming**: Matches actual DELTA50 XYZ file naming convention
4. **Production scaling**: Use existing scaling factors (benchmark-specific derived in Phase 10)

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed molecule naming convention**
- **Found during:** Task 1 implementation
- **Issue:** Plan used `molecule_XX` but actual files are `compound_XX`
- **Fix:** Updated all references to use `compound_XX` pattern
- **Files modified:** runner.py, test_benchmark.py

## Verification Results

```
$ python -m qm_nmr_calc.benchmark status
DELTA50 Benchmark Status
========================
Completed: 0/200 (0.0%)
Remaining: 200

$ pytest tests/test_benchmark.py -v
16 passed in 0.87s
```

## Next Phase Readiness

**Ready for 08-03:** Error metrics calculation module
- Runner produces shifts.json with calculated shifts
- Experimental shifts available for comparison
- MAE/RMSE calculation can compare calculated vs experimental
