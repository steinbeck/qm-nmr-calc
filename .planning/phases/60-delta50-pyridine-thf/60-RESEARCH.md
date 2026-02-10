# Phase 60: DELTA50 Pyridine + THF - Research

**Researched:** 2026-02-10
**Domain:** Quantum chemistry benchmark calculations with NWChem
**Confidence:** HIGH

## Summary

Phase 60 runs the established DELTA50 benchmark pipeline for 2 new solvents: pyridine and THF. This phase leverages the exact infrastructure built in v2.8 (Phases 54-58) and extended in Phase 59.

The benchmark pipeline has been proven with 200 calculations across 4 solvents (methanol, water, acetone, benzene) in v2.8. Phase 59 just added pyridine and THF to the supported solvent list with proper COSMO name mappings. All infrastructure is ready - this is purely an execution phase.

**Primary recommendation:** Run the existing benchmark CLI with `--solvents Pyridine THF --functionals B3LYP --headless` flag. No code changes needed. The system will execute 100 calculations (50 molecules × 2 solvents) using B3LYP/6-311+G(2d,p) with COSMO solvation, storing results in the standard structure at `data/benchmark/results/{compound_name}/{functional}_{solvent}/`.

## Standard Stack

This phase uses the existing benchmark infrastructure - no new libraries required.

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| NWChem | 7.2.0+ | Quantum chemistry calculations | Only DFT engine with COSMO solvation support in this project |
| rdkit | 2025.9.3+ | Molecule handling | Standard cheminformatics library |
| tqdm | 4.67.1+ | Progress tracking | Already used in runner.py for benchmark progress |
| orjson | 3.11.5+ | Fast JSON I/O | Already used for shifts.json output |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| pandas | 2.3.3+ | Result aggregation | Post-execution summary generation |
| matplotlib | 3.10.0+ | Visualization | Phase 61 (analysis) needs this |
| statsmodels | 0.14.0+ | OLS regression | Phase 61 (scaling factor derivation) |

### Command Line Interface
The CLI is already complete from Phase 59:

```bash
python -m qm_nmr_calc.benchmark run \
  --solvents Pyridine THF \
  --functionals B3LYP \
  --headless \
  --processes 4
```

**Installation:** No new dependencies. Project already has all required packages in pyproject.toml.

## Architecture Patterns

### Benchmark Execution Pattern (Already Implemented)

The benchmark runner follows a proven architecture from v2.8:

```
1. Task Matrix Build (runner.py:build_task_matrix)
   - 50 molecules × 2 solvents × 1 functional = 100 tasks

2. Resume Support (runner.py:is_task_complete)
   - Checks for shifts.json marker files
   - Skips completed calculations automatically

3. Sequential Execution with Monitoring (runner.py:run_benchmark)
   - Status updates to status.json (atomic writes)
   - Progress logging to benchmark.log
   - JSONL append-only progress tracking
   - Graceful stop via STOP file
   - Failure threshold (>5 unique molecules fail → pause)

4. Result Storage (runner.py:run_single_calculation)
   - Directory: data/benchmark/results/{molecule_id}/{functional}_{solvent}/
   - Files:
     - shifts.json (h1_shifts, c13_shifts, shielding_data)
     - nwchem.out (full NWChem output)
     - molecule_geom.xyz (optimized geometry)
     - other NWChem artifacts
```

**Example directory structure after execution:**
```
data/benchmark/results/
├── compound_01/
│   ├── B3LYP_Pyridine/
│   │   ├── shifts.json
│   │   ├── nwchem.out
│   │   └── molecule_geom.xyz
│   └── B3LYP_THF/
│       ├── shifts.json
│       ├── nwchem.out
│       └── molecule_geom.xyz
├── compound_02/
│   ├── B3LYP_Pyridine/
│   └── B3LYP_THF/
...
├── benchmark.log
├── status.json
└── progress.jsonl
```

### COSMO Solvation Integration (Already Implemented)

Phase 59 added pyridine and THF support:

**Solvent name flow:**
```
CLI: "Pyridine" (Title-case)
  ↓ runner.py:234 - solvent.lower()
"pyridine" (lowercase)
  ↓ input_gen.py:44 - _validate_solvent()
Validates against SUPPORTED_SOLVENTS ✓
  ↓ input_gen.py:88 - _get_cosmo_solvent_name()
"pyridine" (passthrough, no mapping needed)
  ↓ input_gen.py:92 - COSMO block generation
NWChem input: "solvent pyridine"
```

**NWChem COSMO parameters (verified from official docs):**
- Pyridine: dielectric constant = 12.978
- THF: dielectric constant = 7.4257
- Both are in NWChem's built-in COSMO parameter database
- No custom parameters needed

### Calculation Settings (Already Configured)

From `runner.py:BENCHMARK_PRESETS`:

**B3LYP preset:**
```python
{
    "functional": "b3lyp",
    "basis_set": "6-31G*",              # Geometry optimization
    "nmr_basis_set": "6-311+G(2d,p)",   # NMR shielding calculation
    "max_iter": 150
}
```

**Two-stage calculation workflow:**
1. Geometry optimization with 6-31G* basis (skip_optimization=True for DELTA50)
2. NMR shielding with 6-311+G(2d,p) basis

**Why skip_optimization=True:** DELTA50 molecules have pre-optimized geometries in `data/benchmark/delta50/molecules/*.xyz`. Phase uses these directly to match DELTA50 benchmark protocol.

### Anti-Patterns to Avoid

**DO NOT:**
- Re-run already-complete calculations (use `--resume` flag, which is default)
- Run without `--headless` for long calculations (tqdm overhead adds up)
- Modify `runner.py:SOLVENTS` list (new solvents are opt-in via CLI only)
- Skip monitoring status.json (contains ETA and failure tracking)

## Don't Hand-Roll

Phase 60 requires zero new code. All problems already solved:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Benchmark task queue | Custom job scheduler | runner.py:build_task_matrix + sequential loop | Already handles resume, monitoring, graceful stop |
| Progress tracking | Custom logger | status.json + tqdm + benchmark.log | Atomic status updates, ETA calculation, rolling average |
| NWChem input generation | String templates | input_gen.py functions | COSMO mapping, validated solvents, ISiCLE-compatible format |
| Result validation | Manual checks | is_task_complete() checking shifts.json | Idempotent, allows safe resume |
| Solvent name mapping | Hardcoded strings | SUPPORTED_SOLVENTS + _get_cosmo_solvent_name() | Centralized, tested, handles acetonitrile→acetntrl special case |

**Key insight:** v2.8 ran 200 calculations successfully. Infrastructure is production-ready. Phase 60 is pure execution.

## Common Pitfalls

### Pitfall 1: Running Without Headless Mode
**What goes wrong:** Long-running benchmarks (100 calculations × ~400s each = ~11 hours) time out or lose connection if run interactively.

**Why it happens:** tqdm creates terminal output that can cause buffer issues over long periods. SSH sessions can disconnect.

**How to avoid:**
```bash
# Correct: headless with nohup
nohup python -m qm_nmr_calc.benchmark run \
  --solvents Pyridine THF \
  --functionals B3LYP \
  --headless \
  > /dev/null 2>&1 &

# Monitor progress
python -m qm_nmr_calc.benchmark status
```

**Warning signs:** Terminal hangs, SSH disconnects, "Broken pipe" errors.

### Pitfall 2: Not Monitoring Failures
**What goes wrong:** If >5 unique molecules fail, the runner pauses automatically. Unmonitored runs can stall for hours before human notices.

**Why it happens:** NWChem can fail due to SCF convergence issues, memory limits, or COSMO parameter problems (though pyridine/THF are well-tested).

**How to avoid:**
- Check `status.json` regularly: `python -m qm_nmr_calc.benchmark status`
- Look for `"state": "paused"` or growing `"failures"` array
- Monitor `benchmark.log` for error patterns

**Warning signs:** Progress stalls, `estimated_remaining_hours` not decreasing, `"failed": N` where N > 0.

### Pitfall 3: Assuming All Molecules Will Succeed
**What goes wrong:** Some DELTA50 molecules are challenging (large, flexible, polar). Expecting 100/100 success rate leads to disappointment.

**Why it happens:** DFT convergence is not guaranteed. COSMO adds numerical complexity.

**How to avoid:**
- Accept 95%+ success rate as excellent (≤2-3 failures out of 100)
- Review failures in status.json to identify patterns
- Document any systematic failures (e.g., all failures are same compound)

**Warning signs:** Many failures for same molecule across solvents suggests geometry issue. Failures specific to one solvent suggests COSMO parameter issue.

### Pitfall 4: Not Using Resume Flag
**What goes wrong:** Re-running from scratch after a partial run wastes 4-8 hours of computation.

**Why it happens:** Forgetting that `--resume` is the default, or accidentally using `--no-resume`.

**How to avoid:**
- Default behavior is resume (checks for shifts.json existence)
- Only use `--no-resume` if explicitly wanting to overwrite results
- Verify resume with: `ls data/benchmark/results/compound_01/B3LYP_Pyridine/shifts.json`

**Warning signs:** Seeing "Running 100 benchmark calculations" when you expected fewer (should say "Running X calculations" where X < 100 if some already complete).

### Pitfall 5: Running Both WP04 and B3LYP Together
**What goes wrong:** Phase 60 requirements only specify B3LYP. Adding WP04 doubles execution time (100 → 200 tasks).

**Why it happens:** Seeing WP04 in `BENCHMARK_PRESETS` and assuming it should be run.

**How to avoid:**
- Phase 60 success criteria: "All 50 pyridine DELTA50 molecules calculate successfully with COSMO solvation" (singular solvent test, not cross-functional)
- Run B3LYP only: `--functionals B3LYP`
- WP04 would be a separate phase (not in v2.9 requirements)

**Warning signs:** Task count is 200 instead of 100, ETA is ~22 hours instead of ~11 hours.

## Code Examples

No code changes needed. All patterns already implemented.

### Running the Benchmark (Production Command)

```bash
# Activate virtual environment
source .venv/bin/activate

# Run benchmark in background with headless mode
nohup python -m qm_nmr_calc.benchmark run \
  --solvents Pyridine THF \
  --functionals B3LYP \
  --headless \
  --processes 4 \
  > /dev/null 2>&1 &

# Save process ID for later
echo $! > benchmark.pid
```

**Source:** `src/qm_nmr_calc/benchmark/__main__.py:240-289` (CLI implementation)

### Monitoring Progress

```bash
# Check overall status
python -m qm_nmr_calc.benchmark status

# Tail the log file
tail -f data/benchmark/results/benchmark.log

# Check status.json programmatically
cat data/benchmark/results/status.json | jq '.'
```

**Source:** `src/qm_nmr_calc.benchmark/__main__.py:147-165` (status command)

### Graceful Stop

```bash
# Request stop (completes current calculation, then stops)
python -m qm_nmr_calc.benchmark stop

# Check if acknowledged
python -m qm_nmr_calc.benchmark status
# Should show: "State: STOPPED"
```

**Source:** `src/qm_nmr_calc/benchmark/runner.py:48-62` (stop mechanism)

### Verifying Completion

```bash
# Count completed calculations
find data/benchmark/results -name "shifts.json" -path "*/B3LYP_Pyridine/*" | wc -l
# Should be 50 (or close to it)

find data/benchmark/results -name "shifts.json" -path "*/B3LYP_THF/*" | wc -l
# Should be 50 (or close to it)

# Generate summary CSV
python -m qm_nmr_calc.benchmark summary
# Creates: data/benchmark/results/summary.csv
```

**Source:** `src/qm_nmr_calc/benchmark/runner.py:489-536` (aggregate_results)

### Example shifts.json Output

After successful calculation, each result directory contains:

```json
{
  "molecule_id": "compound_01",
  "functional": "B3LYP",
  "solvent": "Pyridine",
  "h1_shifts": [],  // Empty until scaling factors exist (Phase 61)
  "c13_shifts": [], // Empty until scaling factors exist
  "shielding_data": {
    "index": [1, 2, 3, ...],
    "atom": ["H", "H", "C", ...],
    "shielding": [29.123, 28.456, 156.789, ...]
  }
}
```

**Note:** `h1_shifts` and `c13_shifts` will be empty because no scaling factors exist for pyridine/THF yet. That's expected - Phase 61 derives those factors. The critical data is `shielding_data` which is always populated on successful NWChem completion.

**Source:** `src/qm_nmr_calc/benchmark/runner.py:258-268` (shifts.json format)

## State of the Art

### Current Approach (v2.9)
- **DELTA50 benchmark protocol:** Use pre-optimized geometries, B3LYP/6-311+G(2d,p), COSMO solvation
- **Scaling factors:** OLS linear regression on calculated shieldings vs experimental shifts
- **Quality gates:** R² > 0.99, 1H MAE < 0.2 ppm, 13C MAE < 3.0 ppm

### Historical Context
| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Phase 11.2: CHCl3 + DMSO only | Phase 54-58: Added Methanol, Water, Acetone, Benzene | v2.8 (Feb 2026) | 200 successful calculations, 4 new solvent factors |
| Manual solvent addition | CLI-driven opt-in solvents | Phase 59 (v2.9, Feb 2026) | Can run any supported solvent without code changes |
| Hardcoded acetonitrile | COSMO_NAME_MAP layer | Phase 59 (v2.9, Feb 2026) | Handles NWChem's "acetntrl" abbreviation transparently |

### Phase 60 Position in Workflow

```
Phase 59: Extend CLI + COSMO mapping ✓ COMPLETE
    ↓
Phase 60: Run 100 calculations (Pyridine + THF) ← YOU ARE HERE
    ↓
Phase 61: Derive scaling factors from Phase 60 results
    ↓
Phase 62: Run 100 calculations (Toluene + DCM)
    ↓
Phase 63: Derive scaling factors from Phase 62 results
```

**Deprecation note:** None. Infrastructure is current and production-ready.

## Open Questions

### 1. Expected Failure Rate for Pyridine/THF
**What we know:** v2.8 achieved ~95% success rate for methanol/water/acetone/benzene (no reported systematic failures in Phase 55 verification).

**What's unclear:** Pyridine is aromatic nitrogen-containing, THF is cyclic ether. Both have been calculated before (not exotic), but specific failure modes with DELTA50 molecules unknown.

**Recommendation:**
- Accept 95-100% success (0-5 failures out of 100) as success criteria met
- Document any failures for debugging, but don't block Phase 60 completion
- If >5 failures: investigate whether issue is COSMO parameters (unlikely - verified in NWChem docs) or specific molecule geometries

### 2. Calculation Time Variance
**What we know:** v2.8 average calculation time was 400.2 seconds (from status.json). Range observed: 186-610 seconds per calculation.

**What's unclear:** Whether pyridine/THF will have similar timings. Different solvents have different dielectric constants, affecting SCF convergence.

**Recommendation:**
- Estimate: 100 calculations × 400s = ~11 hours total runtime
- Add 20% buffer for variance: ~13 hours
- Monitor first 5-10 calculations to validate timing assumptions
- If avg_calc_time_seconds in status.json exceeds 600s, investigate convergence issues

### 3. Scaling Factor Quality Prediction
**What we know:** v2.8 factors had R² > 0.995 for all 4 new solvents, MAE well within quality gates.

**What's unclear:** Whether pyridine/THF will maintain same quality (Phase 61 question, not Phase 60).

**Recommendation:**
- Phase 60 only needs successful calculations (shielding_data populated)
- Phase 61 (analysis) will reveal factor quality
- If Phase 61 factors fail quality gates, may need to investigate systematic errors in Phase 60 calculations

## Sources

### Primary (HIGH confidence)
- NWChem COSMO Solvation Model documentation - https://nwchemgit.github.io/COSMO-Solvation-Model.html (verified pyridine and THF in solvent parameter table)
- Project source code:
  - `src/qm_nmr_calc/benchmark/runner.py` (benchmark execution engine)
  - `src/qm_nmr_calc/benchmark/__main__.py` (CLI with Pyridine/THF in choices)
  - `src/qm_nmr_calc/nwchem/input_gen.py` (COSMO name mapping)
- Phase 59 verification report: `.planning/phases/59-benchmark-infrastructure/59-VERIFICATION.md` (confirmed Pyridine/THF support)
- Phase 56 verification report: `.planning/phases/56-scaling-factor-derivation/56-VERIFICATION.md` (confirmed v2.8 success metrics)

### Secondary (MEDIUM confidence)
- Recent Improvements to the NWChem COSMO Module - https://pubs.acs.org/doi/10.1021/acs.jctc.5c01368 (academic validation of COSMO implementation)

### Tertiary (LOW confidence)
None. All findings verified through code inspection or official documentation.

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - All libraries already in project, no new dependencies
- Architecture: HIGH - v2.8 ran 200 calculations successfully, proven infrastructure
- Execution procedure: HIGH - CLI tested in Phase 59, COSMO parameters verified in NWChem docs
- Timing estimates: MEDIUM - Based on v2.8 averages, but solvent-specific variance unknown
- Pitfalls: HIGH - Based on v2.8 execution experience and code analysis
- Expected failure modes: MEDIUM - Extrapolating from v2.8, no Phase 60-specific data yet

**Research date:** 2026-02-10
**Valid until:** 2026-03-10 (30 days - stable infrastructure, no fast-moving dependencies)

**Research methodology:**
1. Analyzed benchmark infrastructure source code (runner.py, __main__.py, input_gen.py)
2. Verified Phase 59 additions (Pyridine/THF in SUPPORTED_SOLVENTS and CLI choices)
3. Confirmed NWChem COSMO parameter support via official documentation
4. Reviewed Phase 56 verification for v2.8 success metrics and quality gates
5. Examined benchmark.log and status.json for timing and failure patterns
6. Cross-referenced Phase 59 VERIFICATION.md for solvent name mapping chain

**Key finding:** Phase 60 is pure execution - zero code changes required. All infrastructure proven in v2.8 (200 calculations). Phase 59 added solvent support. Success depends entirely on NWChem execution, not implementation correctness.
