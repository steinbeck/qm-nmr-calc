# Phase 9: Benchmark Calculations - Research

**Researched:** 2026-01-22
**Domain:** Long-running batch execution, headless process management, status monitoring
**Confidence:** HIGH

## Summary

Phase 9 executes 200 QM calculations (50 DELTA50 molecules x 2 functionals x 2 solvents) using the benchmark infrastructure built in Phase 8. The research focuses on three key areas: (1) headless execution to survive session disconnects, (2) status file patterns for monitoring from any session, and (3) graceful pause/stop mechanisms.

The existing `benchmark/runner.py` provides the calculation loop with resume support via `shifts.json` marker files. This phase extends it with detached execution, richer status monitoring, and human-readable reporting.

**Primary recommendation:** Use `nohup` for process detachment with a JSON status file that the runner updates after each calculation. Implement graceful stop via marker file check (`STOP` file) between calculations, not signal handlers (simpler and safer for long-running NWChem processes).

## Standard Stack

The established libraries/tools for this domain:

### Core (Already in Project)
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| tqdm | 4.67+ | Progress display during interactive runs | Already used in runner.py |
| orjson | 3.11+ | Fast JSON read/write with OPT_INDENT_2 | Project standard for status.json files |
| datetime | stdlib | Timestamps and duration calculation | Standard library |
| signal | stdlib | Signal handlers for graceful shutdown | Standard library |

### New for Headless Execution
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| nohup | system | Detach from terminal, survive SSH disconnect | Default for headless runs |
| pathlib | stdlib | File operations for status/stop marker files | Already project standard |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| nohup | tmux/screen | tmux allows reattach but adds dependency; nohup simpler for fire-and-forget |
| nohup | systemd service | systemd adds complexity for one-time batch run |
| marker file stop | SIGTERM handler | Signal handlers in Python subprocess chains can be tricky; marker file is simpler |

**Installation:**
No new dependencies needed - all tools are system utilities or already in pyproject.toml.

## Architecture Patterns

### Recommended Status File Structure

The status.json file should be written atomically after each calculation:

```
data/benchmark/results/
  status.json           # Global benchmark status (updated after each calc)
  progress.jsonl        # Existing append-only progress log
  compound_01/
    B3LYP_CHCl3/
      shifts.json       # Completion marker (existing)
      shielding.out
    B3LYP_DMSO/
      ...
  compound_02/
    ...
```

### Pattern 1: Status File Schema

**What:** JSON file tracking overall benchmark progress, readable by CLI or any tool
**When to use:** After every calculation completion, read on-demand by status command
**Example:**
```json
{
  "run_id": "2026-01-22T14:30:00",
  "started_at": "2026-01-22T14:30:00Z",
  "updated_at": "2026-01-22T15:45:00Z",
  "state": "running",
  "total_tasks": 200,
  "completed": 45,
  "failed": 2,
  "current_task": {
    "molecule_id": "compound_12",
    "functional": "B3LYP",
    "solvent": "CHCl3",
    "started_at": "2026-01-22T15:42:00Z"
  },
  "failures": [
    {
      "molecule_id": "compound_08",
      "functional": "WP04",
      "solvent": "DMSO",
      "error": "NWChem failed with exit code 1",
      "timestamp": "2026-01-22T15:20:00Z"
    }
  ],
  "estimated_remaining_hours": 12.5,
  "avg_calc_time_seconds": 180.5,
  "pid": 12345
}
```

### Pattern 2: Marker File Graceful Stop

**What:** Check for `STOP` marker file between calculations to pause cleanly
**When to use:** User wants to stop without killing mid-calculation
**Example:**
```python
STOP_FILE = results_dir / "STOP"

def should_stop() -> bool:
    """Check if graceful stop requested."""
    return STOP_FILE.exists()

# In main loop:
for task in tasks:
    if should_stop():
        logger.info("Stop requested, finishing current batch")
        status["state"] = "stopped"
        write_status()
        break

    run_single_calculation(task)
    update_status()
```

### Pattern 3: Atomic Status Write with orjson

**What:** Write status atomically to prevent corruption on read
**When to use:** Always when updating status.json
**Example:**
```python
import orjson
from pathlib import Path

def write_status(status: dict, status_file: Path) -> None:
    """Write status atomically to avoid read-during-write issues."""
    temp_file = status_file.with_suffix('.json.tmp')
    temp_file.write_bytes(
        orjson.dumps(status, option=orjson.OPT_INDENT_2)
    )
    temp_file.rename(status_file)  # Atomic on POSIX
```

### Pattern 4: Headless Execution Command

**What:** Launch benchmark as detached process
**When to use:** For the full 200-calculation run
**Example:**
```bash
# From project root:
nohup python -m qm_nmr_calc.benchmark run --processes 8 > benchmark.log 2>&1 &
echo $! > benchmark.pid

# Later, to check:
python -m qm_nmr_calc.benchmark status

# To stop gracefully:
touch data/benchmark/results/STOP
# Wait for current calculation to finish, then runner exits cleanly
```

### Pattern 5: Pilot Run Before Full Execution

**What:** Run first 5 molecules (20 calculations) to validate pipeline
**When to use:** Before committing to the full 200-calculation run
**Example:**
```bash
# Pilot run with first 5 molecules
python -m qm_nmr_calc.benchmark run \
  --molecules compound_01 compound_02 compound_03 compound_04 compound_05

# Verify results look reasonable
python -m qm_nmr_calc.benchmark status

# Check a few shifts.json files manually
cat data/benchmark/results/compound_01/B3LYP_CHCl3/shifts.json | python -m json.tool
```

### Anti-Patterns to Avoid

- **Signal handler complexity:** Don't try to catch SIGTERM in complex subprocess chains; NWChem processes complicate this. Use marker file instead.
- **Frequent status writes:** Don't write status.json on every loop iteration (e.g., inside tqdm); write once per completed calculation.
- **Tqdm in headless mode:** Tqdm progress bars are for interactive use; in nohup mode, use logging to file instead.
- **Polling status file too fast:** CLI status command should read once and exit, not poll in a loop.

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Process detachment | Custom daemon code | nohup | nohup handles SIGHUP correctly |
| Atomic file writes | Manual temp+rename | pathlib.rename() | rename() is atomic on POSIX |
| JSON serialization | json.dumps() | orjson.dumps() | Project standard, faster, OPT_INDENT_2 |
| Resume tracking | Custom state machine | Existing shifts.json marker | Already implemented in runner.py |
| Progress logging | Custom format | JSONL progress.jsonl | Already implemented in runner.py |

**Key insight:** The existing Phase 8 infrastructure handles most of the hard problems (resume, progress tracking, error handling). Phase 9 only needs to add: (1) richer status.json, (2) headless launch pattern, (3) graceful stop mechanism.

## Common Pitfalls

### Pitfall 1: Tqdm Breaks in Headless Mode

**What goes wrong:** Tqdm writes control characters to log files, making them hard to parse
**Why it happens:** Tqdm detects non-TTY but still outputs update lines
**How to avoid:** Disable tqdm or use `disable=True` when running headless; use logging for progress instead
**Warning signs:** Garbled output in benchmark.log file

### Pitfall 2: Status File Corruption on Read

**What goes wrong:** Reading status.json while runner is writing produces partial/invalid JSON
**Why it happens:** Non-atomic writes; reader catches mid-write
**How to avoid:** Atomic write pattern (temp file + rename); reader should handle JSONDecodeError gracefully
**Warning signs:** Intermittent "invalid JSON" errors in status command

### Pitfall 3: Forgetting to Log to File in Headless Mode

**What goes wrong:** All logging goes to /dev/null or gets lost
**Why it happens:** nohup redirects stdout but Python logging needs explicit file handler
**How to avoid:** Configure logging to file in benchmark runner, not just basicConfig
**Warning signs:** benchmark.log is empty or missing expected messages

### Pitfall 4: PID File Goes Stale

**What goes wrong:** benchmark.pid contains PID of finished process; misleading for monitoring
**Why it happens:** PID file not cleaned up on normal exit
**How to avoid:** Remove PID file on clean exit (in finally block); status command checks if PID is running
**Warning signs:** `ps -p $(cat benchmark.pid)` shows different process or nothing

### Pitfall 5: Failure Threshold Not Enforced

**What goes wrong:** Run continues with >5 failures (10% threshold) without pausing
**Why it happens:** No check for cumulative failure count in loop
**How to avoid:** Check `len(failures)` after each calculation; pause if > threshold
**Warning signs:** Summary shows many failures that could have been caught earlier

### Pitfall 6: Time Estimates Based on Early Calculations Only

**What goes wrong:** ETA wildly inaccurate because early (small) molecules skew average
**Why it happens:** Molecule size varies significantly in DELTA50 set
**How to avoid:** Use rolling average of last N calculations, not overall average
**Warning signs:** ETA says "2 hours" but run takes 20 hours

## Code Examples

Verified patterns from codebase analysis:

### Existing Runner Pattern (Phase 8)
```python
# Source: src/qm_nmr_calc/benchmark/runner.py
def run_benchmark(resume: bool = True, ...) -> list[BenchmarkResult]:
    results_dir = get_results_dir()
    progress_file = results_dir / "progress.jsonl"

    tasks = build_task_matrix(molecules, functionals, solvents)

    if resume:
        tasks = [t for t in tasks if not is_task_complete(results_dir, t)]

    results = []
    for task in tqdm(tasks, desc="DELTA50 Benchmark"):
        result = run_single_calculation(...)
        results.append(result)
        append_to_jsonl(progress_file, {...})

    return results
```

### Existing Status Pattern (Project Standard)
```python
# Source: src/qm_nmr_calc/storage.py
def _write_status(job_id: str, status: JobStatus) -> None:
    """Write status to disk as JSON."""
    status_file = get_job_dir(job_id) / "status.json"
    status_file.write_bytes(
        orjson.dumps(status.model_dump(mode="json"), option=orjson.OPT_INDENT_2)
    )
```

### Enhanced Status Write for Benchmark
```python
# New pattern for benchmark status.json
import orjson
from datetime import datetime, timezone
from pathlib import Path

def update_benchmark_status(
    status_file: Path,
    state: str,
    completed: int,
    failed: int,
    total: int,
    current_task: dict | None,
    failures: list[dict],
    started_at: datetime,
    avg_time: float | None,
) -> None:
    """Update global benchmark status file atomically."""
    now = datetime.now(timezone.utc)

    remaining = total - completed - failed
    eta_seconds = remaining * avg_time if avg_time else None
    eta_hours = eta_seconds / 3600 if eta_seconds else None

    status = {
        "state": state,
        "started_at": started_at.isoformat(),
        "updated_at": now.isoformat(),
        "total_tasks": total,
        "completed": completed,
        "failed": failed,
        "current_task": current_task,
        "failures": failures[-10:],  # Keep last 10 for debugging
        "estimated_remaining_hours": round(eta_hours, 2) if eta_hours else None,
        "avg_calc_time_seconds": round(avg_time, 1) if avg_time else None,
    }

    # Atomic write
    temp_file = status_file.with_suffix('.json.tmp')
    temp_file.write_bytes(orjson.dumps(status, option=orjson.OPT_INDENT_2))
    temp_file.rename(status_file)
```

### Graceful Stop Check
```python
# Pattern for checking STOP marker file
def check_stop_requested(results_dir: Path) -> bool:
    """Check if user requested graceful stop via marker file."""
    stop_file = results_dir / "STOP"
    return stop_file.exists()

def clear_stop_file(results_dir: Path) -> None:
    """Remove STOP file after acknowledging stop request."""
    stop_file = results_dir / "STOP"
    if stop_file.exists():
        stop_file.unlink()
```

### CLI Status Display
```python
# Pattern for human-readable status output
def format_status(status: dict) -> str:
    """Format status.json for human-readable display."""
    lines = [
        "DELTA50 Benchmark Status",
        "=" * 40,
        f"State: {status['state'].upper()}",
        f"Progress: {status['completed']}/{status['total_tasks']} "
        f"({100*status['completed']/status['total_tasks']:.1f}%)",
        f"Failed: {status['failed']}",
    ]

    if status.get("current_task"):
        ct = status["current_task"]
        lines.append(f"Current: {ct['molecule_id']} / {ct['functional']} / {ct['solvent']}")

    if status.get("estimated_remaining_hours"):
        hours = status["estimated_remaining_hours"]
        if hours < 1:
            lines.append(f"ETA: {int(hours * 60)} minutes")
        else:
            lines.append(f"ETA: {hours:.1f} hours")

    return "\n".join(lines)
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Screen/tmux sessions | nohup + status file | Best practice | Simpler for fire-and-forget batch runs |
| Signal handlers for stop | Marker file check | N/A | Safer with subprocess-heavy workloads |
| json.dumps | orjson with atomic write | Project standard | Faster, safer, consistent |

**Deprecated/outdated:**
- None - the existing Phase 8 patterns are current and well-designed

## Time Estimates

Based on the project's existing calculations and the DELTA50 benchmark literature:

| Molecule Size | Typical NMR Calc Time | Notes |
|--------------|----------------------|-------|
| Small (5-10 heavy atoms) | 2-5 minutes | compound_01-10 likely |
| Medium (10-20 heavy atoms) | 5-15 minutes | Most DELTA50 molecules |
| Larger (20+ heavy atoms) | 15-30 minutes | A few complex molecules |

**Estimated total time for 200 calculations:**
- Conservative: 200 x 10 min avg = 33 hours
- Optimistic: 200 x 5 min avg = 17 hours
- Realistic: 20-30 hours for full run

**Pilot run (20 calculations):** 1-3 hours - good validation before committing

## Open Questions

Things that couldn't be fully resolved:

1. **Exact NWChem timing for DELTA50 molecules**
   - What we know: Small molecules take ~5 min, larger ones ~15-30 min
   - What's unclear: Exact distribution of molecule sizes in DELTA50 set
   - Recommendation: Pilot run will provide accurate timing data

2. **Optimal MPI process count**
   - What we know: CLI default is 4 processes; more cores = faster but diminishing returns
   - What's unclear: Machine-specific optimal value (depends on available cores, memory)
   - Recommendation: Let user specify via `--processes` flag; document typical values (4-8)

## Sources

### Primary (HIGH confidence)
- Codebase analysis: `src/qm_nmr_calc/benchmark/runner.py`, `storage.py`, `models.py`
- Codebase analysis: Existing status.json patterns in job storage

### Secondary (MEDIUM confidence)
- [nohup vs screen vs tmux comparison](https://gist.github.com/MangaD/632e8f5a6649c9b2e30e2e5d3926447b) - Process detachment patterns
- [Python signal handling](https://docs.python.org/3/library/signal.html) - Signal handler reference
- [Signal handling best practices](https://johal.in/signal-handling-in-python-custom-handlers-for-graceful-shutdowns/) - Flag-based approach

### Tertiary (LOW confidence)
- [NWChem NMR benchmark study](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-018-0305-8) - ISiCLE/NWChem timing references (different basis sets than DELTA50)

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - Using existing project patterns and standard Unix tools
- Architecture: HIGH - Clear patterns from existing codebase
- Pitfalls: MEDIUM - Based on general best practices, not project-specific failures
- Time estimates: LOW - Need pilot run data to validate

**Research date:** 2026-01-22
**Valid until:** 2026-03-22 (60 days - stable domain, no external dependencies changing)
