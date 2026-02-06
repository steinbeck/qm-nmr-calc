# Phase 53: Conformer Progress Bug Fix - Research

**Researched:** 2026-02-06
**Domain:** Python FastAPI backend + JavaScript frontend polling
**Confidence:** HIGH

## Summary

The conformer progress tracking displays incorrectly during ensemble NMR calculations. The status bar shows "0/2" throughout processing instead of updating to "1/2", "2/2" as conformers complete. This is caused by a missing `update_job_status()` call in the progress callback that updates the step label but never persists the ensemble conformer status changes to disk.

The backend correctly updates conformer status fields in memory (changing from "pending" → "optimizing" → "optimized" → "nmr_running" → "nmr_complete"), but these changes are never written to `status.json` until the final step completion. The frontend polls `/api/v1/jobs/{job_id}` which reads from `status.json`, so it never sees the intermediate conformer status updates.

**Primary recommendation:** Add `update_job_status(job_id, conformer_ensemble=ensemble)` call inside the `on_progress()` callback in `tasks.py:run_ensemble_nmr_task()` to persist conformer status changes to disk after each conformer completes.

## Standard Stack

This is a bug fix in existing code. No new libraries required.

### Current Technology Stack

| Component | Technology | Purpose |
|-----------|-----------|---------|
| Backend | FastAPI + Pydantic | REST API with type-safe models |
| Task Queue | Huey | Background task processing |
| Storage | orjson + filesystem | Job status persistence as JSON |
| Frontend Polling | JavaScript fetch() | Status updates via polling |
| Progress Tracking | In-memory mutation + disk writes | Conformer status updates |

## Architecture Patterns

### Current Progress Tracking Architecture

```
┌─────────────────────────────────────────────────────────────┐
│ run_ensemble_nmr_task() (tasks.py:257-409)                 │
│                                                             │
│  1. Generate conformer ensemble                            │
│  2. update_job_status(conformer_ensemble=ensemble)  ←─── ✓│
│  3. Define on_progress() callback                          │
│  4. Call run_ensemble_dft_and_nmr()                        │
│     │                                                       │
│     └─→ run_conformer_dft_optimization()                   │
│         │  - Updates conformer.status in memory            │
│         │  - Calls progress_callback() after each conf     │
│         │                                                   │
│         └─→ on_progress(step, current, total)              │
│             │  start_step(job_id, step, label)             │
│             │  └─→ Writes current_step_label to disk       │
│             │                                               │
│             ✗ MISSING: update_job_status(                  │
│                   job_id,                                   │
│                   conformer_ensemble=ensemble              │
│                 )                                           │
│                                                             │
│  5. update_job_status(conformer_ensemble=ensemble)  ←─── ✓│
└─────────────────────────────────────────────────────────────┘

Frontend: GET /api/v1/jobs/{job_id} → reads status.json
```

**Key insight:** The `ensemble` object is mutated in-place by `run_conformer_dft_optimization()` and `run_conformer_nmr_calculations()`, but those mutations are never persisted to `status.json` during processing. The frontend only sees the final state when the task completes.

### Pattern: In-Place Mutation Without Persistence

**The Bug Pattern:**
```python
# runner.py:288-338 (run_conformer_dft_optimization)
for conformer in ensemble.conformers:
    conformer.status = "optimizing"  # In-memory mutation
    # ... DFT runs ...
    conformer.status = "optimized"   # In-memory mutation

    if progress_callback:
        progress_callback("optimizing_conformers", current, total)
        # ↑ Callback updates step label but NOT conformer_ensemble
```

**The on_progress callback:**
```python
# tasks.py:329-331
def on_progress(step: str, current: int, total: int):
    """Update job status with progress during conformer processing."""
    start_step(job_id, step, f"{step.replace('_', ' ').title()} ({current}/{total})")
    # ↑ Only updates current_step_label, NOT conformer_ensemble
```

**What's written to disk during processing:**
- ✓ `current_step`: "optimizing_conformers"
- ✓ `current_step_label`: "Optimizing conformers (1/2)"
- ✗ `conformer_ensemble`: Still shows all conformers as "pending"

### Pattern: Status Persistence Points

**Current persistence calls in `run_ensemble_nmr_task()`:**

| Line | What's Persisted | When |
|------|------------------|------|
| 326 | Initial ensemble with all conformers "pending" | After generation |
| 403 | Final ensemble with all statuses complete | After task finishes |

**Missing:** Incremental updates during DFT/NMR processing loops.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Status polling | Custom SSE/WebSocket | Existing 3-second polling | Already works, keep it simple |
| JSON serialization | Custom JSON encoder | orjson (already used) | Fast, handles datetime |
| Progress tracking | Message queue | File-based status.json | Simple, already deployed |

**Key insight:** This is not an architectural problem. The polling mechanism works fine. The bug is simply a missing persistence call.

## Common Pitfalls

### Pitfall 1: Mutating Without Persisting
**What goes wrong:** Object fields are updated in memory but never written to disk.

**Why it happens:** Python objects are mutable by reference. When `run_conformer_dft_optimization()` updates `conformer.status`, it modifies the same object that lives in the task's `ensemble` variable. But that mutation doesn't automatically trigger a disk write.

**How to avoid:** After any in-memory mutation that the frontend needs to see, explicitly call `update_job_status()` to persist changes.

**Warning signs:**
- Frontend shows stale data despite backend logs showing progress
- status.json content doesn't match in-memory state during debugging
- Updates only appear when task fully completes

### Pitfall 2: Over-Writing to Disk
**What goes wrong:** Writing status.json on every loop iteration can cause:
- I/O bottleneck (disk writes are slow)
- File corruption if writes overlap
- Unnecessary CPU for JSON serialization

**Why it happens:** Trying to fix Pitfall 1 by adding `update_job_status()` inside tight loops.

**How to avoid:**
- Write once per conformer completion (not per iteration)
- The current code structure already does this correctly (progress_callback called after each conformer)
- One write per conformer (2-50 conformers typical) is acceptable I/O

**Warning signs:**
- Disk I/O metrics spike during calculations
- Task execution slows down significantly
- Multiple processes trying to write simultaneously

### Pitfall 3: Stale Object References
**What goes wrong:** After calling `update_job_status()`, continuing to use the old `ensemble` object reference.

**Why it happens:** `update_job_status()` creates a new JobStatus object and writes it to disk, but the caller's local variable still points to the old object.

**How to avoid:**
- Option A: Reassign the result: `updated_status = update_job_status(...)`
- Option B: Only update ensemble field, not entire status: `update_job_status(job_id, conformer_ensemble=ensemble)`
- Option C: Mutate in place and persist without recreating (current pattern in runner.py)

**Current code uses Option C:** The `ensemble` object is passed by reference through multiple functions, mutated in place, then persisted with the same reference. This is fine.

## Root Cause Analysis

### The Bug

**File:** `src/qm_nmr_calc/tasks.py`
**Function:** `run_ensemble_nmr_task()`
**Lines:** 329-331

```python
def on_progress(step: str, current: int, total: int):
    """Update job status with progress during conformer processing."""
    start_step(job_id, step, f"{step.replace('_', ' ').title()} ({current}/{total})")
```

**What's missing:** Persisting the mutated `ensemble` object to disk.

**Fix:** Add one line:
```python
def on_progress(step: str, current: int, total: int):
    """Update job status with progress during conformer processing."""
    start_step(job_id, step, f"{step.replace('_', ' ').title()} ({current}/{total})")
    update_job_status(job_id, conformer_ensemble=ensemble)  # ← ADD THIS
```

### Why This Happens

**Call stack when conformer completes:**
1. `run_conformer_dft_optimization()` updates `conformer.status = "optimized"` (line 325)
2. Calls `progress_callback("optimizing_conformers", 1, 2)` (line 338)
3. `on_progress()` in tasks.py receives the callback (line 329)
4. Calls `start_step()` which writes `current_step_label` to disk (line 331)
5. ❌ Never writes the updated `ensemble.conformers[0].status` to disk

**Result:** Frontend polls `/api/v1/jobs/{job_id}`, reads `status.json`, sees:
- ✓ `current_step_label`: "Optimizing conformers (1/2)"
- ❌ `conformer_progress[0].status`: "pending" (should be "optimized")

### Frontend Display Logic

**File:** `src/qm_nmr_calc/api/templates/status.html`
**Function:** `updateConformerProgress()`
**Lines:** 299-362

```javascript
function updateConformerProgress(data) {
    if (data.conformer_progress && data.conformer_progress.length > 0) {
        const total = data.conformer_progress.length;

        // Count conformers by status
        const optimized = data.conformer_progress.filter(c =>
            c.status === 'optimized'
        ).length;  // ← Always 0 because statuses never update

        const nmrComplete = data.conformer_progress.filter(c =>
            c.status === 'nmr_complete'
        ).length;  // ← Always 0 until job finishes

        // Calculate completed count based on current step
        let completed;
        if (data.current_step === 'calculating_nmr') {
            completed = nmrComplete;  // Always 0
        } else {
            completed = optimized + nmrComplete;  // Always 0 + 0
        }

        conformerCount.textContent = completed + '/' + total;
        // Displays "0/2" throughout processing ❌
    }
}
```

**The frontend logic is correct.** It properly filters by status and counts completed conformers. The bug is purely that the backend never sends updated statuses.

## Code Examples

### Current (Buggy) Implementation

**runner.py:336-338**
```python
successful.append(conformer)

# Report progress after each conformer completes (success)
if progress_callback:
    progress_callback("optimizing_conformers", len(successful) + len(failed), len(ensemble.conformers))
```

**tasks.py:329-331**
```python
def on_progress(step: str, current: int, total: int):
    """Update job status with progress during conformer processing."""
    start_step(job_id, step, f"{step.replace('_', ' ').title()} ({current}/{total})")
    # Missing: update_job_status(job_id, conformer_ensemble=ensemble)
```

### Fixed Implementation

**tasks.py:329-332** (add one line)
```python
def on_progress(step: str, current: int, total: int):
    """Update job status with progress during conformer processing."""
    start_step(job_id, step, f"{step.replace('_', ' ').title()} ({current}/{total})")
    update_job_status(job_id, conformer_ensemble=ensemble)  # ← ADD THIS
```

**Why this works:**
1. `start_step()` updates `current_step_label` and writes to disk
2. `update_job_status()` updates `conformer_ensemble` and writes to disk
3. Frontend polls and receives both updated label AND updated conformer statuses
4. `updateConformerProgress()` counts updated statuses correctly
5. Display shows "1/2", "2/2" as conformers complete ✓

### Alternative: Update Inside Progress Callback Call

Another option is to make `progress_callback` accept the ensemble and update it:

**runner.py:336-338**
```python
successful.append(conformer)

# Report progress after each conformer completes (success)
if progress_callback:
    progress_callback("optimizing_conformers", len(successful) + len(failed), len(ensemble.conformers), ensemble)
```

**Verdict:** Don't do this. The callback shouldn't know about the ensemble object. The current architecture (pass counts, persist in tasks.py) is cleaner.

## State of the Art

| Old Approach | Current Approach | Impact |
|--------------|------------------|--------|
| Single-conformer only | Ensemble conformer mode (v2.0) | Progress tracking needed for multi-conformer |
| No progress during processing | Step labels update ("1/2") | Users see some progress, but not conformer completion |
| Immediate completion | 30-90 minute multi-conformer jobs | Users need accurate progress feedback |

**This bug exists because:** Ensemble mode (v2.0) added conformer-level status tracking, but the progress callback only updates step labels, not conformer statuses.

## Recommended Fix

### Change Required

**File:** `src/qm_nmr_calc/tasks.py`
**Function:** `run_ensemble_nmr_task()`
**Line:** 331 (add one line after `start_step()`)

```python
def on_progress(step: str, current: int, total: int):
    """Update job status with progress during conformer processing."""
    start_step(job_id, step, f"{step.replace('_', ' ').title()} ({current}/{total})")
    update_job_status(job_id, conformer_ensemble=ensemble)
```

### Why This Is The Right Fix

1. **Minimal change:** One line addition, no refactoring
2. **Correct scope:** Progress callback has access to `ensemble` via closure
3. **Appropriate frequency:** Called once per conformer (not per iteration)
4. **Consistent pattern:** Other parts of the function already call `update_job_status(job_id, conformer_ensemble=ensemble)` (lines 326, 403)
5. **No side effects:** `update_job_status()` is idempotent and safe to call multiple times

### Testing Strategy

**Verification steps:**

1. Start an ensemble NMR calculation with 2+ conformers
2. Poll `/api/v1/jobs/{job_id}` while job is running
3. Check `conformer_progress` field in response
4. Confirm statuses update: "pending" → "optimizing" → "optimized" → "nmr_running" → "nmr_complete"
5. Check status page displays "1/2", "2/2" as conformers complete
6. Verify progress bar moves from 0% to 50% to 100%

**Edge cases to test:**

- Single conformer (should still work, just updates once)
- Conformer failure (status should show "failed", count should not increase)
- Fast completion (rapid status updates don't cause corruption)

## Open Questions

None. The root cause is confirmed and the fix is straightforward.

## Sources

### Primary (HIGH confidence)
- Source code analysis of `src/qm_nmr_calc/tasks.py:257-409`
- Source code analysis of `src/qm_nmr_calc/nwchem/runner.py:246-578`
- Source code analysis of `src/qm_nmr_calc/storage.py:106-116` (`update_job_status()`)
- Source code analysis of `src/qm_nmr_calc/api/templates/status.html:299-362` (frontend display logic)
- Source code analysis of `src/qm_nmr_calc/api/routers/jobs.py:99-118` (API endpoint conformer_progress building)
- User-reported bug description in `.planning/post-v2.6-problems.md:23-27`

### Secondary (MEDIUM confidence)
- None required (direct code inspection)

## Metadata

**Confidence breakdown:**
- Root cause identification: HIGH - Direct code inspection confirms missing persistence call
- Fix approach: HIGH - One-line fix in well-understood code path
- Side effects: HIGH - No breaking changes, idempotent operation
- Frontend compatibility: HIGH - No frontend changes required, existing logic handles updated statuses

**Research date:** 2026-02-06
**Valid until:** 90+ days (bug fix research, not framework-dependent)
