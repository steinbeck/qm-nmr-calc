# Phase 17: API Integration and Progress Tracking - Research

**Researched:** 2026-01-28
**Domain:** FastAPI API design, Web UI progress tracking, 3Dmol.js conformer visualization
**Confidence:** HIGH

## Summary

This phase extends the existing qm-nmr-calc API and web UI to expose v2.0 ensemble functionality to users. The foundation work (Phases 12-16) has already implemented the computational pipeline: conformer generation, DFT optimization, NMR calculation, and Boltzmann averaging. This phase wires that pipeline to user-facing endpoints and templates.

The primary technical domains are: (1) API schema extension to support ensemble parameters and results, (2) Huey task orchestration for ensemble jobs, (3) progress tracking with per-conformer status updates, (4) web form updates for mode/method selection, and (5) 3Dmol.js conformer selection in the results viewer.

**Primary recommendation:** Extend existing patterns with minimal new concepts. The codebase has established conventions (Pydantic schemas, Huey tasks, Jinja templates, 3Dmol.js viewer) that should be followed rather than introducing new frameworks.

## Standard Stack

The phase uses existing project dependencies with no new libraries required.

### Core (Already in Project)
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| FastAPI | existing | API endpoints | Already handles /jobs, /health |
| Pydantic | existing | Request/response schemas | JobStatusResponse, NMRResultsResponse exist |
| Huey | existing | Background task queue | run_nmr_task pattern established |
| Jinja2 | existing | HTML templating | templates/*.html structure in place |
| 3Dmol.js | 2.5.3 (CDN) | 3D molecule viewer | Already integrated via base.html |
| Pico CSS | 2.x | Styling | Class-less semantic HTML patterns |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| orjson | existing | JSON serialization | Fast status file writes |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Polling | WebSocket | WebSocket adds complexity; polling already works at 3-5s intervals |
| Server-side ETA | Client-side ETA | Server has more context; client-side avoids extra state |

**Installation:**
No new packages required.

## Architecture Patterns

### Recommended Changes

```
src/qm_nmr_calc/
  api/
    schemas.py          # Extend JobSubmitRequest, JobStatusResponse, NMRResultsResponse
    routers/
      jobs.py           # Extend submit_smiles, add geometry endpoints for conformers
      web.py            # Extend form handling for ensemble parameters
  tasks.py              # Add run_ensemble_nmr_task, modify run_nmr_task dispatch
  models.py             # Already has ConformerEnsemble, ConformerData (Phase 12)
  templates/
    submit.html         # Add mode/method dropdowns, advanced options section
    status.html         # Add conformer progress section
    results.html        # Add conformer selector, per-conformer shift table
```

### Pattern 1: Mode-Aware Task Dispatch

**What:** Single submission endpoint dispatches to different Huey tasks based on conformer_mode
**When to use:** Job submission when user selects single vs ensemble mode
**Example:**
```python
# In routers/jobs.py submit_smiles endpoint
if job_status.conformer_mode == "ensemble":
    run_ensemble_nmr_task(job_status.job_id)
else:
    run_nmr_task(job_status.job_id)  # Existing v1.x task
```

### Pattern 2: Hierarchical Step Progress

**What:** Status updates show both pipeline stage and per-conformer progress
**When to use:** Ensemble job status polling
**Example:**
```python
# Extended step labels for ensemble mode
start_step(job_id, "generating_conformers", "Generating conformers")
start_step(job_id, "optimizing_conformers", f"Optimizing conformers (1/{total})")
start_step(job_id, "calculating_nmr", f"Calculating NMR ({completed}/{total})")
start_step(job_id, "averaging_shifts", "Computing Boltzmann average")
```

### Pattern 3: Extended Status Response with Ensemble Metadata

**What:** JobStatusResponse includes ensemble-specific fields when mode=ensemble
**When to use:** Status polling for ensemble jobs
**Example:**
```python
# Additional fields in JobStatusResponse for ensemble mode
class JobStatusResponse(BaseModel):
    # ... existing fields ...
    conformer_mode: str = "single"
    # Ensemble-specific (only present when mode="ensemble")
    conformer_method: Optional[str] = None  # "rdkit_kdg" or "crest"
    conformer_count: Optional[int] = None
    conformer_progress: Optional[list[dict]] = None  # [{id, status}, ...]
    eta_seconds: Optional[int] = None
```

### Pattern 4: Conformer-Aware Geometry Endpoint

**What:** geometry.json endpoint returns multiple geometries for conformer selection
**When to use:** 3D viewer needs to switch between conformers
**Example:**
```python
# Extended /api/v1/jobs/{job_id}/geometry.json
{
    "job_id": "abc123",
    "status": "complete",
    "conformer_mode": "ensemble",
    "conformers": [
        {
            "id": "conf_001",
            "xyz": "...",
            "sdf": "...",
            "energy_kcal": 0.0,
            "population": 0.452,
        },
        ...
    ],
    "h1_assignments": {...},  # Average shifts
    "c13_assignments": {...},
    "per_conformer_shifts": [...]  # Per-conformer detail if needed
}
```

### Anti-Patterns to Avoid
- **Duplicating pipeline logic in tasks.py:** Use run_ensemble_dft_and_nmr from runner.py
- **Blocking conformer loops in API endpoints:** All processing happens in Huey tasks
- **Breaking v1.x compatibility:** Single-conformer mode must work exactly as before
- **Hardcoding CREST availability:** Check health endpoint or detect_crest_available()

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Boltzmann averaging | Custom averaging | `average_ensemble_nmr()` | Already handles edge cases, weight population |
| Conformer filtering | Manual filtering | `apply_post_dft_filter()` | Energy window logic is tricky |
| CREST detection | Inline `which` calls | `detect_crest_available()` | Cached, checks both binaries |
| Step timing | Manual datetime math | `start_step()`, `complete_current_step()` | Existing step tracking system |
| XYZ file handling | String parsing | `load_geometry_file()` | RDKit-backed, validates atoms |

**Key insight:** The computational pipeline (Phases 12-16) is complete. This phase integrates it without modifying the pipeline logic.

## Common Pitfalls

### Pitfall 1: Ensemble/NMR Results List Length Mismatch

**What goes wrong:** `average_ensemble_nmr()` requires `len(conformers) == len(nmr_results)`. After pipeline, ensemble contains ALL conformers (including filtered/failed), but nmr_results only covers successful ones.

**Why it happens:** Some conformers are filtered by post-DFT energy window (status="optimized") or fail DFT/NMR (status="failed").

**How to avoid:** In the Huey task, filter ensemble.conformers to only those with status="nmr_complete" before calling `average_ensemble_nmr`:
```python
nmr_complete_conformers = [c for c in ensemble.conformers if c.status == "nmr_complete"]
# Now len(nmr_complete_conformers) == len(nmr_results)
filtered_ensemble = ensemble.model_copy(update={"conformers": nmr_complete_conformers})
avg_nmr = average_ensemble_nmr(filtered_ensemble, nmr_results)
```

**Warning signs:** IndexError or ValueError from Boltzmann averaging functions.

### Pitfall 2: CREST Fallback Without Warning

**What goes wrong:** User selects CREST method, but CREST isn't installed, and job silently falls back to RDKit.

**Why it happens:** Lenient validation at API layer.

**How to avoid:** Implement the CONTEXT.md decision "Warn and fallback": proceed with RDKit but include warning in job metadata/results. Store `conformer_method_warning: str` in job status.

**Warning signs:** User sees "method: crest" in request but results show RDKit conformer count.

### Pitfall 3: Frontend Expects Immediate Conformer Data

**What goes wrong:** Status.html JavaScript expects conformer progress data immediately, but it's not populated until task starts.

**Why it happens:** Template renders before task begins conformer generation.

**How to avoid:** Check for null/empty before rendering progress. Show "Pending..." state initially.

**Warning signs:** JavaScript errors about undefined properties.

### Pitfall 4: ETA Calculation Without Historical Data

**What goes wrong:** ETA shown as "N/A" because no timing data exists yet.

**Why it happens:** First ensemble job has no prior timings to extrapolate from.

**How to avoid:** Use heuristics based on current step:
- generating_conformers: ~30s for RDKit, ~minutes for CREST
- optimizing_conformers: ~5 min per conformer (varies by molecule size)
- calculating_nmr: ~3 min per conformer
Show "Estimating..." for first job, refine with actual step durations.

**Warning signs:** ETA jumps erratically or shows impossible values.

### Pitfall 5: 3Dmol.js Model/Frame Confusion

**What goes wrong:** Using addModel() for each conformer creates separate models that can't be easily switched.

**Why it happens:** 3Dmol.js has two approaches: frames (animation) vs models (independent).

**How to avoid:** For conformer switching (not animation), use separate models with `removeAllModels()` before adding new one, OR load multi-model file with `addModelsAsFrames` and use `setFrame(index)`.

For this project, recommendation: Load individual XYZ/SDF per conformer (simpler, each has its own geometry.json response).

**Warning signs:** All conformers visible simultaneously, or switching doesn't work.

## Code Examples

Verified patterns from existing codebase:

### Ensemble Task Structure
```python
# Source: tasks.py pattern + nwchem/runner.py
@huey.task()
def run_ensemble_nmr_task(job_id: str) -> dict:
    """Execute ensemble NMR calculation."""
    job_status = load_job_status(job_id)
    preset = PRESETS[PresetName(job_status.input.preset)]

    # Step 1: Generate conformer ensemble
    start_step(job_id, "generating_conformers", "Generating conformers")
    ensemble = generate_conformer_ensemble(
        smiles=job_status.input.smiles,
        job_id=job_id,
        conformer_method=job_status.input.conformer_method or "rdkit_kdg",
        solvent=job_status.input.solvent,
    )
    update_job_status(job_id, conformer_ensemble=ensemble)

    # Step 2-4: DFT -> NMR -> Boltzmann (uses run_ensemble_dft_and_nmr)
    # ... with step tracking updates ...

    # Step 5: Average results
    start_step(job_id, "averaging_shifts", "Computing Boltzmann average")
    nmr_complete = [c for c in ensemble.conformers if c.status == "nmr_complete"]
    filtered_ensemble = ensemble.model_copy(update={"conformers": nmr_complete})
    avg_nmr = average_ensemble_nmr(filtered_ensemble, nmr_results)

    update_job_status(job_id, nmr_results=avg_nmr, conformer_ensemble=ensemble)
```

### Extended JobSubmitRequest
```python
# Source: schemas.py existing pattern
class JobSubmitRequest(BaseModel):
    smiles: str
    # ... existing fields ...
    conformer_mode: str = Field(
        default="ensemble",  # CONTEXT.md: default to ensemble
        description="Conformational sampling mode: 'single' or 'ensemble'",
    )
    conformer_method: Optional[str] = Field(
        None,  # None means use default (rdkit_kdg)
        description="Conformer generation method: 'rdkit_kdg' or 'crest'",
    )
```

### Form with Visible Dropdowns and Advanced Section
```html
<!-- Source: Pico CSS patterns + CONTEXT.md decisions -->
<fieldset>
    <legend>Conformer Mode</legend>

    <label for="conformer_mode">Calculation Mode</label>
    <select id="conformer_mode" name="conformer_mode">
        <option value="ensemble" selected>Ensemble (multi-conformer)</option>
        <option value="single">Single conformer (v1.x behavior)</option>
    </select>

    <label for="conformer_method">Conformer Method</label>
    <select id="conformer_method" name="conformer_method">
        <option value="rdkit_kdg" selected>RDKit KDG (fast)</option>
        <option value="crest" {% if not crest_available %}disabled{% endif %}>
            CREST/xTB (thorough){% if not crest_available %} - not available{% endif %}
        </option>
    </select>
</fieldset>

<details>
    <summary>Advanced Options</summary>
    <fieldset>
        <label for="pre_dft_energy_window">Pre-DFT Energy Window (kcal/mol)</label>
        <input type="number" id="pre_dft_energy_window" name="pre_dft_energy_window"
               value="6.0" min="1" max="20" step="0.5">
        <!-- ... more advanced options ... -->
    </fieldset>
</details>
```

### Conformer Progress in Status Template
```javascript
// Source: status.html existing pattern + extension
function updateConformerProgress(data) {
    if (data.conformer_mode !== 'ensemble') return;

    const progressDiv = document.getElementById('conformer-progress');
    progressDiv.style.display = 'block';

    if (data.conformer_progress) {
        const completed = data.conformer_progress.filter(c =>
            c.status === 'nmr_complete' || c.status === 'optimized'
        ).length;
        const total = data.conformer_progress.length;

        progressDiv.innerHTML = `
            <p>Conformers: ${completed}/${total}</p>
            ${data.eta_seconds ? `<p>ETA: ~${formatEta(data.eta_seconds)}</p>` : ''}
        `;
    }
}
```

### Conformer Selector in Results Viewer
```javascript
// Source: 3Dmol.js docs + results.html pattern
let currentConformerIndex = 0;
let conformerData = null;

async function loadConformerData() {
    const response = await fetch('/api/v1/jobs/' + JOB_ID + '/geometry.json');
    conformerData = await response.json();

    if (conformerData.conformer_mode === 'ensemble' && conformerData.conformers) {
        populateConformerSelector(conformerData.conformers);
    }
    displayConformer(0);
}

function populateConformerSelector(conformers) {
    const selector = document.getElementById('conformer-selector');
    selector.innerHTML = conformers.map((c, i) =>
        `<option value="${i}">${c.id}: ${c.energy_kcal.toFixed(2)} kcal/mol (${(c.population * 100).toFixed(1)}%)</option>`
    ).join('');
    selector.addEventListener('change', e => displayConformer(parseInt(e.target.value)));
}

function displayConformer(index) {
    const conf = conformerData.conformers[index];
    viewer.removeAllModels();
    viewer.removeAllLabels();

    if (conf.sdf) {
        viewer.addModel(conf.sdf, 'sdf');
    } else {
        viewer.addModel(conf.xyz, 'xyz', {assignBonds: true});
    }

    viewer.setStyle({}, {stick: {radius: 0.12}, sphere: {scale: 0.25}});
    addShiftLabels(viewer.getModel(), conformerData.h1_assignments, conformerData.c13_assignments);
    viewer.zoomTo();
    viewer.render();
}
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Single-conformer only | Ensemble with Boltzmann | v2.0 (current) | Much more accurate for flexible molecules |
| ETKDGv3 embedding | KDG embedding | v2.0 | Avoids crystal structure bias |
| No CREST option | Optional CREST | v2.0 | Better conformer ranking for complex molecules |

**Deprecated/outdated:**
- None relevant to this phase. All foundation code (12-16) is current.

## Open Questions

Things that couldn't be fully resolved:

1. **ETA Algorithm Precision**
   - What we know: Can calculate based on step durations
   - What's unclear: Molecule-size scaling factors (larger molecules take proportionally longer)
   - Recommendation: Start with linear extrapolation, refine with logged data. Show "Estimated" prefix.

2. **Per-Conformer Shifts Display Layout**
   - What we know: CONTEXT.md says "side-by-side table"
   - What's unclear: How many conformers fit horizontally (5? 10? scrollable?)
   - Recommendation: Limit to top 5 by population, horizontal scroll for more

3. **Retry Pre-fill Parameters**
   - What we know: CONTEXT.md says "pre-filled form" for retries
   - What's unclear: How to pass failed job's parameters to new submission
   - Recommendation: Query string params or localStorage; simplest is link with query params

## Sources

### Primary (HIGH confidence)
- Project codebase: schemas.py, tasks.py, runner.py, templates/*.html
- Project STATE.md and 17-CONTEXT.md for decisions
- 3Dmol.js GLViewer documentation: https://3dmol.csb.pitt.edu/doc/GLViewer.html
- Pico CSS docs: https://picocss.com/docs/dropdown, https://picocss.com/docs/forms/select

### Secondary (MEDIUM confidence)
- FastAPI polling patterns: https://medium.com/@bhagyarana80/serving-long-running-jobs-with-fastapi-using-webhooks-and-task-polling-860bb0d3e0f9
- Progress tracking UX: https://userguiding.com/blog/progress-trackers-and-indicators

### Tertiary (LOW confidence)
- ETA calculation heuristics (based on general best practices, no chemistry-specific sources)

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - Using existing project libraries
- Architecture: HIGH - Extending established patterns
- Pitfalls: HIGH - Based on STATE.md blockers and codebase analysis
- Code examples: HIGH - Derived from existing codebase patterns
- ETA algorithm: MEDIUM - Heuristic approach, needs refinement

**Research date:** 2026-01-28
**Valid until:** 2026-02-28 (30 days - stable domain, no expected library changes)
