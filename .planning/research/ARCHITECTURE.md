# Architecture Research: Conformational Sampling for NMR

**Domain:** Conformational Sampling for NMR Predictions
**Researched:** 2026-01-26
**Confidence:** HIGH

## Standard Architecture Pattern

### Current Single-Conformer Architecture

```
User Request (SMILES)
         ↓
   FastAPI Endpoint
         ↓
   Create Job Directory (data/jobs/{job_id}/)
         ├─ status.json
         ├─ output/
         └─ logs/
         ↓
   Queue Huey Task (run_nmr_task)
         ↓
   Huey Worker (single process)
         ├─ RDKit 3D embed (ETKDG, deterministic seed)
         ├─ NWChem geometry optimization (B3LYP + COSMO)
         ├─ NWChem NMR shielding (B3LYP + COSMO)
         ├─ Parse shielding output
         ├─ Apply DELTA50 scaling factors
         └─ Generate visualizations
         ↓
   Update status.json (complete)
         ↓
   User polls /jobs/{job_id} for results
```

### Proposed Ensemble Architecture

```
User Request (SMILES, conformer_mode='ensemble')
         ↓
   FastAPI Endpoint
         ↓
   Create Job Directory (data/jobs/{job_id}/)
         ├─ status.json (ensemble_mode=True, conformer_count=N)
         ├─ output/
         │    ├─ initial.xyz (lowest-energy RDKit/CREST conformer)
         │    ├─ conformers/
         │    │    ├─ conf_0001.xyz (input)
         │    │    ├─ conf_0002.xyz (input)
         │    │    └─ ...
         │    ├─ optimized/
         │    │    ├─ conf_0001.xyz (DFT-optimized)
         │    │    ├─ conf_0002.xyz (DFT-optimized)
         │    │    └─ ...
         │    └─ optimized.xyz (lowest-energy conformer for visualization)
         ├─ logs/
         └─ scratch/
              ├─ conf_0001/
              │    ├─ optimize.nw
              │    ├─ optimize.out
              │    ├─ shielding.nw
              │    └─ shielding.out
              ├─ conf_0002/
              └─ ...
         ↓
   Queue Huey Task (run_ensemble_nmr_task)
         ↓
   Huey Worker (single process, sequential conformer processing)
         │
         ├─ [STEP 1: Conformer Generation]
         │    ├─ If mode='rdkit': RDKit ETKDG/KDG → MMFF opt → save to conformers/
         │    └─ If mode='crest': CREST/xTB → save to conformers/
         │    Update: status.json (conformer_count=N, current_step='conformer_generation')
         │
         ├─ [STEP 2: Energy Filtering (Pre-DFT)]
         │    ├─ Filter by MMFF/GFN2-xTB energy window (default: 6 kcal/mol)
         │    └─ Keep M conformers
         │    Update: status.json (conformers_filtered=M)
         │
         ├─ [STEP 3: DFT Optimization Loop]
         │    For each conformer i in 1..M (sequential):
         │         ├─ NWChem geometry optimization (B3LYP + COSMO)
         │         ├─ Parse DFT energy from output
         │         ├─ Save optimized geometry to optimized/conf_{i:04d}.xyz
         │         └─ Update: status.json (conformers_optimized=i/M)
         │
         ├─ [STEP 4: Energy Filtering (Post-DFT)]
         │    ├─ Rank by DFT energy
         │    ├─ Filter by DFT energy window (default: 3 kcal/mol)
         │    └─ Keep P conformers for NMR
         │    Update: status.json (conformers_for_nmr=P)
         │
         ├─ [STEP 5: NMR Shielding Loop]
         │    For each conformer i in surviving set (sequential):
         │         ├─ NWChem NMR shielding (B3LYP + COSMO)
         │         ├─ Parse shielding output
         │         └─ Update: status.json (conformers_nmr_complete=i/P)
         │
         ├─ [STEP 6: Boltzmann Averaging]
         │    ├─ Calculate Boltzmann weights from DFT energies
         │    ├─ Compute weighted-average shieldings per atom
         │    ├─ Apply DELTA50 scaling factors
         │    └─ Generate averaged NMR results
         │
         └─ [STEP 7: Visualization]
              ├─ Generate spectrum plots (1H, 13C)
              ├─ Generate annotated structure (lowest-energy conformer)
              └─ Save conformer_ensemble.json (all conformer data)
         ↓
   Update status.json (complete, nmr_results=averaged_shifts)
         ↓
   User polls /jobs/{job_id} for results
```

## Component Responsibilities

### Existing Components

| Component | Current Responsibility | Changes for Ensemble |
|-----------|----------------------|---------------------|
| `tasks.py` | Single Huey task `run_nmr_task` | Add `run_ensemble_nmr_task` with conformer loop |
| `storage.py` | Job directory creation, status updates | Add ensemble-specific directory structure helpers |
| `models.py` | `JobStatus`, `JobInput`, `NMRResults` | Add `ConformerData`, extend `JobStatus` with ensemble fields |
| `api/schemas.py` | API request/response models | Add `conformer_mode` field to `JobSubmitRequest` |
| `nwchem/runner.py` | Single calc orchestration (`run_calculation`) | Minor changes: accept conformer_id for file naming |
| `shifts.py` | Shielding → shift conversion | No changes (works per-atom, reusable) |

### New Components

| Component | Responsibility | Implementation Notes |
|-----------|---------------|---------------------|
| `conformers.py` | Conformer generation (RDKit, CREST) | Two functions: `generate_rdkit()`, `generate_crest()` |
| `averaging.py` | Boltzmann weighting and averaging | `boltzmann_weights()`, `average_shieldings()` |
| `ensemble.py` | Ensemble job orchestration logic | Helper functions for filtering, energy parsing, etc. |

## Data Flow

### Ensemble Job Flow (Detailed)

```python
# STEP 1: API receives request
POST /jobs
{
    "smiles": "CCCC",
    "solvent": "chcl3",
    "preset": "production",
    "conformer_mode": "ensemble",  # NEW FIELD
    "conformer_method": "rdkit",   # NEW FIELD (optional, default: "rdkit")
    "max_conformers": 50,          # NEW FIELD (optional, default: 50)
    "energy_window_pre": 6.0,      # NEW FIELD (optional, default: 6.0 kcal/mol)
    "energy_window_post": 3.0      # NEW FIELD (optional, default: 3.0 kcal/mol)
}

# STEP 2: Create job directory with ensemble metadata
storage.create_job_directory() → JobStatus
    status.json:
    {
        "job_id": "abc123def456",
        "status": "queued",
        "input": {
            "conformer_mode": "ensemble",
            "conformer_method": "rdkit",
            "max_conformers": 50,
            "energy_window_pre": 6.0,
            "energy_window_post": 3.0
        },
        "ensemble_data": {
            "conformers_generated": 0,
            "conformers_filtered_pre": 0,
            "conformers_optimized": 0,
            "conformers_filtered_post": 0,
            "conformers_nmr_complete": 0,
            "conformer_energies": [],
            "conformer_weights": []
        }
    }

# STEP 3: Queue Huey task
huey.enqueue(run_ensemble_nmr_task(job_id="abc123def456"))

# STEP 4: Huey worker executes task
def run_ensemble_nmr_task(job_id: str) -> dict:
    job_status = load_job_status(job_id)
    job_dir = get_job_dir(job_id)

    # STEP 4a: Generate conformers
    start_step(job_id, "conformer_generation", "Generating conformers")

    if job_status.input.conformer_method == "rdkit":
        conformers = generate_rdkit_conformers(
            smiles=job_status.input.smiles,
            num_confs=job_status.input.max_conformers,
            energy_window=job_status.input.energy_window_pre
        )
        # conformers = [(conf_id, energy, xyz_block), ...]
    elif job_status.input.conformer_method == "crest":
        conformers = generate_crest_conformers(
            smiles=job_status.input.smiles,
            solvent=job_status.input.solvent,
            energy_window=job_status.input.energy_window_pre
        )

    # Save conformers to disk
    conformers_dir = job_dir / "output" / "conformers"
    conformers_dir.mkdir(exist_ok=True)
    for i, (conf_id, energy, xyz) in enumerate(conformers):
        conf_file = conformers_dir / f"conf_{i+1:04d}.xyz"
        conf_file.write_text(xyz)

    update_job_status(job_id, ensemble_data={
        "conformers_generated": len(conformers),
        "conformers_filtered_pre": len(conformers)
    })

    # STEP 4b: DFT optimization loop
    start_step(job_id, "dft_optimization", f"Optimizing {len(conformers)} conformers")

    dft_results = []
    for i, (conf_id, mmff_energy, xyz) in enumerate(conformers):
        # Create conformer-specific scratch directory
        conf_scratch = job_dir / "scratch" / f"conf_{i+1:04d}"
        conf_scratch.mkdir(parents=True, exist_ok=True)

        # Run NWChem optimization
        result = run_calculation(
            smiles=None,  # Already have XYZ
            job_dir=conf_scratch,
            preset=preset,
            solvent=job_status.input.solvent,
            skip_optimization=False,
            geometry_input=xyz,
            skip_nmr=True,  # Only optimize, skip NMR for now
        )

        # Parse DFT energy from output
        dft_energy = parse_dft_energy(result['optimization_output'])

        # Save optimized geometry
        opt_dir = job_dir / "output" / "optimized"
        opt_dir.mkdir(exist_ok=True)
        opt_file = opt_dir / f"conf_{i+1:04d}.xyz"
        opt_file.write_text(result['geometry_file'].read_text())

        dft_results.append({
            "conf_id": i+1,
            "mmff_energy": mmff_energy,
            "dft_energy": dft_energy,
            "geometry_file": opt_file
        })

        # Update progress
        update_job_status(job_id, ensemble_data={
            "conformers_optimized": i+1
        })

    # STEP 4c: Post-DFT energy filtering
    start_step(job_id, "post_dft_filtering", "Filtering by DFT energies")

    # Rank by DFT energy
    dft_results_sorted = sorted(dft_results, key=lambda x: x['dft_energy'])
    min_energy = dft_results_sorted[0]['dft_energy']

    # Filter by energy window
    filtered = [
        r for r in dft_results_sorted
        if (r['dft_energy'] - min_energy) <= job_status.input.energy_window_post
    ]

    update_job_status(job_id, ensemble_data={
        "conformers_filtered_post": len(filtered)
    })

    # STEP 4d: NMR shielding loop
    start_step(job_id, "nmr_shielding", f"Computing NMR for {len(filtered)} conformers")

    shielding_results = []
    for i, conf_data in enumerate(filtered):
        conf_scratch = job_dir / "scratch" / f"conf_{conf_data['conf_id']:04d}"

        # Run NMR shielding calculation
        result = run_calculation(
            smiles=None,
            job_dir=conf_scratch,
            preset=preset,
            solvent=job_status.input.solvent,
            skip_optimization=True,
            geometry_file=conf_data['geometry_file'],
        )

        shielding_results.append({
            "conf_id": conf_data['conf_id'],
            "dft_energy": conf_data['dft_energy'],
            "shielding_data": result['shielding_data']
        })

        update_job_status(job_id, ensemble_data={
            "conformers_nmr_complete": i+1
        })

    # STEP 4e: Boltzmann averaging
    start_step(job_id, "boltzmann_averaging", "Computing weighted average")

    # Extract energies and calculate weights
    energies = [r['dft_energy'] for r in shielding_results]
    weights = boltzmann_weights(energies, temperature=298.15)

    # Average shieldings per atom
    averaged_shielding = average_shieldings(
        shielding_data_list=[r['shielding_data'] for r in shielding_results],
        weights=weights
    )

    # Convert to shifts
    shifts = shielding_to_shift(
        shielding_data=averaged_shielding,
        functional=preset['functional'].upper(),
        basis_set=preset['nmr_basis_set'],
        solvent=job_status.input.solvent
    )

    # Build NMRResults
    nmr_results = NMRResults(...)

    # Save ensemble metadata
    ensemble_file = job_dir / "output" / "conformer_ensemble.json"
    ensemble_file.write_bytes(orjson.dumps({
        "conformers": shielding_results,
        "energies": energies,
        "weights": weights.tolist(),
        "lowest_energy_conformer": dft_results_sorted[0]['conf_id']
    }))

    # STEP 4f: Visualization
    start_step(job_id, "post_processing", "Generating visualizations")

    # Use lowest-energy conformer for structure visualization
    lowest_conf_file = job_dir / "output" / "optimized" / f"conf_{dft_results_sorted[0]['conf_id']:04d}.xyz"
    optimized_xyz = job_dir / "output" / "optimized.xyz"
    optimized_xyz.write_text(lowest_conf_file.read_text())

    generate_spectrum_plot(...)
    generate_annotated_structure(...)

    # Update final status
    complete_current_step(job_id)
    update_job_status(
        job_id,
        nmr_results=nmr_results,
        optimized_geometry_file=str(optimized_xyz),
        ensemble_data={
            "conformer_energies": energies,
            "conformer_weights": weights.tolist()
        }
    )

    return {"success": True, "job_id": job_id}
```

## Job Directory Structure

### Single-Conformer (Current)

```
data/jobs/{job_id}/
├── status.json
├── output/
│   ├── initial.xyz           # RDKit-generated initial geometry
│   ├── optimized.xyz          # DFT-optimized geometry
│   ├── nmr_results.json       # Shifts and metadata
│   ├── spectrum_1H.svg        # 1H NMR spectrum plot
│   ├── spectrum_13C.svg       # 13C NMR spectrum plot
│   └── structure_annotated.svg  # Molecule with shift labels
├── logs/
└── scratch/
    ├── optimize.nw
    ├── optimize.out
    ├── optimize.log
    ├── shielding.nw
    ├── shielding.out
    └── shielding.log
```

### Ensemble Mode (New)

```
data/jobs/{job_id}/
├── status.json                # Extended with ensemble_data fields
├── output/
│   ├── initial.xyz            # Lowest-energy initial conformer (for quick preview)
│   ├── conformers/            # Generated conformers (pre-optimization)
│   │   ├── conf_0001.xyz
│   │   ├── conf_0002.xyz
│   │   └── ...
│   ├── optimized/             # DFT-optimized conformers
│   │   ├── conf_0001.xyz
│   │   ├── conf_0002.xyz
│   │   └── ...
│   ├── optimized.xyz          # Lowest-energy conformer (for visualization)
│   ├── nmr_results.json       # Boltzmann-averaged shifts
│   ├── conformer_ensemble.json  # Full ensemble data (energies, weights, per-conformer shifts)
│   ├── spectrum_1H.svg        # 1H NMR spectrum (averaged)
│   ├── spectrum_13C.svg       # 13C NMR spectrum (averaged)
│   └── structure_annotated.svg  # Lowest-energy conformer with averaged shift labels
├── logs/
└── scratch/
    ├── conf_0001/
    │   ├── optimize.nw
    │   ├── optimize.out
    │   ├── optimize.log
    │   ├── shielding.nw
    │   ├── shielding.out
    │   └── shielding.log
    ├── conf_0002/
    │   └── ...
    └── ...
```

## API Schema Changes

### JobSubmitRequest (api/schemas.py)

```python
class JobSubmitRequest(BaseModel):
    """Request body for SMILES job submission."""

    smiles: str = Field(...)
    name: Optional[str] = Field(None, max_length=100)
    preset: Literal["draft", "production"] = Field(default="production")
    solvent: str = Field(...)
    notification_email: Optional[EmailStr] = Field(None)

    # NEW FIELDS for conformational sampling
    conformer_mode: Literal["single", "ensemble"] = Field(
        default="single",
        description="Single conformer or ensemble averaging"
    )
    conformer_method: Literal["rdkit", "crest"] = Field(
        default="rdkit",
        description="Method for conformer generation (if ensemble mode)"
    )
    max_conformers: int = Field(
        default=50,
        ge=1,
        le=200,
        description="Maximum conformers to generate"
    )
    energy_window_pre: float = Field(
        default=6.0,
        gt=0,
        description="Pre-DFT energy window in kcal/mol"
    )
    energy_window_post: float = Field(
        default=3.0,
        gt=0,
        description="Post-DFT energy window in kcal/mol"
    )
```

### JobInput (models.py)

```python
class JobInput(BaseModel):
    """Input parameters for a calculation job."""

    model_config = ConfigDict(strict=True)

    smiles: str
    name: Optional[str] = None
    preset: Literal["draft", "production"] = "production"
    solvent: str
    notification_email: Optional[str] = None

    # NEW FIELDS
    conformer_mode: Literal["single", "ensemble"] = "single"
    conformer_method: Literal["rdkit", "crest"] = "rdkit"
    max_conformers: int = 50
    energy_window_pre: float = 6.0
    energy_window_post: float = 3.0
```

### EnsembleData (models.py - NEW)

```python
class ConformerData(BaseModel):
    """Data for a single conformer in an ensemble."""

    conf_id: int
    mmff_energy: Optional[float] = None  # Pre-DFT energy (MMFF or GFN2-xTB)
    dft_energy: float  # DFT energy in Hartree
    boltzmann_weight: float  # Population weight
    geometry_file: Optional[str] = None  # Relative path to XYZ file

class EnsembleData(BaseModel):
    """Ensemble calculation metadata."""

    conformers_generated: int = 0
    conformers_filtered_pre: int = 0
    conformers_optimized: int = 0
    conformers_filtered_post: int = 0
    conformers_nmr_complete: int = 0
    conformer_details: list[ConformerData] = []
```

### JobStatus (models.py)

```python
class JobStatus(BaseModel):
    """Complete status of a calculation job."""

    # ... existing fields ...

    # NEW FIELD
    ensemble_data: Optional[EnsembleData] = None
```

### JobStatusResponse (api/schemas.py)

```python
class JobStatusResponse(BaseModel):
    """Response for job status queries."""

    # ... existing fields ...

    # NEW FIELDS
    conformer_mode: str = Field(default="single")
    ensemble_data: Optional[dict] = Field(
        None,
        description="Ensemble progress (if conformer_mode='ensemble')"
    )
```

## Integration Points

### Modified Files

| File | Changes | Reason |
|------|---------|--------|
| `models.py` | Add `ConformerData`, `EnsembleData`, extend `JobInput` and `JobStatus` | Store ensemble metadata |
| `api/schemas.py` | Add ensemble fields to `JobSubmitRequest` and `JobStatusResponse` | API for ensemble requests |
| `storage.py` | Add `create_conformer_directory()`, `save_conformer_ensemble()` | File management |
| `tasks.py` | Add `run_ensemble_nmr_task()` | Huey task for ensemble jobs |
| `nwchem/runner.py` | Minor: add `conformer_id` param for scratch directory naming | Separate conformer workspaces |
| `api/routers/jobs.py` | Route ensemble jobs to new task | Job submission logic |

### New Files

| File | Purpose | Key Functions |
|------|---------|---------------|
| `conformers.py` | Conformer generation | `generate_rdkit_conformers()`, `generate_crest_conformers()` |
| `averaging.py` | Boltzmann weighting | `boltzmann_weights()`, `average_shieldings()` |
| `ensemble.py` | Ensemble orchestration helpers | `filter_by_energy()`, `parse_conformer_energy()` |

## Answers to Specific Questions

### 1. Should each conformer be a separate Huey task, or one task per ensemble?

**Recommendation: One Huey task per ensemble (entire conformer set processed in single task)**

**Rationale:**
- Single Huey worker on single VM (no parallelization benefit from separate tasks)
- Simpler progress tracking (one job_id → one status.json)
- Easier error handling (entire ensemble succeeds or fails together)
- Boltzmann averaging requires all conformer energies → must wait for all anyway
- Sequential processing is deterministic and easier to debug

**Alternative (rejected):** Separate tasks per conformer
- Would require parent-child task coordination
- More complex status tracking (N tasks → 1 result)
- No performance benefit without multiple workers
- Complicates partial failure handling

### 2. How does job progress work when there are N conformers to process?

**Progress Tracking Strategy:**

```python
# status.json structure
{
    "status": "running",
    "current_step": "dft_optimization",
    "current_step_label": "Optimizing conformer 15/30",
    "ensemble_data": {
        "conformers_generated": 30,
        "conformers_filtered_pre": 30,
        "conformers_optimized": 15,      # Updated after each conformer
        "conformers_filtered_post": 0,
        "conformers_nmr_complete": 0
    }
}
```

**UI displays progress as:**
- "Generating conformers... (30 found)"
- "Optimizing conformers... (15/30 complete, 50%)"
- "Computing NMR shielding... (8/20 complete, 40%)"
- "Averaging results..."

**Implementation:**
- After each conformer DFT optimization: `update_job_status(job_id, ensemble_data={"conformers_optimized": i})`
- After each conformer NMR: `update_job_status(job_id, ensemble_data={"conformers_nmr_complete": i})`
- Use `start_step()` to update `current_step_label` with detailed progress

### 3. Where does conformer generation happen (in Huey worker or before queuing)?

**Recommendation: Inside Huey worker (after queuing)**

**Rationale:**
- Conformer generation can take 30-60 minutes for CREST
- Blocking FastAPI endpoint for this long is unacceptable
- RDKit generation is fast (~seconds) but still should be async for consistency
- Keeps API endpoint responsive (returns job_id immediately)

**Flow:**
```python
# FastAPI endpoint (jobs.py)
@router.post("/jobs")
def submit_job(request: JobSubmitRequest):
    # Quick validation only
    validate_smiles(request.smiles)

    # Create job directory (queued status)
    job_status = create_job_directory(...)

    # Queue task (returns immediately)
    if request.conformer_mode == "ensemble":
        run_ensemble_nmr_task(job_status.job_id)
    else:
        run_nmr_task(job_status.job_id)

    return {"job_id": job_status.job_id}

# Huey worker (tasks.py)
@huey.task()
def run_ensemble_nmr_task(job_id: str):
    # Conformer generation happens here
    start_step(job_id, "conformer_generation", "Generating conformers")
    conformers = generate_rdkit_conformers(...)
    # ... rest of pipeline
```

### 4. How does the job directory structure change for multi-conformer jobs?

See "Job Directory Structure" section above. Key changes:
- Add `output/conformers/` for initial geometries
- Add `output/optimized/` for DFT-optimized geometries
- Add `output/conformer_ensemble.json` for full ensemble metadata
- Change `scratch/` to `scratch/conf_{id:04d}/` subdirectories

### 5. What's the data flow for: generate conformers → optimize each → NMR each → average?

See "Data Flow > Ensemble Job Flow (Detailed)" section above for complete code-level flow.

**Summary:**
1. Generate conformers → save to `output/conformers/conf_XXXX.xyz`
2. Filter by pre-DFT energy → reduce N to M conformers
3. Loop: optimize each → parse DFT energy → save to `output/optimized/conf_XXXX.xyz`
4. Filter by post-DFT energy → reduce M to P conformers
5. Loop: NMR shielding for each → parse shielding data
6. Boltzmann average shieldings → convert to shifts → save results

### 6. How to handle partial failures (some conformers fail, others succeed)?

**Strategy: Graceful Degradation**

```python
# In DFT optimization loop
dft_results = []
failed_conformers = []

for i, conformer in enumerate(conformers):
    try:
        result = run_calculation(...)
        dft_results.append(result)
    except RuntimeError as e:
        # Log failure but continue
        failed_conformers.append({
            "conf_id": i,
            "error": str(e),
            "stage": "optimization"
        })
        continue

# Continue if at least 1 conformer succeeded
if len(dft_results) == 0:
    raise RuntimeError("All conformers failed DFT optimization")

# Save failure info to status
update_job_status(job_id, ensemble_data={
    "failed_conformers": failed_conformers
})
```

**Decision Rules:**
- If ALL conformers fail at any stage → job fails, return error
- If SOME conformers fail:
  - Optimization stage: Continue with successful ones, warn user
  - NMR stage: Continue with successful ones, warn user
  - If <10% of conformers survive → warn user (results may be unreliable)
- Always save `failed_conformers` list to `ensemble_data` for debugging

### 7. How does this affect the existing job metadata JSON schema?

See "API Schema Changes" section above. Key additions:

**JobInput:**
- `conformer_mode` (default: "single")
- `conformer_method` (default: "rdkit")
- `max_conformers` (default: 50)
- `energy_window_pre` (default: 6.0)
- `energy_window_post` (default: 3.0)

**JobStatus:**
- `ensemble_data` (optional, only present if `conformer_mode == "ensemble"`)

**Backwards Compatibility:**
- All new fields have defaults (single-conformer mode is default)
- Old jobs without these fields still load (Pydantic defaults apply)
- Old API clients continue to work (new fields optional)

## Build Order Recommendation

### Phase 1: Core Ensemble Infrastructure (v2.0.1)

1. **Models and storage** (foundation)
   - Add `ConformerData`, `EnsembleData` to `models.py`
   - Extend `JobInput` with ensemble fields
   - Add directory helpers to `storage.py`

2. **RDKit conformer generation** (simplest path)
   - Create `conformers.py` with `generate_rdkit_conformers()`
   - Use KDG (not ETKDG) for solution-phase sampling
   - Implement MMFF optimization and energy filtering

3. **Boltzmann averaging** (core algorithm)
   - Create `averaging.py` with `boltzmann_weights()` and `average_shieldings()`
   - Unit tests with known test cases

4. **Ensemble orchestration** (integration)
   - Add `run_ensemble_nmr_task()` to `tasks.py`
   - Sequential conformer processing loop
   - Progress tracking updates

5. **API integration** (user-facing)
   - Add ensemble fields to `api/schemas.py`
   - Route ensemble jobs in `api/routers/jobs.py`
   - Update status endpoint to show ensemble progress

### Phase 2: CREST Integration (v2.1)

1. **CREST detection and execution**
   - Add `generate_crest_conformers()` to `conformers.py`
   - Check for `crest` and `xtb` binaries at startup
   - Fallback to RDKit if not available

2. **CREST output parsing**
   - Parse multi-XYZ ensemble file
   - Extract GFN2-xTB energies

### Phase 3: Optimization and Resilience (v2.2)

1. **Partial failure handling**
   - Try-catch per conformer with graceful degradation
   - Save `failed_conformers` metadata

2. **Energy window tuning**
   - Make pre/post DFT windows user-configurable
   - Add adaptive filtering (e.g., keep top 95% Boltzmann population)

3. **Performance optimizations**
   - Caching of conformer ensembles
   - Optional: Multi-worker support (future, if needed)

## Key Architectural Decisions

### Decision 1: Single Task vs Multiple Tasks

**Chosen:** Single Huey task per ensemble

**Justification:**
- Single worker environment (no parallelization benefit)
- Simpler orchestration
- Atomic success/failure semantics
- Easier progress tracking

### Decision 2: Sequential vs Parallel Conformer Processing

**Chosen:** Sequential (for now)

**Justification:**
- Single worker with 4 MPI processes already uses all cores
- Running multiple NWChem instances in parallel would compete for resources
- Sequential is simpler and deterministic
- Future: Could parallelize if multiple workers available

### Decision 3: Conformer Generation Location

**Chosen:** Inside Huey worker (async)

**Justification:**
- CREST can take 30-60 minutes
- Cannot block FastAPI endpoint
- Consistent with existing async architecture

### Decision 4: Job Directory Structure

**Chosen:** Separate `conformers/` and `optimized/` subdirectories

**Justification:**
- Clear separation of pre-DFT and post-DFT geometries
- Easy to debug by inspecting intermediate files
- Matches ISiCLE's approach (proven pattern)

### Decision 5: Progress Granularity

**Chosen:** Update after each conformer (fine-grained)

**Justification:**
- User feedback for long-running jobs (hours)
- Shows which conformer is currently processing
- Minimal overhead (disk I/O is fast for JSON)

## Sources

**Existing Codebase (HIGH confidence):**
- `/home/chris/develop/qm-nmr-calc/src/qm_nmr_calc/tasks.py` - Current Huey task structure
- `/home/chris/develop/qm-nmr-calc/src/qm_nmr_calc/storage.py` - Job directory management
- `/home/chris/develop/qm-nmr-calc/src/qm_nmr_calc/models.py` - Data models
- `/home/chris/develop/qm-nmr-calc/src/qm_nmr_calc/api/schemas.py` - API schemas
- `/home/chris/develop/qm-nmr-calc/src/qm_nmr_calc/nwchem/runner.py` - NWChem execution
- `/home/chris/develop/qm-nmr-calc/src/qm_nmr_calc/queue.py` - Huey configuration
- `/home/chris/develop/qm-nmr-calc/src/qm_nmr_calc/shifts.py` - Shielding conversion

**Domain Research (HIGH confidence):**
- `/home/chris/develop/qm-nmr-calc/references/conformational_sampling_nmr_analysis.md` - Conformational sampling workflow and best practices

**Architectural Patterns (MEDIUM confidence):**
- Huey task queue documentation - async job processing patterns
- ISiCLE workflow - conformer ensemble processing in computational chemistry

---

*Architecture research for: Conformational Sampling NMR*
*Researched: 2026-01-26*
*Next: Use this architecture to structure roadmap phases for v2.0 milestone*
