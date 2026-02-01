# Architecture: NMReData Export Integration

**Project:** qm-nmr-calc
**Research Date:** 2026-02-01
**Focus:** NMReData export architecture for existing NMR prediction app
**Confidence:** HIGH

## Executive Summary

NMReData is a standardized format for reporting NMR chemical shift assignments embedded in SDF files using tagged metadata. Integration into the existing FastAPI application requires adding NMReData generation to the post-processing pipeline and exposing a new download endpoint. The format extends standard SDF files with NMREDATA_ prefixed tags containing shift assignments, metadata, and solvent information.

**Recommendation:** Generate NMReData files on-demand (not pre-generated) using RDKit's SDF writing capabilities with custom tag injection. This minimizes storage overhead and ensures the export only happens when users request it.

## NMReData Format Overview

### What is NMReData?

NMReData is the standardized file format for NMR data relevant to structural characterization of small molecules. It extends the Structure Data Format (.sdf) by adding tagged metadata with the prefix "NMREDATA_" to include signal assignments, chemical shifts, couplings, 2D correlations, and links to spectra.

**Key characteristics:**
- Based on SDF/MOL file format (RDKit-compatible)
- Chemical structure with 3D coordinates in MOL block
- Tagged metadata fields for NMR data
- Chemical shifts as scalar values (not ranges) in ppm
- Atom numbering matches SDF structure numbering

**Official sources:**
- NMReData Initiative: https://github.com/NMReDATAInitiative
- Specification: https://nmredata.org/wiki/NMReDATA_tag_format
- Academic paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC6226248/

### File Structure

```
[MOL block - chemical structure with 3D coordinates]
[atom count, bond count, etc.]
[atom lines: element, x, y, z, ...]
[bond lines: atom1, atom2, bondtype, ...]
M  END

> <NMREDATA_VERSION>
1.1

> <NMREDATA_LEVEL>
0

> <NMREDATA_SOLVENT>
CDCl3

> <NMREDATA_TEMPERATURE>
298.0

> <NMREDATA_ASSIGNMENT>
Ha, 7.2453, 5
Hb, 3.4521, 8,9,10
C1, 128.45, 1

$$$$
```

### Required Tags

| Tag | Description | Example | Notes |
|-----|-------------|---------|-------|
| `NMREDATA_VERSION` | Format version | `1.1` | Mandatory |
| `NMREDATA_LEVEL` | Complexity indicator | `0` | 0=basic assignment, 3=full coupling |
| `NMREDATA_SOLVENT` | NMR solvent | `CDCl3` | Mandatory |
| `NMREDATA_ASSIGNMENT` | Shift assignments | `Ha, 7.2453, 5` | Core data |
| `NMREDATA_TEMPERATURE` | Temperature in K | `298.0` | Strongly recommended |

### Optional Tags (Relevant to QM-NMR-Calc)

| Tag | Description | Our Use |
|-----|-------------|---------|
| `NMREDATA_FORMULA` | Molecular formula | Can derive from SMILES |
| `NMREDATA_SMILES` | SMILES string | Already stored in job input |
| `NMREDATA_ID` | DOI/identifier | Could use job_id |
| `NMREDATA_INCHI` | InChI string | Can generate from RDKit |

### NMREDATA_ASSIGNMENT Tag Syntax

**Format:** `Label, ChemicalShift, AtomIndices`

**Rules:**
1. Chemical shift in ppm with 4 decimal places: `7.2453` not `7.2`
2. Only single scalar values (no ranges)
3. Atom indices are 1-based (matches SDF numbering convention)
4. Multiple atoms for equivalent positions: `Ha, 1.2345, 4,5,6`
5. Separator: comma-space (`, `)
6. Unknown shifts: use `777.777`

**Examples:**
```
H1, 7.2453, 5
H2, 3.4521, 8
Methyl, 1.2340, 10,11,12
C1, 128.4567, 1
C2, 45.6789, 2
```

**Atom numbering correspondence:**
- SDF atom numbering: 1-based (line order in MOL block)
- RDKit atom indexing: 0-based (internal representation)
- NWChem atom indices: 1-based (our job results use this)
- **NMReData requirement:** 1-based (matches SDF)

**Conversion for our system:**
- Our `AtomShift.index` is already 1-based (NWChem convention)
- RDKit SDF will write atoms in order (0-based becomes 1-based in file)
- **CRITICAL:** Verify atom ordering matches between SMILES->RDKit->SDF and NWChem input geometry

## Current Architecture Analysis

### Existing Job Result Flow

```
[NMR Calculation Task] → run_nmr_task() or run_ensemble_nmr_task()
    ↓
1. NWChem calculation completes
    ↓
2. Parse shieldings from NWChem output
    ↓
3. Convert shieldings to shifts using DELTA50 scaling
    ↓
4. Build NMRResults model (AtomShift objects)
    ↓
5. Save nmr_results.json
    ↓
6. Generate PNG/SVG spectrum plots (1H, 13C)
    ↓
7. Generate PNG/SVG annotated structure
    ↓
8. Save optimized.xyz geometry
    ↓
9. Update job status with nmr_results
    ↓
[Job Complete - status="complete"]
```

### Current File Storage

**Job directory structure:**
```
data/jobs/{job_id}/
├── status.json           # JobStatus model with NMRResults
├── output/
│   ├── nmr_results.json  # AtomShift data
│   ├── optimized.xyz     # DFT-optimized geometry
│   ├── initial.xyz       # RDKit initial geometry
│   ├── spectrum_1H.svg
│   ├── spectrum_1H.png
│   ├── spectrum_13C.svg
│   ├── spectrum_13C.png
│   ├── structure_annotated.svg
│   └── structure_annotated.png
├── logs/
└── scratch/
```

### Current Download Endpoints

**Existing FastAPI routes:**
```
GET /api/v1/jobs/{job_id}/results          → JSON (NMRResults)
GET /api/v1/jobs/{job_id}/geometry         → XYZ file
GET /api/v1/jobs/{job_id}/geometry.sdf     → SDF file (on-demand generated)
GET /api/v1/jobs/{job_id}/geometry.json    → JSON (for 3D viewer)
GET /api/v1/jobs/{job_id}/spectrum/1h.png  → PNG
GET /api/v1/jobs/{job_id}/spectrum/1h.svg  → SVG
GET /api/v1/jobs/{job_id}/spectrum/13c.png → PNG
GET /api/v1/jobs/{job_id}/spectrum/13c.svg → SVG
GET /api/v1/jobs/{job_id}/structure.png    → PNG
GET /api/v1/jobs/{job_id}/structure.svg    → SVG
GET /api/v1/jobs/{job_id}/output           → ZIP (raw NWChem files)
```

### Existing SDF Generation Pattern

From `/api/v1/jobs/{job_id}/geometry.sdf`:

```python
# Read XYZ coordinates
xyz_content = geometry_file.read_text()
xyz_lines = xyz_content.strip().split("\n")

# Parse XYZ: first line is atom count, second is comment, rest are coords
coords = []
for line in xyz_lines[2:]:
    parts = line.split()
    if len(parts) >= 4:
        coords.append((float(parts[1]), float(parts[2]), float(parts[3])))

# Create mol from SMILES and set coordinates
mol = Chem.MolFromSmiles(job_status.input.smiles)
mol = Chem.AddHs(mol)

# Create conformer with coordinates
conf = Chem.Conformer(mol.GetNumAtoms())
for i, (x, y, z) in enumerate(coords):
    if i < mol.GetNumAtoms():
        conf.SetAtomPosition(i, (x, y, z))
mol.AddConformer(conf, assignId=True)

# Generate SDF content
sdf_content = Chem.MolToMolBlock(mol)
```

**This pattern already exists and works!** We can extend it to add NMReData tags.

## Integration Architecture

### Design Decision: On-Demand vs Pre-Generated

**Option A: Pre-generate during post-processing**
- Generate NMReData file in `tasks.py` alongside PNG/SVG outputs
- Store as `output/nmredata.sdf` (single mode) or `output/nmredata_ensemble.sdf` (ensemble mode)
- Endpoint serves static file

**Option B: Generate on-demand (RECOMMENDED)**
- No file storage
- Generate NMReData SDF when endpoint is called
- Similar to current `/geometry.sdf` endpoint

**Recommendation: Option B (On-Demand)**

**Rationale:**
1. **Storage efficiency:** NMReData export is a niche use case. Most users won't download it.
2. **Consistency:** Existing `/geometry.sdf` already uses on-demand generation pattern
3. **Flexibility:** Easy to update format or add fields without regenerating old jobs
4. **Maintainability:** Single source of truth (nmr_results.json), NMReData is a view
5. **Performance:** SDF generation is fast (<100ms), not a bottleneck

**Trade-off:**
- Slight latency on download (acceptable for file export)
- Repeated generation if downloaded multiple times (rare in practice)

### Proposed Endpoint

**Route:** `GET /api/v1/jobs/{job_id}/nmredata.sdf`

**Response:**
- Media type: `chemical/x-mdl-sdfile`
- Filename: `{job_id}_nmredata.sdf`
- Status codes:
  - 200: Success
  - 404: Job not found
  - 409: Job not complete

**Naming convention:**
- Single conformer: `{job_id}_nmredata.sdf`
- Ensemble (averaged): `{job_id}_nmredata_ensemble.sdf`
- Individual conformers: `{job_id}_nmredata_conf_001.sdf` (optional future feature)

### Data Flow for NMReData Generation

```
[User Request: GET /jobs/{id}/nmredata.sdf]
    ↓
1. Load job status
    ↓
2. Validate job complete
    ↓
3. Load nmr_results from job status
    ↓
4. Load optimized.xyz geometry
    ↓
5. Create RDKit Mol from SMILES
    ↓
6. Add explicit hydrogens
    ↓
7. Set 3D coordinates from XYZ
    ↓
8. Generate base SDF (MolToMolBlock)
    ↓
9. Inject NMReData tags
    - NMREDATA_VERSION
    - NMREDATA_LEVEL
    - NMREDATA_SOLVENT
    - NMREDATA_TEMPERATURE
    - NMREDATA_ASSIGNMENT (from h1_shifts + c13_shifts)
    - NMREDATA_FORMULA
    - NMREDATA_SMILES
    ↓
10. Return as Response with proper headers
```

### Component Integration Points

#### 1. API Router: `src/qm_nmr_calc/api/routers/jobs.py`

**New endpoint:**
```python
@router.get(
    "/{job_id}/nmredata.sdf",
    response_class=Response,
    responses={
        200: {
            "description": "NMReData file with chemical shift assignments",
            "content": {"chemical/x-mdl-sdfile": {}}
        },
        404: {"model": ProblemDetail, "description": "Job not found"},
        409: {"model": ProblemDetail, "description": "Job not complete"},
    },
)
async def download_nmredata(job_id: str):
    """Download NMR results in NMReData format (SDF with assignments)."""
    # Implementation in new module
```

#### 2. New Module: `src/qm_nmr_calc/nmredata.py`

**Purpose:** NMReData file generation logic

**Key function:**
```python
def generate_nmredata_sdf(
    smiles: str,
    geometry_xyz: str,
    h1_shifts: list[AtomShift],
    c13_shifts: list[AtomShift],
    solvent: str,
    temperature_k: float = 298.0,
    metadata: dict | None = None,
) -> str:
    """Generate NMReData-formatted SDF file.

    Args:
        smiles: SMILES string (from job input)
        geometry_xyz: XYZ coordinate file content
        h1_shifts: 1H chemical shifts with atom indices
        c13_shifts: 13C chemical shifts with atom indices
        solvent: NMR solvent name
        temperature_k: Temperature in Kelvin
        metadata: Optional dict with formula, version info, etc.

    Returns:
        NMReData SDF file content as string
    """
```

**Implementation strategy:**
1. Reuse existing `_xyz_to_sdf()` pattern from `jobs.py`
2. Generate base MOL block with RDKit
3. Parse and append NMReData tags
4. Use RDKit's `SetProp()` to add SDF tags (cleaner than string manipulation)

#### 3. RDKit Property Injection Pattern

**RDKit supports SDF tag injection via `SetProp()`:**

```python
from rdkit import Chem

mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)

# Set 3D coordinates
# ... (existing pattern)

# Inject NMReData tags as molecule properties
mol.SetProp("NMREDATA_VERSION", "1.1")
mol.SetProp("NMREDATA_LEVEL", "0")
mol.SetProp("NMREDATA_SOLVENT", solvent)
mol.SetProp("NMREDATA_TEMPERATURE", str(temperature_k))

# Build assignment tag content
assignment_lines = []
for shift in h1_shifts:
    # Format: Label, ChemicalShift, AtomIndex
    label = f"H{shift.index}"
    assignment_lines.append(f"{label}, {shift.shift:.4f}, {shift.index}")
for shift in c13_shifts:
    label = f"C{shift.index}"
    assignment_lines.append(f"{label}, {shift.shift:.4f}, {shift.index}")

assignment_content = "\n".join(assignment_lines)
mol.SetProp("NMREDATA_ASSIGNMENT", assignment_content)

# Generate SDF
sdf_content = Chem.MolToMolBlock(mol)
```

**RDKit automatically formats properties as SDF tags!**

#### 4. Atom Ordering Verification

**CRITICAL CONCERN:** Atom ordering must match between:
- NWChem calculation (produces atom indices)
- RDKit Mol object (internal atom order)
- SDF file (atom line order)

**Current approach verification needed:**

1. **Job submission:** User provides SMILES
2. **Initial geometry:** `_generate_initial_xyz()` uses RDKit embedding
3. **NWChem input:** Geometry generated from initial XYZ (order preserved)
4. **NWChem output:** Atom indices are 1-based, order matches input
5. **Result storage:** `AtomShift.index` stores NWChem 1-based indices
6. **SDF generation:** `MolFromSmiles()` creates Mol, `AddHs()` adds hydrogens

**Potential issue:** Does `MolFromSmiles()` → `AddHs()` produce the same atom order as the initial XYZ generation?

**Mitigation strategies:**

**Strategy A: Trust canonical SMILES ordering (RISKY)**
- Assume RDKit produces consistent atom ordering from SMILES
- **Risk:** Hydrogen placement may vary between RDKit versions or parameters

**Strategy B: Verify atom mapping (RECOMMENDED for Phase 1)**
- Add validation that RDKit Mol atom count matches NWChem
- Document assumption that ordering is consistent
- Flag for future enhancement if issues arise

**Strategy C: Store atom mapping explicitly (FUTURE)**
- During initial geometry generation, store RDKit → NWChem atom mapping
- Use mapping to reorder atoms in NMReData generation
- Requires changing JobStatus model (out of scope for this milestone)

**Recommendation:** Use Strategy B with documentation. The existing `/geometry.sdf` endpoint already makes this assumption, so NMReData inherits the same risk profile. Flag for validation testing.

### Ensemble Mode Considerations

**Question:** How to handle ensemble calculations with Boltzmann-weighted shifts?

**Options:**

**Option 1: Single averaged NMReData file (RECOMMENDED)**
- Use Boltzmann-averaged shifts (already computed)
- Geometry: lowest-energy conformer (representative structure)
- Tags indicate ensemble averaging in metadata
- Simple, matches user expectation

**Option 2: Multiple NMReData files**
- One file per conformer with individual shifts
- Endpoint: `/jobs/{id}/nmredata/conformers/{conf_id}.sdf`
- Complex, requires ensemble-aware NMReData consumers

**Option 3: Multi-conformer SDF (advanced)**
- Single SDF with multiple conformers
- NMReData spec supports this but tooling may not
- High complexity

**Recommendation:** Option 1 for MVP. Use Boltzmann-averaged shifts with lowest-energy conformer geometry. Add optional metadata tag to indicate ensemble calculation.

**Custom tag for ensemble metadata:**
```
> <NMREDATA_COMMENT>
Shifts are Boltzmann-weighted averages from ensemble of 15 conformers at 298.15 K. Geometry represents lowest-energy conformer (DFT-optimized, B3LYP/6-31G*).
```

## Implementation Plan

### Phase 1: Core NMReData Generation (MVP)

**Files to create/modify:**

1. **Create:** `src/qm_nmr_calc/nmredata.py`
   - `generate_nmredata_sdf()` function
   - Tag formatting helpers
   - Single-conformer support only

2. **Modify:** `src/qm_nmr_calc/api/routers/jobs.py`
   - Add `GET /{job_id}/nmredata.sdf` endpoint
   - Use existing validation patterns (job complete check)
   - Call `generate_nmredata_sdf()` from new module

3. **Modify:** `docs/API.md` (if exists)
   - Document new endpoint
   - Provide NMReData format reference links

**Testing requirements:**

1. **Unit tests:** `tests/test_nmredata.py`
   - Tag formatting
   - Atom index handling
   - Multi-atom assignment (equivalent hydrogens)

2. **Integration tests:** `tests/api/test_nmredata_endpoint.py`
   - Complete job returns valid SDF
   - Incomplete job returns 409
   - Non-existent job returns 404
   - Validate SDF parseable by RDKit

3. **Validation tests:**
   - Parse generated NMReData with external tool (if available)
   - Verify atom numbering matches between shifts and structure

**Success criteria:**
- Endpoint returns valid SDF file
- NMReData tags present and correctly formatted
- Chemical shifts match nmr_results.json
- File parseable by RDKit and other SDF tools

### Phase 2: Ensemble Support (Post-MVP)

**Enhancements:**

1. **Ensemble metadata tag:**
   - Add `NMREDATA_COMMENT` with ensemble info
   - Include conformer count, temperature, energy window

2. **Representative geometry selection:**
   - Use lowest-energy conformer geometry
   - Document this choice in comment tag

3. **Optional individual conformer export:**
   - New endpoint: `/jobs/{id}/nmredata/conformers/{conf_id}.sdf`
   - Individual shifts (not averaged)
   - Requires parsing per-conformer NMR results (currently Boltzmann-averaged)

**Blocker:** Current ensemble implementation only stores averaged shifts, not per-conformer shifts. Individual conformer export requires architectural change to preserve per-conformer NMR results.

### Phase 3: Advanced Features (Future)

**Potential enhancements:**

1. **Coupling constants (J-values):**
   - Add `NMREDATA_J` tag
   - Requires computing J from NWChem output (not currently parsed)
   - Level 1-3 NMReData compliance

2. **Multiplicity information:**
   - Add `NMREDATA_1D_1H` spectrum tags with multiplicity
   - Requires analyzing J-coupling patterns

3. **2D correlation data:**
   - Add `NMREDATA_2D_*` tags
   - Requires HSQC/HMBC data (not computed by QM-NMR-Calc)

4. **Atom mapping validation:**
   - Store RDKit → NWChem atom mapping during initial geometry generation
   - Use mapping to ensure correct atom indexing in NMReData

5. **Format version support:**
   - Support NMReData version 2.0+ features
   - Metadata tags for authors, institutions

## Technical Specifications

### Dependencies

**No new dependencies required:**
- RDKit: Already used for SDF generation
- Standard library: String formatting, file I/O

**Optional (for validation testing):**
- External NMReData parser (Java tools from NMReDATAInitiative)
- Not required for runtime, only for test validation

### Performance Considerations

**SDF generation cost:**
- `MolFromSmiles()`: ~1ms
- `AddHs()`: ~1ms
- Coordinate assignment: <1ms
- `MolToMolBlock()`: ~1ms
- Tag formatting: <1ms

**Total:** <5ms per request (negligible)

**Comparison:** Existing `/geometry.sdf` has same cost profile and performs well in production.

### Error Handling

**Failure modes:**

1. **Invalid SMILES (should never happen):**
   - Job submission already validates SMILES
   - If RDKit fails: Return 500 with error detail

2. **Atom count mismatch:**
   - RDKit Mol atom count != shift count
   - Log warning, include all available shifts
   - Return SDF with partial data + comment

3. **Missing geometry file:**
   - Return 404 "Geometry file not found"
   - Should not happen for complete jobs

4. **Tag formatting error:**
   - Catch exception, log details
   - Return 500 with sanitized error message

### Security Considerations

**No new attack surface:**
- Read-only operation (no file writes)
- Input validation already done during job submission
- SMILES and geometry data already trusted (came from our system)

**Potential concerns:**
- Large molecules: SDF generation scales with atom count (acceptable, existing `/geometry.sdf` has same profile)
- Malformed XYZ: Caught by existing XYZ parser in `/geometry.sdf`

## Alternative Architectures Considered

### Alternative 1: Pre-Generate During Post-Processing

**Approach:** Add NMReData generation to `tasks.py` post-processing step

**Pros:**
- Faster downloads (static file serving)
- All exports generated at once (consistency)

**Cons:**
- Storage overhead for rarely-used format
- Inflexible (can't update format without reprocessing)
- Breaks existing pattern (only visualizations are pre-generated)

**Verdict:** Rejected. On-demand generation is more maintainable and follows existing patterns.

### Alternative 2: Separate NMReData Export Service

**Approach:** Microservice for NMReData generation

**Pros:**
- Decoupled from main API
- Could support multiple input formats

**Cons:**
- Over-engineered for simple file generation
- Adds deployment complexity
- Network latency between services

**Verdict:** Rejected. Not justified for this use case. Single endpoint is sufficient.

### Alternative 3: Client-Side Generation

**Approach:** Return nmr_results.json and let client generate NMReData

**Pros:**
- No server-side implementation needed
- Maximum flexibility for clients

**Cons:**
- Pushes complexity to every client
- No canonical NMReData export
- Requires JavaScript NMReData library (doesn't exist)

**Verdict:** Rejected. Server-side generation provides better UX.

## Open Questions and Future Research

### 1. Atom Ordering Validation

**Question:** Does `MolFromSmiles()` → `AddHs()` produce consistent atom ordering with NWChem geometry?

**Research needed:**
- Empirical testing with multiple molecules
- Compare RDKit atom indices with NWChem indices
- Document any discrepancies

**Mitigation:** Phase 1 assumes consistency (same as existing `/geometry.sdf`). Flag for validation in Phase 2.

### 2. NMReData Consumer Software

**Question:** What software can read/validate our NMReData files?

**Research needed:**
- Test with Mestrelab Mnova (commercial NMR software)
- Test with ChemInfo NMReData viewer (https://www.cheminfo.org/flavor/nmredata/viewer/)
- Test with NMReDATAInitiative Java tools

**Action:** Phase 1 validation testing should include external parser.

### 3. Equivalent Atoms (Symmetry)

**Question:** How to handle chemically equivalent atoms (e.g., methyl CH3)?

**Current approach:** Assign same shift to all three hydrogens with individual labels:
```
H8, 1.2340, 8
H9, 1.2340, 9
H10, 1.2340, 10
```

**Alternative:** Single label for all equivalent atoms:
```
Methyl_H, 1.2340, 8,9,10
```

**Research needed:**
- What's the NMReData best practice?
- Does it affect downstream software compatibility?

**Decision for Phase 1:** Use individual labels (conservative approach). Revisit in Phase 2 based on validation feedback.

### 4. Solvent Name Mapping

**Question:** Do our NWChem COSMO solvent names match NMReData conventions?

**Our solvents:** `water`, `dmso`, `chloroform`, etc. (from solvents.py)
**NMReData convention:** Typically uses deuterated solvents: `D2O`, `DMSO-d6`, `CDCl3`

**Research needed:**
- Map our solvent names to NMReData conventions
- Add mapping table in `nmredata.py`

**Example mapping:**
```python
SOLVENT_MAPPING = {
    "water": "D2O",
    "dmso": "DMSO-d6",
    "chloroform": "CDCl3",
    "acetonitrile": "CD3CN",
    "methanol": "CD3OD",
}
```

**Action:** Phase 1 implementation should include solvent name mapping.

## Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| Atom ordering mismatch | Low | High | Validation testing, document assumptions |
| Solvent name incompatibility | Medium | Low | Add mapping table |
| NMReData spec changes | Low | Medium | Pin to version 1.1, monitor spec updates |
| Performance degradation | Very Low | Low | SDF generation is fast, same as existing endpoint |
| RDKit version incompatibility | Low | Medium | Pin RDKit version in pyproject.toml |

**Overall risk: LOW**

Most risks are mitigated by following existing patterns and adding validation testing.

## Success Metrics

**Phase 1 (MVP):**
- Endpoint returns valid NMReData SDF for single-conformer jobs
- All required tags present and correctly formatted
- Chemical shifts match nmr_results.json (within floating-point precision)
- File parseable by RDKit `SDMolSupplier`
- 100% test coverage for nmredata.py module
- API documentation updated

**Phase 2 (Ensemble):**
- Ensemble jobs return Boltzmann-averaged NMReData
- Ensemble metadata tag includes conformer count and temperature
- Geometry represents lowest-energy conformer

**Phase 3 (Validation):**
- NMReData validated with external parser (Mnova or ChemInfo)
- Atom ordering verified with test suite
- No user-reported issues with atom indexing

## References

**NMReData Specification:**
- [Official website](https://nmredata.org/wiki/Main_Page)
- [Tag format](https://nmredata.org/wiki/NMReDATA_tag_format)
- [GitHub](https://github.com/NMReDATAInitiative)

**Academic Sources:**
- [NMReData paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC6226248/)
- [PubMed](https://pubmed.ncbi.nlm.nih.gov/29656574/)

**Software:**
- [NMReData Java tools](https://github.com/NMReDATAInitiative/javatools)
- [ChemInfo viewer](https://www.cheminfo.org/flavor/nmredata/viewer/)

**RDKit Documentation:**
- [Getting Started](https://www.rdkit.org/docs/GettingStartedInPython.html)
- [SDWriter](https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html)

**Existing Codebase:**
- Current SDF generation: `src/qm_nmr_calc/api/routers/jobs.py:664-745`
- Job models: `src/qm_nmr_calc/models.py`
- Storage utilities: `src/qm_nmr_calc/storage.py`
