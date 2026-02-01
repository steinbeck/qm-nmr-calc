# NMReData Export Stack Research

**Project:** qm-nmr-calc v2.3 NMReData Export
**Researched:** 2026-02-01
**Confidence:** HIGH

## Executive Summary

**NMReData is an SDF extension, not a separate format.** Add custom tags with `NMREDATA_` prefix to existing SDF files. No dedicated Python library exists. Use RDKit's existing `SDWriter` + `SetProp()` to write NMReData-compliant files directly.

**Key finding:** The project already has all required dependencies (RDKit). Implementation is straightforward tag formatting, not library integration.

## Format Specification

### What is NMReData?

NMReData (NMR electronic DATA initiative) is a **standard to report NMR assignment and parameters** for organic compounds. Published in Magnetic Resonance in Chemistry (2018).

**Technical structure:**
- Extension of standard SDF (Structure Data Format)
- SDF = MOL block + custom tags
- NMReData = SDF + tags prefixed with `NMREDATA_`
- Version 1.1 is current stable, version 2.0 adds extensions

**Official source:** https://nmredata.org (NMReDATAInitiative GitHub organization)

### Core Tags for Chemical Shift Predictions

Based on official examples from [NMReDATAInitiative/Examples-of-NMR-records](https://github.com/NMReDATAInitiative/Examples-of-NMR-records):

#### Required Tags

| Tag | Purpose | Format |
|-----|---------|--------|
| `NMREDATA_VERSION` | Format version | `1.1\` |
| `NMREDATA_ASSIGNMENT` | Links signals to atoms with shifts | Multi-line, see below |

#### Highly Recommended Tags

| Tag | Purpose | Example |
|-----|---------|---------|
| `NMREDATA_SOLVENT` | Solvent used | `CDCl3\` or `(CD3)2CO\` |
| `NMREDATA_LEVEL` | Assignment confidence | `0` (predicted) or `1` (experimental) |
| `NMREDATA_1D_1H` | 1H spectrum metadata | Larmor frequency + signal list |
| `NMREDATA_1D_13C` | 13C spectrum metadata | Larmor frequency + signal list |

#### Optional Tags

| Tag | Purpose |
|-----|---------|
| `NMREDATA_SMILES` | Canonical SMILES |
| `NMREDATA_INCHI` | InChI identifier |
| `NMREDATA_ID` | Database reference |

### NMREDATA_ASSIGNMENT Format

**Structure:** `label, chemical_shift, atom_number[, atom_number, ...]`

**Rules:**
- Chemical shift must be single float (no ranges allowed)
- Use `777.777` if chemical shift unknown
- Atom numbers are 1-indexed (matching MOL file numbering)
- Multiple atoms can be assigned to same signal (equivalent atoms)
- Lines end with backslash `\`

**Example from official repository:**
```
> <NMREDATA_ASSIGNMENT>
s0, 0.89, 20, 21, 22\
s1, 1.37, 34\
s2, 1.42, 32\
s8, 3.84, 23\
s9, 4.18, 18\
Interchangeable=s14, s11\
Interchangeable=s2, s4\
```

**Interchangeable groups:** Optional lines indicating assignment ambiguity (e.g., ortho/meta protons).

### NMREDATA_1D_1H and NMREDATA_1D_13C Format

**Structure:**
```
Spectrum_Location=path/to/spectrum\
Larmor=frequency_in_MHz\
chemical_shift, L=signal_label\
chemical_shift, L=signal_label\
...
```

**Example:**
```
> <NMREDATA_1D_1H>
Spectrum_Location=molecule/10027826\
Larmor=400.0\
0.89, L=s0\
1.37, L=s1\
3.84, L=s8\
```

### Complete File Example

From https://raw.githubusercontent.com/NMReDATAInitiative/Examples-of-NMR-records/master/ambiguous_level_1/10027836.nmredata.sdf:

```
[MOL block with atom coordinates and bonds]
M  END
> <NMREDATA_1D_1H>
Spectrum_Location=molecule/10027826\
Larmor=400.0\
0.89, L=s0\
1.37, L=s1\

> <NMREDATA_VERSION>
1.1\

> <NMREDATA_SOLVENT>
(CD3)2CO\

> <NMREDATA_ASSIGNMENT>
s0, 0.89, 20, 21, 22\
s1, 1.37, 34\
s2, 1.42, 32\

> <NMREDATA_LEVEL>
1\

$$$$
```

**Key observations:**
1. Standard SDF structure (MOL block, `M  END`, then tags)
2. Each tag starts with `> <TAGNAME>`
3. Multi-line tag values end each line with `\`
4. File ends with `$$$$` (SDF record separator)

## Python Library Landscape

### Dedicated NMReData Libraries: None

**Finding:** No Python libraries specifically for NMReData generation exist as of 2026-02-01.

**Evidence:**
- GitHub search: `org:NMReDATAInitiative language:python` returns 0 repositories
- NMReDATAInitiative official repos: Java tools, MATLAB tools, JavaScript visualizers
- djeanner/NMReDATA repo: 88.4% MATLAB, 3.3% C, 3.0% Java, 0% Python
- PyPI search (via web): No packages named `nmredata`, `pynmredata`, `nmredata-py`
- Related NMR libraries (nmrglue, nmrstarlib) handle spectral data, not NMReData format

**Why no library needed:** NMReData is tag formatting convention, not complex binary format. RDKit already provides SDF writing with custom tags.

### Existing Stack: RDKit (Already Installed)

**What qm-nmr-calc already has:**
- RDKit for MOL/SDF parsing and writing
- FastAPI for REST endpoints
- File system storage for results

**What RDKit provides for NMReData:**
```python
from rdkit import Chem

# Read or create molecule
mol = Chem.MolFromSmiles("CCO")

# Set properties (tags) on molecule
mol.SetProp("NMREDATA_VERSION", "1.1")
mol.SetProp("NMREDATA_SOLVENT", "CDCl3")
mol.SetProp("NMREDATA_ASSIGNMENT", "s0, 3.65, 1, 2\ns1, 1.18, 3, 4, 5")

# Write to SDF file
writer = Chem.SDWriter("output.nmredata.sdf")
writer.write(mol)
writer.close()
```

**Confidence:** HIGH — Verified from [RDKit documentation](https://www.rdkit.org/docs/GettingStartedInPython.html) and [RDKit blog on properties](https://greglandrum.github.io/rdkit-blog/posts/2025-01-24-property-tutorial.html).

### Alternative: Manual SDF Writing

Not recommended. RDKit's `SDWriter` handles:
- MOL block formatting (atom coordinates, bonds)
- V2000/V3000 format details
- Property tag escaping
- Record separators

Writing SDF manually is error-prone and duplicates RDKit functionality.

## Implementation Approach

### Recommended Strategy: Direct RDKit Integration

**Phase 1: Core NMReData generation**
1. Add `nmredata.py` module to `src/qm_nmr_calc/`
2. Implement `generate_nmredata_sdf(mol, predictions, metadata)` function
3. Format NMREDATA_ASSIGNMENT from atom-shift mappings
4. Format NMREDATA_1D_1H and NMREDATA_1D_13C from predicted shifts
5. Set tags via `mol.SetProp()`
6. Write via `Chem.SDWriter()`

**Phase 2: API integration**
1. Add `/api/jobs/{job_id}/nmredata` endpoint
2. Load completed job results
3. Call `generate_nmredata_sdf()`
4. Return file download response

**Phase 3: Web UI integration**
1. Add download button on results page
2. Link to new API endpoint
3. Style consistently with existing PNG/XYZ downloads

### Data Mapping (qm-nmr-calc → NMReData)

**Available in qm-nmr-calc results:**
- Optimized 3D geometry (from NWChem) → MOL block coordinates
- Predicted 1H shifts with atom indices → NMREDATA_1D_1H + NMREDATA_ASSIGNMENT
- Predicted 13C shifts with atom indices → NMREDATA_1D_13C + NMREDATA_ASSIGNMENT
- Solvent (CHCl3/DMSO/vacuum) → NMREDATA_SOLVENT
- Molecule object (RDKit Mol) → already available

**Mapping table:**

| qm-nmr-calc Data | NMReData Tag | Transformation |
|------------------|--------------|----------------|
| `results['1H_shifts']` | `NMREDATA_1D_1H` | Format as `shift, L=s{idx}\` per line |
| `results['13C_shifts']` | `NMREDATA_1D_13C` | Format as `shift, L=s{idx}\` per line |
| `results['1H_shifts']` | `NMREDATA_ASSIGNMENT` | Format as `s{idx}, shift, atom_num\` |
| `results['13C_shifts']` | `NMREDATA_ASSIGNMENT` | Append to 1H assignments |
| `metadata['solvent']` | `NMREDATA_SOLVENT` | Map: CHCl3→CDCl3, DMSO→(CD3)2SO, vacuum→None |
| `optimized_geometry` | MOL block | Use `Chem.MolToMolBlock(mol)` |
| Constant: "1.1" | `NMREDATA_VERSION` | Hardcoded |
| Constant: "0" | `NMREDATA_LEVEL` | Predicted (0), not experimental (1) |

### Signal Labeling Convention

**For qm-nmr-calc predictions:**
- Label format: `h{atom_num}` for 1H, `c{atom_num}` for 13C
- Example: Atom 1 (carbon) → `c1`, attached hydrogens → `h2`, `h3`, `h4`
- Avoids generic `s0`, `s1` from experimental spectra
- Clearer link between prediction and structure

**Rationale:** Predicted shifts have unambiguous atom assignments (from QM calculation), unlike experimental spectra where assignment is the challenge.

### Code Structure

```python
# src/qm_nmr_calc/nmredata.py

def generate_nmredata_sdf(
    mol: Chem.Mol,
    h1_shifts: dict[int, float],
    c13_shifts: dict[int, float],
    solvent: str,
    larmor_h1: float = 400.0,
    larmor_c13: float = 100.0,
) -> str:
    """
    Generate NMReData-compliant SDF file content.

    Args:
        mol: RDKit molecule with 3D coordinates
        h1_shifts: {atom_idx: chemical_shift} for 1H
        c13_shifts: {atom_idx: chemical_shift} for 13C
        solvent: "CHCl3", "DMSO", or "vacuum"
        larmor_h1: 1H Larmor frequency in MHz
        larmor_c13: 13C Larmor frequency in MHz

    Returns:
        Complete NMReData SDF file as string
    """
    # Set version and level
    mol.SetProp("NMREDATA_VERSION", "1.1")
    mol.SetProp("NMREDATA_LEVEL", "0")  # Predicted

    # Set solvent
    solvent_map = {
        "CHCl3": "CDCl3",
        "DMSO": "(CD3)2SO",
        "vacuum": "vacuum"
    }
    if solvent in solvent_map:
        mol.SetProp("NMREDATA_SOLVENT", solvent_map[solvent])

    # Format 1H spectrum
    h1_lines = [f"Larmor={larmor_h1}"]
    for atom_idx, shift in sorted(h1_shifts.items()):
        h1_lines.append(f"{shift:.2f}, L=h{atom_idx}")
    mol.SetProp("NMREDATA_1D_1H", "\\\n".join(h1_lines) + "\\")

    # Format 13C spectrum
    c13_lines = [f"Larmor={larmor_c13}"]
    for atom_idx, shift in sorted(c13_shifts.items()):
        c13_lines.append(f"{shift:.2f}, L=c{atom_idx}")
    mol.SetProp("NMREDATA_1D_13C", "\\\n".join(c13_lines) + "\\")

    # Format assignments
    assignment_lines = []
    for atom_idx, shift in sorted(h1_shifts.items()):
        # atom_idx is 0-indexed in RDKit, need 1-indexed for MOL
        assignment_lines.append(f"h{atom_idx}, {shift:.2f}, {atom_idx + 1}")
    for atom_idx, shift in sorted(c13_shifts.items()):
        assignment_lines.append(f"c{atom_idx}, {shift:.2f}, {atom_idx + 1}")
    mol.SetProp("NMREDATA_ASSIGNMENT", "\\\n".join(assignment_lines) + "\\")

    # Write to SDF
    sio = StringIO()
    writer = Chem.SDWriter(sio)
    writer.write(mol)
    writer.close()

    return sio.getvalue()
```

### Integration Points

**API endpoint (`src/qm_nmr_calc/api/routers/jobs.py`):**
```python
from qm_nmr_calc.nmredata import generate_nmredata_sdf

@router.get("/jobs/{job_id}/nmredata")
async def get_nmredata(job_id: str):
    """Download NMReData SDF file for completed job."""
    results = load_job_results(job_id)
    mol = load_molecule(job_id)

    sdf_content = generate_nmredata_sdf(
        mol=mol,
        h1_shifts=results["1H_shifts"],
        c13_shifts=results["13C_shifts"],
        solvent=results["metadata"]["solvent"],
    )

    return Response(
        content=sdf_content,
        media_type="chemical/x-mdl-sdfile",
        headers={"Content-Disposition": f"attachment; filename={job_id}.nmredata.sdf"}
    )
```

**Web UI button (`templates/results.html`):**
```html
<a href="/api/jobs/{{ job_id }}/nmredata"
   download="{{ job_id }}.nmredata.sdf"
   class="download-button">
   Download NMReData
</a>
```

## Dependencies

### Already Installed (No Changes)

| Package | Version | Used For |
|---------|---------|----------|
| RDKit | Latest | MOL/SDF I/O, molecule manipulation |
| FastAPI | Latest | REST API endpoints |

### Not Needed

- No additional Python packages required
- No Java/MATLAB tools needed (those are for NMR software integration)
- No external services or APIs

## Testing Strategy

**Unit tests:**
1. Test `generate_nmredata_sdf()` with simple molecules
2. Verify tag formatting (backslashes, line structure)
3. Validate SDF parsing (RDKit can re-read output)
4. Check atom numbering (0-indexed → 1-indexed conversion)

**Integration tests:**
1. Complete job → download NMReData → parse with RDKit
2. Verify all shifts present in ASSIGNMENT tag
3. Verify solvent mapping
4. Verify version and level tags

**End-to-end test:**
1. Submit ethanol (known structure)
2. Wait for completion
3. Download NMReData file
4. Parse and validate structure + shifts

## Validation Against Official Examples

**Reference files:** https://github.com/NMReDATAInitiative/Examples-of-NMR-records

**Validation checklist:**
- [ ] Version tag: `1.1\` (current stable)
- [ ] Assignment format: `label, shift, atom_num\`
- [ ] Backslash line endings in multi-line tags
- [ ] 1-indexed atom numbers (MOL convention)
- [ ] Level 0 for predicted data
- [ ] Solvent in standard format (CDCl3, not CHCl3)
- [ ] File ends with `$$$$`
- [ ] Parseable by RDKit `Chem.SDMolSupplier()`

## Alternative Approaches Considered

### 1. Use Java Tools from NMReDATAInitiative

**Pros:** Official implementation, guaranteed compliance
**Cons:**
- Requires JRE installation
- Subprocess complexity
- No Python API
- Designed for NMR software integration, not prediction output

**Verdict:** Overkill for our use case. We're generating, not validating experimental data.

### 2. Manual SDF String Formatting

**Pros:** No dependencies
**Cons:**
- Duplicates RDKit functionality
- Error-prone (MOL format is complex)
- Harder to maintain
- Need to handle V2000 vs V3000 formats
- Need to handle special characters in tags

**Verdict:** Don't reinvent RDKit's wheel.

### 3. Contribute Python Library to NMReDATAInitiative

**Pros:** Community benefit
**Cons:**
- Scope creep for v2.3 milestone
- Library design requires broader use cases
- Validation/parsing also needed for full library
- Our use case is generation-only

**Verdict:** Consider for future contribution, not for this milestone.

## Risks and Mitigations

| Risk | Impact | Mitigation |
|------|--------|------------|
| Tag format changes in future NMReData versions | Medium | Use version tag, document against 1.1 spec |
| Atom numbering mismatch (0-indexed vs 1-indexed) | High | Explicit +1 conversion, unit test with known examples |
| Multi-line tag formatting errors | Medium | Test against official examples, use RDKit re-parse |
| Solvent name mapping incomplete | Low | Document supported solvents, error on unknown |
| Conformer ensemble averaging | Medium | Use Boltzmann-weighted shifts (already computed) |

## Open Questions

1. **Conformer ensembles:** Should we include conformer populations in metadata?
   - Suggestion: Add custom tag `NMREDATA_QM_CONFORMERS` with count and populations
   - Not standard but informative for predicted data

2. **Spectrum location:** What to put in `Spectrum_Location`?
   - Suggestion: Use job ID (e.g., `job/abc123`) since we don't have raw spectra
   - Or omit entirely (not applicable for predictions)

3. **Interchangeable groups:** Should we detect equivalent atoms?
   - Suggestion: Defer to future enhancement, requires symmetry analysis
   - Low value for QM predictions (assignments are unambiguous)

## Success Criteria

Implementation complete when:
- [ ] `generate_nmredata_sdf()` implemented and tested
- [ ] API endpoint `/jobs/{job_id}/nmredata` working
- [ ] Download button on results page
- [ ] Output validates against official NMReData examples
- [ ] Documentation updated (README, API docs)
- [ ] Unit tests pass
- [ ] Integration test: ethanol job → download → parse

## Timeline Estimate

Based on implementation complexity:

| Task | Effort | Dependencies |
|------|--------|--------------|
| Core NMReData generation function | 4 hours | RDKit familiarity |
| API endpoint | 2 hours | FastAPI patterns |
| Web UI button | 1 hour | Template structure |
| Unit tests | 3 hours | Test framework |
| Integration tests | 2 hours | Test job fixtures |
| Documentation | 2 hours | — |
| **Total** | **14 hours** | **~2 working days** |

Low complexity because:
- No new dependencies
- Straightforward tag formatting
- Existing patterns for downloads (PNG, XYZ)
- Well-documented format

## Sources

### Official NMReData Resources
- [NMReData tag format specification](https://nmredata.org/wiki/NMReDATA_tag_format)
- [NMReDATA Initiative main site](https://www.nmredata.org/)
- [NMReDATA Initiative GitHub](https://github.com/NMReDATAInitiative)
- [Examples repository](https://github.com/NMReDATAInitiative/Examples-of-NMR-records)
- [NMReDATA publication (ResearchGate)](https://www.researchgate.net/publication/324529794_NMReDATA_a_standard_to_report_the_NMR_assignment_and_parameters_of_organic_compounds)

### RDKit Documentation
- [Getting Started with RDKit](https://www.rdkit.org/docs/GettingStartedInPython.html)
- [RDKit properties tutorial](https://greglandrum.github.io/rdkit-blog/posts/2025-01-24-property-tutorial.html)
- [SDWriter documentation](https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html)

### Python NMR Tools (for context)
- [nmrglue documentation](https://www.nmrglue.com/)
- [nmrstarlib publication](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1580-5)

**Research confidence:** HIGH
- Format specification from official sources
- Concrete examples from official repository
- RDKit functionality verified from current documentation
- No Python library confirmed via multiple search strategies
