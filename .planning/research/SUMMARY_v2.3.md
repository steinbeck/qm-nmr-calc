# v2.3 NMReData Export Research Summary

**Project:** qm-nmr-calc v2.3 NMReData Export
**Domain:** Quantum chemistry NMR prediction with standardized file format export
**Researched:** 2026-02-01
**Confidence:** HIGH

## Executive Summary

NMReData is a standardized format for reporting NMR chemical shift assignments by extending SDF (Structure Data Format) files with custom tags prefixed by `NMREDATA_`. The format explicitly supports both experimental and computational/predicted data, making it ideal for exporting qm-nmr-calc's DFT-predicted chemical shifts. Implementation is straightforward because NMReData is a tag formatting convention, not a complex binary format.

The recommended approach is to generate NMReData files on-demand (not pre-generated) using RDKit's existing SDF writing capabilities with property injection via `SetProp()`. This requires no new dependencies since RDKit is already installed and used for molecular structure handling. The core implementation involves formatting predicted 1H and 13C shifts into the `NMREDATA_ASSIGNMENT` tag with proper atom numbering (1-indexed, matching SDF convention), adding required metadata tags (version, solvent, temperature), and writing via RDKit's `SDWriter`.

The primary risk is atom numbering off-by-one errors (RDKit uses 0-indexed atoms, SDF/NMReData uses 1-indexed). This is mitigated through explicit conversion and unit testing. For ensemble calculations, export Boltzmann-weighted average shifts with the lowest-energy conformer geometry, as the format has limited standardized support for multi-conformer data. Implementation complexity is low (estimated 14 hours) with no new dependencies and follows existing download endpoint patterns.

## Key Findings

### Recommended Stack

No new dependencies required. The project already has all necessary tools to generate NMReData-compliant SDF files.

**Core technologies:**
- **RDKit (already installed)**: SDF/MOL generation and property injection — handles MOL block formatting, atom coordinates, and custom tag writing via `SetProp()` and `SDWriter`
- **FastAPI (already installed)**: REST endpoint creation — provides file download response with proper media type `chemical/x-mdl-sdfile`

**Critical finding:** No Python library exists specifically for NMReData generation. The NMReDATA Initiative provides Java/MATLAB tools for NMR software integration but no Python API. RDKit's standard SDF capabilities are sufficient because NMReData is a tag formatting convention on top of SDF.

### Expected Features

NMReData version 1.1 is the current stable format. Version 2.0 adds explicit 3D structure support (recommended for qm-nmr-calc since we export optimized geometries).

**Must have (table stakes):**
- MOL block with 3D coordinates (from optimized geometry XYZ)
- `NMREDATA_VERSION` tag (set to "1.1" or "2.0")
- `NMREDATA_LEVEL` tag (set to "0" for predicted data with no assignment ambiguity)
- `NMREDATA_SOLVENT` tag (map our codes: CHCl3→CDCl3, DMSO→(CD3)2SO, vacuum→vacuum)
- `NMREDATA_TEMPERATURE` tag (298.15 K default or ensemble temperature)
- `NMREDATA_ASSIGNMENT` tag (core data: atom index → chemical shift mapping)

**Should have (competitive):**
- `NMREDATA_FORMULA` (molecular formula from RDKit)
- `NMREDATA_SMILES` (from job input SMILES)
- `NMREDATA_INCHI` (from RDKit for better structure matching)
- `NMREDATA_ID` with provenance metadata (software version, calculation method, basis set, scaling factors, job ID)
- `NMREDATA_1D_1H` and `NMREDATA_1D_13C` pseudo-spectrum tags (enables visualization in NMReData-compatible tools)

**Defer (v2+):**
- J-coupling constants (`NMREDATA_J` tag) — not computed by qm-nmr-calc
- 2D correlation data (`NMREDATA_2D_*` tags) — experimental technique only
- Per-conformer export for ensembles — format has limited multi-conformer support, use Boltzmann-averaged shifts for MVP
- Multiplicity information — requires J-coupling analysis
- Symmetry detection for equivalent atoms — low value for QM predictions (assignments are unambiguous)

### Architecture Approach

Integration follows the existing on-demand file generation pattern used by the `/geometry.sdf` endpoint. Generate NMReData SDF files when requested rather than pre-generating during post-processing to minimize storage overhead and maintain flexibility for format updates.

**Major components:**

1. **New module (`src/qm_nmr_calc/nmredata.py`)** — Core NMReData generation logic
   - `generate_nmredata_sdf()` function takes molecule, geometry, shifts, solvent, metadata
   - Reuses existing XYZ→SDF conversion pattern from jobs.py
   - Uses RDKit's `SetProp()` to inject NMReData tags as molecule properties
   - Returns complete SDF file content as string

2. **API endpoint (`src/qm_nmr_calc/api/routers/jobs.py`)** — Download interface
   - Add `GET /api/v1/jobs/{job_id}/nmredata.sdf` endpoint
   - Validates job completion status (return 409 if incomplete, 404 if not found)
   - Loads nmr_results and optimized geometry from job storage
   - Calls `generate_nmredata_sdf()` and returns as file download
   - Media type: `chemical/x-mdl-sdfile`, filename: `{job_id}_nmredata.sdf`

3. **Tag formatting logic** — Proper NMReData compliance
   - NMREDATA_ASSIGNMENT format: `label, chemical_shift, atom_number\` (one per line)
   - Use 1-indexed atom numbers (convert from RDKit's 0-indexed internal representation)
   - Label format: `h{atom_idx}` for 1H, `c{atom_idx}` for 13C (clearer than generic s0, s1)
   - Separator: ", " (comma + space, required by specification)
   - Chemical shifts with 4+ decimal places (e.g., 7.2453, not 7.24)

**Data flow:**
1. User requests `/jobs/{id}/nmredata.sdf`
2. Load job status and validate completion
3. Load nmr_results (h1_shifts, c13_shifts) and optimized.xyz
4. Create RDKit Mol from SMILES, add explicit hydrogens
5. Set 3D coordinates from XYZ file
6. Inject NMReData tags via `mol.SetProp()`
7. Write to SDF using `Chem.SDWriter()`
8. Return as file download response

**Ensemble handling:** Use Boltzmann-weighted average shifts (already computed) with lowest-energy conformer geometry. Add `NMREDATA_COMMENT` tag indicating ensemble calculation with conformer count and temperature. Individual per-conformer export deferred to future enhancement.

### Critical Pitfalls

1. **Atom numbering off-by-one errors** — RDKit uses 0-indexed atoms internally, but SDF/MOL files and NMReData use 1-indexed atom numbering. Failing to convert causes all atom assignments to be wrong. **Prevention:** Explicit +1 conversion when writing atom indices to `NMREDATA_ASSIGNMENT` tag. Unit test that first atom is index 1, not 0. Integration test: export and re-import, verify assignments match original.

2. **Separator format inconsistency** — NMReData specification requires ", " (comma + space) as separator in tags. Using "," (comma only) or other variations causes parsing failures in downstream tools. **Prevention:** Use constant `NMREDATA_SEP = ", "` throughout code. Test exported files with regex check: `r",(?! )"` should find zero matches. Validate with NMReData Java parser if available.

3. **Implicit hydrogen atom references** — NMReData requires explicit atom references but typically molecules have implicit hydrogens. Must decide whether to export with explicit H atoms in MOL block or use "H{atom_number}" notation for implicit hydrogens. **Prevention:** Decide on approach (explicit H recommended: use `Chem.AddHs(mol)` before SDF generation). Document choice. Test with molecules having both aliphatic C-H and heteroatom-H bonds.

4. **Conformer ensemble representation not standardized** — NMReData format (versions 1.0-2.0) has limited support for multiple conformers. No standard tags for Boltzmann populations or per-conformer shifts. **Prevention:** Use Strategy A (recommended): Export Boltzmann-weighted average shifts only, single structure (lowest energy conformer), add comment noting ensemble averaging. This preserves simplicity and matches most NMReData use cases.

5. **Missing calculation method metadata** — NMReData format lacks standard tags to distinguish predicted from experimental data or document calculation methodology. Files without metadata may confuse users. **Prevention:** Use `NMREDATA_ID` tag with custom key-value pairs for provenance (software version, DFT method, basis set, scaling factors, job ID). Add comment fields indicating predicted data source.

## Implications for Roadmap

Based on research, suggested phase structure for v2.3 milestone:

### Phase 1: Core NMReData Generation Module
**Rationale:** Foundational functionality must be implemented first before API integration. Low-risk implementation since it only involves data formatting, no external dependencies or complex architecture changes.

**Delivers:**
- `src/qm_nmr_calc/nmredata.py` module with `generate_nmredata_sdf()` function
- Tag formatting logic for all required NMReData fields
- XYZ to MOL conversion with 3D coordinates
- Atom numbering conversion (0-indexed → 1-indexed)
- Solvent name mapping (CHCl3→CDCl3, etc.)

**Addresses features:**
- MOL block generation with optimized geometry
- NMREDATA_VERSION, NMREDATA_LEVEL, NMREDATA_SOLVENT, NMREDATA_TEMPERATURE tags
- NMREDATA_ASSIGNMENT tag with proper formatting
- NMREDATA_FORMULA and NMREDATA_SMILES tags

**Avoids pitfalls:**
- Atom numbering off-by-one (explicit conversion and unit tests)
- Separator inconsistency (use constant from start)
- Implicit hydrogen handling (use explicit H atoms via `AddHs()`)

**Research flag:** Standard implementation, no additional research needed (well-documented format specification and RDKit patterns)

---

### Phase 2: API Endpoint Integration
**Rationale:** API endpoint implementation depends on core generation module from Phase 1. Follows existing download endpoint patterns (`/geometry.sdf`, `/spectrum/*.png`).

**Delivers:**
- `GET /api/v1/jobs/{job_id}/nmredata.sdf` endpoint in jobs.py router
- Job validation logic (404 if not found, 409 if incomplete)
- File download response with proper media type and headers
- Integration with existing job storage and results loading

**Uses stack:**
- FastAPI routing and response classes
- Existing job validation patterns
- Core module from Phase 1

**Implements architecture:**
- On-demand generation (no pre-generated files)
- Same pattern as existing `/geometry.sdf` endpoint
- Stateless file generation

**Avoids pitfalls:**
- No new error modes (reuses existing job loading and validation)
- Performance acceptable (SDF generation <5ms, similar to existing endpoint)

**Research flag:** Standard implementation, no additional research needed (follows existing patterns)

---

### Phase 3: Testing and Validation
**Rationale:** Comprehensive testing required to catch atom numbering errors, format compliance issues, and edge cases. Must validate against NMReData specification and test round-trip parsing.

**Delivers:**
- Unit tests for `nmredata.py` (tag formatting, atom indexing, solvent mapping)
- Integration tests for API endpoint (complete/incomplete jobs, file parsing)
- Validation tests (parse exported file with RDKit, verify atom assignments)
- Test data with known molecules (ethanol, benzene, complex structures)

**Addresses pitfalls:**
- Atom numbering verification (first atom is 1, not 0)
- Format compliance (parseable by RDKit SDMolSupplier)
- Separator consistency (regex validation)
- Missing VERSION tag detection

**Research flag:** May need external validation tools (NMReData Java parser from NMReDATAInitiative) to verify full compliance beyond RDKit parsing

---

### Phase 4: Optional Enhancements
**Rationale:** Post-MVP enhancements add value but aren't required for basic functionality. Can be implemented incrementally without breaking existing exports.

**Delivers:**
- `NMREDATA_ID` tag with full provenance metadata (software, method, basis set, scaling factors, job ID)
- `NMREDATA_INCHI` tag for better structure matching
- `NMREDATA_1D_1H` and `NMREDATA_1D_13C` pseudo-spectrum tags for visualization tool compatibility
- Ensemble metadata comment (conformer count, temperature, Boltzmann weighting)
- Documentation and examples

**Addresses features:**
- Provenance tracking and reproducibility
- Visualization in NMReData-compatible tools (Mestrelab Mnova, ChemInfo viewer)
- Clear distinction between predicted and experimental data

**Research flag:** Standard implementation for most features. May need research for visualization tool compatibility testing (access to Mnova or ChemInfo viewer)

---

### Phase Ordering Rationale

- **Phase 1 before Phase 2:** Core generation logic must be implemented and tested before API integration. Allows unit testing of formatting without HTTP layer complexity.

- **Phase 3 overlaps with Phase 1-2:** Testing can begin during core development. Integration tests require Phase 2 completion but unit tests run against Phase 1 components.

- **Phase 4 deferred:** Enhancements add value but aren't blocking for basic NMReData export. Can be implemented after MVP is validated. Allows iterative improvement based on user feedback.

- **Dependencies:** Phase 2 depends on Phase 1. Phase 3 tests Phase 1 and Phase 2. Phase 4 is independent (enhancements only).

- **Risk mitigation:** Critical pitfalls (atom numbering, separators) addressed in Phase 1 with immediate testing in Phase 3. Ensemble representation complexity deferred via Strategy A (averaged shifts only).

### Research Flags

Phases likely needing deeper research during planning:
- **Phase 3 (Testing):** May need access to external NMReData validation tools (Java parser from NMReDATAInitiative) or chemistry software (Mestrelab Mnova, ChemInfo viewer) to verify full format compliance beyond RDKit parsing. Research how to obtain/run these tools for CI integration.

- **Phase 4 (Enhancements):** If implementing visualization tool compatibility, research which software to target (Mnova vs ChemInfo vs others) and what features they require in NMReData files. May need examples of "pseudo-spectrum" tags for predicted data.

Phases with standard patterns (skip research-phase):
- **Phase 1 (Core Module):** Well-documented format specification, established RDKit patterns, no unknowns
- **Phase 2 (API Endpoint):** Follows existing download endpoint pattern, no new architectural decisions

## Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| Stack | HIGH | RDKit capabilities verified from current documentation and existing codebase usage. No new dependencies needed. |
| Features | HIGH | NMReData specification version 1.1/2.0 stable and well-documented with official examples. Required tags clearly specified. |
| Architecture | HIGH | On-demand generation pattern already proven in `/geometry.sdf` endpoint. Implementation approach straightforward. |
| Pitfalls | MEDIUM | Critical pitfalls identified from format spec and community discussions. Atom numbering validated from SDF spec. Conformer ensemble limitation clear but workaround straightforward. Validation tooling limited (no widely-adopted Python validator). |

**Overall confidence:** HIGH

Implementation is low-risk with no new dependencies and follows existing patterns. Format specification is stable and well-documented with official examples from NMReDATAInitiative repository.

### Gaps to Address

- **Atom ordering consistency:** Current assumption is that `MolFromSmiles()` → `AddHs()` produces same atom ordering as initial XYZ generation used by NWChem. This assumption is inherited from existing `/geometry.sdf` endpoint. Validate empirically with test cases during Phase 3. If discrepancies found, may need to store explicit atom mapping during initial geometry generation (architectural change for future milestone).

- **Validation tooling availability:** Limited Python-based NMReData validators available. NMReDATAInitiative provides Java tools but require JRE installation. During Phase 3, determine if Java validator should be added to CI pipeline or if RDKit round-trip parsing is sufficient. Document validation approach.

- **Ensemble conformer representation:** NMReData format has limited standardized support for multi-conformer data. Strategy A (Boltzmann-averaged shifts with single structure) recommended for MVP but may need enhancement if users request per-conformer export. Monitor NMReDATAInitiative GitHub for conformer tag proposals in future versions.

- **Visualization tool testing:** Phase 4 enhancements assume compatibility with Mestrelab Mnova and ChemInfo viewer. Actual testing may require software licenses or online services. Determine during Phase 4 planning whether this validation is feasible or if format compliance based on specification is sufficient.

## Sources

### Primary (HIGH confidence)
- [NMReData Tag Format Specification v1.1](https://nmredata.org/wiki/NMReDATA_tag_format) — Required tags, ASSIGNMENT format, separator rules
- [NMReData Tag Format v2.0](https://nmredata.org/wiki/NMReDATA_tag_format_2.0) — 3D structure support
- [NMReData Examples Repository](https://github.com/NMReDATAInitiative/Examples-of-NMR-records) — Official example files including computational data
- [NMReData Publication (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6226248/) — Format specification and design rationale
- [RDKit Getting Started](https://www.rdkit.org/docs/GettingStartedInPython.html) — SetProp() and SDWriter documentation
- [RDKit Properties Tutorial](https://greglandrum.github.io/rdkit-blog/posts/2025-01-24-property-tutorial.html) — Custom property injection in SDF files

### Secondary (MEDIUM confidence)
- [NMReData Object Structure](https://nmredata.org/wiki/Nmredata_object_structure) — Multi-conformer discussions
- [NMReData Future Version](https://nmredata.org/wiki/Future_version) — Conformer support proposals
- [SDF File Format Guidance](https://www.nonlinear.com/progenesis/sdf-studio/v0.9/faq/sdf-file-format-guidance.aspx) — 1-indexed atom numbering convention
- [Chemical Table File (Wikipedia)](https://en.wikipedia.org/wiki/Chemical_table_file) — SDF/MOL format background
- Existing qm-nmr-calc codebase — `/geometry.sdf` endpoint pattern in `src/qm_nmr_calc/api/routers/jobs.py`

### Tertiary (LOW confidence)
- [ChemInfo NMReData Viewer](https://www.cheminfo.org/flavor/nmredata/viewer/) — Potential validation tool, accessibility unclear
- [NMReDATAInitiative Java Tools](https://github.com/NMReDATAInitiative/javatools) — Validation parser, requires JRE setup

---
*Research completed: 2026-02-01*
*Ready for roadmap: yes*
