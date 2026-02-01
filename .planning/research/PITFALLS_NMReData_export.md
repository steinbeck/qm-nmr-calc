# NMReData Export Implementation Pitfalls

**Domain:** NMReData export for predicted NMR data
**Researched:** 2026-02-01
**Confidence:** MEDIUM (based on format specification, community discussions, and web sources)

## Executive Summary

NMReData export for PREDICTED (not experimental) NMR data has unique challenges not commonly addressed in the literature, which focuses primarily on experimental spectra. The critical pitfalls involve atom numbering conversion (0-indexed to 1-indexed), hydrogen atom handling (implicit vs explicit), conformer ensemble representation (currently limited format support), and proper metadata to distinguish predicted from experimental data.

## Critical Pitfalls

Mistakes that cause invalid files, data loss, or major rewrites.

### Pitfall 1: Atom Numbering Off-By-One Errors

**What goes wrong:** RDKit uses 0-indexed atom numbering internally, but SDF/MOL files and NMReData use 1-indexed atom numbering. Failing to convert causes all atom assignments to be wrong.

**Why it happens:**
- Different indexing conventions between Python/RDKit (0-indexed) and chemistry file formats (1-indexed)
- Bond blocks, isotope entries, and NMREDATA_ASSIGNMENT tags all use 1-indexed atom numbers
- Easy to miss during implementation when data "looks correct" but is shifted by one

**Consequences:**
- All NMR assignments point to wrong atoms
- Chemical shift for H1 appears on H2, etc.
- Validation tools may not catch this (syntactically valid but chemically wrong)
- Silent data corruption that may not be noticed until manual inspection

**Prevention:**
```python
# WRONG - using RDKit atom index directly
rdkit_atom_idx = 0  # First atom
nmredata_line = f"1H, {shift}, {rdkit_atom_idx}"  # Points to atom 0 (invalid)

# RIGHT - convert to 1-indexed
rdkit_atom_idx = 0  # First atom
sdf_atom_num = rdkit_atom_idx + 1
nmredata_line = f"1H, {shift}, {sdf_atom_num}"  # Points to atom 1 (correct)
```

**Detection:**
- Unit test: Create molecule with known atom numbering, export NMReData, verify first atom is "1" not "0"
- Integration test: Export and re-import, verify assignments match original
- Visual inspection: Check that 3D coordinates in MOL block align with assignment indices

**Phase impact:** Implementation phase (Phase 3) - this must be caught in initial development, not post-deployment

---

### Pitfall 2: Implicit Hydrogen Atom References

**What goes wrong:** NMReData requires explicit references to hydrogen atoms using notation like "H3" (hydrogen on atom 3) when hydrogens are implicit in the structure. Failing to handle this correctly causes missing or incorrect 1H assignments.

**Why it happens:**
- RDKit stores molecules with implicit hydrogens by default
- SDF files typically don't include explicit hydrogen atoms
- NMReData needs to reference hydrogens that don't exist as separate atoms in the MOL block
- The "H3" notation is specific to NMReData format

**Consequences:**
- 1H chemical shifts cannot be assigned (no atom to reference)
- File appears valid but contains no 1H data
- Confusion about which hydrogen on a given heavy atom (if multiple hydrogens present)

**Prevention:**
- Use NMReData "H{atom_number}" notation for implicit hydrogens
- For 13C shifts: reference heavy atom number directly
- For 1H shifts on implicit H: use "H{heavy_atom_number}" format
- Document whether your export includes explicit H atoms or uses implicit notation

**Example:**
```
# For implicit H (typical case):
# Molecule: CH3-CH2-OH (no H atoms in MOL block)
# Atom 1 = C (methyl)
# Atom 2 = C (methylene)
# Atom 3 = O

NMREDATA_ASSIGNMENT:
H1, 1.23, H1    # 1H shift for hydrogen on carbon 1
H2, 3.45, H2    # 1H shift for hydrogen on carbon 2
C1, 14.5, 1     # 13C shift for carbon 1 (direct reference)
C2, 63.2, 2     # 13C shift for carbon 2
```

**Detection:**
- Test with molecule that has both aliphatic C-H and heteroatom-H bonds
- Verify 1H assignments use "H{n}" format when hydrogens are implicit
- Check that number of 1H assignments matches number of hydrogen atoms

**Phase impact:** Implementation phase (Phase 3) - format specification must be understood correctly from the start

---

### Pitfall 3: Separator Format Inconsistency

**What goes wrong:** NMReData specification uses ", " (comma + space) as separator in tags, but implementations may use "," (comma only) or other variations. This causes parsing failures in downstream tools.

**Why it happens:**
- Developer discussion (Bruker's Pavel Kessler) noted challenges with comma separators
- Space after comma is currently required but "may become optional in future versions"
- Different implementations make different assumptions
- Easy to overlook whitespace in format specification

**Consequences:**
- Files fail to parse in some NMReData readers
- Interoperability issues between software platforms
- Need to maintain compatibility with multiple parser versions

**Prevention:**
- Use ", " (comma + space) consistently as current specification requires
- Use a constant/variable in code for separator: `NMREDATA_SEP = ", "`
- Allows easy switching if format changes to comma-only
- Test exported files with available NMReData parsers (Java tools from NMReDATAInitiative)

**Example:**
```python
# GOOD - use constant
NMREDATA_SEP = ", "
assignment = f"H1{NMREDATA_SEP}{shift}{NMREDATA_SEP}{atom_num}"

# BAD - hardcoded inconsistently
assignment1 = f"H1, {shift}, {atom_num}"  # comma-space
assignment2 = f"H2,{shift},{atom_num}"    # comma-only (inconsistent)
```

**Detection:**
- String matching test: verify all commas in output are followed by exactly one space
- Regex check: `r",(?! )"` should find zero matches in output
- Parser validation: run exported file through NMReData Java tools if available

**Phase impact:** Implementation phase (Phase 3) - establish pattern early to avoid inconsistency

---

### Pitfall 4: Missing or Incorrect VERSION Tag

**What goes wrong:** The NMREDATA_VERSION tag is required and specifies format version (1.0, 1.1, or 2.0). Missing or incorrect version causes validation failures and incorrect parsing.

**Why it happens:**
- Not clearly documented as mandatory in all sources
- Different versions have different tag support
- Easy to forget during initial implementation

**Consequences:**
- Parsers may reject file entirely
- File parsed with wrong version expectations (missing features or incorrect interpretation)
- No clear indication of format capabilities to downstream tools

**Prevention:**
- Always include NMREDATA_VERSION tag
- Use version 1.1 (current stable) or 2.0 (if using 3D structures)
- Version 2.0 required if including 3D molecular geometry
- Document version choice in code and user-facing docs

**Example:**
```
> <NMREDATA_VERSION>
1.1

# Or for 3D structure support:
> <NMREDATA_VERSION>
2.0
```

**Detection:**
- File validation: Check that VERSION tag exists
- Version test: Verify version matches features used (2.0 if 3D coords present)
- Parser test: Feed to NMReData tools and check for version warnings

**Phase impact:** Implementation phase (Phase 3) - easy fix if caught early, painful if discovered late

---

### Pitfall 5: Conformer Ensemble Representation Not Supported

**What goes wrong:** NMReData format (versions 1.0-2.0) has limited support for multiple conformers. Current versions don't have standard tags for Boltzmann populations or per-conformer shifts. Trying to export ensemble data without a clear strategy causes data loss.

**Why it happens:**
- qm-nmr-calc uses Boltzmann-weighted ensemble averaging
- NMReData was designed primarily for experimental data (single structure)
- Future versions discuss conformer support but not yet standardized
- No clear community consensus on representing predicted ensemble data

**Consequences:**
- Loss of per-conformer shift information
- No way to represent Boltzmann populations
- Inability to validate ensemble calculation methodology
- Users cannot reproduce weighted averaging

**Current format limitations:**
- Version 2.0 discusses 3D structures but mentions "conformations.sdf" as separate file (non-standard)
- Main .sdf file contains single 3D structure (typically lowest energy conformer)
- No standard tags for: conformer count, energies, populations, per-conformer shifts

**Prevention:**
- **Strategy A (Recommended):** Export Boltzmann-weighted average shifts only
  - Single structure in MOL block (lowest energy conformer OR flat 2D structure)
  - NMREDATA_ASSIGNMENT contains weighted-average chemical shifts
  - Add comment field noting "Boltzmann-weighted average over N conformers"
  - Preserves simplicity, matches most NMReData use cases

- **Strategy B (If detailed provenance needed):** Custom tags with documentation
  - Use standard NMReData for averaged data
  - Add custom tags: `NMREDATA_CONFORMER_COUNT`, `NMREDATA_CONFORMER_POPULATIONS`
  - Document custom tags in README/spec
  - Note: May not be parseable by standard tools

- **Strategy C (Future-proofing):** Wait for version 2.x conformer support
  - Monitor NMReDATAInitiative GitHub for conformer tag proposals
  - Implement when standardized
  - Risk: May not happen on project timeline

**Detection:**
- Design decision in planning phase: Choose strategy before implementation
- Documentation: Clearly state which strategy used and why
- Test: Verify ensemble calculations produce single set of averaged shifts for export

**Phase impact:** Planning phase (Phase 1-2) - architectural decision affects entire implementation

---

## Moderate Pitfalls

Mistakes that cause delays or technical debt but are fixable.

### Pitfall 6: No Indication of Predicted vs Experimental Data

**What goes wrong:** NMReData format doesn't have a standard tag to distinguish predicted/calculated data from experimental data. Files lack metadata about calculation methodology.

**Why it happens:**
- Format designed primarily for experimental spectra
- NMREDATA_LEVEL tag indicates assignment ambiguity, not data source
- No standard tag for "calculation method" or "data source"

**Consequences:**
- Users may confuse predicted shifts with experimental measurements
- No provenance information (DFT method, basis set, scaling factors)
- Difficult to compare predicted vs experimental when both exist
- Lost opportunity to document methodology

**Prevention:**
- Use comment fields (text after ";") to add metadata
- Consider custom tag like `NMREDATA_CALCULATION_METHOD` with documentation
- Include in comments: DFT functional, basis set, solvent model, scaling factors
- Set NMREDATA_LEVEL=0 (complete assignment, no ambiguities) since predictions are deterministic

**Example:**
```
> <NMREDATA_VERSION>
2.0

> <NMREDATA_LEVEL>
0; Predicted data (no assignment ambiguities)

> <NMREDATA_CALCULATION_METHOD>
B3LYP/6-311+G(2d,p) with COSMO solvation, DELTA50-derived scaling factors

> <NMREDATA_ASSIGNMENT>
H1, 1.234, H1; Boltzmann-weighted average over 8 conformers
C1, 14.56, 1; Boltzmann-weighted average over 8 conformers
```

**Detection:**
- Code review: Check that methodology is documented somewhere in output
- User testing: Can users tell this is predicted data?

**Phase impact:** Implementation phase (Phase 3) - add metadata from the start

---

### Pitfall 7: Scaling Factor Application Timing

**What goes wrong:** Applying scaling factors after NMReData export vs before affects whether exported shifts are "raw" shieldings or "scaled" chemical shifts. Inconsistency causes confusion.

**Why it happens:**
- qm-nmr-calc applies DELTA50-derived scaling factors
- Unclear whether NMReData should contain shieldings or chemical shifts
- Different use cases prefer different data

**Consequences:**
- Exported shifts don't match web UI display
- Users can't reproduce scaling from raw data
- Mixing shieldings and shifts in same file

**Prevention:**
- Export chemical shifts (after scaling), not raw shieldings
- NMReData ASSIGNMENT tag contains final predicted chemical shifts
- Document scaling factors in comments or custom tag
- If raw data needed, provide separate download (not NMReData)

**Rationale:**
- NMReData designed for "chemical shift values" (experimental-like)
- Scaling factors convert shieldings → shifts (closer to experimental values)
- Exporting scaled shifts matches user expectations

**Detection:**
- Compare exported shifts to web UI display values (should match)
- Unit test: Verify exported values equal `scaled_shifts`, not `raw_shieldings`

**Phase impact:** Implementation phase (Phase 3) - clarify before implementing export logic

---

### Pitfall 8: Line Continuation Backslash Handling

**What goes wrong:** NMReData tags can span multiple lines using backslash ("\\") before end-of-line. Forgetting backslashes causes truncated data; adding them when not needed causes parser errors.

**Why it happens:**
- Backslash is optional for single-line tags
- Mandatory for multi-line tags
- Easy to add incorrectly or forget when needed

**Consequences:**
- Data truncated (parser stops reading at first newline)
- Parser errors if backslash used incorrectly
- Inconsistent file formatting

**Prevention:**
- For single-line tags: Don't use backslash
- For multi-line tags: Add "\\\\" before each newline (except last)
- Most tags will be single-line for typical molecules

**Example:**
```
# Single line (no backslash):
> <NMREDATA_ASSIGNMENT>
H1, 1.23, H1, C1, 14.5, 1

# Multi-line (backslash required):
> <NMREDATA_ASSIGNMENT>
H1, 1.23, H1\\
C1, 14.5, 1\\
C2, 20.3, 2
```

**Detection:**
- Parser test: Verify multi-line tags parse correctly
- String test: Count newlines and backslashes in tag values

**Phase impact:** Implementation phase (Phase 3) - handle during tag formatting

---

## Minor Pitfalls

Mistakes that cause annoyance but are easily fixable.

### Pitfall 9: Invalid Tag Names

**What goes wrong:** NMReData tag names must follow SDF format restrictions: cannot contain hyphen (-), period (.), <, >, =, %, or space. Must begin with alpha character.

**Why it happens:**
- Natural to use hyphens in tag names (e.g., "NMREDATA_J-COUPLING")
- Not obvious from casual reading of spec

**Consequences:**
- SDF parsers reject file
- Tags silently ignored

**Prevention:**
- Use underscores instead of hyphens: `NMREDATA_J` not `NMREDATA_J-COUPLING`
- Validate tag names against regex: `^[A-Za-z][A-Za-z0-9_]*$`
- Standard tags are already compliant; issue only for custom tags

**Detection:**
- Tag validation function at file write time
- Parser test with strict SDF validator

**Phase impact:** Implementation phase (Phase 3) - easy to catch with validation

---

### Pitfall 10: Comment Handling with Semicolons

**What goes wrong:** Everything after ";" in NMReData tag values is treated as a comment and ignored by parsers. Data containing semicolons gets truncated.

**Why it happens:**
- Semicolon is comment delimiter in NMReData
- Chemical names or SMILES strings might contain semicolons
- Not obvious that data will be truncated

**Consequences:**
- Data loss (everything after ; is ignored)
- Unexpected truncation in labels or metadata

**Prevention:**
- Avoid semicolons in label names and data fields
- Use semicolons only for actual comments
- If semicolons needed in data, escape or quote (check spec for escaping rules)

**Example:**
```
# BAD - semicolon truncates data:
H1, 1.23, H1, source=experimental; solvent=CDCl3
# Parser sees: "H1, 1.23, H1, source=experimental"
# Loses: " solvent=CDCl3"

# GOOD - semicolon marks actual comment:
H1, 1.23, H1; Boltzmann-weighted average
# Parser sees: "H1, 1.23, H1"
# Comment preserved for human readers
```

**Detection:**
- Code review: Check that semicolons only appear before comments
- Test: Parse exported file and verify all data fields intact

**Phase impact:** Implementation phase (Phase 3) - easy fix if noticed

---

## Phase-Specific Warnings

| Phase Topic | Likely Pitfall | Mitigation |
|-------------|---------------|------------|
| Planning/Architecture | Conformer ensemble representation | Choose Strategy A (averaged shifts) unless strong need for per-conformer data |
| SDF generation | Atom numbering off-by-one | Unit test: first atom is index 1, not 0 |
| SDF generation | Implicit vs explicit hydrogens | Decide on implicit (use "H{n}" notation) and test |
| Tag formatting | Separator inconsistency | Use constant `NMREDATA_SEP = ", "` from start |
| Tag content | Missing VERSION tag | Add VERSION tag in template/boilerplate |
| Tag content | No calculation metadata | Add comments or custom tags for provenance |
| Validation | No validation tool available | Manual testing + Java parser from NMReDATAInitiative if possible |
| Testing | Silent data corruption | Integration test: export and re-import, verify atom assignments |

---

## Validation Strategy

### Official Validation Tools

**Status:** LIMITED availability

**What exists:**
- Java tools in `NMReDATAInitiative/javatools` repository
- Example files in `NMReDATAInitiative/Examples-of-NMR-records`
- No widely-adopted standalone validator tool
- Initiative plans "validation service" (not yet available)

**What this means:**
- Validation is primarily manual
- Testing against example files is critical
- Parser errors are your main feedback mechanism

### Recommended Validation Approach

**Level 1: Format Validation**
- SDF syntax (use RDKit to read back exported file)
- Required tags present (VERSION, ASSIGNMENT)
- Tag names valid (no special characters)
- Separator format consistent (comma-space)

**Level 2: Semantic Validation**
- Atom numbers in valid range (1 to num_atoms)
- Hydrogen references valid ("H{n}" where n is valid heavy atom)
- Chemical shift values reasonable (not 999999 or NaN)
- Number of assignments matches expected (one per NMR-active nucleus)

**Level 3: Round-Trip Validation**
- Export NMReData file
- Parse with SDF reader
- Verify atom assignments preserved
- Check chemical shifts match original

**Level 4: Manual Inspection**
- Visual check in text editor
- Compare to example files from NMReDATAInitiative
- Test in chemistry software that supports NMReData (MNova if available)

### Test Data Strategy

- Start with small molecules (5-10 heavy atoms)
- Include heteroatoms (N, O) to test different cases
- Test conformer averaging vs single conformer
- Test all supported solvents
- Use molecules from DELTA50 benchmark (known reference)

---

## qm-nmr-calc Specific Concerns

Based on project context from PROJECT.md:

### Concern 1: RDKit Atom Numbering (0-indexed → 1-indexed)

**Issue:** RDKit uses 0-indexed atoms internally. SDF/NMReData uses 1-indexed.

**Where to convert:**
- When writing atom numbers to NMREDATA_ASSIGNMENT tag
- When writing bond indices to MOL block (already handled by RDKit's SDF writer)
- When mapping shift results (keyed by RDKit atom idx) to SDF atom numbers

**Test case:**
```python
# mol with 5 heavy atoms (RDKit idx 0-4)
# shift_dict = {0: 1.23, 1: 3.45, 2: 5.67, 3: 7.89, 4: 2.34}
# NMReData should reference atoms 1-5, not 0-4
assert "H1, 1.23, H1" in nmredata_assignment  # RDKit atom 0 → SDF atom 1
```

---

### Concern 2: Ensemble vs Single Conformer Modes

**Issue:** qm-nmr-calc supports both modes. NMReData export should handle both.

**Strategy:**
- Single conformer: Use single structure in MOL block, direct shift values
- Ensemble: Use averaged shifts, note in comments, include single 3D structure (lowest energy conformer)

**Metadata to include:**
- Single: "Single conformer calculation"
- Ensemble: "Boltzmann-weighted average over {n} conformers"

---

### Concern 3: Solvent-Specific Scaling Factors

**Issue:** DELTA50 scaling factors vary by solvent. Exported shifts reflect solvent choice.

**Metadata to include:**
- Solvent name (CHCl3, DMSO, vacuum)
- Scaling factors used (slope, intercept)
- Expected MAE (1H: 0.12 ppm, 13C: 2.0 ppm)

**Example comment:**
```
> <NMREDATA_CALCULATION_METHOD>
B3LYP/6-311+G(2d,p) with COSMO solvation (CHCl3), DELTA50-derived scaling factors (1H: slope=0.98, intercept=0.12; 13C: slope=1.02, intercept=-2.5), Expected MAE: 1H 0.12 ppm, 13C 2.0 ppm
```

---

### Concern 4: 3D Structure Export

**Issue:** qm-nmr-calc returns optimized 3D geometries. NMReData v2.0 supports 3D structures in MOL block.

**Recommendation:**
- Use VERSION 2.0 (required for 3D)
- Include optimized geometry (DFT-optimized)
- For ensemble: use lowest energy conformer OR Boltzmann-weighted geometry centroid
- Document which structure is included

---

### Concern 5: Missing J-Coupling Data

**Issue:** qm-nmr-calc predicts chemical shifts only, not J-coupling constants. NMReData has NMREDATA_J tag for couplings.

**Recommendation:**
- Omit NMREDATA_J tag (it's optional)
- If future versions add J-coupling: add tag at that time
- No validation issues from missing optional tags

---

## Implementation Checklist

Before shipping NMReData export:

**Format Compliance:**
- [ ] NMREDATA_VERSION tag present and correct (2.0 recommended)
- [ ] Tag names use valid characters (no hyphens, periods, etc.)
- [ ] Separator is ", " (comma + space) throughout
- [ ] Comments use ";" delimiter correctly
- [ ] Backslash continuation used correctly for multi-line tags

**Atom Numbering:**
- [ ] Atom numbers are 1-indexed (not 0-indexed)
- [ ] Hydrogen references use "H{n}" format for implicit H
- [ ] All atom numbers in valid range (1 to num_atoms)
- [ ] Unit test verifies first atom is 1

**Data Content:**
- [ ] Chemical shifts (not shieldings) exported
- [ ] DELTA50 scaling factors applied before export
- [ ] Boltzmann averaging applied (if ensemble mode)
- [ ] 3D structure included (if VERSION 2.0)

**Metadata:**
- [ ] Calculation method documented (functional, basis set, solvent)
- [ ] Scaling factors documented
- [ ] Ensemble metadata (conformer count, populations) if applicable
- [ ] Clear indication this is PREDICTED data (not experimental)

**Testing:**
- [ ] Round-trip test: export and re-import, verify assignments
- [ ] Visual inspection: compare to NMReDATAInitiative examples
- [ ] Parser test: read with Java tools if available
- [ ] Integration test: verify matches web UI output

**Documentation:**
- [ ] README documents NMReData export feature
- [ ] API docs describe endpoint and response format
- [ ] Known limitations documented (e.g., no J-coupling, averaged ensemble)

---

## Confidence Assessment

| Finding | Confidence | Source |
|---------|-----------|--------|
| Atom numbering 1-indexed | HIGH | NMReData wiki, SDF spec (multiple sources) |
| Separator format (comma-space) | HIGH | NMReData tag format spec, developer discussion |
| Implicit H notation "H{n}" | HIGH | NMReData tag format spec |
| VERSION tag required | MEDIUM | Parser documentation (not explicit in all sources) |
| Conformer ensemble limitations | MEDIUM | Future version discussions, no current standard |
| No predicted vs experimental tag | MEDIUM | Format spec review (absence of such tag) |
| Validation tool availability | HIGH | NMReDATAInitiative GitHub repositories |

**Overall confidence: MEDIUM**

Confidence limited by:
- Official specification spread across wiki pages (not single authoritative document)
- Limited tooling for validation
- Format designed for experimental data (predicted data use case less documented)
- No official examples of predicted/calculated data exports
- Conformer ensemble representation not standardized

**Recommendations to increase confidence:**
1. Test with Java tools from NMReDATAInitiative/javatools
2. Compare exports to files in NMReDATAInitiative/Examples-of-NMR-records
3. Reach out to NMReData community (GitHub, mailing list) for predicted data guidance
4. Start with simple cases, expand gradually

---

## Sources

### Official NMReData Documentation
- [NMReData Main Page](https://nmredata.org/wiki/Main_Page)
- [NMReData Tag Format Specification](https://nmredata.org/wiki/NMReDATA_tag_format)
- [NMReData Tag Format 2.0](https://nmredata.org/wiki/NMReDATA_tag_format_2.0)
- [NMReData Object Structure](https://nmredata.org/wiki/Nmredata_object_structure)
- [NMReData Future Version Discussions](https://nmredata.org/wiki/Future_version)

### GitHub Resources
- [NMReDATAInitiative Organization](https://github.com/NMReDATAInitiative)
- [NMReDATAInitiative Main Website Repo](https://nmredatainitiative.github.io/Main_Page/)
- [cheminfo/nmredata Issues](https://github.com/cheminfo/nmredata/issues)

### Publications
- [NMReData Standard Publication (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6226248/)
- [NMReData Standard Publication (PubMed)](https://pubmed.ncbi.nlm.nih.gov/29656574/)
- [NMReData Tools and Applications (PubMed)](https://pubmed.ncbi.nlm.nih.gov/33729627/)
- [NMReData PDF Specification](https://kups.ub.uni-koeln.de/8163/4/MRC_main_V25.pdf)

### SDF/MOL Format
- [SDF File Format Guidance](https://www.nonlinear.com/progenesis/sdf-studio/v0.9/faq/sdf-file-format-guidance.aspx)
- [Chemical Table File (Wikipedia)](https://en.wikipedia.org/wiki/Chemical_table_file)
- [MOL Format V2000 Tutorial](https://tutorials.technology/cheatsheets/mol_format.html)

### RDKit Documentation
- [RDKit Getting Started Guide](https://www.rdkit.org/docs/GettingStartedInPython.html)
- [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html)
- [RDKit Issue 2362: Atom Parity from SDF with Implicit H](https://github.com/rdkit/rdkit/issues/2362)

### Computational NMR and Boltzmann Weighting
- [Nature Protocols: NMR Chemical Shift Calculation Guide](https://www.nature.com/articles/nprot.2014.042)
- [Automated Framework for NMR Chemical Shift Calculations (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6755567/)
- [Evaluating Error in Boltzmann Weighting (Corin Wagen)](https://corinwagen.github.io/public/blog/20221228_boltzmann_error.html)

---

**Research completed:** 2026-02-01
**Next step:** Use these findings to inform NMReData export roadmap planning
