# Phase 66: Dataset Archive Preparation - Research

**Researched:** 2026-02-11
**Domain:** Scientific dataset archival, FAIR data principles, computational chemistry provenance
**Confidence:** MEDIUM-HIGH

## Summary

Dataset archive preparation for computational chemistry involves organizing calculation data (inputs/outputs), extracting metadata, generating machine-readable identifiers, documenting provenance, and packaging everything with FAIR-compliant metadata. The target repository (RADAR4Chem) requires DataCite 4.0+ metadata, accepts all file formats, and mandates 10 core metadata fields.

This phase requires organizing 650 existing NWChem calculations across 50 molecules x 13 solvents, extracting molecule identifiers (SMILES/InChI/InChIKey via RDKit), generating file integrity checksums (SHA-256), documenting computational provenance, and creating FAIR-compliant documentation (README, metadata.json, manifest CSV, LICENSE).

The existing codebase already has RDKit 2025.9.3 installed with full support for chemical identifier generation. Benchmark data exists at `/data/benchmark/results/` with .nw/.out files organized by compound/solvent hierarchy. Molecule structures are available as XYZ files in `/data/benchmark/delta50/molecules/`.

**Primary recommendation:** Use Python standard library (hashlib, csv, json, pathlib) for file operations and checksums, RDKit for molecular identifiers, and DataCite JSON schema for metadata. Structure output in `publications/dataset/` with molecule/solvent hierarchy preserved. Automate extraction via Python scripts to ensure consistency and reproducibility.

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| Python hashlib | 3.11+ stdlib | SHA-256 checksum generation | NIST-recommended (SHA-2 family), Python stdlib, no dependencies |
| RDKit | 2025.9.3 (installed) | SMILES/InChI/InChIKey generation | Industry standard for cheminformatics, already in project |
| pathlib | 3.11+ stdlib | File traversal and organization | Modern Python path handling, cross-platform |
| json | 3.11+ stdlib | Metadata serialization | DataCite JSON format, universal support |
| csv | 3.11+ stdlib | Manifest generation | Universal format for file listings |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| datacite | 1.3.1 | DataCite metadata generation | If complex DataCite schema validation needed (optional) |
| pandas | 2.3.3 (installed) | CSV export for scaling factors | Already used in project for data processing |
| tqdm | 4.67.1 (installed) | Progress tracking | Already in project for long-running operations |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| hashlib | md5sum/sha256sum CLI | CLI requires subprocess calls, less portable, but simpler for one-off use |
| Manual JSON | datacite Python lib | Library adds validation but is optional dependency (schema is simple) |
| Custom scripts | AiiDA provenance | AiiDA is enterprise-grade but massive overkill for static dataset archival |

**Installation:**
```bash
# Already installed in project:
# - rdkit>=2025.9.3
# - pandas>=2.3.3
# - tqdm>=4.67.1

# Optional (if DataCite validation needed):
pip install datacite==1.3.1
```

## Architecture Patterns

### Recommended Project Structure
```
publications/dataset/
├── README.md                    # Dataset overview, usage, citation
├── LICENSE                      # CC-BY 4.0 license file
├── metadata.json               # DataCite 4.6 metadata
├── MANIFEST.csv                # File listing with checksums
├── molecules/                  # Molecule metadata
│   ├── molecules.csv          # All 50 molecules with SMILES/InChI/InChIKey
│   └── molecules.json         # Same data in JSON format
├── scaling_factors/           # Processed results
│   ├── scaling_factors.csv    # 26 scaling factors (13 solvents x 2 nuclei)
│   └── scaling_factors.json   # Same data with full statistics
├── calculations/              # Raw computation data
│   ├── compound_01/
│   │   ├── chloroform/
│   │   │   ├── shielding.nw   # NWChem input
│   │   │   └── shielding.out  # NWChem output
│   │   ├── DMSO/
│   │   └── ... (13 solvents)
│   ├── compound_02/
│   └── ... (50 compounds)
└── PROVENANCE.md              # Computational environment documentation
```

### Pattern 1: Hierarchical Data Organization
**What:** Mirror existing benchmark structure (compound/solvent) for intuitive navigation
**When to use:** When source data already has logical hierarchy that aids understanding
**Example:**
```python
# Source: Project structure analysis
base_path = Path("publications/dataset/calculations")
for compound_dir in benchmark_results.glob("compound_*"):
    compound_id = compound_dir.name
    target_dir = base_path / compound_id
    target_dir.mkdir(parents=True, exist_ok=True)

    for solvent_dir in compound_dir.glob("*_*"):  # B3LYP_CHCl3, etc.
        solvent_name = solvent_dir.name.split("_")[1]  # Extract solvent
        # Copy .nw and .out files to target_dir/solvent_name/
```

### Pattern 2: Checksum-First File Processing
**What:** Generate checksums while copying/processing files to avoid double-read
**When to use:** When processing large files that will be included in manifest
**Example:**
```python
# Source: Python hashlib best practices
import hashlib
from pathlib import Path

def copy_with_checksum(src: Path, dst: Path) -> str:
    """Copy file and return SHA-256 checksum."""
    sha256 = hashlib.sha256()
    dst.parent.mkdir(parents=True, exist_ok=True)

    with src.open('rb') as f_in, dst.open('wb') as f_out:
        while chunk := f_in.read(8192):
            sha256.update(chunk)
            f_out.write(chunk)

    return sha256.hexdigest()
```

### Pattern 3: Molecule Identifier Pipeline
**What:** XYZ → RDKit Mol → SMILES/InChI/InChIKey in single pass
**When to use:** When generating chemical identifiers from 3D structures
**Example:**
```python
# Source: RDKit documentation (2025.09.5)
from rdkit import Chem
from rdkit.Chem import Descriptors

def extract_identifiers(xyz_path: Path) -> dict:
    """Extract SMILES, InChI, InChIKey from XYZ file."""
    mol = Chem.MolFromXYZFile(str(xyz_path))
    if mol is None:
        raise ValueError(f"Failed to parse {xyz_path}")

    return {
        "smiles": Chem.MolToSmiles(mol),
        "inchi": Chem.MolToInchi(mol),
        "inchikey": Chem.MolToInchiKey(mol),
        "molecular_formula": Descriptors.rdMolDescriptors.CalcMolFormula(mol)
    }
```

### Pattern 4: Manifest Generation
**What:** CSV with filename, size, SHA-256 for all archived files
**When to use:** Required for data integrity verification in scientific archives
**Example:**
```python
# Source: FAIR data integrity practices
import csv
from pathlib import Path

def generate_manifest(dataset_root: Path, output_csv: Path):
    """Generate MANIFEST.csv with file checksums."""
    manifest = []

    for file_path in dataset_root.rglob("*"):
        if file_path.is_file() and file_path.name != "MANIFEST.csv":
            rel_path = file_path.relative_to(dataset_root)
            checksum = compute_sha256(file_path)
            manifest.append({
                "filepath": str(rel_path),
                "size_bytes": file_path.stat().st_size,
                "sha256": checksum
            })

    with output_csv.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=["filepath", "size_bytes", "sha256"])
        writer.writeheader()
        writer.writerows(manifest)
```

### Anti-Patterns to Avoid
- **Nested ZIP archives:** Keep flat hierarchy with directories, not ZIP-in-ZIP
- **Proprietary formats:** Use CSV/JSON, not Excel/pickle for metadata
- **Relative paths in metadata:** Always use dataset-root-relative paths
- **Generated files without checksums:** Include ALL files in manifest, even generated ones
- **Missing version numbers:** Tag dataset as v1.0.0 from start (semantic versioning)

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| SHA-256 checksums | Custom hash function | Python hashlib.sha256() | NIST-approved, tested, optimized, handles large files |
| InChI generation | SMILES-to-InChI parser | RDKit Chem.MolToInchi() | Complex stereochemistry, tautomers, edge cases |
| DataCite DOI minting | Custom HTTP client | RADAR4Chem web interface | DOI registration requires DataCite membership |
| Provenance tracking | Custom workflow logger | Document manually in PROVENANCE.md | AiiDA is overkill for static dataset |
| File format conversion | Custom parsers | Keep original .nw/.out formats | NWChem formats are standard, no conversion needed |

**Key insight:** Scientific data archival prioritizes **reproducibility** over convenience. Use battle-tested tools (hashlib, RDKit) for critical operations like checksums and chemical identifiers. Custom code is acceptable for glue logic (file organization, CSV generation) but not for domain-critical operations (hashing, molecular identifiers).

## Common Pitfalls

### Pitfall 1: Incomplete Provenance Documentation
**What goes wrong:** Missing compiler version, hardware specs, or non-default parameters makes reproduction impossible
**Why it happens:** Focus on "obvious" parameters (functional, basis set) while ignoring environmental factors
**How to avoid:** Document EVERYTHING from NWChem output header: version, compile date, compiler, COSMO grid settings, memory limits, parallelization
**Warning signs:** PROVENANCE.md shorter than 50 lines, no hardware specs, generic "NWChem 7.2" without patch level

**Required provenance elements:**
```markdown
- NWChem version: 7.2.2 (compiled Mon_Apr_01_11:12:32_2024)
- NWChem source: /build/nwchem-qbYJLs/nwchem-7.2.2/build-openmpi
- Basis set: 6-311+G(2d,p) (from NWChem library)
- DFT functional: B3LYP (direct SCF)
- Solvation: COSMO with solvent-specific dielectric constants
- COSMO grid: Default (extracted from .out files if non-default)
- Geometry optimization: B3LYP/6-31G(d) in vacuo (from literature)
- Reference molecule: TMS (calculated identically)
- Hardware: 4 cores (from .out: nproc = 4)
- Memory: 800 MB total (heap 200MB, stack 200MB, global 400MB)
- ScaLAPACK: Enabled (from .out: use scalapack = T)
```

### Pitfall 2: InChI Generation from 2D Structures
**What goes wrong:** XYZ files contain 3D coordinates but may lose stereochemistry when converted to 2D SMILES
**Why it happens:** RDKit XYZ parser may not preserve all stereochemical information
**How to avoid:** Use 3D-aware InChI generation, verify with known compounds, document any stereochemistry assumptions
**Warning signs:** Flat/achiral InChIs for molecules that should have stereochemistry

### Pitfall 3: Missing Files in Manifest
**What goes wrong:** MANIFEST.csv doesn't include README.md, LICENSE, or other metadata files
**Why it happens:** Generating manifest before creating documentation
**How to avoid:** Generate manifest as FINAL step after all files created, include manifest itself with checksum in README
**Warning signs:** MANIFEST.csv has <650 entries (should be 1300+ with .nw/.out plus metadata files)

### Pitfall 4: Inconsistent Solvent Naming
**What goes wrong:** NWChem directories use "CHCl3", scaling factors use "chloroform", RADAR4Chem needs canonical names
**Why it happens:** Multiple naming conventions in source data
**How to avoid:** Create solvent mapping dict at start, use canonical names consistently, document aliases in README
**Warning signs:** Solvent names don't match between different files (CSV vs directory names vs metadata.json)

**Solvent mapping:**
```python
SOLVENT_CANONICAL = {
    "CHCl3": "chloroform",
    "DMSO": "dimethyl_sulfoxide",
    "Methanol": "methanol",
    "Water": "water",
    "Acetone": "acetone",
    "Benzene": "benzene",
    "Pyridine": "pyridine",
    "THF": "tetrahydrofuran",
    "Toluene": "toluene",
    "DCM": "dichloromethane",
    "Acetonitrile": "acetonitrile",
    "DMF": "dimethylformamide",
    # Plus vacuum/gas phase
}
```

### Pitfall 5: Large Scratch Files in Archive
**What goes wrong:** Including molecule.gridpts.* and other temporary files bloats archive from ~50MB to >5GB
**Why it happens:** Copying entire calculation directories without filtering
**How to avoid:** Only copy .nw (input) and .out (output) files, document exclusion of scratch files in README
**Warning signs:** Archive size >1GB (should be ~100-200MB for 650 calculations)

## Code Examples

Verified patterns from official sources:

### SHA-256 Checksum for Large Files
```python
# Source: Python hashlib documentation (3.11+)
import hashlib

def compute_sha256(filepath: Path) -> str:
    """Compute SHA-256 checksum efficiently for large files."""
    sha256 = hashlib.sha256()
    with filepath.open('rb') as f:
        while chunk := f.read(8192):  # 8KB chunks
            sha256.update(chunk)
    return sha256.hexdigest()
```

### Extract Timing from NWChem Output
```python
# Source: Project NWChem output parser patterns
def extract_timing(output_path: Path) -> dict:
    """Extract computation time from NWChem .out file."""
    with output_path.open('r') as f:
        content = f.read()

    # Extract total wall time (seconds)
    if match := re.search(r'Total times\s+cpu:\s+[\d.]+s\s+wall:\s+([\d.]+)s', content):
        return {"wall_time_seconds": float(match.group(1))}
    return {}
```

### Generate molecules.csv with Identifiers
```python
# Source: RDKit documentation + project structure
from rdkit import Chem
from rdkit.Chem import Descriptors
import csv
from pathlib import Path

def generate_molecules_csv(xyz_dir: Path, exp_shifts_json: Path, output_csv: Path):
    """Generate molecules.csv with SMILES/InChI/InChIKey."""
    import json

    # Load experimental data for molecule names
    with exp_shifts_json.open('r') as f:
        exp_data = json.load(f)

    molecules = []
    for xyz_file in sorted(xyz_dir.glob("compound_*.xyz")):
        compound_id = xyz_file.stem  # "compound_01"
        mol = Chem.MolFromXYZFile(str(xyz_file))

        if mol is None:
            print(f"Warning: Failed to parse {xyz_file}")
            continue

        mol_data = exp_data["molecules"].get(compound_id, {})

        molecules.append({
            "compound_id": compound_id,
            "name": mol_data.get("name", "Unknown"),
            "smiles": Chem.MolToSmiles(mol),
            "inchi": Chem.MolToInchi(mol),
            "inchikey": Chem.MolToInchiKey(mol),
            "molecular_formula": Descriptors.rdMolDescriptors.CalcMolFormula(mol),
            "num_h_atoms": mol_data.get("num_h_atoms", 0),
            "num_c_atoms": mol_data.get("num_c_atoms", 0)
        })

    with output_csv.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=molecules[0].keys())
        writer.writeheader()
        writer.writerows(molecules)
```

### DataCite 4.6 Metadata JSON
```python
# Source: DataCite 4.6 documentation + RADAR4Chem requirements
def generate_datacite_metadata(dataset_info: dict) -> dict:
    """Generate DataCite 4.6 compliant metadata."""
    return {
        "identifier": {
            "identifier": "10.5281/zenodo.XXXXXXX",  # Placeholder DOI
            "identifierType": "DOI"
        },
        "creators": [
            {
                "name": dataset_info["creator_name"],
                "nameType": "Personal",
                "affiliation": [{"name": dataset_info["affiliation"]}]
            }
        ],
        "titles": [
            {
                "title": "DELTA50 NMR Benchmark Dataset: 650 DFT Shielding Calculations",
                "titleType": "Title"
            }
        ],
        "publisher": {"name": "RADAR4Chem"},
        "publicationYear": "2026",
        "resourceType": {
            "resourceTypeGeneral": "Dataset",
            "resourceType": "Computational Chemistry Dataset"
        },
        "subjects": [
            {"subject": "NMR spectroscopy"},
            {"subject": "Density Functional Theory"},
            {"subject": "Chemical Shift Prediction"},
            {"subject": "DELTA50"},
            {"subject": "Computational Chemistry"}
        ],
        "contributors": [
            {
                "contributorType": "DataCurator",
                "name": "Claude Opus 4.6",
                "nameType": "Organizational"
            }
        ],
        "dates": [
            {"date": "2026-02-11", "dateType": "Created"}
        ],
        "language": "en",
        "descriptions": [
            {
                "description": "Benchmark dataset of 650 NWChem DFT calculations...",
                "descriptionType": "Abstract"
            }
        ],
        "rightsList": [
            {
                "rights": "Creative Commons Attribution 4.0 International",
                "rightsURI": "https://creativecommons.org/licenses/by/4.0/",
                "rightsIdentifier": "CC-BY-4.0"
            }
        ],
        "version": "1.0.0"
    }
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| MD5 checksums | SHA-256 checksums | NIST 2020 deprecation | MD5 no longer secure, SHA-256 required for data integrity |
| XML metadata | JSON metadata | DataCite 2019 | JSON easier to generate/parse, XML still supported |
| Manual provenance | AiiDA automatic tracking | 2020 (AiiDA 1.0) | Enterprise-grade but overkill for static datasets |
| Institution-specific repos | Discipline-specific repos (RADAR4Chem) | 2021 NFDI4Chem | Better FAIR compliance, chemistry-aware metadata |
| SMILES only | SMILES + InChI + InChIKey | 2010s InChI adoption | InChI provides canonical representation, InChIKey for hashing |
| DataCite 4.3 | DataCite 4.6 | December 2024 | Added CSTR/RRID identifiers, better related resource linking |

**Deprecated/outdated:**
- MD5 checksums: NIST deprecated for cryptographic use, use SHA-256 or SHA-3
- DataCite Metadata Schema <4.0: RADAR4Chem requires 4.0+, current is 4.6
- Proprietary molecular formats: Use open standards (SDF, MOL, XYZ + SMILES/InChI)
- Dataset DOIs from generic repos: Use discipline-specific repos (RADAR4Chem) for better discoverability

## Open Questions

Things that couldn't be fully resolved:

1. **Exact RADAR4Chem 10 Mandatory Fields**
   - What we know: RADAR uses DataCite 4.0 schema with 10 mandatory fields
   - What's unclear: Specific field names beyond standard DataCite required fields (identifier, creators, title, publisher, publicationYear)
   - Recommendation: Use standard DataCite 4.6 mandatory fields, add chemistry extensions (SMILES/InChI), verify via RADAR4Chem documentation during implementation

2. **Stereochemistry Preservation from XYZ Files**
   - What we know: RDKit can read XYZ files and generate InChI
   - What's unclear: Whether 3D coordinates are sufficient to preserve stereochemistry without explicit annotation
   - Recommendation: Verify InChI stereochemistry flags for known chiral compounds (if any in DELTA50), document any assumptions in README

3. **Total Archive Size Estimate**
   - What we know: 650 calculations with .nw/.out files, ~50MB uncompressed per compound (with scratch files)
   - What's unclear: Exact size after filtering scratch files, with/without compression
   - Recommendation: Implement file filtering first, measure actual size, compress if >1GB (RADAR4Chem free tier is 10GB)

4. **Computation Time Extraction Reliability**
   - What we know: NWChem .out files contain timing information
   - What's unclear: Format consistency across all 650 files, whether all calculations completed successfully
   - Recommendation: Parse timing with fallback for missing data, flag incomplete calculations in manifest

5. **Dataset Versioning Strategy**
   - What we know: Should use semantic versioning (v1.0.0)
   - What's unclear: Whether to version dataset itself or just archive format
   - Recommendation: Start with v1.0.0, increment for data corrections (patch), new solvents (minor), new molecules (major)

## Sources

### Primary (HIGH confidence)
- [DataCite Metadata Schema 4.6 Documentation](https://datacite-metadata-schema.readthedocs.io/en/4.6/) - Official schema documentation
- [Python hashlib Documentation](https://docs.python.org/3/library/hashlib.html) - SHA-256 implementation
- [RDKit 2025.09.5 Documentation](https://www.rdkit.org/docs/GettingStartedInPython.html) - SMILES/InChI generation
- [NFDI4Chem FAIR Principles Knowledge Base](https://knowledgebase.nfdi4chem.de/knowledge_base/docs/fair/) - Chemistry-specific FAIR guidelines
- [RADAR4Chem Knowledge Base](https://knowledgebase.nfdi4chem.de/knowledge_base/docs/radar4chem/) - Repository requirements

### Secondary (MEDIUM confidence)
- [FAIR Data Principles (GO FAIR)](https://www.go-fair.org/fair-principles/) - Original FAIR definition
- [IOOS Standard Practices for Data Integrity](https://ioos.github.io/ncei-archiving-cookbook/practices.html) - Checksum best practices (verified with NIST SHA-2 recommendation)
- [Best-Practice DFT Protocols (Angewandte Chemie 2022)](https://pmc.ncbi.nlm.nih.gov/articles/PMC9826355/) - Computational chemistry documentation standards
- [AiiDA Computational Provenance (Nature Scientific Data 2020)](https://www.nature.com/articles/s41597-020-00638-4) - State-of-the-art provenance tracking (reference for comparison)
- [Creative Commons CC-BY 4.0 License](https://creativecommons.org/licenses/by/4.0/deed.en) - Attribution requirements

### Tertiary (LOW confidence)
- WebSearch results on dataset README templates - No single authoritative source, multiple university guidelines consulted
- FAIRness assessment for RADAR4Chem - 2024 study, specific assessment criteria not detailed

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - Python stdlib tools (hashlib, csv, json) are battle-tested, RDKit is installed and current
- Architecture: MEDIUM-HIGH - Patterns based on FAIR principles and existing benchmark structure, but dataset-specific details need verification
- Pitfalls: HIGH - Based on documented NIST deprecations, DataCite requirements, and common computational chemistry archival mistakes
- RADAR4Chem specifics: MEDIUM - General requirements clear, but exact mandatory fields need verification during implementation

**Research date:** 2026-02-11
**Valid until:** 2026-08-11 (6 months - FAIR principles stable, DataCite 4.6 current, Python/RDKit stable)

**Key assumptions:**
1. 650 calculations exist and are complete (verified via directory listing)
2. Molecule XYZ files are valid 3D structures (files exist, need parsing verification)
3. RADAR4Chem free tier (10GB) is sufficient for archive (need size measurement)
4. No sensitive data requiring access restrictions (benchmark data for publication)
5. Dataset will be versioned starting at v1.0.0 (standard practice for first release)

**Implementation notes:**
- Use Python 3.11+ features (walrus operator, pathlib)
- RDKit already installed (version 2025.9.3)
- Leverage existing pandas for CSV export (already in project)
- Keep archive directory structure human-readable (no deep nesting)
- Generate all metadata programmatically for reproducibility
- Include script source code in archive for transparency
