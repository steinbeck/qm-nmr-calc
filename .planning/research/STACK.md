# Technology Stack: DELTA50 Dataset Publication and Papers

**Project:** qm-nmr-calc v3.0 - Dataset and Publication Release
**Milestone:** Publishing DELTA50 benchmark dataset and two papers
**Researched:** 2026-02-11
**Confidence:** HIGH

## Executive Summary

This milestone requires tools for three distinct deliverables: dataset packaging for Radar4Chem repository, a data descriptor paper for Nature Scientific Data, and an application note for Journal of Cheminformatics. The stack leverages existing Python/RDKit infrastructure while adding metadata standards (DataCite), data packaging (Frictionless Data), LaTeX tooling (Springer Nature template), and visualization enhancements (matplotlib/seaborn for publication figures).

**Key principle:** Build on existing validated stack. Do NOT re-implement NMR calculation pipeline - only add publication and packaging layers.

---

## Core Publication Tools

### Dataset Metadata Standards
| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| DataCite Metadata Schema | 4.6 (Dec 2024) | Repository metadata for Radar4Chem | NFDI4Chem standard, based on DataCite 4.0 with 10 mandatory + 13 optional fields. Latest version released Dec 2024. |
| JCAMP-DX | 5.0+ | NMR spectroscopy data format | IUPAC standard for NMR data exchange, human-readable metadata, accepted by nmrXiv and NFDI4Chem repositories. |
| JSON | Native | Structured metadata export | DataCite XML to JSON mapping supported by DataCite REST API, easier to generate programmatically from Python. |

**Rationale:** Radar4Chem uses RADAR Metadata Schema based on DataCite 4.0. Current version is 4.6 (released Dec 2024). JSON format preferred over XML for Python generation. JCAMP-DX is IUPAC standard for NMR data - optional but adds domain-specific value.

### Data Packaging
| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| Frictionless Data Package | 5.18.1+ | Dataset packaging standard | FAIR-compliant packaging format, Python library available (frictionless-py), used by Zenodo and data repositories. |
| pandas | 2.3.3+ (existing) | Tabular data processing | Already in stack, excellent CSV/Parquet export for benchmark tables. |
| Parquet | N/A (format) | Archival data format | Columnar format preserves types, better than CSV for computational data, smaller files, PyArrow integration. |

**Rationale:** Frictionless Data Package provides FAIR-compliant packaging with datapackage.json descriptor. Parquet preserves data types (unlike CSV) and is widely supported in computational chemistry. Both integrate well with existing pandas workflow.

### LaTeX and Document Preparation
| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| Springer Nature LaTeX Template | Dec 2024 | Both journal submissions | Official template for all Springer Nature journals (Nature Scientific Data, Journal of Cheminformatics). Available on Overleaf and direct download. |
| Overleaf | Online | Collaborative LaTeX editing | Official integration with Springer Nature, simplifies submission, handles compilation with TexLive 2021. |
| BibTeX | Native | References management | Required by both journals, embed directly in .TEX file (do NOT use separate .bib/.bbl for Scientific Data). |

**Rationale:** Both target journals are Springer Nature publications. Use the unified Springer Nature LaTeX template (Dec 2024 version). Nature Scientific Data explicitly states "do not use legacy templates" - use official Dec 2024 version. Journal of Cheminformatics compiles with pdfLaTeX and TexLive 2021. Overleaf provides official Springer Nature integration.

**IMPORTANT:** Nature Scientific Data does NOT recommend templates and prefers simple heading structure. If using LaTeX, embed references in .TEX file, NOT separate .bib files.

---

## Data Visualization and Figures

### Publication-Quality Figures
| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| matplotlib | 3.10.0+ (existing) | Scientific plotting | Already in stack, produces publication-quality figures, 20+ years development, industry standard. |
| seaborn | Latest (add) | Statistical visualizations | Built on matplotlib, ideal for correlation heatmaps and statistical distributions in NMR benchmark data. |
| RDKit | 2025.9.3+ (existing) | Molecular structure images | Already in stack, generates publication-quality 2D structures, integrates with matplotlib. |

**Rationale:** Leverage existing matplotlib/RDKit stack. Add seaborn for statistical visualizations (correlation matrices for scaling factors, distribution plots for chemical shift errors). All three integrate seamlessly.

### Figure Requirements (Both Journals)
| Requirement | Specification | Tools |
|-------------|---------------|-------|
| Resolution | 300 DPI minimum (color/photos), 600 DPI (line art) | matplotlib savefig(dpi=600) |
| Format | TIFF or layered Photoshop (preferred), PNG/PDF acceptable | matplotlib supports TIFF, PNG, PDF |
| Width | 80-180mm | matplotlib figsize parameter |
| Maximum figures | 8 for Data Descriptors (Scientific Data) | Plan figure count carefully |

**Implementation:**
```python
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(6.3, 4))  # 6.3 inches = 160mm
# ... plotting code ...
fig.savefig('figure1.tiff', dpi=600, format='tiff', bbox_inches='tight')
```

---

## Data Export and Processing

### Existing Stack (DO NOT CHANGE)
| Technology | Version | Purpose | Notes |
|------------|---------|---------|-------|
| FastAPI | 0.128.0+ | Web service (not used for publication) | Existing infrastructure |
| RDKit | 2025.9.3+ | Molecular structure handling | Use for SMILES, structure images |
| pandas | 2.3.3+ | Data processing | Use for benchmark tables, CSV export |
| scipy | 1.17.0+ | Statistical analysis | Use for regression statistics |
| statsmodels | 0.14.0+ | OLS regression | Use for scaling factor derivation |

### New Additions for Publication
| Library | Version | Purpose | Installation |
|---------|---------|---------|-------------|
| frictionless | 5.18.1+ | Data package creation | `uv add frictionless` |
| seaborn | Latest | Statistical plots | `uv add seaborn` |
| pyarrow | Latest | Parquet file support | `uv add pyarrow` (pandas integration) |

**Rationale:** Minimal additions to existing stack. Frictionless for FAIR packaging, seaborn for publication plots, pyarrow for Parquet export (pandas integration).

---

## Repository Submission Requirements

### Radar4Chem (NFDI4Chem Repository)
| Requirement | Specification | Implementation |
|-------------|---------------|----------------|
| Submission method | Manual upload via web interface | No API available, user leads NFDI4Chem so has access |
| Metadata schema | RADAR Metadata Schema v9.2 (based on DataCite 4.0) | Generate JSON, fill web form manually |
| Mandatory fields | 10 required (title, creator, publication year, resource type, etc.) | Script to extract from existing data |
| Optional fields | 13 recommended (description, keywords, funding, etc.) | Include all available metadata |
| Identifiers | ORCID for authors, CrossRef Funder Registry for funding | User has ORCID, specify grant info |
| Storage limit | 10 GB per project (free tier) | DELTA50 dataset << 10 GB (650 calculations) |
| Preservation | Minimum 25 years | Automatic with DOI assignment |
| License | User-defined (CC-BY recommended) | Choose at upload time |
| DOI | Assigned by DataCite on publication | Cite in papers |

**Workflow:**
1. Export data to Frictionless Data Package (datapackage.json + files)
2. Generate DataCite 4.6 metadata JSON
3. Upload to Radar4Chem web interface
4. Fill mandatory metadata fields
5. Add optional fields (description, keywords, funding)
6. Choose CC-BY-4.0 license
7. Publish and receive DOI
8. Use DOI in Data Availability statements in papers

### Data Package Structure
```
DELTA50-benchmark/
├── datapackage.json          # Frictionless metadata
├── datacite.json             # DataCite 4.6 metadata
├── README.md                 # Human-readable description
├── LICENSE                   # CC-BY-4.0
├── data/
│   ├── benchmark_results.parquet  # Processed results
│   ├── scaling_factors.csv   # 26 OLS parameter sets
│   └── molecule_metadata.csv # DELTA50 molecule info
├── raw_calculations/
│   ├── molecule_001/
│   │   ├── acetone/
│   │   │   ├── input.nw      # NWChem input
│   │   │   └── output.out    # NWChem output
│   │   └── ... (13 solvents)
│   └── ... (50 molecules)
└── scripts/
    ├── extract_shielding.py  # Tensor extraction
    └── derive_scaling.py     # OLS regression
```

---

## Journal-Specific Requirements

### Nature Scientific Data - Data Descriptor
| Requirement | Specification | Tools |
|-------------|---------------|-------|
| Article type | Data Descriptor | Max 8 figures, Methods + Data Records + Technical Validation sections |
| Format | Word or LaTeX | LaTeX preferred, use Springer Nature template (Dec 2024) |
| References | Embedded in .TEX file | Do NOT use separate .bib/.bbl files |
| Figure limit | Maximum 8 figures | Choose most important visualizations |
| Data availability | Public repository with DOI | Radar4Chem DOI (step 1 before submission) |
| Submission | ScholarOne Manuscripts | Online submission system |
| Peer review | Standard peer review | Optional for data descriptors per NFDI4Chem |
| APC | Check current rates | Open access journal |

**Sections required:**
1. **Background & Summary** - DELTA50 benchmark context, NMR prediction validation needs
2. **Methods** - NWChem DFT parameters, COSMO solvation, OLS scaling derivation
3. **Data Records** - Repository location (Radar4Chem DOI), file formats, folder structure
4. **Technical Validation** - R² values, comparison to literature, coverage of chemical space
5. **Usage Notes** - How to use scaling factors, software requirements
6. **Code Availability** - Scripts for tensor extraction and scaling derivation

### Journal of Cheminformatics - Application Note
| Requirement | Specification | Tools |
|-------------|---------------|-------|
| Article type | Software/Application Note | Focus on tool utility, less emphasis on validation |
| Format | LaTeX (required) | Springer Nature template, compiled with pdfLaTeX/TexLive 2021 |
| Submission | BMC/Springer Nature platform | Editable source files required (.TEX + figures) |
| Software paper | Source code + documentation | GitHub repository (existing) + dataset DOI |
| Figures | Publication quality | 300+ DPI, TIFF/PNG/PDF |
| APC | £1690 GBP / $2390 USD / €1990 EUR | Open access fee (check funding/waivers) |
| Peer review | Standard | Expect 2-3 reviewers, 4-8 week turnaround |

**Sections recommended:**
1. **Abstract** - Tool purpose, dataset, availability
2. **Background** - NMR prediction landscape, DELTA50 gap
3. **Implementation** - FastAPI/RDKit/NWChem stack (existing milestone deliverables)
4. **Results and Discussion** - Benchmark performance, scaling factor accuracy
5. **Conclusions** - Tool utility, dataset availability
6. **Availability and Requirements** - GitHub + Radar4Chem DOI, Python 3.11+, Docker

**Focus:** This paper describes the qm-nmr-calc TOOL using the DELTA50 dataset as validation, whereas Scientific Data paper describes the DATASET itself.

---

## File Format Standards

### Primary Formats
| Format | Purpose | Rationale |
|--------|---------|-----------|
| Parquet | Benchmark results table | Preserves types, efficient, widely supported |
| CSV | Scaling factors, molecule metadata | Simple, human-readable, Excel-compatible |
| .nw | NWChem input files | Native format, archival (already have 650) |
| .out | NWChem output files | Native format, archival (already have 650) |
| JSON | Processed shielding tensors | Structured, Python-native, easy parsing |
| JCAMP-DX | NMR spectra (optional) | IUPAC standard, adds domain value |

### FAIR Compliance
| Principle | Implementation |
|-----------|----------------|
| Findable | DOI from Radar4Chem, DataCite metadata, keywords |
| Accessible | Open repository, CC-BY license, no login required |
| Interoperable | Standard formats (Parquet, CSV, JCAMP-DX), clear schema |
| Reusable | MIT license for code, CC-BY for data, documentation |

---

## Metadata Standards Summary

### DataCite 4.6 Mandatory Fields (10)
1. Identifier (DOI from Radar4Chem)
2. Creator (authors with ORCID)
3. Title ("DELTA50 Benchmark Dataset for NMR Quantum Chemical Predictions")
4. Publisher (Radar4Chem / NFDI4Chem)
5. Publication Year (2026)
6. Resource Type (Dataset)
7. Subject (Chemistry, Computational Chemistry, NMR Spectroscopy)
8. Contributor (if applicable)
9. Date (creation, publication)
10. Language (English)

### DataCite 4.6 Recommended Optional Fields (13)
- Description (abstract of dataset)
- Funding Reference (grant information)
- Related Identifier (papers, code repository)
- Rights (CC-BY-4.0)
- Version (1.0.0)
- GeoLocation (not applicable)
- Format (Parquet, CSV, NWChem)
- Size (file sizes, number of calculations)

### MIChI Compliance (NFDI4Chem Metadata Standard)
MIChI (Minimum Information for Chemical Investigations) is NFDI4Chem's developing metadata standard for chemistry datasets. Include:
- Method type (DFT quantum chemistry)
- Instrument/software (NWChem 7.0+)
- Parameters (B3LYP/6-311+G(2d,p), COSMO solvation)
- Sample information (50 molecules, 13 solvents)
- Measured properties (isotropic shielding tensors, 1H and 13C)

---

## Alternatives Considered

| Category | Recommended | Alternative | Why Not |
|----------|-------------|-------------|---------|
| Repository | Radar4Chem | Zenodo, Figshare | User leads NFDI4Chem, domain-specific repository preferred |
| Data packaging | Frictionless Data | Custom JSON | Frictionless is FAIR standard, widely supported |
| Archival format | Parquet | HDF5 | Parquet better for tabular data, HDF5 for multidimensional arrays |
| LaTeX | Springer Nature template | Custom template | Both journals are Springer Nature, use official template |
| Metadata | DataCite 4.6 | Dublin Core | DataCite is requirement for Radar4Chem |
| NMR format | JCAMP-DX | nmrML, NMReDATA | JCAMP-DX is established IUPAC standard, simpler for DFT data |
| Figures | matplotlib + seaborn | Plotly | matplotlib is publication standard, Plotly is interactive/web |

### Why NOT nmrML or NMReDATA?
- **nmrML:** XML-based, complex schema designed for raw instrument data. DELTA50 has computed data (DFT), not experimental spectra.
- **NMReDATA:** Designed for experimental NMR with chemical shift assignments and 2D correlations. Our dataset has isotropic shielding tensors (computational), not assigned spectra.
- **JCAMP-DX:** Simpler, text-based, IUPAC standard, adequate for tabular NMR data. Can include computed shielding as data blocks.

**Decision:** Use JCAMP-DX if including spectral-like data, otherwise stick with Parquet/CSV for computational results.

---

## Installation and Setup

### New Dependencies
```bash
# Add publication tools to existing environment
uv add frictionless seaborn pyarrow
```

### LaTeX Setup
**Option 1: Overleaf (Recommended)**
1. Go to https://www.overleaf.com/latex/templates/springer-nature-latex-template/
2. Open Springer Nature LaTeX Template
3. Start writing directly in browser
4. Submission: Download source files for journal submission

**Option 2: Local LaTeX**
```bash
# Install TexLive 2021+ (required by Journal of Cheminformatics)
sudo apt-get install texlive-full  # Linux
# or download from https://www.tug.org/texlive/

# Download template
wget https://www.springernature.com/gp/authors/campaigns/latex-author-support
# Unzip and follow template README
```

### Verification
```bash
# Test frictionless
python -c "import frictionless; print(frictionless.__version__)"

# Test seaborn
python -c "import seaborn; print(seaborn.__version__)"

# Test parquet
python -c "import pandas as pd; import pyarrow; df = pd.DataFrame({'a': [1,2]}); df.to_parquet('test.parquet'); print('OK')"
```

---

## Implementation Workflow

### Phase 1: Dataset Packaging (Week 1-2)
1. Export processed benchmark results to Parquet
2. Export scaling factors to CSV
3. Organize raw NWChem files (already have 650)
4. Create Frictionless datapackage.json
5. Generate DataCite 4.6 metadata JSON
6. Write README.md with dataset description
7. Test local package structure

### Phase 2: Repository Submission (Week 2)
1. Upload to Radar4Chem web interface
2. Fill mandatory metadata fields (10 required)
3. Add optional fields (description, keywords, funding)
4. Choose CC-BY-4.0 license
5. Peer review (optional, user discretion)
6. Publish and receive DOI
7. Verify DOI resolves correctly

### Phase 3: Scientific Data Paper (Week 3-4)
1. Set up Overleaf with Springer Nature template
2. Write Background & Summary
3. Document Methods (NWChem parameters, COSMO, OLS)
4. Describe Data Records (repository structure, formats)
5. Present Technical Validation (R² plots, error distributions)
6. Create figures with matplotlib/seaborn (max 8)
7. Add Code Availability section
8. Submit to Scientific Data

### Phase 4: Journal of Cheminformatics Paper (Week 5-6)
1. Adapt Springer Nature template for JChem
2. Write application note focusing on TOOL utility
3. Reference DELTA50 dataset DOI
4. Describe qm-nmr-calc implementation
5. Show benchmark results and validation
6. Create publication figures
7. Submit to Journal of Cheminformatics

---

## Quality Checklist

Before submission, verify:

**Dataset (Radar4Chem):**
- [ ] All 650 calculations included (.nw input + .out output)
- [ ] Processed results in Parquet format
- [ ] Scaling factors CSV with all 26 parameter sets
- [ ] datapackage.json validates with Frictionless
- [ ] datacite.json includes all 10 mandatory fields
- [ ] README.md explains dataset structure
- [ ] LICENSE file included (CC-BY-4.0)
- [ ] DOI obtained from Radar4Chem

**Scientific Data Paper:**
- [ ] Uses Springer Nature LaTeX template (Dec 2024)
- [ ] References embedded in .TEX file (no separate .bib)
- [ ] Maximum 8 figures
- [ ] All figures 300+ DPI
- [ ] Data Records section cites Radar4Chem DOI
- [ ] Technical Validation shows R² > 0.99 for all sets
- [ ] Code Availability section complete

**Journal of Cheminformatics Paper:**
- [ ] Uses Springer Nature LaTeX template
- [ ] Editable source files (.TEX + separate figures)
- [ ] Figures 300+ DPI (color), 600+ DPI (line art)
- [ ] Software availability section (GitHub + Docker)
- [ ] Dataset DOI cited from Radar4Chem
- [ ] APC funding confirmed or waiver requested

---

## Cost and Funding Considerations

| Item | Cost | Notes |
|------|------|-------|
| Radar4Chem | Free | 10 GB limit (sufficient for DELTA50) |
| Nature Scientific Data APC | Check current | Open access journal, APC required |
| Journal of Cheminformatics APC | £1690 / $2390 / €1990 | Open access, request waiver if needed |
| LaTeX (Overleaf) | Free | Free tier sufficient for 2 papers |
| Python packages | Free | All open source |

**Total estimated:** $2390-4000 USD (both APCs), may be lower with institutional agreements or waivers.

---

## Timeline Estimates

| Task | Duration | Dependencies |
|------|----------|--------------|
| Dataset packaging | 1-2 weeks | Existing calculations complete |
| Radar4Chem upload | 1-2 days | Package ready |
| DOI assignment | 1-3 days | Upload complete |
| Scientific Data draft | 2-3 weeks | DOI obtained |
| Scientific Data review | 4-8 weeks | Submission complete |
| JChem draft | 2-3 weeks | Can overlap with ScData review |
| JChem review | 4-8 weeks | Submission complete |

**Critical path:** Radar4Chem DOI must be obtained BEFORE submitting either paper (required for Data Availability statements).

---

## Sources and References

**NFDI4Chem and Radar4Chem:**
- [NFDI4Chem Knowledge Base - FAIR Data Principles](https://knowledgebase.nfdi4chem.de/knowledge_base/docs/fair/)
- [NFDI4Chem Knowledge Base - RADAR4Chem](https://knowledgebase.nfdi4chem.de/knowledge_base/docs/radar4chem/)
- [NFDI4Chem Knowledge Base - Metadata / MIChI](https://knowledgebase.nfdi4chem.de/knowledge_base/docs/metadata/)
- [RADAR4Chem About Page](https://radar.products.fiz-karlsruhe.de/en/radarabout/radar4chem)

**DataCite Metadata Schema:**
- [DataCite Metadata Schema 4.6 Documentation](https://datacite-metadata-schema.readthedocs.io/en/4.6/)
- [DataCite Schema Versions](https://schema.datacite.org/versions.html)

**Nature Scientific Data:**
- [Nature Scientific Data Submission Guidelines](https://www.nature.com/sdata/submission-guidelines)
- [Nature Scientific Data Data Policies](https://www.nature.com/sdata/policies/data-policies)

**Journal of Cheminformatics:**
- [Journal of Cheminformatics Submission Guidelines](https://jcheminf.biomedcentral.com/submission-guidelines)
- [Journal of Cheminformatics Fees and Funding](https://jcheminf.biomedcentral.com/submission-guidelines/fees-and-funding)

**LaTeX Templates:**
- [Springer Nature LaTeX Author Support](https://www.springernature.com/gp/authors/campaigns/latex-author-support)
- [Springer Nature LaTeX Template on Overleaf](https://www.overleaf.com/latex/templates/springer-nature-latex-template/myxmhdsbzkyd)

**Data Formats and Standards:**
- [JCAMP-DX IUPAC Standard](https://iupac.org/what-we-do/digital-standards/jcamp-dx/)
- [nmrML Publication](https://pubs.acs.org/doi/10.1021/acs.analchem.7b02795)
- [NMReDATA Initiative](https://www.nmredata.org/)
- [Frictionless Data Framework](https://framework.frictionlessdata.io/)

**Visualization and Figures:**
- [RDKit Getting Started](https://www.rdkit.org/docs/GettingStartedInPython.html)
- [Matplotlib Documentation](https://matplotlib.org/)
- [Seaborn Documentation](https://seaborn.pydata.org/)

**Python Tools:**
- [ChemML Documentation](https://hachmannlab.github.io/chemml/)
- [Awesome Python Chemistry](https://github.com/lmmentel/awesome-python-chemistry)

**FAIR Principles:**
- [GO FAIR - FAIR Principles](https://www.go-fair.org/fair-principles/)
- [The FAIR Guiding Principles - Scientific Data](https://www.nature.com/articles/sdata201618)

---

## Confidence Assessment

| Area | Confidence | Verification |
|------|------------|--------------|
| Radar4Chem requirements | HIGH | Official NFDI4Chem documentation |
| DataCite metadata | HIGH | Official DataCite 4.6 schema (Dec 2024) |
| LaTeX templates | HIGH | Official Springer Nature template (Dec 2024) |
| Journal requirements | HIGH | Official submission guidelines (both journals) |
| Python tools | HIGH | Existing pyproject.toml + official package docs |
| JCAMP-DX format | MEDIUM | IUPAC standard but unsure of applicability to DFT data |
| Figure requirements | HIGH | Multiple journal guideline sources |
| Timeline estimates | MEDIUM | Based on typical journal turnaround, not specific to these journals |

**Overall confidence: HIGH** - All critical technologies verified with official documentation, versions current as of Feb 2026.
