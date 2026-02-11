# Feature Landscape: Publishing DELTA50 Benchmark Dataset and NMR Tool

**Domain:** Computational chemistry dataset publication and scientific software publication
**Researched:** 2026-02-11
**Confidence:** MEDIUM (based on verified publication guidelines and multiple examples)

## Overview

This research examines expected components for three distinct deliverables:
1. **Institutional repository dataset** (RADAR4Chem/NFDI4Chem)
2. **Nature Scientific Data data descriptor** (dataset publication)
3. **Journal of Cheminformatics software article** (web tool publication)

Each has distinct requirements but they are interdependent and should reference each other.

---

## Deliverable 1: Institutional Repository Dataset (RADAR4Chem)

### Table Stakes

Features users/reviewers expect. Missing these makes the dataset feel incomplete or unprofessional.

| Feature | Why Expected | Complexity | Notes |
|---------|--------------|------------|-------|
| **Persistent DOI** | Required by FAIR principles, enables citation | Low | Auto-assigned by RADAR4Chem |
| **README.md** | Entry point for understanding dataset | Low | Plain text documentation |
| **Molecule structures** | Core data for benchmark | Low | XYZ format (already generated) |
| **Raw calculation files** | Required for reproducibility | Low | .nw (input) and .out (output) files |
| **Processed data (JSON/CSV)** | Machine-readable extracted data | Medium | Shielding tensors, energies |
| **Scaling factors** | Core research output | Low | Already in scaling_factors.json |
| **Experimental reference data** | Ground truth for benchmarking | Medium | 1H and 13C experimental shifts |
| **Metadata (author ORCID)** | Required by RADAR4Chem | Low | Already have ORCID |
| **License declaration** | Required for reuse | Low | Suggest CC-BY 4.0 |
| **Version number** | Dataset versioning | Low | Start with v1.0.0 |
| **File organization** | Logical directory structure | Medium | See recommendations below |

**Recommended file organization:**
```
DELTA50-NMR-Benchmark/
├── README.md                       # Overview, citation, usage
├── LICENSE.txt                     # CC-BY 4.0 or similar
├── METADATA.json                   # Machine-readable metadata
├── molecules/                      # Molecule structures
│   ├── structures.sdf              # All 50 molecules in SDF
│   └── structures.xyz              # All 50 in XYZ format
├── calculations/                   # Raw NWChem files
│   ├── molecule_001/
│   │   ├── chloroform.nw
│   │   ├── chloroform.out
│   │   ├── dmso.nw
│   │   └── ...                     # All 13 solvents
│   └── ...                         # All 50 molecules
├── processed_data/                 # Extracted results
│   ├── shielding_tensors.json      # All calculated shieldings
│   ├── shielding_tensors.csv       # Same in CSV format
│   ├── experimental_shifts.csv     # Experimental 1H/13C shifts
│   └── molecular_properties.csv    # MW, formula, SMILES, InChI
├── scaling_factors/                # Regression results
│   ├── scaling_factors.json        # All 26 factor sets
│   ├── scaling_factors.csv         # Tabular format
│   └── validation_metrics.csv      # R², MAE, RMSD per set
└── supplementary/                  # Additional materials
    ├── computational_methods.md    # Detailed calculation protocol
    ├── solvent_parameters.csv      # COSMO parameters used
    └── references.bib              # Bibliography
```

### Differentiators

Features that make this dataset stand out and encourage reuse.

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| **Multiple solvents (13)** | Most benchmarks use gas phase or 1-2 solvents | Low | Already completed |
| **High-quality experimental data** | DELTA50 carefully curated for accuracy | Low | Already validated |
| **Multiple factor sets (26)** | Covers 1H/13C × 13 solvents | Low | Comprehensive coverage |
| **Quality metrics** | R² > 0.99, MAE < 0.15 ppm (1H) | Medium | Statistical validation |
| **InChI/SMILES identifiers** | Links to other chemistry databases | Medium | Enables cross-referencing |
| **Computation time logs** | Helps users estimate resources | Low | Extract from .out files |
| **Provenance tracking** | Software versions, calculation dates | Medium | Reproducibility metadata |
| **Example use cases** | Jupyter notebooks showing usage | High | Tutorial for users |
| **Machine-readable all formats** | JSON + CSV for all data | Medium | Broad tool compatibility |
| **Checksum/validation** | MD5/SHA256 for file integrity | Low | Data integrity verification |

### Anti-Features

Features to explicitly NOT include. Common mistakes in this domain.

| Anti-Feature | Why Avoid | What to Do Instead |
|--------------|-----------|-------------------|
| **Proprietary formats only** | Excludes users without specific software | Always provide open formats (CSV, JSON, XYZ) |
| **Binary-only output** | Not archival, format may become obsolete | Include text-based outputs (.out files are text) |
| **Undocumented acronyms** | Creates barriers to reuse | Define all terms in README |
| **Mixed versions/inconsistent calculations** | Undermines benchmark validity | Document single NWChem version used consistently |
| **Raw spectra images only** | Not machine-readable | Provide numerical data in addition |
| **Excessive file compression** | 10GB limit not a concern for this dataset | Use simple .zip if needed, not multi-part archives |
| **Custom/exotic file formats** | Reduces accessibility | Stick to standard chemistry formats |
| **Incomplete molecule identifiers** | Hard to cross-reference | Always include SMILES, InChI, name |

---

## Deliverable 2: Nature Scientific Data - Data Descriptor

### Table Stakes

Required sections and components per Scientific Data guidelines.

| Feature | Why Expected | Complexity | Notes |
|---------|--------------|------------|-------|
| **Background & Summary** | Motivation and reuse value | Medium | 2-3 paragraphs, no new findings |
| **Methods section** | Detailed calculation protocol | Medium | DFT level, basis set, COSMO params |
| **Data Records section** | What/where data is hosted | Medium | Describe RADAR4Chem deposit |
| **Technical Validation** | Quality checks and statistics | High | Critical for benchmark credibility |
| **Usage Notes** | How to use the data | Medium | File formats, software requirements |
| **Code Availability** | Links to analysis scripts | Medium | GitHub repo or Zenodo |
| **Data Availability** | Repository DOI | Low | Link to RADAR4Chem DOI |
| **Abstract (< 170 words)** | Concise dataset summary | Low | No new scientific findings |
| **Figures (max 8)** | Visual summaries | High | See figure list below |
| **Tables (no limit but concise)** | Key data summaries | Medium | See table list below |
| **References** | Prior work, methods | Low | Cite NWChem, DELTA50 paper, etc. |
| **ORCID for all authors** | Author identification | Low | Required by Nature |
| **Competing interests statement** | COI disclosure | Low | Standard boilerplate |
| **Author contributions (CRediT)** | Who did what | Low | CRediT taxonomy |

**Recommended figures (5-6 total):**
1. **Molecule diversity chart** - Chemical space coverage (functional groups, ring sizes)
2. **Calculation workflow diagram** - NWChem → processing → scaling factors
3. **Accuracy heatmap** - MAE by solvent and nucleus (1H vs 13C)
4. **Correlation plots** - Experimental vs predicted (2-4 representative solvents)
5. **Statistical validation** - R² distribution across all 26 factor sets
6. **Solvent comparison** - How accuracy varies by solvent (violin plot or box plot)

**Recommended tables (3-4 total):**
1. **Dataset summary** - 50 molecules × 13 solvents = 650 calculations, completion rate
2. **Computational details** - NWChem version, DFT functional, basis set, COSMO parameters
3. **Scaling factor summary** - Slope, intercept, R², MAE, RMSD for all 26 sets (may go to supplementary)
4. **Data repository contents** - File counts, sizes, formats

**Supplementary materials (online only):**
- Full molecule list with SMILES, InChI, molecular properties (CSV)
- Complete scaling factor table (all 26 sets)
- Individual molecule validation plots
- Detailed solvent parameters
- NWChem input file template with all parameters
- Quality assurance checks performed

### Differentiators

What makes this data descriptor stand out vs generic submissions.

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| **Comprehensive solvent coverage** | 13 solvents vs typical 1-2 | Low | Unique selling point |
| **Benchmark-quality validation** | R² > 0.99 for all sets | Medium | Exceptional accuracy |
| **Full computational provenance** | Every calculation fully documented | Medium | Reproducibility gold standard |
| **Cross-solvent analysis** | Compare performance across conditions | High | Novel insight |
| **Error analysis by functional group** | Identify systematic biases | High | Practical guidance for users |
| **Reuse examples** | Show 3-5 use cases | High | Jupyter notebooks in supplement |
| **Community-standard molecules** | DELTA50 well-known in NMR community | Low | Builds on trusted foundation |
| **Open science exemplar** | All data, code, methods open | Low | Aligns with Nature values |
| **Machine learning ready** | Structured data for ML training | Medium | Attracts ML community |
| **Interoperability** | Links to PubChem, ChEMBL IDs | Medium | Database cross-referencing |

### Anti-Features

What NOT to include in a data descriptor.

| Anti-Feature | Why Avoid | What to Do Instead |
|--------------|-----------|-------------------|
| **Novel scientific hypotheses** | Data descriptors are NOT research papers | Focus on dataset creation, quality, reusability |
| **Extensive analysis/interpretation** | Out of scope for data descriptor | Save for companion research paper |
| **Comparison with other methods** | This is analysis, not description | Mention briefly in context only |
| **> 8 figures** | Journal guideline violation | Move extras to supplementary |
| **Unvalidated data** | Technical validation is mandatory | Include quality metrics, outlier checks |
| **Vague methods** | "Standard procedures" insufficient | Full parameter lists required |
| **Conclusions about chemistry** | Data descriptors don't draw conclusions | Describe potential uses neutrally |
| **Marketing language** | "Revolutionary", "paradigm-shifting" | Stick to factual description |

---

## Deliverable 3: Journal of Cheminformatics - Software Article

### Table Stakes

Required components for a software article.

| Feature | Why Expected | Complexity | Notes |
|---------|--------------|------------|-------|
| **Implementation section** | Architecture and technical details | Medium | FastAPI, Docker, tech stack |
| **Availability statement** | Where to access the tool | Low | URL + GitHub + Docker Hub |
| **Open source license** | Required by journal | Low | Already MIT (verify) |
| **Source code repository** | GitHub with full code | Low | Already exists |
| **Installation instructions** | How to deploy locally | Medium | Docker compose + manual |
| **API documentation** | Endpoint descriptions | Medium | OpenAPI/Swagger (already have?) |
| **Example usage** | Demonstrate core functions | Medium | cURL, Python, web UI examples |
| **Software dependencies** | All libraries with versions | Low | requirements.txt, Dockerfile |
| **Abstract** | Purpose and functionality | Low | Not hypothesis-driven |
| **Introduction** | Problem context | Medium | Why this tool is needed |
| **Availability of data** | Link to DELTA50 dataset | Low | Cross-reference data descriptor |
| **Test data** | Molecules for demonstrating tool | Low | Subset from DELTA50 |
| **Reproducible examples** | Users can verify functionality | Medium | Step-by-step tutorials |

**Recommended sections:**
1. **Background** - NMR prediction landscape, why DELTA50 scaling factors
2. **Implementation** - FastAPI backend, architecture, scaling factor application
3. **Operation** - Web interface, API endpoints, inputs/outputs
4. **Use Cases** - 3-5 practical examples with results
5. **Availability and Requirements** - OS, Docker, cloud deployment
6. **Conclusions** - Summary, future work

**Recommended figures (4-6 total):**
1. **Architecture diagram** - Frontend, API, calculation engine, data flow
2. **Web interface screenshot** - Annotated showing key features
3. **API workflow** - Input → processing → output with example JSON
4. **Use case example** - Molecule structure + prediction results
5. **Performance metrics** - Response time, throughput (if applicable)
6. **Comparison with web form output** - Show both interfaces

**Recommended tables (2-3 total):**
1. **API endpoints** - Route, method, parameters, response format
2. **Software stack** - Component, technology, version, purpose
3. **Example predictions** - Show accuracy on test molecules

**Supplementary materials:**
- Full API reference (if not in main text)
- Tutorial notebook (Jupyter)
- Docker deployment guide
- Additional use case examples
- Video tutorial (if available)

### Differentiators

What makes this software article stand out.

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| **Multi-solvent prediction** | 13 solvents vs typical gas phase | Low | Rare in web tools |
| **Docker deployment** | One-command local deployment | Low | Reproducibility + ease |
| **RESTful API + Web UI** | Serves both programmatic and interactive users | Medium | Dual interface |
| **OpenAPI compliance** | Standards-compliant documentation | Low | Auto-generated from FastAPI |
| **Benchmark-backed** | Predictions from validated dataset | Low | Credibility from data descriptor |
| **Production deployment** | Live tool users can try immediately | Low | Not just prototype |
| **Cloud-ready** | Deployment instructions for cloud | Medium | Scalability demonstration |
| **Example notebooks** | Jupyter tutorials for API usage | Medium | Lowers adoption barrier |
| **Performance benchmarks** | Response time, concurrent users | Medium | Shows production quality |
| **Integration examples** | Using API from R, Python, JavaScript | High | Broad accessibility |
| **Automated testing** | CI/CD with test coverage | Medium | Software quality signal |
| **FAIR-aligned outputs** | Includes InChI, SMILES in responses | Medium | Interoperability |

### Anti-Features

What NOT to include in a software article.

| Anti-Feature | Why Avoid | What to Do Instead |
|--------------|-----------|-------------------|
| **Closed source** | Journal requires open source | MIT/Apache/GPL license |
| **Proprietary dependencies** | Limits reproducibility | Use only OSS libraries |
| **"Available upon request"** | Violates journal policy | Public GitHub + Docker Hub |
| **Windows-only** | Limits accessibility | Cross-platform (Docker solves this) |
| **No test data** | Reviewers can't verify claims | Include example molecules |
| **Vague installation** | Frustrates users | Step-by-step with error handling |
| **Undocumented API** | Poor usability | Full OpenAPI spec required |
| **Anonymous deployment** | Can't review functionality | Must be accessible to reviewers |
| **Novel chemistry claims** | Out of scope for software paper | Focus on implementation, not chemistry discovery |
| **GUI-only** | No programmatic access | Must have API for automation |

---

## Feature Dependencies Between Deliverables

The three deliverables are interdependent and should reference each other:

```
Institutional Repository Dataset
         ↓
    (provides DOI)
         ↓
Nature Scientific Data Descriptor ←──────┐
         ↓                               │
   (describes dataset)              (uses data)
         ↓                               │
         └──→ Journal Cheminformatics ──┘
              Software Article
              (implements tool using
               scaling factors)
```

**Critical dependencies:**
1. **Repository → Data Descriptor**: Data descriptor MUST cite repository DOI in Data Records section
2. **Repository → Software Article**: Software uses scaling factors from repository
3. **Data Descriptor ↔ Software Article**: Should cite each other as companion papers
4. **Timing**: Publish repository first (get DOI), then papers can cite it

**Cross-referencing requirements:**
- Data descriptor "Data Availability" section → Repository DOI
- Software article "Availability of data" section → Repository DOI
- Software article "Background" section → Cite data descriptor paper
- Data descriptor "Usage Notes" section → Mention web tool as example use

---

## Publication Sequence Recommendation

Based on dependencies and typical review timelines:

### Phase 1: Dataset Preparation and Repository Submission
**Timeline: Weeks 1-2**
1. Organize files in RADAR4Chem structure
2. Generate all CSV/JSON files
3. Write comprehensive README
4. Create METADATA.json
5. Submit to RADAR4Chem
6. Obtain DOI (may take 1-2 weeks)

**Deliverable:** Repository DOI

### Phase 2: Parallel Paper Writing
**Timeline: Weeks 3-6** (while repository is under review)
1. Write Scientific Data data descriptor draft
2. Write Journal Cheminformatics software article draft
3. Generate all figures and tables
4. Prepare supplementary materials

**Deliverable:** Two manuscript drafts

### Phase 3: Sequential Submission
**Timeline: Weeks 7+**
1. Submit Scientific Data data descriptor (includes repository DOI)
2. Wait for initial review/revisions (~4-8 weeks)
3. Submit Journal Cheminformatics after data descriptor acceptance (can cite it)

**Rationale for sequence:**
- Data descriptor establishes dataset credibility
- Software article builds on published dataset
- Both can cite the repository immediately
- Software article can cite data descriptor if timed right

---

## MVP Recommendations

For each deliverable, prioritize:

### Repository (MVP)
1. Core data files (structures, raw calculations, processed JSON)
2. Scaling factors with basic metrics
3. Clear README with citation
4. Standard file organization
5. **Defer:** Example Jupyter notebooks, checksum validation

### Data Descriptor (MVP)
1. All required sections (Background, Methods, Data Records, Technical Validation)
2. 4-5 core figures (diversity, workflow, accuracy, correlation)
3. 3 essential tables (summary, computational details, data contents)
4. Basic supplementary (full molecule list, scaling factors)
5. **Defer:** Cross-solvent analysis, error by functional group, ML examples

### Software Article (MVP)
1. Implementation and Operation sections with architecture
2. Web UI and API documentation
3. 3-5 practical examples
4. Live deployment + Docker Hub
5. 3-4 figures (architecture, UI, workflow)
6. **Defer:** Performance benchmarks, integration examples, video tutorial

---

## Quality Metrics by Deliverable

### Repository Quality Indicators
- [ ] Files organized logically with clear naming
- [ ] README explains usage in < 5 minutes
- [ ] All data in multiple formats (JSON + CSV)
- [ ] Persistent identifiers for molecules (InChI, SMILES)
- [ ] License clearly stated
- [ ] Calculation methods fully documented

### Data Descriptor Quality Indicators
- [ ] Technical validation demonstrates R² > 0.99
- [ ] Methods section enables full reproduction
- [ ] Figures are publication-quality (300+ DPI)
- [ ] Data Records section precisely describes repository contents
- [ ] Usage Notes provide clear guidance
- [ ] No novel chemistry claims (descriptive only)

### Software Article Quality Indicators
- [ ] Tool is live and accessible to reviewers
- [ ] Installation works in < 10 minutes (Docker)
- [ ] API is fully documented (OpenAPI)
- [ ] Examples are reproducible
- [ ] Source code is clean and commented
- [ ] Test suite demonstrates functionality

---

## Common Pitfalls Across Deliverables

### Metadata Consistency
**Pitfall:** Molecule identifiers differ between deliverables
**Prevention:** Generate master molecule table, use everywhere

### Version Control
**Pitfall:** Dataset version in repository doesn't match paper
**Prevention:** Freeze dataset before paper submission, assign version

### DOI Timing
**Pitfall:** Papers submitted before repository DOI issued
**Prevention:** Get repository DOI first, then write papers

### Overpromising
**Pitfall:** Claiming dataset/tool does more than it actually does
**Prevention:** Accurate, factual descriptions only

### Reproducibility Gaps
**Pitfall:** Missing software versions, parameters, or dependencies
**Prevention:** Document everything, include environment files

### Format Inaccessibility
**Pitfall:** Only providing data in hard-to-parse formats
**Prevention:** Always include CSV and JSON alongside any specialized formats

---

## Confidence Assessment

| Area | Confidence | Source Quality | Notes |
|------|------------|----------------|-------|
| **Scientific Data requirements** | MEDIUM | Official guidelines + examples | WebSearch verified with multiple publications |
| **RADAR4Chem requirements** | MEDIUM | Official documentation | 10GB limit, metadata requirements clear |
| **Journal Cheminformatics** | MEDIUM | Official software guidelines | No "application note" type, using "Software" article type |
| **Computational chemistry standards** | HIGH | Multiple recent publications | DELTA50, NMRexp, quantum benchmarks reviewed |
| **FAIR principles** | HIGH | Official GO-FAIR + NFDI4Chem | Well-established standards |
| **File formats** | HIGH | Community standards | XYZ, SDF, CSV, JSON are standard |
| **API best practices** | MEDIUM | OpenAPI spec + examples | FastAPI auto-generates compliant docs |

**Low confidence areas requiring validation:**
- Exact figure count preferences for Scientific Data (guideline is "max 8" but optimal unknown)
- Journal Cheminformatics acceptance criteria for web tools vs downloadable software
- RADAR4Chem review timeline and acceptance rate

---

## Sources

### Publication Guidelines
- [Submission Guidelines | Scientific Data](https://www.nature.com/sdata/submission-guidelines)
- [The Data Descriptor – making your data reusable](https://blogs.nature.com/scientificdata/2013/09/19/the-data-descriptor-making-your-data-reusable/)
- [Software | Journal of Cheminformatics](https://jcheminf.biomedcentral.com/submission-guidelines/preparing-your-manuscript/software)
- [RADAR4Chem | NFDI4Chem Knowledge Base](https://knowledgebase.nfdi4chem.de/knowledge_base/docs/radar4chem/)
- [Publish Research Data in a FAIR way in RADAR4Chem](https://www.nfdi4chem.de/publish-research-data-in-a-fair-way-in-radar4chem/)

### Dataset Examples
- [NMRexp: A database of 3.3 million experimental NMR spectra | Scientific Data](https://www.nature.com/articles/s41597-025-06245-5)
- [Quantum chemical benchmark databases of gold-standard dimer interaction energies | Scientific Data](https://www.nature.com/articles/s41597-021-00833-x)
- [DELTA50: A Highly Accurate Database of Experimental 1H and 13C NMR Chemical Shifts Applied to DFT Benchmarking](https://www.mdpi.com/1420-3049/28/6/2449)
- [IR-NMR multimodal computational spectra dataset for 177K patent-extracted organic molecules | Scientific Data](https://www.nature.com/articles/s41597-025-05729-8)
- [The 100-protein NMR spectra dataset | Scientific Data](https://www.nature.com/articles/s41597-023-02879-5)

### Best Practices
- [Best Practices for Constructing, Preparing, and Evaluating Protein-Ligand Binding Affinity Benchmarks](https://pmc.ncbi.nlm.nih.gov/articles/PMC9662604/)
- [Good Practices in Database Generation for Benchmarking Density Functional Theory](https://wires.onlinelibrary.wiley.com/doi/10.1002/wcms.1737)
- [The FAIR Guiding Principles | Scientific Data](https://www.nature.com/articles/sdata201618)
- [FAIR Data Principles | NFDI4Chem Knowledge Base](https://knowledgebase.nfdi4chem.de/knowledge_base/docs/fair/)
- [Reproducible Research in Computational Chemistry of Materials](https://pubs.acs.org/doi/10.1021/acs.chemmater.7b00799)

### Technical Standards
- [Data Format Standard | NFDI4Chem Knowledge Base](https://knowledgebase.nfdi4chem.de/knowledge_base/docs/format_standards/)
- [Open chemistry: RESTful web APIs, JSON, NWChem](https://pmc.ncbi.nlm.nih.gov/articles/PMC5662523/)
- [Guidelines for authors submitting code & software](https://www.nature.com/documents/GuidelinesCodePublication.pdf)
- [Ten simple rules for writing a paper about scientific software](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008390)
- [Cheminformatics Microservice | Journal of Cheminformatics](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-023-00762-4)

### Accuracy Metrics
- [Accurate Prediction of NMR Chemical Shifts: DFT + GNN](https://pmc.ncbi.nlm.nih.gov/articles/PMC11209944/)
- [DFT/NMR Approach for Configuration Assignment](https://pubs.acs.org/doi/10.1021/acs.joc.9b03129)
- [1H and 13C NMR scaling factors for common solvents](https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.23638)
