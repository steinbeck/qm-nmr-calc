# v3.0 Publication & Dataset Release - Research Summary

**Project:** qm-nmr-calc v3.0 - DELTA50 Dataset Publication
**Domain:** Computational chemistry dataset publication with companion papers
**Researched:** 2026-02-11
**Confidence:** HIGH

## Executive Summary

The v3.0 milestone requires publishing a computational chemistry benchmark dataset (DELTA50: 650 NWChem calculations across 50 molecules and 13 solvents) alongside two companion papers that establish its scientific value and demonstrate its application. This is a publication-focused release rather than a software development effort, requiring expertise in FAIR data principles, academic publishing workflows, and scientific metadata standards.

The recommended approach follows a dataset-first strategy: prepare and deposit the dataset to Radar4Chem to obtain a reserved DOI before writing either paper. This enables both manuscripts to cite the dataset during peer review while keeping the data embargoed until acceptance. The critical path involves three sequential deliverables: (1) dataset archive with FAIR-compliant metadata, (2) Nature Scientific Data descriptor establishing the dataset's quality and reuse value, and (3) Journal of Cheminformatics application note demonstrating the qm-nmr-calc tool using the validated dataset. Total timeline: 12-18 weeks from dataset preparation to manuscript submission.

The primary risk is incomplete computational provenance metadata that looks sufficient initially but causes reproducibility failures years later when other researchers attempt to use the dataset. Prevention requires documenting every computational parameter explicitly (NWChem version including patch level, COSMO grid settings, all non-default parameters, compiler versions) and including representative input/output files. Additional risks include FAIR compliance failures reducing dataset discoverability, insufficient technical validation causing manuscript rejection, and inadequate benchmarking undermining software impact claims.

## Key Findings

### Recommended Stack

The v3.0 milestone builds on the existing validated Python/RDKit/NWChem stack by adding publication-specific tools: metadata standards (DataCite 4.6), data packaging (Frictionless Data), LaTeX infrastructure (Springer Nature template), and enhanced visualization (seaborn for publication-quality statistical plots). The stack is minimal by design—no reimplementation of calculation pipelines, only adding packaging and documentation layers.

**Core technologies:**
- **DataCite Metadata Schema 4.6** (Dec 2024): Radar4Chem's required standard — 10 mandatory fields + 13 optional for chemistry datasets, JSON format preferred for Python generation
- **Frictionless Data Package 5.18.1+**: FAIR-compliant dataset packaging with datapackage.json descriptor, integrates with pandas for CSV/Parquet export
- **Springer Nature LaTeX Template** (Dec 2024): Unified template for both journals (Nature Scientific Data and Journal of Cheminformatics), Overleaf integration available
- **Parquet format**: Archival format for benchmark results preserving data types (better than CSV), PyArrow integration with pandas
- **seaborn (latest)**: Statistical visualizations for publication (correlation heatmaps, distribution plots), built on existing matplotlib stack
- **Radar4Chem repository**: NFDI4Chem domain repository, 10GB free tier sufficient, 25-year preservation guarantee, manual upload only (no API)

**Critical version/format decisions:**
- DataCite 4.6 (Dec 2024) NOT legacy versions—explicit NFDI4Chem requirement
- JSON metadata generation NOT XML—easier Python integration despite XML being canonical
- Springer Nature template (Dec 2024) for BOTH papers—Nature Scientific Data explicitly warns against legacy templates
- Parquet for processed results, CSV for scaling factors—balance between type preservation and human readability
- JCAMP-DX for NMR data optional—IUPAC standard adds domain value but uncertain applicability to DFT computational data

**Installation:** `uv add frictionless seaborn pyarrow` (3 packages added to existing environment)

### Expected Features

Three distinct deliverables with interdependent requirements. Dataset DOI must be obtained first to enable paper citations during peer review.

**Deliverable 1: Radar4Chem Dataset (Table Stakes)**
- Persistent DOI with 10 mandatory DataCite fields + chemistry-specific extensions
- README.md with dataset overview, citation, usage instructions
- 650 NWChem calculations organized by molecule/solvent hierarchy (input .nw files + output .out files)
- Processed data in multiple formats: Parquet (benchmark results), CSV (scaling factors, molecule metadata), JSON (shielding tensors)
- SMILES/InChI/InChIKey identifiers for all 50 molecules (enables cross-database linking, required for FAIR F4)
- Scaling factors for 26 parameter sets (1H/13C × 13 solvents) with validation metrics
- LICENSE.txt (CC-BY 4.0 recommended) and version number (v1.0.0)
- Complete computational provenance: exact NWChem version, all parameters, COSMO grid settings, dielectric constants

**Deliverable 1: Radar4Chem Dataset (Differentiators)**
- 13 explicit solvents (rare—most benchmarks use gas phase or 1-2 solvents)
- High-quality validation: R² > 0.99 for all parameter sets, MAE < 0.15 ppm for 1H
- Machine learning ready: structured data formats, comprehensive molecular descriptors
- Computation time logs extracted from .out files (helps users estimate resources)
- Example Jupyter notebooks showing dataset usage (tutorial for users)
- Checksums (MD5/SHA256) for file integrity verification

**Deliverable 2: Nature Scientific Data Descriptor (Table Stakes)**
- Background & Summary (2-3 paragraphs): motivation and reuse value
- Methods section (3-4 pages): detailed calculation protocol enabling full reproduction (DFT functional, basis set, COSMO parameters with explicit grid settings, SCF convergence criteria)
- Data Records section: precise repository description with Radar4Chem DOI, file formats, folder structure
- Technical Validation (2-3 pages): quantitative metrics (RMSE, R², MAE), statistical analysis, independent recalculation verification, outlier analysis, comparison with literature benchmarks
- Usage Notes: how to apply scaling factors, software requirements, chemical applicability
- Code Availability statement: GitHub + Zenodo DOI for processing scripts
- Data Availability statement: Radar4Chem DOI with CC-BY license specification
- 5-6 figures maximum: molecule diversity chart, workflow diagram, accuracy heatmaps, correlation plots (experimental vs predicted), statistical validation
- ORCID for all authors, competing interests, CRediT author contributions

**Deliverable 2: Scientific Data Descriptor (Differentiators)**
- Comprehensive solvent coverage analysis (how accuracy varies across 13 solvents)
- Benchmark-quality validation with exceptional accuracy (R² > 0.99 consistently)
- Full computational provenance for reproducibility gold standard
- Error analysis by functional group (practical guidance: which molecular features predict poorly?)
- Open science exemplar (data + code + methods fully open)

**Deliverable 3: Journal of Cheminformatics Application Note (Table Stakes)**
- Implementation section: FastAPI architecture, Docker deployment, technology stack
- Availability statement: Live URL + GitHub + Docker Hub (open source MIT license verified)
- API documentation: OpenAPI/Swagger auto-generated from FastAPI
- Example usage: cURL commands, Python client code, web UI demonstrations
- Installation instructions: Docker compose (one-command deployment) + manual setup
- Abstract and Introduction: problem context (why DELTA50 scaling factors needed?)
- Test data: subset of DELTA50 molecules for demonstrating functionality
- 4-6 figures: architecture diagram, web interface screenshot, API workflow, use case examples

**Deliverable 3: Application Note (Differentiators)**
- Multi-solvent prediction (13 solvents vs typical gas phase)
- Dual interface: RESTful API + interactive web UI (serves both programmatic and human users)
- Docker deployment (one-command local installation, reproducibility)
- OpenAPI compliance (standards-compliant auto-generated documentation)
- Benchmark-backed predictions (credibility from validated Scientific Data dataset)
- Production deployment ready (live tool users can try immediately, not just prototype)

**Anti-Features (Explicitly Avoid)**
- Dataset: Proprietary formats only, undocumented acronyms, mixed calculation versions, "available upon request"
- Data Descriptor: Novel scientific hypotheses (out of scope), >8 figures, vague methods using "standard settings," marketing language
- Application Note: Closed source, Windows-only, no test data, undocumented API, >2 page limit

**Publication Sequence (Critical Dependency Chain)**
1. Dataset preparation → Radar4Chem deposit → Reserved DOI (required for papers)
2. Parallel manuscript writing (both cite reserved DOI from start)
3. Scientific Data submission (establishes dataset credibility)
4. Journal Cheminformatics submission (can cite "submitted" descriptor, update to published after acceptance)
5. Coordinate acceptance → publish dataset (DOI goes live) → update application note citation

### Architecture Approach

Three-deliverable structure with dataset-first build order to obtain DOI before writing papers. Key architectural decision is embargo coordination: dataset uploaded and DOI reserved immediately, but data remains private during peer review (4-8 months), then published when both papers are accepted. This maintains novelty while enabling consistent citations.

**Directory structure:**
```
qm-nmr-calc/
├── data/benchmark/delta50/          # Existing source data
└── publications/                    # NEW: All publication deliverables
    ├── dataset/                     # Deliverable 1: Radar4Chem archive
    │   ├── README.md
    │   ├── molecules/               # 50 XYZ structures
    │   ├── calculations/            # 650 NWChem files (solvent/molecule)
    │   ├── results/                 # Parquet, CSV, JSON processed data
    │   └── metadata.json            # DataCite 4.6 + chemistry extensions
    ├── data-descriptor/             # Deliverable 2: Scientific Data paper
    │   ├── manuscript.tex
    │   ├── figures/                 # Development: organized subdirs
    │   └── submission/              # Flattened for journal submission
    └── application-note/            # Deliverable 3: J Cheminf paper
        ├── manuscript.tex
        ├── figures/
        └── submission/
```

**Major components:**
1. **Dataset Package** — Organizes 650 calculations + metadata for archival, produces Radar4Chem archive (tar.gz) + reserved DOI, consumed by both papers
2. **Data Descriptor** — Describes dataset creation/validation/reuse, inputs dataset DOI and methodology, produces Scientific Data manuscript for research community
3. **Application Note** — Describes qm-nmr-calc tool demonstrated with dataset, inputs dataset DOI + data descriptor citation + tool architecture, produces J Cheminf manuscript for tool users
4. **Build System** — Generates submission packages, validates cross-references (dataset DOI consistent across papers, shared bibliography), produces submission-ready archives for journal portals

**Critical architectural patterns:**
- **DOI Reservation Before Writing** (Pattern 1): Reserve dataset DOI immediately after archive preparation, enables consistent citation during peer review, avoids last-minute updates
- **Submission Package Generation** (Pattern 2): Automate flattening directory structure for journal submission (LaTeX systems cannot process nested directories)
- **Cross-Reference Validation** (Pattern 3): Automated checks that dataset DOI appears identically in both papers before submission
- **Embargo Coordination** (Pattern 4): Dataset embargoed throughout peer review (12-18 weeks), published only when both papers accepted, reviewers access via private links

**Anti-patterns to avoid:**
- Writing papers before dataset DOI exists (creates coordination debt)
- Nested subdirectories in submission package (causes LaTeX system failures)
- Publishing dataset before papers ready (loses coordination timing, reviewers can't request dataset changes)
- Independent bibliography files (citation inconsistency across papers)
- Dataset archive without README (violates FAIR reusability)

**Open questions requiring Phase 1 resolution:**
- Does Radar4Chem support DOI reservation for embargoed datasets? (Confidence: LOW—not in docs, fallback: use Zenodo/Figshare which explicitly support this)
- Exact Radar4Chem metadata.json schema requirements? (Confidence: LOW—found general guidance not complete schema, fallback: use Dublin Core standard)

### Critical Pitfalls

**Pitfall 1: Incomplete Computational Provenance Metadata** (CRITICAL—causes long-term reproducibility failures)
- **What goes wrong:** Dataset deposited with generic computational description ("B3LYP/6-311+G(2d,p), COSMO") but missing critical parameters. Other researchers attempt reproduction, get different results (RMSE differs by 0.1-0.2 ppm), dataset credibility damaged.
- **Why dangerous:** Looks fine during manuscript preparation but causes failures years later. NWChem numerical results depend on: exact patch version (7.2.0 vs 7.2.2), COSMO grid parameters (minbem/maxbem/ificos—defaults produce large errors), dielectric constants per solvent, compiler version (GCC vs Intel affects precision), BLAS/LAPACK libraries.
- **Prevention:** Include complete NWChem input files for 5-10 representative molecules, document ALL non-default parameters, provide complete output files enabling validation, create computational methods checklist (software exact version including patch, functional+basis, solvation with all parameters, convergence criteria, computational environment with compiler/MPI/libraries), test cross-environment reproducibility.
- **Detection:** Methods draft uses "standard settings" or "default parameters," collaborators can't answer parameter questions without checking old scripts, no input files version-controlled.

**Pitfall 2: FAIR Principle F1/F4 Violations—Missing Persistent Identifier Links and Chemistry Metadata** (CRITICAL—reduces discoverability)
- **What goes wrong:** Dataset has DOI, papers cite DOI, but metadata doesn't create bidirectional RelatedIdentifier links. Dataset repository doesn't list papers, papers don't appear in DataCite relationships. Separately: metadata completes minimal DataCite fields but skips chemistry-specific extensions (molecular structures, technique ontology terms), dataset invisible to ChemSpider/PubChem/chemistry aggregators.
- **Prevention:** Update dataset metadata when papers accepted (add RelatedIdentifier with IsCitedBy relationship type), include formal dataset citation in papers' References section, configure Crossref deposits with relationship tags. Complete ALL Radar4Chem chemistry-specific fields (SMILES/InChI for 50 molecules, technique tags using ontology terms, software structured entry), include SDF/CSV structure files, follow RADAR Metadata Schema v9.2 chemistry extensions, validate with NFDI4Chem FAIR checklist.
- **Detection:** Planning to "just mention DOI in text," no plan to update dataset metadata post-publication, metadata preview shows only generic DataCite fields (no chemistry extensions), Radar4Chem form completed in <30 minutes (likely skipped optional fields).

**Pitfall 3: Nature Scientific Data—Insufficient Technical Validation Section** (CRITICAL—causes rejection/major revision)
- **What goes wrong:** Technical Validation describes "we checked the data" qualitatively but lacks quantitative validation metrics, statistical analysis, independent verification. Reviewers flag as insufficient, major revisions required, 3-6 month delay.
- **What reviewers expect:** (1) Quantitative metrics: RMSE/MAE/R² for all 650 calculations stratified by nucleus and solvent, statistical distributions, outlier detection; (2) Cross-validation: independent validation set separate from OLS scaling training, comparison with alternative DFT methods; (3) Outlier analysis: statistical tests, chemical explanations for failures; (4) Reproducibility verification: independent recalculation of 5-10 molecules by different team member yielding identical results.
- **Prevention:** Perform quantitative validation during data generation (not after), calculate per-solvent and per-nucleus metrics, include statistical visualizations (correlation plots, residual histograms, error box plots), compare with literature NMR prediction benchmarks (cite specific accuracy comparisons), document independent recalculation verification.
- **Detection:** Technical Validation <500 words (too brief), no quantitative metrics (just "good agreement"), no figures/tables in validation section, no literature comparison, phrases like "carefully checked" without HOW.

**Pitfall 4: Journal of Cheminformatics—Inadequate Benchmarking** (CRITICAL—demonstrates insufficient advantage)
- **What goes wrong:** Application note presents tool but compares against 1-2 convenient baselines instead of comprehensive state-of-the-art survey. Reviewers reject: "does not demonstrate significant advance over previously published software." 6+ month delay to perform additional benchmarks.
- **What reviewers expect:** (1) State-of-the-art comparison: multiple DFT functionals, recent ML approaches, commercial/free tools (ChemDraw, nmrshiftdb2); (2) Multiple benchmark datasets: DELTA50 plus 2-3 independent literature datasets, diverse chemical classes; (3) Statistical rigor: RMSE with confidence intervals, significance tests, per-molecule error analysis, outlier analysis; (4) Computational cost comparison: prediction time, scalability, hardware requirements.
- **Prevention:** Survey NMR prediction papers (last 5 years) to identify current best methods, design benchmark matrix (methods × datasets × metrics), perform fair comparisons (same molecules for all methods, document when competitors fail), include accuracy-cost tradeoff discussion.
- **Detection:** Benchmark table <3 comparison methods, no recent ML approaches, testing only on own dataset, no failure case discussion, no cost comparison.

**Pitfall 5: J. Cheminformatics—Underselling Tool Impact and Novelty** (MODERATE—reduces impact if published)
- **What goes wrong:** Application note focuses on technical implementation ("React web interface") instead of scientific value proposition ("enables NMR structure elucidation for non-experts"). Reviewers assess as "incremental contribution, limited impact."
- **Prevention:** Lead with scientific value proposition (what problem solved? target users? capabilities enabled?), include concrete use cases (structure elucidation workflow, diastereomer assignment, solvent selection), quantify impact metrics (time savings: 48 CPU-hours → 30 seconds; accuracy improvement: 95% correct vs 70% free alternatives; usability: no computational expertise required), highlight novelty explicitly (first benchmark with explicit solvent effects, first web-accessible DFT-accuracy predictor).
- **Detection:** Introduction focuses on implementation details not scientific problem, no user personas or use cases, comparison table lists features (# solvents) not outcomes (accuracy improvement), writing sounds like technical documentation not scientific paper.

## Implications for Roadmap

Based on research, the v3.0 milestone naturally divides into 6 sequential phases with specific synchronization points. The critical path follows the dataset-first strategy with parallel manuscript development.

### Phase 1: Dataset Archive Preparation
**Rationale:** Must prepare complete dataset archive before DOI reservation. This phase organizes existing computational outputs (650 NWChem calculations from v2.x milestones) into FAIR-compliant structure, generates all metadata, and validates completeness. Blocking dependency for all subsequent phases.

**Delivers:**
- publications/dataset/ directory with DELTA50-benchmark archive structure
- README.md and METHODOLOGY.md documentation
- All 650 NWChem .nw inputs and .out outputs organized by molecule/solvent hierarchy
- Processed results exported: Parquet (benchmark results), CSV (scaling factors, molecule metadata), JSON (shielding tensors)
- metadata.json with DataCite 4.6 mandatory fields + chemistry extensions (SMILES/InChI for all molecules)
- Computational provenance checklist: exact NWChem version, all parameters documented, COSMO grid settings explicit, environment details (compiler, MPI, libraries)
- Quality validation: checksums for all files, file manifest CSV

**Addresses features:**
- Dataset table stakes: persistent identifier-ready, complete calculation files, processed data formats, molecular identifiers
- Dataset differentiators: computation time logs, quality metrics calculation

**Avoids pitfalls:**
- Pitfall 1 (incomplete provenance): Implements computational methods checklist with all parameters
- Pitfall 1.6 (file naming): Establishes systematic naming convention from start
- Pitfall 1.4 (versioning): Includes version number (v1.0.0) in all files

**Estimated duration:** 1-2 weeks

### Phase 2: Repository DOI Reservation
**Rationale:** Obtain reserved DOI before manuscript writing to enable consistent citation during peer review. Resolves Open Question 1 (Radar4Chem embargo capability) with fallback to Zenodo/Figshare if needed.

**Delivers:**
- Dataset uploaded to Radar4Chem (or fallback repository) marked as embargoed
- Reserved DOI documented in publications/RESERVED-DOI.txt
- Complete metadata submitted: 10 mandatory DataCite fields + 13 optional + chemistry-specific extensions
- Private reviewer access link (for paper peer review)
- Verification: DOI resolves to embargo page with metadata

**Addresses features:**
- Dataset table stakes: persistent DOI obtained
- FAIR F1 compliance: identifier assigned early

**Avoids pitfalls:**
- Anti-pattern: Writing papers before DOI exists
- Pitfall 1.2 (FAIR F1): Establishes DOI infrastructure for bidirectional linking
- Pitfall 1.3 (FAIR F4): Completes all chemistry-specific metadata fields in repository

**Research flags:** HIGH—Needs resolution of Radar4Chem embargo capability (Open Question 1). If unsupported, requires quick pivot to Zenodo/Figshare.

**Estimated duration:** 1-3 days (upload) + 1-3 days (DOI assignment)

### Phase 3: Data Descriptor Manuscript (Nature Scientific Data)
**Rationale:** Write data descriptor first to establish dataset credibility. Can proceed in parallel with Phase 4 (application note) after Phase 2 completes. Data descriptor is foundational—application note will cite it.

**Delivers:**
- publications/data-descriptor/manuscript.tex with all required sections:
  - Background & Summary (2-3 paragraphs)
  - Methods (3-4 pages) with complete reproducibility detail
  - Data Records (1-2 pages) describing Radar4Chem deposit
  - Technical Validation (2-3 pages) with quantitative metrics
  - Usage Notes, Code Availability, Data Availability statements
- 5-6 publication-quality figures (300+ DPI): molecule diversity, workflow diagram, accuracy heatmaps, correlation plots, statistical validation
- 3-4 tables: dataset summary, computational details, repository contents
- Supplementary materials: full molecule list, complete scaling factor table, validation plots, NWChem input template
- Submission package: flattened directory structure (Pattern 2) ready for journal upload
- Cross-reference validation: dataset DOI cited correctly, all references formatted

**Addresses features:**
- Descriptor table stakes: all required sections, figures within limit (5-6 of max 8), ORCID for authors, Data Availability statement
- Descriptor differentiators: comprehensive solvent analysis, benchmark-quality validation, full provenance, error analysis by functional group

**Avoids pitfalls:**
- Pitfall 2.1 (insufficient validation): Implements comprehensive Technical Validation with quantitative metrics, independent recalculation, outlier analysis, literature comparison
- Pitfall 2.2 (methods detail): Documents complete computational workflow with step-by-step protocol, flowchart, all processing scripts
- Pitfall 2.3 (DAS format): Uses correct Data Availability Statement format with all required elements
- Pitfall 2.4 (code availability): Deposits processing scripts to Zenodo, provides Code Availability statement
- Pitfall 2.5 (figure quality): Follows Nature figure guidelines (300+ DPI, colorblind-friendly, sized for column width)

**Uses stack:**
- Springer Nature LaTeX Template (Dec 2024) via Overleaf
- seaborn for statistical visualizations
- matplotlib for publication figures (300 DPI TIFF export)
- RDKit for molecular structure images

**Research flags:** MEDIUM—Phase-specific research on "standard validation analyses for NMR datasets" recommended. Review 5-10 recent Scientific Data computational chemistry descriptors for Technical Validation patterns.

**Estimated duration:** 2-3 weeks

### Phase 4: Application Note Manuscript (Journal of Cheminformatics)
**Rationale:** Write application note in parallel with data descriptor (both can proceed after Phase 2 DOI reservation). Application note demonstrates qm-nmr-calc tool using DELTA50 dataset as validation, cites data descriptor as "submitted" initially.

**Delivers:**
- publications/application-note/manuscript.tex with sections:
  - Abstract (<170 words)
  - Background (problem context: why DELTA50 scaling needed?)
  - Implementation (FastAPI architecture, Docker, tech stack)
  - Operation (web UI, API endpoints, inputs/outputs)
  - Use Cases (3-5 practical examples: structure elucidation, diastereomer assignment, solvent selection)
  - Availability and Requirements (URL, GitHub, Docker Hub, system requirements)
  - Conclusions
- 4-6 figures: architecture diagram, web interface screenshot, API workflow, use case examples, performance metrics
- 2-3 tables: API endpoints, software stack, example predictions with accuracy metrics
- Comprehensive benchmarking section: comparison with 4-6 state-of-the-art methods (DFT functionals, ML approaches, commercial/free tools), multiple datasets (DELTA50 + 2-3 independent), statistical rigor (RMSE with confidence intervals, significance tests), computational cost comparison
- Submission package: flattened structure ready for journal upload
- Cross-reference validation: dataset DOI cited, data descriptor cited as "submitted to Scientific Data" (update to published reference after acceptance)

**Addresses features:**
- Application note table stakes: implementation section, availability statement, open source license, API documentation, example usage, test data
- Application note differentiators: multi-solvent prediction, dual interface (API + web UI), Docker deployment, OpenAPI compliance, benchmark-backed, production-ready

**Avoids pitfalls:**
- Pitfall 3.1 (inadequate benchmarking): Implements comprehensive benchmark matrix with state-of-the-art comparisons, multiple datasets, statistical rigor, cost analysis
- Pitfall 3.2 (underselling impact): Leads with scientific value proposition, includes concrete use cases, quantifies impact metrics, highlights novelty explicitly
- Pitfall 3.3 (page limits): Checks J Cheminf requirements before writing, writes to length from start, moves details to supplement
- Pitfall 3.4 (tool naming): Verifies "DELTA50" name is distinctive and searchable
- Pitfall 4.2 (citation consistency): Uses identical dataset citation as data descriptor

**Uses stack:**
- Springer Nature LaTeX Template (Dec 2024)
- Existing qm-nmr-calc architecture documentation
- Shared bibliography with data descriptor

**Research flags:** MEDIUM—Phase-specific research on "current state-of-the-art in NMR prediction" recommended for comprehensive benchmarking. Survey recent J Cheminf application notes for structure/style expectations.

**Estimated duration:** 2-3 weeks (can overlap with Phase 3)

### Phase 5: Parallel Manuscript Submission and Peer Review
**Rationale:** Submit both papers in parallel after manuscripts complete. Data descriptor can be submitted immediately; application note cites descriptor as "submitted." Dataset remains embargoed during peer review (4-8 weeks per paper). Coordinate revisions and track acceptance timing.

**Delivers:**
- Scientific Data submission confirmation
- Journal Cheminformatics submission confirmation
- Tracking spreadsheet: paper status, revision deadlines, reviewer comments
- Coordinated revision responses (ensure cross-references remain valid)
- Acceptance notifications for both papers (timing may differ by weeks/months)

**Addresses features:**
- Publication sequence critical dependency: both papers cite embargoed dataset DOI during review
- Embargo coordination: reviewers access dataset via private link, data remains unpublished

**Avoids pitfalls:**
- Pitfall 4.1 (metadata update timing): Sets calendar reminders to update dataset metadata when each paper accepts
- Anti-pattern: Publishing dataset before papers ready

**Estimated duration:** 4-8 weeks per paper (typical academic review cycle), likely sequential acceptances

### Phase 6: Coordinated Publication and Metadata Finalization
**Rationale:** When both papers are accepted, publish the dataset (DOI goes live), update data descriptor citation in application note from "submitted" to published reference, finalize bidirectional metadata linking. This is the final synchronization point ensuring all three deliverables reference each other correctly.

**Delivers:**
- Dataset published on Radar4Chem (embargo lifted, DOI active)
- Dataset metadata updated with RelatedIdentifier entries for both papers (IsCitedBy relationships)
- Application note citation updated: data descriptor reference changes from "submitted" to volume/pages/DOI
- Verification: Crossref and DataCite pages show bidirectional links (dataset lists papers, papers list dataset)
- Archive preprints (arXiv or institutional repository) for open access
- Announcement: project website updated, mailing lists notified, social media (if applicable)

**Addresses features:**
- Cross-publication management: bidirectional DOI linking complete
- FAIR F1 compliance: persistent identifiers fully integrated
- Long-term maintenance plan: institutional contact email, annual review calendar reminder set

**Avoids pitfalls:**
- Pitfall 1.2 (FAIR F1): Implements bidirectional RelatedIdentifier updates within 1 week of each acceptance
- Pitfall 4.1 (metadata timing): Executes planned metadata updates immediately upon acceptance
- Pitfall 4.2 (citation consistency): Verifies identical dataset citations across papers
- Pitfall 4.3 (maintenance plan): Documents dataset maintenance responsibility, sets annual review

**Estimated duration:** 1-2 weeks

### Phase Ordering Rationale

- **Dataset-first strategy (Phases 1-2 before 3-4):** Research consistently recommends obtaining dataset DOI before manuscript writing to enable consistent citation and avoid coordination debt. Radar4Chem embargo capability (Open Question 1) enables this strategy—if embargo not supported, fallback to Zenodo/Figshare which explicitly allow DOI reservation.

- **Parallel manuscript development (Phases 3-4):** Data descriptor and application note can be written simultaneously after DOI reservation. They address different audiences (dataset users vs tool users) with minimal overlap. Application note cites descriptor as "submitted" initially, updated to published after acceptance.

- **Sequential publication (Phase 6 after Phase 5):** Dataset remains embargoed throughout peer review to maintain novelty for papers and allow reviewers to request changes. Publishing dataset only after both papers accept ensures coordinated release and clean cross-reference metadata.

- **Pitfall avoidance through ordering:** Phase 1 addresses computational provenance (Pitfall 1) by creating checklist before any submission. Phase 2 addresses FAIR compliance (Pitfalls 1.2, 1.3) by completing metadata early. Phases 3-4 address manuscript quality (Pitfalls 2.1, 2.2, 3.1, 3.2) through comprehensive documentation and benchmarking. Phase 6 addresses cross-publication management (Pitfalls 4.1, 4.2) through coordinated metadata updates.

### Research Flags

**Phases needing deeper research during planning:**
- **Phase 2 (Repository DOI Reservation):** HIGH priority—must resolve Radar4Chem embargo capability (Open Question 1) before proceeding. If Radar4Chem doesn't support embargo, requires rapid pivot to Zenodo/Figshare. Suggest direct contact with Radar4Chem support or test with dummy upload. Estimated research time: 1-2 days.

- **Phase 3 (Data Descriptor Writing):** MEDIUM priority—recommend phase-specific research on "standard Technical Validation analyses for NMR computational datasets." Review 5-10 recent Nature Scientific Data descriptors for computational chemistry to identify validation patterns reviewers expect. Estimated research time: 1 day.

- **Phase 4 (Application Note Writing):** MEDIUM priority—recommend phase-specific research on "current state-of-the-art in NMR chemical shift prediction" to ensure comprehensive benchmarking (Pitfall 3.1). Survey recent papers (2021-2026) to identify best-performing DFT functionals and emerging ML approaches. Estimated research time: 1 day.

**Phases with standard patterns (skip research-phase):**
- **Phase 1 (Dataset Archive Preparation):** Well-documented FAIR data practices, established file organization patterns, standard metadata schemas. No additional research needed.

- **Phase 5 (Manuscript Submission):** Standard academic publishing workflow, journal submission systems well-documented. No additional research needed.

- **Phase 6 (Publication Coordination):** Administrative/procedural task following documented Crossref/DataCite metadata update procedures. No additional research needed.

## Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| Stack | HIGH | All critical technologies verified with official documentation: DataCite 4.6 schema (Dec 2024), Springer Nature template (Dec 2024), Frictionless Data spec, Radar4Chem requirements from NFDI4Chem knowledge base. Python library versions current as of Feb 2026. |
| Features | MEDIUM | Based on official journal submission guidelines (Nature Scientific Data, Journal of Cheminformatics) and multiple published examples. Exact figure count preferences and review criteria somewhat uncertain—guidelines provide ranges (e.g., "max 8 figures") but not optimal counts. |
| Architecture | MEDIUM | Dataset-first build order and embargo strategy inferred from multiple sources (Zenodo/Figshare documentation, general publishing best practices) but not explicitly documented for Radar4Chem. LOW confidence on Radar4Chem embargo capability (Open Question 1)—not found in documentation, requires Phase 2 resolution. Directory structure and submission packaging well-established patterns. |
| Pitfalls | HIGH | Computational provenance requirements from published reproducibility studies (Chemistry of Materials, J. Chem. Info. Model.). FAIR compliance failures from official NFDI4Chem and DataCite documentation. Manuscript pitfalls from journal submission guidelines and reviewer guidance documents. Real-world evidence provided for critical pitfalls. |

**Overall confidence:** HIGH

Research based on official sources (DataCite 4.6 Dec 2024 schema, NFDI4Chem Knowledge Base, Nature/Springer submission guidelines, published reproducibility studies) with high-confidence technical recommendations. Medium-confidence areas (journal review criteria, Radar4Chem embargo) are procedural rather than technical and can be resolved during execution through direct communication with repositories/journals.

### Gaps to Address

**HIGH priority (blocking Phase 2):**
- **Radar4Chem embargo capability:** Does Radar4Chem support reserving DOIs for embargoed datasets? Not found in documentation. User leads NFDI4Chem infrastructure—suggest direct contact with Radar4Chem support in Phase 2. Fallback: Zenodo/Figshare confirmed to support embargo + DOI reservation.
- **Resolution approach:** Test with dummy dataset upload in early Phase 2. If embargo not supported, pivot to Zenodo/Figshare within 1 day (both meet NFDI4Chem FAIR requirements, domain-appropriate for computational chemistry).

**MEDIUM priority (informative for Phase 3-4):**
- **Scientific Data figure count:** Is 8 figures a hard limit or recommendation? Documentation says "recommended, not mandate" but actual reviewer expectations unclear.
- **Resolution approach:** Review 10 recent Scientific Data computational chemistry descriptors (2024-2026) to identify typical figure counts. Plan 5-6 figures to stay safely under limit with buffer for reviewer-requested additions.

- **Journal Cheminformatics article type:** Does "Software" article type fit better than "Research" for application note? Guidelines mention Software type but requirements unclear.
- **Resolution approach:** Review recent J Cheminf software articles for qm-nmr-calc-similar web tools. Confirm article type with journal during submission process.

**LOW priority (optimization, non-blocking):**
- **Benchmark sufficiency thresholds:** How many comparison methods constitute "adequate" benchmarking for J Cheminformatics? Guidelines require "demonstration of advantage" but don't specify minimum comparisons.
- **Resolution approach:** Aim for 4-6 comparison methods (multiple DFT functionals + ML approach + free tool) based on recent J Cheminf software papers. Reviewers will request additional comparisons if insufficient.

- **JCAMP-DX applicability:** Is JCAMP-DX format appropriate for DFT computational data vs experimental NMR? IUPAC standard but designed for instrument data.
- **Resolution approach:** Optional format—include if straightforward to generate from existing data, otherwise skip. Parquet + CSV sufficient for computational dataset. Consult NFDI4Chem format preferences in Phase 1.

## Sources

### Primary (HIGH confidence)

**NFDI4Chem and Radar4Chem:**
- [NFDI4Chem FAIR Data Principles](https://knowledgebase.nfdi4chem.de/knowledge_base/docs/fair/) — FAIR compliance requirements
- [RADAR4Chem Documentation](https://knowledgebase.nfdi4chem.de/knowledge_base/docs/radar4chem/) — Repository requirements, metadata schema v9.2
- [NFDI4Chem Metadata / MIChI](https://knowledgebase.nfdi4chem.de/knowledge_base/docs/metadata/) — Chemistry-specific metadata extensions
- [RADAR4Chem About Page](https://radar.products.fiz-karlsruhe.de/en/radarabout/radar4chem) — Storage limits, preservation policy

**DataCite and Metadata Standards:**
- [DataCite Metadata Schema 4.6 Documentation](https://datacite-metadata-schema.readthedocs.io/en/4.6/) — Official schema (Dec 2024)
- [DataCite Schema Versions](https://schema.datacite.org/versions.html) — Version history
- [Crossref Relationships Documentation](https://www.crossref.org/documentation/schema-library/markup-guide-metadata-segments/relationships/) — RelatedIdentifier usage
- [DataCite RelatedIdentifier Property](https://support.datacite.org/docs/datacite-metadata-schema) — Bidirectional linking

**Journal Submission Guidelines:**
- [Nature Scientific Data Submission Guidelines](https://www.nature.com/sdata/submission-guidelines) — Required sections, figure limits
- [Nature Scientific Data Data Policies](https://www.nature.com/sdata/policies/data-policies) — Data Availability Statement requirements
- [Nature Scientific Data - For Referees](https://www.nature.com/sdata/policies/for-referees) — Reviewer evaluation criteria
- [Journal of Cheminformatics Submission Guidelines](https://jcheminf.biomedcentral.com/submission-guidelines) — Software article requirements
- [Journal of Cheminformatics Fees and Funding](https://jcheminf.biomedcentral.com/submission-guidelines/fees-and-funding) — APC information

**LaTeX and Publication Tools:**
- [Springer Nature LaTeX Author Support](https://www.springernature.com/gp/authors/campaigns/latex-author-support) — Official template (Dec 2024)
- [Springer Nature LaTeX Template on Overleaf](https://www.overleaf.com/latex/templates/springer-nature-latex-template/myxmhdsbzkyd) — Cloud LaTeX environment
- [Frictionless Data Framework](https://framework.frictionlessdata.io/) — Data packaging specification

**Computational Reproducibility:**
- [Reproducible Research in Computational Chemistry of Materials](https://pubs.acs.org/doi/10.1021/acs.chemmater.7b00799) — Provenance metadata requirements, real-world reproducibility failures
- [Method and Data Sharing and Reproducibility](https://pubs.acs.org/doi/10.1021/acs.jcim.0c01389) — J. Chem. Info. Model. editorial on reproducibility standards
- [NWChem COSMO Solvation Model Documentation](https://nwchemgit.github.io/COSMO-Solvation-Model.html) — Parameter requirements (minbem, maxbem, ificos)

### Secondary (MEDIUM confidence)

**Dataset Publication Examples:**
- [NMRexp: A database of 3.3 million experimental NMR spectra | Scientific Data](https://www.nature.com/articles/s41597-025-06245-5) — Recent NMR dataset descriptor example
- [DELTA50: A Highly Accurate Database of Experimental 1H and 13C NMR Chemical Shifts Applied to DFT Benchmarking](https://www.mdpi.com/1420-3049/28/6/2449) — Original DELTA50 benchmark paper
- [Quantum chemical benchmark databases of gold-standard dimer interaction energies | Scientific Data](https://www.nature.com/articles/s41597-021-00833-x) — Computational chemistry dataset example

**Benchmarking Best Practices:**
- [Best Practices for Constructing, Preparing, and Evaluating Protein-Ligand Binding Affinity Benchmarks](https://pmc.ncbi.nlm.nih.gov/articles/PMC9662604/) — Benchmarking methodology
- [MoleculeNet Benchmarking Framework](https://pubs.rsc.org/en/content/articlehtml/2018/sc/c7sc02664a) — Cheminformatics benchmarking standard

**DOI Reservation and Versioning:**
- [Zenodo DOI Reservation](https://help.zenodo.org/docs/deposit/describe-records/reserve-doi/) — Confirmed embargo + DOI reservation support
- [Figshare DOI Reservation](https://info.figshare.com/user-guide/how-to-reserve-a-doi/) — Alternative repository with embargo
- [DataCite Versioning Guidelines](https://support.datacite.org/docs/versioning) — Dataset version management

### Tertiary (LOW confidence—needs validation)

**Radar4Chem Embargo Capability:** Not found in documentation. Requires Phase 2 verification through direct contact or test upload.

**JCAMP-DX for Computational Data:** IUPAC standard documented for experimental NMR, applicability to DFT computational data uncertain. Optional format—skip if ambiguous.

**Nature Scientific Data Figure Limit Enforcement:** Guidelines state "max 8 figures recommended" but whether this is strictly enforced or flexible for complex datasets unclear from documentation alone.

---
*Research completed: 2026-02-11*
*Ready for roadmap: yes*
