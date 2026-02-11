# Requirements: qm-nmr-calc v3.0

**Defined:** 2026-02-11
**Core Value:** Reliable async NMR predictions with full control over calculation parameters -- submit a molecule, get back accurate 1H/13C shifts without babysitting long-running calculations.

## v3.0 Requirements

Requirements for publication and dataset release. Three interlinked deliverables: Radar4Chem dataset with DOI, Nature Scientific Data descriptor, and J. Cheminformatics application note. Each maps to roadmap phases.

### Dataset Preparation

- [ ] **DATA-01**: 650 NWChem calculations organized by molecule/solvent hierarchy with .nw inputs and .out outputs
- [ ] **DATA-02**: Processed data exports: CSV (scaling factors, molecule metadata), JSON (shielding tensors)
- [ ] **DATA-03**: README.md with dataset overview, citation format, usage instructions
- [ ] **DATA-04**: Complete computational provenance documented (NWChem version, COSMO params, compiler, all non-default settings)
- [ ] **DATA-05**: DataCite 4.6 metadata with chemistry extensions (SMILES/InChI/InChIKey for all 50 molecules)
- [ ] **DATA-06**: CC-BY 4.0 license file and version number (v1.0.0)
- [ ] **DATA-07**: Computation time logs extracted from .out files
- [ ] **DATA-08**: Checksums (SHA256) for file integrity verification
- [ ] **DATA-09**: File manifest CSV listing all files with sizes and checksums

### Repository Upload

- [ ] **REPO-01**: Dataset uploaded to Radar4Chem (or Zenodo/Figshare fallback if embargo not supported)
- [ ] **REPO-02**: Reserved DOI obtained and documented
- [ ] **REPO-03**: All DataCite mandatory fields + chemistry-specific extensions completed in repository metadata
- [ ] **REPO-04**: Dataset published immediately (open science -- no embargo)

### Data Descriptor (Nature Scientific Data)

- [ ] **DESC-01**: Background & Summary section (motivation, dataset scope, reuse value)
- [ ] **DESC-02**: Methods section with complete reproducibility detail (DFT functional, basis set, COSMO params with grid settings, SCF criteria, processing scripts)
- [ ] **DESC-03**: Data Records section describing Radar4Chem deposit (DOI, file formats, folder structure)
- [ ] **DESC-04**: Technical Validation section with quantitative metrics (RMSE, R², MAE per solvent/nucleus, outlier analysis, literature comparison)
- [ ] **DESC-05**: Usage Notes section (how to apply scaling factors, software requirements, chemical applicability)
- [ ] **DESC-06**: Data Availability and Code Availability statements with DOIs
- [ ] **DESC-07**: 5-6 publication-quality figures (molecule diversity, workflow, accuracy heatmaps, correlation plots, statistical validation)
- [ ] **DESC-08**: Cross-solvent accuracy analysis showing how performance varies across 13 solvents
- [ ] **DESC-09**: Error analysis by functional group (which molecular features predict poorly)
- [ ] **DESC-10**: Submission package with flattened directory structure ready for journal upload

### Application Note (J. Cheminformatics)

- [ ] **APP-01**: Background section (problem context, why DELTA50 scaling factors needed)
- [ ] **APP-02**: Implementation section (FastAPI architecture, Docker deployment, tech stack)
- [ ] **APP-03**: Operation section (web UI, API endpoints, inputs/outputs with examples)
- [ ] **APP-04**: Use cases section (3-5 practical examples: structure elucidation, diastereomer assignment, solvent selection)
- [ ] **APP-05**: Availability and Requirements section (GitHub URL, Docker, system requirements, license)
- [ ] **APP-06**: Benchmarking against 4-6 state-of-the-art methods (DFT functionals, ML approaches, free tools)
- [ ] **APP-07**: Computational cost comparison (prediction time, scalability, hardware requirements)
- [ ] **APP-08**: 4-6 publication-quality figures (architecture diagram, web UI screenshot, API workflow, use case examples)
- [ ] **APP-09**: Submission package with flattened directory structure ready for journal upload

### AI Writing Validation

- [ ] **VALID-01**: Both manuscripts pass AI writing pattern audit (no formulaic transitions, no superficial prose, no repetitive paragraph structure)
- [ ] **VALID-02**: All citations verified as real (no phantom references, correct DOIs/metadata)
- [ ] **VALID-03**: Technical terminology used correctly and consistently across both papers
- [ ] **VALID-04**: Both manuscripts have distinct authorial voice with domain-specific insights and critical analysis
- [ ] **VALID-05**: Sentence structure varies naturally (no rigid patterns, varied paragraph openings)

## Future Requirements

Deferred to post-acceptance workflow (not in this milestone).

### Post-Acceptance

- **POST-01**: Update dataset metadata with RelatedIdentifier entries when papers published
- **POST-02**: Update application note citation from "submitted" to published reference
- **POST-03**: Verify bidirectional DOI linking between dataset and papers

## Out of Scope

| Feature | Reason |
|---------|--------|
| Peer review tracking | Milestone ends at submission-ready; review is human-driven timeline |
| Example Jupyter notebooks | Nice-to-have but significant scope; defer to post-publication |
| JCAMP-DX format export | IUPAC standard designed for experimental data, unclear applicability to DFT results |
| Parquet format export | CSV/JSON sufficient for dataset; Parquet adds dependency for marginal benefit |
| Preprint posting | Decision for user after submission, not a code/content deliverable |
| New software features | v3.0 is publication-only, no changes to qm-nmr-calc tool |

## Traceability

| Requirement | Phase | Status |
|-------------|-------|--------|
| DATA-01 | Phase 66 | Pending |
| DATA-02 | Phase 66 | Pending |
| DATA-03 | Phase 66 | Pending |
| DATA-04 | Phase 66 | Pending |
| DATA-05 | Phase 66 | Pending |
| DATA-06 | Phase 66 | Pending |
| DATA-07 | Phase 66 | Pending |
| DATA-08 | Phase 66 | Pending |
| DATA-09 | Phase 66 | Pending |
| REPO-01 | Phase 67 | Pending |
| REPO-02 | Phase 67 | Pending |
| REPO-03 | Phase 67 | Pending |
| REPO-04 | Phase 67 | Pending |
| DESC-01 | Phase 68 | Pending |
| DESC-02 | Phase 68 | Pending |
| DESC-03 | Phase 68 | Pending |
| DESC-04 | Phase 68 | Pending |
| DESC-05 | Phase 68 | Pending |
| DESC-06 | Phase 68 | Pending |
| DESC-07 | Phase 68 | Pending |
| DESC-08 | Phase 68 | Pending |
| DESC-09 | Phase 68 | Pending |
| DESC-10 | Phase 68 | Pending |
| APP-01 | Phase 69 | Pending |
| APP-02 | Phase 69 | Pending |
| APP-03 | Phase 69 | Pending |
| APP-04 | Phase 69 | Pending |
| APP-05 | Phase 69 | Pending |
| APP-06 | Phase 69 | Pending |
| APP-07 | Phase 69 | Pending |
| APP-08 | Phase 69 | Pending |
| APP-09 | Phase 69 | Pending |
| VALID-01 | Phase 70 | Pending |
| VALID-02 | Phase 70 | Pending |
| VALID-03 | Phase 70 | Pending |
| VALID-04 | Phase 70 | Pending |
| VALID-05 | Phase 70 | Pending |

**Coverage:**
- v3.0 requirements: 37 total
- Mapped to phases: 37
- Unmapped: 0

**Coverage validation:** 100% ✓

---
*Requirements defined: 2026-02-11*
*Last updated: 2026-02-11 after roadmap creation*
