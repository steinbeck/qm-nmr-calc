# Domain Pitfalls: Dataset Publication and Academic Papers

**Domain:** Computational Chemistry Dataset Publication + Scientific Papers
**Focus:** DELTA50 benchmark dataset release and two papers (Nature Scientific Data + J. Cheminformatics)
**Researched:** 2026-02-11
**Overall confidence:** MEDIUM-HIGH

## Executive Summary

Publishing computational chemistry datasets and writing data descriptor/application note papers involves unique challenges that differ from traditional research publications. The primary pitfalls cluster around: (1) incomplete computational metadata preventing reproducibility, (2) FAIR compliance failures that reduce dataset discoverability, (3) journal-specific formatting violations causing desk rejection, and (4) inadequate benchmarking/comparison undermining software impact.

This research identifies 19 critical pitfalls across four dimensions: dataset preparation, data descriptor writing, application note writing, and cross-publication management. Each pitfall includes detection warning signs, prevention strategies, and phase-specific guidance.

**Key insight:** The most dangerous pitfall is incomplete computational provenance metadata (Pitfall 1.1), which appears fine during manuscript preparation but causes reproducibility failures years later when other researchers attempt to use the dataset.

---

## Critical Pitfalls

Mistakes that cause rejection, major revisions, or long-term reproducibility failures.

### Pitfall 1.1: Incomplete Computational Provenance Metadata

**What goes wrong:** Dataset deposited without complete computational metadata (exact software versions, all calculation parameters, grid settings, convergence criteria). Looks fine initially but becomes irreproducible when researchers attempt to use the dataset.

**Why it happens:** Researchers assume "NWChem 7.2.0, B3LYP/6-311+G(2d,p), COSMO" is sufficient. It's not. Missing: exact NWChem version (7.2.0 vs 7.2.2), COSMO grid parameters (minbem, maxbem, ificos), dielectric constant for each solvent, SCF convergence thresholds, numerical integration grid settings.

**Consequences:**
- Other researchers cannot reproduce calculations (reported RMSE: 0.1-0.2 ppm for 1H differs from dataset)
- Dataset citation rate drops (researchers move to better-documented alternatives)
- FAIR F2 principle violated (insufficient metadata richness)
- Scientific Data reviewers flag as "insufficient technical detail" requiring major revision

**Real-world evidence:** Studies report that "for most papers, there was little to provide help to a researcher willing to reproduce the calculations: molecular and crystal structures discussed were not provided except as snapshots in figures; input files for calculations were not provided" (source: Reproducible Research in Computational Chemistry of Materials, Chemistry of Materials).

**Specific to NWChem/COSMO:** NWChem default COSMO tessellation parameters (minbem=2, maxbem=3, ificos=0) produce very large errors. Researchers using defaults will not reproduce your results if you used finer grids (minbem=3, maxbem=4, ificos=1).

**Prevention:**
1. **Include complete NWChem input files** for representative molecules (all 50 molecules if possible, minimum 5 diverse examples)
2. **Document all non-default parameters:**
   - COSMO grid settings: `cosmo; minbem 3; maxbem 4; ificos 1; end`
   - Dielectric constant for each solvent (not just "COSMO with CHCl3")
   - SCF convergence thresholds if modified from defaults
   - DFT grid settings (if non-default)
3. **Provide complete output files** for representative calculations (allows validation of reported energies/shifts)
4. **Create computational methods checklist:**
   - Software name and EXACT version (include minor/patch versions)
   - Functional and basis set
   - Solvation model with ALL parameters
   - Geometry optimization settings (if geometries were optimized)
   - Property calculation settings (NMR shielding tensor calculation parameters)
   - SCF/optimization convergence criteria
   - Computational environment (compiler versions if relevant for numerical reproducibility)

**Detection warning signs:**
- Methods section draft uses phrases like "standard settings" or "default parameters"
- Collaborator asks "which COSMO grid did we use?" and you have to check old scripts
- Input files not version-controlled or organized systematically
- Different team members give different answers about computational parameters

**Phase-specific guidance:**
- **Phase 1 (Dataset preparation):** Create input file templates with explicit parameters, version control all scripts
- **Phase 2 (Repository preparation):** Include input/output files for 5-10 representative molecules, create README explaining all parameters
- **Phase 3 (Data descriptor writing):** Dedicate full paragraph to each computational method component (basis set, functional, solvation model, convergence)
- **Phase 4 (Review response):** Be prepared for reviewer requests for additional computational details

**References:**
- [Reproducible Research in Computational Chemistry of Materials](https://pubs.acs.org/doi/10.1021/acs.chemmater.7b00799)
- [COSMO Solvation Model - NWChem Documentation](https://nwchemgit.github.io/COSMO-Solvation-Model.html)

---

### Pitfall 1.2: FAIR Principle F1 Violation - Missing Persistent Identifier Links

**What goes wrong:** Dataset receives DOI, papers cite DOI, but metadata doesn't create bidirectional links. Dataset repository doesn't know about papers, papers don't appear in DataCite/Crossref relationships.

**Why it happens:** Authors treat DOI citation like traditional bibliography (one-way reference). Modern scholarly infrastructure requires structured metadata relationships using RelatedIdentifier fields.

**Consequences:**
- Dataset discovery reduced (researchers finding Paper A don't discover Paper B or dataset)
- Citation tracking broken (dataset DOI doesn't show paper citations)
- FAIR F3 principle violated ("Metadata clearly and explicitly include the identifier of the data they describe")
- Impossible to track dataset reuse across publications

**Real-world evidence:** "The RelatedIdentifier property should be used to create links between versions" and "It is recommended to tag the identifier for the dataset with the type attribute to help identify data citations" (source: Crossref Data Citation Deposit Guide).

**Prevention:**
1. **Update dataset metadata when papers accepted:**
   - Add RelatedIdentifier entries in DataCite metadata for both papers
   - Use relationship type: `IsCitedBy` (for papers citing dataset)
   - Include DOI for both Scientific Data descriptor AND J. Cheminformatics application note
2. **Include structured data citations in papers:**
   - Formal citation in References section (not just URL in text)
   - Use DataCite citation format: `Author, Year, Title, Repository, DOI`
   - Example: `Smith, J. (2026). DELTA50 NMR Benchmark Dataset. Radar4Chem. https://doi.org/10.radar/xyz`
3. **Configure Crossref deposits correctly:**
   - Mark dataset citation with `<dataset>` tag in article XML
   - Include relationship type in metadata deposit
4. **Create explicit bidirectional references:**
   - Dataset README links to papers with DOIs
   - Papers' Data Availability Statement links to dataset with DOI
   - Both Crossref and DataCite metadata include RelatedIdentifier fields

**Detection warning signs:**
- Planning to "just mention the DOI in the text"
- No plan to update dataset metadata after papers publish
- Don't know how to add RelatedIdentifier entries to DataCite metadata
- Repository submission form doesn't ask about related publications

**Phase-specific guidance:**
- **Phase 1 (Dataset preparation):** Reserve DOI early, plan for metadata updates
- **Phase 2 (Repository deposit):** Leave RelatedIdentifier fields empty but documented as "to be updated upon publication"
- **Phase 3 (Paper submission):** Include formal dataset citation in References section
- **Phase 4 (After paper acceptance):** UPDATE dataset metadata with paper DOIs within 1 week of publication

**References:**
- [Crossref Relationships Documentation](https://www.crossref.org/documentation/schema-library/markup-guide-metadata-segments/relationships/)
- [DataCite RelatedIdentifier Property](https://support.datacite.org/docs/datacite-metadata-schema)
- [Data Citation Best Practices](https://www.crossref.org/documentation/reference-linking/data-and-software-citation-deposit-guide/)

---

### Pitfall 1.3: FAIR Principle F4 Violation - Insufficient Machine-Readable Metadata

**What goes wrong:** Dataset discoverable by humans (has DOI, listed in Radar4Chem) but not by machines. Missing structured chemistry-specific metadata required by NFDI4Chem/Chem-DCAT-AP extensions.

**Why it happens:** Researchers complete the minimal DataCite metadata (title, creator, description) but skip domain-specific extensions. Radar4Chem accepts deposit, but dataset lacks chemistry-specific fields (reaction types, molecular formula ranges, measurement techniques).

**Consequences:**
- Dataset invisible to chemistry-specific search tools
- Cannot be discovered via ChemSpider, PubChem, chemistry repository aggregators
- Fails NFDI4Chem compliance requirements
- Reduced dataset reuse and citation

**Real-world evidence:** "Domain-specific repositories extract metadata from analytical data files, with descriptive DataCite metadata as a minimum requirement" and "NFDI4Chem and NFDI4Cat demonstrate how to add domain-specific metadata to a dataset with Chem-DCAT-AP, an extension to DCAT-AP" (source: NFDI4Chem Knowledge Base - For Infrastructure Providers).

**NFDI4Chem-specific requirements:**
- Molecular structure representations (SMILES, InChI, InChIKey for all 50 molecules)
- Measurement technique ontology terms (NMR spectroscopy, quantum chemical calculations)
- Chemical class/taxonomy
- Software used (structured, not free text)
- Standards compliance (MIChI - Minimum Information about a Chemistry Investigation)

**Prevention:**
1. **Use Radar4Chem's chemistry-specific metadata fields:**
   - Molecular structures: Provide SMILES strings for all 50 molecules in structured metadata
   - Technique: Tag as "NMR spectroscopy" AND "density functional theory calculations" using ontology terms
   - Software: Structured entry for "NWChem 7.2.0" (not free text)
   - Domain keywords: Use ChEBI/PubChem ontology terms, not just free text
2. **Include chemistry-specific files in standard formats:**
   - SDF file with all 50 molecular structures
   - CSV file with InChI, InChIKey, SMILES for each molecule
   - Machine-readable calculation metadata (JSON or XML)
3. **Follow RADAR Metadata Schema v9.2 requirements:**
   - 10 mandatory fields (core DataCite)
   - Chemistry-specific optional fields (don't skip these for chemistry data!)
   - TS4NFDI standardized terminology integration
4. **Validate metadata completeness:**
   - Use NFDI4Chem FAIR checklist
   - Verify all 50 molecules have InChIKeys (enables cross-database linking)
   - Check that metadata exports as valid Chem-DCAT-AP JSON-LD

**Detection warning signs:**
- Radar4Chem submission form completed in <30 minutes (likely skipped optional chemistry fields)
- No molecular structure files included (only PDF/text descriptions)
- Metadata preview shows only generic DataCite fields, no chemistry-specific extensions
- No ontology terms selected (only free-text keywords)

**Phase-specific guidance:**
- **Phase 1 (Dataset preparation):** Generate SMILES/InChI for all molecules, organize in CSV
- **Phase 2 (Repository deposit):** Complete ALL Radar4Chem chemistry-specific fields, include SDF/CSV structure files
- **Phase 3 (Validation):** Test dataset discoverability via NFDI4Chem search portal (not just Radar4Chem)
- **Phase 4 (Post-publication):** Monitor dataset appearance in chemistry aggregators (ChemSpider, PubChem)

**References:**
- [NFDI4Chem FAIR Data Principles](https://knowledgebase.nfdi4chem.de/knowledge_base/docs/fair/)
- [RADAR4Chem Metadata Schema v9.2](https://knowledgebase.nfdi4chem.de/knowledge_base/docs/radar4chem/)
- [Chem-DCAT-AP Extension](https://knowledgebase.nfdi4chem.de/knowledge_base/docs/metadata/)

---

### Pitfall 2.1: Nature Scientific Data - Insufficient Technical Validation Section

**What goes wrong:** Technical Validation section describes "we checked the data" but doesn't provide quantitative validation metrics, statistical analysis, or independent verification methods. Reviewers flag as insufficient, requiring major revision.

**Why it happens:** Authors confuse "validation" with "quality control." Reviewers want to know HOW you verified data quality (statistical tests, cross-validation, comparison with independent datasets), not just THAT you checked it.

**Consequences:**
- Manuscript rejected or major revisions required
- Extends publication timeline by 3-6 months
- Reviewers may request additional validation experiments

**Real-world evidence:** "Referees are asked to assess whether the methods and any data processing steps are described in sufficient detail to allow others to reproduce these steps" (source: Nature Scientific Data - For Referees). Scientific Data reviewers specifically evaluate whether Technical Validation convincingly supports data quality.

**What reviewers expect in Technical Validation:**
1. **Quantitative metrics:**
   - RMSE between DFT-calculated and experimental shifts (for your 650 calculations)
   - Statistical distribution analysis (mean absolute error, standard deviation, outlier detection)
   - Correlation coefficients (R²) for each nucleus (1H vs 13C)
2. **Cross-validation:**
   - Compare subset of molecules with literature experimental values
   - Independent validation set (molecules not used in OLS scaling factor determination)
   - Comparison with alternative DFT methods (shows B3LYP/6-311+G(2d,p) is appropriate)
3. **Outlier analysis:**
   - Identify and explain outliers (molecules where prediction fails)
   - Statistical tests for outlier detection
   - Chemical explanation for outliers (e.g., hydrogen bonding, conformational averaging)
4. **Reproducibility verification:**
   - Independent recalculation of 5-10 molecules by different team member
   - Verification that same results obtained from provided input files

**Prevention:**
1. **Perform quantitative validation analysis:**
   - Calculate RMSE, MAE, R² for full dataset
   - Stratify by solvent (do some solvents have higher errors?)
   - Stratify by chemical shift range (aromatic vs aliphatic)
   - Identify and analyze outliers (>2σ deviations)
2. **Include statistical visualizations:**
   - Predicted vs experimental correlation plots
   - Residual distribution histograms
   - Per-solvent error box plots
3. **Compare with literature benchmarks:**
   - How does DELTA50 accuracy compare to previous NMR prediction benchmarks?
   - Cite specific benchmarking studies (e.g., "Our RMSE of X ppm compares favorably with Y ppm reported by [citation]")
4. **Demonstrate reproducibility:**
   - Different team member recalculates 10 molecules, reports identical results
   - Document as: "Independent recalculation of 10 randomly selected molecules by a second researcher yielded identical chemical shifts to within 0.01 ppm"

**Specific to DELTA50:**
- Report separate validation for 1H vs 13C (different accuracy expectations)
- Analyze OLS scaling factor stability across solvents
- Validate that molecules chosen avoid concentration-dependent aggregation (critical for experimental comparison)
- Compare with existing NMR prediction benchmarks (NMRShiftDB, other DFT studies)

**Detection warning signs:**
- Technical Validation section <500 words (likely insufficient)
- No quantitative metrics (just qualitative descriptions like "good agreement")
- No figures/tables in Technical Validation section
- No comparison with literature or alternative methods
- Phrases like "data were carefully checked" without describing HOW

**Phase-specific guidance:**
- **Phase 1 (Data generation):** Calculate validation metrics as data is generated (don't wait until manuscript writing)
- **Phase 2 (Dataset finalization):** Perform independent recalculation, statistical outlier analysis
- **Phase 3 (Data descriptor writing):** Dedicate 2-3 pages to Technical Validation with tables/figures
- **Phase 4 (Pre-submission):** Have domain expert review Technical Validation section specifically

**References:**
- [Nature Scientific Data Submission Guidelines](https://www.nature.com/sdata/submission-guidelines)
- [Nature Scientific Data - For Referees](https://www.nature.com/sdata/policies/for-referees)
- [NMR Prediction Benchmarking Best Practices](https://www.mdpi.com/1420-3049/28/6/2449)

---

### Pitfall 2.2: Nature Scientific Data - Methods Section Lacks Reproducibility Detail

**What goes wrong:** Methods section describes computational workflow but omits details needed for exact reproduction. Reviewers cannot verify that described methods would produce the deposited dataset.

**Why it happens:** Authors describe methods at "publication level" (suitable for journal article) instead of "reproduction level" (suitable for data descriptor). Scientific Data requires more detail than typical papers.

**Consequences:**
- Reviewers flag as "methods not described in sufficient detail to allow reproduction"
- Major revisions required
- May need to regenerate data if undocumented steps cannot be verified

**Real-world evidence:** "Data Descriptors contain details of how datasets were created (Methods), what they contain (Data Records), and how they were checked and validated (Technical Validation). The Methods section should provide sufficient technical detail to ensure that others can understand and reproduce the data collection and processing steps" (source: Scientific Data Submission Guidelines).

**Missing details that commonly cause rejection:**
1. **Computational workflow:**
   - How were input files generated? (manual? scripted? template-based?)
   - Processing pipeline order (geometry optimization → NMR calculation? or direct NMR on provided geometry?)
   - Data extraction method (how were chemical shifts extracted from NWChem output files?)
   - Quality control steps (how were failed calculations detected and handled?)
2. **Data processing:**
   - OLS scaling factor calculation methodology (exact equations used)
   - Reference compound handling (TMS: calculated or experimental reference?)
   - Averaging procedures (if multiple conformers calculated, how averaged?)
3. **Software environment:**
   - Operating system and version
   - Compiler used to build NWChem (affects numerical results)
   - Parallel execution details (single-core? MPI? affects reproducibility)
4. **File formats and organization:**
   - How are files named and organized?
   - What file format stores final chemical shifts? (CSV? JSON? database?)
   - Data processing scripts (provided? which language?)

**Prevention:**
1. **Document complete computational workflow:**
   - Step-by-step protocol from molecular structure → final chemical shift
   - Include script/code snippets for critical steps
   - Specify order of operations
   - Describe quality control at each step
2. **Provide all processing scripts:**
   - Deposit scripts in dataset repository (not just in paper SI)
   - Include: input file generation, output parsing, statistical analysis, OLS scaling
   - Document script dependencies (Python version, required libraries)
3. **Specify software environment completely:**
   - NWChem version: 7.2.0 (include build date if multiple 7.2.0 releases exist)
   - Compiler: GCC 11.2.0 (or Intel Fortran version)
   - MPI implementation: OpenMPI 4.1.1 (if parallel calculations)
   - Operating system: Ubuntu 22.04 LTS
4. **Describe OLS scaling procedure explicitly:**
   - Reference dataset used for OLS fitting (is it separate from DELTA50 or subset?)
   - Statistical method (ordinary least squares regression)
   - Equations: σ_scaled = (σ_raw - intercept) / slope
   - Per-solvent or global scaling? (separate OLS parameters for each solvent?)
   - Report OLS parameters (slope, intercept, R²) in paper

**Specific to DELTA50:**
- Describe how 50 molecules were selected (inclusion/exclusion criteria)
- Explain why 13 solvents chosen (covers polarity range? common lab solvents?)
- Document how experimental reference data was obtained
- Explain aggregation-avoiding criteria

**Detection warning signs:**
- Methods section <1500 words (likely too brief for data descriptor)
- Workflow described in paragraph form (no numbered steps or flowchart)
- No code/scripts mentioned
- References computational details to "standard protocols" without explicit description
- Collaborator asks "how did we do X?" and you have to check old lab notebooks

**Phase-specific guidance:**
- **Phase 1 (Data generation):** Document workflow in lab notebook AS YOU GO (don't reconstruct from memory)
- **Phase 2 (Workflow documentation):** Create flowchart, write step-by-step protocol, test protocol on new molecule
- **Phase 3 (Data descriptor writing):** Dedicate 3-4 pages to Methods with workflow diagram
- **Phase 4 (Pre-submission):** Have researcher unfamiliar with project read Methods and identify ambiguities

**References:**
- [Nature Scientific Data Submission Guidelines - Methods Section](https://www.nature.com/sdata/submission-guidelines)
- [Editorial: Method and Data Sharing and Reproducibility](https://pubs.acs.org/doi/10.1021/acs.jcim.0c01389)

---

### Pitfall 3.1: J. Cheminformatics Application Note - Inadequate Benchmarking

**What goes wrong:** Application note presents NMR prediction tool but compares against limited or outdated baselines. Reviewers reject for insufficient demonstration of advantage over existing tools.

**Why it happens:** Authors benchmark against 1-2 convenient comparisons instead of comprehensive state-of-the-art survey. Common mistake: comparing against older tools or limited benchmark sets.

**Consequences:**
- Manuscript rejected: "does not demonstrate significant advance over previously published software"
- Extends publication timeline by 6+ months (need to perform additional benchmarks)
- Reduced impact: readers assume tool not competitive if benchmarking incomplete

**Real-world evidence:** "Software articles should describe a tool likely to be of broad utility that represents a significant advance over previously published software, usually demonstrated by direct comparison with available related software" (source: BMC Bioinformatics Software Article Guidelines). NMR benchmarking studies show "conflicting results—for example, BMK was found to be the best performing functional by one study but the worst by another" highlighting need for standardized benchmarks (source: Computational NMR Prediction Review).

**Inadequate benchmarking patterns:**
1. **Cherry-picked comparisons:** Only comparing against tools that your method outperforms
2. **Outdated baselines:** Comparing against 5-10 year old methods while ignoring recent ML approaches
3. **Single benchmark dataset:** Only testing on your DELTA50 dataset (circular validation)
4. **Limited chemical diversity:** Testing only on small organic molecules, ignoring heterocycles, organometallics
5. **No statistical significance testing:** Reporting point estimates without confidence intervals

**Expected by reviewers:**
1. **Comparison with state-of-the-art:**
   - DFT methods: Multiple functionals (B3LYP, ωB97X-D, M06-2X, etc.)
   - ML methods: Recently published neural network approaches (if claiming novelty)
   - Commercial tools: ChemDraw NMR predictor, ACD/Labs (if applicable)
   - Free tools: nmrshiftdb2, HOSE code methods
2. **Multiple benchmark datasets:**
   - DELTA50 (your own)
   - Independent literature datasets (e.g., NMRShiftDB subsets)
   - External validation sets from other publications
   - Diverse chemical classes (aromatics, aliphatics, heterocycles)
3. **Statistical rigor:**
   - RMSE with 95% confidence intervals
   - Statistical significance tests (is improvement significant? paired t-test)
   - Per-molecule error analysis (not just aggregate statistics)
   - Outlier analysis (which molecules does your method fail on? do competitors also fail?)
4. **Computational cost comparison:**
   - Prediction time per molecule
   - Scalability to large molecules
   - Hardware requirements

**Prevention:**
1. **Survey existing tools comprehensively:**
   - Literature review: NMR prediction papers from last 5 years
   - Identify current state-of-the-art (highest accuracy, most cited recent methods)
   - Include both DFT and ML approaches
2. **Design benchmark matrix:**
   - Rows: Methods (your tool + 4-6 competitors)
   - Columns: Metrics (RMSE 1H, RMSE 13C, time, cost)
   - Test on multiple datasets (your DELTA50 + 2-3 external)
3. **Perform fair comparisons:**
   - Same molecules for all methods
   - Document when competitor cannot handle certain molecule types
   - Report both successes and failures
4. **Include computational cost analysis:**
   - Your tool: X CPU-hours per molecule
   - DFT baseline: Y CPU-hours per molecule
   - ML alternative: Z seconds per molecule
   - Accuracy-cost tradeoff discussion

**Specific to DELTA50/NMR tool:**
- Compare against: DFT (multiple functionals), HOSE codes, ML (graph neural networks)
- Test on: DELTA50, NMRShiftDB validation set, literature experimental datasets
- Metrics: RMSE (1H, 13C separately), DP4+ probability accuracy, outlier rate
- Chemical diversity: Include challenging cases (quaternary carbons, heteroatoms, aromatics)

**Detection warning signs:**
- Benchmark table has <3 comparison methods
- No comparison with recently published ML approaches
- Only testing on your own dataset
- No discussion of when your method fails
- No computational cost comparison

**Phase-specific guidance:**
- **Phase 1 (Tool development):** Track state-of-the-art continuously, maintain benchmark dataset collection
- **Phase 2 (Benchmarking):** Design comprehensive comparison, run all methods on same datasets
- **Phase 3 (Application note writing):** Dedicate full section to benchmarking with comparison tables
- **Phase 4 (Pre-submission):** Have cheminformatics expert review benchmark design for completeness

**References:**
- [Journal of Cheminformatics Software Guidelines](https://jcheminf.biomedcentral.com/submission-guidelines/preparing-your-manuscript/software)
- [MoleculeNet Benchmarking Framework](https://pubs.rsc.org/en/content/articlehtml/2018/sc/c7sc02664a)
- [NMR Prediction Benchmarking](https://www.mdpi.com/1420-3049/28/6/2449)

---

### Pitfall 3.2: J. Cheminformatics Application Note - Underselling Tool Impact and Novelty

**What goes wrong:** Application note presents tool functionality but fails to articulate WHY it matters, WHO should use it, and WHAT problems it solves better than existing tools. Reviewers assess as "incremental contribution, limited impact."

**Why it happens:** Authors focus on technical implementation details ("we built a web interface with React") instead of scientific value proposition ("enables NMR structure elucidation for non-expert users who cannot run DFT calculations"). Common in software papers written by developers rather than domain scientists.

**Consequences:**
- Desk rejection or rejection after review ("does not demonstrate significant advance")
- Low citation rate even if published (users don't understand when to use tool)
- Missed opportunity for high-impact publication

**Real-world evidence:** Common pattern in software underselling: "developers overlook the contribution when they work on technical pieces but do not explain the business objectives or outcomes of the projects" (source: general software documentation best practices). This applies to academic software: explain scientific objectives, not just technical features.

**Underselling patterns:**
1. **Feature-focused instead of problem-focused:**
   - Bad: "Tool provides web interface for NMR prediction"
   - Good: "Tool democratizes NMR prediction for synthetic chemists without computational expertise"
2. **Missing use case scenarios:**
   - Bad: "Supports 50 molecules and 13 solvents"
   - Good: "Enables rapid diastereomer discrimination during total synthesis campaigns"
3. **No impact quantification:**
   - Bad: "Predictions are accurate"
   - Good: "Predictions achieve 95% correct structural assignment vs 70% for existing free tools"
4. **Ignoring unique advantages:**
   - Bad: Generic description of features
   - Good: Highlighting what ONLY your tool can do (e.g., DELTA50 is first benchmark with explicit solvent coverage)

**Prevention:**
1. **Lead with scientific value proposition:**
   - First paragraph: What scientific problem does this solve?
   - Who is the target user? (computational chemists? synthetic chemists? spectroscopists?)
   - What can users do with your tool that they couldn't do before?
2. **Include concrete use cases:**
   - Use case 1: Structure elucidation (unknown compound → candidate structures → NMR prediction → best match)
   - Use case 2: Diastereomer assignment (synthetic product → predict NMR for each diastereomer → match to experiment)
   - Use case 3: Solvent selection (predict NMR in 13 solvents → choose solvent for maximum chemical shift separation)
3. **Quantify impact with metrics:**
   - Accessibility: "Reduces NMR prediction from 48 CPU-hours to 30 seconds"
   - Accuracy: "Achieves X% correct assignment vs Y% for free alternatives"
   - Usability: "Requires only molecular structure input, no computational expertise"
   - Coverage: "First benchmark with explicit solvent effects"
4. **Compare against user's current workflow:**
   - Without tool: "User must run DFT calculation (48h, requires HPC access, expertise)"
   - With tool: "User uploads structure, receives prediction in 30s via web interface"
5. **Highlight novelty explicitly:**
   - Technical novelty: "First tool to combine X with Y"
   - Scientific novelty: "First benchmark to include explicit solvent effects"
   - Usability novelty: "First web-accessible NMR predictor with DFT accuracy"

**Specific to DELTA50 application note:**
- **Problem:** Existing NMR databases lack solvent-explicit predictions, limiting utility for synthetic chemists
- **Solution:** DELTA50 provides solvent-specific predictions (13 solvents) with quantified accuracy
- **Impact:** Enables solvent selection for structure elucidation, diastereomer assignment in specific solvents
- **Novelty:** First benchmark systematically covering solvent effects with DFT calculations
- **Users:** Synthetic chemists (structure assignment), computational chemists (method benchmarking), spectroscopists (assignment validation)

**Detection warning signs:**
- Introduction focuses on implementation details instead of scientific problem
- No user personas or use cases described
- Comparison table focuses on features (# of solvents) instead of outcomes (accuracy improvement)
- No discussion of who should use tool and when
- Writing sounds like technical documentation instead of scientific paper

**Phase-specific guidance:**
- **Phase 1 (Tool design):** Define target users and use cases BEFORE building features
- **Phase 2 (Validation):** Quantify impact metrics (time savings, accuracy improvement, user success rate)
- **Phase 3 (Application note writing):** Lead with problem/solution, dedicate section to use cases
- **Phase 4 (Pre-submission):** Have non-expert chemist read and explain back what problem tool solves (if they can't, rewrite)

**References:**
- [Journal of Cheminformatics Software Article Guidelines](https://jcheminf.biomedcentral.com/submission-guidelines/preparing-your-manuscript/software)
- [Best Practices in Machine Learning for Chemistry](https://www.nature.com/articles/s41557-021-00716-z)

---

## Moderate Pitfalls

Mistakes that cause delays, revisions, or reduce impact but don't typically cause rejection.

### Pitfall 1.4: Dataset Versioning Strategy Not Planned

**What goes wrong:** Dataset published with DOI, then errors discovered or new molecules added. No versioning plan, resulting in either: (1) silent updates that break reproducibility, or (2) inability to update without confusing citation landscape.

**Why it happens:** Authors assume dataset is "done" at publication. Real-world: errors discovered, new molecules calculated, methods refined.

**Consequences:**
- Reproducibility failures: "I downloaded DELTA50 in 2027 and got different results than 2026 paper"
- Citation confusion: Some papers cite v1, others cite v2, results not comparable
- Violates FAIR principles (versioning enables repeatability)

**Real-world evidence:** "The Iris dataset from the UCI ML Repository was discovered to have multiple widely-publicized versions with differing measurements for certain observations, resulting in non-comparable reported performances of classification models across published papers" (source: Versioning Data Is About More than Revisions, Data Science Journal).

**Prevention:**
1. **Define versioning policy BEFORE publication:**
   - Major version (v2.0): New molecules added, methodology changed
   - Minor version (v1.1): Errors corrected, no molecule additions
   - Patch version (v1.0.1): Typos in metadata, no data changes
2. **Use DataCite versioning:**
   - New major version gets new DOI
   - Use RelatedIdentifier with `IsNewVersionOf` relationship
   - Update metadata of old version with `IsPreviousVersionOf`
3. **Document version history:**
   - CHANGELOG.md in repository
   - Version number in all dataset files
   - Date stamp in metadata
4. **Plan for corrections:**
   - If error discovered: Issue new minor version within 1 month
   - Announce on project website/mailing list
   - Update papers' errata if errors affect published conclusions

**Detection warning signs:**
- No versioning discussion in dataset documentation
- Dataset files not version-stamped
- No plan for handling discovered errors

**Phase-specific guidance:**
- **Phase 1 (Dataset preparation):** Include version number in all files (v1.0.0)
- **Phase 2 (Repository deposit):** Document versioning policy in README
- **Phase 3 (Post-publication):** Monitor for errors, have rapid correction process

**References:**
- [DataCite Versioning Guidelines](https://support.datacite.org/docs/versioning)
- [Benchmark Data Versioning Best Practices](https://arxiv.org/html/2410.24100v1)

---

### Pitfall 1.5: Missing Computational Environment Documentation

**What goes wrong:** Dataset includes input files but not environment details. Researchers attempt reproduction but get slightly different results due to compiler differences, library versions, or numerical precision variations.

**Why it happens:** Authors assume "same software version = same results." Not true for compiled codes like NWChem where compiler, BLAS/LAPACK libraries, and MPI implementation affect numerical results.

**Consequences:**
- Numerical reproducibility failures (differences in 3rd-4th decimal place)
- Users report "bug" in dataset (actually environment differences)
- Reduces dataset trustworthiness

**Prevention:**
1. **Document complete computational environment:**
   - OS: Ubuntu 22.04 LTS
   - Compiler: GCC 11.2.0 or Intel Fortran 2021.4
   - MPI: OpenMPI 4.1.1
   - BLAS/LAPACK: OpenBLAS 0.3.20 or Intel MKL 2021.4
   - Python (for analysis scripts): 3.9.7
2. **Provide environment reproduction tools:**
   - Dockerfile or Singularity container definition
   - Conda environment.yml
   - Module files for HPC systems
3. **Document numerical precision expectations:**
   - "Chemical shifts reproducible to ±0.01 ppm with specified environment"
   - "Different compilers may introduce variations <0.05 ppm"
4. **Test cross-environment reproducibility:**
   - Run same calculation on different systems
   - Document observed variations
   - Report in Technical Validation section

**Detection warning signs:**
- Dataset documentation mentions software but not compiler/libraries
- No container or environment file provided
- Haven't tested reproduction on different systems

**Phase-specific guidance:**
- **Phase 1 (Data generation):** Record environment details as calculations run
- **Phase 2 (Documentation):** Create Dockerfile, test on clean system
- **Phase 3 (Publication):** Include environment documentation in Technical Validation

**References:**
- [Reproducible Research in Computational Chemistry](https://pubs.acs.org/doi/10.1021/acs.chemmater.7b00799)
- [Achieving Reproducibility in Molecular Dynamics](https://pubs.acs.org/doi/10.1021/acs.jced.5c00010)

---

### Pitfall 2.3: Data Availability Statement Incomplete or Incorrect Format

**What goes wrong:** Data Availability Statement (DAS) is vague ("data available upon request") or missing required elements (DOI, access conditions). Nature journal policies require specific DAS format.

**Why it happens:** Authors copy DAS from other papers without checking current journal requirements. Nature portfolio updated DAS requirements in recent years.

**Consequences:**
- Manuscript returned for revision (minor but delays publication)
- Dataset not discoverable if DAS lacks DOI
- Violates journal data policy

**Real-world evidence:** "The Data Availability Statement must list the name of the repository or repositories as well as digital object identifiers (DOIs), accession numbers or codes, or other persistent identifiers for all relevant data" (source: Springer Nature Data Availability Statements Policy).

**Required elements in DAS:**
1. **Repository name:** "Radar4Chem" (not generic "publicly available")
2. **Persistent identifier:** Full DOI (https://doi.org/10.radar/xyz), not just URL
3. **Access conditions:** "Publicly accessible under CC-BY 4.0 license"
4. **Date deposited:** (if relevant for embargo)

**Correct DAS format:**
```
Data Availability
The DELTA50 NMR benchmark dataset is publicly available in the Radar4Chem repository (https://doi.org/10.radar/xyz) under a Creative Commons Attribution 4.0 International (CC-BY 4.0) license. The dataset includes NWChem input files, output files, and calculated NMR chemical shifts for 50 molecules in 13 solvents. Analysis scripts are available at GitHub (https://github.com/user/delta50-analysis).
```

**Incorrect patterns:**
- "Data available upon reasonable request" (violates open data policy)
- "Data in supplementary materials" (should be in proper repository)
- Only URL, no DOI (not persistent identifier)
- No license information (ambiguous reuse rights)

**Prevention:**
1. **Check journal's DAS policy:** Nature portfolio has specific requirements
2. **Include all required elements:** Repository name, DOI, license, access conditions
3. **Deposit data BEFORE manuscript submission:** Need DOI for DAS
4. **Verify DOI resolves:** Test DOI link before submission

**Detection warning signs:**
- Writing DAS before dataset deposited (no DOI yet)
- Copying DAS from old papers without checking current format
- DAS <50 words (likely too brief)

**Phase-specific guidance:**
- **Phase 2 (Repository deposit):** Get DOI from Radar4Chem
- **Phase 3 (Manuscript writing):** Draft DAS with all required elements
- **Phase 4 (Pre-submission):** Verify DOI resolves, test from clean browser

**References:**
- [Springer Nature Data Availability Statements](https://www.springernature.com/gp/authors/research-data-policy/data-availability-statements)
- [Nature Portfolio Reporting Standards](https://www.nature.com/nature-portfolio/editorial-policies/reporting-standards)

---

### Pitfall 2.4: Missing Software/Code Availability Section

**What goes wrong:** Data descriptor includes scripts for data processing but doesn't provide code availability statement or repository. Reviewers flag as incomplete, reproducibility reduced.

**Why it happens:** Authors focus on data, forget that processing scripts are essential for reproducibility. Nature journals increasingly require code availability.

**Consequences:**
- Reviewers request code deposit (delays publication)
- Reproducibility compromised without scripts
- Violates Scientific Data's reproducibility standards

**Prevention:**
1. **Deposit all code in persistent repository:**
   - GitHub + Zenodo DOI (preferred for long-term archival)
   - Include: input file generation scripts, output parsing, statistical analysis
2. **Provide Code Availability statement:**
   ```
   Code Availability
   Scripts for NWChem input file generation, output parsing, and statistical analysis are available at https://github.com/user/delta50-scripts with archived version at Zenodo (https://doi.org/10.5281/zenodo.xyz). Scripts require Python 3.9+ with dependencies listed in requirements.txt.
   ```
3. **Document code:**
   - README with usage examples
   - requirements.txt or environment.yml
   - Example runs with test data
4. **Test code on clean system:**
   - Have collaborator run scripts without your help
   - Verify all dependencies documented

**Detection warning signs:**
- Methods mention "scripts" but no repository linked
- Code only in personal GitHub (not archived)
- No installation/usage documentation

**Phase-specific guidance:**
- **Phase 1 (Data generation):** Version control scripts from day 1
- **Phase 2 (Code preparation):** Clean up scripts, add documentation, create Zenodo archive
- **Phase 3 (Manuscript writing):** Include Code Availability statement

**References:**
- [Nature Portfolio Code Availability](https://www.nature.com/nature-portfolio/editorial-policies/reporting-standards)
- [Scientific Data Submission Guidelines](https://www.nature.com/sdata/submission-guidelines)

---

### Pitfall 3.3: Application Note Exceeds Page Limits or Format Violations

**What goes wrong:** Application note exceeds 2-page limit (common for Journal of Cheminformatics application notes) or violates formatting requirements. Requires shortening/reformatting before review.

**Why it happens:** Authors write standard-length paper, then realize application note has strict limits. Cutting content is difficult after writing.

**Consequences:**
- Desk rejection (if significantly over limit)
- Delays submission by 1-2 weeks while reformatting
- Important content cut to meet limits (reduces quality)

**Prevention:**
1. **Check journal requirements BEFORE writing:**
   - Page limit (typically 2 pages for application notes)
   - Figure limit (often 1-2 figures)
   - Word count (if specified)
   - Section structure (some journals require specific sections)
2. **Write to length from start:**
   - Outline with word budgets per section
   - Introduction: 300 words
   - Implementation: 400 words
   - Validation: 500 words
   - Conclusion: 200 words
3. **Move details to supplementary materials:**
   - Full benchmark tables → supplement
   - Detailed methodology → supplement
   - Main text: summary only
4. **Use figures efficiently:**
   - Multi-panel figures (pack information densely)
   - Reduce white space
   - Supplement for additional figures

**Journal of Cheminformatics application note requirements:**
- **Length:** Up to 2 pages (approximately 1000-1500 words)
- **Sections:** No required structure, but typically Introduction, Implementation, Validation
- **Figures:** 1-2 recommended
- **Availability:** Software must be freely available to non-commercial users for 2 years minimum

**Detection warning signs:**
- Starting to write before checking page limits
- Draft is 4+ pages (will need major cutting)
- No plan for supplementary materials

**Phase-specific guidance:**
- **Phase 1 (Pre-writing):** Read journal's application note requirements in detail
- **Phase 2 (Writing):** Write concise version first, expand if space available
- **Phase 3 (Pre-submission):** Check page count in journal template (not Word default)

**References:**
- [Journal of Cheminformatics Software Article Guidelines](https://jcheminf.biomedcentral.com/submission-guidelines/preparing-your-manuscript/software)
- [BMC Bioinformatics Application Note Guidelines](https://bmcbioinformatics.biomedcentral.com/submission-guidelines/preparing-your-manuscript/software-article)

---

### Pitfall 4.1: Dataset DOI and Paper DOIs Not Cross-Referenced in Time

**What goes wrong:** Dataset published with DOI, papers cite dataset, but metadata not updated bidirectionally before papers publish. Results in broken citation network.

**Why it happens:** Timing mismatch. Dataset DOI obtained during deposit (month 1), papers publish later (month 6-12). Authors forget to update dataset metadata when papers publish.

**Consequences:**
- Citation tracking broken (DataCite doesn't show papers citing dataset)
- Reduced discoverability (researchers finding dataset don't find papers, vice versa)
- Metrics undercount impact

**Prevention:**
1. **Create timeline with metadata update milestones:**
   - Month 1: Dataset deposit, DOI reserved
   - Month 3: Papers submitted (include dataset DOI in manuscript)
   - Month 6-12: Paper 1 accepted → UPDATE dataset metadata within 1 week
   - Month 6-12: Paper 2 accepted → UPDATE dataset metadata within 1 week
2. **Set calendar reminders:**
   - When each paper accepted: "Update Radar4Chem metadata"
   - 1 month after each paper publishes: "Verify Crossref/DataCite linkage"
3. **Document metadata update procedure:**
   - How to update Radar4Chem metadata (log in, edit metadata, add RelatedIdentifier)
   - Required fields: DOI of paper, relationship type (IsCitedBy), title of paper
4. **Verify bidirectional linking:**
   - Check DataCite page for dataset: Does it list both papers as citations?
   - Check Crossref page for papers: Do they show dataset as related resource?

**Detection warning signs:**
- No documented procedure for updating dataset metadata
- No calendar reminders set for metadata updates
- Dataset metadata hasn't been touched since initial deposit (months ago)

**Phase-specific guidance:**
- **Phase 1 (Dataset deposit):** Note in project plan "update metadata when papers publish"
- **Phase 2 (Paper submission):** Confirm dataset metadata allows updates
- **Phase 3 (Paper acceptance):** Update dataset metadata within 1 week
- **Phase 4 (Verification):** Check Crossref/DataCite pages 1 month after publication

**References:**
- [Crossref Relationships Documentation](https://www.crossref.org/documentation/schema-library/markup-guide-metadata-segments/relationships/)
- [DataCite RelatedIdentifier](https://support.datacite.org/docs/datacite-metadata-schema)

---

### Pitfall 4.2: Inconsistent Dataset Citation Format Across Papers

**What goes wrong:** Scientific Data paper cites dataset with one format, J. Cheminformatics paper uses different format. Creates confusion, reduces citation tracking effectiveness.

**Why it happens:** Authors write papers months apart, don't coordinate citation format. Each paper follows different journal's reference style.

**Consequences:**
- Citation parsers may not recognize as same dataset
- Looks unprofessional (inconsistent across related publications)
- Reduces citation tracking effectiveness

**Prevention:**
1. **Establish canonical dataset citation format:**
   - Follow DataCite recommended format
   - Example: `LastName, FirstName (Year). Title. Repository. DOI`
   - Specific: `Smith, J.; Jones, M. (2026). DELTA50: NMR Benchmark Dataset. Radar4Chem. https://doi.org/10.radar/delta50`
2. **Use identical citation in both papers:**
   - Copy exact citation between papers
   - Only difference: journal reference style (numbered vs author-date)
   - Content of citation identical
3. **Verify during manuscript preparation:**
   - Before submission: Compare dataset citations in both papers
   - Ensure DOI, title, authors, year identical

**Specific to DELTA50 publications:**
- Same author list on dataset as papers? (or different?)
- If different: Decide on dataset author list early
- Use dataset DOI in both papers (not repository URL)
- Include dataset version number if applicable

**Detection warning signs:**
- Writing papers separately without cross-checking
- Different collaborators writing different papers (communication gap)
- Dataset citation from memory instead of copying

**Phase-specific guidance:**
- **Phase 1 (Dataset deposit):** Establish canonical citation format
- **Phase 2 (Paper writing):** Use identical dataset citation in both manuscripts
- **Phase 3 (Pre-submission):** Compare dataset citations across papers

---

## Minor Pitfalls

Mistakes that cause annoyance or suboptimal outcomes but are easily fixable.

### Pitfall 1.6: Unclear File Naming and Organization

**What goes wrong:** Dataset files organized idiosyncratically (e.g., `mol_001_calc_final_v2.out`). Users cannot find data for specific molecules or solvents without extensive documentation reading.

**Why it happens:** Files named during calculation workflow, not reorganized for publication. Logical to original researcher, confusing to others.

**Prevention:**
1. **Use systematic file naming convention:**
   - Pattern: `{molecule_ID}_{solvent}_{calc_type}.{ext}`
   - Example: `MOL001_CDCl3_nmr.out` (clear: molecule 1, chloroform, NMR calculation)
2. **Organize in directory structure:**
   ```
   DELTA50/
     input_files/
       MOL001_CDCl3.nw
       MOL001_DMSO.nw
       ...
     output_files/
       MOL001_CDCl3.out
       MOL001_DMSO.out
       ...
     results/
       chemical_shifts.csv
       molecular_structures.sdf
     scripts/
       parse_output.py
       ...
     README.md
   ```
3. **Provide file manifest:**
   - CSV listing all files with descriptions
   - Columns: filename, molecule_ID, solvent, calculation_type, description
4. **Document naming convention in README:**
   - Explain pattern
   - Provide examples
   - Include lookup table (molecule_ID → IUPAC name)

**Detection warning signs:**
- File names include "final", "v2", "revised" (indicates ad-hoc naming)
- Cannot identify molecule/solvent from filename alone
- Files scattered in flat directory

**Phase-specific guidance:**
- **Phase 1 (Data generation):** Use systematic naming from start
- **Phase 2 (Publication prep):** Rename/reorganize files if needed
- **Phase 3 (Documentation):** Create file manifest, document convention

---

### Pitfall 2.5: Figures Not Optimized for Data Descriptor Format

**What goes wrong:** Figures designed for presentation/thesis, not publication. Text too small, colors not colorblind-friendly, resolution too low.

**Why it happens:** Authors reuse existing figures instead of redesigning for journal.

**Prevention:**
1. **Follow Nature Scientific Data figure guidelines:**
   - Minimum 300 DPI resolution
   - Size figures for column width (85 mm) or full width (170 mm)
   - Font size ≥7 pt in final figure
   - Colorblind-friendly color schemes (use ColorBrewer)
2. **Design for black & white printing:**
   - Use different line styles (solid, dashed, dotted) not just colors
   - Symbols distinguishable in grayscale
3. **High information density:**
   - Multi-panel figures (a, b, c subpanels)
   - Minimize white space
   - Clear, informative captions

**Detection warning signs:**
- Figures copied from presentations
- Text illegible when printed at column width
- Red/green color scheme (colorblind unfriendly)

**Phase-specific guidance:**
- **Phase 3 (Figure creation):** Design in journal template dimensions from start
- **Phase 4 (Pre-submission):** Print figures in grayscale to check readability

---

### Pitfall 3.4: Software Name/Branding Not Memorable or Searchable

**What goes wrong:** Tool named generically ("NMR Predictor") or as awkward acronym ("QCNMRPREDSFV2"). Users cannot find tool via search, cannot remember name.

**Why it happens:** Authors choose first name that comes to mind, don't consider searchability/memorability.

**Prevention:**
1. **Choose distinctive, searchable name:**
   - Good: "DELTA50" (distinctive, short, searchable)
   - Bad: "NMR Tool" (generic, unsearchable)
   - Bad: "QCNMRPRED" (unmemorable acronym)
2. **Check name availability:**
   - Google search: Does name return unrelated results?
   - Domain name available? (even if not using, signals availability)
   - PyPI/CRAN/NPM: Name available? (if distributing as package)
3. **Trademark check:**
   - Commercial tools with similar names? (avoid confusion)
4. **Make name meaningful:**
   - DELTA50: "Dataset of Experimental and Linear-scaled Theoretical Assignments, 50 molecules"
   - Not required, but helpful if natural

**Detection warning signs:**
- Tool name is generic noun ("Calculator", "Predictor")
- Acronym requires 5+ word expansion
- Google search for name returns unrelated content on first page

**Phase-specific guidance:**
- **Phase 1 (Early development):** Choose name, register domain/GitHub org
- **Phase 2 (Pre-publication):** Verify name uniqueness, searchability

---

### Pitfall 4.3: No Plan for Dataset Long-Term Maintenance

**What goes wrong:** Dataset published, then authors move institutions/leave academia. No one maintains dataset if issues discovered. Repository may delete after inactivity.

**Why it happens:** Authors assume "published = done forever." Reality: contact emails become invalid, URLs break, issues reported with no response.

**Prevention:**
1. **Use institutional contact, not personal:**
   - Contact email: lab group email or institutional address (not gmail)
   - Or: Create project email that forwards to multiple people
2. **Plan for handoff:**
   - Document who will maintain dataset if PI leaves
   - Train junior lab member on dataset structure
3. **Set annual review:**
   - Calendar reminder: Review dataset page annually
   - Check: DOI still resolves? Contact email works? Issues reported?
4. **Use stable repository:**
   - Radar4Chem has institutional backing (less likely to disappear than personal website)
   - Avoid personal Dropbox/Google Drive

**Detection warning signs:**
- Contact email is personal gmail
- Only one person knows dataset structure
- No documentation for future maintainers

**Phase-specific guidance:**
- **Phase 1 (Dataset deposit):** Use institutional contact
- **Phase 2 (Post-publication):** Set annual review reminder
- **Phase 3 (PI transition):** Handoff dataset maintenance responsibility

---

## Summary Table: Pitfalls by Phase

| Phase | Critical Pitfalls | Moderate Pitfalls | Minor Pitfalls |
|-------|-------------------|-------------------|----------------|
| **Dataset Preparation** | 1.1 (Incomplete provenance), 1.2 (FAIR F1), 1.3 (FAIR F4) | 1.4 (No versioning plan), 1.5 (Environment docs) | 1.6 (File naming) |
| **Data Descriptor Writing** | 2.1 (Insufficient validation), 2.2 (Methods detail) | 2.3 (DAS format), 2.4 (Code availability) | 2.5 (Figure optimization) |
| **Application Note Writing** | 3.1 (Inadequate benchmarking), 3.2 (Underselling impact) | 3.3 (Page limit violations) | 3.4 (Poor naming) |
| **Cross-Publication Management** | - | 4.1 (Metadata update timing), 4.2 (Citation inconsistency) | 4.3 (No maintenance plan) |

---

## Roadmap Implications

### Phase 1: Dataset Preparation and Repository Deposit
**Primary risk areas:**
- Computational provenance metadata (Pitfall 1.1)
- FAIR compliance (Pitfalls 1.2, 1.3)
- File organization (Pitfall 1.6)

**Recommended approach:**
1. Create input file templates with explicit parameters (addresses 1.1)
2. Generate SMILES/InChI for all molecules early (addresses 1.3)
3. Establish systematic file naming convention (addresses 1.6)
4. Reserve DOI early (addresses 1.2)

**Deep research flags:** None. Standard practices well-documented.

---

### Phase 2: Data Descriptor Writing (Nature Scientific Data)
**Primary risk areas:**
- Technical Validation section (Pitfall 2.1)
- Methods reproducibility detail (Pitfall 2.2)
- Data/code availability statements (Pitfalls 2.3, 2.4)

**Recommended approach:**
1. Perform quantitative validation analysis during data generation (not after)
2. Document workflow as step-by-step protocol with flowchart
3. Deposit code in Zenodo before manuscript submission
4. Draft DAS/CAS with all required elements

**Deep research flags:**
- Phase-specific research on "what validation analyses are standard for NMR datasets?"
- Review 5-10 recent Scientific Data descriptors for computational chemistry datasets

---

### Phase 3: Application Note Writing (J. Cheminformatics)
**Primary risk areas:**
- Benchmarking comprehensiveness (Pitfall 3.1)
- Impact articulation (Pitfall 3.2)
- Page limits (Pitfall 3.3)

**Recommended approach:**
1. Design benchmark matrix early (methods × datasets × metrics)
2. Write to page limit from start (outline with word budgets)
3. Lead with scientific value proposition (not technical features)
4. Include concrete use case scenarios

**Deep research flags:**
- Phase-specific research on "current state-of-the-art in NMR prediction" (for benchmarking)
- Survey recent J. Cheminformatics application notes for structure/style

---

### Phase 4: Cross-Publication Management
**Primary risk areas:**
- Metadata update timing (Pitfall 4.1)
- Citation consistency (Pitfall 4.2)

**Recommended approach:**
1. Create timeline with metadata update milestones
2. Set calendar reminders for metadata updates when papers publish
3. Establish canonical dataset citation format used in both papers
4. Verify bidirectional linking 1 month after publication

**Deep research flags:** None. Procedural/administrative task.

---

## Confidence Assessment

| Pitfall Category | Confidence | Evidence Base |
|------------------|------------|---------------|
| FAIR compliance issues | HIGH | Official NFDI4Chem, DataCite documentation + published literature |
| Computational metadata requirements | HIGH | Published reproducibility studies, NWChem documentation |
| Nature Scientific Data formatting | MEDIUM | Official submission guidelines (accessible) but specific rejection examples limited (not publicly shared) |
| J. Cheminformatics requirements | MEDIUM | Official guidelines + published articles, but limited information on common rejection reasons |
| Cross-referencing mechanics | MEDIUM-HIGH | Crossref/DataCite documentation clear, but implementation examples limited |

**Sources relied upon:**
- NFDI4Chem Knowledge Base (FAIR, Radar4Chem requirements): HIGH confidence
- Published reproducibility studies (Chemistry of Materials, J. Chem. Info. Model.): HIGH confidence
- Journal submission guidelines (Nature Scientific Data, J. Cheminformatics): MEDIUM-HIGH confidence (official but sometimes vague on rejection reasons)
- WebSearch results on benchmarking practices: MEDIUM confidence (good overview but would benefit from official benchmarking framework documentation)

**Low confidence areas requiring validation:**
1. Specific radar4chem metadata schema v9.2 field requirements (found general description but not complete field list)
2. Nature Scientific Data exact Technical Validation expectations (guidelines describe what to include but not quantitative expectations)
3. J. Cheminformatics benchmark sufficiency thresholds (how many comparisons is "enough"?)

**Recommendations:**
- Consult NFDI4Chem directly for Radar4Chem submission checklist (user is NFDI4Chem infrastructure lead, has direct access)
- Review 10-15 recent Scientific Data descriptors for computational chemistry datasets (see what passed peer review)
- Contact J. Cheminformatics editors for application note expectations if unclear

---

## Sources

### FAIR Principles and NFDI4Chem
- [NFDI4Chem FAIR Data Principles](https://knowledgebase.nfdi4chem.de/knowledge_base/docs/fair/)
- [RADAR4Chem Metadata Requirements](https://knowledgebase.nfdi4chem.de/knowledge_base/docs/radar4chem/)
- [NFDI4Chem Infrastructure Providers Guide](https://knowledgebase.nfdi4chem.de/knowledge_base/docs/publishing_standards_infrastructure/)
- [GO FAIR Principles](https://www.go-fair.org/fair-principles/)
- [FAIR Guiding Principles - Nature Scientific Data](https://www.nature.com/articles/sdata201618)

### Computational Chemistry Reproducibility
- [Reproducible Research in Computational Chemistry of Materials - Chemistry of Materials](https://pubs.acs.org/doi/10.1021/acs.chemmater.7b00799)
- [Method and Data Sharing and Reproducibility - J. Chem. Info. Model.](https://pubs.acs.org/doi/10.1021/acs.jcim.0c01389)
- [The Role of Metadata in Reproducible Computational Research - arXiv](https://arxiv.org/pdf/2006.08589)
- [Best Practices in Machine Learning for Chemistry - Nature Chemistry](https://www.nature.com/articles/s41557-021-00716-z)

### NMR Prediction Benchmarking
- [DELTA50 Dataset Paper - Molecules](https://www.mdpi.com/1420-3049/28/6/2449)
- [MoleculeNet Benchmark Framework - Chemical Science](https://pubs.rsc.org/en/content/articlehtml/2018/sc/c7sc02664a)
- [Toward Unified NMR Prediction Benchmark - Nature Computational Science](https://www.nature.com/articles/s43588-025-00783-z)
- [Real-time NMR Prediction with 3D GNN - Chemical Science](https://pubs.rsc.org/en/content/articlehtml/2021/sc/d1sc03343c)

### Journal Guidelines
- [Nature Scientific Data Submission Guidelines](https://www.nature.com/sdata/submission-guidelines)
- [Nature Scientific Data - For Referees](https://www.nature.com/sdata/policies/for-referees)
- [Nature Scientific Data Data Policies](https://www.nature.com/sdata/policies/data-policies)
- [Journal of Cheminformatics Software Article Guidelines](https://jcheminf.biomedcentral.com/submission-guidelines/preparing-your-manuscript/software)
- [BMC Bioinformatics Software Article Guidelines](https://bmcbioinformatics.biomedcentral.com/submission-guidelines/preparing-your-manuscript/software-article)

### Data Citation and Versioning
- [Crossref Relationships Documentation](https://www.crossref.org/documentation/schema-library/markup-guide-metadata-segments/relationships/)
- [Crossref Data Citation Guide](https://www.crossref.org/documentation/reference-linking/data-and-software-citation-deposit-guide/)
- [DataCite Metadata Schema](https://support.datacite.org/docs/datacite-metadata-schema)
- [DataCite Versioning Guidelines](https://support.datacite.org/docs/versioning)
- [Versioning Data Is About More than Revisions - Data Science Journal](https://datascience.codata.org/articles/10.5334/dsj-2021-012/)
- [Benchmark Data Repositories Best Practices - arXiv](https://arxiv.org/html/2410.24100v1)

### Data Availability and Reporting Standards
- [Springer Nature Data Availability Statements](https://www.springernature.com/gp/authors/research-data-policy/data-availability-statements)
- [Nature Portfolio Reporting Standards](https://www.nature.com/nature-portfolio/editorial-policies/reporting-standards)
- [PLOS Data Availability](https://journals.plos.org/plosone/s/data-availability)

### NWChem and COSMO Solvation
- [NWChem COSMO Solvation Model Documentation](https://nwchemgit.github.io/COSMO-Solvation-Model.html)
- [COSMO Solvation Model - Wikipedia](https://en.wikipedia.org/wiki/COSMO_solvation_model)

### Chemistry Reproducibility Crisis
- [Chemistry's Reproducibility Crisis - Chemistry World](https://www.chemistryworld.com/news/chemistrys-reproducibility-crisis-that-youve-probably-never-heard-of/4011693.article)

---

**End of PITFALLS_dataset_publication.md**
