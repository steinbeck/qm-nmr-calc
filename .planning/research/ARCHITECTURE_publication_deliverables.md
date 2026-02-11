# Architecture Patterns: DELTA50 Publication & Paper Deliverables

**Domain:** Computational chemistry dataset publication with companion papers
**Researched:** 2026-02-11
**Confidence:** MEDIUM (verified with official sources for publication workflows, LOW for Radar4Chem specifics)

## Executive Summary

Publishing the DELTA50 benchmark requires coordinating three deliverables with interdependencies: (1) a dataset archive on Radar4Chem, (2) a data descriptor in Scientific Data (Nature), and (3) an application note in Journal of Cheminformatics. The critical architectural decision is the build order: **dataset DOI must be reserved first**, enabling both papers to cite it during peer review, but the dataset can remain embargoed until papers are accepted.

The recommended architecture separates concerns into three parallel workstreams that converge at specific synchronization points, minimizing blocking dependencies while maintaining cross-reference integrity.

## Recommended Architecture

### Three-Deliverable Structure

```
qm-nmr-calc/
├── data/benchmark/delta50/          # Source data (existing)
│   ├── molecules/                   # 50 XYZ structures
│   ├── experimental_shifts.json     # Reference data
│   ├── scaling_factors.json         # Computed results
│   └── SCALING-FACTORS.md          # Analysis documentation
│
├── publications/                    # NEW: All publication deliverables
│   ├── dataset/                     # Deliverable 1: Radar4Chem archive
│   │   ├── README.md               # Dataset documentation
│   │   ├── METHODOLOGY.md          # Computational methods
│   │   ├── molecules/              # 50 structures (copied)
│   │   ├── calculations/           # NEW: NWChem I/O organized by solvent
│   │   ├── results/                # Scaling factors, plots, analysis
│   │   └── metadata.json           # Radar4Chem required metadata
│   │
│   ├── data-descriptor/            # Deliverable 2: Scientific Data paper
│   │   ├── manuscript.tex          # Main LaTeX source
│   │   ├── manuscript.bib          # Bibliography (shared references)
│   │   ├── figures/                # Plots and diagrams
│   │   ├── tables/                 # Data tables
│   │   └── submission/             # Generated submission package
│   │
│   └── application-note/           # Deliverable 3: J Cheminf paper
│       ├── manuscript.tex          # Main LaTeX source
│       ├── manuscript.bib          # Bibliography (shared references)
│       ├── figures/                # Tool screenshots, architecture
│       ├── tables/                 # Performance metrics
│       └── submission/             # Generated submission package
│
└── .planning/                       # Project management (existing)
```

### Component Boundaries

| Component | Responsibility | Inputs | Outputs | Consumers |
|-----------|---------------|--------|---------|-----------|
| **Dataset Package** | Organize 650 NWChem calculations + metadata for archival | data/benchmark/delta50/ | Radar4Chem archive (tar.gz) + reserved DOI | Both papers |
| **Data Descriptor** | Describe dataset creation, validation, reuse value | Dataset DOI, methodology | Scientific Data manuscript | Research community |
| **Application Note** | Describe qm-nmr-calc tool, demonstrate with dataset | Dataset DOI, data descriptor citation, tool architecture | J Cheminf manuscript | Tool users |
| **Build System** | Generate submission packages, validate cross-refs | All manuscripts + dataset | Submission-ready archives | Journal portals |

### Data Flow

```
Phase 1: Dataset Preparation
data/benchmark/delta50/ → publications/dataset/ → Radar4Chem upload (embargoed) → Reserved DOI

Phase 2: Parallel Manuscript Development
Reserved DOI → data-descriptor/manuscript.tex (cite dataset)
Reserved DOI → application-note/manuscript.tex (cite dataset)

Phase 3: Submission Coordination
data-descriptor/ → Scientific Data submission
application-note/ → J Cheminf submission
(Dataset remains embargoed during peer review)

Phase 4: Publication Coordination
Both papers accepted → Publish dataset (DOI goes live)
Data descriptor published → Application note cites descriptor + dataset
```

## Critical Architectural Decisions

### Decision 1: Dataset-First Build Order

**Rationale:** Both papers must cite the dataset DOI. Zenodo and Figshare allow reserving DOIs before publication, but Radar4Chem capabilities are unclear (LOW confidence).

**Options:**
1. **Reserve DOI first** (recommended if Radar4Chem supports this)
   - Upload dataset to Radar4Chem, mark as embargoed
   - Reserve DOI immediately
   - Both papers cite reserved DOI during writing
   - Publish dataset when both papers are accepted

2. **Use Zenodo/Figshare instead** (fallback if Radar4Chem doesn't support embargo)
   - These platforms explicitly support DOI reservation
   - Upload embargoed dataset, get reserved DOI
   - Cite in papers during peer review
   - Publish when ready

3. **Dataset last** (NOT recommended)
   - Write papers without dataset DOI
   - Add DOI during final revisions
   - Creates timing risk if papers accept at different times

**Chosen:** Option 1 with Option 2 as fallback. Verify Radar4Chem embargo capabilities early in Phase 1.

### Decision 2: Paper Submission Order

**Rationale:** The application note references the data descriptor, creating a dependency. However, both can be submitted in parallel if the data descriptor citation is included as "submitted" or "in review."

**Options:**
1. **Sequential:** Data descriptor first, wait for acceptance, then submit application note
   - Pro: Application note can cite published descriptor
   - Con: Delays application note by 3-6 months (peer review timeline)

2. **Parallel with provisional citation:** Submit both simultaneously
   - Pro: Faster overall timeline
   - Con: Application note must update descriptor citation if status changes
   - Mitigation: Use "submitted to Scientific Data" citation initially

3. **Application note first:** Submit tool paper before dataset paper
   - Con: Backwards dependency (tool paper should reference data descriptor)
   - NOT recommended

**Chosen:** Option 2 (parallel submission). The application note manuscript.tex includes a citation like:
```latex
\cite{AuthorYear-submitted}  % Dataset descriptor (submitted to Scientific Data)
```
Updated to journal reference once descriptor is accepted.

### Decision 3: Directory Structure for Manuscripts

**Rationale:** LaTeX submission systems (particularly Elsevier's Editorial Manager) cannot process folders with directory structures. All files must be in a flat directory for submission.

**Architecture:**
- **Source structure:** Use subdirectories for development (figures/, tables/)
- **Submission structure:** Flatten before submission (all files in submission/)
- **Build system:** Script to copy and flatten files for journal submission

**Example:**
```bash
# Development structure (human-friendly)
data-descriptor/
├── manuscript.tex
├── manuscript.bib
├── figures/
│   ├── fig1-workflow.pdf
│   └── fig2-validation.pdf
└── tables/
    └── table1-statistics.csv

# Submission structure (journal-friendly)
data-descriptor/submission/
├── manuscript.tex
├── manuscript.bib
├── fig1-workflow.pdf
├── fig2-validation.pdf
└── table1-statistics.csv  # Converted to LaTeX table
```

### Decision 4: Dataset Archive Organization

**Rationale:** 650 NWChem calculations (13 solvents × 50 molecules) with .nw inputs and .out outputs are large but compressible. Zenodo's 50GB limit (100GB exceptional) should be sufficient. Radar4Chem accepts all file types but structure is unclear.

**Recommended organization:**
```
delta50-dataset/
├── README.md                    # Dataset overview, citation info
├── METHODOLOGY.md               # Computational methods in detail
├── LICENSE.txt                  # CC-BY-4.0 or CC0
├── metadata.json                # Machine-readable metadata
├── molecules/
│   ├── compound_01.xyz          # 50 structures
│   └── ...
├── calculations/
│   ├── vacuum/
│   │   ├── compound_01/
│   │   │   ├── optimize.nw      # Geometry optimization input
│   │   │   ├── optimize.out     # Output
│   │   │   ├── shielding.nw     # NMR calculation input
│   │   │   └── shielding.out    # Output with shielding tensors
│   │   └── compound_02/
│   │       └── ...
│   ├── CHCl3/
│   │   └── ...                  # Same structure per solvent
│   ├── DMSO/
│   ├── Methanol/
│   └── ...                      # 13 solvents total
└── results/
    ├── scaling_factors.json     # Per-solvent scaling parameters
    ├── experimental_shifts.json # Reference data
    ├── statistics.csv           # MAE, RMSD, R² per solvent
    └── plots/
        ├── correlation-CHCl3-1H.pdf
        └── ...
```

**Size estimation:**
- 650 calculations × ~2 files (avg .nw + .out) × ~50KB avg = ~65MB compressed
- Figures/plots: ~5MB
- Total: ~70MB (well under limits)

**Archive format:** tar.gz for repository upload, organized as above.

## Integration Points

### Integration 1: Dataset DOI → Both Papers

**What:** Both manuscripts cite the dataset DOI in Data Availability sections.

**When:** After dataset DOI is reserved (Phase 1 complete).

**How:**
- Data descriptor: "The DELTA50 dataset is available at https://doi.org/10.XXXX/YYYY [DOI reserved]"
- Application note: "Benchmark calculations use the DELTA50 dataset (https://doi.org/10.XXXX/YYYY)"

**Validation:** Check that both .bib files include identical dataset reference.

### Integration 2: Data Descriptor Citation → Application Note

**What:** Application note cites the data descriptor to justify dataset quality.

**When:** After data descriptor is submitted (Phase 2).

**How:**
```bibtex
@article{AuthorYear-descriptor,
  title={The DELTA50 Dataset: 650 NWChem Calculations for NMR Chemical Shift Prediction},
  author={Author Names},
  journal={Scientific Data},
  note={Submitted},  % Update to volume/pages when published
  year={2026}
}
```

**Validation:** Update citation status as descriptor progresses through peer review.

### Integration 3: Application Note → Tool Repository

**What:** Application note describes qm-nmr-calc architecture from docs/architecture.md.

**When:** During manuscript writing (Phase 2).

**How:** Reference existing documentation, add citation to code repository.

**Validation:** Ensure architecture described in paper matches current tool implementation.

### Integration 4: Shared Bibliography

**What:** Both papers cite overlapping references (NWChem, DFT methods, DELTA50 original paper).

**When:** Throughout manuscript writing.

**How:** Maintain a shared `publications/shared-references.bib` file, symlinked or copied into each manuscript directory.

**Validation:** Use biblatex or bibtex to ensure consistent citation styles across papers.

## Patterns to Follow

### Pattern 1: DOI Reservation Before Writing

**What:** Reserve dataset DOI before writing manuscripts that cite it.

**When:** Immediately after dataset archive is prepared.

**Why:** Enables consistent citation during peer review, avoids last-minute updates.

**Example:**
```bash
# Phase 1: Prepare dataset
cd publications/dataset
tar -czf delta50-dataset.tar.gz delta50-dataset/

# Upload to Radar4Chem (or Zenodo), mark as embargoed
# Receive reserved DOI: 10.XXXX/delta50-nmr-benchmark

# Phase 2: Write papers with DOI
echo "Dataset DOI: 10.XXXX/delta50-nmr-benchmark" >> ../RESERVED-DOI.txt
# Both manuscripts cite this DOI from the start
```

### Pattern 2: Submission Package Generation

**What:** Automate flattening directory structure for journal submission.

**When:** Before submitting to journals.

**Why:** Ensures compliance with LaTeX submission system requirements.

**Example script:**
```bash
#!/bin/bash
# publications/data-descriptor/build-submission.sh

OUTDIR="submission"
rm -rf "$OUTDIR"
mkdir -p "$OUTDIR"

# Copy main files
cp manuscript.tex "$OUTDIR/"
cp manuscript.bib "$OUTDIR/"

# Flatten figures
cp figures/*.pdf "$OUTDIR/"

# Flatten tables (convert CSV to LaTeX if needed)
cp tables/*.tex "$OUTDIR/"

# Create submission archive
cd "$OUTDIR"
zip ../data-descriptor-submission.zip *
cd ..
echo "Submission package: data-descriptor-submission.zip"
```

### Pattern 3: Cross-Reference Validation

**What:** Automated checks that all cross-references between deliverables are valid.

**When:** Before each submission.

**Why:** Prevents inconsistent citations across papers.

**Example validation:**
```python
# publications/validate-crossrefs.py
import re

def extract_doi(file_path):
    """Extract DOI citations from LaTeX file."""
    with open(file_path) as f:
        content = f.read()
    return re.findall(r'doi\.org/(10\.\S+)', content)

# Check both papers cite the same dataset DOI
descriptor_dois = extract_doi('data-descriptor/manuscript.tex')
appnote_dois = extract_doi('application-note/manuscript.tex')
dataset_doi = open('RESERVED-DOI.txt').read().strip().split(': ')[1]

assert dataset_doi in descriptor_dois, "Data descriptor missing dataset DOI"
assert dataset_doi in appnote_dois, "Application note missing dataset DOI"
print("✓ All cross-references valid")
```

### Pattern 4: Embargo Coordination

**What:** Keep dataset embargoed until both papers are accepted.

**When:** Throughout peer review process.

**Why:** Maintains novelty for papers, allows reviewers to access data via private links.

**Timeline:**
```
Week 0:  Upload embargoed dataset → Reserve DOI
Week 1:  Submit both manuscripts (cite reserved DOI)
Week 12: Data descriptor accepted
Week 18: Application note accepted
Week 19: Publish dataset (DOI goes live)
Week 20: Both papers published with live dataset links
```

**Coordination:** Track paper status in shared spreadsheet:
```
Deliverable          | Status      | DOI/Citation
---------------------|-------------|------------------
Dataset (Radar4Chem) | Embargoed   | 10.XXXX/delta50 (reserved)
Data Descriptor      | In Review   | Submitted to Sci Data
Application Note     | In Review   | Submitted to J Cheminf
```

## Anti-Patterns to Avoid

### Anti-Pattern 1: Writing Papers Before Dataset DOI Exists

**What:** Write manuscripts with placeholder "DOI: TBD" that must be updated later.

**Why bad:**
- Requires updating both manuscripts if DOI comes late
- May require re-submission if DOI changes after initial review
- Creates coordination debt

**Instead:** Reserve DOI first (Pattern 1), write papers with final DOI from the start.

### Anti-Pattern 2: Nested Subdirectories in Submission Package

**What:** Submit LaTeX manuscripts with figures in `figures/` and tables in `tables/` subdirectories.

**Why bad:**
- Elsevier Editorial Manager (and similar systems) cannot process nested directories
- Causes submission failures
- Requires manual flattening at last minute

**Instead:** Use Pattern 2 (build script) to flatten automatically before submission.

### Anti-Pattern 3: Publishing Dataset Before Papers Are Ready

**What:** Publish dataset immediately after upload to get DOI, before papers are written.

**Why bad:**
- Dataset becomes public before papers establish its scientific value
- Reviewers may request dataset changes, but it's already published
- Loses coordination with paper publication timing

**Instead:** Use embargo feature (Pattern 4) to reserve DOI while keeping data private during peer review.

### Anti-Pattern 4: Independent Bibliography Files

**What:** Maintain separate .bib files for each paper with duplicate entries.

**Why bad:**
- Inconsistent citation formatting across papers
- Duplicate maintenance effort when references change
- Risk of citing different versions of same paper

**Instead:** Use shared bibliography (Integration 4) with consistent entries.

### Anti-Pattern 5: Dataset Archive Without README

**What:** Upload tar.gz of raw NWChem files without documentation explaining structure.

**Why bad:**
- Dataset is not reusable without understanding file organization
- Violates FAIR principles (Findable, Accessible, Interoperable, Reusable)
- Data descriptor paper must compensate for missing dataset documentation

**Instead:** Include comprehensive README.md in dataset archive (Decision 4 structure).

## Scalability Considerations

| Concern | At Current Scale (650 calcs) | At 10× Scale (6500 calcs) | At 100× Scale (65K calcs) |
|---------|------------------------------|---------------------------|---------------------------|
| **Dataset Size** | ~70MB compressed | ~700MB (under Zenodo 50GB limit) | ~7GB (may need institutional repo) |
| **File Organization** | Solvent/molecule hierarchy works | Same structure, more solvents/molecules | May need sharding (e.g., batches of 1000) |
| **Metadata** | Single metadata.json | Single JSON still manageable | May need database + JSON export |
| **Paper Length** | Fits Scientific Data format | May need supplementary methods | Definitely needs supplementary info |
| **Build Time** | Seconds to flatten submission | ~1 minute to tar.gz dataset | May need parallel compression |

**Recommendation:** Current architecture scales to 10× without changes. At 100× scale, consider:
- Institutional repository instead of Zenodo/Radar4Chem
- Database-backed metadata instead of flat JSON
- Supplementary information for data descriptor instead of main text

## Open Questions and Research Gaps

### HIGH Priority (blocking Phase 1)

1. **Radar4Chem embargo capability:** Does Radar4Chem support reserving DOIs for embargoed datasets?
   - **Confidence:** LOW (not found in documentation)
   - **Fallback:** Use Zenodo or Figshare (confirmed to support embargo + DOI reservation)
   - **Resolution:** Contact Radar4Chem support or test with dummy upload in Phase 1

2. **Radar4Chem metadata requirements:** What fields are required in metadata.json?
   - **Confidence:** LOW (found general guidance but not schema)
   - **Fallback:** Use Dublin Core metadata standard (widely accepted)
   - **Resolution:** Consult Radar4Chem Quickstart Guide or support

### MEDIUM Priority (informative for Phase 2)

3. **Scientific Data figure limit:** Is 8 figures a hard limit or guideline?
   - **Confidence:** MEDIUM (documentation says "recommended, not mandate")
   - **Impact:** May need to consolidate plots into multi-panel figures
   - **Resolution:** Review recent data descriptors for typical figure count

4. **J Cheminf article type:** Does "Software" article type fit better than "Research"?
   - **Confidence:** MEDIUM (found Software type in guidelines)
   - **Impact:** Software articles may have different length/scope requirements
   - **Resolution:** Review J Cheminf submission guidelines in detail during Phase 2

### LOW Priority (optimization)

5. **Shared reference management:** Should publications/ use a tool like Zotero or just shared .bib?
   - **Confidence:** HIGH (either works, tooling is preference)
   - **Impact:** Minimal, affects workflow not architecture
   - **Resolution:** Decide based on team preference in Phase 2

## Build Order Summary

**Critical Path:**
```
1. Prepare dataset archive (publications/dataset/)
2. Upload to Radar4Chem (or fallback to Zenodo)
3. Reserve DOI (embargoed dataset)
4. Write data descriptor manuscript (cite reserved DOI)
5. Write application note manuscript (cite reserved DOI + submitted descriptor)
6. Submit both papers in parallel
7. Iterate on peer review (dataset remains embargoed)
8. Both papers accepted → Publish dataset (DOI goes live)
9. Data descriptor published → Update application note citation
10. Application note published → Complete
```

**Parallelization opportunities:**
- Steps 4-5: Write both manuscripts simultaneously (different authors or time-boxed)
- Steps 6-7: Both papers in peer review concurrently
- Build scripts (Pattern 2) and validation (Pattern 3) can be developed anytime

**Blocking dependencies:**
- Step 3 blocks Step 4 and Step 5 (both papers need DOI)
- Step 4 (submission) informs Step 5 (application note can cite "submitted" descriptor)
- Step 8 blocks final publication (dataset DOI must go live before papers cite it as published)

## Recommended Phase Structure for Roadmap

Based on this architecture, suggest the following phase breakdown:

1. **Phase: Dataset Archive Preparation**
   - Create publications/dataset/ structure
   - Collect 650 NWChem calculations from data/jobs/ or rerun if needed
   - Organize into solvent/molecule hierarchy
   - Write README.md and METHODOLOGY.md
   - Generate metadata.json

2. **Phase: DOI Reservation**
   - Research Radar4Chem embargo (resolve Open Question 1)
   - Upload embargoed dataset
   - Reserve DOI
   - Document DOI in publications/RESERVED-DOI.txt

3. **Phase: Data Descriptor Manuscript**
   - Set up publications/data-descriptor/ LaTeX project
   - Write Methods, Data Records, Technical Validation sections
   - Create figures (correlation plots, validation statistics)
   - Create submission package (Pattern 2)
   - Submit to Scientific Data

4. **Phase: Application Note Manuscript**
   - Set up publications/application-note/ LaTeX project
   - Describe qm-nmr-calc architecture and features
   - Demonstrate with DELTA50 dataset
   - Create figures (architecture diagram, performance metrics)
   - Create submission package (Pattern 2)
   - Submit to Journal of Cheminformatics

5. **Phase: Peer Review Coordination**
   - Track status of both papers
   - Respond to reviewer comments
   - Update cross-references if needed (e.g., descriptor gets DOI)
   - Coordinate acceptance timing

6. **Phase: Publication**
   - Publish dataset when both papers are accepted
   - Update application note citation when descriptor is published
   - Finalize both papers
   - Archive preprints (arXiv or institutional repo)

**Phase dependencies:**
- Phases 3-4 can run in parallel after Phase 2 completes
- Phase 5 blocks Phase 6
- Phase 6 is final synchronization point

## Sources

**Dataset Publication:**
- [RADAR4Chem overview](https://radar.products.fiz-karlsruhe.de/en/radarabout/radar4chem)
- [RADAR4Chem documentation](https://chemotion.net/docs/eln/interfaces/radar)
- [Zenodo DOI reservation](https://help.zenodo.org/docs/deposit/describe-records/reserve-doi/)
- [Figshare DOI reservation](https://info.figshare.com/user-guide/how-to-reserve-a-doi/)
- [Figshare publishing workflow](https://help.figshare.com/article/publishing-a-dataset-at-the-same-time-as-the-associated-paper)

**Data Descriptor Guidelines:**
- [Scientific Data submission guidelines](https://www.nature.com/sdata/submission-guidelines)
- [Scientific Data data policies](https://www.nature.com/sdata/policies/data-policies)
- [Scientific Data for referees](https://www.nature.com/sdata/policies/for-referees)
- [Data citation roadmap](https://www.nature.com/articles/sdata2018259)

**Journal of Cheminformatics:**
- [J Cheminf submission guidelines](https://jcheminf.biomedcentral.com/submission-guidelines)
- [J Cheminf software article type](https://jcheminf.biomedcentral.com/submission-guidelines/preparing-your-manuscript/software)

**Computational Chemistry Datasets:**
- [Managing computational chemistry big data](https://pubs.acs.org/doi/10.1021/ci500593j)
- [Digital data repositories in chemistry](https://pubs.acs.org/doi/10.1021/ci500302p)
- [NMRexp 3.3M spectra dataset](https://www.nature.com/articles/s41597-025-06245-5)
- [DELTA50 dataset paper](https://www.mdpi.com/1420-3049/28/6/2449)

**LaTeX Manuscript Structure:**
- [Elsevier LaTeX instructions](https://www.elsevier.com/researcher/author/policies-and-guidelines/latex-instructions)
- [Springer Nature LaTeX support](https://www.springernature.com/gp/authors/campaigns/latex-author-support)

**Confidence Assessment:**
- DOI reservation workflow: HIGH (verified with Zenodo/Figshare official docs)
- Scientific Data format: HIGH (official submission guidelines)
- Radar4Chem specifics: LOW (general info only, lacking technical details)
- Build order strategy: MEDIUM (inferred from multiple sources, not explicitly documented)
- Dataset organization: MEDIUM (based on similar datasets, not domain standard)
