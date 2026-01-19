# Features Research

**Domain:** Async NMR QM Calculation Web Service
**Researched:** 2026-01-19
**Confidence:** MEDIUM (verified against multiple authoritative sources)

## Similar Tools

The NMR prediction landscape splits into two distinct categories:

### Database/ML-Based Prediction Services (Fast, Approximate)

| Tool | Type | Strengths | Limitations |
|------|------|-----------|-------------|
| [nmrdb.org](https://www.nmrdb.org/) | Web, Free | 1H/13C/2D prediction, neural network based, instant results | No QM accuracy, limited to trained chemical space |
| [nmrshiftdb2](https://nmrshiftdb.nmr.uni-koeln.de/) | Web, Free | Open database, 1H/13C/other nuclei, searchable | HOSE-code based, database-limited accuracy |
| [NMRium](https://www.nmrium.com/) | Web, Free | Processing + prediction, excellent UI, drag-drop, publication export | Prediction is ML-based, not QM |
| [PROSPRE](https://prospre.ca/) | Web, Free | ML-based 1H, multiple solvents (D2O, CDCl3, DMSO, CD3OD) | 1H only, no 13C |
| [ACD/Labs NMR Predictors](https://www.acdlabs.com/products/spectrus-platform/nmr-predictors/) | Commercial | HOSE + neural network, extensive database, many nuclei | Expensive, desktop software |
| [Mnova NMRPredict](https://mestrelab.com/main-product/nmr-predict) | Commercial | Ensemble prediction, excellent accuracy, 1H/13C/15N/19F/31P | Commercial license required |

### QM-Based Calculation Tools (Slow, Accurate)

| Tool | Type | Strengths | Limitations |
|------|------|-----------|-------------|
| [Gaussian](https://gaussian.com/nmr/) | Desktop, Commercial | Industry standard, GIAO method, extensive documentation | Expensive, requires local compute |
| [NWChem](https://nwchemgit.github.io/) | Desktop, Open Source | Free, DFT NMR shieldings, ISiCLE integration | Requires expertise, no web UI |
| [ORCA](https://www.faccts.de/orca/) | Desktop, Free (academic) | Good 19F NMR, scalable | Command-line only |
| [WebMO](https://www.webmo.net/) | Web Interface | Web UI for multiple QM packages including NWChem, job management | Not NMR-specialized, general QC interface |

**Key Insight:** There is no publicly available web service that provides QM-level NMR predictions with a clean REST API. This is the gap qm-nmr-calc fills.

### ISiCLE (Our Calculation Engine)

[ISiCLE](https://github.com/pnnl/isicle) from PNNL provides:
- Automated DFT NMR chemical shift calculations via NWChem
- Conformer generation and Boltzmann weighting
- Solvent effects via COSMO model
- Support for multiple DFT methods (B3LYP, BLYP, etc.)
- Validated accuracy: MAE 0.33 ppm (1H), 3.93 ppm (13C) with B3LYP/cc-pVDZ

Sources: [ISiCLE GitHub](https://github.com/pnnl/isicle), [ISiCLE NMR Paper](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-018-0305-8)

---

## Table Stakes

Features users **expect** from any NMR calculation service. Missing these makes the product feel incomplete.

### Molecular Input

| Feature | Why Expected | Complexity | Notes |
|---------|--------------|------------|-------|
| **SMILES string input** | Universal, copy-paste from anywhere | Low | Canonical SMILES preferred |
| **Structure file upload** (SDF/MOL) | Standard exchange format, preserves 3D if available | Low | RDKit handles parsing |
| **Structure drawing widget** | Users without SMILES need visual input | Medium | JSME or Ketcher embed, not custom |
| **Input validation** | Chemically impossible structures waste compute | Low | Valence checks, aromaticity handling |
| **Example molecules** | Helps users test without own data | Low | Include 2-3 well-known structures |

**Recommendation:** Start with SMILES + file upload. Add drawing widget in v1 if time permits, otherwise v1.1.

### Calculation Options

| Feature | Why Expected | Complexity | Notes |
|---------|--------------|------------|-------|
| **Preset levels** (Fast/Standard/Accurate) | Users don't know DFT parameters | Low | Map to basis set/functional combos |
| **Nuclei selection** (1H, 13C) | Users may only need one | Low | Default to both |
| **Solvent selection** | Shifts depend on solvent | Low | Chloroform, DMSO, water common |
| **Job naming** | Users track multiple calculations | Low | Optional field |

**Recommended Presets:**
- **Draft/Fast:** B3LYP/6-31G(d) - minutes, ~0.5 ppm 1H accuracy
- **Standard:** B3LYP/cc-pVDZ - tens of minutes, ~0.33 ppm 1H accuracy (ISiCLE validated)
- **Publication:** B3LYP/cc-pVTZ or pcS-2 - hours, highest accuracy

Sources: [ISiCLE Paper](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-018-0305-8), [ResearchGate DFT Discussion](https://www.researchgate.net/post/What-is-the-best-DFT-functional-to-perform-NMR-calculations-in-Gaussian)

### Job Management

| Feature | Why Expected | Complexity | Notes |
|---------|--------------|------------|-------|
| **Job ID returned immediately** | Async pattern standard | Low | UUID or similar |
| **Status endpoint** | Poll for completion | Low | queued/running/completed/failed |
| **Job list/history** | Track past calculations | Low | Filesystem-backed for v1 |
| **Estimated completion time** | Set expectations | Medium | Rough estimate from molecule size |

### Results Delivery

| Feature | Why Expected | Complexity | Notes |
|---------|--------------|------------|-------|
| **Chemical shifts as JSON** | Programmatic access, atom-indexed | Low | Core deliverable |
| **Atom assignment** | Know which shift belongs to which atom | Medium | ISiCLE provides this |
| **Optimized geometry** | Users need final structure | Low | XYZ or SDF format |
| **Raw output files** | Experts want NWChem logs | Low | Zip download |
| **Downloadable results** | Users need local copies | Low | Multiple formats |

### Notifications

| Feature | Why Expected | Complexity | Notes |
|---------|--------------|------------|-------|
| **Email on completion** | Don't need to poll constantly | Medium | Optional, requires email service |
| **Email on failure** | Know when to investigate | Medium | Include error summary |

Sources: [AWS Async API Patterns](https://aws.amazon.com/blogs/architecture/managing-asynchronous-workflows-with-a-rest-api/), [Azure Async Request-Reply](https://learn.microsoft.com/en-us/azure/architecture/patterns/async-request-reply)

---

## Differentiators

Features that would **stand out**. Not expected, but create competitive advantage and delight.

### Enhanced Visualization

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| **Simulated spectrum plot** | Visual interpretation, publication-ready | Medium | Lorentzian line shapes, configurable |
| **Annotated structure image** | See shifts on molecule | Medium | RDKit drawing + label overlay |
| **Interactive spectrum viewer** | Zoom, peak selection, hover details | High | Consider NMRium embed or SpeckTackle |
| **Structure-spectrum linking** | Click peak, highlight atom | High | v2+ feature, powerful for assignment |

**Recommendation:** Static spectrum plot + annotated structure for v1. Interactive viewer for v1.1+.

### Advanced Calculation Features

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| **Conformer ensemble averaging** | More accurate for flexible molecules | Medium | ISiCLE supports this natively |
| **Multiple theory level comparison** | See accuracy vs speed tradeoff | Medium | Run same molecule at multiple levels |
| **J-coupling calculation** | Complete NMR picture | High | Computationally expensive |
| **Scaling factors** | Empirical correction for systematic error | Low | Literature values available |
| **Custom DFT parameters** | Expert users need fine control | Low | Expose basis/functional/solvation |

### API Quality

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| **OpenAPI/Swagger docs** | Easy integration | Low | Auto-generate from FastAPI |
| **Webhook callbacks** | No polling needed | Medium | POST to user URL on completion |
| **Batch submission** | Submit multiple molecules | Medium | API-first, UI can follow |
| **Rate limiting with clear errors** | Fair use, predictable behavior | Low | Standard API practice |

### User Experience

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| **Progress indicators** | Know calculation is working | Medium | ISiCLE stages: opt -> NMR -> postprocess |
| **Calculation queue position** | Set timing expectations | Low | Simple counter |
| **Shareable result links** | Collaboration, reproducibility | Low | Unique URLs per job |
| **Result comparison** | Compare different methods/molecules | Medium | Side-by-side view |

Sources: [WebMO](https://www.webmo.net/), [NMRium Features](https://www.nmrium.com/features), [SpeckTackle](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-015-0065-7)

---

## Anti-Features

Features to **deliberately NOT build** for v1. Common in the space but add complexity without proportional value.

### Spectrum Processing

| Anti-Feature | Why Avoid | What to Do Instead |
|--------------|-----------|-------------------|
| **FID processing** | We predict, not process experimental data | Return simulated spectrum only |
| **Phase/baseline correction** | Processing feature, not prediction | Out of scope |
| **Peak picking on uploaded spectra** | Different problem domain | Focus on prediction |

**Rationale:** NMRium and Mnova already do processing well. Our value is QM prediction.

### Calculation Scope Creep

| Anti-Feature | Why Avoid | What to Do Instead |
|--------------|-----------|-------------------|
| **Other nuclei (15N, 31P, 19F)** | Different parameterization, less validated | 1H/13C only for v1, expand later |
| **Protein NMR** | Completely different scale and methods | Small molecule focus |
| **Dynamics/relaxation** | Different calculation type | Static shieldings only |
| **Reaction mechanism NMR** | Much more complex workflow | Single molecule only |

**Rationale:** Scope discipline. ISiCLE is validated for 1H/13C small molecules.

### Infrastructure Complexity

| Anti-Feature | Why Avoid | What to Do Instead |
|--------------|-----------|-------------------|
| **Real-time WebSocket updates** | Overkill for minute+ calculations | Polling + email sufficient |
| **User accounts for v1** | Adds auth complexity | Job-based access, single user |
| **Cloud auto-scaling** | Single VM constraint | Fixed capacity, queue management |
| **Database** | Adds operational burden | Filesystem storage |

**Rationale:** Single VM deployment, keep simple.

### UI Scope

| Anti-Feature | Why Avoid | What to Do Instead |
|--------------|-----------|-------------------|
| **Mobile-optimized UI** | Desktop users, complex visualizations | Responsive but desktop-first |
| **Custom spectrum editor** | Processing, not our domain | Static/read-only visualization |
| **Molecule library/favorites** | Requires user accounts | Export/import jobs |

Sources: [PROJECT.md out of scope](file:///home/christoph_steinbeck_gmail_com/develop/qm-nmr-calc/.planning/PROJECT.md)

---

## User Experience Patterns

How similar tools handle common workflows.

### Job Submission Flow

**Best practices observed:**

1. **Progressive disclosure**
   - Simple by default (SMILES + preset)
   - Advanced options collapsed/hidden
   - NMRium, WebMO both use this pattern

2. **Immediate feedback**
   - Validate input before submission
   - Show molecule preview from SMILES
   - PROSPRE shows structure drawing after SMILES paste

3. **Clear confirmation**
   - Return job ID immediately
   - Show expected wait time
   - Provide status URL

**Recommended flow for qm-nmr-calc:**
```
[Input SMILES/file] -> [Preview molecule] -> [Select preset] -> [Submit]
     |                       |                    |
     v                       v                    v
  Validation            2D structure          Optional: advanced
  feedback              drawing shown         params expandable
                                                    |
                                                    v
                                           [Job ID + Status URL]
```

### Status Checking Patterns

**Polling (most common):**
- Return `202 Accepted` with `Location` header to status endpoint
- Status endpoint returns `{status, progress, result_url}`
- Client polls every N seconds

**Webhook (advanced):**
- Client provides callback URL at submission
- Server POSTs to callback on completion
- Fallback to polling always available

**Email (supplementary):**
- Optional email field at submission
- Send on completion with result link
- Doesn't replace API status

**Recommended:** Polling as primary, email as optional, webhook as v1.1 enhancement.

Sources: [Azure Async Pattern](https://learn.microsoft.com/en-us/azure/architecture/patterns/async-request-reply), [AWS Async Workflows](https://aws.amazon.com/blogs/architecture/managing-asynchronous-workflows-with-a-rest-api/)

### Results Presentation

**Best practices observed:**

1. **Summary first, details on demand**
   - Key results (shifts) visible immediately
   - Raw files downloadable but not prominent
   - NMRium shows spectrum first, data tables secondary

2. **Multiple download formats**
   - JSON for programmatic use
   - CSV/text for spreadsheet users
   - Images (PNG/SVG) for publications
   - Raw files (zip) for experts

3. **Visual + tabular**
   - Spectrum plot for quick interpretation
   - Table of shifts for precise values
   - PROSPRE returns table, nmrdb.org returns spectrum + table

**Recommended results page structure:**
```
+---------------------------+
|  Job Info (name, params)  |
+---------------------------+
|                           |
|   [Spectrum Plot]         |
|                           |
+-------------+-------------+
| [Annotated  | [Shifts     |
|  Structure] |  Table]     |
+-------------+-------------+
| Downloads:                |
| [JSON] [CSV] [PNG] [ZIP]  |
+---------------------------+
```

### Input Methods Comparison

| Tool | SMILES | InChI | File Upload | Drawing | Drag-Drop |
|------|--------|-------|-------------|---------|-----------|
| PROSPRE | Yes | Yes | SDF | MarvinJS | No |
| nmrdb.org | Yes | No | No | JSME | No |
| NMRium | No | No | Yes | No | Yes |
| WebMO | No | No | Yes | Built-in 3D | No |

**Recommendation:** SMILES + SDF file upload for v1. Drawing widget (JSME or Ketcher) for v1.1.

---

## Feature Dependencies

```
Input Validation
     |
     v
Molecule Parsing (RDKit)
     |
     +---> 2D Structure Preview
     |
     v
ISiCLE Calculation
     |
     +---> Geometry Optimization
     |           |
     |           v
     |     NMR Shielding Calculation
     |           |
     +---> Raw Output Files
     |
     v
Post-processing
     |
     +---> JSON Shifts (with atom assignment)
     +---> Spectrum Generation
     +---> Annotated Structure Image
     |
     v
Results Storage
     |
     +---> API Access
     +---> Web UI Display
     +---> Download Endpoints
```

---

## MVP Recommendation

For v1.0, prioritize in this order:

### Must Have (Week 1-2 focus)
1. SMILES + MOL/SDF file input with validation
2. Three calculation presets (Fast/Standard/Publication)
3. Job submission returning job ID
4. Status polling endpoint
5. JSON results with atom-assigned shifts
6. Raw file download (zip)

### Should Have (Week 3-4 focus)
7. Optimized geometry download (XYZ/SDF)
8. Simulated spectrum plot (static PNG)
9. Annotated structure image
10. Email notification (optional field)
11. Job history/list endpoint
12. Clean web UI for all above

### Nice to Have (if time permits)
13. Advanced parameter exposure (basis/functional/solvent)
14. Molecule preview before submission
15. Estimated completion time
16. Queue position indicator

### Defer to v1.1+
- Interactive spectrum viewer
- Structure drawing widget
- Webhook callbacks
- Batch submission
- Result comparison view
- Additional nuclei (15N, 19F, 31P)

---

## Sources

### Primary (HIGH confidence)
- [ISiCLE GitHub](https://github.com/pnnl/isicle) - Calculation engine documentation
- [ISiCLE NMR Paper](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-018-0305-8) - Accuracy benchmarks
- [NWChem Documentation](https://nwchemgit.github.io/Properties.html) - DFT NMR methods

### Secondary (MEDIUM confidence)
- [NMRium Features](https://www.nmrium.com/features) - UI/UX patterns
- [PROSPRE](https://prospre.ca/) - Input/output patterns
- [nmrdb.org](https://www.nmrdb.org/) - Prediction UI patterns
- [WebMO](https://www.webmo.net/) - Job management patterns
- [Azure Async Pattern](https://learn.microsoft.com/en-us/azure/architecture/patterns/async-request-reply) - API design

### Tertiary (LOW confidence - community sources)
- [ResearchGate DFT Discussion](https://www.researchgate.net/post/What-is-the-best-DFT-functional-to-perform-NMR-calculations-in-Gaussian) - Functional/basis set choices
