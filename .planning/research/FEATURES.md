# Feature Research: Conformational Sampling for NMR Predictions

**Domain:** Conformational Sampling for NMR Predictions
**Researched:** 2026-01-26
**Confidence:** HIGH

## Feature Landscape

Research across ISiCLE, CENSO, Spartan, and academic literature reveals clear patterns in how conformational sampling features are implemented for NMR prediction. The landscape divides into table stakes (features users expect from any ensemble NMR tool), differentiators (features that provide competitive advantage), and anti-features (commonly requested but problematic features to avoid).

### Table Stakes (Users Expect These)

| Feature | Why Expected | Complexity | Notes |
|---------|--------------|------------|-------|
| Conformer generation method choice | Different methods trade speed vs accuracy; users need control | Medium | RDKit (fast, drug-like molecules) vs CREST (slow, high accuracy). ISiCLE and CENSO both support multiple backends. |
| Energy window filtering | Exclude high-energy irrelevant conformers; standard practice | Low | Typical: 3-6 kcal/mol. Two-stage filtering (pre-DFT wide, post-DFT tight) is common. |
| Boltzmann-weighted average shifts | Core purpose of ensemble calculations; non-negotiable | Medium | Weight by DFT energies from optimization step. Formula: δ_avg = Σ p_i × δ_i where p_i = exp(-E_i/RT) / Σ exp(-E_j/RT) |
| Temperature parameter | Boltzmann weighting is temperature-dependent; users need control | Low | Default 298.15 K. Essential for experimental condition matching. |
| Single-conformer fallback mode | Rigid molecules don't need ensemble; wastes compute | Low | Users explicitly choose single vs ensemble. Avoids forcing ensemble on benzene. |
| Per-atom shift assignments | Users need atom-by-atom predictions for structure elucidation | Low | Already in v1.x. Ensemble mode returns weighted average per atom. |
| Lowest-energy geometry return | Users expect the "most likely" structure | Low | Standard output. Important for 3D visualization. |

### Differentiators (Competitive Advantage)

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| Automatic CREST detection | Zero-config high-accuracy mode when xtb/crest installed | Low | Most tools require manual backend selection. Auto-detection provides better UX. |
| Weighted-average-only API | Simpler, cleaner API focused on experimental comparison | Low | Spartan/CENSO show per-conformer detail. Weighted average only is sufficient for structure validation. Decision: v2.0 ships this way. |
| Population % metadata | Shows ensemble composition without full conformer dump | Medium | Response includes: conformer count, energy range, top 3 conformer populations. Useful transparency without overwhelming detail. |
| Two-stage energy filtering | Reduces DFT calculations without losing important conformers | Medium | Pre-DFT: 6 kcal/mol (broad). Post-DFT: 3 kcal/mol or 95% cumulative population. Saves compute on MMFF/xTB artifacts. |
| CREST optional (not required) | App works without system deps; graceful fallback to RDKit | Medium | Most tools require all deps. RDKit-only mode enables usage without xtb/crest installation. |
| Scaling factors per solvent | Existing v1.1 infrastructure extends naturally to ensemble | Low | Competitive advantage from v1.1 carries forward. Ensemble doesn't change this. |
| Interactive 3D lowest conformer | Shows most populated structure with shift labels | Low | Extends existing 3Dmol.js viewer. Ensemble jobs show lowest-energy conformer, not arbitrary. |

### Anti-Features (Commonly Requested, Often Problematic)

| Feature | Why Requested | Why Problematic | Alternative |
|---------|---------------|-----------------|-------------|
| Automatic rigidity detection | "Let the software decide if I need ensemble" | Rotatable bond count is crude (5-8 bonds = maybe flexible?). nConf20 descriptor better but adds complexity. Misclassification wastes compute or gives poor results. | **User explicitly chooses** single vs ensemble mode. Users know their molecules better than heuristics. |
| Per-conformer shift detail in API | "I want to see all individual conformers" | Overwhelming data (20 conformers × 30 atoms = 600 values). Rarely actionable—users compare to experiment (one spectrum). Complicates API design. | **Population metadata only**: conformer count, energy range, top populations. Weighted average is what matters. |
| Real-time progress updates | "Show me conformer 5/20 processing" | WebSocket complexity for marginal UX gain. Most jobs complete faster than user checks. Polling status works fine for v1.x. | **Standard polling**: Status endpoint returns "generating_conformers", "optimizing_conformers", "calculating_nmr", "complete". Sufficient granularity. |
| DP4+ probability scoring | "Give me confidence in structure assignment" | DP4+ requires experimental spectrum input. Out of scope for prediction-only tool. Adds statistical complexity. | **Future consideration (v3+)**: Structure validation mode with experimental input. V2.0 focuses on prediction. |
| Conformer geometry downloads | "Export all conformer XYZ files" | Large file downloads (20 conformers × 1-3 KB each). Storage overhead. Rarely needed—users want shifts, not geometries. | **Lowest-energy geometry only**: Standard for single-conformer jobs. Sufficient for visualization and downstream use. |
| Automatic method selection by molecule size | "Use CREST for small, RDKit for large" | Size is poor proxy for flexibility. Cyclohexane (18 atoms) needs ensemble. Rigid steroids (30+ atoms) don't. Logic would be wrong often. | **User choice of method**: RDKit (fast) or CREST (accurate). User knows flexibility better than atom count heuristic. |
| MMFF energy pre-filtering | "Filter conformers by MMFF before DFT" | MMFF energies anticorrelate with DFT (Spearman ρ ~ -0.1 to -0.45). Can discard important conformers. | **GFN2-xTB pre-filtering only**: If CREST used, xTB energies filter well (ρ ~ 0.39-0.47). If RDKit used, go straight to DFT. |

## Feature Dependencies

```
Existing v1.1 Features (Dependencies)
├── SMILES/SDF input handling → Conformer generation needs starting geometry
├── NWChem DFT optimization → Each conformer optimized with existing code
├── NWChem GIAO NMR → Each conformer's shieldings calculated with existing code
├── Scaling factors (DELTA50) → Applied to weighted-average shifts (not per-conformer)
├── Job queue (Huey) → Ensemble jobs queue N×DFT+N×NMR subtasks
├── Status polling → Extended with ensemble-specific states
└── 3Dmol.js viewer → Shows lowest-energy conformer geometry

New v2.0 Features (Internal Dependencies)
├── Conformer generation (RDKit KDG)
│   └── Requires: RDKit (existing dep), SMILES or starting geometry
├── Conformer generation (CREST/xTB)
│   ├── Requires: xtb + crest binaries on PATH (optional)
│   └── Fallback: RDKit if not detected
├── Energy window filtering
│   ├── Pre-DFT: Filter RDKit MMFF (crude) or CREST xTB (good) energies
│   └── Post-DFT: Filter by DFT energies from optimization step
├── Multi-conformer DFT optimization
│   └── Requires: Energy-filtered conformer list, existing NWChem opt logic
├── Multi-conformer NMR calculation
│   └── Requires: Optimized conformer geometries + DFT energies, existing NWChem NMR logic
├── Boltzmann weighting
│   ├── Requires: DFT energies (from optimization), shifts (from NMR calc), temperature parameter
│   └── Formula: Standard Boltzmann averaging over per-conformer shifts
└── Population metadata
    └── Requires: Boltzmann weights, conformer count, energy range
```

## MVP Definition

### Launch With (v2.0)

**Core conformational sampling workflow:**
1. User chooses single-conformer (v1.x behavior) or ensemble mode
2. User chooses RDKit (fast) or CREST (accurate) if ensemble mode
3. User sets energy window (default 6 kcal/mol pre-DFT, 3 kcal/mol post-DFT)
4. User sets temperature (default 298.15 K)
5. System generates conformers, filters by energy window
6. System runs DFT optimization on each surviving conformer
7. System re-filters by DFT energies (tighter window or 95% cumulative population)
8. System runs NMR calculation on each surviving conformer
9. System Boltzmann-weights shifts by DFT energies
10. API returns weighted-average shifts only (not per-conformer detail)

**Metadata in response:**
- Conformer count (initial, post-pre-filter, post-DFT-filter, final NMR count)
- Energy range (min/max relative energies in kcal/mol)
- Top 3 conformer populations (as percentages)
- Method used (RDKit KDG or CREST+xTB)
- Temperature used for Boltzmann weighting

**3D viewer:**
- Shows lowest-energy conformer geometry with shift labels (same as v1.x single-conformer mode)

**Status polling granularity:**
- "generating_conformers" (RDKit or CREST running)
- "filtering_conformers" (energy window filtering)
- "optimizing_conformers" (DFT geometry optimizations in progress, X/N complete)
- "calculating_nmr" (NMR shielding calculations in progress, X/N complete)
- "averaging_shifts" (Boltzmann weighting)
- "complete"
- "failed"

**Error handling:**
- CREST not found → graceful fallback to RDKit with warning in response metadata
- Zero conformers after filtering → error, suggest wider energy window
- DFT optimization fails on all conformers → error
- DFT optimization fails on some conformers → warning, continue with survivors

### Add After Validation (v2.x)

**Post-launch enhancements** (after v2.0 user feedback):

1. **Adaptive energy window** (v2.1)
   - Automatically tighten post-DFT filter based on 95% cumulative Boltzmann population
   - Current v2.0 plan: user sets window or uses default 3 kcal/mol
   - Enhancement: calculate cumulative population, keep conformers until 95% threshold

2. **Parallel conformer processing** (v2.1)
   - v2.0: sequential DFT optimizations (simpler implementation)
   - v2.1: parallel Huey tasks for independent conformer DFT jobs
   - Benefit: 20 conformers processed in ~15 min instead of 5 hours (for 15 min/conformer)

3. **Conformer caching** (v2.2)
   - Cache RDKit or CREST conformer ensembles by (SMILES, method, energy_window, temperature)
   - Resubmit same molecule with different NMR settings → skip conformer generation
   - Significant speedup for parameter sweeps or solvent comparison

4. **Rotatable bond count display** (v2.2)
   - Show in web UI: "This molecule has 6 rotatable bonds (moderately flexible)"
   - Informational only, not used for automatic decisions
   - Helps users decide single vs ensemble mode

5. **Entropy corrections** (v2.3)
   - Add GFN2-xTB quasi-RRHO entropy corrections to DFT energies
   - Improves Boltzmann weighting accuracy (free energy instead of electronic energy)
   - Requires xtb frequency calculation on DFT-optimized geometries

### Future Consideration (v3+)

**Major features requiring architectural changes or new scope:**

1. **DP4+ structure validation** (v3.0)
   - Accept experimental NMR spectrum as input
   - Calculate DP4+ probability for candidate structures
   - Requires: experimental spectrum parser, statistical framework, conformer comparison logic
   - Use case: Structure elucidation with confidence metrics

2. **Per-conformer detail export** (v3.0)
   - Optional API endpoint: `/jobs/{job_id}/conformers` returns full conformer detail
   - Includes: per-conformer geometries, energies, shifts, populations
   - Use case: Advanced users debugging ensemble or publishing conformational analysis
   - Not in v2.0: Adds API complexity, niche use case

3. **Real-time progress streaming** (v3.1)
   - WebSocket-based progress updates during conformer processing
   - Shows: current conformer ID, optimization convergence, energy updates
   - Requires: WebSocket infrastructure, more complex state management
   - Benefit: Better UX for long-running jobs (hours), but polling sufficient for v2.0

4. **Multi-temperature ensemble** (v3.1)
   - Calculate Boltzmann weights at multiple temperatures (e.g., 273 K, 298 K, 323 K)
   - Compare temperature-dependent conformer populations
   - Use case: Variable-temperature NMR experiments

5. **Automatic rigidity heuristics** (v3.2)
   - Calculate nConf20 or nTABS descriptor for flexibility quantification
   - Suggest (not enforce) ensemble vs single-conformer mode
   - Requires: Conformational complexity library, descriptor calculation
   - v2.0 decision: User explicit choice is sufficient

6. **Solvent-specific conformer generation** (v3.2)
   - Run CREST with ALPB solvent model matching NMR solvent
   - Currently: CREST runs in vacuum (or generic ALPB), DFT uses COSMO
   - Benefit: Solvent-dependent conformer populations (e.g., polar solvents stabilize different conformers)

## Real-World Tool Comparison

### ISiCLE (PNNL, Open Source)

**Conformer approach:**
- Molecular dynamics simulations for conformer generation
- Boltzmann-weighted by Gibbs free energy
- NWChem backend (same as qm-nmr-calc)
- Demonstrated on 80 methylcyclohexane conformers with excellent results

**What they got right:**
- Automated workflow (minimal user input)
- NWChem integration for NMR
- COSMO implicit solvent
- Boltzmann weighting by free energy

**What's different from v2.0 plan:**
- ISiCLE uses MD for conformers (slower, more thorough)
- qm-nmr-calc uses RDKit or CREST (faster, sufficient for drug-like molecules)
- ISiCLE focuses on automation; qm-nmr-calc gives user control (method choice, energy window)

### CENSO (Grimme Lab, Open Source)

**Conformer approach:**
- CREST metadynamics for exhaustive conformer search
- Multi-stage filtering: GFN2-xTB → DFT single-point → high-level DFT
- Command-line tool with Python API
- Generates input for ANMR (spectrum simulation)

**What they got right:**
- Best-in-class conformer sampling (CREST)
- Intelligent energy filtering (multi-stage)
- Transparent intermediate results

**What's different from v2.0 plan:**
- CENSO is command-line; qm-nmr-calc is web service
- CENSO requires xtb/crest; qm-nmr-calc has RDKit fallback
- CENSO uses ANMR for spectra; qm-nmr-calc directly returns shifts
- CENSO is expert tool; qm-nmr-calc targets broader users

### Spartan (Commercial, GUI)

**Conformer approach:**
- Force field conformational search (MMFF, SYBYL)
- Quantum mechanical refinement (DFT or semi-empirical)
- Boltzmann-weighted NMR spectra (1D and 2D COSY)
- GUI-based workflow with visualization

**What they got right:**
- Integrated workflow (search → optimize → NMR → average in one click)
- Boltzmann-weighted 1D and 2D NMR spectra
- Visual conformer distribution display
- Multi-step recipes (MM + QM)

**What's different from v2.0 plan:**
- Spartan is GUI-only; qm-nmr-calc is API-first with web UI
- Spartan shows per-conformer detail; qm-nmr-calc weighted average only (simplicity)
- Spartan commercial ($1000s); qm-nmr-calc open source

## Implementation Priorities by User Impact

### High-Impact, Must-Have (v2.0)
1. **Conformer generation choice** (RDKit or CREST) → Users need speed/accuracy tradeoff
2. **Boltzmann-weighted average** → Core feature, non-negotiable
3. **Energy window filtering** → Reduces compute cost without losing accuracy
4. **Single-conformer fallback** → Don't force ensemble on rigid molecules
5. **Temperature parameter** → Experimental condition matching

### Medium-Impact, Should-Have (v2.0)
6. **Population metadata** → Transparency without overwhelming detail
7. **Lowest-energy geometry** → Standard output, important for visualization
8. **Status granularity** → Users want to know progress on long jobs
9. **Graceful CREST fallback** → Works without system deps

### Low-Impact, Nice-to-Have (v2.x)
10. **Adaptive energy window** → Optimization, not critical
11. **Parallel conformer processing** → Speed improvement, not functionality
12. **Conformer caching** → Optimization for parameter sweeps
13. **Rotatable bond display** → Informational, not actionable

### Future/Niche (v3+)
14. **DP4+ probability** → Requires experimental spectrum input (scope change)
15. **Per-conformer detail export** → Advanced users only
16. **Real-time progress streaming** → Marginal UX improvement over polling
17. **Multi-temperature ensemble** → Variable-temperature NMR experiments (niche)
18. **Automatic rigidity heuristics** → User choice sufficient
19. **Solvent-specific conformers** → Accuracy improvement, added complexity

## Sources

### Academic Literature
- [ISiCLE NMR Framework](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-018-0305-8) — Yesiltepe et al. 2018, automated NWChem-based NMR with conformer Boltzmann weighting
- [Fully Automated Quantum NMR](https://onlinelibrary.wiley.com/doi/full/10.1002/anie.201708266) — Grimme et al. 2017, automated quantum chemistry NMR workflow
- [Conformational Search Importance for NMR](https://pmc.ncbi.nlm.nih.gov/articles/PMC9694776/) — Cuadrado et al. 2022, force field comparison showing MMFF poor correlation with DFT
- [DP4+ Structure Assignment](https://pubs.acs.org/doi/abs/10.1021/ja105035r) — Smith & Goodman 2010, DP4 probability for stereochemical assignment
- [DP4-AI Automation](https://pmc.ncbi.nlm.nih.gov/articles/PMC8152620/) — Automated NMR data analysis from spectrometer to structure

### Software Documentation
- [CENSO GitHub](https://github.com/grimme-lab/CENSO) — Commandline ENergetic SOrting conformer analysis tool
- [CENSO NMR Calculation Docs](https://xtb-docs.readthedocs.io/en/latest/CENSO_docs/censo_nmr.html) — Documentation for NMR spectrum calculation workflow
- [ISiCLE GitHub](https://github.com/pnnl/isicle) — In silico chemical library engine for property prediction
- [Spartan Conformational Search FAQ](https://downloads.wavefun.com/FAQ/Conformational_Searching.html) — Conformational search methods and energy profile questions
- [Generation of Gaussian Input from Spartan Conformers](https://protocols.scienceexchange.com/protocols/generation-of-gaussian-09-input-files-for-the-computation-of-1h-and-13c-nmr-chemical-shifts-of-structures-from-a-spartan-14-conformational-search) — Protocol for NMR calculation workflow

### Methodological Reviews
- [NMR Ensemble Determination](https://www.nature.com/articles/s41570-023-00494-x) — Ensemble determination by NMR data deconvolution
- [Molecular Flexibility Beyond Rotatable Bonds](https://pubs.acs.org/doi/10.1021/acs.jcim.6b00565) — nConf20 descriptor for conformational flexibility quantification
- [Torsion Angular Bin Strings](https://pubs.acs.org/doi/10.1021/acs.jcim.4c01513) — TABS method for molecular flexibility quantification
- [CREST Conformational Space Exploration](https://pubs.aip.org/aip/jcp/article/160/11/114110/3278084/CREST-A-program-for-the-exploration-of-low-energy) — Automated low-energy molecular chemical space exploration

### Computational Best Practices
- [Automated NMR Workflows for Non-Experts](https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/10.1002/mrc.5540) — Best practices for automated NMR data processing
- [Low-Cost Methods for Conformer Screening](https://pubs.acs.org/doi/10.1021/acs.jpca.4c01407) — GFN-FF, GFN2-xTB, DFTB3, HF-3c comparison
- [GEOM Energy-Annotated Conformers](https://www.nature.com/articles/s41597-022-01288-4) — Large-scale conformer dataset with energy annotations

---
*Feature research for: Conformational Sampling NMR*
*Research confidence: HIGH — verified with authoritative sources (academic literature, official software documentation)*
*Researched: 2026-01-26*
