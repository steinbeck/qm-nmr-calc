# Roadmap: qm-nmr-calc

## Milestones

- ✅ **v1.0 MVP** - Phases 1-6 (shipped 2026-01-20)
- ✅ **v1.1 Accurate Chemical Shifts** - Phases 7-11.2 (shipped 2026-01-25)
- ✅ **v2.0 Conformational Sampling** - Phases 12-17 (shipped 2026-01-28)
- ✅ **v2.0.1 Conformer Pre-selection** - Phase 24 (shipped 2026-01-30)
- ⏸️ **v2.1 UI Redesign** - Phases 18-23 (paused at Phase 21)

## Overview

v2.1 modernizes the web interface with bento grid layouts and glassmorphism effects. Pure CSS approach using native backdrop-filter and CSS Grid provides sophisticated visual design without build tools or frameworks. Delivers professional scientific aesthetic with clear component separation while maintaining accessibility and mobile performance.

The roadmap follows a foundation-first strategy: establish design system and reusable components (Phase 18), then apply to each page in complexity order (Results most complex, Status simplest), wrap with responsive implementation, and validate with accessibility testing. Critical decisions implemented early: design tokens prevent inconsistency, glass opacity ranges ensure WCAG contrast compliance, mobile performance budgets prevent GPU-intensive effects on devices.

## Phases

<details>
<summary>✅ v1.0 MVP (Phases 1-6) - SHIPPED 2026-01-20</summary>

### Phase 1: Project Setup and Foundation
**Goal**: Runnable service with basic architecture
**Plans**: 3 plans
**Status**: Complete

### Phase 2: NWChem Integration
**Goal**: Execute NWChem calculations via Huey
**Plans**: 3 plans
**Status**: Complete

### Phase 3: NMR Calculation Pipeline
**Goal**: End-to-end NMR predictions
**Plans**: 2 plans
**Status**: Complete

### Phase 4: Web Interface
**Goal**: Usable web UI for job submission
**Plans**: 2 plans
**Status**: Complete

### Phase 5: Result Visualization
**Goal**: Visual outputs (spectrum, structure)
**Plans**: 2 plans
**Status**: Complete

### Phase 6: Polish and Documentation
**Goal**: Production-ready v1.0
**Plans**: 2 plans
**Status**: Complete

</details>

<details>
<summary>✅ v1.1 Accurate Chemical Shifts (Phases 7-11.2) - SHIPPED 2026-01-25</summary>

### Phase 7: Remove ISiCLE Dependency
**Goal**: Direct NWChem integration
**Plans**: 2 plans
**Status**: Complete

### Phase 8: COSMO Solvation
**Goal**: Accurate solvent modeling
**Plans**: 2 plans
**Status**: Complete

### Phase 9: DELTA50 Benchmark
**Goal**: NWChem-specific scaling factors
**Plans**: 3 plans
**Status**: Complete

### Phase 10: Enhanced Visualization
**Goal**: Interactive 3D viewer with annotations
**Plans**: 2 plans
**Status**: Complete

### Phase 11: Geometry Input Options
**Goal**: Accept pre-optimized structures
**Plans**: 2 plans
**Status**: Complete

### Phase 11.1: Solvent Expansion Research
**Goal**: Research additional solvents
**Plans**: 1 plan
**Status**: Complete

### Phase 11.2: Vacuum Benchmark
**Goal**: Gas-phase scaling factors
**Plans**: 3 plans
**Status**: Complete

</details>

<details>
<summary>✅ v2.0 Conformational Sampling (Phases 12-17) - SHIPPED 2026-01-28</summary>

### Phase 12: Conformer Data Model and Storage
**Goal**: Foundation for multi-conformer calculations -- data structures, file organization, backward compatibility
**Depends on**: Phase 11.2
**Requirements**: DFT-03, API-04 (partial)
**Success Criteria** (what must be TRUE):
  1. System tracks conformer IDs, energies (with units), and geometries in structured format
  2. Each conformer uses isolated scratch directory to prevent NWChem file conflicts
  3. Canonical atom ordering established and verified across conformer lifecycle
  4. Existing v1.x single-conformer jobs load and run without modification
  5. Job directory structure supports per-conformer outputs (`output/conformers/`, `output/optimized/`)
**Plans**: 3 plans
**Status**: Complete

### Phase 13: RDKit KDG Conformer Generation
**Goal**: Generate conformer ensembles using RDKit distance geometry (no external dependencies)
**Depends on**: Phase 12
**Requirements**: CONF-03, CONF-06, FILT-01
**Success Criteria** (what must be TRUE):
  1. System generates conformers using KDG method (pure distance geometry, no crystal bias)
  2. MMFF optimization and RMSD deduplication produces diverse low-energy set
  3. Pre-DFT energy window filter (default 6 kcal/mol) removes high-energy conformers
  4. Adaptive conformer count scales with rotatable bond count
  5. App works fully without CREST/xTB installed (RDKit-only mode validated)
**Plans**: 3 plans
**Status**: Complete

### Phase 14: Boltzmann Averaging Implementation
**Goal**: Numerically stable Boltzmann weighting and ensemble averaging
**Depends on**: Phase 13
**Requirements**: BOLTZ-01, BOLTZ-02, BOLTZ-03, BOLTZ-04
**Success Criteria** (what must be TRUE):
  1. System calculates Boltzmann weights from DFT energies using exp-normalize trick (no overflow/underflow)
  2. System computes population-weighted average chemical shifts across conformer ensemble
  3. Temperature parameter (default 298.15 K) correctly applied in Boltzmann calculation
  4. Unit tests validate numerics with wide energy ranges (0-20 kcal/mol) and known test cases
  5. Single-conformer edge case (ensemble of 1) returns correct degeneracy (weight = 1.0)
**Plans**: 2 plans
**Status**: Complete

### Phase 15: Multi-Conformer NWChem Integration
**Goal**: Sequential DFT optimization and NMR calculation loop for conformer ensemble
**Depends on**: Phase 14
**Requirements**: DFT-01, DFT-02, DFT-04, DFT-05, FILT-02, FILT-03
**Success Criteria** (what must be TRUE):
  1. System runs DFT geometry optimization on each conformer surviving pre-DFT filter
  2. System extracts DFT energies from optimization step for Boltzmann weighting
  3. Post-DFT energy window filter (default 3 kcal/mol) applied before NMR calculations
  4. System runs NMR shielding calculation on each conformer surviving post-DFT filter
  5. Partial conformer failures handled gracefully (job continues with successful conformers if >50% succeed)
  6. User can configure energy window parameters (pre-DFT and post-DFT thresholds)
**Plans**: 3 plans
**Status**: Complete

### Phase 16: CREST Integration (Optional High-Accuracy Mode)
**Goal**: Production-quality conformer generation using CREST metadynamics (when available)
**Depends on**: Phase 15
**Requirements**: CONF-04, CONF-05
**Success Criteria** (what must be TRUE):
  1. System generates conformers using CREST/xTB when binaries available on PATH
  2. System auto-detects CREST/xTB availability and reports status in API health endpoint
  3. CREST includes ALPB solvation model matching job solvent parameter
  4. CREST timeout (default 7200s) with clear error message prevents hanging on macrocycles (fail-fast, no auto-fallback)
  5. Environment variables (OMP_STACKSIZE, GFORTRAN_UNBUFFERED_ALL) configured for stability
**Plans**: 3 plans
**Status**: Complete

### Phase 17: API Integration and Progress Tracking
**Goal**: User-facing ensemble mode with conformer selection, metadata, and progress updates
**Depends on**: Phase 16
**Requirements**: CONF-01, CONF-02, API-01, API-02, API-03, PROG-01, PROG-02
**Success Criteria** (what must be TRUE):
  1. User can choose between single-conformer mode (v1.x behavior) and ensemble mode when submitting job
  2. User can choose RDKit KDG or CREST as conformer generation method in ensemble mode
  3. API returns Boltzmann-weighted average shifts (not per-conformer detail)
  4. API response includes ensemble metadata (conformer count, energy range, top 3 populations, method used, temperature)
  5. 3D viewer shows lowest-energy conformer geometry with shift labels
  6. Job status includes ensemble-specific states (generating_conformers, optimizing_conformers X/N, calculating_nmr X/N, averaging_shifts)
  7. Web UI displays conformer progress for ensemble jobs
**Plans**: 5 plans
**Status**: Complete

</details>

### v2.1 UI Redesign (In Progress)

**Milestone Goal:** Modern bento grid layout with glassmorphism effects -- visually appealing interface with clear component separation, professional scientific aesthetic, and accessible design.

#### Phase 18: CSS Foundation and Design System
**Goal**: Reusable CSS architecture with design tokens, glass components, and bento grid utilities
**Depends on**: Phase 17
**Requirements**: CSS-01, CSS-02, CSS-03, CSS-04, CSS-05, LAYOUT-06, GLASS-01, GLASS-02, GLASS-03, GLASS-04, GLASS-05, GLASS-06
**Success Criteria** (what must be TRUE):
  1. Design token system defines all colors, spacing, glass effects, and transitions as CSS custom properties
  2. CSS Cascade Layers architecture (reset, base, layout, components, utilities) manages style priority
  3. Reusable .glass-card component with BEM variants delivers glassmorphism without duplication
  4. Bento grid system (CSS Grid with named areas) supports asymmetric card layouts
  5. Safari -webkit-backdrop-filter prefix and @supports fallbacks ensure cross-browser compatibility
  6. Pico CSS fully removed from base template, replaced with custom multi-file CSS
**Plans**: 4 plans

Plans:
- [x] 18-01-PLAN.md — CSS architecture foundation (layers, reset, tokens)
- [x] 18-02-PLAN.md — Base styles and bento grid layout system
- [x] 18-03-PLAN.md — Glass card component with BEM variants
- [x] 18-04-PLAN.md — Template integration and Pico CSS removal
**Status**: Complete

#### Phase 19: Results Page Redesign
**Goal**: Bento grid layout for results page with prominent visualizations and clear data organization
**Depends on**: Phase 18
**Requirements**: RESULTS-01, RESULTS-02, RESULTS-03, RESULTS-04, RESULTS-05, RESULTS-06, RESULTS-07
**Success Criteria** (what must be TRUE):
  1. 3D molecular viewer displayed in hero card position (largest, most prominent)
  2. 1H and 13C spectrum images displayed in medium-sized cards with clear labels
  3. Ensemble metadata (conformer count, populations) displayed in dedicated card for ensemble jobs
  4. Chemical shift tables organized in cards with right-aligned numbers and proper formatting
  5. Download buttons grouped in compact card with clear affordances
  6. Calculation details (solvent, functional, basis set) visible in metadata card
  7. Conformer selector integrated with 3D viewer card (ensemble jobs only)
**Plans**: 3 plans

Plans:
- [x] 19-01-PLAN.md — Create results-page.css with viewer-card and page-specific components
- [x] 19-02-PLAN.md — Update results.html template with bento grid structure
- [x] 19-03-PLAN.md — Visual verification checkpoint
**Status**: Complete

#### Phase 20: Submit Page Redesign
**Goal**: Clean form layout with logical grouping and solid backgrounds for usability
**Depends on**: Phase 18
**Requirements**: SUBMIT-01, SUBMIT-02, SUBMIT-03, SUBMIT-04, SUBMIT-05
**Success Criteria** (what must be TRUE):
  1. Form inputs organized in logical groups with clear visual hierarchy
  2. Form inputs use solid backgrounds (not glassmorphic) for usability and contrast
  3. Molecule preview area provides visual feedback for SMILES/file input
  4. Required fields clearly distinguished from optional fields with visual indicators
  5. Conformer mode and method options presented with clear descriptions
**Plans**: 2 plans

Plans:
- [x] 20-01-PLAN.md — Create submit-page.css with form layout components
- [x] 20-02-PLAN.md — Update template with two-column layout and SmilesDrawer preview
**Status**: Complete

#### Phase 21: Status Page Redesign
**Goal**: Job progress visualization with clear step tracking and status indicators
**Depends on**: Phase 18
**Requirements**: STATUS-01, STATUS-02, STATUS-03, STATUS-04, STATUS-05
**Success Criteria** (what must be TRUE):
  1. Job progress displayed with visual progress indicators (bars, percentages)
  2. Step tracking shows completed, current, and pending steps with clear states
  3. 3D molecule preview displayed during calculation (if geometry available)
  4. Ensemble progress (X/N conformers) clearly visible for ensemble jobs
  5. Error states displayed clearly if job fails with actionable messages
**Plans**: 2 plans

Plans:
- [x] 21-01-PLAN.md — Create status-page.css with step tracker and progress components
- [x] 21-02-PLAN.md — Update template with bento layout and JavaScript for step tracker
**Status**: Complete

#### Phase 22: Responsive and Layout Polish
**Goal**: Mobile-first breakpoints with performance-optimized glass effects and responsive card layouts
**Depends on**: Phase 19, Phase 20, Phase 21
**Requirements**: LAYOUT-01, LAYOUT-02, LAYOUT-03, LAYOUT-04, LAYOUT-05
**Success Criteria** (what must be TRUE):
  1. Page layouts use CSS Grid bento system with asymmetric card arrangements on desktop
  2. Cards support variable sizes (1x1, 2x1, 2x2, full-width) that adapt responsively
  3. Consistent gutters (16-24px) between all cards at all breakpoints
  4. Responsive breakpoints collapse complex grids: desktop asymmetric → tablet 2-3 columns → mobile single column
  5. Cards have smooth hover state transitions (transform, shadow) with reduced motion support
**Plans**: 2 plans

Plans:
- [x] 22-01-PLAN.md — Touch-safe hover, GPU shadow animation, reduced motion for cards
- [x] 22-02-PLAN.md — Status page reduced motion and visual verification
**Status**: Complete

#### Phase 23: Accessibility and Testing
**Goal**: WCAG compliance, cross-browser validation, and performance optimization
**Depends on**: Phase 22
**Requirements**: A11Y-01, A11Y-02, A11Y-03, A11Y-04, A11Y-05
**Success Criteria** (what must be TRUE):
  1. All text on glass backgrounds meets WCAG 4.5:1 contrast ratio minimum (validated with tools)
  2. Keyboard navigation works for all interactive elements with visible focus indicators
  3. Reduced motion support via prefers-reduced-motion media query disables animations
  4. Mobile performance optimized (2-3 glass elements max, 6px blur vs 10px desktop, Lighthouse 90+ score)
  5. Cross-browser testing passes on Chrome, Firefox, Safari, Edge with webkit prefix verification
**Plans**: 2 plans

Plans:
- [ ] 23-01-PLAN.md — Add focus tokens and comprehensive :focus-visible rules
- [ ] 23-02-PLAN.md — Validate accessibility, performance, and cross-browser compliance

<details>
<summary>✅ v2.0.1 Conformer Pre-selection Hotfix - SHIPPED 2026-01-30</summary>

**Hotfix Goal:** Reduce conformers sent to expensive DFT optimization from ~40 to ~8 using RMSD clustering and xTB energy ranking.

### Phase 24: Conformer Pre-selection with xTB
**Goal**: Efficient conformer filtering using RMSD clustering and xTB semi-empirical ranking
**Depends on**: Phase 17 (v2.0 complete)
**Requirements**: Performance optimization for ensemble calculations
**Success Criteria** (what must be TRUE):
  1. RMSD-based Butina clustering reduces redundant conformers (40 → ~8-12 clusters)
  2. xTB (GFN2) single-point energy ranking available when xTB binary detected
  3. System selects ~8 diverse, low-energy conformers for DFT optimization
  4. Fallback to MMFF ranking when xTB not available (RDKit-only mode works)
  5. Flexible molecules (hexanol, decane) complete in <3 hours instead of 10+ hours
**Plans**: 3 plans

Plans:
- [x] 24-01-PLAN.md — RMSD clustering with Butina algorithm (RDKit)
- [x] 24-02-PLAN.md — xTB integration for energy ranking
- [x] 24-03-PLAN.md — Pipeline integration and testing
**Status**: Complete

</details>

## Progress

**Execution Order:**
Phases execute in numeric order: 1 -> 6 (v1.0) -> 7 -> 11.2 (v1.1) -> 12 -> 17 (v2.0) -> 18 -> 23 (v2.1)

| Phase | Milestone | Plans Complete | Status | Completed |
|-------|-----------|----------------|--------|-----------|
| 1. Project Setup | v1.0 | 3/3 | Complete | 2026-01-20 |
| 2. NWChem Integration | v1.0 | 3/3 | Complete | 2026-01-20 |
| 3. NMR Pipeline | v1.0 | 2/2 | Complete | 2026-01-20 |
| 4. Web Interface | v1.0 | 2/2 | Complete | 2026-01-20 |
| 5. Visualization | v1.0 | 2/2 | Complete | 2026-01-20 |
| 6. Polish | v1.0 | 2/2 | Complete | 2026-01-20 |
| 7. Direct NWChem | v1.1 | 2/2 | Complete | 2026-01-25 |
| 8. COSMO Solvation | v1.1 | 2/2 | Complete | 2026-01-25 |
| 9. DELTA50 Benchmark | v1.1 | 3/3 | Complete | 2026-01-25 |
| 10. 3D Visualization | v1.1 | 2/2 | Complete | 2026-01-25 |
| 11. Geometry Input | v1.1 | 2/2 | Complete | 2026-01-25 |
| 11.1 Solvent Research | v1.1 | 1/1 | Complete | 2026-01-25 |
| 11.2 Vacuum Benchmark | v1.1 | 3/3 | Complete | 2026-01-25 |
| 12. Data Model | v2.0 | 3/3 | Complete | 2026-01-26 |
| 13. RDKit Generation | v2.0 | 3/3 | Complete | 2026-01-27 |
| 14. Boltzmann Averaging | v2.0 | 2/2 | Complete | 2026-01-27 |
| 15. NWChem Integration | v2.0 | 3/3 | Complete | 2026-01-27 |
| 16. CREST Integration | v2.0 | 3/3 | Complete | 2026-01-28 |
| 17. API Integration | v2.0 | 5/5 | Complete | 2026-01-28 |
| 18. CSS Foundation | v2.1 | 4/4 | Complete | 2026-01-29 |
| 19. Results Redesign | v2.1 | 3/3 | Complete | 2026-01-29 |
| 20. Submit Redesign | v2.1 | 2/2 | Complete | 2026-01-29 |
| 21. Status Redesign | v2.1 | 2/2 | Complete | 2026-01-31 |
| 22. Responsive Polish | v2.1 | 2/2 | Complete | 2026-01-31 |
| 23. Accessibility | v2.1 | 0/2 | Not started | - |
| 24. Conformer Pre-selection | v2.0.1 | 3/3 | Complete | 2026-01-30 |

---

*Last updated: 2026-01-31 after Phase 23 planning complete*
