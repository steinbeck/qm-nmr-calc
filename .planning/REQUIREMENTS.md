# Requirements: qm-nmr-calc v2.1

**Defined:** 2026-01-29
**Core Value:** Reliable async NMR predictions with full control over calculation parameters -- submit a molecule, get back accurate 1H/13C shifts without babysitting long-running calculations.

## v2.1 Requirements

Requirements for v2.1 UI Redesign. Each maps to roadmap phases.

### Layout System

- [ ] **LAYOUT-01**: Page layouts use CSS Grid bento grid system with asymmetric card arrangements
- [ ] **LAYOUT-02**: Cards support variable sizes (1x1, 2x1, 2x2, full-width)
- [ ] **LAYOUT-03**: Consistent gutters (16-24px) between all cards
- [ ] **LAYOUT-04**: Responsive breakpoints adapt layout for desktop, tablet, and mobile
- [ ] **LAYOUT-05**: Cards have smooth hover state transitions (transform, shadow)
- [ ] **LAYOUT-06**: Design tokens (CSS custom properties) for colors, spacing, and effects

### Glassmorphism

- [ ] **GLASS-01**: Glass cards use backdrop-filter blur effect (8-12px)
- [ ] **GLASS-02**: Cards have semi-transparent backgrounds (85-95% opacity for accessibility)
- [ ] **GLASS-03**: Subtle border highlights on glass cards (1px light border)
- [ ] **GLASS-04**: Safari compatibility with -webkit-backdrop-filter prefix
- [ ] **GLASS-05**: Layered depth with multiple glass intensity levels
- [ ] **GLASS-06**: Subtle box shadows for card elevation

### Results Page

- [ ] **RESULTS-01**: 3D molecular viewer in hero card position (largest, prominent)
- [ ] **RESULTS-02**: 1H and 13C spectrum images in medium cards
- [ ] **RESULTS-03**: Ensemble metadata displayed in dedicated card (for ensemble jobs)
- [ ] **RESULTS-04**: Chemical shift tables in organized card layout
- [ ] **RESULTS-05**: Download buttons grouped in compact card
- [ ] **RESULTS-06**: Calculation details (solvent, functional, basis set) visible
- [ ] **RESULTS-07**: Conformer selector integrated with 3D viewer card

### Submit Page

- [ ] **SUBMIT-01**: Clean form layout with logical grouping
- [ ] **SUBMIT-02**: Form inputs use solid backgrounds (not glass) for usability
- [ ] **SUBMIT-03**: Molecule preview area for visual feedback
- [ ] **SUBMIT-04**: Clear visual hierarchy for required vs optional fields
- [ ] **SUBMIT-05**: Conformer mode and method options clearly presented

### Status Page

- [ ] **STATUS-01**: Job progress displayed with visual indicators
- [ ] **STATUS-02**: Step tracking shows completed and current steps
- [ ] **STATUS-03**: 3D molecule preview during calculation
- [ ] **STATUS-04**: Ensemble progress (X/N conformers) clearly visible
- [ ] **STATUS-05**: Error states displayed clearly if job fails

### Accessibility & Performance

- [ ] **A11Y-01**: All text meets WCAG 4.5:1 contrast ratio minimum
- [ ] **A11Y-02**: Keyboard navigation works for all interactive elements
- [ ] **A11Y-03**: Reduced motion support via prefers-reduced-motion media query
- [ ] **A11Y-04**: Mobile performance optimized (reduced blur, limited glass elements)
- [ ] **A11Y-05**: Focus indicators visible for keyboard users

### CSS Architecture

- [ ] **CSS-01**: Pico CSS replaced with custom CSS framework
- [ ] **CSS-02**: CSS Cascade Layers for style priority management
- [ ] **CSS-03**: BEM naming convention for component classes
- [ ] **CSS-04**: Multi-file CSS organization (no build step required)
- [ ] **CSS-05**: Base template updated with new stylesheet loading

## v2.x Requirements

Deferred to future release. Tracked but not in current roadmap.

### Enhanced Interactivity

- **INTERACT-01**: Card expansion to focus mode (full-screen view)
- **INTERACT-02**: Drag-and-drop card reordering
- **INTERACT-03**: User preference persistence (card states)

### Dark Mode

- **DARK-01**: Dark mode color scheme
- **DARK-02**: System preference detection (prefers-color-scheme)
- **DARK-03**: Manual toggle between light and dark modes

## Out of Scope

Explicitly excluded. Documented to prevent scope creep.

| Feature | Reason |
|---------|--------|
| JavaScript framework (React, Vue) | Keeps stack simple, Jinja2 templates work well |
| CSS preprocessor (Sass, Less) | Modern CSS has variables and nesting, no build step needed |
| Animation library | CSS transitions sufficient for hover effects |
| Dark mode | Adds complexity, defer to v2.2 |
| Card drag-and-drop | Nice-to-have, not essential for v2.1 |

## Traceability

Which phases cover which requirements. Updated during roadmap creation.

| Requirement | Phase | Status |
|-------------|-------|--------|
| LAYOUT-01 | TBD | Pending |
| LAYOUT-02 | TBD | Pending |
| LAYOUT-03 | TBD | Pending |
| LAYOUT-04 | TBD | Pending |
| LAYOUT-05 | TBD | Pending |
| LAYOUT-06 | TBD | Pending |
| GLASS-01 | TBD | Pending |
| GLASS-02 | TBD | Pending |
| GLASS-03 | TBD | Pending |
| GLASS-04 | TBD | Pending |
| GLASS-05 | TBD | Pending |
| GLASS-06 | TBD | Pending |
| RESULTS-01 | TBD | Pending |
| RESULTS-02 | TBD | Pending |
| RESULTS-03 | TBD | Pending |
| RESULTS-04 | TBD | Pending |
| RESULTS-05 | TBD | Pending |
| RESULTS-06 | TBD | Pending |
| RESULTS-07 | TBD | Pending |
| SUBMIT-01 | TBD | Pending |
| SUBMIT-02 | TBD | Pending |
| SUBMIT-03 | TBD | Pending |
| SUBMIT-04 | TBD | Pending |
| SUBMIT-05 | TBD | Pending |
| STATUS-01 | TBD | Pending |
| STATUS-02 | TBD | Pending |
| STATUS-03 | TBD | Pending |
| STATUS-04 | TBD | Pending |
| STATUS-05 | TBD | Pending |
| A11Y-01 | TBD | Pending |
| A11Y-02 | TBD | Pending |
| A11Y-03 | TBD | Pending |
| A11Y-04 | TBD | Pending |
| A11Y-05 | TBD | Pending |
| CSS-01 | TBD | Pending |
| CSS-02 | TBD | Pending |
| CSS-03 | TBD | Pending |
| CSS-04 | TBD | Pending |
| CSS-05 | TBD | Pending |

**Coverage:**
- v2.1 requirements: 35 total
- Mapped to phases: 0 (pending roadmap)
- Unmapped: 35

---
*Requirements defined: 2026-01-29*
*Last updated: 2026-01-29 after initial definition*
