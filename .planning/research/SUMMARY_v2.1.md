# v2.1 UI Redesign Research Summary

**Project:** qm-nmr-calc v2.1 UI Redesign
**Domain:** Scientific web application - NMR prediction results visualization
**Researched:** 2026-01-29
**Overall Confidence:** HIGH

## Executive Summary

The v2.1 UI redesign adds bento grid layouts and glassmorphism effects to three pages: Submit (form), Status (progress tracking), and Results (spectra, 3D viewer, tables). Research reveals this is achievable using pure CSS Grid with native backdrop-filter, requiring no build tools or frameworks. This aligns perfectly with the project's existing "no build step" philosophy (FastAPI + Jinja2 + vanilla JS + CDN-loaded libraries).

**Recommended approach:** Use CSS Grid with grid-template-areas for bento layouts, native backdrop-filter for glassmorphism, and CSS custom properties for design tokens. Multi-file CSS architecture using CSS Cascade Layers provides maintainability without preprocessors. BEM naming prevents specificity conflicts. Replace Pico CSS entirely with custom styles.

**Key risks:** Accessibility (WCAG contrast violations on glass surfaces), mobile performance (backdrop-filter is GPU-intensive), and WebGL conflicts (3Dmol.js canvas rendering with glass overlays). All risks are mitigated with specific prevention strategies: higher opacity backgrounds for contrast, reduced blur on mobile, and isolation of canvas elements from glass effects.

## Key Findings

### Recommended Stack

Pure CSS with modern native features eliminates complexity while delivering sophisticated visual design. No build tools, no frameworks, no preprocessors needed.

**Core technologies:**
- **CSS Grid with grid-template-areas**: Bento layouts with named regions — semantic, flexible, easy to maintain
- **Native backdrop-filter**: Glassmorphism with 92% browser support (96% with fallbacks) — no polyfills needed
- **CSS Custom Properties**: Design tokens and theming — replaces preprocessor variables
- **CSS Cascade Layers**: Architectural organization — explicit priority without specificity wars
- **BEM naming convention**: Component boundaries — prevents naming conflicts

**Browser compatibility:**
- CSS Grid: 96%+ support (since 2017)
- backdrop-filter: 92% native, 96% with fallbacks (requires -webkit- prefix for Safari)
- Container queries: 90%+ (optional enhancement)

**Key decision:** Replace Pico CSS entirely. Custom CSS is simpler than overriding framework defaults for asymmetric bento grids and precise glassmorphism control.

### Expected Features

#### Must Have (Table Stakes)

From bento grid and scientific UI patterns:

- **Visual hierarchy via card sizing**: 3-4 size variations (hero 2x2, large 2x1, medium 1x1, small)
- **Clear section boundaries**: Consistent padding, borders, shadows on glassmorphic cards
- **Responsive stacking**: Desktop complex grid → tablet 2-3 columns → mobile single column
- **Large visualization cards**: Spectra and 3D viewer need 50-60% viewport for usability
- **Data table formatting**: Right-aligned numbers, monospace fonts, proper scientific notation
- **Units and labels**: All scientific data shows units (ppm, eV, Hz) and axis labels
- **Metadata accessibility**: Calculation parameters in dedicated card, not buried
- **Loading states**: Per-card skeleton loaders for async API calls

#### Should Have (Differentiators)

Modern interactivity and usability enhancements:

- **Interactive hover states**: Subtle lift (4px translateY) with shadow enhancement
- **Micro-animations**: Card entry animations, smooth transitions (200-300ms ease-out)
- **Focus mode**: Expand single card to fullscreen for detailed spectra inspection
- **Dark mode support**: Light mode is primary (research target), but dark mode reduces eye strain
- **Accessibility enhancements**: Keyboard navigation, ARIA labels, screen reader support
- **Contextual tooltips**: Explain technical terms inline for less experienced users
- **Download affordances**: Clear export buttons on relevant cards (spectra, tables, 3D models)

#### Defer (Anti-Features / v2+)

Features that degrade usability or are inappropriate:

- **More than 9 visible cards**: Overwhelming, reduces scanability — use progressive disclosure instead
- **Auto-playing animations**: Distracting for scientific work — use hover/click triggers only
- **Bento grid on Submit page**: Linear workflow suits traditional forms, not asymmetric grid
- **Excessive color variation**: Scientific data already colorful — neutral card backgrounds only
- **Image-only spectra**: Users need CSV/JSON downloads, not just pretty pictures
- **Autoplay 3D animations**: Users want control — static pose with play button

### Architecture Approach

Multi-file CSS architecture provides maintainability and debugging clarity. HTTP/2 makes multiple files performant without concatenation.

**File structure:**
```
static/css/
├── layers.css              # Layer order declaration
├── reset.css               # Minimal reset
├── tokens.css              # Design tokens (colors, spacing, glass effects)
├── base.css                # Element defaults
├── layout.css              # Bento grid system
├── components/
│   ├── glassmorphic-card.css
│   ├── form.css
│   ├── status-indicator.css
│   ├── molecule-viewer.css
│   └── spectrum-display.css
└── utilities.css           # Override layer
```

**Major components:**
1. **Bento Grid System** (layout.css) — CSS Grid with responsive custom properties, grid-auto-flow: dense for packing, span utilities for asymmetry
2. **Glassmorphic Card** (components/glassmorphic-card.css) — Reusable .glass-card class with backdrop-filter, BEM elements for header/body/footer, modifiers for variants
3. **Design Token System** (tokens.css) — CSS custom properties for glass effects, spacing, colors, transitions, responsive breakpoints

**Key patterns:**
- CSS Cascade Layers (@layer reset, base, layout, components, utilities) for explicit priority
- BEM naming (Block__Element--Modifier) for component boundaries
- Responsive tokens (custom properties update at breakpoints)
- Progressive enhancement (@supports for backdrop-filter fallbacks)

**Bento grid sizing for NMR results:**
- Hero (2x2): Primary 1H spectrum or 3D viewer
- Large (2x1): Secondary 13C spectrum or chemical shift table
- Medium (1x1): Metadata, downloads, conformer selector
- Small (1x0.5): Quick stats, status badges

### Critical Pitfalls

Top risks identified from research with prevention strategies:

1. **WCAG Contrast Violations on Glass Surfaces**
   - **Risk:** Translucent backgrounds fail 4.5:1 contrast ratio for text
   - **Prevention:** Use 85-95% opacity (not 10-30%), add text shadows, test with WCAG tools, avoid glass on form inputs
   - **Impact:** Accessibility violations, unreadable text, legal compliance issues
   - **Phase concern:** Phase 2 (Results page) and Phase 3 (Submit page)

2. **Mobile Performance Degradation**
   - **Risk:** backdrop-filter is GPU-intensive, causes battery drain and stuttering
   - **Prevention:** Limit to 2-3 glass elements per viewport on mobile, reduce blur to 6px (vs 10px desktop), never animate backdrop-filter
   - **Impact:** Janky scrolling, device overheating, poor UX on tablets/phones
   - **Phase concern:** All phases, especially Phase 4 (Responsive) and Phase 6 (Testing)

3. **Z-Index and Stacking Context Chaos**
   - **Risk:** backdrop-filter creates new stacking context, breaks dropdown/modal z-index
   - **Prevention:** Document z-index scale (0-99 base, 100-199 glass, 200-299 dropdowns, 300+ modals), use portal pattern for modals
   - **Impact:** Dropdowns hidden behind glass, modals don't overlay properly
   - **Phase concern:** Phase 3 (Submit page forms) and Phase 5 (Interactivity)

4. **3Dmol.js Canvas Rendering Conflicts**
   - **Risk:** WebGL canvas + backdrop-filter compositing causes flicker, context loss
   - **Prevention:** Never apply backdrop-filter over canvas, use isolation: isolate on canvas container, test rotation/zoom with glass visible
   - **Impact:** 3D viewer flickers, performance drops, visual glitches
   - **Phase concern:** Phase 2 (Results page with 3D viewer)

5. **Bento Grid Content Overflow**
   - **Risk:** Long SMILES strings or molecule names break grid layout
   - **Prevention:** Set min-width: 0 on grid items, use text-overflow: ellipsis pattern, add word-break: break-all for SMILES
   - **Impact:** Grid misalignment, horizontal scrollbars, responsive breakpoints fail
   - **Phase concern:** Phase 2 (Results page) and Phase 4 (Responsive)

## Implications for Roadmap

Based on research, suggested phase structure for v2.1 UI redesign:

### Phase 1: CSS Foundation and Design System
**Rationale:** Establish foundation before UI work. Define design tokens, glass effect classes, bento grid system. Remove Pico CSS.

**Delivers:**
- Complete design token system (colors, spacing, glass effects, transitions)
- CSS Cascade Layers architecture (layers.css + reset.css + base.css)
- Reusable .glass-card component with variants
- Bento grid layout utilities (layout.css)
- Safari webkit prefix for backdrop-filter
- @supports fallbacks for older browsers

**Addresses:**
- Table stakes: consistent spacing, clear boundaries
- Pitfall 1: Define glass opacity ranges (85-95%) for contrast
- Pitfall 8: Include -webkit-backdrop-filter prefix
- Pitfall 10: Add solid background fallbacks

**Research needed:** None — well-documented CSS patterns

---

### Phase 2: Results Page Redesign
**Rationale:** Most complex page (spectra, 3D viewer, tables, metadata). Tests all architectural decisions early.

**Delivers:**
- Bento grid layout for Results page (12-column grid with named areas)
- 1H/13C spectrum cards (hero size, 2x2 or 2x1)
- 3D molecule viewer card (hero size, 2x2) with canvas isolation
- Chemical shift table cards (large size, 2x1, proper number formatting)
- Metadata card (medium size, 1x1)
- Download affordances (medium size, 1x1)

**Addresses:**
- Features: large visualization cards, data table formatting, units/labels
- Pitfall 4: Test with long SMILES strings and molecule names
- Pitfall 7: Verify 3Dmol.js canvas rendering with glass cards (z-index, performance)
- Pitfall 2: Performance test on mobile device

**Uses:**
- CSS Grid with grid-template-areas
- Glassmorphic card components
- Design tokens from Phase 1

**Research needed:** None — standard implementation of established patterns

---

### Phase 3: Submit Page Redesign
**Rationale:** Form UX requires solid backgrounds (no glass on inputs). Simpler layout than Results page.

**Delivers:**
- Bento grid layout for Submit page (simplified 2-column form layout)
- Form inputs with solid backgrounds (not glassmorphic)
- Strong focus indicators (3px outline, 2px offset)
- File upload dropzone with visual feedback
- Validation message display with proper contrast
- Clear disabled states (opacity 0.6, cursor not-allowed)

**Addresses:**
- Features: scannable layout, responsive stacking
- Pitfall 6: Form inputs with solid backgrounds, keyboard navigation
- Pitfall 1: Ensure form label contrast meets WCAG
- Pitfall 3: Test file upload preview z-index

**Anti-pattern avoided:** NOT using glassmorphism on form inputs (usability issue)

**Research needed:** None — standard form patterns with BEM components

---

### Phase 4: Status Page Redesign
**Rationale:** Simplest page (progress indicators, status cards). Low complexity.

**Delivers:**
- Status cards in auto-fit grid (grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)))
- Progress indicators with proper contrast on glass
- Loading states visible on glass backgrounds
- Polling update handling (avoid unnecessary re-renders)

**Addresses:**
- Features: loading states, clear section boundaries
- Pitfall 1: Verify progress indicator contrast

**Research needed:** None — straightforward grid layout

---

### Phase 5: Responsive Implementation
**Rationale:** Mobile-first breakpoints with reduced glass effects for performance.

**Delivers:**
- Mobile breakpoint (< 768px): single column, reduced blur (6px), 1-2 glass elements max
- Tablet breakpoint (768-1024px): 2-3 columns, simplified spans
- Desktop breakpoint (1024px+): full asymmetric bento grid
- Responsive design tokens (custom properties update at breakpoints)
- @media (hover: hover) for touch device optimization
- @media (prefers-reduced-motion) for accessibility

**Addresses:**
- Table stakes: responsive stacking
- Pitfall 2: Mobile performance optimization (reduced blur, fewer glass elements)
- Pitfall 5: Bento grid collapse behavior

**Research needed:** None — standard responsive patterns

---

### Phase 6: Accessibility and Polish
**Rationale:** Final pass for WCAG compliance, cross-browser testing, performance audit.

**Delivers:**
- WCAG contrast validation on all text elements (WebAIM analyzer)
- Keyboard navigation testing (tab order, focus indicators)
- Screen reader testing (ARIA labels, semantic HTML)
- Safari testing (webkit prefix verification)
- Performance audit (Lighthouse 90+ score)
- Battery drain test (30 minutes continuous use on Android mid-range)
- Cross-browser testing (Chrome, Firefox, Safari, Edge, mobile)

**Addresses:**
- Features: accessibility enhancements
- Pitfall 1: WCAG contrast violations (final verification)
- Pitfall 2: Mobile performance (audit and optimize)
- Pitfall 8: Safari webkit prefix verification

**Research needed:** None — standard testing protocols

---

### Phase Ordering Rationale

**Dependencies:**
- Phase 1 provides foundation for all subsequent phases (design system must exist first)
- Phase 2 tests most complex scenarios early (3Dmol.js canvas, multiple card types)
- Phases 3-4 can proceed in parallel after Phase 2
- Phase 5 wraps layout decisions from Phases 2-4
- Phase 6 final validation (can't test until implementation complete)

**Risk mitigation:**
- Testing 3Dmol.js conflicts early (Phase 2) avoids late-stage architecture changes
- Mobile performance considerations built into Phase 5 (not afterthought)
- Accessibility validation in Phase 6 (final pass, not ignored)

**Grouping logic:**
- Phase 1: Foundation (CSS architecture)
- Phases 2-4: Page-specific implementations (Results, Submit, Status)
- Phase 5: Cross-cutting concern (responsiveness)
- Phase 6: Quality assurance (testing, polish)

### Research Flags

**No phases need deeper research.** All patterns are well-documented:
- CSS Grid bento layouts: HIGH confidence from MDN, multiple 2026 tutorials
- Glassmorphism implementation: HIGH confidence from MDN, browser compatibility data
- Form accessibility: HIGH confidence from WCAG guidelines, established patterns
- Performance optimization: HIGH confidence from browser documentation, community best practices

**Standard patterns throughout:**
- CSS Grid with grid-template-areas (MDN official docs)
- BEM naming convention (getbem.com official spec)
- CSS Cascade Layers (Smashing Magazine, CSS-Tricks)
- backdrop-filter with fallbacks (MDN, Can I Use data)

## Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| Stack | HIGH | Pure CSS approach with 96%+ browser support, official MDN documentation |
| Features | MEDIUM-HIGH | Strong consensus on bento grid patterns, less specific guidance on scientific visualization |
| Architecture | HIGH | CSS Cascade Layers and BEM are well-documented, multiple authoritative sources |
| Pitfalls | HIGH | Glassmorphism accessibility and performance issues extensively documented (Nielsen Norman, Axess Lab, Can I Use) |

**Overall confidence:** HIGH

### Gaps to Address

Minor gaps that need validation during implementation:

1. **3Dmol.js canvas + glassmorphism interaction**: Theoretical understanding is solid (WebGL compositing, stacking contexts), but specific behavior needs testing in Phase 2. Mitigation strategy is clear (isolation, z-index management).

2. **Mobile performance thresholds**: Research recommends 2-3 glass elements max on mobile with 6px blur, but exact performance impact depends on device. Phase 6 testing will validate and adjust if needed.

3. **Scientific data display on glass backgrounds**: General guidance exists for data tables and visualization, but specific contrast requirements for NMR spectra annotations need testing. Phase 2 will validate with actual data.

All gaps have clear validation points in roadmap phases. No unknowns that block planning.

## Sources

### Stack Research (STACK_v2.1_ui_redesign.md)

**CSS Grid and Bento Layouts (HIGH confidence):**
- MDN Grid Template Areas (official docs)
- Speckyboy CSS Bento Grids (2026 tutorial)
- Codemotion Bento Layout Tutorial (step-by-step implementation)
- iamsteve Bento Layout Guide (practical examples)
- FreeCodeCamp Bento Grids (educational resource)

**Glassmorphism (HIGH confidence):**
- MDN backdrop-filter (official specification)
- Job Huntley Glassmorphism Guide 2026 (current trends)
- Inverness Design Studio Glassmorphism 2026 (design patterns)
- UX Pilot Glassmorphism Best Practices (usability focus)
- Glass UI Generator (interactive tool for values)

**Browser Support (HIGH confidence):**
- Can I Use backdrop-filter (96% with fallbacks)
- LambdaTest CSS Backdrop Filter (cross-browser compatibility)

**Performance and Animation (HIGH confidence):**
- Josh Comeau CSS Transitions (educational authority)
- design.dev CSS Transitions Guide (best practices)
- WebPeak Animation Trends 2026 (current standards)

**Accessibility (HIGH confidence):**
- Axess Lab Glassmorphism Accessibility (authoritative testing)
- New Target Glassmorphism Accessibility (WCAG focus)
- WebAIM Contrast Guidelines (official WCAG resource)

### Features Research (FEATURES_bento_grid_ui.md)

**Bento Grid Patterns (HIGH confidence):**
- Mockuuups Studio Best Examples 2026 (curated collection)
- SaaSFrame 43 UI Examples (practical patterns)
- SaaSFrame Bento Layout Trend (design analysis)
- Bentogrids.com Collection (gallery resource)
- UX Girl Bento Grid Patterns (design principles)

**Responsive Implementation (HIGH confidence):**
- Frontend Mentor Bento Grid (CSS layout example)
- Baltech AI Dashboard Bento Grids (data visualization focus)
- Tailwind CSS Bento Grids (official component library)

**Dashboard Patterns (MEDIUM confidence):**
- Muzli Dashboard Design 2026 (design inspiration)
- Pencil & Paper Dashboard UX Patterns (UX analysis)
- DesignRush Dashboard Principles (best practices)
- UXPin Dashboard Design 2025 (design guide)

**Interactive Patterns (MEDIUM confidence):**
- WriterDock UI Trends 2026 (trend analysis)
- Awwwards Interactive Bento Grid (award-winning examples)
- Figma Bento Animation Community (design resources)

**Scientific Visualization (MEDIUM confidence):**
- Medium Chemical Modelling Dashboard (molecular data focus)
- Mol* Molecular Viewer (WebGL viewer reference)
- Tableau Visual Best Practices (data visualization authority)

### Architecture Research (ARCHITECTURE.md)

**CSS Organization (HIGH confidence):**
- MDN Organizing CSS (official guide)
- Acodez CSS Best Practices 2026 (current standards)
- Nick Paolini Modern CSS Toolkit 2026 (comprehensive guide)

**CSS Cascade Layers (HIGH confidence):**
- CSS-Tricks Cascade Layers Guide (authoritative tutorial)
- Smashing Magazine Cascade Layers (detailed implementation)
- MDN Cascade Layers (official specification)

**BEM Naming (HIGH confidence):**
- getbem.com (official BEM specification)
- Sparkbox BEM by Example (practical guide)
- GeeksforGeeks BEM Convention (educational resource)

**Glassmorphism Implementation (HIGH confidence):**
- Kinetools Modern UI Glassmorphism 2026 (complete guide)
- Josh Comeau backdrop-filter Tutorial (advanced techniques)
- Bento Grid Generator (interactive tool)

### Pitfalls Research (PITFALLS_v2.1_UI_redesign.md)

**Glassmorphism Accessibility (HIGH confidence):**
- Axess Lab Glassmorphism Accessibility (comprehensive testing)
- UXPilot Glassmorphism Features (12 best practices)
- Nielsen Norman Group Glassmorphism (usability authority)
- New Target Glassmorphism Accessibility (WCAG compliance)

**Performance and Browser Compatibility (HIGH confidence):**
- Can I Use backdrop-filter (verified 96% support)
- MDN backdrop-filter (official specification)
- GitHub Xen-HTML Issue #219 (battery drain documentation)
- LambdaTest Compatibility Score (cross-browser data)
- Smashing Magazine Image Effects Performance (benchmark study)

**Bento Grid Implementation (MEDIUM confidence):**
- LogRocket Bento Grid UX (design analysis)
- Mockuuups Studio Examples (2026 patterns)
- ibelick Bento Grid Layouts (tutorial)

**CSS Grid Text Overflow (HIGH confidence):**
- TheLinuxCode Text Overflow Fix (troubleshooting guide)
- Mozilla Bugzilla Grid Cell Overflow (official bug tracker)
- FrontendTools CSS Grid 2025-2026 (modern techniques)

**WebGL and Canvas Performance (MEDIUM confidence):**
- Dev3lop WebGL vs Canvas Benchmarks (performance study)
- MDN Optimizing Canvas (official guide)
- Mozilla Bugzilla CSS Blur Performance (browser bug documentation)

---

**Research completed:** 2026-01-29
**Ready for roadmap:** Yes
**Recommended next step:** Proceed to requirements definition and roadmap creation
