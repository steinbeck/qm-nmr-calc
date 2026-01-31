# Phase 22: Responsive and Layout Polish - Research

**Researched:** 2026-01-31
**Domain:** CSS Grid responsive design, mobile-first breakpoints, accessibility-first animations
**Confidence:** HIGH

## Summary

Phase 22 focuses on responsive design patterns, mobile-first CSS breakpoints, and performance-optimized animations for the existing bento grid layout system. The project has already established a solid foundation with CSS Grid bento layouts, design tokens, and glass card components. This phase refines responsive behavior, ensures smooth cross-device experiences, and adds accessible hover interactions.

The research confirms that the existing approach (pure CSS, CSS Grid with `grid-template-areas`, mobile-first breakpoints at 768px and 1024px) aligns with 2026 best practices. Modern techniques include container queries (90%+ browser support for size queries), `clamp()` for fluid scaling, and `@media (hover: hover)` for touch-safe hover states.

Key findings show that CSS Grid's `repeat(auto-fill, minmax())` pattern can reduce media query complexity, `gap` property should use responsive scaling (via custom properties or `clamp()`), and GPU-accelerated properties (`transform`, `opacity`) should be preferred over CPU-bound properties (`box-shadow`) for animations.

**Primary recommendation:** Enhance existing bento grid system with fluid gap scaling, implement `@media (hover: hover)` for touch-safe interactive cards, and ensure all animations respect `prefers-reduced-motion` with opacity-based alternatives.

## Standard Stack

The established libraries/tools for this domain:

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| CSS Grid | Native CSS | Layout system for asymmetric bento grids | Built-in browser feature, no dependencies |
| CSS Custom Properties | Native CSS | Design tokens for responsive values | Enables dynamic theming and responsive scaling |
| CSS Cascade Layers | Native CSS | Architecture organization | Manages specificity without !important |
| Media Queries Level 4 | Native CSS | Responsive breakpoints and feature detection | Standard for `hover: hover`, `prefers-reduced-motion` |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| CSS Container Queries | Native CSS (Baseline 2023) | Component-responsive sizing | When components need to adapt to container size, not viewport |
| clamp() function | Native CSS | Fluid typography and spacing | For gap, padding, font-size that scales smoothly |
| prefers-reduced-motion | Native CSS | Accessibility-safe animations | All animations and transitions (required for WCAG 2.3.3) |
| @media (hover: hover) | Media Queries Level 4 | Touch-safe hover states | Prevents sticky hover on touch devices |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| CSS Grid | Flexbox | Grid better for 2D layouts; Flexbox for 1D alignment within grid cells |
| Media Queries | Container Queries | Media queries for page layout; Container queries for reusable components |
| Fixed breakpoints | `auto-fill` + `minmax()` | Auto-responsive grids reduce media queries but sacrifice explicit control |

**Installation:**
```bash
# No installation needed - all native CSS features
# Verify browser support: https://caniuse.com
```

## Architecture Patterns

### Recommended Project Structure
```
src/qm_nmr_calc/api/static/css/
├── tokens.css           # Design tokens with responsive overrides
├── layout.css           # Bento grid system with breakpoints
├── utilities.css        # Responsive utility classes
└── components/
    └── glass-card.css   # Card component with hover states
```

### Pattern 1: Mobile-First Responsive Grid
**What:** Start with single-column mobile layout, progressively enhance for tablet and desktop using `min-width` media queries.

**When to use:** Page-level layouts where screen size dictates structure changes.

**Example:**
```css
/* Source: MDN CSS Grid Layout Guide */
/* Mobile: Single column (default) */
.bento-grid {
  display: grid;
  grid-template-columns: 1fr;
  gap: var(--bento-gap);
}

/* Tablet: 3 columns */
@media (min-width: 768px) {
  .bento-grid {
    grid-template-columns: repeat(3, 1fr);
  }
}

/* Desktop: 6 columns */
@media (min-width: 1024px) {
  .bento-grid {
    grid-template-columns: repeat(6, 1fr);
  }
}
```

### Pattern 2: Fluid Responsive Gap with clamp()
**What:** Use `clamp()` to scale gap values smoothly between minimum and maximum sizes based on viewport width.

**When to use:** When you want spacing to adapt fluidly without hard breakpoint jumps.

**Example:**
```css
/* Source: Web.dev CSS clamp() guide */
:root {
  /* Scales from 12px (mobile) to 24px (desktop) */
  --bento-gap: clamp(0.75rem, 2.5vw, 1.5rem);
}

/* Alternative: Responsive custom property overrides */
:root {
  --bento-gap: 0.75rem; /* 12px mobile */
}

@media (min-width: 768px) {
  :root {
    --bento-gap: 1rem; /* 16px tablet */
  }
}

@media (min-width: 1024px) {
  :root {
    --bento-gap: 1.5rem; /* 24px desktop */
  }
}
```

### Pattern 3: Touch-Safe Hover States
**What:** Apply hover effects only on devices that support hover (mice, trackpads), not touch devices.

**When to use:** All interactive hover animations to prevent sticky hover states on mobile.

**Example:**
```css
/* Source: MDN Media Queries Level 4 */
/* No hover styles in base CSS */
.glass-card--interactive {
  transition: transform var(--transition-base),
              box-shadow var(--transition-base);
}

/* Apply hover only on hover-capable devices */
@media (hover: hover) and (pointer: fine) {
  .glass-card--interactive:hover {
    transform: translateY(-4px);
    box-shadow: 0 12px 48px hsl(0 0% 0% / 0.2);
  }
}

/* Alternative: Touch interaction with :active */
@media (hover: none) {
  .glass-card--interactive:active {
    transform: scale(0.98);
  }
}
```

### Pattern 4: GPU-Accelerated Animations
**What:** Animate only `transform` and `opacity` properties, which are GPU-accelerated. Avoid animating `box-shadow` directly.

**When to use:** All transitions and animations for smooth 60fps performance.

**Example:**
```css
/* Source: Performance best practices research */
/* GOOD: GPU-accelerated properties */
.card {
  transition: transform 250ms ease, opacity 250ms ease;
}

.card:hover {
  transform: translateY(-4px);
  opacity: 0.95;
}

/* BAD: CPU-bound, causes repaints */
.card {
  transition: box-shadow 250ms ease; /* Avoid! */
}

/* ALTERNATIVE: Pseudo-element shadow with opacity animation */
.card {
  position: relative;
}

.card::after {
  content: '';
  position: absolute;
  inset: 0;
  box-shadow: 0 12px 48px hsl(0 0% 0% / 0.2);
  opacity: 0;
  transition: opacity 250ms ease;
  pointer-events: none;
  z-index: -1;
}

.card:hover::after {
  opacity: 1;
}
```

### Pattern 5: Reduced Motion Accessible Animations
**What:** Provide alternative animations for users with `prefers-reduced-motion` setting. Replace transform-based motion with opacity changes.

**When to use:** All animations (required for WCAG 2.3.3 compliance).

**Example:**
```css
/* Source: MDN prefers-reduced-motion guide */
/* Default: Motion animation */
.status-badge[data-status="running"]::before {
  animation: pulse 1.5s ease-in-out infinite;
}

@keyframes pulse {
  0%, 100% { transform: scale(1); opacity: 1; }
  50% { transform: scale(1.1); opacity: 0.8; }
}

/* Reduced motion: Replace scale with opacity fade */
@media (prefers-reduced-motion: reduce) {
  .status-badge[data-status="running"]::before {
    animation: dissolve 1.5s ease-in-out infinite;
  }

  @keyframes dissolve {
    0%, 100% { opacity: 1; }
    50% { opacity: 0.5; }
  }
}
```

### Pattern 6: Asymmetric Bento Grid with Variable Spans
**What:** Use `grid-column: span N` and `grid-row: span N` to create cards of different sizes (1x1, 2x1, 2x2, full-width).

**When to use:** Dashboard-style layouts with visual hierarchy and asymmetric card arrangements.

**Example:**
```css
/* Source: Existing layout.css + bento grid research */
.bento-grid {
  display: grid;
  grid-template-columns: repeat(6, 1fr);
  gap: var(--bento-gap);
}

/* Variable card sizes */
.bento-grid__item--span-2 {
  grid-column: span 2; /* 1x1 card (2 of 6 columns) */
}

.bento-grid__item--span-3 {
  grid-column: span 3; /* 2x1 card (3 of 6 columns) */
}

.bento-grid__item--feature {
  grid-column: span 2; /* 2x2 card (2 cols x 2 rows) */
  grid-row: span 2;
}

.bento-grid__item--span-6 {
  grid-column: span 6; /* Full-width banner */
}

/* Responsive collapse */
@media (max-width: 767px) {
  .bento-grid {
    grid-template-columns: 1fr; /* Single column */
  }

  /* Reset all spans to single column */
  .bento-grid__item--span-2,
  .bento-grid__item--span-3,
  .bento-grid__item--span-6,
  .bento-grid__item--feature {
    grid-column: 1;
    grid-row: auto;
  }
}
```

### Pattern 7: Container Queries for Component Responsiveness (Optional Enhancement)
**What:** Use container queries for reusable components that adapt based on their container size, not viewport.

**When to use:** When components need to be truly reusable across different layout contexts (sidebar vs main area).

**Example:**
```css
/* Source: MDN Container Queries guide */
/* Define container */
.card-container {
  container-type: inline-size;
  container-name: card-wrapper;
}

/* Component responds to container width, not viewport */
@container card-wrapper (min-width: 400px) {
  .glass-card__title {
    font-size: var(--text-2xl);
  }

  .glass-card__body {
    display: grid;
    grid-template-columns: 1fr 1fr;
    gap: var(--space-md);
  }
}

@container card-wrapper (max-width: 399px) {
  .glass-card__title {
    font-size: var(--text-lg);
  }

  .glass-card__body {
    display: block;
  }
}
```

**Note:** Container queries are optional for this phase. They're best used when components need to adapt independently of viewport. For page-level layouts, media queries remain the standard.

### Anti-Patterns to Avoid
- **Desktop-first media queries:** Use `max-width` sparingly; prefer mobile-first `min-width` for smaller initial CSS payload.
- **Animating box-shadow directly:** CPU-bound, causes repaints. Use pseudo-element with opacity instead.
- **Sticky hover on touch:** Always wrap hover styles in `@media (hover: hover)` to prevent sticky states on mobile.
- **Disabling all animations for reduced motion:** Replace motion with opacity/fade instead of removing animations entirely.
- **Fixed pixel values everywhere:** Use custom properties and `clamp()` for responsive scaling.
- **grid-auto-flow: dense without consideration:** Can reorder items visually, breaking tab order and accessibility.

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Responsive breakpoints | Custom JavaScript resize listeners | CSS Media Queries | Built-in, performant, standard approach |
| Fluid typography/spacing | Manual calculations per breakpoint | `clamp(min, preferred, max)` | Single-line fluid scaling without breakpoints |
| Touch vs hover detection | JavaScript touch event detection | `@media (hover: hover)` | CSS-only, reliable, no JavaScript needed |
| Animation performance | Animating all properties | `transform` + `opacity` only | GPU-accelerated, smooth 60fps |
| Reduced motion support | JavaScript user preference detection | `@media (prefers-reduced-motion)` | Built-in accessibility, respects OS setting |
| Component responsiveness | Multiple media queries per component | Container Queries (optional) | Component adapts to container, not viewport |
| Grid gaps/gutters | Margins on grid items | CSS Grid `gap` property | Cleaner, no edge case margins |
| Responsive grids | Many media queries for columns | `repeat(auto-fill, minmax())` | Self-responsive, fewer breakpoints needed |

**Key insight:** Modern CSS provides native solutions for responsive design that are more performant, maintainable, and accessible than JavaScript-based alternatives. The existing codebase correctly uses pure CSS without frameworks or preprocessors.

## Common Pitfalls

### Pitfall 1: Sticky Hover States on Touch Devices
**What goes wrong:** Hover styles remain "stuck" after tapping on mobile devices, creating confusing UX where cards stay in hover state.

**Why it happens:** Touch events trigger `:hover` pseudo-class but don't have a natural "unhover" action. The hover state persists until another element is tapped.

**How to avoid:** Wrap all hover styles in `@media (hover: hover) and (pointer: fine)` to apply only on devices with hover capability (mice, trackpads).

**Warning signs:**
- Testing on mobile shows cards staying in hover state after tap
- Users report "stuck" or "highlighted" cards on touch devices
- Hover effects appear during scrolling on mobile

### Pitfall 2: Animating box-shadow Causes Performance Issues
**What goes wrong:** Animating `box-shadow` property directly causes janky animations, especially on mobile devices.

**Why it happens:** `box-shadow` is CPU-bound and triggers repaints on every frame. It's not GPU-accelerated like `transform` and `opacity`.

**How to avoid:** Use pseudo-element with `box-shadow` and animate its `opacity` instead, or use `filter: drop-shadow()` which may be GPU-accelerated.

**Warning signs:**
- Animations drop below 60fps in DevTools Performance panel
- Janky hover transitions on mobile devices
- Performance warnings in browser DevTools

### Pitfall 3: Removing All Animations for Reduced Motion
**What goes wrong:** Using `@media (prefers-reduced-motion: reduce)` to disable all animations removes helpful visual feedback.

**Why it happens:** Misunderstanding the purpose of `prefers-reduced-motion` as "no animations" instead of "reduced motion."

**How to avoid:** Replace motion-based animations (scale, translate, rotate) with opacity-based animations. Keep visual feedback, just make it less motion-intensive.

**Warning signs:**
- Users with reduced motion see no loading indicators
- No visual feedback for interactive elements
- Accessibility testing shows lack of state changes

### Pitfall 4: Hard Breakpoint Jumps in Spacing
**What goes wrong:** Gap values jump abruptly at breakpoints (12px → 24px), creating jarring layout shifts.

**Why it happens:** Using fixed pixel values in media queries instead of fluid scaling.

**How to avoid:** Use `clamp()` function for gap values: `gap: clamp(0.75rem, 2.5vw, 1.5rem)` scales smoothly without breakpoints.

**Warning signs:**
- Visual "pop" when resizing browser across breakpoints
- Layout feels jumpy rather than smooth
- Spacing looks wrong at edge of breakpoint ranges

### Pitfall 5: grid-auto-flow: dense Breaks Accessibility
**What goes wrong:** Using `grid-auto-flow: dense` to fill gaps causes visual order to differ from source order, breaking keyboard navigation.

**Why it happens:** Dense packing algorithm reorders items visually for optimal space usage, but tab order follows DOM order.

**How to avoid:** Only use `dense` when visual order matches source order, or when items are purely decorative. For interactive content, preserve source order.

**Warning signs:**
- Tab order jumps around visually when using keyboard
- Screen readers navigate in unexpected order
- WCAG 1.3.2 (Meaningful Sequence) failure

### Pitfall 6: Container Queries Without container-type
**What goes wrong:** Container queries don't apply because no ancestor has `container-type` declared.

**Why it happens:** Container queries require explicit containment context, unlike media queries which work globally.

**How to avoid:** Always declare `container-type: inline-size` (or `size`) on the parent container before using `@container` queries.

**Warning signs:**
- Container query styles never apply
- Browser DevTools show no matching container
- Component doesn't adapt to container width

### Pitfall 7: Over-Reliance on Container Queries
**What goes wrong:** Using container queries for everything when media queries are more appropriate, leading to performance overhead.

**Why it happens:** Excitement about new technology without understanding when to use it.

**How to avoid:** Use media queries for page-level layouts; container queries for reusable components that need context-independence. Don't replace all media queries.

**Warning signs:**
- Every component has container query styles
- Performance issues on lower-end devices
- Code is more complex than needed

## Code Examples

Verified patterns from official sources:

### Example 1: Complete Responsive Bento Grid
```css
/* Source: Existing layout.css + MDN Grid Layout patterns */
/* Mobile-first responsive grid system */

/* Mobile: Single column */
.bento-grid {
  display: grid;
  grid-template-columns: 1fr;
  gap: var(--bento-gap);
  width: 100%;
}

.bento-grid__item {
  min-height: var(--bento-min-height);
  display: flex;
  flex-direction: column;
}

/* Tablet: 3 columns (768px - 1024px) */
@media (min-width: 768px) {
  .bento-grid {
    grid-template-columns: repeat(3, 1fr);
  }

  .bento-grid__item--span-2 {
    grid-column: span 2;
  }

  .bento-grid__item--span-3 {
    grid-column: span 3; /* Full width on tablet */
  }

  .bento-grid__item--feature {
    grid-column: span 2;
    grid-row: span 2;
  }
}

/* Desktop: 6 columns (1024px+) */
@media (min-width: 1024px) {
  .bento-grid {
    grid-template-columns: repeat(6, 1fr);
  }

  .bento-grid__item--span-2 {
    grid-column: span 2;
  }

  .bento-grid__item--span-3 {
    grid-column: span 3;
  }

  .bento-grid__item--span-4 {
    grid-column: span 4;
  }

  .bento-grid__item--span-6 {
    grid-column: span 6;
  }
}
```

### Example 2: Responsive Design Tokens
```css
/* Source: Existing tokens.css + clamp() best practices */
:root {
  /* Spacing scale - mobile-first */
  --space-xs: 0.25rem;    /* 4px */
  --space-sm: 0.5rem;     /* 8px */
  --space-md: 1rem;       /* 16px */
  --space-lg: 1.5rem;     /* 24px */
  --space-xl: 2rem;       /* 32px */
  --space-2xl: 3rem;      /* 48px */

  /* Bento grid - responsive gaps */
  --bento-gap: 0.75rem;   /* 12px mobile */
  --bento-min-height: 120px;

  /* Glass effects - mobile performance */
  --glass-blur-light: blur(6px);
  --glass-blur-medium: blur(8px);
  --glass-blur-heavy: blur(10px);

  /* Transitions */
  --transition-fast: 150ms ease;
  --transition-base: 250ms ease;
  --transition-slow: 350ms ease;
}

/* Tablet overrides */
@media (min-width: 768px) {
  :root {
    --bento-gap: 1rem;  /* 16px tablet */
    --glass-blur-light: blur(8px);
    --glass-blur-medium: blur(12px);
  }
}

/* Desktop overrides */
@media (min-width: 1024px) {
  :root {
    --bento-gap: 1.5rem;  /* 24px desktop */
    --glass-blur-medium: blur(12px);
    --glass-blur-heavy: blur(16px);
  }
}

/* Alternative: Fluid gap with clamp() */
:root {
  /* Scales smoothly from 12px (mobile) to 24px (desktop) */
  --bento-gap-fluid: clamp(0.75rem, 2.5vw, 1.5rem);
}
```

### Example 3: Touch-Safe Interactive Cards
```css
/* Source: MDN Media Queries + existing glass-card.css */
.glass-card--interactive {
  /* Base styles for all devices */
  background: var(--glass-bg-light);
  backdrop-filter: var(--glass-blur-medium);
  -webkit-backdrop-filter: blur(12px);
  border-radius: var(--glass-radius);
  padding: var(--space-lg);

  /* Performance-optimized transitions */
  transition: transform var(--transition-base),
              box-shadow var(--transition-base);
}

/* Hover effects ONLY on devices with hover capability */
@media (hover: hover) and (pointer: fine) {
  .glass-card--interactive:hover {
    transform: translateY(-4px);
    box-shadow: 0 12px 48px hsl(0 0% 0% / 0.2);
    cursor: pointer;
  }

  .glass-card--interactive:active {
    transform: translateY(-2px);
  }
}

/* Touch devices: Subtle press feedback */
@media (hover: none) {
  .glass-card--interactive:active {
    transform: scale(0.98);
    opacity: 0.95;
  }
}

/* Reduced motion: Remove transform, keep shadow change */
@media (prefers-reduced-motion: reduce) {
  .glass-card--interactive {
    transition: box-shadow var(--transition-base);
  }

  @media (hover: hover) {
    .glass-card--interactive:hover {
      transform: none;
      box-shadow: 0 12px 48px hsl(0 0% 0% / 0.2);
    }
  }
}
```

### Example 4: GPU-Accelerated Shadow Animation
```css
/* Source: Performance best practices research */
/* BAD: Direct box-shadow animation (CPU-bound) */
.card-bad {
  box-shadow: 0 4px 16px hsl(0 0% 0% / 0.1);
  transition: box-shadow 250ms ease;
}

.card-bad:hover {
  box-shadow: 0 12px 48px hsl(0 0% 0% / 0.2);
}

/* GOOD: Pseudo-element shadow with opacity animation (GPU-accelerated) */
.card-good {
  position: relative;
  box-shadow: 0 4px 16px hsl(0 0% 0% / 0.1);
}

.card-good::after {
  content: '';
  position: absolute;
  inset: 0;
  box-shadow: 0 12px 48px hsl(0 0% 0% / 0.2);
  opacity: 0;
  transition: opacity 250ms ease;
  pointer-events: none;
  z-index: -1;
  border-radius: inherit;
}

.card-good:hover::after {
  opacity: 1;
}
```

### Example 5: Accessible Animations with Reduced Motion
```css
/* Source: MDN prefers-reduced-motion + existing status-page.css */
/* Default: Pulsing animation with scale */
.status-badge[data-status="running"]::before {
  content: '';
  width: 8px;
  height: 8px;
  border-radius: 50%;
  background: currentColor;
  animation: pulse 1.5s ease-in-out infinite;
}

@keyframes pulse {
  0%, 100% {
    transform: scale(1);
    opacity: 1;
  }
  50% {
    transform: scale(1.1);
    opacity: 0.8;
  }
}

/* Reduced motion: Replace scale with opacity fade */
@media (prefers-reduced-motion: reduce) {
  .status-badge[data-status="running"]::before {
    animation: dissolve 2s ease-in-out infinite;
  }

  @keyframes dissolve {
    0%, 100% {
      opacity: 1;
    }
    50% {
      opacity: 0.5;
    }
  }
}
```

### Example 6: Responsive Step Tracker
```css
/* Source: Existing status-page.css */
/* Desktop: Horizontal layout */
.step-tracker__list {
  display: flex;
  justify-content: space-between;
  list-style: none;
  padding: 0;
  margin: 0;
  position: relative;
}

/* Connecting line */
.step-tracker__list::before {
  content: '';
  position: absolute;
  top: 12px;
  left: 24px;
  right: 24px;
  height: 2px;
  background: var(--color-border);
}

.step-tracker__item {
  display: flex;
  flex-direction: column;
  align-items: center;
  gap: var(--space-xs);
  text-align: center;
  min-width: 80px;
}

/* Mobile: Vertical layout */
@media (max-width: 767px) {
  .step-tracker__list {
    flex-direction: column;
    align-items: flex-start;
    gap: var(--space-md);
  }

  /* Reposition connecting line for vertical layout */
  .step-tracker__list::before {
    top: 0;
    bottom: 0;
    left: 11px;
    right: auto;
    width: 2px;
    height: auto;
  }

  .step-tracker__item {
    flex-direction: row;
    text-align: left;
    gap: var(--space-sm);
    min-width: auto;
  }
}
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Desktop-first max-width queries | Mobile-first min-width queries | ~2015 | Smaller initial CSS payload, progressive enhancement |
| Fixed pixel breakpoints | Content-driven breakpoints + clamp() | 2023-2024 | Smoother responsive scaling, fewer hard jumps |
| JavaScript resize listeners | Media queries + ResizeObserver | 2020 | Better performance, declarative CSS |
| Media queries only | Media queries + Container queries | 2023 (Baseline) | Component-level responsiveness, better reusability |
| :hover without media query | @media (hover: hover) | 2021-2022 | Touch-safe hover states, no sticky hover |
| All animations disabled for reduced motion | Motion replaced with opacity | 2019-2020 | Better accessibility, maintains visual feedback |
| Animating all properties | Transform + opacity only | Ongoing | GPU acceleration, 60fps performance |
| calc() for fluid sizing | clamp(min, preferred, max) | 2020 | Simpler syntax, built-in constraints |

**Deprecated/outdated:**
- **Desktop-first @media (max-width)**: Mobile-first is standard; max-width creates larger initial CSS
- **JavaScript-based responsive behavior**: CSS media queries and container queries are performant and declarative
- **Polyfills for CSS Grid**: 95%+ browser support; no need for Autoprefixer grid polyfills
- **Hover without hover capability check**: Always wrap in `@media (hover: hover)` for touch devices
- **Removing all animations for reduced motion**: Replace motion with opacity instead

## Open Questions

Things that couldn't be fully resolved:

1. **Container queries adoption timeline**
   - What we know: 90%+ browser support for size queries (Baseline 2023); style queries partially supported
   - What's unclear: Whether style queries for custom properties will be practical in this project's timeframe (Firefox support pending)
   - Recommendation: Use container size queries if component reusability is critical, but media queries are sufficient for page-level layouts. Container queries are optional enhancement, not required for this phase.

2. **Optimal gap scaling approach: clamp() vs media queries**
   - What we know: `clamp()` provides fluid scaling; media queries provide explicit control
   - What's unclear: Whether smooth fluid scaling or explicit breakpoint control is better for this design system
   - Recommendation: Start with media query approach (already established in tokens.css), optionally add `clamp()` if design benefits from smoother transitions. Test both approaches with real content.

3. **Performance budget for mobile glass effects**
   - What we know: Prior decision is "2-3 glass elements max" on mobile with reduced blur (6px vs 10px desktop)
   - What's unclear: Exact device testing shows this performs at 60fps on target devices (needs validation)
   - Recommendation: Validate performance on low-end Android devices (target: 60fps with 2-3 backdrop-filter elements). Use DevTools Performance panel to measure.

## Sources

### Primary (HIGH confidence)
- [MDN: Realizing common layouts using grids](https://developer.mozilla.org/en-US/docs/Web/CSS/CSS_grid_layout/Realizing_common_layouts_using_grids) - CSS Grid responsive patterns
- [MDN: prefers-reduced-motion](https://developer.mozilla.org/en-US/docs/Web/CSS/Reference/At-rules/@media/prefers-reduced-motion) - Accessibility guidelines
- [MDN: CSS container queries](https://developer.mozilla.org/en-US/docs/Web/CSS/Guides/Containment/Container_queries) - Container query usage and limitations
- [W3C WCAG: C39 prefers-reduced-motion technique](https://www.w3.org/WAI/WCAG21/Techniques/css/C39) - Official accessibility guidance
- Existing codebase: layout.css, tokens.css, glass-card.css, status-page.css, submit-page.css

### Secondary (MEDIUM confidence)
- [LogRocket: Container queries in 2026](https://blog.logrocket.com/container-queries-2026/) - Current state and limitations
- [Can I Use: CSS Container Queries (Size)](https://caniuse.com/css-container-queries) - Browser support data
- [Web.dev: CSS min(), max(), and clamp()](https://web.dev/articles/min-max-clamp) - Fluid sizing techniques
- [DEV Community: Responsive Design Breakpoints 2025 Playbook](https://dev.to/gerryleonugroho/responsive-design-breakpoints-2025-playbook-53ih) - Modern breakpoint strategies
- [BrowserStack: Responsive Design Breakpoints](https://www.browserstack.com/guide/responsive-design-breakpoints) - Industry standards
- [iamsteve: Build a bento layout with CSS grid](https://iamsteve.me/blog/bento-layout-css-grid) - Bento grid implementation patterns
- [CSS-Tricks: prefers-reduced-motion](https://css-tricks.com/almanac/rules/m/media/prefers-reduced-motion/) - Implementation examples
- [Tobias Ahlin: How to animate box-shadow](https://tobiasahlin.com/blog/how-to-animate-box-shadow/) - Performance optimization patterns

### Tertiary (LOW confidence - WebSearch findings)
- [Medium: Handle :hover CSS on mobile touch screen](https://arturocreates.medium.com/handle-hover-css-on-mobile-touch-screen-69142ea79fe7) - Touch device patterns
- [DEV Community: Costly CSS Properties](https://dev.to/leduc1901/costly-css-properties-and-how-to-optimize-them-3bmd) - Performance considerations
- [Medium: How box-shadow impacts performance](https://medium.com/@emadfanaeian/how-box-shadow-and-transition-impact-performance-in-hybrid-mobile-apps-b5973087a4b8) - Mobile performance

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - Native CSS features with excellent documentation and browser support
- Architecture patterns: HIGH - Verified with MDN official docs, existing codebase follows best practices
- Performance optimization: HIGH - GPU acceleration patterns well-documented, verified with multiple sources
- Accessibility (reduced motion): HIGH - WCAG official guidance, MDN documentation
- Container queries: MEDIUM - Good browser support but style queries partially supported; optional for this phase
- Pitfalls: HIGH - Common issues well-documented across multiple authoritative sources

**Research date:** 2026-01-31
**Valid until:** 2026-02-28 (30 days - CSS standards are stable)
