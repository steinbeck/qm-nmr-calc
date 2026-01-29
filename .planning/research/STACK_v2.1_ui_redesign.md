# Stack Research: Bento Grid UI with Glassmorphism

**Project:** qm-nmr-calc v2.1 UI Redesign
**Researched:** 2026-01-29
**Confidence:** HIGH

## Executive Summary

For the bento grid layout with glassmorphism effects, **use pure CSS Grid with modern native CSS features**. No build tools, no frameworks, no preprocessors. The existing FastAPI/Jinja2/vanilla JS stack supports this perfectly.

**Key approach:**
- CSS Grid with `grid-template-areas` for named layout regions
- Native `backdrop-filter: blur()` for glassmorphism (92% browser support)
- CSS custom properties for theming and reusability
- Transform/opacity for GPU-accelerated animations
- Progressive enhancement with `@supports` fallbacks

This maintains the project's "no build step" philosophy while delivering modern visual design.

---

## Recommended CSS Approach

### Core Stack Decision

**Use:** Pure CSS Grid + CSS Custom Properties + Native backdrop-filter
**Why:**
- Zero build complexity (aligns with existing 3Dmol.js CDN approach)
- Excellent browser support (all modern browsers since 2023-2024)
- GPU-accelerated performance
- Maintainable without preprocessor knowledge

**Replace:** Pico CSS (currently used)
**Why:**
- Bento grids need custom asymmetric layouts (Pico's semantic containers fight this)
- Glassmorphism requires precise control over backgrounds and blur
- Custom CSS is simpler than overriding framework defaults

---

## CSS Grid for Bento Layouts

### Grid-Template-Areas (Recommended)

Use **named grid areas** for semantic, maintainable layouts.

```css
.results-layout {
  display: grid;
  grid-template-columns: repeat(12, 1fr);
  grid-template-rows: auto auto auto;
  gap: 1rem;
  grid-template-areas:
    "viewer viewer viewer viewer viewer spectra spectra spectra spectra spectra spectra spectra"
    "viewer viewer viewer viewer viewer h1-table h1-table h1-table c13-table c13-table c13-table c13-table"
    "meta meta meta meta meta h1-table h1-table h1-table c13-table c13-table c13-table c13-table";
}

.viewer-card { grid-area: viewer; }
.spectra-card { grid-area: spectra; }
.h1-table-card { grid-area: h1-table; }
.c13-table-card { grid-area: c13-table; }
.metadata-card { grid-area: meta; }
```

**Why this approach:**
- **Semantic**: Named areas are self-documenting
- **Flexible**: Easy to rearrange layout by changing template string
- **Responsive**: Different `grid-template-areas` for mobile/tablet/desktop
- **Maintainable**: No span calculations, just descriptive names

**Source:** [MDN Grid Template Areas](https://developer.mozilla.org/en-US/docs/Web/CSS/Guides/Grid_layout/Grid_template_areas), [iamsteve Bento Layout Tutorial](https://iamsteve.me/blog/bento-layout-css-grid)

### Grid-Auto-Flow: Dense

Use `grid-auto-flow: dense` to fill gaps with smaller items.

```css
.bento-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
  grid-auto-rows: minmax(150px, auto);
  gap: 1rem;
  grid-auto-flow: dense; /* Magic property for tight packing */
}

.bento-item-large {
  grid-column: span 2;
  grid-row: span 2;
}

.bento-item-wide {
  grid-column: span 2;
}

.bento-item-tall {
  grid-row: span 2;
}
```

**Why `dense`:**
- Fills holes in the grid with smaller items
- Creates the "tetris" bento aesthetic
- Works automatically without manual positioning

**Note:** Dense packing can disrupt document order for screen readers. Use `tabindex` or reorder DOM if keyboard navigation is critical.

**Source:** [Speckyboy CSS Bento Grids](https://speckyboy.com/css-bento-grid-layouts/), [Codemotion Bento Layout](https://www.codemotion.com/magazine/frontend/lets-create-a-bento-box-design-layout-using-modern-css/)

### Responsive Breakpoints

Redefine grid structure at breakpoints, not individual items.

```css
/* Mobile: stacked */
@media (max-width: 768px) {
  .results-layout {
    grid-template-columns: 1fr;
    grid-template-areas:
      "viewer"
      "spectra"
      "h1-table"
      "c13-table"
      "meta";
  }
}

/* Tablet: 2-column */
@media (min-width: 769px) and (max-width: 1024px) {
  .results-layout {
    grid-template-columns: repeat(6, 1fr);
    grid-template-areas:
      "viewer viewer viewer spectra spectra spectra"
      "h1-table h1-table h1-table c13-table c13-table c13-table"
      "meta meta meta meta meta meta";
  }
}
```

**Why this works:**
- Single source of truth per breakpoint (the `grid-template-areas` string)
- Items don't need individual media queries
- Easier to reason about layout changes

**Source:** [FreeCodeCamp Bento Grids](https://www.freecodecamp.org/news/bento-grids-in-web-design/)

### Container Queries (Optional Enhancement)

Container queries are production-ready (baseline support since 2023). Use for self-contained responsive cards.

```css
.card {
  container-type: inline-size;
  container-name: card;
}

@container card (min-width: 400px) {
  .card-content {
    display: grid;
    grid-template-columns: 1fr 1fr;
  }
}
```

**When to use:**
- Cards that appear in multiple contexts (sidebar vs main area)
- Component-level responsiveness independent of viewport

**When NOT to use:**
- Page-level layouts (use media queries)
- If container size depends on content (creates circular dependency)

**Source:** [Container Queries in 2026](https://blog.logrocket.com/container-queries-2026/), [MDN Container Queries](https://developer.mozilla.org/en-US/docs/Web/CSS/Guides/Containment/Container_queries)

---

## Glassmorphism Implementation

### Core Properties

Glassmorphism requires **four properties working together:**

```css
.glass-card {
  /* 1. Semi-transparent background (REQUIRED) */
  background-color: rgba(255, 255, 255, 0.15);

  /* 2. Backdrop blur (the "frosted glass" effect) */
  backdrop-filter: blur(10px) brightness(1.05);

  /* 3. Subtle border for definition */
  border: 1px solid rgba(255, 255, 255, 0.3);

  /* 4. Soft shadow for depth */
  box-shadow: 0 8px 32px rgba(0, 0, 0, 0.1);

  /* Optional: rounded corners for softer aesthetic */
  border-radius: 12px;
}
```

**Source:** [Job Huntley Glassmorphism Guide](https://www.jobhuntley.com/blog/web-design-trends-for-2026-the-rise-of-glassmorphism-and-how-to-achieve-it-with-css), [MDN backdrop-filter](https://developer.mozilla.org/en-US/docs/Web/CSS/backdrop-filter)

### Light Mode Recommendations

For light backgrounds, use **lower alpha values** to maintain contrast.

```css
/* Light mode: Higher transparency (lower alpha) */
.glass-card-light {
  background-color: rgba(255, 255, 255, 0.1);  /* 10% opacity */
  backdrop-filter: blur(8px) saturate(120%);
  border: 1px solid rgba(255, 255, 255, 0.25);
  box-shadow: 0 4px 16px rgba(0, 0, 0, 0.08);
}
```

**Alpha value guidelines for light mode:**
- `0.1 - 0.15`: Over colorful backgrounds or images
- `0.15 - 0.25`: Over gradient backgrounds
- `0.25 - 0.4`: Over solid light backgrounds (white/gray)

**Why lower alpha:** Light mode has less inherent contrast. Too much opacity washes out the blur effect and makes text hard to read.

**Source:** [Inverness Glassmorphism 2026](https://invernessdesignstudio.com/glassmorphism-what-it-is-and-how-to-use-it-in-2026), [UX Pilot Glassmorphism Best Practices](https://uxpilot.ai/blogs/glassmorphism-ui)

### Backdrop-Filter Values

**Blur radius:**
- `5px - 10px`: Subtle glass effect, better performance
- `10px - 20px`: Standard glassmorphism (recommended)
- `20px+`: Heavy blur, performance cost, use sparingly

**Combining filters:**
```css
/* Enhance depth with multiple filters */
backdrop-filter: blur(10px) saturate(150%) brightness(1.1);
```

- `saturate()`: Boost colors behind glass (110%-150%)
- `brightness()`: Lighten/darken backdrop (0.9-1.1)
- `contrast()`: Increase definition (100%-120%)

**Performance note:** Each filter function adds GPU cost. Keep to 1-2 filters per element.

**Source:** [MDN backdrop-filter](https://developer.mozilla.org/en-US/docs/Web/CSS/backdrop-filter)

### Browser Support & Fallback

**Browser support: 92% (as of 2026)**

| Browser | Version | Support |
|---------|---------|---------|
| Chrome | 76+ | Full |
| Firefox | 103+ | Full |
| Safari | 9+ | Full |
| Edge | 79+ | Full |
| iOS Safari | 9+ | Full |

**Fallback strategy:**

```css
.glass-card {
  /* Fallback: solid background */
  background-color: rgba(255, 255, 255, 0.9);

  /* Progressive enhancement */
  backdrop-filter: blur(10px);
}

/* Higher opacity fallback if backdrop-filter unsupported */
@supports not (backdrop-filter: blur(1px)) {
  .glass-card {
    background-color: rgba(255, 255, 255, 0.95);
  }
}
```

**Why this works:**
- Modern browsers get glassmorphism
- Older browsers get solid semi-transparent background
- Functionality preserved in both cases

**Source:** [Can I Use backdrop-filter](https://caniuse.com/css-backdrop-filter), [MDN backdrop-filter](https://developer.mozilla.org/en-US/docs/Web/CSS/backdrop-filter)

### Accessibility Considerations

**CRITICAL: Glassmorphism often fails WCAG contrast requirements.**

**WCAG 2.2 minimums:**
- Normal text: 4.5:1 contrast ratio
- Large text (18px+ or 14px+ bold): 3:1 contrast ratio
- UI components: 3:1 contrast ratio

**Solutions:**

1. **Text shadows for legibility:**
```css
.glass-card-text {
  color: #1a1a1a;
  text-shadow: 0 1px 2px rgba(255, 255, 255, 0.8);
}
```

2. **Darker inner background for text areas:**
```css
.glass-card-content {
  background: rgba(255, 255, 255, 0.6); /* Higher opacity under text */
  padding: 1rem;
  border-radius: 8px;
}
```

3. **Use glassmorphism for container cards, solid backgrounds for text content:**
```css
.glass-outer {
  background: rgba(255, 255, 255, 0.15);
  backdrop-filter: blur(10px);
  padding: 1rem;
}

.glass-outer .content-area {
  background: rgba(255, 255, 255, 0.95); /* Nearly solid for text */
  padding: 1rem;
}
```

4. **Respect user preferences:**
```css
@media (prefers-reduced-transparency: reduce) {
  .glass-card {
    backdrop-filter: none;
    background-color: rgba(255, 255, 255, 0.95);
  }
}
```

**Source:** [Axess Lab Glassmorphism Accessibility](https://axesslab.com/glassmorphism-meets-accessibility-can-frosted-glass-be-inclusive/), [New Target Glassmorphism Accessibility](https://www.newtarget.com/web-insights-blog/glassmorphism/), [WebAIM Contrast](https://webaim.org/articles/contrast/)

---

## Animations and Transitions

### Performance-First Approach

**Animate ONLY transform and opacity** (GPU-accelerated, no reflow).

```css
.glass-card {
  transition: transform 0.2s ease-out, opacity 0.2s ease-out;
}

.glass-card:hover {
  transform: translateY(-4px) scale(1.01);
  opacity: 1;
}
```

**Why:**
- `transform` and `opacity` are GPU-accelerated
- No layout recalculation (reflow)
- Smooth 60fps performance

**NEVER animate:** `width`, `height`, `top`, `left`, `margin`, `padding` (causes reflow, janky).

**Source:** [Josh Comeau CSS Transitions](https://www.joshwcomeau.com/animation/css-transitions/), [design.dev CSS Transitions Guide](https://design.dev/guides/css-transitions/)

### Hover Effects

**Card lift on hover:**
```css
.glass-card {
  transform: translateZ(0); /* Create stacking context for GPU */
  transition: transform 0.2s ease-out, box-shadow 0.2s ease-out;
}

.glass-card:hover {
  transform: translateY(-4px);
  box-shadow: 0 12px 40px rgba(0, 0, 0, 0.15);
}
```

**Backdrop blur intensity change (use sparingly, GPU-intensive):**
```css
.glass-card {
  backdrop-filter: blur(8px);
  transition: backdrop-filter 0.3s ease;
}

.glass-card:hover {
  backdrop-filter: blur(12px);
}
```

**WARNING:** Animating `backdrop-filter` is GPU-intensive. Use only on key interactive elements, not all cards.

**Source:** [MDN backdrop-filter](https://developer.mozilla.org/en-US/docs/Web/CSS/backdrop-filter), [CSS-Tricks Transitions](https://css-tricks.com/almanac/properties/t/transition/)

### Timing Functions

**Best practices for UI transitions:**
- **0.2-0.3s duration**: Feels responsive without being abrupt
- **ease-out**: Things entering or growing (most UI interactions)
- **ease-in**: Things leaving or shrinking
- **ease-in-out**: Symmetrical movements

```css
/* Quick, snappy hover enter */
.card {
  transition: transform 0.2s ease-out;
}

/* Slightly slower, relaxed hover exit */
.card:not(:hover) {
  transition: transform 0.3s ease-in;
}
```

**Source:** [design.dev CSS Transitions](https://design.dev/guides/css-transitions/)

### Mobile Considerations

**Disable hover effects on touch devices:**

```css
@media (hover: hover) {
  .glass-card:hover {
    transform: translateY(-4px);
  }
}
```

**Why:** Touch devices don't have hover states. This prevents awkward tap-to-hover behavior.

**Source:** [WebPeak Animation Trends 2026](https://webpeak.org/blog/css-js-animation-trends/)

---

## CSS Custom Properties (Variables)

Use custom properties for **theming and DRY principles**.

```css
:root {
  /* Glassmorphism theme */
  --glass-bg-light: rgba(255, 255, 255, 0.15);
  --glass-border-light: rgba(255, 255, 255, 0.3);
  --glass-shadow: 0 8px 32px rgba(0, 0, 0, 0.1);
  --glass-blur: blur(10px);
  --glass-blur-heavy: blur(16px);

  /* Spacing (match bento grid gaps) */
  --bento-gap: 1rem;
  --bento-gap-lg: 1.5rem;

  /* Border radius */
  --radius-sm: 8px;
  --radius-md: 12px;
  --radius-lg: 16px;

  /* Transitions */
  --transition-fast: 0.2s ease-out;
  --transition-medium: 0.3s ease-out;
}

.glass-card {
  background-color: var(--glass-bg-light);
  backdrop-filter: var(--glass-blur);
  border: 1px solid var(--glass-border-light);
  box-shadow: var(--glass-shadow);
  border-radius: var(--radius-md);
  transition: transform var(--transition-fast);
}
```

**Why custom properties:**
- Consistency across components
- Easy theme adjustments (change one value, update all cards)
- No preprocessor needed
- Scoped updates (can override at component level)

---

## What NOT to Use

### 1. CSS Frameworks (Tailwind, Bootstrap, etc.)

**Why avoid:**
- Requires build step (Node.js, PostCSS)
- Adds bundle size and complexity
- Bento layouts need custom CSS anyway
- Conflicts with existing "no build" philosophy (3Dmol.js from CDN)

**Exception:** If project already has a build pipeline, Tailwind's utility classes can speed up glassmorphism (`backdrop-blur-lg`, `bg-white/10`). But given current stack, **not recommended**.

### 2. CSS Preprocessors (Sass, Less)

**Why avoid:**
- Modern CSS custom properties replace most use cases
- Adds build step
- Nesting is coming to native CSS (some browser support already)
- Grid and backdrop-filter work identically in preprocessors

**Use CSS custom properties instead.**

### 3. CSS-in-JS (Styled-components, Emotion)

**Why avoid:**
- Requires JavaScript framework (React, Vue)
- Adds runtime overhead
- Jinja2 templates are server-rendered
- Overkill for styling-only changes

**Project uses Jinja2 + vanilla JS, not a JS framework.**

### 4. JavaScript Grid Libraries (Masonry.js, Packery)

**Why avoid:**
- CSS Grid with `grid-auto-flow: dense` achieves bento packing natively
- No JS dependency
- Better performance (layout in CSS, not JS)
- More maintainable

**Use CSS Grid's native features instead.**

**Source:** [DEV Community Responsive Bento Grid](https://dev.to/velox-web/how-to-build-a-responsive-bento-grid-with-tailwind-css-no-masonryjs-3f2c)

### 5. Heavy Animation Libraries (GSAP, Anime.js)

**Why avoid:**
- CSS transitions handle hover/state changes perfectly
- Adds bundle size
- Overkill for simple card animations

**Use CSS transitions. Reserve JS animation for complex timing (none needed here).**

---

## Browser Compatibility Summary

| Feature | Support | Fallback Needed? |
|---------|---------|------------------|
| CSS Grid | 96%+ (since 2017) | No |
| `grid-template-areas` | 96%+ (since 2017) | No |
| `grid-auto-flow: dense` | 96%+ (since 2017) | No |
| `backdrop-filter` | 92% (since 2023) | Yes (solid bg) |
| CSS Custom Properties | 97%+ (since 2016) | No |
| Container Queries | 90%+ (baseline 2023) | Optional enhancement |
| `@supports` | 97%+ (since 2013) | No |

**Target browsers (2026):**
- Chrome/Edge 100+
- Firefox 100+
- Safari 15+
- iOS Safari 15+

All key features are supported. Only `backdrop-filter` needs a fallback (solid background).

---

## Implementation Strategy

### Phase 1: Remove Pico CSS, Add Base Styles

1. Remove Pico CSS `<link>` from `base.html`
2. Create `static/css/main.css` with:
   - CSS reset/normalize
   - Custom properties (colors, spacing, radius, transitions)
   - Base typography
   - Utility classes

### Phase 2: Glassmorphism Card Component

Create reusable `.glass-card` class:
```css
.glass-card {
  background-color: rgba(255, 255, 255, 0.15);
  backdrop-filter: blur(10px);
  border: 1px solid rgba(255, 255, 255, 0.3);
  box-shadow: 0 8px 32px rgba(0, 0, 0, 0.1);
  border-radius: 12px;
  padding: 1.5rem;
  transition: transform 0.2s ease-out, box-shadow 0.2s ease-out;
}

@supports not (backdrop-filter: blur(1px)) {
  .glass-card {
    background-color: rgba(255, 255, 255, 0.95);
  }
}

@media (hover: hover) {
  .glass-card:hover {
    transform: translateY(-4px);
    box-shadow: 0 12px 40px rgba(0, 0, 0, 0.15);
  }
}
```

### Phase 3: Bento Grid Layouts

Create page-specific layouts:
- `results.html`: 12-column grid with named areas (viewer, spectra, tables, metadata)
- `submit.html`: Simpler 2-column form layout
- `status.html`: Status cards in auto-fit grid

### Phase 4: Responsive Breakpoints

Add media queries for mobile/tablet, redefining `grid-template-areas` at each breakpoint.

### Phase 5: Accessibility Pass

1. Test contrast ratios with WebAIM checker
2. Add text shadows or darker backgrounds where needed
3. Test keyboard navigation (tab order with `grid-auto-flow: dense`)
4. Add `prefers-reduced-transparency` support

---

## File Structure

```
static/
├── css/
│   ├── main.css          # Base styles, variables, resets
│   ├── components.css    # Glass cards, buttons, forms
│   ├── layouts.css       # Bento grid layouts per page
│   └── responsive.css    # Media queries
└── js/
    └── [existing files]
```

**Alternative (single file):**
```
static/
└── css/
    └── main.css  # All styles in one file (simpler for small project)
```

Single file is fine for this project size. Split if CSS exceeds ~800 lines.

---

## Sources

### Bento Grid Layouts
- [CSS Bento Grid Layouts](https://speckyboy.com/css-bento-grid-layouts/)
- [Build a Bento Layout with CSS Grid](https://iamsteve.me/blog/bento-layout-css-grid)
- [Bento Box Design with Modern CSS](https://www.codemotion.com/magazine/frontend/lets-create-a-bento-box-design-layout-using-modern-css/)
- [Master CSS Grid for Bento UI](https://medium.com/@lilskyjuicebytes/design-to-code-1-brewbolt-bento-ui-with-html-css-a128f64ebceb)
- [FreeCodeCamp Bento Grids](https://www.freecodecamp.org/news/bento-grids-in-web-design/)
- [MDN Grid Template Areas](https://developer.mozilla.org/en-US/docs/Web/CSS/Guides/Grid_layout/Grid_template_areas)

### Glassmorphism
- [Glassmorphism in 2026](https://invernessdesignstudio.com/glassmorphism-what-it-is-and-how-to-use-it-in-2026)
- [Web Design Trends 2026: Glassmorphism](https://www.jobhuntley.com/blog/web-design-trends-for-2026-the-rise-of-glassmorphism-and-how-to-achieve-it-with-css)
- [UX Pilot Glassmorphism UI](https://uxpilot.ai/blogs/glassmorphism-ui)
- [MDN backdrop-filter](https://developer.mozilla.org/en-US/docs/Web/CSS/backdrop-filter)
- [Glass UI Generator](https://ui.glass/generator/)
- [Create Glassmorphic UI with CSS](https://blog.openreplay.com/create-glassmorphic-ui-css/)

### Browser Support
- [Can I Use backdrop-filter](https://caniuse.com/css-backdrop-filter)
- [LambdaTest CSS Backdrop Filter](https://www.lambdatest.com/web-technologies/css-backdrop-filter)

### Animations & Performance
- [Josh Comeau CSS Transitions](https://www.joshwcomeau.com/animation/css-transitions/)
- [design.dev CSS Transitions Guide](https://design.dev/guides/css-transitions/)
- [CSS-Tricks Transitions](https://css-tricks.com/almanac/properties/t/transition/)
- [WebPeak CSS/JS Animation Trends](https://webpeak.org/blog/css-js-animation-trends/)

### Accessibility
- [Axess Lab Glassmorphism Accessibility](https://axesslab.com/glassmorphism-meets-accessibility-can-frosted-glass-be-inclusive/)
- [New Target Glassmorphism Accessibility](https://www.newtarget.com/web-insights-blog/glassmorphism/)
- [Nielsen Norman Glassmorphism](https://www.nngroup.com/articles/glassmorphism/)
- [WebAIM Contrast](https://webaim.org/articles/contrast/)

### Container Queries
- [Container Queries in 2026](https://blog.logrocket.com/container-queries-2026/)
- [MDN Container Queries](https://developer.mozilla.org/en-US/docs/Web/CSS/Guides/Containment/Container_queries)
- [Josh Comeau Container Queries](https://www.joshwcomeau.com/css/container-queries-introduction/)
