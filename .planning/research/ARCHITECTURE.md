# Architecture Research: CSS Organization for Bento Grid + Glassmorphism

**Project:** qm-nmr-calc UI Redesign
**Researched:** 2026-01-29
**Confidence:** HIGH

## Executive Summary

Modern CSS in 2026 provides native features that eliminate the need for build tools while maintaining excellent organization and maintainability. For this Jinja2 template-based app, a multi-file CSS architecture using CSS Cascade Layers, custom properties, and BEM naming provides the ideal balance of maintainability and simplicity.

**Key findings:**
- CSS Cascade Layers provide architectural structure without preprocessors
- CSS custom properties handle theming and design tokens natively
- Multiple CSS files can be loaded without build steps
- BEM naming prevents specificity conflicts
- Native CSS nesting is available in all modern browsers

## Recommended File Structure

### Directory Organization

```
src/qm_nmr_calc/api/static/css/
├── layers.css              # Layer order declaration (loaded first)
├── reset.css               # Minimal reset/normalization
├── tokens.css              # Design tokens (custom properties)
├── base.css                # Base element styles
├── layout.css              # Layout systems (bento grid)
├── components/
│   ├── glassmorphic-card.css
│   ├── form.css
│   ├── status-indicator.css
│   ├── molecule-viewer.css
│   └── spectrum-display.css
└── utilities.css           # Utility classes
```

### Loading Order in base.html

```html
<head>
    <!-- 1. Layer declarations (must be first) -->
    <link rel="stylesheet" href="{{ url_for('static', path='/css/layers.css') }}">

    <!-- 2. Reset layer -->
    <link rel="stylesheet" href="{{ url_for('static', path='/css/reset.css') }}">

    <!-- 3. Design tokens -->
    <link rel="stylesheet" href="{{ url_for('static', path='/css/tokens.css') }}">

    <!-- 4. Base and layout -->
    <link rel="stylesheet" href="{{ url_for('static', path='/css/base.css') }}">
    <link rel="stylesheet" href="{{ url_for('static', path='/css/layout.css') }}">

    <!-- 5. Components -->
    <link rel="stylesheet" href="{{ url_for('static', path='/css/components/glassmorphic-card.css') }}">
    <link rel="stylesheet" href="{{ url_for('static', path='/css/components/form.css') }}">
    <link rel="stylesheet" href="{{ url_for('static', path='/css/components/status-indicator.css') }}">
    <link rel="stylesheet" href="{{ url_for('static', path='/css/components/molecule-viewer.css') }}">
    <link rel="stylesheet" href="{{ url_for('static', path='/css/components/spectrum-display.css') }}">

    <!-- 6. Utilities (highest priority) -->
    <link rel="stylesheet" href="{{ url_for('static', path='/css/utilities.css') }}">
</head>
```

**Rationale:** Multiple small files improve maintainability and debugging. Cascade order is explicit through loading sequence. No concatenation or minification needed for development or production (HTTP/2 handles multiple files efficiently).

## CSS Cascade Layers Implementation

### layers.css (Load First)

```css
/**
 * CSS Cascade Layers Declaration
 *
 * Defines layer priority order. Must be loaded before any other stylesheets.
 * Priority: reset (lowest) → base → layout → components → utilities (highest)
 */

@layer reset, base, layout, components, utilities;
```

**Why:** Explicit layer declaration at the top provides single-location control over cascade priority. Un-layered styles (if any) automatically get highest priority, enabling incremental migration.

### Wrapping Styles in Layers

Each file wraps its contents in the appropriate layer:

```css
/* reset.css */
@layer reset {
    *, *::before, *::after {
        box-sizing: border-box;
    }

    body {
        margin: 0;
        line-height: 1.5;
    }
}
```

```css
/* base.css */
@layer base {
    body {
        font-family: var(--font-sans);
        color: var(--color-text);
        background: var(--color-bg);
    }

    h1, h2, h3 {
        line-height: 1.2;
    }
}
```

```css
/* components/glassmorphic-card.css */
@layer components {
    .glass-card {
        background: var(--glass-bg);
        backdrop-filter: var(--glass-blur);
        border: var(--glass-border);
        border-radius: var(--glass-radius);
    }
}
```

**Benefit:** Layer priority is automatic. Component styles never accidentally override utilities. No !important needed.

## CSS Custom Properties (Design Tokens)

### tokens.css Structure

```css
/**
 * Design Tokens
 *
 * All visual design decisions as CSS custom properties.
 * Organized by category for maintainability.
 */

:root {
    /* === Color System === */

    /* Base colors */
    --color-white: hsl(0 0% 100%);
    --color-gray-50: hsl(220 13% 98%);
    --color-gray-100: hsl(220 13% 95%);
    --color-gray-200: hsl(220 13% 91%);
    --color-gray-300: hsl(220 13% 83%);
    --color-gray-400: hsl(220 13% 69%);
    --color-gray-500: hsl(220 13% 50%);
    --color-gray-600: hsl(220 13% 38%);
    --color-gray-700: hsl(220 13% 25%);
    --color-gray-800: hsl(220 13% 15%);
    --color-gray-900: hsl(220 13% 9%);

    /* Brand colors */
    --color-primary: hsl(217 91% 60%);      /* Blue */
    --color-primary-hover: hsl(217 91% 50%);
    --color-success: hsl(142 76% 36%);      /* Green */
    --color-warning: hsl(45 93% 47%);       /* Yellow */
    --color-error: hsl(0 84% 60%);          /* Red */

    /* Semantic colors */
    --color-bg: var(--color-gray-50);
    --color-text: var(--color-gray-900);
    --color-text-muted: var(--color-gray-600);
    --color-border: var(--color-gray-300);

    /* === Glassmorphism === */

    /* Background blur values */
    --glass-blur-light: blur(8px);          /* Subtle frost */
    --glass-blur-medium: blur(12px);        /* Standard glass */
    --glass-blur-heavy: blur(16px);         /* Deep blur */

    /* Glass backgrounds (semi-transparent) */
    --glass-bg-light: hsl(0 0% 100% / 0.7);
    --glass-bg-medium: hsl(0 0% 100% / 0.5);
    --glass-bg-subtle: hsl(0 0% 100% / 0.9);

    /* Glass borders */
    --glass-border-light: 1px solid hsl(0 0% 100% / 0.3);
    --glass-border-medium: 1px solid hsl(0 0% 100% / 0.5);

    /* Combined glass effect tokens */
    --glass-blur: var(--glass-blur-medium);
    --glass-bg: var(--glass-bg-light);
    --glass-border: var(--glass-border-light);
    --glass-shadow: 0 8px 32px hsl(0 0% 0% / 0.1);
    --glass-radius: var(--radius-lg);

    /* === Bento Grid === */

    --bento-cols: 6;                        /* Desktop columns */
    --bento-rows: 4;                        /* Desktop rows */
    --bento-gap: 1rem;                      /* Gap between cards */
    --bento-min-height: 120px;              /* Minimum card height */

    /* === Spacing === */

    --space-xs: 0.25rem;    /* 4px */
    --space-sm: 0.5rem;     /* 8px */
    --space-md: 1rem;       /* 16px */
    --space-lg: 1.5rem;     /* 24px */
    --space-xl: 2rem;       /* 32px */
    --space-2xl: 3rem;      /* 48px */
    --space-3xl: 4rem;      /* 64px */

    /* === Border Radius === */

    --radius-sm: 0.25rem;   /* 4px */
    --radius-md: 0.5rem;    /* 8px */
    --radius-lg: 1rem;      /* 16px */
    --radius-xl: 1.5rem;    /* 24px */
    --radius-full: 9999px;  /* Pill shape */

    /* === Typography === */

    --font-sans: system-ui, -apple-system, 'Segoe UI', Roboto, sans-serif;
    --font-mono: 'Cascadia Code', 'Fira Code', 'Monaco', monospace;

    --text-xs: 0.75rem;     /* 12px */
    --text-sm: 0.875rem;    /* 14px */
    --text-base: 1rem;      /* 16px */
    --text-lg: 1.125rem;    /* 18px */
    --text-xl: 1.25rem;     /* 20px */
    --text-2xl: 1.5rem;     /* 24px */
    --text-3xl: 1.875rem;   /* 30px */
    --text-4xl: 2.25rem;    /* 36px */

    /* === Transitions === */

    --transition-fast: 150ms ease;
    --transition-base: 250ms ease;
    --transition-slow: 350ms ease;

    /* === Z-index === */

    --z-base: 0;
    --z-dropdown: 100;
    --z-sticky: 200;
    --z-modal: 300;
    --z-toast: 400;
}

/* Responsive token adjustments */
@media (max-width: 768px) {
    :root {
        --bento-cols: 2;
        --bento-rows: auto;
        --bento-gap: 0.75rem;
        --glass-blur: blur(6px);  /* Reduced for mobile performance */
    }
}

@media (max-width: 480px) {
    :root {
        --bento-cols: 1;
        --space-lg: 1rem;
        --space-xl: 1.5rem;
    }
}
```

**Why custom properties:**
- Single source of truth for design decisions
- Runtime theming without preprocessors
- Easy to override for responsive breakpoints
- Semantic naming (--glass-blur instead of hardcoded blur(12px))
- Composable (--glass-bg can reference --color-white)

## Component Naming Convention: BEM

### BEM Structure

**Block__Element--Modifier**

- Block: standalone component (.glass-card)
- Element: part of block (.glass-card__header)
- Modifier: variant or state (.glass-card--featured)

### Example: Glassmorphic Card

```css
/* components/glassmorphic-card.css */
@layer components {
    /* Block */
    .glass-card {
        background: var(--glass-bg);
        backdrop-filter: var(--glass-blur);
        -webkit-backdrop-filter: var(--glass-blur);  /* Safari support */
        border: var(--glass-border);
        border-radius: var(--glass-radius);
        box-shadow: var(--glass-shadow);
        padding: var(--space-lg);
        transition: transform var(--transition-base);
    }

    /* Elements */
    .glass-card__header {
        margin-bottom: var(--space-md);
        padding-bottom: var(--space-md);
        border-bottom: 1px solid var(--color-border);
    }

    .glass-card__title {
        font-size: var(--text-2xl);
        font-weight: 600;
        margin: 0;
    }

    .glass-card__body {
        color: var(--color-text);
    }

    .glass-card__footer {
        margin-top: var(--space-lg);
        padding-top: var(--space-md);
        border-top: 1px solid var(--color-border);
    }

    /* Modifiers */
    .glass-card--featured {
        background: var(--glass-bg-medium);
        border-width: 2px;
        box-shadow: 0 12px 48px hsl(0 0% 0% / 0.15);
    }

    .glass-card--subtle {
        background: var(--glass-bg-subtle);
        backdrop-filter: var(--glass-blur-light);
        -webkit-backdrop-filter: var(--glass-blur-light);
    }

    .glass-card--interactive:hover {
        transform: translateY(-4px);
        box-shadow: 0 12px 48px hsl(0 0% 0% / 0.15);
    }
}
```

### HTML Usage in Jinja2 Templates

```html
{% extends "base.html" %}

{% block content %}
<div class="glass-card glass-card--featured">
    <header class="glass-card__header">
        <h2 class="glass-card__title">{{ molecule.name }}</h2>
    </header>
    <div class="glass-card__body">
        <p>Chemical formula: {{ molecule.formula }}</p>
    </div>
    <footer class="glass-card__footer">
        <a href="{{ url }}" role="button">View Results</a>
    </footer>
</div>
{% endblock %}
```

**Benefits:**
- No naming conflicts across components
- Clear component boundaries in HTML
- Easy to search (grep for "glass-card")
- Self-documenting (name indicates hierarchy)

## Bento Grid Layout System

### layout.css

```css
/**
 * Bento Grid Layout
 *
 * Asymmetric grid system for card-based layouts.
 * Uses CSS Grid with custom properties for responsive control.
 */

@layer layout {
    .bento-grid {
        display: grid;
        grid-template-columns: repeat(var(--bento-cols), 1fr);
        gap: var(--bento-gap);
        width: 100%;
        padding: var(--space-lg);
    }

    /* Grid cell spanning utilities */
    .bento-grid__item {
        min-height: var(--bento-min-height);
        display: flex;
        flex-direction: column;
    }

    /* Span patterns for asymmetry */
    .bento-grid__item--span-2 {
        grid-column: span 2;
    }

    .bento-grid__item--span-3 {
        grid-column: span 3;
    }

    .bento-grid__item--span-row-2 {
        grid-row: span 2;
    }

    .bento-grid__item--span-row-3 {
        grid-row: span 3;
    }

    /* Feature grid item (2x2) */
    .bento-grid__item--feature {
        grid-column: span 2;
        grid-row: span 2;
    }

    /* Full-width item */
    .bento-grid__item--full {
        grid-column: 1 / -1;
    }

    /* Auto-flow dense for optimal packing */
    .bento-grid--dense {
        grid-auto-flow: dense;
    }

    /* Responsive adjustments handled by token media queries */
}
```

### Usage Example

```html
<div class="bento-grid bento-grid--dense">
    <!-- Featured molecule viewer (2x2) -->
    <div class="bento-grid__item bento-grid__item--feature">
        <div class="glass-card">
            <div id="molecule-viewer"></div>
        </div>
    </div>

    <!-- 1H spectrum (2x1) -->
    <div class="bento-grid__item bento-grid__item--span-2">
        <div class="glass-card">
            <img src="{{ h1_spectrum_url }}" alt="1H NMR Spectrum">
        </div>
    </div>

    <!-- 13C spectrum (2x1) -->
    <div class="bento-grid__item bento-grid__item--span-2">
        <div class="glass-card">
            <img src="{{ c13_spectrum_url }}" alt="13C NMR Spectrum">
        </div>
    </div>

    <!-- Metadata card (1x1) -->
    <div class="bento-grid__item">
        <div class="glass-card glass-card--subtle">
            <dl class="metadata-list">
                <dt>Solvent</dt>
                <dd>{{ job.solvent }}</dd>
            </dl>
        </div>
    </div>

    <!-- Downloads card (1x1) -->
    <div class="bento-grid__item">
        <div class="glass-card glass-card--subtle">
            <h3>Downloads</h3>
            <a href="{{ sdf_url }}" download>Download SDF</a>
        </div>
    </div>
</div>
```

**Key features:**
- Responsive via custom property breakpoints
- Dense packing with `grid-auto-flow: dense`
- Flexible spanning without complex calculations
- Works with glassmorphic cards naturally

## Template Integration Patterns

### Page-Specific CSS (Optional)

For page-specific styles, use inline `<style>` blocks in template heads:

```html
{% extends "base.html" %}

{% block head %}
<style>
    /* Results page specific layout adjustments */
    .results-layout {
        max-width: 1400px;
        margin: 0 auto;
    }
</style>
{% endblock %}
```

**When to use:**
- Page has unique layout needs
- Style won't be reused elsewhere
- Keeps component CSS files focused

### Dynamic Classes from Backend

Pass CSS modifier classes from FastAPI views:

```python
# api/routes.py
@router.get("/status/{job_id}")
async def get_status(job_id: str):
    job = get_job(job_id)

    # Map job status to CSS modifier
    status_modifier_map = {
        "pending": "status-indicator--pending",
        "running": "status-indicator--running",
        "complete": "status-indicator--complete",
        "failed": "status-indicator--error",
    }

    return templates.TemplateResponse("status.html", {
        "request": request,
        "job": job,
        "status_class": status_modifier_map[job.status],
    })
```

```html
<!-- status.html -->
<div class="status-indicator {{ status_class }}">
    <span class="status-indicator__dot"></span>
    <span class="status-indicator__label">{{ job.status }}</span>
</div>
```

**Benefits:**
- Backend controls styling logic
- No JS needed for status classes
- Type-safe with Python enums

## Implementation Order

### Phase 1: Foundation (Days 1-2)

1. Create directory structure
2. Write layers.css (layer declarations)
3. Write tokens.css (all design tokens)
4. Write reset.css (minimal reset)
5. Write base.css (element defaults)
6. Update base.html with new link tags

**Deliverable:** Design tokens established, base styles in place

### Phase 2: Layout System (Day 3)

1. Write layout.css (bento grid system)
2. Test grid on simple HTML mockup
3. Adjust responsive breakpoints in tokens.css
4. Document grid span patterns

**Deliverable:** Working bento grid with responsive behavior

### Phase 3: Glassmorphic Components (Days 4-5)

1. Write components/glassmorphic-card.css
2. Test backdrop-filter browser support
3. Add -webkit- prefixes for Safari
4. Verify blur performance on mobile

**Deliverable:** Reusable glass card component

### Phase 4: Page-Specific Components (Days 6-8)

Create component CSS files in parallel:
- components/form.css (submit page forms)
- components/status-indicator.css (status page)
- components/molecule-viewer.css (3D viewer wrapper)
- components/spectrum-display.css (results page spectra)

**Deliverable:** All page components styled

### Phase 5: Template Integration (Days 9-10)

1. Migrate submit.html to bento + glass
2. Migrate status.html to bento + glass
3. Migrate results.html to bento + glass
4. Remove Pico CSS CDN link
5. Write utilities.css (final overrides)

**Deliverable:** All templates using new design system

### Phase 6: Polish (Days 11-12)

1. Add hover states and transitions
2. Verify responsive behavior on mobile
3. Test cross-browser (Chrome, Firefox, Safari)
4. Optimize glassmorphism for performance
5. Document component usage in comments

**Deliverable:** Production-ready UI

## Anti-Patterns to Avoid

### 1. Inline Styles in HTML

**Don't:**
```html
<div style="backdrop-filter: blur(12px); background: rgba(255,255,255,0.7);">
```

**Do:**
```html
<div class="glass-card">
```

**Why:** Inline styles defeat maintainability. Design changes require HTML edits across all templates.

### 2. Overly-Specific Selectors

**Don't:**
```css
body main.container article.glass-card div.content p {
    color: var(--color-text);
}
```

**Do:**
```css
.glass-card__content {
    color: var(--color-text);
}
```

**Why:** High specificity makes overrides difficult and defeats layer architecture.

### 3. Hardcoded Values

**Don't:**
```css
.card {
    backdrop-filter: blur(12px);
    background: rgba(255, 255, 255, 0.7);
    border-radius: 16px;
}
```

**Do:**
```css
.glass-card {
    backdrop-filter: var(--glass-blur);
    background: var(--glass-bg);
    border-radius: var(--glass-radius);
}
```

**Why:** Design tokens enable consistent changes. Hardcoded values create maintenance burden.

### 4. !important for Cascade Control

**Don't:**
```css
.my-component {
    color: red !important;
}
```

**Do:**
```css
/* Use layers for priority */
@layer utilities {
    .text-error {
        color: var(--color-error);
    }
}
```

**Why:** !important breaks cascade architecture. Layers provide explicit priority without hacks.

### 5. Single Monolithic CSS File

**Don't:**
```
custom.css (2000+ lines)
```

**Do:**
```
layers.css
tokens.css
base.css
layout.css
components/*.css
utilities.css
```

**Why:** Small focused files are easier to maintain, debug, and understand. HTTP/2 makes multiple files performant.

## Browser Compatibility Notes

### backdrop-filter Support (2026)

**Supported:**
- Chrome 76+ (2019)
- Safari 9+ (2015, with -webkit- prefix)
- Firefox 103+ (2022)
- Edge 79+ (2020)

**Fallback for older browsers:**
```css
.glass-card {
    /* Fallback for no backdrop-filter support */
    background: var(--color-white);

    /* Progressive enhancement */
    @supports (backdrop-filter: blur(1px)) {
        background: var(--glass-bg);
        backdrop-filter: var(--glass-blur);
        -webkit-backdrop-filter: var(--glass-blur);
    }
}
```

### CSS Cascade Layers Support (2026)

**Supported:**
- Chrome 99+ (March 2022)
- Safari 15.4+ (March 2022)
- Firefox 97+ (February 2022)
- Edge 99+ (March 2022)

**Fallback:** Not needed. Layers degrade gracefully to normal cascade rules in older browsers.

### Native CSS Nesting Support (2026)

**Supported:**
- Chrome 112+ (April 2023)
- Safari 16.5+ (May 2023)
- Firefox 117+ (August 2023)

**Decision:** Use nesting sparingly. It's supported but not required. BEM reduces nesting need.

## Performance Considerations

### Glassmorphism Performance

**Backdrop-filter is GPU-accelerated but expensive:**
- Limit to 8-12 glass elements per page
- Use lighter blur on mobile (--glass-blur: blur(6px))
- Avoid animating backdrop-filter (animate transform instead)

**Optimization pattern:**
```css
.glass-card {
    /* Static blur, good performance */
    backdrop-filter: var(--glass-blur);

    /* Animate transform, not blur */
    transition: transform var(--transition-base);
}

.glass-card--interactive:hover {
    transform: translateY(-4px);
    /* Don't animate backdrop-filter! */
}
```

### Multiple CSS Files

**HTTP/2 advantage:** Multiple small files are fine. Browser parallelizes downloads.

**No concatenation needed:** Development and production use same files. Reduces build complexity.

**If optimization desired later:** CDN can handle minification without changing development workflow.

## Confidence Assessment

| Area | Confidence | Source |
|------|------------|--------|
| CSS Layers architecture | HIGH | MDN, CSS-Tricks (official docs) |
| BEM naming convention | HIGH | getbem.com (official spec) |
| Glassmorphism tokens | HIGH | Multiple 2026 guides, generator tools |
| Bento grid patterns | HIGH | Codemotion tutorial, multiple examples |
| Browser support | HIGH | MDN compatibility tables |
| Performance considerations | MEDIUM | Community best practices, no official benchmarks |

## Sources

**CSS Organization:**
- [Organizing your CSS - MDN](https://developer.mozilla.org/en-US/docs/Learn/CSS/Building_blocks/Organizing)
- [CSS Best Practices 2026 - Acodez](https://acodez.in/css-best-practices/)
- [The Modern CSS Toolkit 2026 - Nick Paolini](https://www.nickpaolini.com/blog/modern-css-toolkit-2026)

**CSS Cascade Layers:**
- [Cascade Layers Guide - CSS-Tricks](https://css-tricks.com/css-cascade-layers/)
- [Getting Started With CSS Cascade Layers - Smashing Magazine](https://www.smashingmagazine.com/2022/01/introduction-css-cascade-layers/)
- [Cascade layers - MDN](https://developer.mozilla.org/en-US/docs/Learn_web_development/Core/Styling_basics/Cascade_layers)

**BEM Naming:**
- [BEM Naming Convention - getbem.com](https://getbem.com/naming/)
- [BEM by Example - Sparkbox](https://sparkbox.com/foundry/bem_by_example)
- [Understanding BEM Methodology - GeeksforGeeks](https://www.geeksforgeeks.org/css/understanding-the-css-bem-convention/)

**Glassmorphism:**
- [How to Create Modern UI with Glassmorphism Effects 2026 - Kinetools](https://medium.com/@Kinetools/how-to-create-modern-ui-with-glassmorphism-effects-a-complete-2026-guide-2b1d71856542)
- [Glassmorphism in 2026 - Inverness Design Studio](https://invernessdesignstudio.com/glassmorphism-what-it-is-and-how-to-use-it-in-2026)
- [backdrop-filter - MDN](https://developer.mozilla.org/en-US/docs/Web/CSS/Reference/Properties/backdrop-filter)
- [Next-level frosted glass with backdrop-filter - Josh Comeau](https://www.joshwcomeau.com/css/backdrop-filter/)

**Bento Grids:**
- [Let's Create a Bento Box Design Layout - Codemotion](https://www.codemotion.com/magazine/frontend/lets-create-a-bento-box-design-layout-using-modern-css/)
- [Bento Grid with CSS Grid - iamsteve](https://iamsteve.me/blog/bento-layout-css-grid)
- [Bento Grid Generator](https://www.bentogridgenerator.com/)

**Design Tokens:**
- [Glass UI Generator](https://ui.glass/generator/)
- [CSS Glassmorphism Generator - 10015 Tools](https://10015.io/tools/css-glassmorphism-generator)
