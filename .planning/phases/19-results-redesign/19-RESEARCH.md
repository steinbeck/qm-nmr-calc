# Phase 19: Results Page Redesign - Research

**Researched:** 2026-01-29
**Domain:** Bento grid layout application, 3Dmol.js integration, glassmorphism cards, responsive design
**Confidence:** HIGH

## Summary

This phase applies the CSS foundation from Phase 18 to the results page, transforming the existing linear layout into a bento grid with glassmorphic cards. The results page has the most complex content mix: interactive 3D viewer (WebGL), spectrum images, data tables, metadata cards, and downloads.

Research confirms the Phase 18 CSS components are sufficient for this redesign. The bento grid system provides the necessary span utilities (2-col, 3-col, feature 2x2), glass cards provide visual containment, and responsive breakpoints handle mobile gracefully. The primary concern is 3Dmol.js WebGL canvas inside glass cards -- research shows this works but requires explicit height, `position: relative`, and potentially `isolation: isolate` on the container to avoid compositing issues with backdrop-filter.

The Jinja2 conditional templating already used in the current results.html handles ensemble vs single-conformer display. This pattern continues in the redesigned layout with conditional card visibility based on `conformer_mode`.

**Primary recommendation:** Use a 6-column bento grid with the 3D viewer as a span-4 feature card (not 2x2), spectra as span-2 cards side by side, and supporting content in span-2 or span-3 cards. Apply glass cards to all content sections except the 3D viewer container itself (which uses solid white for WebGL performance).

## Standard Stack

### Core (Already Available from Phase 18)

| Component | File | Purpose | Why Standard |
|-----------|------|---------|--------------|
| Bento Grid | `layout.css` | 6-column responsive grid with span utilities | Established in Phase 18, provides all needed span values |
| Glass Card | `glass-card.css` | Glassmorphic container component | BEM structure with header/body/footer elements |
| Design Tokens | `tokens.css` | Consistent spacing, colors, typography | Already includes all needed values |
| Utilities | `utilities.css` | Text alignment, display helpers | Provides `text-right`, `text-mono` |

### Supporting (Already Available)

| Component | File | Purpose | When to Use |
|-----------|------|---------|-------------|
| Legacy Styles | `legacy.css` | Modal, viewer container, legend styles | Keep for image modal, 3D viewer sizing |
| Base Styles | `base.css` | Table, form defaults | Table styling for shift data |

### New CSS Needed

| Component | Purpose | Why Needed |
|-----------|---------|------------|
| `results-page.css` | Page-specific overrides and viewer card | 3D viewer needs special treatment (solid bg, explicit height) |

**Installation:** No new dependencies. All CSS is already in place from Phase 18.

## Architecture Patterns

### Recommended Bento Grid Layout

```
Desktop (6 columns):
+---------------------------+---------------+---------------+
|                           |   1H Spectrum |  13C Spectrum |
|      3D Molecular Viewer  |    (span-2)   |    (span-2)   |
|         (span-4)          +---------------+---------------+
|                           |  Ensemble Metadata (span-2)   |
+---------------------------+---+---------------------------+
|    Calculation Details    |   | 1H Shift Table | 13C Shift Table |
|         (span-2)          |   |    (span-2)    |    (span-2)    |
+---------------------------+---+---------------+---------------+
|       Downloads (span-6 full width)                       |
+-----------------------------------------------------------+
```

**Rationale for span-4 (not 2x2 feature):**
- 3D viewer needs horizontal space more than vertical (molecule rotation)
- Allows spectra cards to align right side at same row start
- Easier responsive collapse (span-4 becomes full width on tablet/mobile)

### Responsive Collapse Strategy

**Tablet (3 columns):**
- 3D viewer: span-3 (full width)
- Spectra: span-3 each (stack vertically)
- All other cards: span-3 (full width)

**Mobile (1 column):**
- All cards: full width, natural stacking order
- Order is preserved from HTML source

### Pattern 1: 3D Viewer Card (Special Case)

**What:** WebGL canvas needs different treatment than glass cards
**When to use:** Any interactive 3D/WebGL content

**Why special:**
- `backdrop-filter` on parent can cause compositing performance issues
- Canvas needs explicit height (can't rely on aspect-ratio)
- 3Dmol.js requires `position: relative` on container

**Example:**
```html
<div class="bento-grid__item bento-grid__item--span-4 bento-grid__item--span-row-2">
    <div class="viewer-card">
        <header class="viewer-card__header">
            <h3 class="viewer-card__title">Interactive 3D Structure</h3>
            <!-- Conformer selector (ensemble only) -->
            {% if job.conformer_mode == 'ensemble' %}
            <select id="conformer-selector" class="viewer-card__selector">
                <!-- Populated by JS -->
            </select>
            {% endif %}
        </header>
        <div id="viewer-container-3d" class="viewer-card__canvas">
            <!-- 3Dmol.js renders here -->
        </div>
        <footer class="viewer-card__footer">
            <div class="viewer-legend">...</div>
        </footer>
    </div>
</div>
```

```css
/* results-page.css */
@layer components {
    .viewer-card {
        /* Solid background (NOT glass) for WebGL performance */
        background: var(--color-white);
        border: 1px solid var(--color-border);
        border-radius: var(--glass-radius);
        box-shadow: var(--glass-shadow);
        padding: var(--space-lg);
        display: flex;
        flex-direction: column;
        height: 100%;
    }

    .viewer-card__canvas {
        flex: 1;
        min-height: 350px;  /* Explicit minimum */
        position: relative;  /* Required by 3Dmol.js */
        isolation: isolate;  /* Prevent compositing issues */
        background: white;
        border-radius: var(--space-xs);
    }
}
```

**Source:** [3Dmol.js Tutorial](https://3dmol.csb.pitt.edu/doc/tutorial-code.html), [GitHub Issue #446](https://github.com/3dmol/3Dmol.js/issues/446)

### Pattern 2: Spectrum Image Cards

**What:** Medium-sized cards for spectrum images with click-to-enlarge
**When to use:** 1H and 13C spectrum display

**Example:**
```html
<div class="bento-grid__item bento-grid__item--span-2">
    <article class="glass-card glass-card--interactive">
        <header class="glass-card__header">
            <h3 class="glass-card__title"><sup>1</sup>H NMR Spectrum</h3>
        </header>
        <div class="glass-card__body">
            <figure class="spectrum-figure">
                <img src="/api/v1/jobs/{{ job.job_id }}/spectrum/1h.png"
                     alt="1H NMR spectrum"
                     class="spectrum-figure__image"
                     onclick="openModal(this.src)">
            </figure>
        </div>
    </article>
</div>
```

```css
@layer components {
    .spectrum-figure {
        margin: 0;
        text-align: center;
    }

    .spectrum-figure__image {
        max-width: 100%;
        height: auto;
        max-height: 280px;
        object-fit: contain;
        border-radius: var(--space-xs);
        cursor: pointer;
        transition: transform var(--transition-base), box-shadow var(--transition-base);
    }

    @media (hover: hover) {
        .spectrum-figure__image:hover {
            transform: scale(1.02);
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
        }
    }
}
```

### Pattern 3: Chemical Shift Tables in Glass Cards

**What:** Right-aligned numeric tables for shift data
**When to use:** 1H and 13C chemical shift listings

**Example:**
```html
<div class="bento-grid__item bento-grid__item--span-2">
    <article class="glass-card">
        <header class="glass-card__header">
            <h3 class="glass-card__title"><sup>1</sup>H Chemical Shifts</h3>
        </header>
        <div class="glass-card__body">
            <table class="shift-table">
                <thead>
                    <tr>
                        <th>Atom</th>
                        <th class="text-right">Shift (ppm)</th>
                    </tr>
                </thead>
                <tbody id="h1-shifts-body">
                    <!-- Populated by JS -->
                </tbody>
            </table>
        </div>
    </article>
</div>
```

```css
@layer components {
    .shift-table {
        width: 100%;
        border-collapse: collapse;
    }

    .shift-table th,
    .shift-table td {
        padding: var(--space-sm) var(--space-md);
        text-align: left;
        border-bottom: 1px solid var(--color-border);
    }

    .shift-table th {
        font-weight: 600;
        color: var(--color-text-muted);
        font-size: var(--text-sm);
    }

    /* Right-align numeric columns with tabular figures */
    .shift-table td:last-child,
    .shift-table th:last-child {
        text-align: right;
        font-variant-numeric: tabular-nums;
    }

    .shift-table td:last-child {
        font-family: var(--font-mono);
    }
}
```

**Source:** [MDN font-variant-numeric](https://developer.mozilla.org/en-US/docs/Web/CSS/Reference/Properties/font-variant-numeric)

### Pattern 4: Ensemble Metadata Card (Conditional)

**What:** Card showing conformer count, populations, generation method
**When to use:** Only for ensemble jobs (multi-conformer)

**Example:**
```html
{% if job.conformer_mode == 'ensemble' %}
<div class="bento-grid__item bento-grid__item--span-2">
    <article class="glass-card glass-card--featured" id="ensemble-metadata">
        <header class="glass-card__header">
            <h3 class="glass-card__title">Ensemble Summary</h3>
        </header>
        <div class="glass-card__body">
            <dl class="metadata-list">
                <div class="metadata-list__item">
                    <dt>Conformers</dt>
                    <dd id="meta-conformer-count">-</dd>
                </div>
                <div class="metadata-list__item">
                    <dt>Method</dt>
                    <dd id="meta-method">-</dd>
                </div>
                <!-- More metadata... -->
            </dl>
            <!-- Top conformers table -->
            <h4>Top Contributors</h4>
            <table class="conformer-table" id="top-conformers-table">
                <!-- Populated by JS -->
            </table>
        </div>
    </article>
</div>
{% endif %}
```

### Pattern 5: Conformer Selector in Viewer Card

**What:** Dropdown integrated into viewer card header
**When to use:** Ensemble jobs only, allows switching displayed geometry

**Example:**
```html
<header class="viewer-card__header">
    <h3 class="viewer-card__title">Interactive 3D Structure</h3>
    {% if job.conformer_mode == 'ensemble' %}
    <div class="viewer-card__controls">
        <label for="conformer-selector" class="sr-only">View Conformer</label>
        <select id="conformer-selector" class="conformer-select">
            <!-- Populated by JS -->
        </select>
        <span id="current-conformer-info" class="conformer-info"></span>
    </div>
    {% endif %}
</header>
```

```css
@layer components {
    .viewer-card__header {
        display: flex;
        justify-content: space-between;
        align-items: center;
        gap: var(--space-md);
        margin-bottom: var(--space-md);
        flex-wrap: wrap;
    }

    .viewer-card__controls {
        display: flex;
        align-items: center;
        gap: var(--space-sm);
    }

    .conformer-select {
        /* Solid background for form inputs (accessibility) */
        background: var(--color-white);
        border: 1px solid var(--color-border);
        border-radius: var(--space-xs);
        padding: var(--space-xs) var(--space-sm);
        font-size: var(--text-sm);
        min-width: 200px;
    }

    .conformer-info {
        font-size: var(--text-sm);
        color: var(--color-text-muted);
        font-style: italic;
    }
}
```

### Pattern 6: Downloads Card (Compact)

**What:** Full-width card with button grid for downloads
**When to use:** Bottom of results page

**Example:**
```html
<div class="bento-grid__item bento-grid__item--span-6">
    <article class="glass-card glass-card--subtle">
        <header class="glass-card__header">
            <h3 class="glass-card__title">Downloads</h3>
        </header>
        <div class="glass-card__body">
            <div class="download-grid">
                <a href="/api/v1/jobs/{{ job.job_id }}/geometry" class="download-btn">
                    <span class="download-btn__icon">...</span>
                    <span class="download-btn__label">Geometry (XYZ)</span>
                </a>
                <!-- More download buttons -->
            </div>
        </div>
    </article>
</div>
```

### Anti-Patterns to Avoid

- **Glass card on WebGL container:** Don't apply backdrop-filter directly on the 3Dmol.js container -- use solid white background for performance
- **Fixed pixel heights on grid items:** Use min-height + flex-grow so cards adapt to content
- **Inline styles for conditional display:** Use Jinja2 `{% if %}` blocks to exclude entire cards from HTML
- **Duplicating responsive breakpoints:** Let bento grid handle responsive -- don't add custom breakpoints in results-page.css

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Grid layout | Custom grid CSS | Phase 18 `layout.css` bento grid | Already tested, responsive breakpoints included |
| Card styling | New card styles | Phase 18 `glass-card.css` | BEM structure, modifiers, fallbacks all present |
| Table alignment | Manual padding | `font-variant-numeric: tabular-nums` + `text-align: right` | Native CSS feature, cleaner code |
| Conditional cards | JS show/hide | Jinja2 `{% if %}` blocks | Cards not in DOM if not needed, better performance |
| Image modal | New modal system | Existing `legacy.css` dialog styles | Already working, tested with Pico dialog |

**Key insight:** Phase 18 did the hard work. This phase is about applying existing components to new content, not creating new components.

## Common Pitfalls

### Pitfall 1: 3Dmol.js Container Has No Height

**What goes wrong:** 3D viewer appears as collapsed (0 height) or invisible.

**Why it happens:** 3Dmol.js requires the container div to have explicit dimensions. CSS Grid can give width via span, but height needs to be set explicitly.

**How to avoid:**
- Set `min-height: 350px` on the canvas container
- Use `flex: 1` on canvas inside flex column to fill available space
- Ensure `position: relative` on container (3Dmol.js requirement)

**Warning signs:**
- Viewer div has 0 computed height
- 3Dmol.js logs "Unable to create viewer" or similar

**Source:** [3Dmol.js Tutorial](https://3dmol.csb.pitt.edu/doc/tutorial-code.html)

### Pitfall 2: Backdrop-filter Performance with WebGL

**What goes wrong:** Sluggish rotation/zoom in 3D viewer, dropped frames, high GPU usage.

**Why it happens:** backdrop-filter causes continuous GPU compositing. When combined with WebGL canvas rendering, GPU contention occurs.

**How to avoid:**
- Use solid background on viewer card (NOT glass-card class)
- Add `isolation: isolate` to canvas container
- Limit glass elements near WebGL content

**Warning signs:**
- Frame rate drops during rotation
- Browser DevTools shows high GPU usage
- Viewer stutters when scrolling page

**Source:** [Shadcn UI backdrop-filter issue](https://github.com/shadcn-ui/ui/issues/327), [Mozilla Bug 1718471](https://bugzilla.mozilla.org/show_bug.cgi?id=1718471)

### Pitfall 3: Spectra Images Stretched or Cropped

**What goes wrong:** Spectrum PNGs don't fit card dimensions properly.

**Why it happens:** Fixed aspect ratio images in flexible containers need `object-fit` and max dimensions.

**How to avoid:**
- Use `object-fit: contain` on images
- Set `max-height` to prevent overflow
- Allow `height: auto` for natural scaling

**Warning signs:**
- Images appear distorted
- Scrollbars appear on card
- Image aspect ratio looks wrong

### Pitfall 4: Ensemble Cards Appear for Single-Conformer Jobs

**What goes wrong:** Empty "Ensemble Summary" card visible on single-conformer results.

**Why it happens:** JavaScript populates cards dynamically, leaving empty containers if condition not in template.

**How to avoid:**
- Use Jinja2 `{% if job.conformer_mode == 'ensemble' %}` to exclude cards from DOM entirely
- Don't rely on JS to hide/show entire cards
- JS only updates content within cards that exist

**Warning signs:**
- Empty cards with headers but no content
- "Conformer Selector" visible on single-conformer job
- Layout has unexpected gaps

### Pitfall 5: Mobile Layout Order Wrong

**What goes wrong:** Important content (3D viewer) appears below less important content on mobile.

**Why it happens:** CSS Grid flattens to source order on mobile. If HTML order doesn't match importance, mobile suffers.

**How to avoid:**
- Order HTML by importance: 3D viewer first, then spectra, then tables
- Don't rely on CSS Grid placement for visual hierarchy
- Test mobile layout early

**Warning signs:**
- Downloads appear before results on mobile
- Scrolling past metadata to see spectra

## Code Examples

### Complete Results Page HTML Structure

```html
{% extends "base.html" %}

{% block title %}Results: {{ job.input_name or job.job_id }} | QM NMR Calculator{% endblock %}

{% block head %}
<script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.5.3/3Dmol-min.js"></script>
<link rel="stylesheet" href="{{ url_for('static', path='/css/results-page.css') }}">
{% endblock %}

{% block content %}
<div class="page-wrapper">
    <header class="results-header">
        <h1>NMR Calculation Results</h1>
        <p class="muted">{{ job.input_name or job.job_id }}</p>
    </header>

    <div class="bento-grid bento-grid--dense">
        <!-- 3D Viewer (Hero Position) -->
        <div class="bento-grid__item bento-grid__item--span-4 bento-grid__item--span-row-2">
            <div class="viewer-card">
                <!-- Content... -->
            </div>
        </div>

        <!-- 1H Spectrum -->
        <div class="bento-grid__item bento-grid__item--span-2">
            <article class="glass-card glass-card--interactive">
                <!-- Content... -->
            </article>
        </div>

        <!-- 13C Spectrum -->
        <div class="bento-grid__item bento-grid__item--span-2">
            <article class="glass-card glass-card--interactive">
                <!-- Content... -->
            </article>
        </div>

        <!-- Ensemble Metadata (conditional) -->
        {% if job.conformer_mode == 'ensemble' %}
        <div class="bento-grid__item bento-grid__item--span-2">
            <article class="glass-card glass-card--featured" id="ensemble-metadata">
                <!-- Content... -->
            </article>
        </div>
        {% endif %}

        <!-- Calculation Details -->
        <div class="bento-grid__item bento-grid__item--span-2">
            <article class="glass-card">
                <!-- Content... -->
            </article>
        </div>

        <!-- 1H Shift Table -->
        <div class="bento-grid__item bento-grid__item--span-2">
            <article class="glass-card">
                <!-- Content... -->
            </article>
        </div>

        <!-- 13C Shift Table -->
        <div class="bento-grid__item bento-grid__item--span-2">
            <article class="glass-card">
                <!-- Content... -->
            </article>
        </div>

        <!-- Downloads (Full Width) -->
        <div class="bento-grid__item bento-grid__item--span-6">
            <article class="glass-card glass-card--subtle">
                <!-- Content... -->
            </article>
        </div>
    </div>

    <p class="back-link"><a href="/">Submit another job</a></p>
</div>

<!-- Image Modal (reuse from legacy) -->
<dialog id="image-modal">
    <!-- Same as current -->
</dialog>
{% endblock %}
```

### results-page.css Component CSS

```css
/* results-page.css */
/* Page-specific components for results page bento layout */

@layer components {
    /* ===== Results Header ===== */
    .results-header {
        margin-bottom: var(--space-xl);
    }

    .results-header h1 {
        margin: 0;
        font-size: var(--text-2xl);
    }

    /* ===== Viewer Card (Solid Background) ===== */
    .viewer-card {
        background: var(--color-white);
        border: 1px solid var(--color-border);
        border-radius: var(--glass-radius);
        box-shadow: var(--glass-shadow);
        padding: var(--space-lg);
        display: flex;
        flex-direction: column;
        height: 100%;
    }

    .viewer-card__header {
        display: flex;
        justify-content: space-between;
        align-items: center;
        gap: var(--space-md);
        margin-bottom: var(--space-md);
        padding-bottom: var(--space-md);
        border-bottom: 1px solid var(--color-border);
        flex-wrap: wrap;
    }

    .viewer-card__title {
        font-size: var(--text-xl);
        font-weight: 600;
        margin: 0;
    }

    .viewer-card__controls {
        display: flex;
        align-items: center;
        gap: var(--space-sm);
    }

    .viewer-card__canvas {
        flex: 1;
        min-height: 350px;
        position: relative;
        isolation: isolate;
        background: white;
        border-radius: var(--space-xs);
        border: 1px solid var(--color-gray-200);
    }

    .viewer-card__footer {
        margin-top: var(--space-md);
        padding-top: var(--space-md);
        border-top: 1px solid var(--color-border);
    }

    /* ===== Conformer Selector ===== */
    .conformer-select {
        background: var(--color-white);
        border: 1px solid var(--color-border);
        border-radius: var(--space-xs);
        padding: var(--space-xs) var(--space-sm);
        font-size: var(--text-sm);
        min-width: 200px;
        cursor: pointer;
    }

    .conformer-select:focus {
        outline: 2px solid var(--color-primary);
        outline-offset: 2px;
    }

    .conformer-info {
        font-size: var(--text-sm);
        color: var(--color-text-muted);
        font-style: italic;
    }

    /* ===== Spectrum Figure ===== */
    .spectrum-figure {
        margin: 0;
        text-align: center;
    }

    .spectrum-figure__image {
        max-width: 100%;
        height: auto;
        max-height: 280px;
        object-fit: contain;
        border-radius: var(--space-xs);
        cursor: pointer;
        transition: transform var(--transition-base), box-shadow var(--transition-base);
    }

    @media (hover: hover) {
        .spectrum-figure__image:hover {
            transform: scale(1.02);
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
        }
    }

    /* ===== Shift Tables ===== */
    .shift-table {
        width: 100%;
        border-collapse: collapse;
    }

    .shift-table th,
    .shift-table td {
        padding: var(--space-sm) var(--space-md);
        text-align: left;
        border-bottom: 1px solid var(--color-border);
    }

    .shift-table th {
        font-weight: 600;
        color: var(--color-text-muted);
        font-size: var(--text-sm);
        text-transform: uppercase;
        letter-spacing: 0.5px;
    }

    .shift-table td:last-child,
    .shift-table th:last-child {
        text-align: right;
        font-variant-numeric: tabular-nums;
    }

    .shift-table td:last-child {
        font-family: var(--font-mono);
    }

    .shift-table tbody tr:last-child td {
        border-bottom: none;
    }

    /* ===== Metadata List ===== */
    .metadata-list {
        display: grid;
        grid-template-columns: repeat(2, 1fr);
        gap: var(--space-sm);
    }

    .metadata-list__item {
        display: flex;
        flex-direction: column;
    }

    .metadata-list__item dt {
        font-size: var(--text-sm);
        color: var(--color-text-muted);
        margin-bottom: var(--space-xs);
    }

    .metadata-list__item dd {
        margin: 0;
        font-weight: 500;
    }

    /* ===== Download Grid ===== */
    .download-grid {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(140px, 1fr));
        gap: var(--space-sm);
    }

    .download-btn {
        display: flex;
        flex-direction: column;
        align-items: center;
        gap: var(--space-xs);
        padding: var(--space-md);
        background: var(--color-gray-50);
        border: 1px solid var(--color-border);
        border-radius: var(--space-sm);
        text-decoration: none;
        color: var(--color-text);
        transition: background var(--transition-fast), border-color var(--transition-fast);
    }

    .download-btn:hover {
        background: var(--color-gray-100);
        border-color: var(--color-primary);
    }

    .download-btn__label {
        font-size: var(--text-sm);
        text-align: center;
    }

    /* ===== Ensemble Note ===== */
    .ensemble-note {
        padding: var(--space-md);
        background: var(--color-primary-light);
        border-left: 4px solid var(--color-primary);
        border-radius: var(--space-xs);
        font-size: var(--text-sm);
        margin-top: var(--space-md);
    }

    /* ===== Conformer Table ===== */
    .conformer-table {
        width: 100%;
        border-collapse: collapse;
        margin-top: var(--space-sm);
        font-size: var(--text-sm);
    }

    .conformer-table th,
    .conformer-table td {
        padding: var(--space-xs) var(--space-sm);
        text-align: left;
        border-bottom: 1px solid var(--color-border);
    }

    .conformer-table th {
        font-weight: 600;
        color: var(--color-text-muted);
    }

    .conformer-table td:nth-child(2),
    .conformer-table td:nth-child(3) {
        text-align: right;
        font-variant-numeric: tabular-nums;
    }

    /* ===== Back Link ===== */
    .back-link {
        margin-top: var(--space-xl);
        text-align: center;
    }
}
```

### JavaScript Changes Required

The existing JavaScript in results.html will work with minimal changes:

1. **Selector updates:** Change element IDs/classes to match new HTML structure
2. **Conformer selector:** Already works, just needs CSS styling
3. **Image modal:** Reuse existing `openModal()` function
4. **3Dmol.js:** Container ID stays `viewer-container-3d`, no changes needed

```javascript
// Existing JS remains functional
// Only update references if element IDs change

// Example: Updating a selector reference
document.getElementById('conformer-selector-container')
// becomes just checking if selector exists (Jinja2 conditionally renders)
const selectorContainer = document.getElementById('conformer-selector');
if (selectorContainer) {
    // Ensemble mode
}
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `<article>` semantic blocks | BEM glass-card components | Phase 18 | Consistent styling, clear component boundaries |
| 2-column spectra grid | 6-column bento grid | Phase 18 | Flexible asymmetric layouts |
| Inline style display:none | Jinja2 {% if %} blocks | Best practice | Cleaner DOM, no hidden elements |
| Fixed height viewers | Flex + min-height | Modern CSS | Responsive, content-adaptive |

**Deprecated/outdated:**
- Pico CSS classes (`.secondary`, `role="button"`) -- replace with glass-card and custom button styles
- `article.results-container` max-width -- bento grid handles layout
- Inline styles for conditional display -- use Jinja2 conditionals

## Open Questions

### 1. Exact Card Sizing for Spectra

**What we know:**
- Spectrum PNGs are generated by matplotlib
- Current max-height is 400px in legacy.css
- Glass cards have padding that reduces available space

**What's unclear:**
- Optimal max-height for spectrum images in bento cards
- Whether aspect ratio constraint is needed

**Recommendation:**
- Start with max-height: 280px (accounts for card padding)
- Test with actual spectrum PNGs
- Adjust in implementation if needed

### 2. Dense Grid Packing Behavior

**What we know:**
- `bento-grid--dense` uses `grid-auto-flow: dense`
- Dense packing can reorder visual display vs source order

**What's unclear:**
- Whether dense packing causes accessibility issues (visual order vs DOM order)
- If ensemble metadata card would shift unexpectedly

**Recommendation:**
- Use dense packing initially (helps fill gaps)
- Test with screen reader to verify reading order is logical
- Remove if reordering is problematic

### 3. Mobile Touch Target Size for Download Buttons

**What we know:**
- WCAG recommends 44x44px minimum touch targets
- Current download buttons may be smaller on mobile

**What's unclear:**
- Exact padding needed for WCAG compliance
- Whether button text truncates on small screens

**Recommendation:**
- Set min-height: 44px on download buttons
- Test on mobile device during implementation
- Adjust padding if text wraps poorly

## Sources

### Primary (HIGH confidence)

- [3Dmol.js Tutorial](https://3dmol.csb.pitt.edu/doc/tutorial-code.html) - Official documentation for viewer setup
- [3Dmol.js Issue #446](https://github.com/3dmol/3Dmol.js/issues/446) - Container sizing requirements
- [MDN font-variant-numeric](https://developer.mozilla.org/en-US/docs/Web/CSS/Reference/Properties/font-variant-numeric) - Tabular numbers for tables
- [FastAPI Templates](https://fastapi.tiangolo.com/advanced/templates/) - Jinja2 integration

### Secondary (MEDIUM confidence)

- [FreeCodeCamp Bento Grids](https://www.freecodecamp.org/news/bento-grids-in-web-design/) - Bento grid best practices
- [Codemotion Bento Layout](https://www.codemotion.com/magazine/frontend/lets-create-a-bento-box-design-layout-using-modern-css/) - CSS Grid implementation patterns
- [iamsteve Bento Layout](https://iamsteve.me/blog/bento-layout-css-grid) - Grid template areas usage
- [Shadcn UI backdrop-filter issue](https://github.com/shadcn-ui/ui/issues/327) - Performance with WebGL

### Tertiary (LOW confidence - requires validation)

- [Mozilla Bug 1718471](https://bugzilla.mozilla.org/show_bug.cgi?id=1718471) - Firefox backdrop-filter performance (may be fixed)
- Bento grid sizing recommendations (based on general best practices, not empirically tested with this content)

## Metadata

**Confidence breakdown:**
- Bento grid layout: HIGH - Phase 18 CSS already provides all needed utilities
- 3D viewer integration: HIGH - 3Dmol.js documentation clear on requirements
- Glass card application: HIGH - Phase 18 components ready, just need assembly
- Responsive behavior: MEDIUM - Phase 18 breakpoints exist, but specific card ordering needs testing
- Performance with WebGL: MEDIUM - Known concern, mitigation strategies documented

**Research date:** 2026-01-29
**Valid until:** 2026-04-29 (90 days - stable CSS patterns, unlikely to change)
