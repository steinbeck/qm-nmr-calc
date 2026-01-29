# Phase 21: Status Page Redesign - Research

**Researched:** 2026-01-29
**Domain:** Job progress visualization, CSS step trackers, status indicators, progress bars
**Confidence:** HIGH

## Summary

This phase redesigns the job status page to provide clear visual progress tracking during NMR calculations. The current implementation has basic progress elements but lacks the refined styling and clear visual hierarchy established in Phases 18-20. Research confirms that all required functionality can be achieved with pure CSS using the existing design system, with the native HTML `<progress>` element as the foundation for progress bars.

The key insight is that the existing status.html already has the JavaScript polling infrastructure and data binding in place. The redesign focuses on visual presentation: applying glass-card styling, creating step tracker components with clear state indicators (completed/active/pending), improving the conformer progress display, and adding the 3D molecule preview that already works on other pages.

Research into status indicator patterns (Carbon Design System) confirms that status indicators should use at least two of: color, shape, or symbol for accessibility. The existing status display uses only color, which must be enhanced with icons or shapes. The HTML `<progress>` element requires careful cross-browser styling with vendor prefixes but provides built-in accessibility semantics.

**Primary recommendation:** Apply existing glass-card and bento-grid components to the status page layout. Create a step-tracker component with BEM naming that shows completed (checkmark), active (spinner), and pending (circle) states using CSS-only icons. Use the native `<progress>` element for conformer progress with custom styling. Reuse the viewer-card pattern for the 3D preview (solid white background for WebGL performance).

## Standard Stack

### Core (Already Available from Phases 18-20)

| Component | File | Purpose | Why Standard |
|-----------|------|---------|--------------|
| Glass Card | `glass-card.css` | Container for status sections | BEM structure, accessibility fallbacks |
| Bento Grid | `layout.css` | Page layout structure | Responsive breakpoints established |
| Design Tokens | `tokens.css` | Colors for status states | `--color-success`, `--color-warning`, `--color-error` defined |
| Viewer Card | `results-page.css` | 3D molecule preview | Solid white for WebGL, already tested |
| 3Dmol.js | CDN 2.5.3 | Interactive molecule viewer | Already used in status.html |

### Supporting

| Element | Source | Purpose | When to Use |
|---------|--------|---------|-------------|
| `<progress>` | HTML native | Progress bar semantics | Conformer progress display |
| CSS animations | Native | Spinner, pulse effects | Active step indicator |
| Unicode symbols | Native | Status icons | Checkmark, circle, spinner fallback |

### New CSS Needed

| Component | Purpose | Why Needed |
|-----------|---------|------------|
| `status-page.css` | Step tracker, progress bar styling, status cards | Page-specific components following established patterns |

**Installation:** No new dependencies. All features use existing CSS infrastructure and native HTML.

## Architecture Patterns

### Recommended Bento Grid Layout

```
Desktop (6 columns):
+----------------------------------+---------------------------+
|         Job Summary Card         |     3D Molecule Preview   |
|           (span-3)               |         (span-3)          |
+----------------------------------+---------------------------+
|              Step Progress Tracker (span-6)                  |
+--------------------------------------------------------------+
|            Conformer Progress Section (span-6)               |
|   [Progress bar + conformer table (ensemble only)]           |
+--------------------------------------------------------------+
|              Error Display (span-6, if failed)               |
+--------------------------------------------------------------+
```

**Rationale:**
- Job summary and 3D preview side-by-side (50/50 split)
- Step tracker full width for clear linear progression
- Conformer progress full width with expandable details
- Error section appears below all progress when needed

### Pattern 1: Step Tracker Component (Pure CSS)

**What:** Horizontal timeline showing completed, active, and pending calculation steps
**When to use:** Display the multi-step calculation process

**Example HTML:**
```html
<div class="step-tracker">
    <ol class="step-tracker__list">
        <li class="step-tracker__item step-tracker__item--complete">
            <span class="step-tracker__icon" aria-hidden="true"></span>
            <span class="step-tracker__label">Geometry Optimization</span>
            <span class="step-tracker__duration">45s</span>
        </li>
        <li class="step-tracker__item step-tracker__item--active">
            <span class="step-tracker__icon" aria-hidden="true"></span>
            <span class="step-tracker__label">NMR Shielding</span>
            <span class="step-tracker__duration">Running...</span>
        </li>
        <li class="step-tracker__item step-tracker__item--pending">
            <span class="step-tracker__icon" aria-hidden="true"></span>
            <span class="step-tracker__label">Post-processing</span>
        </li>
    </ol>
</div>
```

**Example CSS:**
```css
@layer components {
    .step-tracker {
        padding: var(--space-lg);
    }

    .step-tracker__list {
        display: flex;
        justify-content: space-between;
        list-style: none;
        padding: 0;
        margin: 0;
        position: relative;
    }

    /* Connecting line between steps */
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
        position: relative;
        z-index: 1;
        flex: 1;
        text-align: center;
    }

    .step-tracker__icon {
        width: 24px;
        height: 24px;
        border-radius: 50%;
        display: flex;
        align-items: center;
        justify-content: center;
        font-size: 14px;
        background: var(--color-white);
        border: 2px solid var(--color-border);
    }

    /* Complete state: Green with checkmark */
    .step-tracker__item--complete .step-tracker__icon {
        background: var(--color-success);
        border-color: var(--color-success);
        color: var(--color-white);
    }

    .step-tracker__item--complete .step-tracker__icon::before {
        content: '\2713';  /* Checkmark */
    }

    /* Active state: Blue with spinner animation */
    .step-tracker__item--active .step-tracker__icon {
        background: var(--color-primary);
        border-color: var(--color-primary);
        color: var(--color-white);
        animation: pulse 1.5s ease-in-out infinite;
    }

    .step-tracker__item--active .step-tracker__icon::before {
        content: '\25CF';  /* Filled circle as spinner fallback */
    }

    /* Pending state: Gray outline */
    .step-tracker__item--pending .step-tracker__icon {
        background: var(--color-white);
        border-color: var(--color-gray-300);
    }

    .step-tracker__item--pending .step-tracker__icon::before {
        content: '\25CB';  /* Empty circle */
        color: var(--color-gray-400);
    }

    .step-tracker__label {
        font-size: var(--text-sm);
        font-weight: 500;
        color: var(--color-text);
    }

    .step-tracker__item--pending .step-tracker__label {
        color: var(--color-text-muted);
    }

    .step-tracker__duration {
        font-size: var(--text-xs);
        color: var(--color-text-muted);
    }

    @keyframes pulse {
        0%, 100% { transform: scale(1); opacity: 1; }
        50% { transform: scale(1.1); opacity: 0.8; }
    }
}

/* Mobile: Stack vertically */
@media (max-width: 767px) {
    .step-tracker__list {
        flex-direction: column;
        align-items: flex-start;
        gap: var(--space-md);
    }

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
    }
}
```

**Source:** Adapted from [CSS Progress Wizard](https://christabor.github.io/css-progress-wizard/), [Progress Tracker by Nigel O'Toole](https://nigelotoole.github.io/progress-tracker/)

### Pattern 2: Conformer Progress Bar (Native Element)

**What:** Styled `<progress>` element for ensemble job conformer completion
**When to use:** Showing X/N conformers complete

**Example HTML:**
```html
<div class="conformer-progress">
    <div class="conformer-progress__header">
        <span class="conformer-progress__stage">NMR Calculations</span>
        <span class="conformer-progress__count">5 / 12 conformers</span>
    </div>
    <progress
        class="conformer-progress__bar"
        value="5"
        max="12"
        aria-label="Conformer progress: 5 of 12 complete">
    </progress>
    <div class="conformer-progress__eta">
        Estimated time remaining: ~8 minutes
    </div>
</div>
```

**Example CSS:**
```css
@layer components {
    .conformer-progress {
        padding: var(--space-md);
    }

    .conformer-progress__header {
        display: flex;
        justify-content: space-between;
        align-items: baseline;
        margin-bottom: var(--space-sm);
    }

    .conformer-progress__stage {
        font-weight: 600;
        color: var(--color-text);
    }

    .conformer-progress__count {
        font-size: var(--text-sm);
        color: var(--color-text-muted);
        font-variant-numeric: tabular-nums;
    }

    /* Progress bar styling - requires vendor prefixes */
    .conformer-progress__bar {
        width: 100%;
        height: 1.5rem;
        appearance: none;
        -webkit-appearance: none;
        border: none;
        border-radius: var(--space-xs);
        background: var(--color-gray-100);
        overflow: hidden;
    }

    /* WebKit (Chrome, Safari, Edge) */
    .conformer-progress__bar::-webkit-progress-bar {
        background: var(--color-gray-100);
        border-radius: var(--space-xs);
    }

    .conformer-progress__bar::-webkit-progress-value {
        background: var(--color-primary);
        border-radius: var(--space-xs);
        transition: width var(--transition-base);
    }

    /* Firefox */
    .conformer-progress__bar::-moz-progress-bar {
        background: var(--color-primary);
        border-radius: var(--space-xs);
    }

    .conformer-progress__eta {
        margin-top: var(--space-sm);
        font-size: var(--text-sm);
        color: var(--color-text-muted);
    }
}
```

**Source:** [MDN progress element](https://developer.mozilla.org/en-US/docs/Web/HTML/Reference/Elements/progress), [Custom CSS Progress Bar](https://nikitahl.com/progress-bar-css)

### Pattern 3: Status Summary Card

**What:** Glass card showing job metadata and current status
**When to use:** Primary job information display

**Example HTML:**
```html
<article class="glass-card">
    <header class="glass-card__header">
        <h2 class="glass-card__title">Job Status</h2>
        <span class="status-badge status-badge--running">Running</span>
    </header>
    <div class="glass-card__body">
        <dl class="metadata-list">
            <div class="metadata-list__item">
                <dt>Job ID</dt>
                <dd><code>a1b2c3d4e5f6</code></dd>
            </div>
            <div class="metadata-list__item">
                <dt>Elapsed Time</dt>
                <dd id="elapsed-time">2m 34s</dd>
            </div>
            <!-- More metadata... -->
        </dl>
    </div>
</article>
```

**Example CSS:**
```css
@layer components {
    .status-badge {
        display: inline-flex;
        align-items: center;
        gap: var(--space-xs);
        padding: var(--space-xs) var(--space-sm);
        border-radius: var(--space-xs);
        font-size: var(--text-sm);
        font-weight: 500;
    }

    /* Status badge variants */
    .status-badge--queued {
        background: var(--color-gray-100);
        color: var(--color-gray-600);
    }

    .status-badge--running {
        background: hsl(195 85% 95%);
        color: var(--color-primary);
    }

    .status-badge--complete {
        background: hsl(142 76% 95%);
        color: var(--color-success);
    }

    .status-badge--failed {
        background: hsl(0 84% 95%);
        color: var(--color-error);
    }

    /* Status badge icon (dot) */
    .status-badge::before {
        content: '';
        width: 8px;
        height: 8px;
        border-radius: 50%;
        background: currentColor;
    }

    .status-badge--running::before {
        animation: pulse 1.5s ease-in-out infinite;
    }
}
```

### Pattern 4: Error Display Card

**What:** Alert-style card for displaying job failure information
**When to use:** When job.status === 'failed'

**Example HTML:**
```html
<article class="glass-card glass-card--error">
    <header class="glass-card__header">
        <h3 class="glass-card__title">
            <span class="error-icon" aria-hidden="true">!</span>
            Job Failed
        </h3>
    </header>
    <div class="glass-card__body">
        <p class="error-message" id="error-message">
            NWChem geometry optimization did not converge after 100 cycles.
        </p>
        <div class="error-actions">
            <a href="/" class="btn">Try Another Molecule</a>
        </div>
    </div>
</article>
```

**Example CSS:**
```css
@layer components {
    .glass-card--error {
        background: hsl(0 84% 97%);
        border-color: hsl(0 84% 80%);
    }

    .glass-card--error .glass-card__header {
        border-bottom-color: hsl(0 84% 80%);
    }

    .glass-card--error .glass-card__title {
        color: var(--color-error);
        display: flex;
        align-items: center;
        gap: var(--space-sm);
    }

    .error-icon {
        width: 24px;
        height: 24px;
        border-radius: 50%;
        background: var(--color-error);
        color: var(--color-white);
        display: flex;
        align-items: center;
        justify-content: center;
        font-weight: bold;
    }

    .error-message {
        color: var(--color-text);
        line-height: 1.6;
        margin-bottom: var(--space-lg);
    }

    .error-actions {
        display: flex;
        gap: var(--space-md);
    }
}
```

### Pattern 5: 3D Molecule Preview (Reuse from Results)

**What:** Interactive 3D viewer showing molecule geometry during calculation
**When to use:** STATUS-03 requirement - display geometry if available

**Reuse existing pattern from `results-page.css`:**
- `viewer-card` with solid white background (not glass)
- `viewer-card__canvas` with `position: relative` and `isolation: isolate`
- `min-height: 350px` but can be reduced to `min-height: 250px` for status page

The existing status.html already initializes 3Dmol.js and loads geometry - just need to apply proper styling.

### Anti-Patterns to Avoid

- **Glass card on WebGL container:** Continue using solid white background for 3D viewer
- **Color-only status indicators:** Always pair color with shape or symbol for accessibility
- **Animating backdrop-filter:** Use transform/opacity animations only
- **Custom progress bar divs:** Use native `<progress>` element for accessibility semantics
- **Hiding conformer progress with CSS:** Use Jinja2 conditionals to exclude from DOM when not ensemble mode

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Progress bar semantics | Custom div with width% | Native `<progress>` element | Built-in ARIA role, screen reader support |
| Step tracker layout | Manual flex calculations | Existing flexbox pattern | Well-tested responsive behavior |
| Status colors | Hardcoded hex values | Design tokens (`--color-success`, etc.) | Consistent with design system |
| 3D viewer container | New styling approach | `viewer-card` from results-page.css | Already solves WebGL/glass conflict |
| Conformer table | New table styling | `conformer-table` from results-page.css | Already has tabular-nums, proper alignment |
| Polling infrastructure | New JS implementation | Existing status.html JavaScript | Already handles all data binding |

**Key insight:** The status page already has working functionality. This phase is about applying existing CSS patterns to improve visual presentation, not building new features.

## Common Pitfalls

### Pitfall 1: Progress Element Styling is Browser-Specific

**What goes wrong:** Progress bar looks different or unstyled across browsers.

**Why it happens:** The `<progress>` element has different internal structure in WebKit vs Firefox. CSS must target vendor-specific pseudo-elements.

**How to avoid:**
- Always set `appearance: none` / `-webkit-appearance: none` first
- Style both `::-webkit-progress-bar` AND `::-webkit-progress-value` for WebKit
- Style `::-moz-progress-bar` for Firefox
- Provide fallback background color on the element itself

**Warning signs:**
- Default browser chrome visible
- Progress bar has no fill color
- Styling works in Chrome but not Firefox/Safari

**Source:** [MDN progress element](https://developer.mozilla.org/en-US/docs/Web/HTML/Reference/Elements/progress), [BrowserStack cross-browser progress bar](https://www.browserstack.com/guide/how-to-create-cross-browser-compatible-html-progress-bar)

### Pitfall 2: Status Indicators Rely Only on Color

**What goes wrong:** Users with color vision deficiency cannot distinguish status states.

**Why it happens:** Default implementation uses red/green/yellow colors without additional visual cues.

**How to avoid:**
- Add shape differentiation: checkmark for complete, circle for pending, spinner for active
- Add text labels alongside badges
- Ensure 3:1 minimum contrast ratio for icons/shapes
- Test with grayscale filter to verify distinguishability

**Warning signs:**
- All status badges look identical in grayscale
- No icons, only color dots
- User confusion about job state

**Source:** [Carbon Design System Status Indicators](https://carbondesignsystem.com/patterns/status-indicator-pattern/), [WCAG 1.4.1 Use of Color](https://www.w3.org/WAI/WCAG21/Understanding/use-of-color.html)

### Pitfall 3: Step Tracker Connecting Line Misaligned

**What goes wrong:** The horizontal line connecting steps doesn't align with step icons, or breaks on mobile.

**Why it happens:** Using absolute positioning without accounting for responsive layout changes.

**How to avoid:**
- Position line relative to step list container, not individual items
- Use `left: 24px; right: 24px;` (half icon width) to start/end at icon centers
- Switch to vertical line on mobile with different absolute positioning
- Set `z-index` on items to ensure icons appear above line

**Warning signs:**
- Line extends beyond first/last step icons
- Line appears above icons instead of behind
- Line doesn't switch to vertical on mobile

### Pitfall 4: Ensemble Progress Hidden but DOM Present

**What goes wrong:** Empty conformer progress section visible on single-conformer jobs.

**Why it happens:** Using JavaScript/CSS to hide section instead of excluding from DOM.

**How to avoid:**
- Use Jinja2 `{% if job.conformer_mode == 'ensemble' %}` to conditionally render
- JavaScript should only update content within existing elements, not show/hide major sections
- Keep the existing conditional pattern from status.html

**Warning signs:**
- Empty card frames visible
- "Conformer Progress" header with no content
- Layout gaps where content should be excluded

### Pitfall 5: Elapsed Time Timer Accumulates Errors

**What goes wrong:** Displayed elapsed time drifts from actual elapsed time over long calculations.

**Why it happens:** Using `setInterval` with incremented seconds instead of calculating from start time.

**How to avoid:**
- Store `CREATED_AT` timestamp once
- On each tick, calculate: `Math.floor((now - CREATED_AT) / 1000)`
- Existing status.html already does this correctly - don't change the pattern

**Warning signs:**
- Elapsed time slightly off after 30+ minutes
- Timer "jumps" when tab is backgrounded

### Pitfall 6: 3D Viewer Initialization Before Container Exists

**What goes wrong:** 3Dmol.js fails to render, console shows "Unable to create viewer".

**Why it happens:** JavaScript runs before DOM is ready, or container has zero dimensions.

**How to avoid:**
- Initialize viewer in `DOMContentLoaded` event (existing pattern)
- Ensure container has explicit dimensions via CSS (min-height)
- Verify container ID matches between HTML and JavaScript

**Warning signs:**
- Empty white box where viewer should be
- Console error about viewer creation
- Viewer works on results page but not status page

## Code Examples

### Complete Status Page HTML Structure

```html
{% extends "base.html" %}

{% block title %}Job Status | QM NMR Calculator{% endblock %}

{% block page_css %}
<link rel="stylesheet" href="{{ url_for('static', path='/css/pages/status-page.css') }}">
{% endblock %}

{% block head %}
<script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.5.3/3Dmol-min.js"></script>
{% endblock %}

{% block content %}
<div class="page-wrapper">
    <header class="status-header">
        <h1>Job Status</h1>
        <span class="status-badge" id="status-badge" data-status="{{ job.status }}">
            {{ job.status }}
        </span>
    </header>

    <div class="bento-grid">
        <!-- Job Summary (span-3) -->
        <div class="bento-grid__item bento-grid__item--span-3">
            <article class="glass-card">
                <header class="glass-card__header">
                    <h2 class="glass-card__title">Job Details</h2>
                </header>
                <div class="glass-card__body">
                    <dl class="metadata-list">
                        <div class="metadata-list__item">
                            <dt>Job ID</dt>
                            <dd><code>{{ job.job_id }}</code></dd>
                        </div>
                        <div class="metadata-list__item">
                            <dt>Elapsed Time</dt>
                            <dd id="elapsed-time">--</dd>
                        </div>
                        {% if job.input_name %}
                        <div class="metadata-list__item">
                            <dt>Molecule</dt>
                            <dd>{{ job.input_name }}</dd>
                        </div>
                        {% endif %}
                        <div class="metadata-list__item">
                            <dt>SMILES</dt>
                            <dd><code class="smiles-code">{{ job.input_smiles }}</code></dd>
                        </div>
                        <div class="metadata-list__item">
                            <dt>Solvent</dt>
                            <dd>{{ job.solvent }}</dd>
                        </div>
                        <div class="metadata-list__item">
                            <dt>Preset</dt>
                            <dd>{{ job.preset }}</dd>
                        </div>
                    </dl>
                </div>
            </article>
        </div>

        <!-- 3D Molecule Preview (span-3) -->
        <div class="bento-grid__item bento-grid__item--span-3">
            <div class="viewer-card viewer-card--compact">
                <header class="viewer-card__header">
                    <h3 class="viewer-card__title">Molecular Structure</h3>
                </header>
                <div id="viewer-container-3d" class="viewer-card__canvas">
                    <!-- 3Dmol.js renders here -->
                </div>
                <footer class="viewer-card__footer">
                    <small class="muted">RDKit-generated geometry</small>
                </footer>
            </div>
        </div>

        <!-- Step Progress Tracker (span-6) -->
        <div class="bento-grid__item bento-grid__item--span-6" id="progress-section">
            <article class="glass-card">
                <header class="glass-card__header">
                    <h3 class="glass-card__title">Calculation Progress</h3>
                </header>
                <div class="glass-card__body">
                    <div class="step-tracker" id="step-tracker">
                        <!-- Populated by JavaScript -->
                    </div>
                </div>
            </article>
        </div>

        <!-- Conformer Progress (span-6) - ONLY FOR ENSEMBLE -->
        {% if job.conformer_mode == 'ensemble' %}
        <div class="bento-grid__item bento-grid__item--span-6">
            <article class="glass-card">
                <header class="glass-card__header">
                    <h3 class="glass-card__title">Conformer Progress</h3>
                </header>
                <div class="glass-card__body">
                    <div class="conformer-progress">
                        <div class="conformer-progress__header">
                            <span class="conformer-progress__stage" id="conformer-stage">--</span>
                            <span class="conformer-progress__count" id="conformer-count">--</span>
                        </div>
                        <progress
                            class="conformer-progress__bar"
                            id="conformer-progress-bar"
                            value="0"
                            max="100"
                            aria-label="Conformer calculation progress">
                        </progress>
                        <div class="conformer-progress__eta" id="eta-display"></div>
                    </div>

                    <details class="conformer-details">
                        <summary>Conformer Details</summary>
                        <table class="conformer-table" id="conformer-table">
                            <thead>
                                <tr>
                                    <th>ID</th>
                                    <th>Status</th>
                                    <th>Energy</th>
                                    <th>Population</th>
                                </tr>
                            </thead>
                            <tbody id="conformer-table-body">
                                <!-- Populated by JavaScript -->
                            </tbody>
                        </table>
                    </details>
                </div>
            </article>
        </div>
        {% endif %}

        <!-- Error Display (span-6) - Hidden by default -->
        <div class="bento-grid__item bento-grid__item--span-6" id="error-section" style="display: none;">
            <article class="glass-card glass-card--error">
                <header class="glass-card__header">
                    <h3 class="glass-card__title">
                        <span class="error-icon" aria-hidden="true">!</span>
                        Job Failed
                    </h3>
                </header>
                <div class="glass-card__body">
                    <p class="error-message" id="error-message"></p>
                    <div class="error-actions">
                        <a href="/" class="btn">Submit New Job</a>
                    </div>
                </div>
            </article>
        </div>
    </div>

    <footer class="status-footer">
        <small id="polling-status">Page refreshes automatically while job is running</small>
        <br>
        <a href="/">Submit another job</a>
    </footer>
</div>
{% endblock %}
```

### status-page.css Component CSS

```css
/*
 * Status Page Components
 *
 * Page-specific styling for the job status/progress page.
 * Uses design tokens from tokens.css and follows Phase 18 patterns.
 *
 * Key components:
 * - .step-tracker: Horizontal progress steps with state indicators
 * - .conformer-progress: Progress bar for ensemble job conformers
 * - .status-badge: Inline badge showing job status
 * - .glass-card--error: Error state card variant
 */

@layer components {

    /* ===== Status Header ===== */
    .status-header {
        display: flex;
        align-items: center;
        justify-content: space-between;
        gap: var(--space-md);
        margin-bottom: var(--space-xl);
        flex-wrap: wrap;
    }

    .status-header h1 {
        margin: 0;
        font-size: var(--text-2xl);
    }


    /* ===== Status Badge ===== */
    .status-badge {
        display: inline-flex;
        align-items: center;
        gap: var(--space-xs);
        padding: var(--space-xs) var(--space-sm);
        border-radius: var(--space-xs);
        font-size: var(--text-sm);
        font-weight: 500;
        text-transform: capitalize;
    }

    .status-badge::before {
        content: '';
        width: 8px;
        height: 8px;
        border-radius: 50%;
        background: currentColor;
    }

    .status-badge[data-status="queued"] {
        background: var(--color-gray-100);
        color: var(--color-gray-600);
    }

    .status-badge[data-status="running"] {
        background: hsl(195 85% 95%);
        color: var(--color-primary);
    }

    .status-badge[data-status="running"]::before {
        animation: pulse 1.5s ease-in-out infinite;
    }

    .status-badge[data-status="complete"] {
        background: hsl(142 76% 95%);
        color: var(--color-success);
    }

    .status-badge[data-status="failed"] {
        background: hsl(0 84% 95%);
        color: var(--color-error);
    }


    /* ===== Step Tracker ===== */
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
        position: relative;
        z-index: 1;
        flex: 1;
        text-align: center;
        min-width: 80px;
    }

    .step-tracker__icon {
        width: 24px;
        height: 24px;
        border-radius: 50%;
        display: flex;
        align-items: center;
        justify-content: center;
        font-size: 12px;
        background: var(--color-white);
        border: 2px solid var(--color-border);
        transition: background var(--transition-fast),
                    border-color var(--transition-fast);
    }

    /* Complete state */
    .step-tracker__item--complete .step-tracker__icon {
        background: var(--color-success);
        border-color: var(--color-success);
        color: var(--color-white);
    }

    .step-tracker__item--complete .step-tracker__icon::before {
        content: '\2713';
    }

    /* Active state */
    .step-tracker__item--active .step-tracker__icon {
        background: var(--color-primary);
        border-color: var(--color-primary);
        color: var(--color-white);
        animation: pulse 1.5s ease-in-out infinite;
    }

    .step-tracker__item--active .step-tracker__icon::before {
        content: '\2022';
        font-size: 16px;
    }

    /* Pending state */
    .step-tracker__item--pending .step-tracker__icon::before {
        content: '\25CB';
        color: var(--color-gray-400);
    }

    .step-tracker__label {
        font-size: var(--text-sm);
        font-weight: 500;
        color: var(--color-text);
        max-width: 100px;
    }

    .step-tracker__item--pending .step-tracker__label {
        color: var(--color-text-muted);
    }

    .step-tracker__duration {
        font-size: var(--text-xs);
        color: var(--color-text-muted);
        font-variant-numeric: tabular-nums;
    }


    /* ===== Conformer Progress ===== */
    .conformer-progress__header {
        display: flex;
        justify-content: space-between;
        align-items: baseline;
        margin-bottom: var(--space-sm);
    }

    .conformer-progress__stage {
        font-weight: 600;
        color: var(--color-text);
    }

    .conformer-progress__count {
        font-size: var(--text-sm);
        color: var(--color-text-muted);
        font-variant-numeric: tabular-nums;
    }

    /* Progress bar styling */
    .conformer-progress__bar {
        width: 100%;
        height: 1.5rem;
        appearance: none;
        -webkit-appearance: none;
        border: none;
        border-radius: var(--space-xs);
        background: var(--color-gray-100);
        overflow: hidden;
    }

    .conformer-progress__bar::-webkit-progress-bar {
        background: var(--color-gray-100);
        border-radius: var(--space-xs);
    }

    .conformer-progress__bar::-webkit-progress-value {
        background: var(--color-primary);
        border-radius: var(--space-xs);
        transition: width var(--transition-base);
    }

    .conformer-progress__bar::-moz-progress-bar {
        background: var(--color-primary);
        border-radius: var(--space-xs);
    }

    .conformer-progress__eta {
        margin-top: var(--space-sm);
        font-size: var(--text-sm);
        color: var(--color-text-muted);
    }

    .conformer-details {
        margin-top: var(--space-lg);
    }

    .conformer-details summary {
        cursor: pointer;
        font-weight: 500;
        color: var(--color-text);
    }


    /* ===== Error Card Variant ===== */
    .glass-card--error {
        background: hsl(0 84% 97%);
        border-color: hsl(0 84% 80%);
    }

    .glass-card--error .glass-card__header {
        border-bottom-color: hsl(0 84% 80%);
    }

    .glass-card--error .glass-card__title {
        color: var(--color-error);
        display: flex;
        align-items: center;
        gap: var(--space-sm);
    }

    .error-icon {
        width: 24px;
        height: 24px;
        border-radius: 50%;
        background: var(--color-error);
        color: var(--color-white);
        display: flex;
        align-items: center;
        justify-content: center;
        font-weight: bold;
        font-size: var(--text-sm);
    }

    .error-message {
        color: var(--color-text);
        line-height: 1.6;
        margin-bottom: var(--space-lg);
    }

    .error-actions {
        display: flex;
        gap: var(--space-md);
    }


    /* ===== Compact Viewer Card Variant ===== */
    .viewer-card--compact .viewer-card__canvas {
        min-height: 250px;
    }


    /* ===== Status Footer ===== */
    .status-footer {
        margin-top: var(--space-xl);
        text-align: center;
    }


    /* ===== Conformer Status Colors in Table ===== */
    .conformer-table tr[data-status="nmr_complete"] td:nth-child(2) {
        color: var(--color-success);
    }

    .conformer-table tr[data-status="failed"] td:nth-child(2) {
        color: var(--color-error);
    }

    .conformer-table tr[data-status="optimizing"] td:nth-child(2),
    .conformer-table tr[data-status="nmr_running"] td:nth-child(2) {
        color: var(--color-warning);
    }


    /* ===== Animation ===== */
    @keyframes pulse {
        0%, 100% { transform: scale(1); opacity: 1; }
        50% { transform: scale(1.1); opacity: 0.8; }
    }
}


/* ===== Mobile: Stack step tracker vertically ===== */
@media (max-width: 767px) {
    @layer components {
        .step-tracker__list {
            flex-direction: column;
            align-items: flex-start;
            gap: var(--space-md);
        }

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

        .step-tracker__label {
            max-width: none;
        }

        .status-header {
            flex-direction: column;
            align-items: flex-start;
        }
    }
}
```

### JavaScript Changes Required

The existing JavaScript in status.html will need updates to:

1. **Populate step tracker:** Build step list from `steps_completed` and `current_step`
2. **Update status badge:** Set `data-status` attribute for CSS styling
3. **Show/hide error section:** Toggle error card visibility on failure

```javascript
// Example: Building step tracker from API response
function updateStepTracker(data) {
    const tracker = document.getElementById('step-tracker');
    if (!tracker) return;

    // Define all possible steps in order
    const allSteps = [
        { id: 'conformer_generation', label: 'Conformer Generation' },
        { id: 'geometry_optimization', label: 'Geometry Optimization' },
        { id: 'nmr_calculation', label: 'NMR Calculation' },
        { id: 'post_processing', label: 'Post-processing' }
    ];

    // Build step list HTML
    const completedIds = data.steps_completed.map(s => s.step);

    let html = '<ol class="step-tracker__list">';

    allSteps.forEach(step => {
        let stateClass = 'step-tracker__item--pending';
        let duration = '';

        if (completedIds.includes(step.id)) {
            stateClass = 'step-tracker__item--complete';
            const completed = data.steps_completed.find(s => s.step === step.id);
            duration = formatDuration(completed.duration_seconds);
        } else if (data.current_step === step.id) {
            stateClass = 'step-tracker__item--active';
            duration = 'Running...';
        }

        html += `
            <li class="step-tracker__item ${stateClass}">
                <span class="step-tracker__icon" aria-hidden="true"></span>
                <span class="step-tracker__label">${step.label}</span>
                ${duration ? `<span class="step-tracker__duration">${duration}</span>` : ''}
            </li>
        `;
    });

    html += '</ol>';
    tracker.innerHTML = html;
}
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Custom div progress bars | Native `<progress>` element | HTML5 (2014) | Built-in accessibility, semantic meaning |
| Color-only status indicators | Color + shape + text | WCAG 2.1 (2018) | Accessible to color-blind users |
| Inline styles for status | CSS custom properties + data attributes | CSS3/modern browsers | Consistent theming, easier maintenance |
| JavaScript show/hide | Jinja2 conditional rendering | Best practice | Cleaner DOM, no flash of unstyled content |

**Deprecated/outdated:**
- Pico CSS inline styles in current status.html - replace with glass-card components
- Color-only status differentiation - add icons/shapes
- Manual class toggling for status - use data attributes with CSS attribute selectors

## Open Questions

### 1. Step Order for Ensemble vs Single Jobs

**What we know:**
- Single conformer: Geometry Optimization -> NMR Calculation -> Post-processing
- Ensemble: Conformer Generation -> (per-conformer optimization) -> (per-conformer NMR) -> Averaging

**What's unclear:**
- Should step tracker show individual conformer steps or aggregate?
- How to represent "3/12 conformers optimized" in step tracker format?

**Recommendation:**
- Keep step tracker simple: show phase-level steps (Conformer Gen, Optimization, NMR, Averaging)
- Use conformer progress bar below for detailed X/N tracking
- Mark step as "active" while any conformers are processing that phase

### 2. ETA Accuracy and Display

**What we know:**
- API returns `eta_seconds` for running jobs
- ETA is calculated from average conformer processing time

**What's unclear:**
- How accurate is ETA for first few conformers (cold start)?
- Should ETA be hidden until enough samples for accuracy?

**Recommendation:**
- Display ETA only after 2+ conformers complete
- Show "Calculating estimate..." initially
- Round ETA to nearest minute for display

### 3. 3D Viewer Update During Calculation

**What we know:**
- Geometry is available early (RDKit-generated before DFT optimization)
- Final optimized geometry differs from initial

**What's unclear:**
- Should viewer update to show DFT-optimized geometry as conformers complete?
- Performance impact of updating viewer during polling?

**Recommendation:**
- Show initial RDKit geometry (current behavior)
- Don't update during calculation (avoid user confusion)
- Note "RDKit-generated geometry" in footer
- Final optimized geometry shown on results page

## Sources

### Primary (HIGH confidence)

- [MDN progress element](https://developer.mozilla.org/en-US/docs/Web/HTML/Reference/Elements/progress) - Native element specification, pseudo-elements, accessibility
- [MDN ARIA progressbar role](https://developer.mozilla.org/en-US/docs/Web/Accessibility/ARIA/Reference/Roles/progressbar_role) - Accessibility requirements
- [Carbon Design System Status Indicators](https://carbondesignsystem.com/patterns/status-indicator-pattern/) - Color, shape, accessibility patterns
- Existing codebase: Phase 18-20 CSS files, status.html, schemas.py - Established patterns

### Secondary (MEDIUM confidence)

- [CSS Progress Wizard](https://christabor.github.io/css-progress-wizard/) - Step tracker flexbox pattern (not maintained but pattern valid)
- [Progress Tracker by Nigel O'Toole](https://nigelotoole.github.io/progress-tracker/) - Step tracker implementation reference
- [BrowserStack cross-browser progress bar](https://www.browserstack.com/guide/how-to-create-cross-browser-compatible-html-progress-bar) - Browser compatibility testing
- [Custom CSS Progress Bar](https://nikitahl.com/progress-bar-css) - Cross-browser styling techniques

### Tertiary (LOW confidence - requires validation)

- [Loading.io CSS spinners](https://loading.io/css/) - Animation patterns (may need performance testing)
- WebSearch results on step wizard libraries (most are unmaintained, patterns extracted)

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - Using established Phase 18-20 components
- Step tracker pattern: HIGH - Flexbox layout, CSS-only icons well-documented
- Progress bar styling: HIGH - MDN documentation comprehensive, cross-browser tested
- Status indicators: HIGH - Carbon Design System provides authoritative guidance
- JavaScript integration: HIGH - Existing status.html code works, minimal changes needed
- Mobile responsiveness: MEDIUM - Patterns exist but need testing with actual content

**Research date:** 2026-01-29
**Valid until:** 2026-04-29 (90 days - CSS patterns stable)
