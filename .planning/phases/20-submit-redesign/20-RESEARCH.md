# Phase 20: Submit Page Redesign - Research

**Researched:** 2026-01-29
**Domain:** HTML Form Design, CSS Form Styling, Molecule Preview, Accessibility
**Confidence:** HIGH

## Summary

Phase 20 redesigns the submit page to improve form usability through logical grouping, solid backgrounds for form inputs, molecule preview feedback, clear required/optional field indicators, and well-described conformer options. The research confirms that the existing CSS infrastructure from Phase 18 provides all necessary foundation components.

The primary challenge is adding real-time molecule preview from SMILES input. SmilesDrawer is the established client-side library for this purpose - it is dependency-free, MIT licensed, and renders SMILES strings to HTML canvas elements with minimal code. The library is available via CDN with well-documented error handling for invalid SMILES.

Form accessibility follows established patterns: fieldsets with legends for logical grouping, solid white backgrounds for all inputs (already implemented in base.css), and required field indication using both visual asterisks AND the `required` attribute plus `aria-required="true"` for screen reader compatibility.

**Primary recommendation:** Create a two-column form layout using CSS Grid with logical fieldset groups, add a SmilesDrawer-powered preview card, and use existing design tokens and glass-card components from Phase 18.

## Standard Stack

The established libraries/tools for this domain:

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| SmilesDrawer | 1.2.0 | Client-side SMILES to 2D structure rendering | Dependency-free, MIT license, fast canvas rendering, built-in error handling |
| CSS Grid | Native | Form layout in two columns | Built into Phase 18 layout.css |
| Fieldset/Legend | HTML5 | Form grouping | Native accessibility, screen reader support |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| Design Tokens | Phase 18 | Consistent spacing, colors, borders | All styling |
| glass-card.css | Phase 18 | Container styling | Molecule preview card |
| base.css form styles | Phase 18 | Input, select, button styling | All form elements |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| SmilesDrawer | RDKit.js | RDKit.js is larger (WASM), more powerful but overkill for preview |
| SmilesDrawer | 3Dmol.js | Already using 3Dmol.js for 3D; 2D preview better for form UX |
| CSS Grid form | Flexbox | Grid provides cleaner two-column alignment |

**CDN URL:**
```html
<script src="https://unpkg.com/smiles-drawer@1.2.0/dist/smiles-drawer.min.js"></script>
```

## Architecture Patterns

### Recommended Form Layout Structure
```
.submit-page
  .page-wrapper
    .page-header
    .submit-layout (CSS Grid: 2 columns on desktop)
      .submit-form (left column)
        fieldset.form-group.form-group--molecule
          legend
          .form-group__field (SMILES input)
          .form-group__field (file upload)
        fieldset.form-group.form-group--settings
          legend
          .form-group__field (solvent)
          .form-group__field (preset)
        fieldset.form-group.form-group--conformer
          legend
          .form-group__field (mode)
          .form-group__field (method)
          details.form-group__advanced
            (max conformers)
        fieldset.form-group.form-group--optional
          legend
          .form-group__field (name)
          .form-group__field (email)
        button[type=submit]
      .molecule-preview (right column, sticky)
        .preview-card (solid background, not glass)
          canvas#molecule-preview
          .preview-card__status
```

### Pattern 1: Form Group Component
**What:** Reusable fieldset styling with BEM naming
**When to use:** All form sections in submit page
**Example:**
```css
/* Source: BEM naming from Phase 18 patterns */
@layer components {
    .form-group {
        border: 1px solid var(--color-border);
        border-radius: 0.5rem;
        padding: var(--space-lg);
        margin-bottom: var(--space-lg);
        background: var(--color-white);
    }

    .form-group legend {
        font-weight: 600;
        font-size: var(--text-lg);
        padding: 0 var(--space-sm);
        color: var(--color-text);
    }

    .form-group__field {
        margin-bottom: var(--space-md);
    }

    .form-group__field:last-child {
        margin-bottom: 0;
    }
}
```

### Pattern 2: Required Field Indicator
**What:** Visual asterisk + programmatic required attribute
**When to use:** All required fields (solvent, conformer_mode)
**Example:**
```html
<!-- Source: WCAG/WAI best practices -->
<label for="solvent">
    Solvent
    <span class="required-indicator" aria-hidden="true">*</span>
</label>
<select id="solvent" name="solvent" required aria-required="true">
```

```css
.required-indicator {
    color: var(--color-error);
    margin-left: var(--space-xs);
}

/* Form instructions above form */
.form-instructions {
    font-size: var(--text-sm);
    color: var(--color-text-muted);
    margin-bottom: var(--space-lg);
}
```

### Pattern 3: Molecule Preview Card
**What:** Solid-background card with canvas for SmilesDrawer
**When to use:** Right column of submit form
**Example:**
```css
/* Solid background like viewer-card from results-page.css */
.preview-card {
    background: var(--color-white);
    border: 1px solid var(--color-border);
    border-radius: var(--glass-radius);
    box-shadow: var(--glass-shadow);
    padding: var(--space-lg);
    position: sticky;
    top: var(--space-lg);
}

.preview-card__canvas {
    width: 100%;
    aspect-ratio: 1;
    background: var(--color-gray-50);
    border: 1px solid var(--color-gray-200);
    border-radius: var(--space-xs);
}

.preview-card__status {
    margin-top: var(--space-md);
    font-size: var(--text-sm);
    text-align: center;
}

.preview-card__status--valid {
    color: var(--color-success);
}

.preview-card__status--invalid {
    color: var(--color-error);
}

.preview-card__status--empty {
    color: var(--color-text-muted);
}
```

### Pattern 4: Two-Column Submit Layout
**What:** CSS Grid responsive layout
**When to use:** Submit page wrapper
**Example:**
```css
.submit-layout {
    display: grid;
    grid-template-columns: 1fr 350px;
    gap: var(--space-xl);
    align-items: start;
}

@media (max-width: 900px) {
    .submit-layout {
        grid-template-columns: 1fr;
    }

    .molecule-preview {
        order: -1; /* Preview first on mobile */
    }
}
```

### Anti-Patterns to Avoid
- **Glass effect on form inputs:** Reduces contrast, fails WCAG. Use solid white backgrounds (already in base.css).
- **Placeholder as label:** Disappears on focus, inaccessible. Use visible labels above inputs.
- **Asterisk without instructions:** Explain "* indicates required" at form top.
- **Nested fieldsets:** Causes erratic screen reader behavior. Use flat fieldset structure.
- **CSS reordering of form fields:** Breaks tab order. Keep visual and DOM order aligned.

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| SMILES validation | Custom regex | SmilesDrawer.parse() error callback | SMILES syntax is complex, semantic validation needed |
| 2D molecule rendering | Canvas drawing | SmilesDrawer | Handles stereochemistry, aromatic rings, layout |
| Form input styling | Custom styles | base.css existing styles | Already WCAG compliant, consistent |
| Required field indication | Just asterisks | asterisk + required + aria-required | Accessibility requires programmatic indication |
| Focus ring styling | Custom focus | base.css :focus-visible | Already implemented with proper offset |

**Key insight:** SmilesDrawer handles the complex chemistry parsing and 2D coordinate generation. Attempting to validate SMILES with regex will miss semantic errors and produce incorrect previews.

## Common Pitfalls

### Pitfall 1: Glass Background on Inputs
**What goes wrong:** Text contrast fails WCAG 4.5:1 ratio when input background is semi-transparent
**Why it happens:** Visual consistency with glassmorphism theme
**How to avoid:** Phase 18 already solved this - base.css sets solid white background on all inputs
**Warning signs:** Text hard to read on colorful backgrounds

### Pitfall 2: SmilesDrawer Canvas Sizing
**What goes wrong:** Canvas renders blurry or at wrong dimensions
**Why it happens:** Canvas element size vs CSS size mismatch, devicePixelRatio not handled
**How to avoid:** Set canvas width/height attributes, not just CSS dimensions
**Warning signs:** Fuzzy molecule rendering, especially on retina displays

```javascript
// Correct canvas sizing
const canvas = document.getElementById('molecule-preview');
const dpr = window.devicePixelRatio || 1;
const rect = canvas.getBoundingClientRect();
canvas.width = rect.width * dpr;
canvas.height = rect.height * dpr;
canvas.getContext('2d').scale(dpr, dpr);
```

### Pitfall 3: File Input Mutual Exclusion
**What goes wrong:** User submits both SMILES and file, or neither
**Why it happens:** Two input methods without clear mutual exclusivity
**How to avoid:** Existing JavaScript in submit.html already handles this - disable other input when one has value
**Warning signs:** Confusing which input will be used

### Pitfall 4: Debounce Missing on SMILES Input
**What goes wrong:** SmilesDrawer.parse() called on every keystroke, laggy UX
**Why it happens:** Eager validation without rate limiting
**How to avoid:** Debounce input handler by 300-500ms
**Warning signs:** Slow typing, canvas flickering

```javascript
// Debounce pattern
let debounceTimer;
smilesInput.addEventListener('input', () => {
    clearTimeout(debounceTimer);
    debounceTimer = setTimeout(() => {
        updateMoleculePreview(smilesInput.value);
    }, 400);
});
```

### Pitfall 5: Conformer Method Description Hidden
**What goes wrong:** Users don't understand RDKit vs CREST difference
**Why it happens:** Just showing option names without explanation
**How to avoid:** Add descriptive help text below select, use small element styling
**Warning signs:** Support questions about which method to choose

## Code Examples

Verified patterns from official sources:

### SmilesDrawer Initialization and Usage
```javascript
// Source: https://github.com/reymond-group/smilesDrawer
const smilesDrawer = new SmilesDrawer.Drawer({
    width: 350,
    height: 350,
    bondThickness: 0.6,
    bondLength: 15,
    terminalCarbons: true,
    explicitHydrogens: false,
    themes: {
        light: {
            C: '#222',
            O: '#e74c3c',
            N: '#3498db',
            S: '#f1c40f',
            P: '#e67e22',
            F: '#1abc9c',
            CL: '#1abc9c',
            BR: '#e74c3c',
            I: '#9b59b6',
            H: '#aaa'
        }
    }
});

function updateMoleculePreview(smiles) {
    const canvas = document.getElementById('molecule-preview');
    const status = document.getElementById('preview-status');

    if (!smiles || smiles.trim() === '') {
        // Clear canvas, show empty state
        const ctx = canvas.getContext('2d');
        ctx.clearRect(0, 0, canvas.width, canvas.height);
        status.textContent = 'Enter a SMILES string to preview';
        status.className = 'preview-card__status preview-card__status--empty';
        return;
    }

    SmilesDrawer.parse(smiles, function(tree) {
        // Valid SMILES - draw molecule
        smilesDrawer.draw(tree, 'molecule-preview', 'light', false);
        status.textContent = 'Valid SMILES';
        status.className = 'preview-card__status preview-card__status--valid';
    }, function(err) {
        // Invalid SMILES - show error
        const ctx = canvas.getContext('2d');
        ctx.clearRect(0, 0, canvas.width, canvas.height);
        status.textContent = 'Invalid SMILES syntax';
        status.className = 'preview-card__status preview-card__status--invalid';
    });
}
```

### Required Field with Accessibility
```html
<!-- Source: W3C WAI ARIA2 technique -->
<p class="form-instructions">
    Fields marked with <span aria-hidden="true">*</span>
    <span class="sr-only">asterisk</span> are required.
</p>

<div class="form-group__field">
    <label for="solvent">
        Solvent
        <span class="required-indicator" aria-hidden="true">*</span>
    </label>
    <select id="solvent" name="solvent" required aria-required="true">
        <option value="" disabled selected>Select a solvent...</option>
        {% for solvent in solvents %}
        <option value="{{ solvent.value }}">{{ solvent.label }}</option>
        {% endfor %}
    </select>
</div>
```

### Fieldset with Descriptive Legend
```html
<!-- Source: W3C WAI form grouping tutorial -->
<fieldset class="form-group form-group--conformer">
    <legend>Conformer Settings</legend>

    <div class="form-group__field">
        <label for="conformer_mode">
            Calculation Mode
            <span class="required-indicator" aria-hidden="true">*</span>
        </label>
        <select id="conformer_mode" name="conformer_mode" required aria-required="true">
            <option value="ensemble" selected>Ensemble (recommended)</option>
            <option value="single">Single conformer</option>
        </select>
        <small class="field-help">
            Ensemble mode generates multiple conformers for more accurate predictions.
        </small>
    </div>

    <div class="form-group__field">
        <label for="conformer_method">Conformer Method</label>
        <select id="conformer_method" name="conformer_method">
            <option value="rdkit_kdg" selected>RDKit KDG (fast)</option>
            <option value="crest" {% if not crest_available %}disabled{% endif %}>
                CREST/xTB (thorough){% if not crest_available %} - not installed{% endif %}
            </option>
        </select>
        <small class="field-help">
            RDKit is fast and always available. CREST provides better sampling for complex molecules.
        </small>
    </div>
</fieldset>
```

### Field Help Text Styling
```css
/* Source: Phase 18 design tokens + WAI patterns */
.field-help {
    display: block;
    margin-top: var(--space-xs);
    font-size: var(--text-sm);
    color: var(--color-text-muted);
    line-height: 1.4;
}
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Pico CSS form defaults | Custom base.css form styles | Phase 18 | Full control, WCAG compliant |
| Single-column form | Two-column with preview | Phase 20 | Better visual feedback |
| No molecule preview | SmilesDrawer real-time preview | Phase 20 | User knows input is valid |
| Asterisk-only required | asterisk + required + aria-required | Phase 20 | Screen reader accessible |

**Deprecated/outdated:**
- Placeholder text as labels: Removed, using visible labels
- Pico CSS fieldset styling: Replaced by custom form-group component

## Open Questions

Things that couldn't be fully resolved:

1. **File upload preview**
   - What we know: SMILES preview via SmilesDrawer works well
   - What's unclear: Whether to show preview for uploaded MOL/SDF files (requires parsing)
   - Recommendation: Start with SMILES preview only. MOL/SDF preview would require additional parsing library (potentially RDKit.js or server-side)

2. **Canvas retina handling**
   - What we know: devicePixelRatio scaling needed for crisp rendering
   - What's unclear: SmilesDrawer's internal handling of high-DPI
   - Recommendation: Test on retina display, apply manual scaling if needed

## Sources

### Primary (HIGH confidence)
- Phase 18 CSS infrastructure (layers.css, tokens.css, base.css, layout.css, components/glass-card.css) - direct codebase review
- [SmilesDrawer GitHub](https://github.com/reymond-group/smilesDrawer) - CDN URL, API documentation
- [SmilesDrawer Documentation](https://smilesdrawer.readthedocs.io/en/latest/) - Usage patterns, options
- [W3C WAI ARIA2 Technique](https://www.w3.org/WAI/WCAG21/Techniques/aria/ARIA2) - Required field accessibility
- [W3C WAI Form Grouping Tutorial](https://www.w3.org/WAI/tutorials/forms/grouping/) - Fieldset/legend patterns

### Secondary (MEDIUM confidence)
- [WebAIM Form Controls](https://webaim.org/techniques/forms/controls) - Required field best practices
- [Deque Required Form Fields](https://www.deque.com/blog/anatomy-of-accessible-forms-required-form-fields/) - Asterisk + aria-required combination
- [TetraLogical Fieldset/Legend](https://tetralogical.com/blog/2025/01/31/foundations-fieldset-and-legend/) - Screen reader behavior

### Tertiary (LOW confidence)
- [Glassmorphism Accessibility](https://www.newtarget.com/web-insights-blog/glassmorphism/) - Solid background recommendation for inputs (validates Phase 18 approach)

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - SmilesDrawer is well-documented, Phase 18 CSS verified in codebase
- Architecture: HIGH - Patterns derived from existing results-page.css and WAI guidelines
- Pitfalls: MEDIUM - Canvas sizing and debounce based on general web dev knowledge

**Research date:** 2026-01-29
**Valid until:** 2026-03-29 (60 days - stable domain, no expected changes)
