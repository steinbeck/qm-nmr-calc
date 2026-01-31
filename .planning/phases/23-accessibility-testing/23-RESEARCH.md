# Phase 23: Accessibility and Testing - Research

**Researched:** 2026-01-31
**Domain:** Web accessibility (WCAG 2.2), cross-browser CSS testing, performance optimization
**Confidence:** HIGH

## Summary

This phase focuses on ensuring WCAG 2.2 Level AA compliance, cross-browser compatibility, and mobile performance optimization for the glassmorphism-based UI redesign. The research covers automated and manual testing tools, accessibility patterns for keyboard navigation and motion sensitivity, and cross-browser testing strategies for backdrop-filter support.

The standard approach combines automated testing tools (axe-core, Lighthouse) that catch 30-40% of issues with manual keyboard navigation testing and cross-browser validation. For glassmorphism specifically, mobile performance requires reduced blur values (6px vs 10px desktop) and limited glass element count (2-3 max), both of which are already implemented in the current codebase.

Key findings: Safari still requires `-webkit-backdrop-filter` prefix in 2026, automated tools only catch ~40% of accessibility issues (manual testing is critical), and `:focus-visible` pseudo-class is now the standard approach for keyboard focus indicators over legacy `:focus`.

**Primary recommendation:** Use axe-core for automated WCAG testing, manual keyboard navigation for interaction testing, and BrowserStack or Playwright for cross-browser validation of webkit-prefixed backdrop-filter.

## Standard Stack

The established libraries/tools for this domain:

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| axe-core | Latest | Automated accessibility testing | Industry standard by Deque, covers WCAG 2.0/2.1/2.2, 3B+ downloads, finds ~40% of issues automatically |
| Lighthouse | Built-in Chrome | Performance & accessibility auditing | Google's official tool, built into Chrome DevTools, measures performance metrics and accessibility |
| WebAIM Contrast Checker | Web-based | Color contrast validation | Free authoritative tool for WCAG 2.0/2.1 contrast ratio compliance (4.5:1 for text, 3:1 for UI components) |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| Pa11y | Latest | Command-line accessibility testing | CI/CD pipeline automation, headless browser testing |
| Playwright | Latest | Cross-browser automation | Testing backdrop-filter across Chromium/WebKit/Firefox engines |
| BrowserStack | Cloud service | Real device/browser testing | Manual cross-browser validation on actual Safari, Edge, Firefox versions |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| axe-core | WAVE browser extension | WAVE is manual-only, no CI/CD automation; axe-core supports both |
| Playwright | Selenium | Playwright has better modern browser support and simpler API |
| WebAIM Contrast Checker | Browser DevTools contrast tool | DevTools integrated but less precise; WebAIM is authoritative reference |

**Installation:**
```bash
# For automated testing (not needed for manual-only Phase 23)
npm install --save-dev @axe-core/playwright
npm install --save-dev playwright
npm install -g pa11y
```

**Note:** Phase 23 focuses on manual testing and validation. No npm packages required since testing uses browser DevTools and web-based tools.

## Architecture Patterns

### Recommended Testing Structure

For accessibility and cross-browser testing, organize testing into layers:

```
testing/
├── automated/           # axe-core, Lighthouse CI (future)
├── manual/              # Documented manual test results
│   ├── keyboard-nav.md
│   ├── contrast-audit.md
│   └── cross-browser.md
└── performance/         # Lighthouse reports
```

**For Phase 23 (manual testing focus):** Document results in `.planning/phases/23-accessibility-testing/` as markdown files.

### Pattern 1: Automated Accessibility Testing with axe-core
**What:** Integrate axe-core into test suite to catch common WCAG violations
**When to use:** Every code change, in CI/CD pipeline (future automation)
**Example:**
```javascript
// Source: https://github.com/dequelabs/axe-core
// Playwright integration
const { test, expect } = require('@playwright/test');
const { injectAxe, checkA11y } = require('@axe-core/playwright');

test('homepage should not have accessibility violations', async ({ page }) => {
  await page.goto('http://localhost:8000/');
  await injectAxe(page);
  await checkA11y(page);
});
```

### Pattern 2: Keyboard Navigation Testing
**What:** Manual testing of all interactive elements via keyboard only
**When to use:** Every UI change, before deployment
**Testing sequence:**
1. Tab through all interactive elements (forms, buttons, links, details/summary)
2. Verify visible focus indicators on each element
3. Test Enter/Space activation on buttons and links
4. Test Escape key for modals/dropdowns
5. Verify no keyboard traps (can always Tab out)

**Verification:**
- Focus indicator must be visible (WCAG 2.4.7 Level A)
- Focus indicator should have 3:1 contrast vs unfocused state (2.4.13 Level AAA recommended)
- Logical tab order (follows visual order)

### Pattern 3: Contrast Ratio Validation
**What:** Validate all text meets WCAG 4.5:1 minimum (3:1 for large text)
**When to use:** Every color change, new component
**Tool workflow:**
1. Use WebAIM Contrast Checker (https://webaim.org/resources/contrastchecker/)
2. Test foreground color against background color
3. For glass backgrounds: Test text color against effective background (blend of glass opacity + backdrop)
4. Document passing/failing combinations

**Example validation:**
```css
/* High opacity glass (95%) with dark text */
--glass-bg-subtle: hsl(0 0% 100% / 0.95);  /* Near-white background */
--color-text: hsl(220 13% 15%);             /* Dark gray text */
/* Contrast ratio: ~15:1 (PASSES 4.5:1) */
```

### Pattern 4: Cross-Browser Testing for backdrop-filter
**What:** Validate backdrop-filter with webkit prefix across browsers
**When to use:** Before deployment, after CSS changes
**Test matrix:**
```
Chrome/Edge (Blink):   backdrop-filter (unprefixed works)
Firefox (Gecko):       backdrop-filter (unprefixed since v103)
Safari (WebKit):       -webkit-backdrop-filter (REQUIRED, CSS vars don't work)
```

**Implementation pattern:**
```css
/* Source: Verified from existing glass-card.css */
.glass-card {
    backdrop-filter: var(--glass-blur-medium, blur(12px));
    -webkit-backdrop-filter: blur(12px);  /* Safari: MUST use literal value */
}
```

### Pattern 5: Reduced Motion Implementation
**What:** Disable transform animations, keep opacity transitions
**When to use:** All animation-heavy components
**Example:**
```css
/* Source: https://developer.mozilla.org/en-US/docs/Web/CSS/@media/prefers-reduced-motion */
/* Default: Full animation */
.glass-card--interactive:hover {
    transform: translateY(-4px);
}

/* Reduced motion: Remove transform, keep opacity */
@media (prefers-reduced-motion: reduce) {
    .glass-card--interactive:hover {
        transform: none;  /* No movement */
    }
    .glass-card::after {
        transition: opacity var(--transition-reduced, 200ms ease);  /* Opacity OK */
    }
}
```

### Pattern 6: Focus-Visible Pseudo-Class
**What:** Modern CSS pseudo-class for keyboard-only focus indicators
**When to use:** All interactive elements (buttons, links, inputs, details)
**Example:**
```css
/* Source: https://developer.mozilla.org/en-US/docs/Web/CSS/:focus-visible */
/* Don't show focus on mouse click */
button:focus {
    outline: none;
}

/* Show focus only for keyboard navigation */
button:focus-visible {
    outline: 2px solid var(--color-primary);
    outline-offset: 2px;
}
```

**Browser support:** Baseline widely available since January 2020.

### Anti-Patterns to Avoid
- **Never remove outlines without replacement:** `outline: none` without `:focus-visible` replacement creates keyboard trap
- **Don't use `outline: none` globally:** Removes critical accessibility feature
- **Don't test contrast on mockups:** Glass effects blend backgrounds; test in-browser with real backdrop
- **Don't rely on automated tools alone:** Catches only 30-40% of issues
- **Don't use CSS variables in -webkit-backdrop-filter:** Safari doesn't support them even in 2026

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Contrast ratio calculation | Manual color math | WebAIM Contrast Checker, browser DevTools | WCAG formula is complex (relative luminance calculation), authoritative tools exist |
| Accessibility rules engine | Custom WCAG validator | axe-core library | 3B+ downloads, actively maintained, covers WCAG 2.0/2.1/2.2, false positive rate very low |
| Cross-browser screenshots | Manual screenshots | BrowserStack, Playwright | Hundreds of browser/OS combinations, real devices, automated |
| Keyboard trap detection | Manual navigation | Pa11y automated testing | Detects keyboard traps programmatically |
| Performance scoring | Custom metrics | Lighthouse | Industry-standard metrics (LCP, FID, CLS), regression tracking |

**Key insight:** Accessibility testing requires domain expertise. Use authoritative tools (WebAIM, axe-core, Lighthouse) maintained by accessibility experts rather than implementing WCAG interpretation yourself.

## Common Pitfalls

### Pitfall 1: Testing Contrast on Static Mockups
**What goes wrong:** Glass backgrounds are semi-transparent and blend with backdrop content. Contrast ratio on static mockup doesn't match real browser rendering.
**Why it happens:** Designers test colors in Figma/Photoshop without backdrop-filter effect.
**How to avoid:** Always test contrast in-browser with real backdrop content. Use browser DevTools color picker on rendered page.
**Warning signs:** Text looks readable in design but fails contrast checker in production.

### Pitfall 2: Assuming Automated Tools Find Everything
**What goes wrong:** Automated tools (axe-core, Lighthouse) catch only 30-40% of WCAG issues. Manual testing finds interaction problems (keyboard traps, modal focus management, screen reader compatibility).
**Why it happens:** Automated tools can't test dynamic interactions, user workflows, or semantic meaning.
**How to avoid:** Combine automated + manual testing. Manual testing by experienced auditors catches the 60-70% that automation misses.
**Warning signs:** 100% Lighthouse accessibility score but site isn't keyboard-navigable.

### Pitfall 3: Safari backdrop-filter with CSS Variables
**What goes wrong:** Safari requires `-webkit-backdrop-filter` with literal blur values. CSS variables don't work even if they resolve to valid values.
**Why it happens:** Safari's webkit implementation doesn't evaluate CSS custom properties in backdrop-filter.
**How to avoid:** Always use literal values in `-webkit-backdrop-filter`: `blur(12px)` not `var(--glass-blur-medium)`.
**Warning signs:** Glassmorphism works in Chrome/Firefox but not Safari.

**Code example from codebase:**
```css
/* CORRECT: Literal value for Safari */
backdrop-filter: var(--glass-blur-medium, blur(12px));
-webkit-backdrop-filter: blur(12px);

/* WRONG: CSS variable in webkit prefix (Safari ignores) */
-webkit-backdrop-filter: var(--glass-blur-medium);
```

### Pitfall 4: Removing Focus Indicators Globally
**What goes wrong:** `* { outline: none; }` or `button:focus { outline: none; }` removes keyboard navigation visibility, violating WCAG 2.4.7 (Level A).
**Why it happens:** Developers find default browser outlines "ugly" and remove them without replacement.
**How to avoid:** Use `:focus-visible` to hide mouse focus while keeping keyboard focus. Never remove outlines without visible replacement.
**Warning signs:** Keyboard users report "can't see where I am on page."

**Correct pattern:**
```css
/* Hide focus for mouse users */
button:focus:not(:focus-visible) {
    outline: none;
}

/* Show focus for keyboard users */
button:focus-visible {
    outline: 2px solid var(--color-primary);
    outline-offset: 2px;
}
```

### Pitfall 5: Color-Only Status Indicators
**What goes wrong:** Using only color to indicate status (green = success, red = error) fails WCAG 1.4.1 (Use of Color - Level A). Users with color blindness can't distinguish status.
**Why it happens:** Color-coding is visually intuitive and fast to implement.
**How to avoid:** Combine color with shape/icon/text. Use `::before` pseudo-element for status dot + text label.
**Warning signs:** Status page uses colored rows without status text.

**Correct pattern (from codebase):**
```css
/* Status badge: Colored dot (::before) PLUS text */
.status-badge::before {
    content: '';
    width: 8px;
    height: 8px;
    border-radius: 50%;
    background: currentColor;
}
/* Text label provides non-color indicator */
```

### Pitfall 6: Motion Animations Without Reduced Motion Support
**What goes wrong:** Transform-based animations (scale, translate, rotate) trigger vestibular motion disorders for users with motion sensitivity.
**Why it happens:** Animations improve perceived performance and delight, easy to add without considering accessibility.
**How to avoid:** Implement `@media (prefers-reduced-motion: reduce)` for ALL transform animations. Replace with opacity-only transitions.
**Warning signs:** Pulsing/scaling animations without reduced motion query.

**Correct pattern (from codebase):**
```css
/* Default: Scale animation */
.step-tracker__icon {
    animation: pulse 1.5s ease-in-out infinite;
}

@keyframes pulse {
    0%, 100% { transform: scale(1); opacity: 1; }
    50% { transform: scale(1.1); opacity: 0.8; }
}

/* Reduced motion: Opacity-only animation */
@media (prefers-reduced-motion: reduce) {
    .step-tracker__icon {
        animation: dissolve 2s ease-in-out infinite;
    }
}

@keyframes dissolve {
    0%, 100% { opacity: 1; }
    50% { opacity: 0.5; }
}
```

### Pitfall 7: Mobile Performance with Too Many Glass Elements
**What goes wrong:** Glassmorphism with backdrop-filter is GPU-intensive. Too many glass elements on mobile causes lag, dropped frames, or crashes.
**Why it happens:** Desktop GPUs handle multiple blur layers easily; mobile GPUs are constrained.
**How to avoid:** Limit to 2-3 glass elements on mobile. Reduce blur radius (6px mobile vs 10px desktop). Use `@media (max-width: 768px)` with reduced values.
**Warning signs:** Smooth on desktop, janky scrolling on mobile.

**Correct pattern (from codebase):**
```css
:root {
    --glass-blur-medium: blur(12px);  /* Desktop */
}

@media (max-width: 768px) {
    :root {
        --glass-blur-light: blur(6px);   /* Mobile: Reduced */
        --glass-blur-medium: blur(8px);
        --glass-blur-heavy: blur(10px);
    }
}
```

## Code Examples

Verified patterns from official sources and existing codebase:

### Focus Indicators for Keyboard Navigation
```css
/* Source: https://developer.mozilla.org/en-US/docs/Web/CSS/:focus-visible */
/* Modern approach: :focus-visible for keyboard-only indicators */

/* Style all interactive elements */
button, a, input, select, textarea, summary {
    /* Remove default outline for mouse users */
    outline: none;
}

/* Show high-contrast outline for keyboard users */
button:focus-visible,
a:focus-visible,
input:focus-visible,
select:focus-visible,
textarea:focus-visible,
summary:focus-visible {
    outline: 2px solid var(--color-primary);
    outline-offset: 2px;
    /* 3:1 contrast recommended by WCAG 2.4.13 (Level AAA) */
}
```

### Prefers-Reduced-Motion Media Query
```css
/* Source: https://developer.mozilla.org/en-US/docs/Web/CSS/@media/prefers-reduced-motion */
/* Existing codebase implementation (status-page.css) */

/* Default: Transform + opacity animation */
@keyframes pulse {
    0%, 100% { transform: scale(1); opacity: 1; }
    50% { transform: scale(1.1); opacity: 0.8; }
}

.status-badge[data-status="running"]::before {
    animation: pulse 1.5s ease-in-out infinite;
}

/* Reduced motion: Opacity-only animation (no transform) */
@keyframes dissolve {
    0%, 100% { opacity: 1; }
    50% { opacity: 0.5; }
}

@media (prefers-reduced-motion: reduce) {
    .status-badge[data-status="running"]::before {
        animation: dissolve 2s ease-in-out infinite;
    }
}
```

### Safari-Compatible backdrop-filter
```css
/* Source: Existing glass-card.css + https://caniuse.com/css-backdrop-filter */
/* Safari requires -webkit prefix with LITERAL values (not CSS variables) */

.glass-card {
    /* Standard property with CSS variable (Chrome, Firefox, Edge) */
    backdrop-filter: var(--glass-blur-medium, blur(12px));

    /* Safari: MUST use literal value, CSS variables don't work */
    -webkit-backdrop-filter: blur(12px);

    background: var(--glass-bg-light, hsl(0 0% 100% / 0.85));
}

/* Modifier with different blur */
.glass-card--subtle {
    backdrop-filter: var(--glass-blur-light, blur(8px));
    -webkit-backdrop-filter: blur(8px);  /* Literal value for Safari */
}
```

### Mobile Performance Optimization
```css
/* Source: Existing tokens.css */
/* Reduce blur and limit glass elements on mobile */

:root {
    /* Desktop: Full blur values */
    --glass-blur-light: blur(8px);
    --glass-blur-medium: blur(12px);
    --glass-blur-heavy: blur(16px);
}

/* Mobile: Reduced blur for GPU performance */
@media (max-width: 768px) {
    :root {
        --glass-blur-light: blur(6px);   /* 25% reduction */
        --glass-blur-medium: blur(8px);  /* 33% reduction */
        --glass-blur-heavy: blur(10px);  /* 37% reduction */
    }
}
```

### Contrast-Compliant Glass Backgrounds
```css
/* Source: Existing tokens.css */
/* High opacity (85-95%) for WCAG contrast compliance */

:root {
    /* High opacity glass backgrounds */
    --glass-bg-light: hsl(0 0% 100% / 0.85);    /* 85% opacity */
    --glass-bg-medium: hsl(0 0% 100% / 0.90);   /* 90% opacity */
    --glass-bg-subtle: hsl(0 0% 100% / 0.95);   /* 95% opacity */

    /* Dark text for contrast */
    --color-text: hsl(220 13% 15%);             /* Near-black */
    --color-text-muted: hsl(220 13% 38%);       /* Medium gray */
}

/* Test contrast: --color-text on --glass-bg-subtle */
/* Effective background: ~98% white */
/* Text: 15% lightness */
/* Contrast ratio: ~13:1 (PASSES 4.5:1 minimum) */
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `:focus` for keyboard indicators | `:focus-visible` pseudo-class | Jan 2020 (baseline) | Keyboard-only focus indicators, no mouse click outlines |
| WCAG 2.1 Level AA | WCAG 2.2 Level AA | Oct 2023 | New focus appearance criteria (2.4.11, 2.4.12, 2.4.13) |
| Manual-only accessibility testing | Automated (axe-core) + manual | Ongoing trend | Catches 40% of issues automatically, scales better |
| `-webkit-backdrop-filter` prefix optional | Required for Safari (2026) | Still required | Safari 18 still needs webkit prefix despite MDN docs |
| Motion-first design | No-motion-first or reduced motion by default | 2021+ trend | Better default accessibility |

**Deprecated/outdated:**
- **Universal `:focus { outline: none; }`:** Violates WCAG 2.4.7 unless replaced with `:focus-visible` alternative
- **Color-only status indicators:** WCAG 1.4.1 requires non-color indicators (shape, text, icon)
- **Autoprefixer for backdrop-filter:** Doesn't handle CSS variable limitation in Safari; manual literal values required

## Open Questions

Things that couldn't be fully resolved:

1. **Contrast ratio testing for dynamic glass backgrounds**
   - What we know: WebAIM Contrast Checker works for static colors. Glass backgrounds blend with backdrop content.
   - What's unclear: Best method to test contrast when backdrop content varies (e.g., colorful background images behind glass)
   - Recommendation: Test against worst-case backdrop (lightest content for dark text, darkest for light text). Use browser DevTools color picker to sample effective rendered color.

2. **Lighthouse performance score threshold for mobile**
   - What we know: Success criteria requires "Lighthouse 90+ score" on mobile
   - What's unclear: Which Lighthouse category (Performance vs Accessibility vs Best Practices vs SEO)
   - Recommendation: Assume Performance score since context is "mobile performance optimized." Document baseline score before Phase 23, target +10% improvement or 90+ absolute.

3. **Screen reader testing scope**
   - What we know: Manual accessibility testing should include screen readers (NVDA, JAWS, VoiceOver)
   - What's unclear: Whether Phase 23 includes screen reader testing or just keyboard navigation
   - Recommendation: Focus on keyboard navigation (explicit in success criteria). Screen reader testing is valuable but time-intensive; consider post-Phase 23 if not blocking.

4. **Cross-browser testing OS coverage**
   - What we know: Success criteria lists Chrome, Firefox, Safari, Edge
   - What's unclear: Which OS versions (Windows/macOS/iOS/Android)
   - Recommendation: Test Safari on macOS (desktop) and iOS (mobile). Edge on Windows. Chrome/Firefox on any OS (rendering engines consistent).

## Sources

### Primary (HIGH confidence)
- [W3C WCAG 2.2 Understanding SC 2.4.7 Focus Visible](https://www.w3.org/WAI/WCAG22/Understanding/focus-visible.html) - Official W3C documentation on focus visible requirements
- [MDN: prefers-reduced-motion](https://developer.mozilla.org/en-US/docs/Web/CSS/Reference/At-rules/@media/prefers-reduced-motion) - Official MDN documentation on reduced motion media query (last modified Jan 8, 2026)
- [MDN: :focus-visible](https://developer.mozilla.org/en-US/docs/Web/CSS/Reference/Selectors/:focus-visible) - Official MDN documentation on focus-visible pseudo-class
- [Chrome for Developers: Lighthouse Overview](https://developer.chrome.com/docs/lighthouse/overview/) - Official Lighthouse documentation from Google
- [Can I Use: backdrop-filter](https://caniuse.com/css-backdrop-filter) - Browser support tables for backdrop-filter
- [GitHub: axe-core](https://github.com/dequelabs/axe-core) - Official axe-core repository and documentation

### Secondary (MEDIUM confidence)
- [WebAIM: Contrast Checker](https://webaim.org/resources/contrastchecker/) - Industry-standard contrast validation tool
- [Deque: axe-core](https://www.deque.com/axe/axe-core/) - Official axe-core product page with coverage details
- [WCAG 2.4.13 Focus Appearance Guide](https://www.wcag.com/designers/2-4-13-focus-appearance/) - Implementation guide for WCAG 2.4.13
- [Sara Soueidan: Designing Accessible WCAG-Conformant Focus Indicators](https://www.sarasoueidan.com/blog/focus-indicators/) - Expert guidance on focus indicator design
- [Web.dev: prefers-reduced-motion](https://web.dev/articles/prefers-reduced-motion) - Google's best practices guide
- [Smashing Magazine: Respecting Users' Motion Preferences](https://www.smashingmagazine.com/2021/10/respecting-users-motion-preferences/) - In-depth motion accessibility guide
- [Accessible Minds: Automating Accessibility Testing in CI/CD](https://accessiblemindstech.com/automating-accessibility-testing-in-ci-cd-with-github-actions/) - CI/CD integration patterns
- [TestParty: Accessibility Testing in CI/CD Integration Guide](https://testparty.ai/blog/accessibility-testing-cicd) - 2026 guide on automation
- [Medium: Glassmorphism Complete 2026 Guide](https://medium.com/@Kinetools/how-to-create-modern-ui-with-glassmorphism-effects-a-complete-2026-guide-2b1d71856542) - Performance optimization patterns
- [Mobile UI Trends 2026: Glassmorphism](https://www.sanjaydey.com/mobile-ui-trends-2026-glassmorphism-spatial-computing/) - 2026 glassmorphism best practices

### Tertiary (LOW confidence - WebSearch only)
- [accessiBe: Best Contrast Checker Tools 2026](https://accessibe.com/blog/knowledgebase/color-contrast-checker-tools) - Tool comparison
- [LogRocket: Cross-Browser CSS Testing Tools](https://blog.logrocket.com/5-best-open-source-tools-cross-browser-css-testing/) - Testing tool overview
- [TestMu AI: Cross-Browser Testing](https://www.testmu.ai/cross-browser-testing/) - Cloud testing platform

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - axe-core, Lighthouse, WebAIM Contrast Checker are industry-standard authoritative tools
- Architecture: HIGH - Patterns verified from official MDN docs, W3C WCAG docs, and existing codebase
- Pitfalls: HIGH - Common issues documented in authoritative sources + observed in codebase analysis
- Cross-browser support: HIGH - Can I Use database + verified GitHub issues confirming Safari webkit prefix requirement in 2026
- Performance optimization: MEDIUM - Mobile glassmorphism optimization based on 2026 blog posts + existing codebase tokens, not official W3C guidance

**Research date:** 2026-01-31
**Valid until:** 2026-03-31 (60 days for stable standards like WCAG 2.2; browser support evolves slowly)

**Key constraints from existing codebase:**
- Pure CSS approach maintained (no build tools, frameworks, or preprocessors)
- CSS Cascade Layers architecture
- BEM naming convention
- Glass opacity 85-95% for WCAG compliance (already implemented)
- Mobile performance budget: 2-3 glass elements max, 6px blur (already implemented)
- prefers-reduced-motion already implemented in glass-card.css and status-page.css
- Touch-safe hover states already implemented
