# Pitfalls Research: Bento Grid and Glassmorphism UI Redesign

**Project:** qm-nmr-calc v2.1 UI Redesign
**Domain:** Scientific web application UI with interactive elements
**Researched:** 2026-01-29
**Overall confidence:** HIGH (glassmorphism), MEDIUM (bento grid specifics)

## Executive Summary

Implementing bento grid layouts with glassmorphism effects for a scientific web application presents unique challenges beyond typical marketing sites. The combination of interactive forms, canvas-based 3D visualization (3Dmol.js), data tables, and file uploads requires careful attention to contrast, performance, z-index management, and responsive breakpoints. Critical risks include accessibility failures (WCAG contrast violations), performance degradation (especially on mobile), and breaking existing interactive functionality.

---

## Critical Pitfalls

These mistakes can cause accessibility violations, major performance issues, or broken functionality requiring significant rework.

### Pitfall 1: WCAG Contrast Violations on Glassmorphic Surfaces

**What goes wrong:**
Glassmorphism relies on transparency and backdrop-filter blur, which naturally reduce contrast. Text on semi-transparent "glass" surfaces often fails to meet WCAG 2.2 requirements (4.5:1 for normal text, 3:1 for large text/UI elements). The problem is worse in light mode where backgrounds can be unpredictable, and when the background behind the glass has color variation (e.g., gradients, images, or colored cards).

**Why it happens:**
- Translucent backgrounds allow varying colors to show through
- Backdrop blur reduces edge sharpness, making low-contrast text harder to read
- Background content shifts as users scroll, causing text readability to vary dynamically
- Light mode designs typically use white/light glass over light backgrounds, reducing contrast

**Consequences:**
- Inaccessible to users with visual impairments, color blindness, or age-related vision decline
- Form labels and data table text become unreadable
- WCAG compliance failure (legal/institutional requirements for scientific tools)
- Cognitive load increases for all users when contrast is marginal

**Prevention:**
1. **Add semi-opaque overlay to glass surfaces**: Pair blur with 10-30% opacity solid white/dark tint
   ```css
   .glass-card {
     background: rgba(255, 255, 255, 0.85); /* NOT just 0.3 */
     backdrop-filter: blur(10px);
   }
   ```
2. **Test contrast with WCAG tools** on every glassmorphic element with text
3. **Use dark text on light glass** (light mode) or light text on dark glass
4. **Provide fallback solid backgrounds** for browsers without backdrop-filter support:
   ```css
   .glass-card {
     background: rgba(255, 255, 255, 0.95); /* Fallback */
     backdrop-filter: blur(10px);
   }
   @supports not (backdrop-filter: blur(10px)) {
     .glass-card {
       background: rgba(255, 255, 255, 1); /* Solid fallback */
     }
   }
   ```
5. **Avoid glassmorphism on critical interactive elements** like form inputs, buttons, error messages

**Detection (warning signs):**
- Squinting required to read text on glass surfaces
- Text readability varies when scrolling
- Contrast checker tools report ratios below 4.5:1
- Stakeholders comment "looks cool but hard to read"

**Which phase should address this:**
- Phase 1 (CSS framework setup): Define glass effect classes with proper opacity
- Phase 2 (Results page redesign): Test contrast on all data tables and shift annotations
- Phase 3 (Submit page): Ensure form labels meet contrast requirements

---

### Pitfall 2: Mobile Performance Degradation and Battery Drain

**What goes wrong:**
`backdrop-filter: blur()` is GPU-intensive and dramatically impacts performance on mobile devices. Multiple glass elements, large blur radii, or animating blur effects cause frame drops, device heating, and rapid battery drain. The problem compounds when glass effects are applied to large areas or many elements simultaneously.

**Why it happens:**
- Browser must render entire scene behind element, apply blur filter, then composite on top
- Blur operations are exponentially more expensive at higher pixel counts and blur radii
- Mobile GPUs have limited power compared to desktop
- Each additional glass layer multiplies rendering cost
- Scrolling with fixed glass elements triggers continuous recomposition

**Consequences:**
- Janky scrolling and animation stutter on mobile
- Device overheating during extended use
- Poor user experience on tablets/phones
- Battery drain makes app unusable for field work
- 3Dmol.js canvas performance degrades when glass elements overlay or are nearby

**Prevention:**
1. **Limit blur radius**: Keep between 6-10px on mobile (vs 10-20px desktop)
   ```css
   .glass-card {
     backdrop-filter: blur(10px);
   }
   @media (max-width: 768px) {
     .glass-card {
       backdrop-filter: blur(6px);
     }
   }
   ```
2. **Reduce glass element count**: Maximum 2-3 glassmorphic elements per viewport on mobile
3. **Never animate backdrop-filter**: Animates poorly even on desktop
4. **Avoid glass on large areas**: Small cards only, not full-page backgrounds
5. **Test on lower-end Android devices** (not just iPhone/iPad)
6. **Use `will-change` sparingly**: Only on elements that actually animate, remove after animation completes
7. **Consider disabling glass on mobile entirely**:
   ```css
   @media (max-width: 768px) and (prefers-reduced-motion: reduce) {
     .glass-card {
       backdrop-filter: none;
       background: rgba(255, 255, 255, 0.95);
     }
   }
   ```

**Detection (warning signs):**
- Frame rate drops during scroll testing on mobile
- Device warm to touch after 5 minutes of use
- 3D molecule viewer lags when glass elements are visible
- Browser developer tools show high GPU usage

**Which phase should address this:**
- Phase 1 (CSS framework): Build mobile-optimized variants with reduced blur
- Phase 4 (Results page): Test 3Dmol.js canvas performance with glass cards
- Phase 6 (Responsive testing): Performance audit on Android mid-range device

---

### Pitfall 3: Z-Index and Stacking Context Chaos

**What goes wrong:**
Applying `backdrop-filter` creates a new stacking context, which breaks carefully constructed z-index hierarchies. Interactive elements (dropdowns, modals, tooltips, file upload previews) suddenly appear behind glass cards, or the 3Dmol.js canvas renders incorrectly relative to glass overlays.

**Why it happens:**
- `backdrop-filter` implicitly creates stacking context (like `transform`, `opacity < 1`, `filter`)
- Multiple nested stacking contexts make z-index values local, not global
- Each glass card becomes a containing block for positioned descendants
- 3Dmol.js canvas elements have their own stacking requirements

**Consequences:**
- Dropdown menus hidden behind glass cards
- Modal dialogs don't properly overlay the page
- File upload previews disappear
- 3D molecule viewer canvas overlaps incorrectly with UI chrome
- Interactive tooltips on NMR shift annotations don't display

**Prevention:**
1. **Document stacking context hierarchy** in CSS comments
2. **Use consistent z-index scale**: Define ranges for different UI layers
   ```css
   /* Z-index scale:
    * 0-99: Base content
    * 100-199: Glass cards
    * 200-299: Dropdowns/popovers
    * 300-399: Modals
    * 400+: Notifications
    */
   .glass-card { z-index: 100; }
   .dropdown-menu { z-index: 200; }
   ```
3. **Test all interactive states**: Hover dropdowns, file pickers, tooltips
4. **Avoid deeply nested glass elements**: Keep hierarchy flat
5. **Use portal/teleport pattern for modals**: Render modals at root level, outside glass card DOM tree
6. **Test 3Dmol.js canvas independently**: Verify it renders above/below glass as intended

**Detection (warning signs):**
- UI elements disappear when hovering/clicking
- Dropdown menus cut off by parent containers
- Modals appear behind page content
- Canvas rendering order looks wrong

**Which phase should address this:**
- Phase 2 (Results page): Test 3Dmol.js canvas z-index with surrounding glass cards
- Phase 3 (Submit page): Verify file upload dropzone and form validation tooltips
- Phase 5 (Status page): Ensure progress indicators display correctly

---

## Moderate Pitfalls

These mistakes cause technical debt, layout bugs, or degraded UX but are recoverable without major rewrites.

### Pitfall 4: Bento Grid Content Overflow and Text Truncation

**What goes wrong:**
CSS Grid items containing overflowing text expand to fit content instead of truncating with ellipsis. Long molecule names, file paths, or SMILES strings break grid layouts. Grid cells become misaligned when some cells have more content than others.

**Why it happens:**
- Grid items have implicit `min-width: auto`, allowing them to grow beyond track size
- Text overflow requires explicit width constraints
- Flex/grid children need `min-width: 0` to honor parent constraints
- `text-overflow: ellipsis` doesn't work without `white-space: nowrap` and `overflow: hidden`

**Consequences:**
- Grid layout breaks on wide content (long SMILES strings, file names)
- Horizontal scrollbars appear unexpectedly
- Cards in bento grid become misaligned
- Responsive breakpoints fail because content forces expansion

**Prevention:**
1. **Set `min-width: 0` on grid items**:
   ```css
   .bento-grid {
     display: grid;
     grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
   }
   .bento-grid > * {
     min-width: 0; /* Allow shrinking below content size */
   }
   ```
2. **Use text truncation pattern for dynamic content**:
   ```css
   .truncate {
     white-space: nowrap;
     overflow: hidden;
     text-overflow: ellipsis;
   }
   ```
3. **Add word-break for long tokens without spaces**:
   ```css
   .smiles-string {
     word-break: break-all; /* Break SMILES mid-token if needed */
   }
   ```
4. **Test with extreme content**: Very long molecule names, 500+ character SMILES

**Detection (warning signs):**
- Grid layout shifts when content changes
- Horizontal scrollbars on grid containers
- Cards of unequal height in asymmetric bento layouts
- Text overflows card boundaries

**Which phase should address this:**
- Phase 2 (Results page): Test with long molecule names and SMILES strings
- Phase 3 (Submit page): Verify file name display in upload widget
- Phase 6 (Responsive): Test content overflow at all breakpoints

---

### Pitfall 5: Responsive Breakpoint Collapse Chaos

**What goes wrong:**
Bento grids with asymmetric layouts (e.g., 1 large card + 3 small cards) break awkwardly at tablet breakpoints. Media queries override grid layouts incorrectly, leaving whitespace gaps or cards that don't collapse to single column on mobile. Glassmorphic effects become cluttered when too many elements stack vertically.

**Why it happens:**
- Asymmetric grid designs don't map cleanly to 2-column tablet or 1-column mobile
- `grid-template-columns` needs different definitions at each breakpoint
- `min-width` in media queries may not apply to larger breakpoints
- Visual hierarchy gets lost when large "hero" cards shrink to same size as detail cards

**Consequences:**
- Awkward layouts at 768px-1024px (tablet range)
- Too much vertical scrolling on mobile
- Visual hierarchy flattened (all cards same size)
- Glassmorphism feels cluttered with many stacked cards

**Prevention:**
1. **Design mobile-first, enhance upward**:
   ```css
   /* Mobile: single column */
   .bento-grid {
     display: grid;
     grid-template-columns: 1fr;
     gap: 1rem;
   }
   /* Tablet: 2 columns */
   @media (min-width: 768px) {
     .bento-grid {
       grid-template-columns: repeat(2, 1fr);
     }
   }
   /* Desktop: asymmetric bento */
   @media (min-width: 1024px) {
     .bento-grid {
       grid-template-columns: 2fr 1fr 1fr;
     }
   }
   ```
2. **Use `grid-auto-flow: dense`** to fill gaps in asymmetric layouts
3. **Test on real tablet devices** (not just browser resize)
4. **Reduce glass element count on mobile**: Some cards become solid backgrounds
5. **Maintain visual hierarchy**: Hero cards remain larger even on tablet

**Detection (warning signs):**
- Large whitespace gaps at 768px width
- Cards don't fit viewport width on mobile
- Excessive vertical scrolling
- Layout feels cramped between breakpoints

**Which phase should address this:**
- Phase 4 (Responsive design): Define breakpoints for all pages
- Phase 6 (Testing): Validate on iPad, Android tablet, mobile devices

---

### Pitfall 6: Form Input and Interactive Element Usability

**What goes wrong:**
Glassmorphic form inputs have poor visual affordance (hard to tell what's clickable), low contrast between input field and background, and unclear focus states. File upload dropzones with glass effects don't provide sufficient visual feedback. Disabled states are indistinguishable from enabled states.

**Why it happens:**
- Translucent backgrounds reduce visual distinction between input and surrounding card
- Placeholder text often too low contrast
- Focus rings get lost against blurred backgrounds
- Glass effect makes borders subtle/invisible

**Consequences:**
- Users don't recognize form inputs
- No clear feedback when clicking/typing
- Accessibility issues for keyboard navigation (focus invisible)
- File uploads feel unresponsive

**Prevention:**
1. **Use solid backgrounds for form inputs**, not glass:
   ```css
   .glass-card input {
     background: white; /* Solid, not transparent */
     border: 1px solid rgba(0, 0, 0, 0.2);
   }
   ```
2. **Strong focus indicators**:
   ```css
   input:focus {
     outline: 3px solid #0066cc;
     outline-offset: 2px;
   }
   ```
3. **Clear disabled states**:
   ```css
   input:disabled {
     opacity: 0.6;
     cursor: not-allowed;
   }
   ```
4. **File upload visual feedback**: Change background on drag-over
5. **Test keyboard navigation**: Tab through all inputs, verify focus visibility

**Detection (warning signs):**
- Inputs hard to locate visually
- Clicking input provides no feedback
- Focus state invisible when tabbing
- File drag-drop feels unresponsive

**Which phase should address this:**
- Phase 3 (Submit page): Design form inputs with solid backgrounds
- Phase 5 (Interactive states): Define focus, hover, disabled states

---

### Pitfall 7: 3Dmol.js Canvas Rendering Conflicts

**What goes wrong:**
Glassmorphic cards overlaying or adjacent to the 3Dmol.js canvas cause rendering glitches. The canvas element may flicker, appear above glass when it should be behind, or experience frame rate drops when glass elements are visible. Backdrop-filter blur over canvas areas causes WebGL context loss on some browsers.

**Why it happens:**
- WebGL canvases use hardware acceleration separate from CSS rendering
- Compositing layers conflict between canvas and backdrop-filter
- GPU resources shared between WebGL and CSS blur operations
- Browser compositing order unpredictable with canvas + glass

**Consequences:**
- 3D molecule viewer flickers or renders incorrectly
- Performance degradation when both canvas and glass are visible
- WebGL context loss errors
- Visual glitches during canvas manipulation (rotation, zoom)

**Prevention:**
1. **Never apply backdrop-filter directly over canvas**:
   ```css
   /* Bad */
   .glass-overlay {
     backdrop-filter: blur(10px);
     position: absolute;
     top: 0; /* Covers canvas */
   }

   /* Good */
   .glass-card.adjacent-to-canvas {
     backdrop-filter: none; /* Disable near canvas */
   }
   ```
2. **Use isolation for canvas containers**:
   ```css
   .canvas-container {
     isolation: isolate;
     z-index: 50;
   }
   ```
3. **Test canvas interactions** with glass visible: Rotation, zoom, resize
4. **Consider disabling glass effects on Results page** where canvas is prominent
5. **Monitor for WebGL context loss**: Add error handlers in 3Dmol.js init

**Detection (warning signs):**
- Canvas flickers during page load
- 3D viewer frame rate drops with glass cards visible
- Console errors about WebGL context
- Visual glitches when rotating molecule

**Which phase should address this:**
- Phase 2 (Results page): Test 3Dmol.js with glass card layouts
- Phase 4 (Integration testing): Verify canvas performance with all effects enabled

---

## Minor Pitfalls

Annoying but fixable issues that don't block launch.

### Pitfall 8: Safari Webkit Prefix Requirement

**What goes wrong:**
`backdrop-filter` doesn't work in Safari without the `-webkit-` prefix. Designs look perfect in Chrome/Firefox but have no glass effect in Safari, falling back to plain transparent backgrounds.

**Why it happens:**
- Safari requires `-webkit-backdrop-filter` for compatibility
- Standard `backdrop-filter` not recognized in older Safari versions
- Easy to forget when developing primarily in Chrome

**Consequences:**
- Glass effect missing in Safari
- Lower visual appeal on macOS/iOS
- Inconsistent experience across browsers

**Prevention:**
1. **Always include both prefixed and unprefixed**:
   ```css
   .glass-card {
     -webkit-backdrop-filter: blur(10px);
     backdrop-filter: blur(10px);
   }
   ```
2. **Use autoprefixer**: Automate with PostCSS build step
3. **Test in Safari early**: Don't wait until end of project

**Detection (warning signs):**
- Glass effects missing in Safari
- Backgrounds fully transparent without blur

**Which phase should address this:**
- Phase 1 (CSS framework): Include prefixes in all glass classes

---

### Pitfall 9: Overly Aggressive Blur Values

**What goes wrong:**
Using `blur(20px)` or higher creates muddy, unclear UIs. Text becomes harder to read, visual hierarchy is lost, and performance suffers disproportionately (blur cost scales exponentially).

**Why it happens:**
- Designers excited by effect, push values too high
- "More blur = more glass effect" assumption

**Consequences:**
- Muddy appearance
- Poor readability
- Worse performance for marginal visual gain

**Prevention:**
1. **Keep blur between 6-15px**: Sweet spot for glass effect without mud
2. **Use 8-10px as default**: Subtle but effective
3. **Reserve 15px+ for hero sections only**: Not repeated elements

**Detection (warning signs):**
- UI feels muddy or unclear
- Text hard to read even with good contrast
- Stakeholders say "too blurry"

**Which phase should address this:**
- Phase 1 (CSS framework): Define standard blur values as CSS variables

---

### Pitfall 10: Forgetting Solid Background Fallback

**What goes wrong:**
Older browsers or browsers with backdrop-filter disabled (e.g., Firefox with WebRender disabled) show fully transparent backgrounds instead of glass, making text unreadable.

**Why it happens:**
- Developers test only in modern Chrome
- Fallback styles not defined for unsupported browsers

**Consequences:**
- Broken UI in older browsers
- Unreadable text on transparent backgrounds

**Prevention:**
1. **Always provide solid fallback**:
   ```css
   .glass-card {
     background: rgba(255, 255, 255, 0.95); /* Fallback */
   }
   @supports (backdrop-filter: blur(10px)) {
     .glass-card {
       background: rgba(255, 255, 255, 0.75); /* With blur */
       backdrop-filter: blur(10px);
     }
   }
   ```
2. **Test with backdrop-filter disabled**: Firefox about:config flag

**Detection (warning signs):**
- Text unreadable in Firefox (older versions)
- Transparent backgrounds without blur

**Which phase should address this:**
- Phase 1 (CSS framework): Include @supports fallbacks

---

## Browser Compatibility Summary

### backdrop-filter Support (2026)

**Excellent support:** 96% global coverage

| Browser | Support | Notes |
|---------|---------|-------|
| Chrome 76+ | Full | Unprefixed |
| Firefox 103+ | Full | Enabled by default since Firefox 103 |
| Safari 9+ | Full | Requires `-webkit-` prefix |
| Edge 17+ | Full | Chromium-based |
| IE 11 | None | No support, fallback required |

**Key compatibility requirements:**
- Always include `-webkit-backdrop-filter` for Safari
- Provide solid background fallback with `@supports`
- Test on Firefox 102 and earlier if targeting older browsers

### CSS Grid Support (2026)

**Near-universal support:** 98% global coverage

| Browser | Support | Notes |
|---------|---------|-------|
| Chrome 57+ | Full | Subgrid since 117 (Sept 2023) |
| Firefox 52+ | Full | Subgrid since 71 (Dec 2019) |
| Safari 10.1+ | Full | Subgrid since 16 (Sept 2022) |
| Edge 16+ | Full | Subgrid since 117 (Sept 2023) |
| IE 11 | Partial | Old spec with `-ms-` prefixes, avoid |

**Key compatibility notes:**
- Subgrid has 97% support in 2026 (safe to use)
- IE 11 supports outdated Grid spec, don't bother with polyfills
- Use `@supports` for progressive enhancement if needed

---

## Performance Guidelines

### Acceptable Use

| Scenario | Blur Radius | Element Count | Device | Notes |
|----------|-------------|---------------|--------|-------|
| Desktop hero card | 10-15px | 1-2 | Desktop | Large cards, not repeated |
| Desktop data cards | 8-10px | 3-5 | Desktop | Moderate use |
| Tablet cards | 6-8px | 2-3 | Tablet | Reduced blur |
| Mobile cards | 6px or none | 1-2 | Mobile | Minimal glass effect |
| Adjacent to canvas | 0px | 0 | All | Avoid near WebGL |

### Performance Budget

**Maximum per viewport:**
- Desktop: 5 glass elements with 10px blur
- Tablet: 3 glass elements with 8px blur
- Mobile: 2 glass elements with 6px blur

**Never:**
- Animate backdrop-filter
- Apply glass to full-page backgrounds
- Use blur > 15px on repeated elements
- Stack more than 2 glass layers

---

## Accessibility Checklist

Before shipping any glass effect:

- [ ] Text contrast ratio >= 4.5:1 (normal text)
- [ ] Large text / UI elements >= 3:1
- [ ] Tested with WCAG color contrast analyzer tool
- [ ] Focus indicators visible on all interactive elements
- [ ] Form inputs have solid backgrounds (not glass)
- [ ] Disabled states visually distinct
- [ ] Works with browser zoom at 200%
- [ ] Screen reader testing completed (glass doesn't affect semantics)
- [ ] Keyboard navigation works (focus visible)
- [ ] Respects `prefers-reduced-motion` (consider disabling blur)

---

## Phase-Specific Warnings

### Phase 1: CSS Framework Setup
- Define glass effect classes with proper fallbacks
- Include Safari webkit prefix
- Set blur value range (6-15px via CSS variables)
- Create solid background fallbacks for all glass styles
- Document z-index scale for stacking contexts

### Phase 2: Results Page Redesign
- Test 3Dmol.js canvas with glass cards (z-index, performance)
- Verify data table contrast on glass backgrounds
- Test NMR shift annotation tooltip visibility
- Validate SMILES string truncation in grid cells
- Performance test on mobile device

### Phase 3: Submit Page Redesign
- Use solid backgrounds for all form inputs (no glass)
- Strong focus indicators for keyboard navigation
- Test file upload dropzone visual feedback
- Verify form validation message visibility
- Test long file name display with text truncation

### Phase 4: Status Page Redesign
- Ensure progress indicators contrast properly
- Test polling updates don't cause glass re-renders
- Verify loading states visible on glass backgrounds

### Phase 5: Responsive Design Implementation
- Define mobile-first breakpoints (320px, 768px, 1024px)
- Reduce glass element count on mobile
- Reduce blur radius on mobile (6px)
- Test tablet layouts (768-1024px range)
- Verify bento grid collapses cleanly to single column

### Phase 6: Cross-Browser and Performance Testing
- Safari testing (webkit prefix)
- Firefox with backdrop-filter disabled
- Performance audit on mid-range Android device
- Battery drain test (30 minutes continuous use)
- Lighthouse accessibility audit
- WCAG contrast validation

---

## Testing Protocol

### Required Tests Before Shipping

**Contrast Testing:**
1. Run WCAG Color Contrast Analyzer on every text element on glass
2. Test with different page backgrounds (white, gray, gradient)
3. Verify contrast at all responsive breakpoints

**Performance Testing:**
1. Lighthouse performance audit (target: 90+ score)
2. Manual testing on Android mid-range device
3. GPU usage monitoring with DevTools
4. Battery drain test: 30 minutes continuous interaction
5. Frame rate testing during scroll (target: 60fps on desktop, 30fps on mobile)

**Interaction Testing:**
1. Tab through all form inputs (focus visible)
2. Test all dropdowns, tooltips, modals (correct z-index)
3. 3Dmol.js canvas interaction (rotation, zoom) with glass visible
4. File upload drag-drop visual feedback
5. Hover states on all interactive elements

**Browser Testing:**
1. Chrome (latest)
2. Firefox (latest + version 102 for fallback test)
3. Safari (latest, verify webkit prefix)
4. Edge (latest)
5. Mobile Safari (iOS)
6. Chrome Mobile (Android)

**Responsive Testing:**
1. Desktop 1920x1080
2. Laptop 1366x768
3. Tablet 768x1024
4. Mobile 375x667
5. Test breakpoint transitions (767px, 768px, 1023px, 1024px)

---

## Recovery Strategies

If you encounter these issues late:

**Contrast violations found:**
- Increase background opacity (0.85-0.95 range)
- Add semi-opaque tint overlay
- Switch to solid backgrounds for affected elements
- Reduce blur radius to improve edge sharpness

**Performance issues found:**
- Reduce blur radius by 30% (10px → 7px)
- Decrease glass element count (remove from less important cards)
- Disable backdrop-filter on mobile entirely
- Use `@media (prefers-reduced-motion: reduce)` to disable blur

**Z-index chaos:**
- Flatten stacking context hierarchy
- Move modals to root level with portal pattern
- Document all stacking contexts in comments
- Use consistent z-index scale

**Bento grid breaks:**
- Simplify to symmetric grid on mobile
- Use single-column layout below 768px
- Test with extreme content lengths
- Add `min-width: 0` to all grid items

---

## Sources

### Glassmorphism and Accessibility
- [Glassmorphism Meets Accessibility: Can Glass Be Inclusive? | Axess Lab](https://axesslab.com/glassmorphism-meets-accessibility-can-frosted-glass-be-inclusive/)
- [12 Glassmorphism UI Features, Best Practices, and Examples | UXPilot](https://uxpilot.ai/blogs/glassmorphism-ui)
- [Glassmorphism: Definition and Best Practices | Nielsen Norman Group](https://www.nngroup.com/articles/glassmorphism/)
- [Dark Glassmorphism: The Aesthetic That Will Define UI in 2026 | Medium](https://medium.com/@developer_89726/dark-glassmorphism-the-aesthetic-that-will-define-ui-in-2026-93aa4153088f)
- [Glassmorphism with Website Accessibility in Mind | New Target](https://www.newtarget.com/web-insights-blog/glassmorphism/)

### Performance and Browser Compatibility
- [CSS Backdrop Filter | Can I use](https://caniuse.com/css-backdrop-filter)
- [backdrop-filter | MDN Web Docs](https://developer.mozilla.org/en-US/docs/Web/CSS/backdrop-filter)
- [backdrop-filter: blur causes battery drain/lag | GitHub Issue](https://github.com/Matchstic/Xen-HTML/issues/219)
- [Cross Browser Compatibility Score of CSS Backdrop Filter | LambdaTest](https://www.lambdatest.com/web-technologies/css-backdrop-filter)
- [Web Image Effects Performance Showdown | Smashing Magazine](https://www.smashingmagazine.com/2016/05/web-image-effects-performance-showdown/)

### Bento Grid Layouts
- [Bite-sized bento grid UX designs | LogRocket Blog](https://blog.logrocket.com/ux-design/bento-grids-ux/)
- [Best Bento Grid Design Examples [2026] | Mockuuups Studio](https://mockuuups.studio/blog/post/best-bento-grid-design-examples/)
- [Creating Bento Grid Layouts | ibelick](https://ibelick.com/blog/create-bento-grid-layouts)
- [Build a bento layout with CSS grid | iamsteve](https://iamsteve.me/blog/bento-layout-css-grid)

### CSS Grid and Text Overflow
- [CSS Text Overflow Ellipsis Not Working | TheLinuxCode](https://thelinuxcode.com/css-text-overflow-ellipsis-not-working-fixed/)
- [Overflow text breaks out of CSS Grid cells | Mozilla Bugzilla](https://bugzilla.mozilla.org/show_bug.cgi?id=1528030)
- [CSS Subgrid Browser Support 2025-2026 | FrontendTools](https://www.frontendtools.tech/blog/mastering-css-grid-2025)
- [Modern CSS Layout Techniques | FrontendTools](https://www.frontendtools.tech/blog/modern-css-layout-techniques-flexbox-grid-subgrid-2025)

### Form Inputs and Interactive Elements
- [Design Best Glassmorphism Login Form Using HTML CSS | 2026](https://abduldev.com/design-a-glassmorphism-login-form-using-html-and-css/)
- [How to Create Modern UI with Glassmorphism Effects | Medium](https://medium.com/@Kinetools/how-to-create-modern-ui-with-glassmorphism-effects-a-complete-2026-guide-2b1d71856542)

### WebGL and Canvas Performance
- [Real-Time Dashboard Performance: WebGL vs Canvas | Dev3lop](https://dev3lop.com/real-time-dashboard-performance-webgl-vs-canvas-rendering-benchmarks/)
- [Optimizing canvas | MDN Web Docs](https://developer.mozilla.org/en-US/docs/Web/API/Canvas_API/Tutorial/Optimizing_canvas)
- [CSS blur filter is order of magnitude slower than Chrome | Mozilla Bugzilla](https://bugzilla.mozilla.org/show_bug.cgi?id=925025)

---

**Confidence Assessment:**

| Area | Confidence | Source Quality |
|------|-----------|----------------|
| Glassmorphism accessibility | HIGH | Multiple authoritative sources (Nielsen Norman, Axess Lab) + MDN |
| backdrop-filter performance | HIGH | Can I Use data, MDN docs, GitHub issues, browser bug trackers |
| Browser compatibility | HIGH | Can I Use verified data (96% support), MDN compatibility tables |
| Bento grid responsive issues | MEDIUM | WebSearch + developer community patterns, less authoritative |
| Canvas + glass conflicts | MEDIUM | Performance research + WebGL best practices, needs project-specific testing |
| Z-index stacking context | HIGH | MDN specification, well-documented CSS behavior |

**Overall Assessment:** Research is comprehensive for glassmorphism pitfalls (accessibility, performance, browser compat). Bento grid research is solid but less critical since grid layout is well-established. Canvas-specific interactions need testing during implementation but general guidance is sound.
