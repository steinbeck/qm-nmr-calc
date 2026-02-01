# Phase 31: Documentation Polish and Cross-References - Research

**Researched:** 2026-02-01
**Domain:** Documentation consistency, cross-referencing, and link validation
**Confidence:** HIGH

## Summary

This phase focuses on the final polish of the v2.2 documentation milestone. The research examined the existing documentation files (README.md, docs/installation.md, docs/usage.md, docs/architecture.md, docs/libraries.md, docs/science.md), identified cross-reference patterns and gaps, and compiled best practices for ensuring consistency.

The existing documentation is already well-structured with a consistent audience-focused approach. Each document targets either academic researchers or developers and uses a clear pattern of introductory paragraphs, tables for structured data, code examples with source attribution, and "Related Documentation" sections. The primary work involves: (1) auditing and standardizing cross-references between documents, (2) ensuring consistent heading structure and formatting, (3) validating all internal and external links, and (4) verifying table of contents accuracy.

**Primary recommendation:** Systematically audit each document against a formatting checklist, add missing cross-references following the existing "Related Documentation" pattern, and validate all links using markdown-link-check before finalizing.

## Standard Stack

The established tools for documentation polish:

### Core
| Tool | Version | Purpose | Why Standard |
|------|---------|---------|--------------|
| markdown-link-check | latest | Link validation | Industry standard, GitHub Action integration |
| Manual review | - | Content consistency | Human judgment for cross-reference relevance |

### Supporting
| Tool | Version | Purpose | When to Use |
|------|---------|---------|-------------|
| linkcheckmd (Python) | 3.0+ | Alternative link checker | If Node.js tools not preferred |
| grep/bash | - | Heading structure audit | Quick consistency checks |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| markdown-link-check | md-dead-link-check | Python-based, pre-commit integration |
| Manual TOC check | doctoc | Automated but may not match custom formatting |

**Installation:**
```bash
# Node.js link checker
npm install -g markdown-link-check

# Or Python alternative
pip install linkcheckmd
```

## Architecture Patterns

### Documentation Structure Pattern

The existing documentation follows a consistent pattern that should be maintained:

```
docs/
├── installation.md   # Setup and prerequisites
├── usage.md          # User-facing workflow guide
├── architecture.md   # Developer internals
├── libraries.md      # Integration reference
└── science.md        # Methodology and theory
```

### Cross-Reference Pattern (Existing)

Each document already has a "Related Documentation" section at the bottom. The pattern is:

```markdown
## Related Documentation

- [Document Name](relative-path.md) - Brief description
- [Section Name](path.md#anchor) - Brief description
```

### Heading Hierarchy Pattern

Documents should follow this consistent structure:

```markdown
# Document Title                    <- H1: One per document
## Major Section                    <- H2: Top-level divisions
### Subsection                      <- H3: Within sections
#### Detail (sparingly)             <- H4: Only when necessary
```

**Key rules:**
- One H1 per document (the document title)
- Never skip heading levels (H2 -> H4 is invalid)
- Maximum depth of H4 (avoid H5/H6)
- Use sentence case for all headings

### Anti-Patterns to Avoid
- **Orphan documents:** Pages with no inbound or outbound links
- **Circular-only references:** Documents that only reference each other with no hub
- **Inconsistent anchor formats:** Mixing `#section-name` and `#section_name`
- **Redundant cross-references:** Linking to the same page multiple times in one section

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Link validation | Custom regex scripts | markdown-link-check | Handles edge cases, HTTP validation |
| TOC generation | Manual list maintenance | Document structure matches TOC | Existing structure is correct |
| Heading extraction | Custom parsing | grep -n "^#" | Simple and reliable |
| Anchor generation | Manual ID assignment | Auto-generated from heading text | GitHub/common markdown processors handle this |

**Key insight:** The existing documentation structure is already good. The polish phase is about verification and minor fixes, not restructuring.

## Common Pitfalls

### Pitfall 1: Broken Relative Links After Restructuring
**What goes wrong:** Moving files breaks relative links to other docs
**Why it happens:** Paths like `./architecture.md` or `../README.md` are fragile
**How to avoid:** Use consistent relative paths from document location, validate with link checker
**Warning signs:** 404 errors when navigating between docs

### Pitfall 2: Inconsistent Heading Levels
**What goes wrong:** Jumping from H2 to H4, or multiple H1s
**Why it happens:** Copy-paste from other documents, focus on visual appearance
**How to avoid:** Audit each document's heading structure systematically
**Warning signs:** Screen readers announce confusing structure, TOC looks wrong

### Pitfall 3: Stale Cross-References
**What goes wrong:** "See section X" references that no longer exist
**Why it happens:** Section headings change but references don't update
**How to avoid:** Search for all internal links and verify targets exist
**Warning signs:** #anchor links that scroll to wrong places or top of page

### Pitfall 4: External Link Rot
**What goes wrong:** External URLs (DOIs, GitHub, official docs) become unavailable
**Why it happens:** Websites restructure, projects move, pages deleted
**How to avoid:** Validate all external links before release
**Warning signs:** 404/301 responses, redirects to different content

### Pitfall 5: Inconsistent Formatting Conventions
**What goes wrong:** Some docs use `**bold**` for emphasis, others use `*italic*`
**Why it happens:** Multiple authoring sessions, no style guide reference
**How to avoid:** Document existing patterns, apply consistently
**Warning signs:** Visual inconsistency when reading sequentially

## Code Examples

Verified patterns from the existing documentation:

### Cross-Reference Pattern (from architecture.md)
```markdown
## Related Documentation

- [README](../README.md) - High-level overview and quick start
- [Installation Guide](installation.md) - System dependencies and setup
- [Usage Guide](usage.md) - Web UI and REST API reference
- [NMR Methodology](science.md) - DP4+, linear scaling, Boltzmann averaging
- [Library Documentation](libraries.md) - RDKit, NWChem, Huey integrations
```

### Section Link Pattern (from usage.md)
```markdown
**Prerequisites:** Complete the [Installation Guide](installation.md) before using this guide.
```

### Internal Section Reference
```markdown
For detailed methodology and scientific references, see the [Science documentation](science.md).
```

### Code Block with Source Attribution (from libraries.md)
```markdown
```python
# Source: src/qm_nmr_calc/validation.py, lines 10-39
from rdkit import Chem, RDLogger
...
```
```

## Formatting Checklist

Based on existing documentation patterns, verify each document against:

### Document Structure
- [ ] Single H1 title at top
- [ ] H2 for major sections
- [ ] No skipped heading levels
- [ ] Target audience stated in intro
- [ ] Related Documentation section at bottom

### Cross-References
- [ ] Links to related docs where context helps
- [ ] Section anchors work correctly
- [ ] No orphan pages (every doc linked from at least one other)
- [ ] README links to all docs in Documentation section

### Formatting Consistency
- [ ] Code blocks use consistent fencing (triple backticks)
- [ ] Language specified for syntax highlighting
- [ ] Tables use consistent alignment
- [ ] Lists use consistent punctuation (period vs no period)
- [ ] Emphasis uses consistent pattern (bold for UI elements, code for technical terms)

### Links
- [ ] All relative links resolve
- [ ] All external links return 200
- [ ] No duplicate link text with different targets
- [ ] DOI links use doi.org format

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Manual link checking | Automated with markdown-link-check | 2020+ | Catches broken links in CI |
| Static TOC | Auto-generated or heading-based | 2018+ | TOC stays in sync |
| Plain URLs | Descriptive link text | Accessibility standards | Better screen reader support |

**Deprecated/outdated:**
- Absolute URLs to same-repo files (use relative paths)
- Auto-numbered headings (use semantic structure)

## Existing Documentation Audit

### Current Cross-Reference Matrix

| Document | Links To | Linked From |
|----------|----------|-------------|
| README.md | installation, usage, architecture, libraries, science | (entry point) |
| installation.md | usage, architecture, science | README, usage |
| usage.md | installation, architecture, science | README, installation |
| architecture.md | README, installation, usage, science, libraries | README, usage |
| libraries.md | architecture, README | README, architecture |
| science.md | README | README, architecture, usage, installation |

### Identified Gaps
1. **science.md** only links back to README (needs Related Documentation section)
2. **libraries.md** has navigation footer but limited in-context links
3. **installation.md** could link to architecture.md for "how it works" context

### Heading Structure Audit

| Document | H1 | H2 | H3 | H4 | Issues |
|----------|----|----|----|----|--------|
| README.md | 1 | 8 | 0 | 0 | Clean |
| installation.md | 1 | 8 | 14 | 0 | Clean |
| usage.md | 1 | 12 | 26 | 0 | Clean |
| architecture.md | 1 | 10 | 16 | 0 | Clean |
| libraries.md | 1 | 6 | 19 | 0 | Clean |
| science.md | 1 | 10 | 16 | 0 | Clean |

All documents follow proper heading hierarchy - no issues found.

## Open Questions

Things that could not be fully resolved:

1. **External Link Validation Scope**
   - What we know: DOI links should work, GitHub URLs should resolve
   - What's unclear: Should we validate NWChem/RDKit doc links (may change frequently)?
   - Recommendation: Validate all external links, but document which are dependencies

2. **Anchor ID Format**
   - What we know: GitHub auto-generates from heading text (lowercase, hyphens)
   - What's unclear: Are there any manual anchor IDs in the docs?
   - Recommendation: Rely on auto-generated anchors, verify they work

## Sources

### Primary (HIGH confidence)
- Direct examination of existing documentation files
- [Google Developer Documentation Style Guide](https://developers.google.com/style) - Formatting standards
- [W3C WAI Headings Tutorial](https://www.w3.org/WAI/tutorials/page-structure/headings/) - Accessibility requirements

### Secondary (MEDIUM confidence)
- [markdown-link-check](https://github.com/tcort/markdown-link-check) - Link validation tool
- [Markdown Guide Basic Syntax](https://www.markdownguide.org/basic-syntax/) - Markdown best practices
- [Write the Docs Style Guides](https://www.writethedocs.org/guide/writing/style-guides/) - Technical writing patterns

### Tertiary (LOW confidence)
- Web search results for 2026 documentation practices (general trends, not verified)

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - Link validation tools are well-established
- Architecture patterns: HIGH - Patterns extracted directly from existing docs
- Pitfalls: HIGH - Common documentation issues, well-documented in industry
- Cross-reference audit: HIGH - Direct file examination

**Research date:** 2026-02-01
**Valid until:** 60+ days (documentation best practices are stable)
