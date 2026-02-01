---
phase: 27-usage-guide
plan: 01
subsystem: documentation
tags: [usage, web-ui, conformers, solvents, presets, workflow]

# Dependency graph
requires:
  - phase: 26-installation-guide
    provides: Installation documentation with validation steps
provides:
  - Web UI workflow documentation (submit, status, results pages)
  - Calculation modes guide (single vs ensemble with decision flowchart)
  - Solvent selection documentation (CHCl3, DMSO, vacuum with COSMO)
  - Preset comparison (draft vs production with accuracy metrics)
  - Tips for best results section
affects: [27-02-usage-guide, 28-architecture-docs, 30-science-methodology]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - Usage documentation structure: Introduction → Web UI Workflow → Concepts → Tips → Related Docs

key-files:
  created: []
  modified:
    - docs/usage.md

key-decisions:
  - "Document web UI as primary interface, REST API in separate plan"
  - "Include decision flowchart for calculation mode selection"
  - "Provide accuracy comparison table between presets"
  - "Cross-reference installation.md as prerequisite"

patterns-established:
  - "Feature descriptions match actual UI field names from templates"
  - "Tables used for settings comparison (options, accuracy, timing)"
  - "Related Documentation section links to all relevant docs"

# Metrics
duration: 2min
completed: 2026-02-01
---

# Phase 27 Plan 01: Web UI Workflow and Calculation Concepts Summary

**368-line usage guide documenting web UI workflow (submit/status/results pages), calculation modes (single vs ensemble with decision guidance), solvent selection (COSMO solvation), and preset tradeoffs (draft vs production with accuracy metrics)**

## Performance

- **Duration:** 2 minutes
- **Started:** 2026-02-01T07:31:39Z
- **Completed:** 2026-02-01T07:33:19Z
- **Tasks:** 2
- **Files modified:** 1

## Accomplishments

- Complete web UI workflow documentation covering all three pages (submit, status, results)
- Calculation modes section with single vs ensemble comparison and decision flowchart
- Solvent selection guide explaining COSMO solvation and when to use each solvent
- Preset comparison with accuracy metrics (MAE values) and timing estimates
- Tips for best results section with practical guidance

## Task Commits

Each task was committed atomically:

1. **Task 1: Web UI Workflow and Calculation Settings Documentation** - `9019b4f` (docs)
   - Web UI workflow: submit, status, results pages
   - Calculation modes: single vs ensemble with decision guide
   - Solvent selection: CHCl3, DMSO, vacuum with COSMO explanation
   - Presets: draft vs production with accuracy comparison

**Plan metadata:** (to be committed)

## Files Created/Modified

- `docs/usage.md` - Comprehensive usage guide (368 lines)
  - Introduction with target audience and prerequisites
  - Web UI Workflow (submit, status, results pages)
  - Calculation Modes (single vs ensemble with decision flowchart)
  - Solvent Selection (COSMO solvation model)
  - Calculation Presets (draft vs production comparison)
  - Tips for Best Results section
  - Related Documentation cross-references

## Decisions Made

**Documentation scope:**
- Focus on web UI workflow in this plan, REST API in 27-02
- Include decision flowchart for calculation mode selection
- Provide accuracy metrics from presets.py (MAE estimates)

**Content from source files:**
- Verified UI elements against submit.html, status.html, results.html
- Solvent descriptions match solvents.py SUPPORTED_SOLVENTS
- Preset configurations match presets.py PRESETS dictionary
- Step names in status page match actual task steps

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None - documentation created based on source code analysis.

## User Setup Required

None - documentation only, no external service configuration.

## Next Phase Readiness

**Ready for Plan 27-02 (REST API Reference):**
- Web UI workflow complete
- Calculation concepts documented (modes, solvents, presets)
- REST API endpoints referenced but not detailed
- Result interpretation guide pending

**Context for API documentation:**
- Same calculation concepts apply to API
- Web UI workflow shows expected user journey
- Endpoint URLs visible in template code

**No blockers or concerns.**

---
*Phase: 27-usage-guide*
*Completed: 2026-02-01*
