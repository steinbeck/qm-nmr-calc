---
phase: 26-installation-guide
plan: 01
subsystem: documentation
tags: [installation, nwchem, uv, crest, xtb, conda, troubleshooting]

# Dependency graph
requires:
  - phase: 25-readme-and-docs-structure
    provides: Documentation directory structure and README with Quick Start
provides:
  - Comprehensive installation documentation (664 lines)
  - System dependencies guide (NWChem, OpenMPI, Python 3.11+)
  - uv package manager setup instructions
  - Optional CREST/xTB configuration for ensemble mode
  - Environment validation steps with health endpoint
  - Troubleshooting section with 8 common issues
affects: [27-usage-guide, 28-architecture-docs, 29-library-docs, 30-science-methodology, 31-contributor-guide]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - Installation documentation structure: Prerequisites → System Deps → Project Setup → Optional Tools → Validation → Troubleshooting → Help

key-files:
  created: []
  modified:
    - docs/installation.md

key-decisions:
  - "Installation guide targets both academic researchers and developers"
  - "Document both conda and manual installation for CREST/xTB"
  - "Include health endpoint validation as primary environment check"
  - "Troubleshooting covers 8 most common installation issues"
  - "Optional tools section clearly states app works without CREST/xTB"

patterns-established:
  - "Documentation includes clear 'Why this matters' explanations for dependencies"
  - "Code blocks are copy-paste ready with expected output"
  - "Use markdown admonitions (> Note:) for important callouts"
  - "Cross-reference related docs (usage.md, architecture.md, science.md)"

# Metrics
duration: 3min
completed: 2026-01-31
---

# Phase 26 Plan 01: Installation Guide Summary

**Comprehensive 664-line installation guide covering NWChem/OpenMPI setup, uv package manager, optional CREST/xTB ensemble tools, environment validation via health endpoint, and troubleshooting for 8 common issues**

## Performance

- **Duration:** 3 minutes
- **Started:** 2026-01-31T22:28:06Z
- **Completed:** 2026-01-31T22:31:01Z
- **Tasks:** 2
- **Files modified:** 1

## Accomplishments

- Complete installation documentation (664 lines) covering all INSTALL-01 through INSTALL-05 requirements
- System dependencies with multi-distro instructions (Ubuntu/Debian, Fedora/RHEL)
- Clear explanation of why each dependency matters (NWChem for DFT, OpenMPI for parallelization, Python 3.11+)
- Optional CREST/xTB setup with both conda (recommended) and manual installation paths
- Step-by-step environment validation using /api/v1/health/ready endpoint
- Comprehensive troubleshooting section covering 8 common installation issues
- Cross-references to usage.md, architecture.md, and science.md for next steps

## Task Commits

Each task was committed atomically:

1. **Task 1: Core Installation Documentation** - `a2b1f25` (docs)
   - System dependencies: NWChem, OpenMPI, Python 3.11+
   - Package manager: uv installation and setup
   - Project setup: git clone and uv sync workflow
   - Requirements INSTALL-01 and INSTALL-02 satisfied

2. **Task 2: Optional Tools, Validation, and Troubleshooting** - `9fc7bb4` (docs)
   - CREST/xTB setup: conda and manual installation
   - Environment validation: step-by-step health checks
   - Troubleshooting: 8 common issues with detailed fixes
   - Getting Help: health diagnostics and issue reporting
   - Requirements INSTALL-03, INSTALL-04, and INSTALL-05 satisfied

**Plan metadata:** (to be committed)

## Files Created/Modified

- `docs/installation.md` - Comprehensive installation guide (664 lines)
  - Prerequisites Overview
  - System Dependencies (NWChem, OpenMPI, Python)
  - Project Setup with uv
  - Optional: CREST and xTB Setup
  - Environment Validation
  - Troubleshooting Common Issues
  - Getting Help

## Decisions Made

**Installation guide structure:**
- Target both academic researchers and developers
- Provide multi-distro instructions (Ubuntu/Debian + Fedora/RHEL)
- Include both recommended (conda) and manual paths for CREST/xTB
- Use health endpoint as primary validation method

**Optional tools documentation:**
- Clearly state app works fully without CREST/xTB (uses RDKit instead)
- Document detect_crest_available() pattern (both binaries required)
- Explain GFN2-xTB force field usage in ensemble generation

**Troubleshooting approach:**
- Problem/Symptom/Fix format for each issue
- Cover most common scenarios based on dependency stack
- Include diagnostic commands (health endpoint, which, versions)

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None - documentation created based on existing codebase patterns and implementation.

## User Setup Required

None - documentation only, no external service configuration.

## Next Phase Readiness

**Ready for Phase 27 (Usage Guide):**
- Installation guide complete with clear next steps
- Cross-references to usage.md established
- Health endpoint documented as validation method
- Troubleshooting covers installation-specific issues

**Context for future documentation phases:**
- Installation patterns established (Prerequisites → Setup → Validation → Troubleshooting)
- Health endpoint serves as diagnostic tool for all docs
- Cross-reference pattern established for related docs

**No blockers or concerns.**

---
*Phase: 26-installation-guide*
*Completed: 2026-01-31*
