# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-02)

**Core value:** Reliable async NMR predictions with full control over calculation parameters -- submit a molecule, get back accurate 1H/13C shifts without babysitting long-running calculations.
**Current focus:** v2.4 Docker Deployment - Phase 40 Documentation

## Current Position

Milestone: v2.4 Docker Deployment
Phase: 40 of 40 (Documentation)
Plan: 01 of 02 complete
Status: In progress
Last activity: 2026-02-03 -- Completed 40-01-PLAN.md (Docker quick start in README)

Progress: [####################] 100% (v1.0-v2.3) | [##################..] 90% (v2.4)

## Performance Metrics

**Velocity:**
- Total plans completed: 95 (v1.0: 16, v1.1: 21, v2.0: 18, v2.0.1: 3, v2.1: 17, v2.2: 10, v2.3: 3, v2.4: 7)
- Average duration: ~7 min
- Total execution time: ~661 min (~11 hours)

**By Milestone:**

| Milestone | Phases | Plans | Duration | Status |
|-----------|--------|-------|----------|--------|
| v1.0 Core NMR Service | 6 | 16 | 2 days | Shipped 2026-01-20 |
| v1.1 Accurate Chemical Shifts | 8 | 21 | 5 days | Shipped 2026-01-25 |
| v2.0 Conformational Sampling | 6 | 18 | ~2 days | Shipped 2026-01-28 |
| v2.0.1 Conformer Pre-selection | 1 | 3 | 18 min | Shipped 2026-01-30 |
| v2.1 UI Redesign | 6 | 17 | ~3 days | Shipped 2026-01-31 |
| v2.2 Documentation | 7 | 10 | 2 days | Shipped 2026-02-01 |
| v2.3 NMReData Export | 3 | 3 | 1 day | Shipped 2026-02-01 |
| v2.4 Docker Deployment | 6 | TBD | - | In progress |

## Accumulated Context

### Decisions

All prior decisions logged in PROJECT.md Key Decisions table.
v2.3 decisions archived to milestones/v2.3-ROADMAP.md.

**v2.4 Decisions:**
- Worker image amd64-only due to CREST/xTB lacking arm64 binaries
- API image supports amd64+arm64 for broader deployment options
- GITHUB_TOKEN authentication for GHCR (not PAT)
- GHA cache with per-image scope to avoid cache eviction
- Docker as primary deployment method in README (not source installation)
- Pre-built GHCR images referenced as default (not build-from-source)

### Roadmap Evolution

- v1.0: 6 phases (1-6), shipped 2026-01-20
- v1.1: 8 phases (7-11.2, including 3 inserted), shipped 2026-01-25
- v2.0: 6 phases (12-17), shipped 2026-01-28
- v2.0.1: 1 phase (24), shipped 2026-01-30
- v2.1: 6 phases (18-23), shipped 2026-01-31
- v2.2: 7 phases (25-31), shipped 2026-02-01
- v2.3: 3 phases (32-34), shipped 2026-02-01
- v2.4: 6 phases (35-40), in progress

### Pending Todos

- Dark mode (deferred from v2.1)
- Enhanced interactivity (deferred from v2.1)
- User accounts and calculation history
- NMReData enhancements: INCHI tag, 1D pseudo-spectrum tags, J-coupling (requires additional QM)

### Blockers/Concerns

**Active:**
None

**Research Flags (from v2.4 research):**
- Phase 35: Resolved - Used Miniconda Python 3.11 for glibc compatibility, MPI configured with OMPI_ALLOW_RUN_AS_ROOT
- Phase 36: Resolved - Added X11 libraries (libxrender1, libxext6, libexpat1) for RDKit drawing
- Phase 37: Resolved - SIGINT for Huey graceful shutdown, 5-min grace period, 512MB shm_size for MPI
- Phase 39: Resolved - Worker amd64-only, API multi-arch; GITHUB_TOKEN for auth

## Session Continuity

Last session: 2026-02-03 21:51 UTC
Stopped at: Completed 40-01-PLAN.md (Docker quick start in README)
Resume file: None
Next: 40-02-PLAN.md (deployment guide)
Tests: All tests passing (356 tests)
Codebase: ~6,400 LOC Python, ~2,450 LOC tests, ~950 LOC templates, ~2,400 LOC CSS, ~4,100 LOC docs
Docker: Worker image 2.1GB, API image ~733MB, Caddy reverse proxy, GHCR publishing workflow
