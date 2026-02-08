# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-07)

**Core value:** Reliable async NMR predictions with full control over calculation parameters -- submit a molecule, get back accurate 1H/13C shifts without babysitting long-running calculations.
**Current focus:** v2.8 Expanded Solvent Support

## Current Position

Milestone: v2.8 Expanded Solvent Support
Phase: 55 of 58 (DELTA50 Benchmark Calculations)
Plan: 1 of 2 in current phase
Status: In progress (55-01 complete, 55-02 pending summary)
Last activity: 2026-02-08 — Completed 55-01-PLAN.md (Methanol + Water benchmarks)

Progress: [#####################.....] 116 plans complete across 11 milestones (v1.0-v2.7) + v2.8 Phases 54-55

## Performance Metrics

**Velocity:**
- Total plans completed: 116 (v1.0: 16, v1.1: 21, v2.0: 18, v2.0.1: 3, v2.1: 17, v2.2: 10, v2.3: 3, v2.4: 8, v2.5: 4, v2.6: 4, v2.7: 9, v2.8: 2)
- Average duration: ~6.4 min (excluding benchmark compute time)
- Total execution time: ~731 min (~12.2 hours) + ~8.5 hours benchmark compute

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
| v2.4 Docker Deployment | 6 | 8 | ~2 hours | Shipped 2026-02-03 |
| v2.5 ARM64 Docker Support | 4 | 4 | ~1 day | Shipped 2026-02-04 |
| v2.6 GCP Spot Deployment | 5 | 4 | ~1 day | Shipped 2026-02-05 |
| v2.7 Automated GCP Deployment | 5 | 9 | ~32 min | Shipped 2026-02-06 |
| v2.8 Expanded Solvent Support | 5 | 2 | ~8.5h compute | In progress |

## Accumulated Context

### Decisions

All prior decisions logged in PROJECT.md Key Decisions table.
v2.7 decisions archived to .planning/milestones/v2.7-ROADMAP.md.

| Decision | Phase | Rationale |
|----------|-------|-----------|
| Only benzene needed in input_gen SUPPORTED_SOLVENTS | 54-01 | Water, acetone, methanol already supported at NWChem layer |
| Title-case solvent names in runner, lowercase for NWChem | 54-01 | Matches existing CHCl3/DMSO convention; runner.lower() handles conversion |
| Added 'direct' to DFT shielding input for CPHF convergence | 55-01 | Required for NMR property calcs with COSMO; without it all calcs fail |
| Set NWChem cwd to scratch dir for file isolation | 55-01 | Prevents molecule.movecs cross-contamination between sequential runs |

### Key Context for v2.8

- NWChem COSMO now knows water, methanol, acetone, benzene by name. All 7 solvents supported.
- Benchmark CLI accepts all 6 solvents: CHCl3, DMSO, Methanol, Water, Acetone, Benzene.
- Benchmark is compute-intensive: 50 molecules x 4 solvents = 200 NWChem calculations (hours of compute).
- All 4 solvents use same experimental shifts from CDCl3 (Grimblat et al. 2023).
- Only B3LYP functional needed (WP04 out of scope).
- Pipeline per solvent: extend CLI -> run benchmark -> analyze -> copy factors -> add to solvents.py/shifts.py.
- NMR shielding + COSMO requires 'direct' in DFT block for CPHF convergence.
- Methanol + Water benchmarks: 100/100 complete, 0 failures. 342 H + 221 C data points per solvent.
- Acetone + Benzene benchmarks also completed concurrently (200 total calculations done).

### Pending Todos

- Dark mode (deferred from v2.1)
- Enhanced interactivity (deferred from v2.1)
- User accounts and calculation history
- NMReData enhancements: INCHI tag, 1D pseudo-spectrum tags, J-coupling
- Phase 48.1 implementation (machine info display - v2.6)

### Blockers/Concerns

**Active:**
None

## Session Continuity

Last session: 2026-02-08T18:32:00Z
Stopped at: Completed 55-01-PLAN.md (Methanol + Water benchmarks)
Resume file: None
Next: Complete 55-02-PLAN.md (Acetone + Benzene verification/summary), then Phase 56
Tests: 430 tests (15 new solvent tests added in Phase 54)
Codebase: ~7,300 LOC Python, ~3,050 LOC tests, ~950 LOC templates, ~2,400 LOC CSS, ~4,800 LOC docs, ~2,170 LOC GCP scripts
