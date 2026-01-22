# Phase 9: Benchmark Calculations - Context

**Gathered:** 2026-01-22
**Status:** Ready for planning

<domain>
## Phase Boundary

Execute 200 QM calculations (50 DELTA50 molecules × 2 functionals × 2 solvents) to generate raw shielding data for scaling factor derivation in Phase 10. The benchmark infrastructure from Phase 8 is already built.

**Calculation Matrix:**
| Functional | Solvent | Count | Purpose |
|------------|---------|-------|---------|
| B3LYP | CHCl3 | 50 | General-purpose baseline (1H + 13C) |
| B3LYP | DMSO | 50 | General-purpose in alternate solvent |
| WP04 | CHCl3 | 50 | 1H-optimized functional |
| WP04 | DMSO | 50 | 1H-optimized in alternate solvent |

</domain>

<decisions>
## Implementation Decisions

### Execution Sequence
- Run by molecule first: all 4 conditions for compound_01, then compound_02, etc.
- No priority between conditions — run all 4 (B3LYP/CHCl3, B3LYP/DMSO, WP04/CHCl3, WP04/DMSO) equally
- **Pilot first:** Run 5 molecules (20 calculations) to validate pipeline before full run
- After pilot validation, continue with remaining 45 molecules

### Failure Handling
- Continue and log on failure — don't stop the entire run
- **10% threshold:** If >5 of 50 molecules fail, pause to investigate
- No auto-retry — document failures for analysis (may indicate molecule-specific issues)
- Preserve full NWChem .out file + error message + timing for all failures

### Progress Monitoring
- **Status file as interface:** Write status.json that can be checked from any session
- CLI command to read and display status (human-friendly output)
- Show: completed/total, current molecule, estimated time remaining
- No push notifications — check status manually
- Generate BENCHMARK-RESULTS.md summary report on completion

### Resource Management
- Run on dedicated/idle machine — use all available resources
- **Sequential execution:** One calculation at a time (NWChem uses all cores internally)
- No time constraints — start and let it run to completion
- **Graceful pause:** Support clean stop after current calculation completes

### Headless Operation
- **Run as detached process:** Must survive terminal/session disconnect (nohup, screen, tmux, or similar)
- Status file is the monitoring interface — independent of running session
- Graceful stop via signal or marker file that runner checks between calculations

### Claude's Discretion
- Exact status.json schema
- Choice of detachment mechanism (nohup vs screen vs tmux)
- Pilot molecule selection (first 5 or representative sample)

</decisions>

<specifics>
## Specific Ideas

- Status file should be readable by both humans (CLI command) and machines (JSON)
- Pilot run validates the full pipeline before committing to 200 calculations
- User will monitor from Claude Code session but runner must continue independently

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 09-benchmark-calculations*
*Context gathered: 2026-01-22*
