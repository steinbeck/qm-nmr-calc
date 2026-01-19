# Phase 1: Foundation - Context

**Gathered:** 2026-01-19
**Status:** Ready for planning

<domain>
## Phase Boundary

Establish working calculation infrastructure with job queue and ISiCLE integration. This phase delivers: job queuing and execution, ISiCLE/NWChem geometry optimization on a test molecule, clear error handling for failed calculations, and job state persistence across restarts.

</domain>

<decisions>
## Implementation Decisions

### Job queue design
- Store comprehensive metadata: status, timestamps, input parameters, molecule info, resource usage (memory, CPU time), and detailed calculation logs
- One directory per job structure: `jobs/abc123/` with status.json, input files, and outputs — self-contained and easy to inspect
- Interrupted jobs (running at process restart) marked as failed with clear reason ("interrupted - process restart") — no automatic resume, user must decide whether to resubmit
- No job timeout — let calculations run until natural completion or failure

### ISiCLE integration
- Thin wrapper approach — pass parameters through to ISiCLE, don't over-abstract
- Fail fast at startup — validate NWChem binary exists and works when service starts
- Store ISiCLE + NWChem versions with each job for reproducibility

### Error handling
- Failed jobs preserve all partial outputs: geometry steps, logs, temp files — everything needed for debugging
- No automatic retries — fail immediately, user decides whether to resubmit
- Error messages show both: user-friendly summary + option to see full technical details (NWChem error codes, file paths, etc.)

### Project structure
- Top-level package: `qm_nmr_calc`
- Job data lives inside project: `./data/jobs/`
- Use modern Python tooling: uv + pyproject.toml

### Claude's Discretion
- NWChem configuration approach (env vars, config file, or both)
- Module organization (flat vs domain-based subpackages) — decide based on expected complexity
- Error categorization scheme (input error vs convergence failure vs system error, etc.)

</decisions>

<specifics>
## Specific Ideas

- Jobs that were running when process died should fail with clear message explaining the interruption — user needs to know why and decide whether to resubmit
- "Fail for a reason" philosophy — automatic resume doesn't make sense for QM calculations that fail, they need debugging

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 01-foundation*
*Context gathered: 2026-01-19*
