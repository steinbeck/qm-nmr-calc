# Phase 4: Results Delivery - Context

**Gathered:** 2026-01-19
**Status:** Ready for planning

<domain>
## Phase Boundary

Users can retrieve all calculation results in multiple formats. This includes JSON with atom-assigned shifts, downloadable optimized geometry files (XYZ/SDF), raw NWChem output files, and opt-in email notifications when jobs complete.

</domain>

<decisions>
## Implementation Decisions

### Claude's Discretion

User has delegated all implementation decisions to Claude. Apply these principles:

**JSON Results:**
- Follow existing API patterns (JobStatusResponse structure)
- Include shifts grouped by nucleus (1H, 13C) with atom indices
- Include calculation metadata (functional, basis set, solvent) for reproducibility

**Download Endpoints:**
- RESTful URL patterns consistent with existing /api/v1/jobs/{job_id} structure
- Separate endpoints for different file types (geometry, raw output)
- Appropriate content-type headers and filename suggestions

**Email Notifications:**
- Opt-in field at submission time (email address)
- Send on job completion (success or failure)
- Include job ID, status, and link to results
- Simple, functional email format

**Error Handling:**
- Consistent with existing RFC 7807 ProblemDetail pattern
- Clear errors for: job not found, job not complete, file not available
- Appropriate HTTP status codes (404 for not found, 409 for not ready)

</decisions>

<specifics>
## Specific Ideas

No specific requirements — open to standard approaches that follow existing project patterns.

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope.

</deferred>

---

*Phase: 04-results-delivery*
*Context gathered: 2026-01-19*
