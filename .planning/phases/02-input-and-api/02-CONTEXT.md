# Phase 2: Input and API - Context

**Gathered:** 2026-01-19
**Status:** Ready for planning

<domain>
## Phase Boundary

REST API for molecule submission and job status polling. Users can submit molecules via SMILES string or MOL/SDF file upload, receive a job ID immediately, and poll for job status. OpenAPI/Swagger documentation served at /docs.

</domain>

<decisions>
## Implementation Decisions

### Molecule Input Formats
- Strict SMILES validation — reject invalid SMILES immediately with specific error
- Accept MOL and SDF file uploads (single molecule only per file)
- Reject multi-molecule SDF files with clear error
- Detailed validation error messages — include specific error location and fix suggestion if possible
- Optional molecule name/label field for user identification (e.g., "Aspirin")

### API Response Design
- Job ID format: 12-character hex string (URL-safe, compact)
- Successful submission returns full JobStatus object (not just job ID)
- HTTP 202 Accepted for successful job submission (async processing)
- Basic 4 job statuses: queued, running, complete, failed
- No progress percentage for running jobs (QM calculations hard to estimate)

### Claude's Discretion
- API error format (RFC 7807 vs simple JSON)
- API versioning (/api/v1 prefix or not)
- Endpoint naming style (resource-based vs action-based)
- File upload endpoint design (same endpoint with content-type or separate)
- Service health endpoint (/health)
- Poll interval hints (retry-after header)
- 404 vs 200 response for non-existent job IDs

</decisions>

<specifics>
## Specific Ideas

No specific requirements — open to standard FastAPI/REST approaches.

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope.

</deferred>

---

*Phase: 02-input-and-api*
*Context gathered: 2026-01-19*
