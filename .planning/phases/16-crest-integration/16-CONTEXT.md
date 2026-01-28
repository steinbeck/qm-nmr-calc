# Phase 16: CREST Integration (Optional High-Accuracy Mode) - Context

**Gathered:** 2026-01-28
**Status:** Ready for planning

<domain>
## Phase Boundary

Add CREST/xTB as an optional conformer generation backend. Auto-detect availability at startup, run CREST with ALPB solvation when requested and available, handle hung processes with a timeout, and feed CREST conformers into the existing DFT+NMR pipeline. The app must continue to work without CREST installed (RDKit-only mode validated in Phases 13-15).

</domain>

<decisions>
## Implementation Decisions

### Detection & availability
- Require both `crest` and `xtb` binaries on PATH -- partial install (xTB only) means CREST mode unavailable
- Detection happens once at startup, result cached -- no per-request PATH scanning
- Health endpoint gets a simple boolean field: `crest_available: true/false` (no version info)
- If user explicitly requests CREST mode but it's not installed, **fail the job immediately** with a clear error message -- no silent fallback

### Timeout & fallback policy
- No timeouts on DFT/NMR calculations -- those run as long as they need
- CREST-specific hang timeout (generous, e.g. 2-4 hours) to catch infinite loops on macrocycles
- When CREST times out: **fail the job** with a message suggesting RDKit mode -- no auto-fallback
- Timeout is fixed server-side (config/env var) -- not user-configurable in API

### Solvation mapping
- Hardcode two ALPB mappings: CHCl3 -> 'chcl3', DMSO -> 'dmso'
- Vacuum/gas-phase jobs: **skip CREST entirely**, force RDKit mode -- CREST solvation is the benefit
- Trust standard ALPB keyword mapping, validate in integration tests
- If a future solvent is unsupported by ALPB: force RDKit mode for that solvent (CREST unavailable)

### CREST output handling
- Skip RMSD deduplication (CREST does this internally), but apply pre-DFT energy window filter
- Store CREST energies in Hartree (matches existing DFT energy convention and ConformerData.energy_unit)
- Split CREST multi-structure XYZ into individual per-conformer XYZ files (matches RDKit path convention)
- Track conformer_mode in job parameters only -- no duplication in ensemble metadata

### Claude's Discretion
- CREST command-line argument construction
- Environment variable setup (OMP_STACKSIZE, GFORTRAN_UNBUFFERED_ALL)
- Exact timeout value within 2-4 hour range
- Multi-XYZ parsing implementation details
- CREST output directory structure and cleanup

</decisions>

<specifics>
## Specific Ideas

- User emphasized: NMR GIAO calcs take as long as they take -- timeouts are only for detecting hung CREST processes, never for legitimate calculation time
- Fail-fast philosophy: explicit CREST request without CREST installed = immediate error, not degraded service
- Vacuum jobs don't benefit from CREST enough to justify running it without solvation

</specifics>

<deferred>
## Deferred Ideas

None -- discussion stayed within phase scope

</deferred>

---

*Phase: 16-crest-integration*
*Context gathered: 2026-01-28*
