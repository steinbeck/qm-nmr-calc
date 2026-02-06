---
phase: 49-config-foundation-and-pricing-query
plan: 02
subsystem: infra
tags: [httpx, cloudprice, spot-pricing, gcp, caching, bash-library]

# Dependency graph
requires:
  - phase: 49-01
    provides: gcp/__init__.py package marker for Python imports
provides:
  - GCP spot pricing query module with CloudPrice.net API integration
  - Hardcoded regional fallback for offline/API-failure scenarios
  - File-based pricing cache with 24-hour TTL
  - Bash pricing library (get_cheapest_region, get_cheapest_zone)
  - CLI script for pricing queries with JSON output
affects: [50-machine-type-selection, 51-deployment-orchestration]

# Tech tracking
tech-stack:
  added: []
  patterns: [three-tier-resolution (cache -> API -> fallback), file-based-ttl-cache, python-bash-bridge-via-json]

key-files:
  created:
    - gcp/query_pricing.py
    - gcp/lib/pricing.sh
    - tests/test_gcp_pricing.py
  modified: []

key-decisions:
  - "CloudPrice.net API queried best-effort with hardcoded fallback always available (LOW API confidence)"
  - "Cache stores API results only, not fallback data (avoids caching stale fallback as 'fresh')"
  - "10 fallback regions across 4 continents ranked by historical spot pricing patterns"

patterns-established:
  - "Three-tier resolution: cache -> API -> hardcoded fallback (never fails)"
  - "Bash library sources Python JSON output via jq for shell integration"
  - "File-based TTL cache: JSON with cached_at timestamp, 24h expiry"

# Metrics
duration: 8min
completed: 2026-02-06
---

# Phase 49 Plan 02: Spot Pricing Query Summary

**CloudPrice.net API integration with httpx, 24h file-based TTL cache, and 10-region hardcoded fallback across 4 continents**

## Performance

- **Duration:** ~8 min
- **Started:** 2026-02-06T10:12:07Z
- **Completed:** 2026-02-06T14:04:29Z
- **Tasks:** 2 (TDD: RED + GREEN)
- **Files created:** 3

## Accomplishments
- Pricing module with three-tier resolution: cache -> CloudPrice.net API -> hardcoded fallback
- 19 tests covering fallback structure, cache TTL, mocked API responses, and integration cascading
- Bash library providing get_cheapest_region(), get_cheapest_zone(), get_pricing_table() for shell scripts
- CLI with --cpu-cores, --ram-gb, --cache-dir, --refresh flags, JSON output parseable by jq

## Task Commits

Each task was committed atomically:

1. **Task 1: Write failing tests for pricing module** - `af25673` (test) - RED phase, 19 tests failing with ImportError
2. **Task 2: Implement pricing module and make tests pass** - `bc44adf` (feat) - GREEN phase, all 19 tests passing

## Files Created/Modified
- `gcp/query_pricing.py` - Python pricing module: API query, cache, fallback, CLI
- `gcp/lib/pricing.sh` - Bash library wrapping Python output for shell consumption
- `tests/test_gcp_pricing.py` - 19 tests: fallback (6), cache (5), API mocked (4), integration (4)

## Decisions Made
- CloudPrice.net API queried best-effort with httpx (10s timeout); hardcoded fallback always available since API reliability is LOW confidence
- Cache stores only API results, not fallback data -- avoids treating stale fallback as "fresh" cached data
- 10 fallback regions covering US (5), Europe (2), Asia (2), South America (1), ranked by historical spot pricing
- API key optional via CLOUDPRICE_API_KEY env var -- API may work without auth for basic queries
- get_ranked_regions() never fails -- always returns a non-empty list

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
- Pre-existing test failures in test_conformer_nwchem.py (ValueError unpacking) and test_nwchem_integration.py (NWChem MPI errors) -- unrelated to this plan, not caused by pricing changes
- CloudPrice.net API returns connection error in current environment (expected per LOW confidence rating in research) -- fallback works correctly

## User Setup Required

None - no external service configuration required. CLOUDPRICE_API_KEY environment variable is optional.

## Next Phase Readiness
- Pricing data available for Phase 50 (Machine Type Selection) to use in region selection
- Bash library ready for Phase 51 (Deployment Orchestration) to source
- Cache directory (gcp/.cache/) auto-created on first API success

---
*Phase: 49-config-foundation-and-pricing-query*
*Completed: 2026-02-06*
