---
phase: 39-cicd-ghcr-publishing
plan: 01
subsystem: infra
tags: [github-actions, ghcr, docker, ci-cd, multi-arch]

# Dependency graph
requires:
  - phase: 35-docker-worker-image
    provides: Dockerfile.worker for containerized worker
  - phase: 36-docker-api-image
    provides: Dockerfile.api for containerized API
provides:
  - GitHub Actions workflow for automated Docker image publishing
  - OCI labels for GHCR repository linking
  - Pre-built image documentation for deployment
affects: [40-deployment-polish, future-releases]

# Tech tracking
tech-stack:
  added: [docker/build-push-action@v6, docker/metadata-action@v5, docker/login-action@v3, docker/setup-qemu-action@v3, docker/setup-buildx-action@v3]
  patterns: [release-tag-triggered-ci, matrix-strategy-multi-image, gha-cache-scoping]

key-files:
  created: [.github/workflows/publish-images.yml]
  modified: [Dockerfile.api, Dockerfile.worker, docker-compose.yml]

key-decisions:
  - "Worker image amd64-only due to CREST/xTB lacking arm64 binaries"
  - "API image supports amd64+arm64 for broader deployment options"
  - "GITHUB_TOKEN authentication instead of PAT for security"
  - "GHA cache with per-image scope to avoid cache eviction"

patterns-established:
  - "Release-tag triggers: v* pattern for semantic versioning"
  - "Matrix strategy: single workflow for multiple images"
  - "OCI labels: source/description/licenses for GHCR linking"

# Metrics
duration: 8min
completed: 2026-02-03
---

# Phase 39 Plan 01: CI/CD + GHCR Publishing Summary

**GitHub Actions workflow for automated Docker image publishing to GHCR on release tags with multi-arch support for API (amd64+arm64) and amd64-only for worker**

## Performance

- **Duration:** 8 min
- **Started:** 2026-02-03T20:10:56Z
- **Completed:** 2026-02-03T20:19:00Z
- **Tasks:** 3
- **Files modified:** 4

## Accomplishments
- GitHub Actions workflow with matrix strategy for dual image builds
- Automatic semver tagging (v2.4.0, v2.4, latest) via metadata-action
- OCI labels for automatic GHCR-to-repository linking
- Documentation for pre-built image deployment option

## Task Commits

Each task was committed atomically:

1. **Task 1: Create GitHub Actions workflow** - `96fc6fa` (feat)
2. **Task 2: Add OCI labels to Dockerfiles** - `9c0dc0f` (feat)
3. **Task 3: Document pre-built image usage** - `da101ab` (docs)

## Files Created/Modified
- `.github/workflows/publish-images.yml` - CI/CD workflow for GHCR publishing (new)
- `Dockerfile.api` - Added OCI labels for GHCR linking
- `Dockerfile.worker` - Added OCI labels for GHCR linking
- `docker-compose.yml` - Added pre-built image deployment documentation

## Decisions Made
- **Worker amd64-only:** CREST and xTB only provide x86_64 pre-built binaries; arm64 would require conda-forge or compilation
- **API multi-arch:** No architecture-specific binaries, QEMU build time acceptable for 733MB image
- **GITHUB_TOKEN over PAT:** Built-in, auto-rotates, no manual secret management required
- **GHA cache with scope:** Separate cache per image (scope=api, scope=worker) to prevent cache thrashing

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None

## User Setup Required

None - workflow uses GITHUB_TOKEN which is automatically available to all GitHub Actions.

To publish images:
1. Create a release tag: `git tag v2.4.0 && git push origin v2.4.0`
2. GitHub Actions will automatically build and push to GHCR
3. Images available at:
   - `ghcr.io/steinbeck/qm-nmr-calc-api:v2.4.0`
   - `ghcr.io/steinbeck/qm-nmr-calc-worker:v2.4.0`

## Next Phase Readiness
- CI/CD infrastructure complete
- Ready for Phase 40: Deployment polish and documentation
- First actual publication will occur when v2.4.0 tag is pushed

---
*Phase: 39-cicd-ghcr-publishing*
*Completed: 2026-02-03*
