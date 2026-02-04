# Phase 43: CI/CD Integration - Research

**Researched:** 2026-02-04
**Domain:** GitHub Actions native ARM64 runners, Docker multi-arch manifests, GHCR publishing
**Confidence:** HIGH

## Summary

This phase updates the existing CI/CD workflow to build ARM64 worker images using GitHub's native `ubuntu-24.04-arm` runner and create multi-architecture manifest lists. The key architecture decision is whether to use parallel native runners (fast) or QEMU emulation (slow). Since QEMU builds of the 2.1GB worker image would timeout or take 30+ minutes, the approach MUST use separate native runners with digest-based artifact passing.

The workflow pattern involves three jobs: (1) build amd64 image on `ubuntu-latest`, (2) build arm64 image on `ubuntu-24.04-arm`, (3) merge job that downloads digests and creates multi-arch manifest using `docker buildx imagetools create`. This is the Docker-documented best practice for native multi-platform builds.

**Primary recommendation:** Use the matrix-based native runner pattern from Docker's official documentation and the [sredevopsorg/multi-arch-docker-github-workflow](https://github.com/sredevopsorg/multi-arch-docker-github-workflow) reference implementation. Build platform-specific images by digest, upload digests as artifacts, then merge into a manifest list.

## Standard Stack

The established tools for multi-arch CI/CD with GitHub Actions:

### Core Actions
| Action | Version | Purpose | Why Standard |
|--------|---------|---------|--------------|
| `docker/setup-buildx-action` | v3.12.0 | Enable buildx builder | Required for multi-platform builds and imagetools |
| `docker/login-action` | v3.7.0 | Registry authentication | Official, supports GHCR with GITHUB_TOKEN |
| `docker/build-push-action` | v6.18.0 | Build and push images | Most maintained, supports push-by-digest |
| `docker/metadata-action` | v5.10.0 | Generate tags/labels | Automatic semver handling |
| `actions/upload-artifact` | v4 (v6.0.0 latest) | Share digests between jobs | Required for parallel builds |
| `actions/download-artifact` | v4 (v7.0.0 latest) | Retrieve digests for merge | Required for manifest creation |

### GitHub-Hosted Runners
| Runner Label | Architecture | Use Case |
|--------------|--------------|----------|
| `ubuntu-latest` | x86_64 (amd64) | Standard amd64 builds |
| `ubuntu-24.04-arm` | aarch64 (arm64) | Native ARM64 builds |
| `ubuntu-22.04-arm` | aarch64 (arm64) | Alternative ARM64 runner |

### Supporting
| Tool | Purpose | When to Use |
|------|---------|-------------|
| `actions/checkout@v4` | Clone repository | Every workflow |
| `docker buildx imagetools create` | Create multi-arch manifest | Merge job |
| GITHUB_TOKEN | Authentication | Built-in, no PAT needed |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Native ARM runners | QEMU emulation | QEMU too slow for 2.1GB worker image (30+ min vs ~5 min) |
| Artifact-based digest passing | Registry cache with tags | Artifacts cleaner, no intermediate tags polluting registry |
| `int128/docker-manifest-create-action` | Raw `buildx imagetools` | Raw command more flexible, official pattern |

## Architecture Patterns

### Recommended Workflow Structure
```
jobs:
  build:          # Matrix job: runs on native runners per platform
    strategy:
      matrix:
        platform: [linux/amd64, linux/arm64]
    runs-on: [runner selected by platform]
    outputs: uploads digest artifact

  merge:          # Single job: creates multi-arch manifest
    needs: [build]
    runs-on: ubuntu-latest
    steps: downloads digests, creates manifest list
```

### Pattern 1: Dynamic Runner Selection
**What:** Select runner based on matrix platform value
**When to use:** Always - matches platform to native runner
**Example:**
```yaml
# Source: https://github.com/sredevopsorg/multi-arch-docker-github-workflow
runs-on: ${{ matrix.platform == 'linux/amd64' && 'ubuntu-latest' || matrix.platform == 'linux/arm64' && 'ubuntu-24.04-arm' }}
```

### Pattern 2: Push by Digest
**What:** Push platform-specific images by digest (not tag), export digest for later
**When to use:** Required for native multi-runner builds
**Example:**
```yaml
# Source: Docker multi-platform documentation
- name: Build and push by digest
  id: build
  uses: docker/build-push-action@v6
  with:
    context: .
    platforms: ${{ matrix.platform }}
    outputs: type=image,name=${{ env.GHCR_IMAGE }},push-by-digest=true,name-canonical=true,push=true
    cache-from: type=gha,scope=${{ github.repository }}-${{ matrix.platform }}
    cache-to: type=gha,scope=${{ github.repository }}-${{ matrix.platform }}

- name: Export digest
  run: |
    mkdir -p /tmp/digests
    digest="${{ steps.build.outputs.digest }}"
    touch "/tmp/digests/${digest#sha256:}"
```

### Pattern 3: Artifact-Based Digest Passing
**What:** Upload digest files as artifacts, download in merge job
**When to use:** Required for cross-job communication
**Example:**
```yaml
# Build job - upload
- uses: actions/upload-artifact@v4
  with:
    name: digests-${{ env.PLATFORM_PAIR }}  # e.g., digests-linux-amd64
    path: /tmp/digests/*
    retention-days: 1

# Merge job - download
- uses: actions/download-artifact@v4
  with:
    path: /tmp/digests
    pattern: digests-*
    merge-multiple: true
```

### Pattern 4: Manifest Creation with imagetools
**What:** Create multi-arch manifest from platform-specific digests
**When to use:** Final merge step
**Example:**
```yaml
# Source: Docker buildx documentation
- name: Create manifest list and push
  working-directory: /tmp/digests
  run: |
    docker buildx imagetools create \
      $(jq -cr '.tags | map("-t " + .) | join(" ")' <<< "$DOCKER_METADATA_OUTPUT_JSON") \
      $(printf '${{ env.GHCR_IMAGE }}@sha256:%s ' *)
```

### Anti-Patterns to Avoid
- **QEMU for large images:** Worker image is 2.1GB - QEMU builds will timeout or take 30+ minutes
- **Single workflow file with all platforms:** Forces QEMU; use matrix with native runners instead
- **Hardcoded digests:** Always use outputs from build step, never hardcode SHA values
- **Same artifact name for different platforms:** v4 artifacts are immutable - use platform-specific names
- **Using docker/setup-qemu-action for native builds:** Not needed when using native ARM runners

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Multi-arch manifest | Manual `docker manifest create` | `docker buildx imagetools create` | Handles attestations, annotations, OCI compliance |
| Cross-job data | Environment variables, cache | `actions/upload-artifact` + `actions/download-artifact` | Official pattern, immutable, reliable |
| Tag generation | Shell scripts parsing git refs | `docker/metadata-action` | Handles edge cases, OCI labels, multiple tag formats |
| Platform name formatting | Manual string replace | `${platform//\//-}` bash pattern | Standard pattern, used in reference workflows |
| Registry auth | Manual docker login | `docker/login-action` | Handles token refresh, supports multiple registries |

**Key insight:** The Docker Actions ecosystem handles the complexity of multi-arch builds. The pattern of push-by-digest + artifact upload + imagetools merge is explicitly documented and battle-tested.

## Common Pitfalls

### Pitfall 1: ARM64 Runner Availability
**What goes wrong:** Workflow fails with "no runners available" for `ubuntu-24.04-arm`
**Why it happens:** ARM64 runners are FREE for public repos only; private repos require Team/Enterprise plan
**How to avoid:**
- Ensure repository is public, OR
- Ensure organization has Team/Enterprise plan for private repos
**Warning signs:** Job queued indefinitely, "no matching runners" error

### Pitfall 2: Artifact Name Collision
**What goes wrong:** "Artifact with name 'digests' already exists" error
**Why it happens:** v4 artifacts are immutable; matrix jobs try to upload to same name
**How to avoid:** Use platform-specific artifact names:
```yaml
name: digests-${{ env.PLATFORM_PAIR }}  # digests-linux-amd64, digests-linux-arm64
```
**Warning signs:** "Artifact already exists" error after first platform builds

### Pitfall 3: Missing `merge-multiple: true`
**What goes wrong:** Merge job only sees digests from one platform
**Why it happens:** download-artifact@v4 creates subdirectories per artifact by default
**How to avoid:** Set `merge-multiple: true` to flatten all digests into single directory:
```yaml
- uses: actions/download-artifact@v4
  with:
    path: /tmp/digests
    pattern: digests-*
    merge-multiple: true
```
**Warning signs:** Only one platform in final manifest

### Pitfall 4: Cache Scope Conflicts
**What goes wrong:** Cache invalidation or collision between platforms
**Why it happens:** Default cache scope is global; both platforms overwrite each other
**How to avoid:** Use platform-specific cache scopes:
```yaml
cache-from: type=gha,scope=${{ github.repository }}-${{ matrix.platform }}
cache-to: type=gha,scope=${{ github.repository }}-${{ matrix.platform }}
```
**Warning signs:** "Cache not found" despite previous successful builds, inconsistent build times

### Pitfall 5: Worker vs API Different Patterns
**What goes wrong:** Trying to apply same multi-arch pattern to both images
**Why it happens:** API image can use QEMU (small, no compilation); worker requires native
**How to avoid:**
- API image: Keep existing single-job QEMU pattern (already works)
- Worker image: New multi-job native pattern (this phase)
**Warning signs:** Unnecessarily complex workflow, slower API builds

### Pitfall 6: Permissions for Package Write
**What goes wrong:** 403 Forbidden when pushing to GHCR
**Why it happens:** Default GITHUB_TOKEN lacks `packages: write`
**How to avoid:** Explicitly set permissions in each job:
```yaml
permissions:
  contents: read
  packages: write
```
**Warning signs:** "denied: permission_denied" in push logs

### Pitfall 7: Metadata JSON Environment Variable
**What goes wrong:** `DOCKER_METADATA_OUTPUT_JSON` empty in merge job
**Why it happens:** Environment variable only set after metadata-action runs in that job
**How to avoid:** Run `docker/metadata-action` in the merge job, not just build jobs
**Warning signs:** Empty tags in `imagetools create` command

## Code Examples

Verified patterns from official sources:

### Complete Multi-Arch Workflow for Worker Image
```yaml
# Source: https://github.com/sredevopsorg/multi-arch-docker-github-workflow
# Adapted for qm-nmr-calc worker image
name: Publish Docker Images

on:
  push:
    tags:
      - 'v*'

env:
  REGISTRY: ghcr.io
  WORKER_IMAGE: ghcr.io/${{ github.repository_owner }}/qm-nmr-calc-worker

jobs:
  # ==========================================================================
  # Build API image (existing pattern - QEMU ok for small image)
  # ==========================================================================
  build-api:
    runs-on: ubuntu-latest
    permissions:
      contents: read
      packages: write
    steps:
      - uses: actions/checkout@v4
      - uses: docker/setup-qemu-action@v3
      - uses: docker/setup-buildx-action@v3
      - uses: docker/login-action@v3
        with:
          registry: ${{ env.REGISTRY }}
          username: ${{ github.actor }}
          password: ${{ secrets.GITHUB_TOKEN }}
      - uses: docker/metadata-action@v5
        id: meta
        with:
          images: ${{ env.REGISTRY }}/${{ github.repository_owner }}/qm-nmr-calc-api
          tags: |
            type=semver,pattern={{version}}
            type=semver,pattern={{major}}.{{minor}}
            type=raw,value=latest
      - uses: docker/build-push-action@v6
        with:
          context: .
          file: Dockerfile.api
          platforms: linux/amd64,linux/arm64
          push: true
          tags: ${{ steps.meta.outputs.tags }}
          labels: ${{ steps.meta.outputs.labels }}
          cache-from: type=gha,scope=api
          cache-to: type=gha,mode=max,scope=api

  # ==========================================================================
  # Build Worker image per-platform on native runners
  # ==========================================================================
  build-worker:
    strategy:
      fail-fast: false
      matrix:
        platform:
          - linux/amd64
          - linux/arm64

    runs-on: ${{ matrix.platform == 'linux/amd64' && 'ubuntu-latest' || 'ubuntu-24.04-arm' }}
    permissions:
      contents: read
      packages: write

    steps:
      - name: Prepare platform pair
        run: |
          platform=${{ matrix.platform }}
          echo "PLATFORM_PAIR=${platform//\//-}" >> $GITHUB_ENV

      - uses: actions/checkout@v4

      - uses: docker/metadata-action@v5
        id: meta
        with:
          images: ${{ env.WORKER_IMAGE }}

      - name: Set up Docker Buildx
        uses: docker/setup-buildx-action@v3
        with:
          platforms: ${{ matrix.platform }}

      - uses: docker/login-action@v3
        with:
          registry: ${{ env.REGISTRY }}
          username: ${{ github.actor }}
          password: ${{ secrets.GITHUB_TOKEN }}

      - name: Build and push by digest
        id: build
        uses: docker/build-push-action@v6
        with:
          context: .
          file: ${{ matrix.platform == 'linux/amd64' && 'Dockerfile.worker' || 'Dockerfile.worker.arm64' }}
          platforms: ${{ matrix.platform }}
          labels: ${{ steps.meta.outputs.labels }}
          outputs: type=image,name=${{ env.WORKER_IMAGE }},push-by-digest=true,name-canonical=true,push=true
          cache-from: type=gha,scope=worker-${{ env.PLATFORM_PAIR }}
          cache-to: type=gha,scope=worker-${{ env.PLATFORM_PAIR }}

      - name: Export digest
        run: |
          mkdir -p /tmp/digests
          digest="${{ steps.build.outputs.digest }}"
          touch "/tmp/digests/${digest#sha256:}"

      - name: Upload digest
        uses: actions/upload-artifact@v4
        with:
          name: digests-${{ env.PLATFORM_PAIR }}
          path: /tmp/digests/*
          if-no-files-found: error
          retention-days: 1

  # ==========================================================================
  # Merge Worker digests into multi-arch manifest
  # ==========================================================================
  merge-worker:
    runs-on: ubuntu-latest
    needs: [build-worker]
    permissions:
      contents: read
      packages: write

    steps:
      - name: Download digests
        uses: actions/download-artifact@v4
        with:
          path: /tmp/digests
          pattern: digests-*
          merge-multiple: true

      - uses: docker/metadata-action@v5
        id: meta
        with:
          images: ${{ env.WORKER_IMAGE }}
          tags: |
            type=semver,pattern={{version}}
            type=semver,pattern={{major}}.{{minor}}
            type=raw,value=latest

      - uses: docker/setup-buildx-action@v3

      - uses: docker/login-action@v3
        with:
          registry: ${{ env.REGISTRY }}
          username: ${{ github.actor }}
          password: ${{ secrets.GITHUB_TOKEN }}

      - name: Create manifest list and push
        working-directory: /tmp/digests
        run: |
          docker buildx imagetools create \
            $(jq -cr '.tags | map("-t " + .) | join(" ")' <<< "$DOCKER_METADATA_OUTPUT_JSON") \
            $(printf '${{ env.WORKER_IMAGE }}@sha256:%s ' *)

      - name: Inspect image
        run: |
          docker buildx imagetools inspect '${{ env.WORKER_IMAGE }}:${{ steps.meta.outputs.version }}'
```

### Selective Dockerfile by Platform
```yaml
# Source: Pattern for using different Dockerfiles per architecture
file: ${{ matrix.platform == 'linux/amd64' && 'Dockerfile.worker' || 'Dockerfile.worker.arm64' }}
```

### Inspecting Multi-Arch Manifest
```bash
# Verify manifest contains both architectures
docker buildx imagetools inspect ghcr.io/owner/qm-nmr-calc-worker:latest
# Should show:
# Name:      ghcr.io/owner/qm-nmr-calc-worker:latest
# MediaType: application/vnd.oci.image.index.v1+json
# Manifests:
#   Platform: linux/amd64
#   Platform: linux/arm64
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| QEMU for all platforms | Native ARM runners | Jan 2025 (public preview) | 10x faster ARM builds |
| Self-hosted ARM runners | `ubuntu-24.04-arm` | Jan 2025 | No infrastructure management |
| `actions/upload-artifact@v3` | `actions/upload-artifact@v4` | Jan 2025 (v3 deprecated) | Immutable artifacts, faster |
| Single-job multi-platform | Multi-job with native runners | 2024+ | Bypasses QEMU limitations |
| `docker manifest create` | `docker buildx imagetools create` | 2022+ | Handles attestations, OCI |

**Deprecated/outdated:**
- `actions/upload-artifact@v3`: Deprecated January 2025, must use v4
- `docker/build-push-action@v5`: Current is v6
- Self-hosted ARM runners: No longer needed for most use cases
- QEMU for production builds: Only acceptable for small images (<500MB)

## Open Questions

Things that couldn't be fully resolved:

1. **Build Time for ARM64 Worker**
   - What we know: Native ARM runners eliminate QEMU overhead
   - What's unclear: Exact build time for 2.1GB worker with conda-forge packages
   - Recommendation: Monitor first build, expect ~10-15 minutes (similar to amd64)

2. **Cache Efficiency Across Architectures**
   - What we know: Each platform gets own cache scope
   - What's unclear: Whether micromamba layer caching works efficiently on ARM
   - Recommendation: Use `mode=max` for cache-to; monitor cache hit rates

3. **Repository Visibility Requirement**
   - What we know: `ubuntu-24.04-arm` is free for PUBLIC repos
   - What's unclear: If repo becomes private, workflow will break
   - Recommendation: Document dependency on public repo status in workflow comments

## Sources

### Primary (HIGH confidence)
- [GitHub Changelog - Linux ARM64 Runners](https://github.blog/changelog/2025-01-16-linux-arm64-hosted-runners-now-available-for-free-in-public-repositories-public-preview/) - Runner availability and labels
- [sredevopsorg/multi-arch-docker-github-workflow](https://github.com/sredevopsorg/multi-arch-docker-github-workflow) - Complete reference implementation
- [GitHub Actions Artifacts v4](https://github.blog/news-insights/product-news/get-started-with-v4-of-github-actions-artifacts/) - Artifact best practices
- GitHub API queries for latest action versions (docker/build-push-action@v6.18.0, docker/metadata-action@v5.10.0, etc.)

### Secondary (MEDIUM confidence)
- [Blacksmith Blog - Multi-Platform ARM64 Images](https://www.blacksmith.sh/blog/building-multi-platform-docker-images-for-arm64-in-github-actions) - Performance comparisons
- [InfoQ - GitHub Actions ARM64 Runners](https://www.infoq.com/news/2025/02/github-actions-linux-arm64/) - Announcement coverage
- [GitHub Community Discussion #148648](https://github.com/orgs/community/discussions/148648) - ARM64 runner details

### Tertiary (LOW confidence)
- Various community blog posts on cache optimization - Used for best practice patterns only

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - Official Docker Actions, verified versions via GitHub API
- Architecture: HIGH - Docker-documented pattern, reference implementation verified
- Pitfalls: HIGH - Well-documented in official sources and community discussions

**Research date:** 2026-02-04
**Valid until:** 2026-05-04 (GitHub Actions ecosystem relatively stable; ARM runners are GA)
