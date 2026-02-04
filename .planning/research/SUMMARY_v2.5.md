# Research Summary: v2.5 ARM64 Docker Support

**Project:** qm-nmr-calc v2.5 - ARM64 Docker Support
**Domain:** Multi-architecture Docker deployment for computational chemistry
**Researched:** 2026-02-04
**Confidence:** HIGH

## Executive Summary

ARM64 native support for the qm-nmr-calc worker container is now **fully achievable** using conda-forge packages. The key insight driving this milestone is that all required computational chemistry packages (NWChem 7.3.1, xTB 6.7.1, CREST 3.0.2) are now available on conda-forge for linux-aarch64 -- a capability that did not exist when v2.4 was developed. This eliminates the fundamental blocker noted in the current Dockerfile.worker comment that "arm64 not supported - CREST/xTB lack arm64 binaries."

The recommended approach is a **separate ARM64 Dockerfile** (`Dockerfile.worker.arm64`) using `mambaorg/micromamba:2.5.0-debian12-slim` as the base image, with all computational chemistry packages installed via conda-forge. This preserves the proven x86_64 approach while enabling native ARM64 execution. GitHub Actions free ARM64 runners (`ubuntu-24.04-arm`, generally available since August 2025) enable fast native builds without QEMU emulation overhead.

The primary risks are SIGILL crashes from x86-specific SIMD instructions and OpenBLAS misconfiguration -- both well-understood issues with clear mitigations. The conda-forge approach sidesteps most pitfalls by providing pre-built, architecture-appropriate binaries. The main trade-off is slightly larger image size (~500MB overhead from conda) compared to the pre-compiled binary approach used for x86_64.

## Key Findings

### Recommended Stack

The ARM64 worker requires a fundamentally different packaging approach than x86_64 due to the absence of pre-compiled ARM64 binaries from upstream projects (grimme-lab/xtb and crest-lab/crest only publish x86_64 releases).

**Core technologies:**
- **mambaorg/micromamba:2.5.0-debian12-slim** (base image): Multi-arch support, 2x faster package resolution than miniconda, Debian base avoids Alpine/musl issues
- **conda-forge packages** (NWChem 7.3.1, xTB 6.7.1, CREST 3.0.2, RDKit 2025.09.5): All verified available for linux-aarch64 as of 2026-02-04
- **GitHub Actions ubuntu-24.04-arm**: Free native ARM64 runner for public repos, avoids 22x QEMU slowdown

**Installation strategy:**
- Heavy scientific packages via conda-forge (native ARM64 binaries)
- Pure-Python packages (huey, jinja2, pydantic) via pip within conda environment
- Single `env-worker-arm64.yaml` defines entire environment

### Expected Features

**Must have (table stakes):**
- TS-1: Native execution on Apple Silicon (M1/M2/M3) without SIGILL crashes
- TS-2: Automatic architecture detection via multi-arch manifest (same docker compose command works everywhere)
- TS-3: Numerically equivalent results to x86_64 (< 0.01 ppm tolerance)
- TS-4: All calculation modes work (single-conformer, CREST ensemble, RDKit-only)
- TS-5: Multi-arch image publishing to GHCR

**Should have (differentiators included in MVP):**
- D-1: AWS Graviton deployment support (15-25% cost savings)
- D-2: Local development without emulation (primary user motivation)

**Defer (post-MVP optimization):**
- D-3: Cache-efficient multi-arch builds
- D-4: Architecture-specific performance tuning (e.g., different OMP_NUM_THREADS defaults)

### Architecture Approach

Create a separate `Dockerfile.worker.arm64` rather than a unified Dockerfile with build args. This separation provides cleaner maintenance, easier debugging, and preserves the proven x86_64 approach without risk of regression. The CI/CD workflow builds each architecture on native runners in parallel, then merges into a single multi-arch manifest using `docker buildx imagetools create`.

**Major components:**
1. **Dockerfile.worker.arm64** (new): micromamba base + conda-forge packages for all computational chemistry dependencies
2. **env-worker-arm64.yaml** (new): Conda environment definition with pinned versions matching x86_64 parity
3. **publish-images.yml** (modified): Separate build jobs for amd64/arm64 + manifest merge job

**Key differences from x86_64:**

| Aspect | x86_64 | ARM64 |
|--------|--------|-------|
| Base image | ghcr.io/nwchemgit/nwchem-dev/amd64 | mambaorg/micromamba:2.5.0-debian12-slim |
| NWChem | Pre-installed in base | conda-forge package |
| xTB/CREST | GitHub release binaries | conda-forge packages |
| Estimated image size | ~2-3GB | ~3-4GB |

### Critical Pitfalls

1. **SIGILL from x86 SIMD instructions** (BLOCKER): Pre-compiled x86_64 binaries use AVX/AVX2 that cannot be emulated. **Prevention:** Use conda-forge packages exclusively for ARM64.

2. **Wrong OpenBLAS/BLAS configuration** (BLOCKER): Intel MKL unavailable on ARM64; generic OpenBLAS kernel is 10-100x slower. **Prevention:** conda-forge NWChem includes properly configured OpenBLAS; verify with `ldd nwchem | grep blas`.

3. **OpenBLAS thread creation failures** (HIGH): Container crashes with "pthread_create failed". **Prevention:** Set `ENV OPENBLAS_NUM_THREADS=4` explicitly in Dockerfile.

4. **Miniconda architecture mismatch** (BLOCKER): x86_64 installer causes Python to run under emulation. **Prevention:** Use micromamba multi-arch image or Miniforge ARM64 installer.

5. **Build time explosion with QEMU** (HIGH for CI): ARM64 build via emulation takes 30-60 minutes vs 3-5 minutes native. **Prevention:** Use `ubuntu-24.04-arm` GitHub Actions runner.

## Implications for Roadmap

Based on research, suggested phase structure:

### Phase 1: ARM64 Dockerfile Creation

**Rationale:** Foundation must exist before CI/CD integration or testing
**Delivers:** `Dockerfile.worker.arm64` and `env-worker-arm64.yaml`
**Addresses:** TS-1 (native execution), TS-4 (all modes work)
**Avoids:** Pitfalls #1 (SIGILL), #2 (BLAS), #3 (threading), #4 (Miniconda)
**Estimated effort:** Small (template provided in STACK.md)

### Phase 2: Local Validation

**Rationale:** Must verify container works before CI investment
**Delivers:** Validated ARM64 container on Apple Silicon
**Implements:** Validation checklist from PITFALLS.md
**Tests:**
- Worker starts without SIGILL
- NWChem DFT calculation completes
- CREST conformer search completes
- xTB ranking works
- Results match x86_64 (< 0.01 ppm)

### Phase 3: CI/CD Integration

**Rationale:** Automation required for sustainable multi-arch publishing
**Delivers:** Updated `publish-images.yml` with ARM64 support
**Addresses:** TS-5 (multi-arch publishing), TS-2 (auto-detection)
**Avoids:** Pitfall #5 (QEMU slowdown) via native runner
**Components:**
- `build-worker-amd64` job (ubuntu-latest)
- `build-worker-arm64` job (ubuntu-24.04-arm)
- `create-worker-manifest` job (merge into single tag)

### Phase 4: Documentation and Release

**Rationale:** Users need awareness of new capability
**Delivers:** v2.5.0 release with ARM64 support
**Addresses:** D-1 (Graviton documentation), D-2 (local dev documentation)
**Includes:**
- README update noting ARM64 support
- Performance expectations documentation
- Any architecture-specific notes

### Phase Ordering Rationale

- **Phase 1 before 2:** Cannot test without Dockerfile
- **Phase 2 before 3:** Must validate locally before investing in CI changes
- **Phase 3 before 4:** Cannot document/release without working CI
- **Parallel potential:** Phases 1-2 are sequential; Phase 3 can begin once Phase 2 validation is partially complete

### Research Flags

**Phases with standard patterns (no additional research needed):**
- **Phase 1:** Dockerfile patterns well-documented; conda-forge packages verified; template in STACK.md
- **Phase 3:** GitHub Actions multi-arch patterns well-documented; reference implementation available

**Phases that may need attention during execution:**
- **Phase 2 (validation):** May discover NWChem basis set path issues specific to conda-forge package. Mitigation: Verify `NWCHEM_BASIS_LIBRARY` environment variable.
- **Phase 4 (documentation):** May need performance benchmarking if significant differences found. Mitigation: Document expectations, not guarantees.

## Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| Stack | HIGH | All packages verified on conda-forge; base image verified on Docker Hub |
| Features | HIGH | Table stakes are standard multi-arch expectations; clearly scoped |
| Architecture | HIGH | Separate Dockerfile approach is proven pattern; NWChem ARM64 image verified |
| Pitfalls | MEDIUM-HIGH | Based on issue trackers and documented problems; some edge cases possible |

**Overall confidence:** HIGH

### Gaps to Address

- **NWChem basis set paths:** conda-forge NWChem may have different `NWCHEM_BASIS_LIBRARY` path than base image. Verify during Phase 2 validation.
- **OpenMPI configuration:** conda-forge NWChem uses OpenMPI; verify `OMPI_ALLOW_RUN_AS_ROOT` settings work correctly.
- **Performance benchmarking:** No quantitative data on ARM64 vs x86_64 calculation times. Document during Phase 2 but accept "equivalent" as sufficient.

## Sources

### Primary (HIGH confidence)
- [conda-forge NWChem](https://anaconda.org/conda-forge/nwchem) - Version 7.3.1, linux-aarch64 verified 2026-02-04
- [conda-forge xTB](https://anaconda.org/conda-forge/xtb) - Version 6.7.1, linux-aarch64 verified 2026-02-04
- [conda-forge CREST](https://anaconda.org/conda-forge/crest) - Version 3.0.2, linux-aarch64 verified 2026-02-04
- [Docker Multi-Platform Documentation](https://docs.docker.com/build/building/multi-platform/) - Buildx patterns
- [GitHub ARM64 Runners GA](https://github.blog/changelog/2025-08-07-arm64-hosted-runners-for-public-repositories-are-now-generally-available/) - ubuntu-24.04-arm availability

### Secondary (MEDIUM confidence)
- [mambaorg/micromamba Docker Hub](https://hub.docker.com/r/mambaorg/micromamba) - Multi-arch manifest verified
- [NWChem Container Documentation](https://nwchemgit.github.io/Containers.html) - ARM64 image availability
- [actuated.com Native ARM64 Performance](https://actuated.com/blog/native-arm64-for-github-actions) - 22x speedup data

### Issue Trackers (pitfall validation)
- [Docker for Mac #7137](https://github.com/docker/for-mac/issues/7137) - Rosetta breaks amd64 images with AVX
- [OpenBLAS #4825](https://github.com/OpenMathLib/OpenBLAS/issues/4825) - aarch64 Docker build threading issues
- [OpenMPI #12600](https://github.com/open-mpi/ompi/issues/12600) - Fortran module path issues

---
*Research completed: 2026-02-04*
*Ready for roadmap: yes*
