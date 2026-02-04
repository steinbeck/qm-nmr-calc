# Feature Landscape: ARM64 Docker Support

**Domain:** Multi-architecture Docker deployment for scientific computing
**Project:** qm-nmr-calc v2.5
**Researched:** 2026-02-04
**Overall Confidence:** HIGH

## Executive Summary

ARM64 Docker support enables qm-nmr-calc to run natively on Apple Silicon Macs (M1/M2/M3) and ARM-based cloud instances (AWS Graviton, Azure Arm VMs). The current x86_64-only worker fails with SIGILL (illegal instruction) on Apple Silicon due to AVX instructions in compiled binaries.

The core enabler is that **all required scientific packages (NWChem 7.3.1, CREST 3.0.2, xTB 6.7.1) are now available on conda-forge for linux-aarch64**. This wasn't true when v2.4 was built -- hence the original amd64-only decision.

## Table Stakes

Features users expect when ARM64 support is claimed. Missing = ARM64 support feels incomplete.

### TS-1: Native Execution on Apple Silicon

| Aspect | Requirement |
|--------|-------------|
| **What** | Worker container runs without QEMU emulation on M1/M2/M3 Macs |
| **Why Expected** | "ARM64 support" meaningless if it still emulates x86_64 |
| **User Experience** | `docker compose up` pulls ARM64 image automatically; calculations run at native speed |
| **Complexity** | Medium (requires new Dockerfile based on conda-forge instead of NWChem base image) |

**Verification:** Run a small NMR calculation end-to-end on Apple Silicon without SIGILL or emulation warnings.

### TS-2: Automatic Architecture Detection

| Aspect | Requirement |
|--------|-------------|
| **What** | Same `docker compose up` command works on x86_64 and ARM64 |
| **Why Expected** | Users shouldn't need to know their architecture or modify compose files |
| **User Experience** | Pull `ghcr.io/steinbeck/qm-nmr-calc-worker:latest` and Docker selects correct platform |
| **Complexity** | Low (multi-arch manifest handles this automatically) |

**How it works:** Multi-architecture images use a manifest list (also called "fat manifest"). When Docker pulls an image, the registry returns the manifest list. Docker automatically selects the correct variant based on the host's architecture (linux/amd64 or linux/arm64).

### TS-3: Correct Numerical Results

| Aspect | Requirement |
|--------|-------------|
| **What** | ARM64 calculations produce numerically equivalent results to x86_64 |
| **Why Expected** | Chemistry doesn't change with CPU architecture |
| **User Experience** | Same molecule, same settings = same chemical shifts (within floating-point tolerance) |
| **Complexity** | Low (verify, not implement -- conda-forge NWChem passes upstream tests) |

**Tolerance:** Expect sub-0.01 ppm differences due to floating-point determinism. Any larger discrepancy indicates a build problem.

### TS-4: All Calculation Modes Work

| Aspect | Requirement |
|--------|-------------|
| **What** | Single-conformer, CREST ensemble, and RDKit-only modes all function on ARM64 |
| **Why Expected** | ARM64 support means full feature parity, not degraded mode |
| **User Experience** | All existing API parameters work identically |
| **Complexity** | Medium (requires NWChem, CREST, and xTB all working) |

**Dependencies verified on conda-forge linux-aarch64:**
- NWChem 7.3.1 (DFT calculations)
- CREST 3.0.2 (conformer generation)
- xTB 6.7.1 (conformer ranking, CREST backend)
- RDKit (molecule handling, conformer generation fallback)

### TS-5: Multi-Arch Image Publishing

| Aspect | Requirement |
|--------|-------------|
| **What** | CI/CD publishes both linux/amd64 and linux/arm64 images to GHCR |
| **Why Expected** | Users need ARM64 images available in registry |
| **User Experience** | Standard `docker pull` works on both architectures |
| **Complexity** | Medium (requires buildx configuration in GitHub Actions) |

**CI Strategy Options:**
1. **QEMU emulation** (slow but simple): Single runner builds both architectures
2. **Native runners** (fast, recommended): Matrix strategy with `ubuntu-latest` (amd64) and `ubuntu-24.04-arm` (arm64)
3. **Parallel build + manifest merge**: Build architectures separately, merge manifest

## Differentiators

Features that add value beyond basic ARM64 support. Not expected, but appreciated.

### D-1: AWS Graviton Deployment Support

| Aspect | Detail |
|--------|--------|
| **What** | Documented and tested deployment to AWS Graviton instances |
| **Value Proposition** | 20-40% cost savings vs x86 instances with equivalent or better performance |
| **User Experience** | Deploy to Graviton (t4g, c7g, etc.) with same docker-compose.yml |
| **Complexity** | Low (ARM64 support inherently enables this) |

**Performance expectations from AWS benchmarks:**
- 15-25% lower compute cost vs equivalent x86 instances
- Web applications: up to 30% faster on Graviton4 vs Graviton3
- Better price-performance ratio for compute-intensive workloads

### D-2: Local Development Without Emulation

| Aspect | Detail |
|--------|--------|
| **What** | Developers on Apple Silicon can run full stack natively during development |
| **Value Proposition** | No more QEMU overhead (85% slower than native reported) |
| **User Experience** | Build and run tests at native speed; `docker compose up` starts quickly |
| **Complexity** | Included with ARM64 support |

**Current pain point:** Developers on M1/M2/M3 Macs cannot run the worker locally. They must either:
- Use a remote x86 server
- Attempt QEMU emulation (SIGILL due to AVX instructions)
- Skip worker testing entirely

### D-3: Cache-Efficient Multi-Arch Builds

| Aspect | Detail |
|--------|--------|
| **What** | CI caches layer builds per-architecture to speed rebuilds |
| **Value Proposition** | Faster CI, lower costs |
| **User Experience** | Subsequent builds significantly faster |
| **Complexity** | Low (buildx cache configuration) |

### D-4: Architecture-Specific Performance Tuning

| Aspect | Detail |
|--------|--------|
| **What** | Different default parallelism settings for ARM64 vs x86_64 |
| **Value Proposition** | Optimal performance on each architecture |
| **User Experience** | Transparent; auto-detected |
| **Complexity** | Low (environment variable defaults) |

**Consideration:** ARM cores often have different characteristics (e.g., efficiency cores on Apple Silicon). May need different `NWCHEM_NPROC` and `OMP_NUM_THREADS` defaults.

## Anti-Features

Things to deliberately NOT build. Common mistakes in multi-arch Docker work.

### AF-1: Separate Image Tags for Architectures

| Anti-Feature | Why Avoid | What to Do Instead |
|--------------|-----------|-------------------|
| `worker:latest-arm64` and `worker:latest-amd64` tags | Forces users to know their architecture; breaks automation | Use multi-arch manifest with single tag (`:latest`, `:v2.5.0`) |

Docker's manifest list exists specifically to solve this. A single tag contains multiple architecture variants. The registry handles selection automatically.

### AF-2: QEMU-Based CI Builds

| Anti-Feature | Why Avoid | What to Do Instead |
|--------------|-----------|-------------------|
| Build ARM64 images via QEMU on x86 runners | 10x slower builds; 30+ minute CI times; unreliable for compute-heavy compilation | Use native ARM64 runners (GitHub `ubuntu-24.04-arm`) |

QEMU emulation for builds is acceptable for simple images but problematic for scientific computing containers that compile native code.

### AF-3: Architecture-Specific Compose Files

| Anti-Feature | Why Avoid | What to Do Instead |
|--------------|-----------|-------------------|
| `docker-compose-arm64.yml` vs `docker-compose.yml` | Confusing; error-prone; maintenance burden | Single `docker-compose.yml` that works everywhere via multi-arch images |

The whole point of multi-arch images is to make architecture transparent to users.

### AF-4: Rosetta/QEMU Fallback Mode

| Anti-Feature | Why Avoid | What to Do Instead |
|--------------|-----------|-------------------|
| Documentation suggesting "use Rosetta if ARM64 doesn't work" | Defeats purpose; 20-85% performance penalty; may still fail on AVX | Native ARM64 binaries only |

The current x86_64 worker uses AVX instructions that cannot be emulated even with Rosetta. A proper ARM64 image is the only solution.

### AF-5: Different Feature Sets Per Architecture

| Anti-Feature | Why Avoid | What to Do Instead |
|--------------|-----------|-------------------|
| "ARM64 supports single-conformer only" | Users expect feature parity; creates support burden | Full feature parity or don't ship |

If CREST or xTB had ARM64 issues, better to delay ARM64 release than ship degraded mode.

### AF-6: Build-Time Architecture Detection in Dockerfile

| Anti-Feature | Why Avoid | What to Do Instead |
|--------------|-----------|-------------------|
| `RUN if [ "$(uname -m)" = "aarch64" ]; then ...` | Breaks buildx cross-compilation; complicates maintenance | Use separate Dockerfiles or buildx `--platform` with multi-stage builds |

Buildx handles architecture-specific builds properly. Manual detection is fragile and defeats the tooling.

## Feature Dependencies

```
TS-5: Multi-Arch Publishing
    |
    +---> TS-2: Auto Architecture Detection
    |         |
    |         +---> TS-1: Native Apple Silicon Execution
    |         |
    |         +---> D-1: Graviton Deployment
    |
    +---> TS-4: All Calculation Modes
              |
              +---> TS-3: Correct Numerical Results
```

**Dependency explanation:**
- Multi-arch publishing (TS-5) is the foundation -- images must exist in registry
- Auto-detection (TS-2) enables transparent user experience
- Native execution (TS-1) and Graviton (D-1) are both enabled by multi-arch images
- Correct results (TS-3) validates that ARM64 NWChem is scientifically sound
- All modes (TS-4) requires all dependencies (NWChem, CREST, xTB) working

## MVP Recommendation

For v2.5 MVP, prioritize:

1. **TS-1: Native Apple Silicon Execution** - Core user need; blocks local development
2. **TS-2: Auto Architecture Detection** - Essential UX; comes "for free" with multi-arch
3. **TS-3: Correct Numerical Results** - Must verify before release
4. **TS-4: All Calculation Modes** - Feature parity required
5. **TS-5: Multi-Arch Publishing** - Enables users to pull images

All table stakes are MVP requirements.

**Defer to post-MVP:**
- D-3 (Cache-efficient builds): Optimization; can iterate
- D-4 (Architecture-specific tuning): Optimization; defaults are fine initially

**Include in MVP:**
- D-1 (Graviton support): Zero extra work if ARM64 works
- D-2 (Local development): Primary motivation for this milestone

## Performance Expectations

### Local Development (Apple Silicon)

| Scenario | Current (x86_64 emulation) | Expected (ARM64 native) |
|----------|---------------------------|------------------------|
| Container start | N/A (SIGILL crash) | < 5 seconds |
| Small molecule NMR | N/A | ~2-5 minutes |
| Conformer search | N/A | Same as x86_64 |

### AWS Graviton Deployment

| Metric | x86_64 Instance | Graviton Instance | Notes |
|--------|-----------------|-------------------|-------|
| Hourly cost | Baseline | 15-25% lower | Instance type dependent |
| Calculation time | Baseline | ~equivalent | May vary by workload |
| Price-performance | Baseline | Better | AWS benchmarks show 40% improvement |

### CI Build Times

| Strategy | Estimated Time | Notes |
|----------|----------------|-------|
| QEMU emulation | 30-60 minutes | Not recommended |
| Native runners (parallel) | 10-15 minutes | Recommended |
| Native runners (sequential) | 15-25 minutes | Acceptable |

## Validation Checklist

Before claiming ARM64 support complete:

- [ ] Worker starts on Apple Silicon without SIGILL
- [ ] Worker starts on Linux aarch64 (docker run test)
- [ ] NWChem runs basic DFT calculation
- [ ] CREST generates conformers
- [ ] xTB ranks conformers
- [ ] Full NMR calculation produces results
- [ ] Results match x86_64 within tolerance (< 0.01 ppm)
- [ ] GHCR shows multi-arch manifest
- [ ] `docker pull` on both architectures succeeds
- [ ] Single `docker-compose.yml` works on both

## Sources

### Docker Multi-Architecture
- [Docker Multi-Platform Builds Documentation](https://docs.docker.com/build/building/multi-platform/) - Official Docker documentation on multi-platform builds
- [Docker Multi-Arch Blog Post](https://www.docker.com/blog/multi-arch-images/) - Building Multi-Arch Images for Arm and x86
- [Docker Manifest Documentation](https://docs.docker.com/reference/cli/docker/manifest/) - How manifest lists work

### GitHub Actions Multi-Arch
- [Blacksmith ARM64 Runners Guide](https://www.blacksmith.sh/blog/building-multi-platform-docker-images-for-arm64-in-github-actions) - Native ARM64 builds in GitHub Actions
- [GitHub Multi-Arch Workflow Example](https://github.com/sredevopsorg/multi-arch-docker-github-workflow) - How to build without QEMU
- [Docker Build Push Action Multi-Platform](https://docs.docker.com/build/ci/github-actions/multi-platform/) - Official GitHub Actions guide

### conda-forge Packages (verified current)
- [NWChem on conda-forge](https://anaconda.org/conda-forge/nwchem) - Version 7.3.1, supports linux-aarch64, osx-arm64
- [CREST on conda-forge](https://anaconda.org/conda-forge/crest) - Version 3.0.2, supports linux-aarch64, osx-arm64
- [xTB on conda-forge](https://anaconda.org/conda-forge/xtb) - Version 6.7.1, supports linux-aarch64, osx-arm64

### AWS Graviton
- [AWS Graviton Documentation](https://aws.amazon.com/ec2/graviton/) - Official EC2 Graviton page
- [AWS Graviton Cost Savings Guide](https://sedai.io/blog/aws-graviton-usage-guide) - Real-world benchmarks and savings
- [AWS Graviton Performance Blog](https://aws.amazon.com/blogs/compute/how-potential-performance-upside-with-aws-graviton-helps-reduce-your-costs-further/) - Performance upside analysis

### Apple Silicon Docker
- [Docker on Apple Silicon Performance Guide](https://oneuptime.com/blog/post/2026-01-16-docker-mac-apple-silicon/view) - Best practices for M1/M2/M3
- [QEMU vs Rosetta Performance Analysis](https://medium.com/@guillem.riera/the-most-performant-docker-setup-on-macos-apple-silicon-m1-m2-m3-for-x64-amd64-compatibility-da5100e2557d) - Performance comparison
