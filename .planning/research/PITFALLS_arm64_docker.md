# Domain Pitfalls: ARM64 Docker for NWChem/xTB/CREST

**Domain:** ARM64 containers for computational chemistry (NWChem/xTB/CREST)
**Researched:** 2026-02-04
**Overall Confidence:** MEDIUM-HIGH

This document catalogs pitfalls specific to building and running ARM64 (aarch64) Docker containers for the qm-nmr-calc computational chemistry stack. These are in addition to the general Docker deployment pitfalls documented in `PITFALLS_v2.4_docker_deployment.md`.

**Context:** The project already hit SIGILL from AVX instructions in x86_64 binaries under Rosetta 2 emulation, demonstrating the real-world impact of these architecture issues.

---

## Critical Pitfalls

Mistakes that cause container failure, SIGILL crashes, or require complete rebuilds.

### Pitfall 1: SIGILL from x86-Specific SIMD Instructions (AVX/AVX2/AVX-512)

**What goes wrong:** Containers built for x86_64 crash with `SIGILL: illegal instruction` when run on ARM64 hosts under Rosetta 2 or QEMU emulation.

**Why it happens:**
- Scientific computing binaries are often compiled with AVX/AVX2/AVX-512 optimizations for x86_64
- Pre-compiled binaries from releases (xTB, CREST) are typically x86_64-only with SIMD optimizations
- Rosetta 2 and QEMU cannot emulate all AVX instructions (AVX support in QEMU is recent and incomplete on ARM64)
- Intel MKL is x86-only and will SIGILL on ARM64

**Consequences:**
- Immediate crash on first calculation attempt
- No useful error message - just "Illegal instruction" or exit code 132
- Works in dev (x86), fails in production (ARM64)
- The current `Dockerfile.worker` downloads `xtb-6.7.1-linux-x86_64.tar.xz` and `crest-gnu-12-ubuntu-latest.tar.xz` which are x86_64 binaries

**Prevention:**
1. **Use architecture-specific base images:**
   ```dockerfile
   # Instead of: FROM ghcr.io/nwchemgit/nwchem-dev/amd64:latest
   # Use: FROM ghcr.io/nwchemgit/nwchem-dev/arm64:latest (for ARM64 builds)
   ```

2. **Download architecture-appropriate binaries:**
   ```dockerfile
   ARG TARGETARCH
   RUN if [ "$TARGETARCH" = "arm64" ]; then \
         # Use conda-forge packages or compile from source
         conda install -c conda-forge xtb crest; \
       else \
         # Use pre-compiled x86_64 binaries
         wget xtb-6.7.1-linux-x86_64.tar.xz ...; \
       fi
   ```

3. **Build from source for ARM64:** xTB and CREST must be compiled with ARM64-appropriate flags (no AVX, use NEON)

4. **Use conda-forge for ARM64 packages:**
   - xTB is available on conda-forge for linux-aarch64
   - CREST is available on conda-forge for linux-aarch64
   - NWChem is available on conda-forge for linux-aarch64

**Detection (warning signs):**
- Exit code 132 (128 + 4 = SIGILL)
- "Illegal instruction" in dmesg or container logs
- Works on x86_64 host, fails on ARM64 host
- Error occurs immediately when binary starts, not during calculation

**Which phase should address:** Phase 1 - This is a BLOCKER. Container will not function without architecture-appropriate binaries.

**Sources:**
- [Docker for Mac Issue #7137 - Rosetta emulation breaks amd64 images](https://github.com/docker/for-mac/issues/7137)
- [Docker for Mac Issue #6620 - QEMU AVX support request](https://github.com/docker/for-mac/issues/6620)
- [conda-forge xTB packages](https://anaconda.org/conda-forge/xtb) - linux-aarch64 builds available

---

### Pitfall 2: Wrong OpenBLAS/BLAS Configuration for ARM64

**What goes wrong:** NWChem calculations produce wrong results, crash, or run extremely slowly on ARM64.

**Why it happens:**
- Intel MKL is not available for ARM64 - must use OpenBLAS
- OpenBLAS requires correct TARGET specification for ARM64 optimization
- Default OpenBLAS may use generic kernels that are 10-100x slower
- MacOS Accelerate framework cannot be used in Linux containers
- NWChem documentation explicitly warns: "Please do not use the Mac OS X default BLAS and LAPACK libraries"

**Consequences:**
- Calculations that take 5 minutes on x86_64 take hours on ARM64
- Subtle numerical differences leading to wrong chemical results
- NWChem crashes with "Error: cannot find basis set" if BLAS is misconfigured

**Prevention:**
1. **Use OpenBLAS with correct ARM64 TARGET:**
   ```bash
   # In Dockerfile or build script
   export TARGET=ARMV8  # Generic ARM64
   # Or for specific chips:
   export TARGET=CORTEXA72  # For AWS Graviton, etc.
   export TARGET=NEOVERSEN2  # For newer ARM servers
   ```

2. **Verify OpenBLAS is linked correctly:**
   ```bash
   ldd /path/to/nwchem | grep -i blas
   # Should show libopenblas, NOT libmkl
   ```

3. **Use official NWChem ARM64 image:** The `ghcr.io/nwchemgit/nwchem-dev/arm64` image is pre-configured with OpenBLAS

4. **Set NWChem environment for OpenBLAS:**
   ```dockerfile
   ENV BLASOPT="-lopenblas"
   ENV LAPACK_LIB="-lopenblas"
   ```

**Detection:**
- `ldd nwchem | grep mkl` shows Intel MKL linked (WILL FAIL on ARM64)
- Calculations are 10-100x slower than expected
- Numerical results differ from x86_64 reference calculations
- "Error: cannot find basis set" on ARM64 but not x86_64

**Which phase should address:** Phase 1 - BLAS configuration is fundamental to NWChem functionality.

**Sources:**
- [NWChem Compiling Documentation](https://nwchemgit.github.io/Compiling-NWChem.html) - BLASOPT configuration
- [OpenBLAS ARMV8 compilation](https://github.com/OpenMathLib/OpenBLAS/issues/2409)

---

### Pitfall 3: OpenBLAS Thread Creation Failures in Containers

**What goes wrong:** Container crashes with "OpenBLAS blas_thread_init: pthread_create failed for thread X of Y: Operation not permitted"

**Why it happens:**
- OpenBLAS tries to create threads based on detected CPU count
- Docker containers may have restricted seccomp profiles or cgroup limits
- ARM64 containers may have different thread limits than x86_64
- OpenBLAS on ARM64 may have SVE-related threading issues

**Consequences:**
- Container fails to start or crashes on first calculation
- Error: "OpenBLAS blas_thread_init: pthread_create failed"
- Error: "Resource temporarily unavailable"

**Prevention:**
1. **Set OPENBLAS_NUM_THREADS explicitly:**
   ```dockerfile
   ENV OPENBLAS_NUM_THREADS=1
   # Or match OMP_NUM_THREADS
   ENV OPENBLAS_NUM_THREADS=${OMP_NUM_THREADS:-4}
   ```

2. **Increase shared memory:**
   ```yaml
   # docker-compose.yml
   services:
     worker:
       shm_size: '512m'
   ```

3. **Set sysctl for overcommit (host-level):**
   ```bash
   sysctl vm.overcommit_memory=1
   ```

4. **Exclude problematic ARM64 kernels in custom OpenBLAS build:**
   ```bash
   # Avoid A64FX kernel which requires SVE
   make DYNAMIC_LIST='ARMV8 CORTEXA53 CORTEXA57 NEOVERSEN1'
   ```

**Detection:**
- Error message explicitly mentions "pthread_create failed"
- Container starts but crashes when first BLAS operation runs
- Works with `OPENBLAS_NUM_THREADS=1` but fails with higher values

**Which phase should address:** Phase 1 - Add OPENBLAS_NUM_THREADS to Dockerfile environment.

**Sources:**
- [OpenBLAS Issue #4825 - Building inside aarch64 docker fails](https://github.com/OpenMathLib/OpenBLAS/issues/4825)
- [Label Studio Issue #3070 - OpenBLAS pthread_create failed](https://github.com/HumanSignal/label-studio/issues/3070)

---

### Pitfall 4: Miniconda Installer Architecture Mismatch

**What goes wrong:** Python environment fails to install or runs in emulation mode, causing massive slowdowns.

**Why it happens:**
- Current Dockerfile hardcodes x86_64 Miniconda installer:
  ```dockerfile
  wget ... Miniconda3-py311_24.11.1-0-Linux-x86_64.sh
  ```
- Running x86_64 Python under QEMU is 5-20x slower
- pip packages may fail to install ARM64 wheels

**Consequences:**
- Python startup takes 10+ seconds instead of <1 second
- pip installs fail with compilation errors (no ARM64 wheels)
- RDKit, NumPy, SciPy all run in emulation (extremely slow)

**Prevention:**
1. **Use architecture-appropriate Miniconda:**
   ```dockerfile
   ARG TARGETARCH
   RUN if [ "$TARGETARCH" = "arm64" ]; then \
         wget https://repo.anaconda.com/miniconda/Miniconda3-py311_24.11.1-0-Linux-aarch64.sh -O /tmp/miniconda.sh; \
       else \
         wget https://repo.anaconda.com/miniconda/Miniconda3-py311_24.11.1-0-Linux-x86_64.sh -O /tmp/miniconda.sh; \
       fi \
       && bash /tmp/miniconda.sh -b -p $CONDA_DIR
   ```

2. **Or use miniforge which auto-detects:**
   ```dockerfile
   RUN curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" \
       && bash Miniforge3-*.sh -b -p $CONDA_DIR
   ```

3. **Verify native execution:**
   ```bash
   file $(which python)
   # Should show: ELF 64-bit LSB executable, ARM aarch64 (for ARM64)
   ```

**Detection:**
- `file $(which python)` shows x86_64 on ARM64 host
- Python operations are extremely slow
- pip shows "no matching distribution found" for packages with only x86_64 wheels

**Which phase should address:** Phase 1 - BLOCKER for ARM64 build.

**Sources:**
- [conda-forge miniforge](https://github.com/conda-forge/miniforge) - Multi-architecture installers
- [Anaconda archive](https://repo.anaconda.com/miniconda/) - Architecture-specific installers

---

## Moderate Pitfalls

Mistakes that cause degraded performance, difficult debugging, or require significant rework.

### Pitfall 5: Build Time Explosion with QEMU Emulation

**What goes wrong:** Building ARM64 images on x86_64 hosts takes hours instead of minutes.

**Why it happens:**
- Docker buildx uses QEMU to emulate ARM64 on x86_64
- Fortran compilation (NWChem, xTB, CREST) is CPU-intensive
- QEMU emulation is 5-20x slower than native
- Scientific computing builds involve heavy compilation

**Consequences:**
- CI/CD pipeline timeouts (builds exceed 6 hours)
- Developer frustration waiting for builds
- Temptation to skip ARM64 testing

**Prevention:**
1. **Use native ARM64 build runners:**
   - GitHub Actions: `runs-on: [self-hosted, linux, ARM64]`
   - AWS: Use Graviton instances for builds
   - Docker Build Cloud: Provides native ARM64 builders

2. **Use multi-stage builds to minimize emulated work:**
   ```dockerfile
   # Stage 1: Native build for Python deps (can be cached)
   FROM --platform=$BUILDPLATFORM python:3.11-slim AS python-deps
   RUN pip download --platform manylinux2014_aarch64 --only-binary :all: rdkit ...

   # Stage 2: Target platform for runtime
   FROM --platform=$TARGETPLATFORM ghcr.io/nwchemgit/nwchem-dev/arm64
   COPY --from=python-deps /wheels /wheels
   ```

3. **Use pre-built conda packages instead of compiling:**
   ```dockerfile
   # Instead of compiling xTB/CREST from source:
   RUN conda install -c conda-forge xtb crest
   ```

4. **Aggressive layer caching:**
   ```dockerfile
   # Cache computational chemistry tool installation separately
   COPY requirements-chemistry.txt .
   RUN conda install --file requirements-chemistry.txt
   # Then copy application code (changes more frequently)
   COPY src/ ./src/
   ```

**Detection:**
- `docker buildx build --platform linux/arm64` takes >30 minutes
- Build logs show heavy compilation activity under QEMU
- Native x86_64 build takes <5 minutes, ARM64 build takes >2 hours

**Which phase should address:** Phase 2 - CI/CD optimization. Phase 1 can use pre-built packages to avoid this.

**Sources:**
- [Docker Multi-Platform Builds](https://docs.docker.com/build/building/multi-platform/) - QEMU limitations
- [Faster Multi-Platform Builds Guide](https://www.docker.com/blog/faster-multi-platform-builds-dockerfile-cross-compilation-guide/)

---

### Pitfall 6: OpenMPI Fortran Module Path Issues on ARM64

**What goes wrong:** NWChem fails to compile or link with "cannot find mpi.mod" or CMake errors about MPI::MPI_Fortran paths.

**Why it happens:**
- OpenMPI installs Fortran `.mod` files in architecture-specific paths
- Cross-compilation or ARM64 builds may reference wrong paths (x86_64 paths)
- Bug reported: "Imported target 'MPI::MPI_Fortran' includes non-existent path '/usr/lib/x86_64-linux-gnu/fortran/...'"

**Consequences:**
- Build fails with missing module errors
- Linking fails with undefined MPI symbols
- Works on x86_64, fails on ARM64

**Prevention:**
1. **Use official pre-built NWChem ARM64 image** (already has correct MPI configuration)

2. **If building from source, verify MPI paths:**
   ```bash
   mpif90 --showme:compile
   mpif90 --showme:link
   # Verify paths exist and are for correct architecture
   ```

3. **Set explicit MPI paths:**
   ```dockerfile
   ENV MPI_INCLUDE_PATH=/usr/lib/aarch64-linux-gnu/openmpi/include
   ENV MPI_LIB_PATH=/usr/lib/aarch64-linux-gnu/openmpi/lib
   ```

4. **Ensure gfortran is installed before OpenMPI:**
   ```dockerfile
   RUN apt-get install -y gfortran && apt-get install -y libopenmpi-dev
   ```

**Detection:**
- CMake error: "MPI::MPI_Fortran includes non-existent path"
- Fortran compilation error: "cannot find mpi.mod"
- Build succeeds on x86_64, fails on ARM64

**Which phase should address:** Phase 1 if building from source; N/A if using pre-built images.

**Sources:**
- [OpenMPI Issue #12600 - Fortran .mod files location](https://github.com/open-mpi/ompi/issues/12600)
- [Debian Bug #1086087 - ARM64 MPI_Fortran path issue](https://www.mail-archive.com/debian-science-maintainers@alioth-lists.debian.net/msg60538.html)

---

### Pitfall 7: Performance Regression with Generic ARM64 Kernels

**What goes wrong:** ARM64 calculations run significantly slower than expected despite native execution.

**Why it happens:**
- OpenBLAS generic ARMV8 kernel lacks optimizations for specific ARM chips
- DGEMM performance degraded in OpenBLAS 0.3.19+ for generic ARMV8
- SVE (Scalable Vector Extension) optimizations may not be available
- NEON optimizations not used if compiled with wrong flags

**Consequences:**
- DFT calculations 2-5x slower than optimal
- Users perceive ARM64 as "slow" when it's just misconfigured
- Cost savings from ARM64 instances negated by longer runtime

**Prevention:**
1. **Specify optimal TARGET for deployment hardware:**
   ```bash
   # For AWS Graviton2 (Neoverse N1):
   export TARGET=NEOVERSEN1

   # For AWS Graviton3 (Neoverse V1):
   export TARGET=NEOVERSEVERSION

   # For Apple Silicon (via Rosetta not recommended, but):
   # Use native macOS build instead of Docker
   ```

2. **Benchmark actual performance:**
   ```python
   # Add benchmark script to validate performance
   # Compare ARM64 vs x86_64 for reference molecule
   ```

3. **Use OpenBLAS version < 0.3.19 if targeting generic ARMV8:**
   ```dockerfile
   # Or use conda-forge which handles version selection
   RUN conda install -c conda-forge openblas
   ```

**Detection:**
- Same molecule takes >2x longer on ARM64 vs comparable x86_64
- `lscpu` shows specific ARM chip, but OpenBLAS uses generic kernel
- Logs show "Using generic ARMV8 kernel"

**Which phase should address:** Phase 2 - Performance optimization after basic functionality.

**Sources:**
- [OpenBLAS Issue #3581 - ARMv8 DGEMM performance](https://github.com/OpenMathLib/OpenBLAS/issues/3581)
- [OpenBLAS ARM64 targets](https://www.openblas.net/Changelog.txt)

---

### Pitfall 8: RDKit ARM64 Wheel Availability

**What goes wrong:** RDKit installation fails on ARM64 with compilation errors.

**Why it happens:**
- RDKit has complex C++ dependencies (Boost, Eigen)
- Pre-built wheels may not exist for linux-aarch64
- Building from source requires extensive build environment

**Consequences:**
- pip install fails with "no matching distribution"
- Build from source fails with missing dependencies
- Multi-hour build times if dependencies need compilation

**Prevention:**
1. **Use conda-forge for RDKit on ARM64:**
   ```dockerfile
   # RDKit has ARM64 builds on conda-forge
   RUN conda install -c conda-forge rdkit
   ```

2. **Verify wheel availability before pip install:**
   ```bash
   pip download --platform manylinux2014_aarch64 --only-binary :all: rdkit
   # If this fails, use conda instead
   ```

3. **Use miniforge base image (conda-forge default channel):**
   ```dockerfile
   FROM condaforge/miniforge3:latest
   RUN conda install rdkit
   ```

**Detection:**
- `pip install rdkit` fails with "no matching distribution found"
- Build attempts show missing Boost/Eigen headers
- Works on x86_64, fails on ARM64

**Which phase should address:** Phase 1 - RDKit is core dependency.

**Sources:**
- [conda-forge RDKit](https://anaconda.org/conda-forge/rdkit) - Multi-platform builds
- [Python Speed - Docker build problems on Mac](https://pythonspeed.com/articles/docker-build-problems-mac/)

---

## Minor Pitfalls

Annoyances that are easily fixable but waste time if not anticipated.

### Pitfall 9: Architecture Detection in Multi-Stage Builds

**What goes wrong:** Build scripts fail because architecture detection happens at wrong stage.

**Why it happens:**
- `uname -m` returns host architecture in some build stages
- `$TARGETARCH` and `$BUILDPLATFORM` behave differently in different contexts
- Shell scripts may not receive correct architecture variable

**Prevention:**
1. **Use Docker BuildKit ARG variables:**
   ```dockerfile
   ARG TARGETPLATFORM
   ARG TARGETARCH
   ARG BUILDPLATFORM

   RUN echo "Building for $TARGETARCH on $BUILDPLATFORM"
   ```

2. **Verify at runtime:**
   ```dockerfile
   RUN arch && uname -m  # Should match TARGETARCH
   ```

**Which phase should address:** Phase 1 - Part of Dockerfile setup.

---

### Pitfall 10: Container Registry Multi-Arch Manifest Issues

**What goes wrong:** Pulling image on ARM64 gets x86_64 variant, or vice versa.

**Why it happens:**
- Multi-arch manifest not properly created with buildx
- Missing `--push` or incorrect manifest list
- Registry doesn't support multi-arch manifests

**Prevention:**
1. **Use buildx with proper manifest creation:**
   ```bash
   docker buildx build --platform linux/amd64,linux/arm64 \
     --tag myimage:latest --push .
   ```

2. **Verify manifest:**
   ```bash
   docker manifest inspect myimage:latest
   # Should show both amd64 and arm64 variants
   ```

3. **Use GitHub Container Registry (ghcr.io) which supports multi-arch**

**Which phase should address:** Phase 2 - CI/CD pipeline setup.

---

## Phase-Specific Warning Summary

| Phase | Topic | Pitfall | Priority |
|-------|-------|---------|----------|
| Phase 1 | Binary Selection | SIGILL from x86 SIMD (#1) | **BLOCKER** |
| Phase 1 | BLAS Configuration | Wrong OpenBLAS for ARM64 (#2) | **BLOCKER** |
| Phase 1 | Python Environment | Miniconda architecture mismatch (#4) | **BLOCKER** |
| Phase 1 | Threading | OpenBLAS thread creation failure (#3) | HIGH |
| Phase 1 | Dependencies | RDKit ARM64 wheels (#8) | HIGH |
| Phase 1 | Build | Fortran/MPI path issues (#6) | MEDIUM (if building from source) |
| Phase 2 | CI/CD | Build time explosion (#5) | HIGH |
| Phase 2 | Performance | Generic ARM64 kernel slowness (#7) | MEDIUM |
| Phase 2 | Deployment | Multi-arch manifest issues (#10) | MEDIUM |
| Phase 2 | Build | Architecture detection (#9) | LOW |

---

## Recommended Approach: Dockerfile.worker.arm64

Given the pitfalls above, the recommended approach for ARM64 support:

### Option A: Separate ARM64 Dockerfile (Recommended)

Create `Dockerfile.worker.arm64` that:
1. Uses `ghcr.io/nwchemgit/nwchem-dev/arm64` as base
2. Uses Miniforge for Python (auto-detects ARM64)
3. Installs xTB/CREST via conda-forge (has ARM64 builds)
4. Installs RDKit via conda-forge
5. Sets OPENBLAS_NUM_THREADS explicitly

```dockerfile
# Dockerfile.worker.arm64
FROM ghcr.io/nwchemgit/nwchem-dev/arm64:latest

# Use Miniforge for ARM64-native Python
RUN curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-aarch64.sh" \
    && bash Miniforge3-Linux-aarch64.sh -b -p /opt/conda \
    && rm Miniforge3-Linux-aarch64.sh

ENV PATH="/opt/conda/bin:$PATH"

# Install computational chemistry tools via conda-forge (has ARM64 builds)
RUN conda install -c conda-forge xtb crest rdkit

# Prevent OpenBLAS thread issues
ENV OPENBLAS_NUM_THREADS=4

# ... rest of configuration
```

### Option B: Multi-Arch Dockerfile with ARG

Single Dockerfile with architecture conditionals:

```dockerfile
ARG TARGETARCH

FROM ghcr.io/nwchemgit/nwchem-dev/${TARGETARCH}:latest

ARG TARGETARCH
RUN if [ "$TARGETARCH" = "arm64" ]; then \
      curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-aarch64.sh"; \
    else \
      curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"; \
    fi \
    && bash Miniforge3-*.sh -b -p /opt/conda \
    && rm Miniforge3-*.sh
```

---

## Validation Checklist for ARM64 Container

Before deploying ARM64 container:

### Binary Verification
- [ ] `file /opt/xtb/bin/xtb` shows "ARM aarch64" (not x86-64)
- [ ] `file /opt/crest/crest` shows "ARM aarch64"
- [ ] `file $(which nwchem)` shows "ARM aarch64"
- [ ] `file $(which python)` shows "ARM aarch64"

### BLAS Verification
- [ ] `ldd $(which nwchem) | grep blas` shows libopenblas (not libmkl)
- [ ] `echo $OPENBLAS_NUM_THREADS` is set
- [ ] Small matrix test completes without error

### Functional Tests
- [ ] Simple xTB optimization completes
- [ ] CREST conformer search completes
- [ ] NWChem DFT calculation completes
- [ ] RDKit molecule parsing works

### Performance Baseline
- [ ] Reference molecule calculation time within 2x of x86_64 time

---

## Sources Summary

### Official Documentation
- [NWChem Containers](https://nwchemgit.github.io/Containers.html) - ARM64 image availability
- [NWChem Compiling Guide](https://nwchemgit.github.io/Compiling-NWChem.html) - BLAS configuration
- [Docker Multi-Platform Builds](https://docs.docker.com/build/building/multi-platform/) - Buildx and QEMU
- [conda-forge ARM64 support](https://conda-forge.org/blog/2020/10/29/macos-arm64/)

### Issue Trackers (Real-World Problems)
- [Docker for Mac #7137 - Rosetta breaks amd64 images](https://github.com/docker/for-mac/issues/7137)
- [Docker for Mac #6620 - QEMU AVX support](https://github.com/docker/for-mac/issues/6620)
- [OpenBLAS #4825 - aarch64 Docker build fails](https://github.com/OpenMathLib/OpenBLAS/issues/4825)
- [OpenBLAS #3581 - ARMv8 DGEMM performance](https://github.com/OpenMathLib/OpenBLAS/issues/3581)
- [OpenMPI #12600 - Fortran module paths](https://github.com/open-mpi/ompi/issues/12600)

### Package Repositories
- [conda-forge xTB](https://anaconda.org/conda-forge/xtb) - linux-aarch64 builds
- [conda-forge CREST](https://anaconda.org/conda-forge/crest) - linux-aarch64 builds
- [conda-forge NWChem](https://anaconda.org/conda-forge/nwchem) - linux-aarch64 builds
- [conda-forge miniforge](https://github.com/conda-forge/miniforge) - ARM64 installer
