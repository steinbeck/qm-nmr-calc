# Domain Pitfalls: Containerizing qm-nmr-calc

**Domain:** Containerized scientific computing (NWChem/CREST/xTB) with task queue
**Researched:** 2026-02-02
**Overall Confidence:** MEDIUM-HIGH

This document catalogs pitfalls specific to containerizing the qm-nmr-calc application, which combines:
- Large computational chemistry binaries (NWChem ~500MB)
- Semi-empirical tools (CREST/xTB) with OpenMP requirements
- SQLite-based task queue (Huey)
- Filesystem-based job storage with GB-scale scratch files
- MPI subprocess execution

---

## Critical Pitfalls

Mistakes that cause application failure, data loss, or require architectural rewrites.

### Pitfall 1: SQLite Locking Failures in Multi-Container Deployments

**What goes wrong:** Huey uses SqliteHuey with `./data/huey.db`. When running API container and worker container(s) accessing the same SQLite database on a shared volume, "database is locked" errors occur randomly under load.

**Why it happens:** SQLite locks the entire database during writes. With multiple processes (API enqueueing, worker dequeuing, signal handlers updating status), lock contention is inevitable. The problem is compounded when:
- Docker volume is mounted over NFS/SMB (network filesystems break POSIX locking)
- Multiple worker processes/threads access concurrently
- Transactions are not short enough

**Consequences:**
- Jobs stuck in "queued" state
- Status updates lost (job completes but status shows "running")
- Worker crashes on repeated lock failures
- Data corruption in edge cases

**Prevention:**
1. **Phase 1 (MVP):** Keep SQLite for single-worker deployments only. Configure `fsync=True` (already done) and ensure WAL mode with busy_timeout:
   ```python
   huey = SqliteHuey('qm-nmr-calc', filename='./data/huey.db',
                     fsync=True, timeout=60)  # 60 second busy_timeout
   ```
2. **Phase 2 (Production):** Migrate to RedisHuey for multi-container deployments. Redis handles concurrent access properly and is Huey's recommended production backend.
3. **Never:** Mount SQLite database on network filesystems (NFS, SMB, CIFS)

**Detection (warning signs):**
- Intermittent `sqlite3.OperationalError: database is locked` in logs
- Jobs disappearing from queue without completion
- Status discrepancies between API response and actual job state

**Which phase should address:** Phase 1 (single-container) can use SQLite with timeout configuration. Phase 2 (multi-container) MUST use Redis.

**Sources:**
- [Database locked error with SqliteHuey - GitHub Issue #445](https://github.com/coleifer/huey/issues/445)
- [SQLite Locking Documentation](https://sqlite.org/lockingv3.html)
- [Huey Storage Guide](https://huey.readthedocs.io/en/latest/guide.html)

---

### Pitfall 2: NWChem Shared Memory Exhaustion

**What goes wrong:** NWChem fails with cryptic MPI errors or segmentation faults when running parallel calculations inside containers.

**Why it happens:** Docker containers have a default `/dev/shm` size of 64MB. NWChem's MPI implementation (and OpenMPI generally) uses shared memory for inter-process communication. With 4 MPI processes doing DFT calculations, 64MB is grossly insufficient.

**Consequences:**
- Immediate segfault on calculation start
- Cryptic MPI errors: `prte_init() failed`
- Calculations that work locally fail in container
- Intermittent failures based on molecule size

**Prevention:**
1. **Always** specify `--shm-size 256m` minimum (or `512m` for larger molecules):
   ```yaml
   # docker-compose.yml
   services:
     worker:
       shm_size: '512m'
   ```
2. Set in Dockerfile comments as mandatory requirement
3. Test with largest expected molecule during container validation

**Detection:**
- Segfaults immediately after "mpirun" in logs
- Error messages containing "shm" or "shared memory"
- Calculations that worked locally fail in container

**Which phase should address:** Phase 1 - must be configured from first container build.

**Sources:**
- [NWChem Containers Wiki](https://github.com/nwchemgit/nwchem/wiki/Containers)
- [NWChem Docker run examples](https://nwchemgit.github.io/Containers.html)

---

### Pitfall 3: MPI Process Binding Failures in Containers

**What goes wrong:** NWChem's `mpirun` commands fail or hang when process binding is enabled (default in OpenMPI 5.x).

**Why it happens:** OpenMPI tries to bind processes to specific CPU cores, but Docker containers report CPU limits differently than the host. OpenMPI 5.x has known issues detecting container CPU restrictions properly. The current code uses `--bind-to none` which is correct, but this can be lost during Dockerfile configuration.

**Consequences:**
- mpirun hangs indefinitely
- Only subset of requested processes start
- Error: "There are not enough slots available in the system"

**Prevention:**
1. **Always** use `--bind-to none` for mpirun in containers (already in code)
2. **Verify** OpenMPI version compatibility - OpenMPI 5.0.x has documented issues
3. Consider using OpenMPI 4.x in container if 5.x causes problems
4. Set `--oversubscribe` if running more processes than container CPU limit

Current code pattern to preserve:
```python
cmd = f"mpirun --bind-to none -n {processes} nwchem {infile} > {outfile}"
```

**Detection:**
- mpirun exits immediately with code 1
- Processes hang after "Starting NWChem"
- stderr contains "not enough slots" or binding errors

**Which phase should address:** Phase 1 - validate MPI configuration in first container.

**Sources:**
- [OpenMPI Docker Issues #12431](https://github.com/open-mpi/ompi/issues/12431)
- [OpenMPI 5.0.x binding issues #12967](https://github.com/open-mpi/ompi/issues/12967)

---

### Pitfall 4: Scratch Directory Disk Exhaustion

**What goes wrong:** NWChem calculations fail mid-run due to "No space left on device" when scratch files exceed container disk limits.

**Why it happens:** A single DFT calculation on a moderate molecule (20-30 atoms) can generate 1-5GB of scratch files. Conformer ensemble calculations multiply this by conformer count. Default Docker container disk space or `/tmp` sizing is often 10GB or less.

**Consequences:**
- Calculations fail mid-run (wasted compute time)
- Partial results with corrupted output files
- Container becomes unresponsive if root filesystem fills
- Job status stuck at "running" forever

**Prevention:**
1. **Mount dedicated scratch volume** with adequate space:
   ```yaml
   volumes:
     - scratch_data:/app/data/jobs  # Persistent, sized appropriately
     - type: tmpfs
       target: /scratch
       tmpfs:
         size: 20G  # For active calculation scratch
   ```
2. **Configure NWChem SCRATCH_DIR** to dedicated volume
3. **Implement scratch cleanup** after calculation completion (already done in `storage.py`)
4. **Monitor disk usage** and alert before exhaustion
5. Size scratch volume based on: `max_concurrent_jobs * max_conformers * 5GB`

**Detection:**
- "No space left on device" in NWChem output
- Calculations that succeed on small molecules fail on larger ones
- `df -h` shows /tmp or container overlay at 100%

**Which phase should address:** Phase 1 - configure volume sizing from initial deployment.

**Sources:**
- [NWChem Scratch Directory Documentation](https://github.com/nwchemgit/nwchem/wiki/Scratch_Dir)
- [TACC I/O Best Practices](https://docs.tacc.utexas.edu/tutorials/managingio/)

---

### Pitfall 5: CREST/xTB OpenMP Stack Overflow

**What goes wrong:** CREST conformer generation crashes with segmentation faults on larger molecules.

**Why it happens:** xTB uses OpenMP for parallelization and requires large stack sizes for molecular calculations. Default container stack limits (8MB) are insufficient for molecules with 50+ atoms.

**Consequences:**
- Silent segfault during CREST run
- Truncated or missing conformer output
- Inconsistent failures (depends on molecule size)
- No useful error message - just "Killed" or exit code 139

**Prevention:**
1. **Set OpenMP environment variables in container:**
   ```dockerfile
   ENV OMP_STACKSIZE=2G
   ENV OMP_NUM_THREADS=4,1
   ENV OMP_MAX_ACTIVE_LEVELS=1
   ENV GFORTRAN_UNBUFFERED_ALL=1
   ```
2. **Set ulimit in container:**
   ```dockerfile
   # In entrypoint script
   ulimit -s unlimited
   ```
3. **Verify settings** are actually applied (check from within running container)

Current code already sets these in `run_crest()`:
```python
env["OMP_STACKSIZE"] = "2G"
env["GFORTRAN_UNBUFFERED_ALL"] = "1"
```

**Detection:**
- Exit code 139 (SIGSEGV) from CREST
- Segfaults only on larger molecules
- `dmesg` shows "Out of memory" for stack

**Which phase should address:** Phase 1 - must be in base Dockerfile.

**Sources:**
- [xTB Setup Documentation](https://xtb-docs.readthedocs.io/en/latest/setup.html)
- [xTB Stack Overflow Issue #191](https://github.com/grimme-lab/xtb/issues/191)

---

## Moderate Pitfalls

Mistakes that cause degraded performance, debugging difficulty, or technical debt.

### Pitfall 6: Job State Lost on Container Restart

**What goes wrong:** Running jobs are lost when worker container restarts, leaving orphaned "running" status in database and no way to recover.

**Why it happens:** The current architecture stores job status as "running" in `status.json` and relies on Huey signals for completion. If worker dies mid-calculation:
- Huey signal handlers (`SIGNAL_INTERRUPTED`) may not fire
- Job status remains "running" forever
- No mechanism to detect orphaned jobs on startup

**Consequences:**
- Jobs stuck in "running" state after restart
- Users wait forever for results that will never come
- Manual intervention required to reset job status
- Scratch files from interrupted jobs accumulate

**Prevention:**
1. **Implement orphan detection on worker startup:**
   ```python
   def detect_orphaned_jobs():
       """Mark jobs still 'running' from previous worker as failed."""
       running_jobs = list_jobs_by_status('running')
       for job_id in running_jobs:
           update_job_status(job_id, status='failed',
                           error_message='Worker restarted - job interrupted')
   ```
2. **Store worker PID with job** to identify true owner
3. **Implement job heartbeat** for long-running calculations
4. **Use proper restart policy** in Docker:
   ```yaml
   restart: unless-stopped
   ```

**Detection:**
- Jobs with "running" status but no active worker process
- Scratch directories for jobs with no active calculation
- User complaints about stuck jobs after deployments

**Which phase should address:** Phase 2 - when implementing production reliability features.

---

### Pitfall 7: Hardcoded Paths Break in Container

**What goes wrong:** Application fails because it expects paths like `./data/jobs` that don't exist in container filesystem.

**Why it happens:** Current code uses relative paths:
```python
DATA_DIR = Path("./data/jobs")
huey = SqliteHuey('qm-nmr-calc', filename='./data/huey.db', ...)
```

These assume working directory context that may differ in container.

**Consequences:**
- FileNotFoundError on startup
- Jobs written to wrong location (container overlay, not volume)
- Data loss on container restart

**Prevention:**
1. **Use environment variables for all data paths:**
   ```python
   DATA_DIR = Path(os.environ.get('QM_DATA_DIR', './data/jobs'))
   HUEY_DB = os.environ.get('QM_HUEY_DB', './data/huey.db')
   ```
2. **Set WORKDIR in Dockerfile:**
   ```dockerfile
   WORKDIR /app
   ```
3. **Validate paths exist on startup:**
   ```python
   def validate_paths():
       DATA_DIR.mkdir(parents=True, exist_ok=True)
   ```
4. **Document required volume mounts** clearly

**Detection:**
- FileNotFoundError or PermissionError on startup
- Data written to unexpected locations
- Volume appears empty despite job completions

**Which phase should address:** Phase 1 - must be configurable from first container.

---

### Pitfall 8: NWChem Binary Not Found / Wrong Version

**What goes wrong:** Container reports "nwchem not found" or uses wrong NWChem version.

**Why it happens:** NWChem is notoriously difficult to compile. Pre-built binaries may:
- Not be in PATH in container context
- Be compiled for wrong architecture (amd64 vs arm64)
- Have missing runtime dependencies (OpenMPI, BLAS)
- Be outdated version without needed features

**Consequences:**
- Application exits on startup (current `validate_nwchem()` behavior)
- Subtle calculation differences with different versions
- Runtime crashes if dependencies missing

**Prevention:**
1. **Use official NWChem Docker images as base:**
   ```dockerfile
   FROM ghcr.io/nwchemgit/nwchem-dev/amd64 AS nwchem
   # Or use pre-built binary
   ```
2. **Validate NWChem on container startup** (already implemented)
3. **Pin specific NWChem version** and document compatibility
4. **Test with actual calculations**, not just `which nwchem`

**Detection:**
- "FATAL: nwchem not found in PATH" on startup
- Calculation results differ from development environment
- Missing library errors in mpirun output

**Which phase should address:** Phase 1 - critical for first working container.

**Sources:**
- [NWChem Docker Images](https://github.com/nwchemgit/nwchem-dockerfiles)
- [Official NWChem Container Registry](https://hub.docker.com/r/nwchemorg/nwchem-dev)

---

### Pitfall 9: Container Image Size Explosion

**What goes wrong:** Docker image becomes 5-10GB, making deployments slow and expensive.

**Why it happens:** Naive Dockerfile includes:
- Full NWChem source + build artifacts
- Development dependencies
- Multiple Python virtual environments
- Unstripped binaries and debug symbols

**Consequences:**
- Slow CI/CD pipelines (10+ minutes to push/pull)
- Expensive container registry storage
- Slow container startup on fresh nodes
- Exceeds some registry size limits

**Prevention:**
1. **Use multi-stage builds:**
   ```dockerfile
   FROM python:3.11-slim AS python-deps
   # Install Python deps only

   FROM ghcr.io/nwchemgit/nwchem-dev/amd64 AS nwchem-base
   # Copy only runtime files from python-deps
   ```
2. **Clean up in same layer as install:**
   ```dockerfile
   RUN apt-get update && apt-get install -y ... \
       && rm -rf /var/lib/apt/lists/*
   ```
3. **Use .dockerignore** to exclude:
   ```
   .git
   .venv
   __pycache__
   tests/
   *.pyc
   data/jobs/*
   ```
4. **Strip unnecessary NWChem modules** if using custom build

**Detection:**
- `docker images` shows >3GB for production image
- CI pipeline timeouts during push
- Deployment time dominated by image pull

**Which phase should address:** Phase 1 - optimize from first Dockerfile.

---

### Pitfall 10: Environment Variable Conflicts

**What goes wrong:** Container's environment variables conflict with NWChem/CREST/xTB requirements.

**Why it happens:** Multiple tools need specific environment settings:
- NWChem: needs `NWCHEM_TOP`, `NWCHEM_BASIS_LIBRARY`
- CREST: needs `OMP_STACKSIZE`, `OMP_NUM_THREADS`
- xTB: needs `XTBPATH`
- OpenMPI: needs `OMPI_MCA_*` settings

Parent images or orchestration systems may set conflicting values.

**Consequences:**
- NWChem can't find basis sets
- CREST uses wrong thread count
- MPI binding conflicts
- Hard-to-debug intermittent failures

**Prevention:**
1. **Explicitly set all required env vars in Dockerfile:**
   ```dockerfile
   ENV NWCHEM_BASIS_LIBRARY=/usr/share/nwchem/libraries/
   ENV OMP_STACKSIZE=2G
   ENV OMP_NUM_THREADS=4,1
   ENV OMPI_MCA_btl_base_warn_component_unused=0
   ```
2. **Override at runtime for flexibility:**
   ```yaml
   environment:
     - OMP_NUM_THREADS=${OMP_NUM_THREADS:-4,1}
   ```
3. **Document all required environment variables**
4. **Test in clean container** (not inheriting local env)

**Detection:**
- "Error: cannot find basis set" from NWChem
- Calculations use wrong CPU count
- Different behavior between local and container

**Which phase should address:** Phase 1 - must be correct in first Dockerfile.

---

## Minor Pitfalls

Annoyances that are easily fixable but waste time if not anticipated.

### Pitfall 11: Container Runs as Root

**What goes wrong:** Files created in mounted volumes are owned by root, causing permission issues on host.

**Why it happens:** Docker containers run as root by default. When writing to mounted volumes, files inherit root ownership.

**Consequences:**
- Cannot delete job files without sudo on host
- Permission denied when app tries to read files created by host
- Security scanners flag container as vulnerable

**Prevention:**
1. **Run as non-root user:**
   ```dockerfile
   RUN useradd -m appuser
   USER appuser
   ```
2. **Match host UID for development:**
   ```yaml
   user: "${UID}:${GID}"
   ```
3. **Set proper permissions on volume directories**

**Which phase should address:** Phase 1 - security best practice.

---

### Pitfall 12: Log Output Lost

**What goes wrong:** Container logs don't capture NWChem/CREST output, making debugging impossible.

**Why it happens:** Current code redirects subprocess output to files:
```python
cmd = f"... nwchem {infile} > {outfile} 2> {logfile}"
```
These files are in scratch directories that may be cleaned up or inside the container.

**Prevention:**
1. **Also capture to container stdout for Docker logging:**
   ```python
   # Tee output to both file and stdout
   cmd = f"... nwchem {infile} 2>&1 | tee {outfile}"
   ```
2. **Mount log directory as volume**
3. **Implement log aggregation** for production
4. **Set appropriate log retention**

**Which phase should address:** Phase 2 - when implementing observability.

---

### Pitfall 13: Time Zone Issues

**What goes wrong:** Job timestamps are wrong or inconsistent.

**Why it happens:** Containers default to UTC. Host may be in different timezone. Current code uses `datetime.utcnow()` which is correct, but display/comparison issues can arise.

**Prevention:**
1. **Always use UTC internally** (already done)
2. **Set TZ in container explicitly:**
   ```dockerfile
   ENV TZ=UTC
   ```
3. **Convert to user timezone only for display**

**Which phase should address:** Phase 1 - simple fix in Dockerfile.

---

## Phase-Specific Warning Summary

| Phase | Topic | Likely Pitfall | Priority |
|-------|-------|---------------|----------|
| Phase 1 | Base Container | NWChem binary not found (#8) | BLOCKER |
| Phase 1 | Base Container | Shared memory exhaustion (#2) | BLOCKER |
| Phase 1 | Base Container | MPI binding failures (#3) | BLOCKER |
| Phase 1 | Base Container | OpenMP stack overflow (#5) | HIGH |
| Phase 1 | Base Container | Hardcoded paths (#7) | HIGH |
| Phase 1 | Base Container | Environment conflicts (#10) | MEDIUM |
| Phase 1 | Base Container | Image size (#9) | MEDIUM |
| Phase 1 | Base Container | Root user (#11) | LOW |
| Phase 1 | Storage | Scratch disk exhaustion (#4) | HIGH |
| Phase 2 | Multi-container | SQLite locking (#1) | BLOCKER |
| Phase 2 | Reliability | Orphaned jobs on restart (#6) | HIGH |
| Phase 2 | Observability | Lost logs (#12) | MEDIUM |

---

## Recommended Testing Checklist

Before declaring container deployment ready:

### Functional Tests
- [ ] NWChem validation passes on startup
- [ ] Small molecule calculation completes successfully
- [ ] Large molecule (50+ atoms) doesn't crash with stack overflow
- [ ] CREST conformer generation completes
- [ ] Multiple concurrent jobs don't corrupt each other
- [ ] Job survives worker container restart (status correctly updated)

### Resource Tests
- [ ] Calculation succeeds with only 512MB shm_size
- [ ] 20GB scratch space handles largest expected molecule
- [ ] Container runs successfully as non-root
- [ ] Image size under 3GB

### Integration Tests
- [ ] API container can enqueue jobs
- [ ] Worker container processes queue
- [ ] Status updates visible through API
- [ ] Volumes persist across container restarts

---

## Sources Summary

### Official Documentation
- [NWChem Containers Wiki](https://github.com/nwchemgit/nwchem/wiki/Containers) - Container configuration, shm_size
- [NWChem Scratch Directory](https://github.com/nwchemgit/nwchem/wiki/Scratch_Dir) - Scratch file management
- [xTB Setup Documentation](https://xtb-docs.readthedocs.io/en/latest/setup.html) - OpenMP configuration
- [Huey Documentation](https://huey.readthedocs.io/en/latest/guide.html) - Storage backend selection

### Issue Trackers (Real-World Problems)
- [Huey SqliteHuey locking #445](https://github.com/coleifer/huey/issues/445)
- [OpenMPI Docker issues #12431](https://github.com/open-mpi/ompi/issues/12431)
- [xTB Stack overflow #191](https://github.com/grimme-lab/xtb/issues/191)

### Best Practices
- [Docker Anti-Patterns - Codefresh](https://codefresh.io/blog/docker-anti-patterns/)
- [Docker Restart Policies](https://docs.docker.com/engine/containers/start-containers-automatically/)
- [SQLite WAL mode for concurrency](https://sqlite.org/lockingv3.html)
