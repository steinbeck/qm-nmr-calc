# Pitfalls Research: Async NMR Calculation Service

**Domain:** Async NMR QM calculation web service with ISiCLE/NWChem
**Researched:** 2026-01-19
**Overall Confidence:** MEDIUM (ISiCLE-specific details LOW, general patterns HIGH)

---

## ISiCLE/NWChem Pitfalls

### Critical: NWChem Memory Exhaustion

**What goes wrong:** NWChem calculations fail with cryptic errors like "ptsalloc: increase memory", "ga_create failed", or "ao_replicated: insufficient memory". These failures can occur mid-calculation, wasting hours of compute time.

**Why it happens:** NWChem has a complex memory model with three regions (stack, heap, global) that partition total memory. The default is 512 MB total, which is insufficient for medium-sized molecules. Memory is allocated *per processor core*, not total.

**Consequences:**
- Calculations abort after significant runtime
- No partial results saved
- Server resources consumed without output

**Warning signs:**
- Calculations for molecules >30 atoms fail unexpectedly
- Log files show memory allocation errors
- Jobs work on small test molecules but fail on real inputs

**Prevention:**
1. Always set explicit memory in NWChem input: `memory total 4 gb` (per core)
2. Use `semidirect memsize N filesize 0` for better memory efficiency
3. Calculate required memory based on molecule size and basis set
4. Implement pre-flight memory estimation before queuing jobs

**Detection:** Monitor NWChem stderr for memory-related warnings early in calculation.

**Phase to address:** Phase 1 (ISiCLE integration) - bake memory configuration into job templates.

**Confidence:** HIGH - [NWChem Memory Wiki](https://github.com/nwchemgit/nwchem/wiki/Memory), [NWChem FAQ](https://nwchemgit.github.io/FAQ.html)

---

### Critical: NWChem Scratch Disk Overflow

**What goes wrong:** Calculations fail because scratch directory fills up. Integral files can be *enormous* - 5-6 GB per benzene calculation with cc-pvtz basis, multiplied by number of cores.

**Why it happens:** NWChem writes temporary integral files to scratch directory during SCF. On parallel jobs, this scales with core count. A single conformer calculation can generate 50+ GB of temp files.

**Consequences:**
- Disk full crashes the entire VM (not just the job)
- Other jobs fail due to shared filesystem
- VM may become unresponsive

**Warning signs:**
- `*.aoints.*` and `*.gridpts.*` files growing rapidly
- Filesystem approaching capacity during calculation
- Jobs that worked before suddenly fail

**Prevention:**
1. Set dedicated scratch directory: `scratch_dir /tmp/nwchem_scratch`
2. Use direct SCF to minimize disk I/O: add `direct` to scf block
3. Implement per-job scratch directories that auto-cleanup
4. Set disk quotas or monitor filesystem usage
5. Add cleanup trap in job wrapper scripts

**Detection:** Monitor scratch directory size during calculations. Alert at 80% capacity.

**Phase to address:** Phase 1 - implement scratch directory management from the start.

**Confidence:** HIGH - [NWChem Scratch_Dir Documentation](https://nwchemgit.github.io/Scratch_Dir.html)

---

### Critical: ISiCLE CREST Conformer Generation Failures

**What goes wrong:** CREST (conformer sampling step) crashes with segmentation faults, especially for larger molecules. The `-nocross` workaround may miss important conformers.

**Why it happens:** CREST's GC (genetic crossing) algorithm can trigger memory issues, particularly in Conda-installed versions. Reported with CREST 2.12.

**Consequences:**
- Entire NMR workflow fails at early stage
- Workarounds may produce incomplete conformer sets
- Unreliable results for certain molecule classes

**Warning signs:**
- SIGSEGV errors in logs
- Crashes "sometime prior to (or at the beginning of) the GC algorithm"
- Works on small molecules, fails on larger ones

**Prevention:**
1. Test ISiCLE installation thoroughly with various molecule sizes
2. Consider `-nocross` flag as fallback (document limitation)
3. Implement conformer count validation (did CREST produce reasonable output?)
4. Have timeout + retry logic for CREST step

**Detection:** Parse CREST output for expected conformer count. Zero conformers = failure.

**Phase to address:** Phase 1 - test matrix for conformer generation.

**Confidence:** MEDIUM - [CREST GitHub Issue #203](https://github.com/crest-lab/crest/issues/203)

---

### Moderate: NWChem AUTOZ Coordinate Failure

**What goes wrong:** NWChem fails to generate internal coordinates for certain molecular geometries, particularly linear molecules with 4+ atoms in a chain.

**Why it happens:** AUTOZ automatic coordinate generation cannot handle strictly linear arrangements.

**Consequences:** Job fails at geometry setup before any calculation begins.

**Prevention:**
1. Use `geometry NOAUTOZ` in input files to skip automatic coordinate generation
2. Use Cartesian coordinates consistently
3. Validate molecular geometry before submission

**Detection:** Check NWChem output for AUTOZ warnings in first 100 lines.

**Phase to address:** Phase 1 - use NOAUTOZ by default in ISiCLE wrapper.

**Confidence:** HIGH - [NWChem FAQ](https://nwchemgit.github.io/FAQ.html)

---

### Moderate: ISiCLE API Instability

**What goes wrong:** ISiCLE's Python API has documented issues including missing attributes (XTBWrapper missing `.xyz`), unexpected TypeErrors in example scripts.

**Why it happens:** ISiCLE development has slowed; API surface not fully hardened.

**Consequences:**
- Wrapper code needs defensive handling
- Version upgrades may break integration
- Error messages may be unhelpful

**Warning signs:**
- AttributeError on ISiCLE objects
- TypeError from ISiCLE functions
- Behavior differs from documentation

**Prevention:**
1. Pin ISiCLE version strictly
2. Build abstraction layer around ISiCLE calls
3. Comprehensive try/catch around all ISiCLE operations
4. Test suite covering full workflow before deployment

**Detection:** Integration tests with real molecules.

**Phase to address:** Phase 1 - abstraction layer around ISiCLE.

**Confidence:** HIGH - [ISiCLE GitHub Issues #47, #14](https://github.com/pnnl/isicle/issues)

---

### Minor: NWChem Basis Function Linear Dependency

**What goes wrong:** SCF convergence fails due to near-linearly-dependent basis functions.

**Why it happens:** Certain basis sets in certain molecular environments create mathematical instability.

**Consequences:** Calculation never converges, runs until timeout.

**Prevention:**
1. Allow NWChem to handle linear dependencies by NOT setting `lindep:n_dep 0`
2. Use well-tested functional/basis combinations for NMR (B3LYP/6-31G* as default)
3. Implement convergence monitoring

**Phase to address:** Phase 1 - sensible defaults in job templates.

**Confidence:** HIGH - [NWChem FAQ](https://nwchemgit.github.io/FAQ.html)

---

## Async Processing Pitfalls

### Critical: FastAPI BackgroundTasks Unsuitable for Long Jobs

**What goes wrong:** Using FastAPI's built-in `BackgroundTasks` for hour-long calculations causes jobs to be lost on server restart, no progress tracking, no retry capability.

**Why it happens:** BackgroundTasks runs in the same event loop as FastAPI. It's designed for quick tasks (email notifications), not multi-hour computations.

**Consequences:**
- Server restart kills all running calculations
- No way to check job progress
- No persistence of job state
- No retry on failure

**Warning signs:**
- Jobs disappear after server redeploy
- Cannot tell if job is running or stuck
- No recovery from transient failures

**Prevention:**
1. Use dedicated task queue (Celery, RQ, Huey, or lighter-weight custom solution)
2. Persist job state to filesystem/database before starting
3. Implement heartbeat mechanism for running jobs
4. Design jobs to be restartable from checkpoints

**Detection:** Test by restarting server during calculation.

**Phase to address:** Phase 1 (Core) - foundational architecture decision.

**Confidence:** HIGH - [FastAPI Background Tasks Documentation](https://fastapi.tiangolo.com/tutorial/background-tasks/), [Leapcell Blog](https://leapcell.io/blog/managing-background-tasks-and-long-running-operations-in-fastapi)

---

### Critical: Non-Idempotent Job Execution

**What goes wrong:** When a job fails and retries, it duplicates work, corrupts results, or enters inconsistent state.

**Why it happens:** Jobs that don't track their own progress can't safely resume. Partial results from failed runs may persist.

**Consequences:**
- Duplicate conformer calculations
- Corrupted result files
- Resource waste from repeated work
- Confused users seeing partial/wrong results

**Warning signs:**
- Same job shows different results on different runs
- Output files have unexpected content after retry
- Disk space grows unexpectedly

**Prevention:**
1. Make each job step idempotent (can safely re-run)
2. Track progress in job metadata file
3. Clean up partial outputs before retry
4. Use atomic file operations (write to temp, rename to final)

**Detection:** Test by killing job mid-execution and restarting.

**Phase to address:** Phase 1 - job state management design.

**Confidence:** HIGH - [Medium Article on Idempotency](https://medium.com/@vivekburman1997/data-engineering-part-1-idempotency-retry-and-recovery-b3631a9b8b6f)

---

### Critical: Zombie Process Accumulation

**What goes wrong:** Failed or timed-out NWChem processes leave zombie processes that consume PIDs and eventually prevent new jobs from running.

**Why it happens:** NWChem with MPI can exit uncleanly (MPI_Abort). Parent process (job runner) may not properly wait() on children.

**Consequences:**
- PID table fills up
- Server becomes unable to start new processes
- System slowdown

**Warning signs:**
- `ps aux | grep Z` shows growing zombie count
- New jobs fail to start with "cannot fork"
- System load appears low but jobs won't run

**Prevention:**
1. Proper subprocess handling: always wait() on child processes
2. Use process groups for NWChem + MPI to enable clean kills
3. Implement periodic zombie reaping
4. Set job timeouts with proper cleanup

**Detection:** Monitor zombie process count; alert if > 0.

**Phase to address:** Phase 1 - subprocess management wrapper.

**Confidence:** HIGH - [NREL ESIFHPC3 Issue #2](https://github.com/NREL/ESIFHPC3/issues/2)

---

### Moderate: No Timeout on Hung Calculations

**What goes wrong:** A calculation that will never converge runs forever, blocking resources and never returning results.

**Why it happens:** SCF can fail to converge. NWChem will keep trying until max iterations (which may be set very high).

**Consequences:**
- Job runs for days with no useful output
- Resources blocked from other jobs
- User waits indefinitely

**Prevention:**
1. Set explicit timeouts per calculation phase
2. Implement wall-clock timeout for entire job
3. NWChem: set reasonable `maxiter` in SCF/DFT blocks
4. Return "calculation did not converge" as valid result

**Detection:** Job runtime exceeding 2x expected time.

**Phase to address:** Phase 1 - timeout configuration.

**Confidence:** HIGH - [Python subprocess docs](https://docs.python.org/3/library/subprocess.html)

---

### Moderate: Concurrent Job Resource Contention

**What goes wrong:** Multiple simultaneous jobs compete for CPU, memory, and disk, causing all jobs to run slowly or fail.

**Why it happens:** Single VM with no job scheduling/limiting allows uncontrolled parallelism.

**Consequences:**
- All jobs slow down dramatically
- Memory pressure causes OOM kills
- Disk I/O saturation

**Prevention:**
1. Implement job concurrency limit (likely 1 for single VM)
2. Queue jobs with fair ordering
3. Reserve resources per job
4. Monitor resource usage per job

**Detection:** CPU at 100%, memory pressure, I/O wait.

**Phase to address:** Phase 1 - job queue with concurrency control.

**Confidence:** HIGH

---

## Input Validation Pitfalls

### Critical: SMILES Parsing Variations

**What goes wrong:** Valid SMILES in one toolkit fails to parse in another, or parses differently. User-submitted SMILES may be syntactically valid but semantically problematic.

**Why it happens:** SMILES specification has ambiguities. RDKit has specific sanitization rules. Aromatic atom notation is case-sensitive.

**Consequences:**
- Jobs fail at molecule loading
- Different results from "equivalent" inputs
- Hard to reproduce calculations

**Warning signs:**
- `Chem.MolFromSmiles()` returns None
- "non-ring atom marked aromatic" errors
- Different atom counts than expected

**Prevention:**
1. Validate SMILES with RDKit immediately on submission
2. Canonicalize SMILES before processing
3. Return clear error messages with specific parsing failures
4. Sanitize with explicit options:
   ```python
   params = Chem.SmilesParserParams()
   params.sanitize = True
   mol = Chem.MolFromSmiles(smiles, params)
   if mol is None:
       # Handle error
   ```

**Detection:** Validate at API boundary, before job enters queue.

**Phase to address:** Phase 1 - input validation layer.

**Confidence:** HIGH - [RDKit Blog on Sanitization](https://greglandrum.github.io/rdkit-blog/posts/2025-06-27-sanitization-and-file-parsing.html)

---

### Moderate: Molecule Size Limits

**What goes wrong:** User submits molecule that would take days/weeks to calculate, or requires more memory than available.

**Why it happens:** No pre-validation of molecule complexity vs. available resources.

**Consequences:**
- Jobs run for excessive time
- Memory exhaustion
- User frustrated waiting

**Prevention:**
1. Set atom count limits (e.g., max 100 atoms for draft, 50 for production quality)
2. Estimate calculation time from molecular properties
3. Warn user before accepting long jobs
4. Implement calculation tiers with different limits

**Detection:** Count atoms/heavy atoms before accepting job.

**Phase to address:** Phase 1 - validation with resource estimation.

**Confidence:** MEDIUM

---

### Minor: Invalid Charge States

**What goes wrong:** User provides SMILES with implicit charge that ISiCLE interprets incorrectly, or forgets to specify charge for charged molecules.

**Why it happens:** SMILES encodes some charge information, but formal charge must often be explicitly set for QM calculations.

**Consequences:**
- Wrong calculation (wrong number of electrons)
- Convergence failures
- Incorrect NMR predictions

**Prevention:**
1. Detect charged groups in input molecule
2. Prompt user to confirm charge state
3. Validate charge consistency between SMILES and explicit charge parameter

**Detection:** Check molecular charge during validation.

**Phase to address:** Phase 2 - enhanced validation.

**Confidence:** MEDIUM

---

## Resource Management Pitfalls

### Critical: Disk Space Exhaustion

**What goes wrong:** Job outputs, scratch files, and logs accumulate until disk is full, crashing all services.

**Why it happens:** Each calculation produces multiple output files. Failed jobs may leave orphaned files. No automatic cleanup.

**Consequences:**
- VM becomes unresponsive
- All jobs fail
- Service downtime

**Warning signs:**
- Filesystem > 80% full
- Old job directories not being cleaned
- Growing `/tmp` or scratch directories

**Prevention:**
1. Job-specific directories with automated cleanup on completion
2. Scratch directory cleanup on job end (success or failure)
3. Retention policy for completed jobs (auto-delete after N days)
4. Disk space monitoring with alerts
5. Trap signals to ensure cleanup on job kill

**Detection:** Monitor disk usage; alert at 70%, critical at 85%.

**Phase to address:** Phase 1 - storage management from day one.

**Confidence:** HIGH

---

### Moderate: Memory Leaks in Long-Running Service

**What goes wrong:** The API server or job worker leaks memory over time, eventually requiring restart.

**Why it happens:** Python garbage collection edge cases. ISiCLE/RDKit object lifecycle issues. Large molecule data kept in memory.

**Consequences:**
- Service becomes slow
- OOM kills lose jobs
- Unplanned restarts

**Prevention:**
1. Run calculations in separate processes (fork per job)
2. Explicit cleanup of large objects
3. Memory monitoring with automatic worker restart at threshold
4. Periodic scheduled restarts during quiet periods

**Detection:** Monitor worker process RSS over time.

**Phase to address:** Phase 2 - operational hardening.

**Confidence:** MEDIUM

---

### Moderate: Single Point of Failure

**What goes wrong:** Single VM means any failure takes down entire service.

**Why it happens:** Intentional simplification for v1 (single VM deployment).

**Consequences:**
- Hardware/network issues cause complete outage
- No capacity for burst load

**Prevention for v1:**
1. Implement graceful degradation (queue jobs even if worker is down)
2. Persist job state to survive restarts
3. Automatic service restart on failure
4. Health check endpoints for monitoring

**Prevention for future:**
- Queue with separate worker VM
- Filesystem storage on network mount

**Phase to address:** Phase 1 (mitigation), post-v1 (full solution).

**Confidence:** HIGH

---

## Prevention Strategies Summary

| Category | Key Strategy | Phase |
|----------|--------------|-------|
| ISiCLE/NWChem | Memory config + scratch management + abstraction layer | 1 |
| Async Processing | Proper task queue (not BackgroundTasks) + idempotent jobs | 1 |
| Input Validation | RDKit validation + size limits at API boundary | 1 |
| Resource Management | Per-job directories + cleanup + monitoring | 1 |

---

## Phase-Specific Warnings

| Phase | Topic | Likely Pitfall | Mitigation |
|-------|-------|---------------|------------|
| Phase 1 | ISiCLE setup | Installation/environment issues | Docker container with pinned deps |
| Phase 1 | Job queue | Wrong tool choice (BackgroundTasks) | Use proper task queue |
| Phase 1 | NWChem integration | Memory/disk failures | Configure limits upfront |
| Phase 2 | Scaling | Resource contention | Strict concurrency limits |
| Phase 2 | Reliability | No retry logic | Idempotent checkpointed jobs |
| Phase 3 | UI | Long-polling timeout | WebSocket or generous polling interval |

---

## Sources

### ISiCLE/NWChem
- [NWChem FAQ](https://nwchemgit.github.io/FAQ.html) - Common errors and solutions
- [NWChem Memory Wiki](https://github.com/nwchemgit/nwchem/wiki/Memory) - Memory configuration
- [NWChem Scratch_Dir](https://nwchemgit.github.io/Scratch_Dir.html) - Scratch directory management
- [ISiCLE GitHub Issues](https://github.com/pnnl/isicle/issues) - Known bugs and limitations
- [CREST Issue #203](https://github.com/crest-lab/crest/issues/203) - SIGSEGV in conformer generation

### Async Processing
- [FastAPI Background Tasks](https://fastapi.tiangolo.com/tutorial/background-tasks/) - Official documentation
- [Managing Long-Running Operations in FastAPI](https://leapcell.io/blog/managing-background-tasks-and-long-running-operations-in-fastapi) - Best practices
- [Idempotency, Retry, and Recovery](https://medium.com/@vivekburman1997/data-engineering-part-1-idempotency-retry-and-recovery-b3631a9b8b6f) - Design patterns
- [Python subprocess docs](https://docs.python.org/3/library/subprocess.html) - Timeout handling

### Input Validation
- [RDKit Sanitization Blog](https://greglandrum.github.io/rdkit-blog/posts/2025-06-27-sanitization-and-file-parsing.html) - SMILES parsing options
- [RDKit Issue #6996](https://github.com/rdkit/rdkit/issues/6996) - Error handling challenges

### Process Management
- [NREL ESIFHPC3 Issue #2](https://github.com/NREL/ESIFHPC3/issues/2) - NWChem job termination
- [Linux Zombie Process Cleanup](https://www.linuxjournal.com/content/how-kill-zombie-processes-linux) - Process management
