# Problems Encountered During v2.6 GCP Spot Deployment

## 1. Docker `nproc` Returns 1 Inside Containers

**The biggest issue.** Inside Docker containers, `nproc` returns 1 even when the container has access to all host CPUs (64 cores). This caused OpenMPI to think there was only 1 slot available, refusing to launch 64 processes with a "not enough slots" error.

**Fix:** Added `--oversubscribe` flag to the `mpirun` command (`f12abb6`).

## 2. CPU Count Detection

Related to the above -- the system needed a way to correctly determine how many cores to use. Python's `os.cpu_count()` works correctly inside Docker (unlike `nproc`), so the fix dynamically sets `NWCHEM_NPROC` based on the VM's actual CPU count (`4593443`).

## 3. Worker Memory Limit Too Low

The default Docker Compose memory limits were insufficient for a high-memory VM (256 GB). NWChem with 64 processes uses ~2 GB per process (~128 GB total). Had to increase the worker container memory limit to 120 GB, applied in two rounds (`e0b86eb`, `9bc8489`).

## 4. HTTPS/Caddy Configuration

Caddy was configured for HTTPS by default, but the deployment used a bare IP address (34.44.198.243) with no domain ownership. Caddy can't issue TLS certs for an IP.

**Fix:** Changed the Caddyfile to serve plain HTTP on port 80.

## 5. Conformer Progress Tracking Bug (Display Only)

The UI status bar doesn't update during conformer processing. It shows "Optimizing conformers (0/2)" even after conf_001 finishes. Conformer statuses stay "pending" even when complete. This is a display-only bug -- calculations proceed correctly, but users get no progress feedback.

**Status:** Still unfixed.

## 6. macOS vs Linux Memory Detection

The machine info footer feature needed to detect RAM, but `free` doesn't exist on macOS. Required platform-specific detection (`5ada354`).

## 7. Missing Production Dependency

`httpx` was used for the machine info feature but wasn't listed in production dependencies, causing import failures on deployment (`8a35cd7`).

## 8. Calculation Time Expectations

The pinene NMR shielding calculation took ~34 minutes per conformer on 64 cores (much longer than initially expected). A 2-conformer job takes ~80-90 minutes total. This wasn't a bug per se, but the job appeared stuck because the expected completion time was way off.

---

**Summary:** The hardest problems were all around the Docker-on-GCP environment -- container cgroup isolation breaking CPU detection, memory limits needing manual bumps, and MPI refusing to oversubscribe slots. The rest were typical deployment papercuts (HTTPS without a domain, missing deps, platform differences).
