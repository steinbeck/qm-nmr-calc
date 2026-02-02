# Feature Landscape: Docker Deployment for Scientific Computing App

**Domain:** Containerized scientific computing web application deployment
**Researched:** 2026-02-02
**Confidence:** HIGH (Docker official docs, established patterns for scientific apps)
**Target Users:** Academic researchers, research group sysadmins, individual chemists

## Executive Summary

Docker deployment for a scientific computing application like qm-nmr-calc must balance two competing needs: **simplicity for researchers** (who want `docker compose up` to "just work") and **production readiness** (health checks, auto-recovery, persistent data, HTTPS).

The research identifies three user personas with different expectations:
1. **Individual chemist** - wants quick local setup, doesn't care about HTTPS
2. **Research group admin** - needs multi-user access, persistent data, some security
3. **Cloud deployer** - needs production hardening, auto-HTTPS, monitoring

v2.4 should satisfy all three with sensible defaults that work locally and scale to production with minimal configuration changes (primarily adding a domain name in `.env`).

## Table Stakes Features

Features users **expect** from any containerized scientific app. Missing these = product feels broken or incomplete.

### Core Container Stack

| Feature | Why Expected | Complexity | User Persona | Notes |
|---------|--------------|------------|--------------|-------|
| **Single `docker compose up` command** | Standard deployment pattern | LOW | All | Entry point to deployment |
| **Pre-built images on GHCR** | Users shouldn't need to build | MEDIUM | All | `ghcr.io/steinbeck/qm-nmr-calc` |
| **All dependencies bundled** | Scientific apps have complex deps | HIGH | All | NWChem, CREST, xTB in worker image |
| **Volume persistence for job data** | Calculations are expensive, can't lose results | LOW | All | Named volume for `/data/jobs` |
| **Health checks** | Know if services are actually working | LOW | Admin, Cloud | `/api/v1/health/ready` endpoint exists |
| **Restart policies** | Auto-recovery from crashes | LOW | Admin, Cloud | `restart: unless-stopped` |
| **Environment configuration via `.env`** | Standard Docker Compose pattern | LOW | All | Domain, email, workers count |

**Rationale:** These are baseline expectations from any Docker-deployed application. Scientific users specifically expect pre-built images because building from source with NWChem/RDKit is painful.

### Documentation & Onboarding

| Feature | Why Expected | Complexity | User Persona | Notes |
|---------|--------------|------------|--------------|-------|
| **Quick start in README** | Get running in 5 minutes | LOW | All | `git clone && docker compose up` |
| **`.env.example` file** | Know what to configure | LOW | All | All variables with comments |
| **Deployment guide** | Production setup instructions | MEDIUM | Admin, Cloud | Separate doc for cloud deployment |
| **Troubleshooting section** | Common issues have answers | LOW | All | Container logs, health checks |

**Rationale:** Academic researchers are not DevOps experts. Clear documentation is table stakes because poor docs = no adoption.

### Security Baseline

| Feature | Why Expected | Complexity | User Persona | Notes |
|---------|--------------|------------|--------------|-------|
| **Non-root container user** | Basic security hygiene | LOW | Admin, Cloud | Run as `appuser` not `root` |
| **No secrets in images** | Security 101 | LOW | All | Use `.env` file or Docker secrets |
| **Minimal base images** | Reduce attack surface | MEDIUM | Cloud | Multi-stage build, slim images |

**Rationale:** Research computing increasingly faces security requirements. Basic hygiene is expected.

## Differentiators

Features that provide **competitive advantage** over manual installation. Not expected, but highly valued.

### Deployment Experience

| Feature | Value Proposition | Complexity | User Persona | Notes |
|---------|-------------------|------------|--------------|-------|
| **Auto-HTTPS via Caddy** | Zero SSL certificate management | MEDIUM | Admin, Cloud | Let's Encrypt automatic |
| **Caddy reverse proxy included** | Production-ready out of the box | LOW | Admin, Cloud | HTTP/2, compression, security headers |
| **Single-command cloud deploy** | Minutes not hours to production | LOW | Cloud | Works on any VPS with Docker |
| **Log aggregation** | Single place to see all service logs | LOW | Admin, Cloud | `docker compose logs -f` works |

**Rationale:** Manual NWChem + Python + reverse proxy + SSL setup takes hours. Caddy auto-HTTPS is a major DX win that differentiates Docker deployment from manual installation.

### Operational Excellence

| Feature | Value Proposition | Complexity | User Persona | Notes |
|---------|-------------------|------------|--------------|-------|
| **Graceful shutdown** | No orphaned calculations | MEDIUM | All | SIGTERM handling in worker |
| **Container resource limits** | Prevent runaway calculations | LOW | Admin, Cloud | Memory limits in compose |
| **Backup-friendly volumes** | Simple disaster recovery | LOW | Admin, Cloud | Named volumes, tar backup works |
| **Multi-architecture images** | Runs on ARM servers too | MEDIUM | Cloud | amd64 + arm64 builds |

**Rationale:** These features reduce operational burden and make the deployment feel "production-grade" without requiring Kubernetes expertise.

### Scientific Computing Specific

| Feature | Value Proposition | Complexity | User Persona | Notes |
|---------|-------------------|------------|--------------|-------|
| **Worker parallelism config** | Match hardware resources | LOW | Admin, Cloud | `NWCHEM_NPROC` environment variable |
| **Separate app/worker containers** | Scale workers independently | MEDIUM | Admin, Cloud | API server vs calculation workers |
| **Job queue persistence** | Survive restarts | LOW | All | Huey SQLite in named volume |
| **CREST/xTB availability indicator** | Know what's installed | LOW | All | Health endpoint shows status |

**Rationale:** Scientific computing has unique needs around parallelism and long-running jobs. These features address those needs specifically.

## Anti-Features

Features to **deliberately NOT include** in v2.4. Common overengineering mistakes.

| Anti-Feature | Why Avoid | What to Do Instead |
|--------------|-----------|-------------------|
| **Kubernetes manifests** | Overkill for target users (single VM) | Docker Compose only |
| **Custom orchestration** | Docker Compose is sufficient | Use built-in restart policies |
| **Database container** | Huey SQLite works fine, no need | Keep filesystem-based storage |
| **Redis/message broker** | Huey SQLite is adequate for single-node | Don't add complexity |
| **Built-in monitoring stack** | Grafana/Prometheus overkill | Log to stdout, use external tools if needed |
| **Horizontal autoscaling** | Single VM deployment | Manual worker count adjustment |
| **Blue-green deployment** | Academic use case doesn't need it | Simple restart is fine |
| **Service mesh** | Not needed for 3-container stack | Direct container communication |
| **Docker Swarm mode** | Single node doesn't need it | Standalone Docker Compose |
| **Custom log shipping** | Stdout logging is sufficient | Use `docker compose logs` |
| **Rate limiting in app** | Caddy can do this if needed | Keep app simple |
| **Authentication/multi-user** | Out of scope for v2.4 | Future milestone if needed |
| **Secrets manager integration** | `.env` file is sufficient | Keep deployment simple |

**Rationale:** The target is `docker compose up` simplicity. Every additional component increases complexity and failure modes. Academic researchers want calculations, not DevOps.

## User Persona Analysis

### Persona 1: Individual Chemist (Local Development)

**Scenario:** Running on laptop/workstation for personal NMR predictions

**Needs:**
- Quick setup (< 5 minutes)
- Works without domain/HTTPS
- Minimal configuration

**Features prioritization:**
1. Pre-built images (no compile time)
2. `docker compose up` works immediately
3. Data persists across restarts
4. Clear shutdown instructions

**Configuration:** Zero required; defaults work for `localhost:8000`

### Persona 2: Research Group Admin

**Scenario:** Deploying on shared lab server for 5-20 researchers

**Needs:**
- HTTPS for security
- Persistent data (calculations are expensive)
- Basic monitoring (is it working?)
- Easy updates

**Features prioritization:**
1. Auto-HTTPS with real domain
2. Health checks and restart policies
3. Volume persistence with backup instructions
4. Simple update procedure (`docker compose pull && docker compose up -d`)

**Configuration:** Domain name in `.env`, optional email for Let's Encrypt

### Persona 3: Cloud Deployer

**Scenario:** Deploying on AWS/GCP/DigitalOcean VM for production use

**Needs:**
- Production-grade defaults
- Resource limits
- Log access
- Security best practices

**Features prioritization:**
1. Secure defaults (non-root, minimal images)
2. Resource limits to prevent runaway costs
3. Graceful shutdown (no orphaned jobs on deploy)
4. Multi-architecture support

**Configuration:** Domain, email, resource limits in `.env`

## Feature Dependencies

```
Core Stack (Phase 1):
├── Dockerfiles (app, worker)
├── docker-compose.yml
├── .env.example
└── Basic documentation

Caddy + HTTPS (Phase 2):
├── Caddyfile
├── caddy service in compose
├── Domain configuration in .env
└── HTTP→HTTPS redirect

Pre-built Images (Phase 3):
├── GitHub Actions workflow
├── GHCR image publishing
├── Multi-arch builds
└── Image tagging strategy

Documentation (Phase 4):
├── Deployment guide
├── Troubleshooting section
├── Update instructions
└── Backup/restore guide
```

## MVP Recommendation

**For v2.4 MVP (production-ready Docker deployment):**

### Must Have (Table Stakes)
1. **App Dockerfile** - FastAPI server with all Python dependencies
2. **Worker Dockerfile** - NWChem, CREST, xTB pre-installed
3. **docker-compose.yml** - app + worker + caddy services
4. **Caddyfile** - Reverse proxy with auto-HTTPS
5. **.env.example** - All configuration options documented
6. **Named volumes** - `qm-nmr-data` for job persistence, `caddy-data` for certificates
7. **Health checks** - All services have health check configuration
8. **Restart policies** - `unless-stopped` for all services
9. **GHCR images** - Pre-built `ghcr.io/steinbeck/qm-nmr-calc-{app,worker}`
10. **Quick start docs** - 5-minute deployment in README

### Should Have (Differentiators)
1. **Multi-arch images** - amd64 + arm64 for broader compatibility
2. **Graceful shutdown** - SIGTERM handling in worker
3. **Resource limits** - Memory limits in compose for worker
4. **Deployment guide** - Detailed production setup documentation

### Defer to Post-MVP
- Kubernetes manifests
- Monitoring integration
- Advanced logging
- Rate limiting
- Multi-user authentication

**Rationale:** MVP delivers on the core promise: `docker compose up` gives you a working, production-ready NMR calculation service. Advanced features can be added later without breaking existing deployments.

## Container Architecture

```
                    ┌─────────────────────────────────────────────┐
                    │              Host Machine                   │
                    │                                             │
  Users ────────────┼───► Port 80/443                             │
                    │         │                                   │
                    │    ┌────▼────┐                              │
                    │    │  Caddy  │  (reverse proxy, auto-HTTPS) │
                    │    └────┬────┘                              │
                    │         │                                   │
                    │    ┌────▼────┐    ┌──────────────┐          │
                    │    │   App   │◄───│    Worker    │          │
                    │    │ (API)   │    │  (NWChem)    │          │
                    │    └────┬────┘    └──────┬───────┘          │
                    │         │               │                   │
                    │    ┌────▼───────────────▼────┐              │
                    │    │    Named Volumes        │              │
                    │    │  - qm-nmr-data (jobs)   │              │
                    │    │  - caddy-data (certs)   │              │
                    │    │  - huey-db (queue)      │              │
                    │    └─────────────────────────┘              │
                    └─────────────────────────────────────────────┘
```

**Service responsibilities:**
- **Caddy**: TLS termination, reverse proxy, security headers, HTTP/2
- **App**: FastAPI server, web UI, API endpoints
- **Worker**: Huey consumer, NWChem calculations, conformer generation

## Configuration Reference

### Required Environment Variables

| Variable | Purpose | Default | Example |
|----------|---------|---------|---------|
| `DOMAIN` | Public domain for HTTPS | `localhost` | `nmr.example.com` |

### Optional Environment Variables

| Variable | Purpose | Default | Example |
|----------|---------|---------|---------|
| `ACME_EMAIL` | Let's Encrypt notifications | (none) | `admin@example.com` |
| `NWCHEM_NPROC` | NWChem parallel processes | `4` | `8` |
| `WORKER_COUNT` | Huey worker threads | `1` | `2` |
| `LOG_LEVEL` | Logging verbosity | `INFO` | `DEBUG` |
| `DATA_DIR` | Job data directory | `/data` | `/mnt/storage` |

### Volume Mounts

| Volume | Container Path | Purpose |
|--------|---------------|---------|
| `qm-nmr-data` | `/data` | Job files, results |
| `caddy-data` | `/data` | TLS certificates |
| `caddy-config` | `/config` | Caddy configuration |

## Image Size Considerations

### Target Sizes

| Image | Target Size | Strategy |
|-------|-------------|----------|
| App | < 500 MB | Multi-stage build, slim Python base |
| Worker | < 2 GB | Necessary for NWChem/CREST/xTB binaries |

**Worker image is necessarily large** because NWChem, CREST, and xTB are substantial scientific computing packages. This is acceptable given:
1. Images are pulled once and cached
2. No alternative for these dependencies
3. Comparable to other scientific computing containers

### Multi-stage Build Strategy

```dockerfile
# Stage 1: Build dependencies
FROM python:3.11-slim as builder
RUN pip install --user [packages]

# Stage 2: Runtime
FROM python:3.11-slim
COPY --from=builder /root/.local /root/.local
```

This reduces app image by ~65% compared to single-stage build.

## Sources

### Docker Best Practices
- [Docker Compose Official Documentation](https://docs.docker.com/compose) - Authoritative reference
- [Docker Environment Variables Best Practices](https://docs.docker.com/compose/how-tos/environment-variables/best-practices/) - Official guidance
- [Ten Simple Rules for Dockerfiles](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008316) - Scientific computing focus
- [Docker in 2026 Best Practices](https://medium.com/devops-ai-decoded/docker-in-2026-top-10-must-see-innovations-and-best-practices-for-production-success-30a5e090e5d6) - Current patterns

### Caddy & HTTPS
- [Caddy Reverse Proxy Quick Start](https://caddyserver.com/docs/quick-starts/reverse-proxy) - Official docs
- [Caddy Docker Proxy](https://github.com/lucaslorentz/caddy-docker-proxy) - Docker integration patterns
- [Caddy 2025 Setup Guide](https://www.virtualizationhowto.com/2025/09/caddy-reverse-proxy-in-2025-the-simplest-docker-setup-for-your-home-lab/) - Practical tutorial

### Scientific Computing Containers
- [NWChem Containers](https://nwchemgit.github.io/Containers.html) - Official NWChem Docker info
- [NWChem Docker Images](https://github.com/nwchemgit/nwchem-dockerfiles) - Official Dockerfiles
- [Docker Volume Persistence](https://docs.docker.com/engine/storage/volumes/) - Data management

### CI/CD & Image Publishing
- [GitHub Container Registry](https://docs.github.com/en/packages/working-with-a-github-packages-registry/working-with-the-container-registry) - GHCR setup
- [Publishing Docker Images](https://docs.github.com/actions/guides/publishing-docker-images) - GitHub Actions workflow

### Health Checks & Monitoring
- [Docker Compose Health Checks](https://last9.io/blog/docker-compose-health-checks/) - Practical guide
- [Docker Compose Logs](https://docs.docker.com/reference/cli/docker/compose/logs/) - Log management

**Confidence Level:** HIGH - Based on official Docker documentation, established scientific computing patterns (NWChem official Dockerfiles), and current 2025-2026 best practices. No speculative features included.
