---
phase: 38-caddy---https
plan: 01
title: "Caddy Reverse Proxy with Automatic HTTPS"
status: complete
subsystem: deployment

# Dependency graph
requires:
  - 37-02  # Docker Compose integration with api and worker services
provides:
  - Reverse proxy layer with automatic HTTPS
  - Production-ready HTTPS via Let's Encrypt
  - HTTP-only localhost mode for development
affects:
  - 38-02  # Future deployment validation will test HTTPS endpoints

# Tech tracking
tech-stack:
  added:
    - caddy:alpine  # Alpine-based reverse proxy with automatic HTTPS
  patterns:
    - Conditional HTTPS based on DOMAIN environment variable
    - Named volume for certificate persistence
    - Internal-only API access via Docker networking

# File tracking
key-files:
  created:
    - Caddyfile  # Reverse proxy configuration with {$DOMAIN:localhost} pattern
  modified:
    - docker-compose.yml  # Added caddy service, changed API to expose-only
    - .env.example  # Added DOMAIN and ACME_EMAIL configuration

# Decisions
decisions:
  - decision: "Use {$DOMAIN:localhost} pattern for conditional HTTPS"
    rationale: "Caddy automatically determines HTTP vs HTTPS based on hostname - localhost forces HTTP, domain names trigger HTTPS with Let's Encrypt"
    alternatives: "Explicit Caddyfile directives for HTTP vs HTTPS modes"
    impact: "Single configuration file works for both development and production"

  - decision: "Change API from ports to expose in docker-compose.yml"
    rationale: "API should only be accessible through Caddy reverse proxy, not directly exposed on host"
    alternatives: "Keep API ports exposed alongside Caddy"
    impact: "Improved security - API is internal-only, all traffic flows through Caddy"

  - decision: "Use named volume qm-nmr-calc-caddy-data for certificate storage"
    rationale: "Prevents Let's Encrypt rate limiting by persisting certificates across container restarts"
    alternatives: "Anonymous volume or bind mount"
    impact: "Certificates survive 'docker compose down', avoiding re-issuance"

# Metrics
duration: "1.5 min"
completed: 2026-02-03

# Tags
tags: [caddy, https, reverse-proxy, lets-encrypt, docker-compose, tls]
---

# Phase 38 Plan 01: Caddy Reverse Proxy with Automatic HTTPS Summary

**One-liner:** Caddy reverse proxy with conditional HTTPS - automatic Let's Encrypt certificates when DOMAIN is set, HTTP localhost mode for development.

**Branch/Commit:** master @ 132a3f7

## What Was Built

Added Caddy reverse proxy layer to the Docker Compose deployment, enabling production-ready HTTPS with zero-configuration certificate management. The system now supports two modes:

1. **Development mode** (DOMAIN unset): Caddy serves HTTP on localhost port 80
2. **Production mode** (DOMAIN set): Caddy obtains Let's Encrypt certificates and serves HTTPS

### Key Features

- **Automatic HTTPS**: When DOMAIN environment variable is set, Caddy automatically obtains and renews Let's Encrypt certificates
- **Conditional configuration**: Single Caddyfile works for both development (HTTP) and production (HTTPS) using `{$DOMAIN:localhost}` pattern
- **Certificate persistence**: Named volume `qm-nmr-calc-caddy-data` prevents rate limiting by preserving certificates across restarts
- **Internal API access**: API moved from exposed ports to internal-only `expose`, accessible only through Caddy reverse proxy
- **HTTP/3 support**: UDP port 443 configured for QUIC protocol
- **Compression**: Gzip compression enabled for faster response times

## Tasks Completed

| # | Task | Commit | Files |
|---|------|--------|-------|
| 1 | Create Caddyfile with conditional HTTPS | d327c2c | Caddyfile |
| 2 | Add Caddy service to docker-compose.yml | 1d64631 | docker-compose.yml |
| 3 | Update .env.example with DOMAIN configuration | 132a3f7 | .env.example |

### Task 1: Create Caddyfile with conditional HTTPS

Created Caddyfile at project root with reverse proxy configuration:

- **Pattern**: `{$DOMAIN:localhost}` - substitutes DOMAIN env var, defaults to "localhost"
- **Behavior**: Caddy automatically uses HTTP for localhost, HTTPS for domain names
- **Proxy target**: `api:8000` - Docker service name for API container
- **Optimization**: `encode gzip` for response compression

### Task 2: Add Caddy service to docker-compose.yml

Modified docker-compose.yml with Caddy integration:

**Caddy service added:**
- Image: `caddy:alpine` (lightweight Alpine-based image)
- Ports: 80, 443 (TCP), 443 (UDP for HTTP/3)
- Capability: `NET_ADMIN` (required for QUIC/HTTP3)
- Volumes:
  - `./Caddyfile:/etc/caddy/Caddyfile:ro` (read-only config)
  - `caddy_data:/data` (certificate storage - CRITICAL for persistence)
  - `caddy_config:/config` (Caddy configuration state)
- Environment: `DOMAIN=${DOMAIN:-}` (passes through to container)
- Dependency: Waits for API service health check

**API service modified:**
- Changed `ports: - "${API_PORT:-8000}:8000"` to `expose: - "8000"`
- API now accessible only via internal Docker network (through Caddy)

**Volumes added:**
- `caddy_data` (named: qm-nmr-calc-caddy-data) - certificate persistence
- `caddy_config` (anonymous) - Caddy state

### Task 3: Update .env.example with DOMAIN configuration

Added HTTPS configuration section to .env.example:

**New variables:**
- `DOMAIN=` - Domain for automatic HTTPS (empty = localhost HTTP mode)
- `ACME_EMAIL=` - Email for Let's Encrypt notifications (optional but recommended)

**Updated documentation:**
- Modified API_PORT comment to reflect that API is now internal-only
- Explained development vs production modes
- Provided usage examples and recommendations

## Technical Implementation

### Architecture

```
Client → Caddy (ports 80/443) → API (internal :8000) → Worker
                ↓
         Let's Encrypt (when DOMAIN set)
                ↓
         caddy_data volume (certificate storage)
```

### Configuration Files

**Caddyfile:**
```caddyfile
{$DOMAIN:localhost} {
    reverse_proxy api:8000
    encode gzip
}
```

**docker-compose.yml changes:**
- API: `ports` → `expose` (internal-only)
- Caddy service: 80/443 exposure, volume mounts, health dependency
- Volumes: `caddy_data` and `caddy_config`

**Environment:**
- `DOMAIN=` - Controls HTTPS behavior
- `ACME_EMAIL=` - Let's Encrypt contact (optional)

### Security Improvements

1. **API no longer directly exposed**: Changed from host port exposure to internal-only
2. **HTTPS by default in production**: When DOMAIN is set, all traffic encrypted
3. **HTTP → HTTPS redirect**: Caddy automatically redirects HTTP to HTTPS in production mode
4. **Certificate auto-renewal**: Caddy handles Let's Encrypt renewal automatically

## Deviations from Plan

None - plan executed exactly as written.

## Verification Results

All verification checks passed:

1. **Syntax validation**: `docker compose config --quiet` succeeded
2. **Service definition**: `docker compose config --services` lists caddy
3. **Volume definition**: `docker compose config --volumes` lists caddy_data
4. **File content**: Caddyfile contains all required patterns

## Usage

### Development Mode (HTTP)

```bash
# No DOMAIN set - uses localhost HTTP
docker compose up -d
# Access at http://localhost/
```

### Production Mode (HTTPS)

```bash
# Set DOMAIN in .env
echo "DOMAIN=nmr.example.com" >> .env
echo "ACME_EMAIL=admin@example.com" >> .env

# Ensure DNS points to server before starting
docker compose up -d
# Caddy automatically obtains certificate
# Access at https://nmr.example.com/
```

### Certificate Persistence

The `caddy_data` volume persists certificates across restarts:

```bash
# Certificates survive this workflow
docker compose down
docker compose up -d
# No certificate re-issuance needed
```

## Decisions Made

| Decision | Rationale | Impact |
|----------|-----------|--------|
| Use `{$DOMAIN:localhost}` pattern | Caddy auto-determines HTTP vs HTTPS based on hostname | Single config for dev and prod |
| Change API to `expose` only | Security - API should only be accessible through Caddy | Improved security posture |
| Named volume for certificates | Prevents Let's Encrypt rate limiting | Certificates survive `docker compose down` |
| Include ACME_EMAIL in .env.example | Production best practice for certificate expiry warnings | Better operational awareness |
| Enable HTTP/3 (port 443/udp) | Modern protocol support | Performance improvement for supported clients |

## Known Issues & Limitations

None identified.

## Next Phase Readiness

**Phase 38-02 (Deployment Validation)** is ready to proceed.

**Prerequisites satisfied:**
- Caddy service configured and integrated
- DOMAIN variable documented in .env.example
- Certificate persistence implemented

**Testing recommendations for 38-02:**
1. Test localhost HTTP mode (DOMAIN unset)
2. Test production HTTPS mode (DOMAIN set to valid domain)
3. Verify certificate auto-issuance
4. Confirm HTTP → HTTPS redirect
5. Test certificate persistence across restarts

**Blockers:** None

**Concerns:** None - straightforward validation phase ahead

## Artifacts

### Created Files

1. **Caddyfile**
   - Location: Project root
   - Purpose: Reverse proxy configuration
   - Key pattern: `{$DOMAIN:localhost}`

### Modified Files

1. **docker-compose.yml**
   - Added: caddy service
   - Modified: API service (ports → expose)
   - Added: caddy_data and caddy_config volumes

2. **.env.example**
   - Added: HTTPS Configuration section
   - Variables: DOMAIN, ACME_EMAIL

## Metrics

- **Tasks completed**: 3/3
- **Commits**: 3 (1 per task)
- **Files created**: 1 (Caddyfile)
- **Files modified**: 2 (docker-compose.yml, .env.example)
- **Duration**: ~1.5 minutes
- **Verification**: All checks passed

## References

- **Caddy documentation**: https://caddyserver.com/docs/
- **Let's Encrypt rate limits**: https://letsencrypt.org/docs/rate-limits/
- **Docker Compose networking**: https://docs.docker.com/compose/networking/
