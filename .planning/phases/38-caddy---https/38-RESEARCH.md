# Phase 38: Caddy + HTTPS - Research

**Researched:** 2026-02-03
**Domain:** Reverse proxy, automatic HTTPS, Docker deployment
**Confidence:** HIGH

## Summary

Caddy is a modern web server that provides automatic HTTPS via Let's Encrypt with zero configuration beyond specifying a domain name. For this phase, we need to add Caddy as a reverse proxy in front of the existing FastAPI container, handling TLS termination and automatic certificate management.

The key challenge is supporting two modes: (1) production mode where a domain is configured and HTTPS is automatic, and (2) development/localhost mode where no domain is set and HTTP-only access is provided. Caddy's environment variable substitution with default values (`{$DOMAIN:localhost}`) elegantly solves this.

The current docker-compose.yml exposes the API on port 8000. With Caddy, we'll change the API to only expose internally (no host port mapping) and have Caddy handle all external traffic on ports 80/443.

**Primary recommendation:** Use Caddyfile with `{$DOMAIN:localhost}` pattern to enable automatic HTTPS when DOMAIN is set, falling back to HTTP-only localhost mode when unset.

## Standard Stack

The established tools for this domain:

### Core
| Tool | Version | Purpose | Why Standard |
|------|---------|---------|--------------|
| Caddy | 2.10.x | Reverse proxy + automatic HTTPS | Only server with automatic HTTPS by default; handles Let's Encrypt without configuration |
| caddy:alpine | 2.10-alpine | Docker image | Official image, small footprint (~45MB) |

### Supporting
| Tool | Version | Purpose | When to Use |
|------|---------|---------|-------------|
| Docker named volumes | N/A | Persist certificates | Always - prevents rate limiting |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Caddy | Nginx + Certbot | More config, manual cert renewal, but more familiar to some |
| Caddy | Traefik | More features for microservices, but more complex for single-app |
| caddy:alpine | caddy:latest | Larger image but same functionality |

**Installation:**
No installation needed - uses official Docker image `caddy:alpine`.

## Architecture Patterns

### Recommended Project Structure
```
project/
├── docker-compose.yml       # Add caddy service
├── Caddyfile                 # Caddy configuration
├── .env.example              # Add DOMAIN variable
└── .env                      # User's configuration
```

### Pattern 1: Environment Variable Domain with Default
**What:** Use `{$DOMAIN:localhost}` to conditionally enable HTTPS
**When to use:** When supporting both production (with domain) and development (localhost)
**Example:**
```caddyfile
# Source: https://caddyexamples.com/examples/using-environment-variable-for-hostname/
{$DOMAIN:localhost} {
    reverse_proxy api:8000
}
```

When `DOMAIN` is unset/empty: serves `localhost` over HTTP (no certificate)
When `DOMAIN=example.com`: serves `example.com` over HTTPS (auto certificate from Let's Encrypt)

### Pattern 2: HTTP Prefix for Explicit HTTP-Only
**What:** Prefix address with `http://` to disable automatic HTTPS
**When to use:** When you explicitly want HTTP-only
**Example:**
```caddyfile
# Source: https://caddyserver.com/docs/caddyfile/concepts
http://localhost {
    reverse_proxy api:8000
}
```

### Pattern 3: Docker Compose Service Configuration
**What:** Standard Caddy service with proper volumes and networking
**When to use:** Always for Docker deployments
**Example:**
```yaml
# Source: https://github.com/caddyserver/caddy-docker
services:
  caddy:
    image: caddy:alpine
    restart: unless-stopped
    ports:
      - "80:80"
      - "443:443"
      - "443:443/udp"  # HTTP/3 (QUIC)
    volumes:
      - ./Caddyfile:/etc/caddy/Caddyfile:ro
      - caddy_data:/data
      - caddy_config:/config
    depends_on:
      - api

volumes:
  caddy_data:
    external: true  # Survives docker-compose down
  caddy_config:
```

### Anti-Patterns to Avoid
- **Anonymous volumes:** Using VOLUME without named mounts causes certificate loss and rate limiting
- **Destroying containers without volumes:** Certificates stored in `/data` must persist
- **Using $HOSTNAME in Docker:** This is reserved by Docker; use a different variable name like `$DOMAIN`
- **Exposing API port to host:** Once Caddy is added, API should only be accessible via internal Docker network

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Certificate management | Manual certbot scripts | Caddy automatic HTTPS | Edge cases with renewal, storage, multi-domain |
| HTTP to HTTPS redirect | Custom redirect rules | Caddy auto-redirect | Enabled by default, handles edge cases |
| Reverse proxy headers | Manual header config | Caddy defaults | X-Forwarded-For/Proto/Host set automatically |
| Health check endpoint | Custom endpoint | Existing /health | Already implemented in API |

**Key insight:** Caddy's value is handling TLS complexity automatically. Any manual configuration defeats the purpose.

## Common Pitfalls

### Pitfall 1: Certificate Rate Limiting
**What goes wrong:** Let's Encrypt returns 429 errors, deployment fails
**Why it happens:** Certificates not persisted, re-issued on every container restart
**How to avoid:** Use named volumes for `/data` and `/config`; mark data volume as `external: true`
**Warning signs:** "too many certificates already issued" errors in Caddy logs

### Pitfall 2: Port Conflicts
**What goes wrong:** Caddy fails to start, "address already in use"
**Why it happens:** Another service (or old API port mapping) using 80/443
**How to avoid:** Remove API's host port mapping; ensure no other services on 80/443
**Warning signs:** Immediate container crash on startup

### Pitfall 3: Network Isolation
**What goes wrong:** Caddy can't reach backend, "connection refused"
**Why it happens:** Services not on same Docker network
**How to avoid:** All services in same compose file share default network; use service name as hostname
**Warning signs:** 502 Bad Gateway errors

### Pitfall 4: Localhost HTTPS Confusion
**What goes wrong:** Browser warnings about untrusted certificates on localhost
**Why it happens:** Caddy uses self-signed certs for localhost, browser doesn't trust them
**How to avoid:** For localhost, use explicit `http://localhost` or accept self-signed cert
**Warning signs:** Browser security warnings in development

### Pitfall 5: Environment Variable Name Collision
**What goes wrong:** Domain configuration doesn't work
**Why it happens:** Using `$HOSTNAME` which Docker reserves
**How to avoid:** Use `$DOMAIN` or another variable name
**Warning signs:** Unexpected container hostname appearing in Caddyfile

## Code Examples

Verified patterns from official sources:

### Complete Caddyfile for qm-nmr-calc
```caddyfile
# Source: https://caddyserver.com/docs/caddyfile/concepts
# When DOMAIN is unset: serves http://localhost (no TLS)
# When DOMAIN is set: serves https://{DOMAIN} (auto TLS)

{$DOMAIN:localhost} {
    # Reverse proxy to FastAPI backend
    reverse_proxy api:8000

    # Enable compression
    encode gzip
}
```

### Complete docker-compose.yml Caddy Service
```yaml
# Source: https://github.com/caddyserver/caddy-docker
services:
  caddy:
    image: caddy:alpine
    restart: unless-stopped
    ports:
      - "80:80"
      - "443:443"
      - "443:443/udp"
    cap_add:
      - NET_ADMIN  # Required for QUIC/HTTP3
    volumes:
      - ./Caddyfile:/etc/caddy/Caddyfile:ro
      - caddy_data:/data
      - caddy_config:/config
    environment:
      - DOMAIN=${DOMAIN:-}
    depends_on:
      api:
        condition: service_healthy

volumes:
  caddy_data:
    name: qm-nmr-calc-caddy-data
  caddy_config:
```

### Modified API Service (remove host port)
```yaml
services:
  api:
    # ... existing config ...
    expose:
      - "8000"  # Internal only, not published to host
    # REMOVE: ports: - "${API_PORT:-8000}:8000"
```

### Updated .env.example
```bash
# Domain for HTTPS certificates (leave empty for localhost HTTP-only)
# When set, Caddy automatically obtains Let's Encrypt certificates
# Example: DOMAIN=nmr.example.com
DOMAIN=
```

### Health Check Integration
```caddyfile
{$DOMAIN:localhost} {
    reverse_proxy api:8000 {
        health_uri /health
        health_interval 30s
    }
    encode gzip
}
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Manual certbot + cron | Caddy automatic HTTPS | Caddy 2.0 (2020) | Zero-config certificates |
| Nginx + Let's Encrypt | Caddy single-binary | Caddy 2.0 | Simpler deployment |
| Individual certs per domain | Wildcard utilization | Caddy 2.8+ (2024) | Fewer cert requests |
| Separate HTTP/3 config | Built-in QUIC | Caddy 2.6+ | Just expose 443/udp |

**Deprecated/outdated:**
- Caddy 1.x: Completely different config format, not compatible
- Manual ACME clients: Caddy handles this automatically

## Open Questions

Things that couldn't be fully resolved:

1. **HTTP/3 (QUIC) Priority**
   - What we know: Requires 443/udp port and NET_ADMIN capability
   - What's unclear: Whether it's worth the complexity for this application
   - Recommendation: Include it (one line), but document it's optional

2. **Email for Let's Encrypt**
   - What we know: Recommended for certificate expiry notifications
   - What's unclear: Whether to make it configurable or hardcode
   - Recommendation: Add optional `ACME_EMAIL` env var

## Sources

### Primary (HIGH confidence)
- [Caddy Official Documentation - Automatic HTTPS](https://caddyserver.com/docs/automatic-https) - Certificate automation behavior
- [Caddy Official Documentation - Caddyfile Concepts](https://caddyserver.com/docs/caddyfile/concepts) - Environment variable syntax
- [Caddy Official Documentation - reverse_proxy](https://caddyserver.com/docs/caddyfile/directives/reverse_proxy) - Proxy configuration
- [Caddy Official Documentation - Global Options](https://caddyserver.com/docs/caddyfile/options) - auto_https control
- [Caddy Official Documentation - Patterns](https://caddyserver.com/docs/caddyfile/patterns) - Common configurations

### Secondary (MEDIUM confidence)
- [Caddy Docker GitHub](https://github.com/caddyserver/caddy-docker) - Docker image patterns
- [Caddy Examples - Environment Variable Hostname](https://caddyexamples.com/examples/using-environment-variable-for-hostname/) - Domain fallback pattern
- [endoflife.date/caddy](https://endoflife.date/caddy) - Version 2.10.2 current

### Tertiary (LOW confidence)
- Community discussions on rate limiting - Confirmed best practices around volumes

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - Official documentation and Docker image
- Architecture: HIGH - Well-documented patterns from official sources
- Pitfalls: HIGH - Multiple sources confirm rate limiting and volume issues

**Research date:** 2026-02-03
**Valid until:** 2026-03-03 (30 days - Caddy is stable)
