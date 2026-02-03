# Docker Deployment Guide

Deploy qm-nmr-calc with Docker for production use. This guide covers cloud VPS setup, HTTPS configuration, and operational tasks.

## Quick Reference

| Task | Command |
|------|---------|
| Start services | `docker compose up -d` |
| Stop services | `docker compose down` |
| View logs | `docker compose logs -f` |
| Restart single service | `docker compose restart api` |
| Rebuild from source | `docker compose up -d --build` |
| Check service health | `docker compose ps` |

## Prerequisites

- Linux server (Ubuntu 22.04+ recommended)
- Docker Engine 24+ and Docker Compose v2
- Domain name (for HTTPS) or localhost access
- 4+ CPU cores, 8+ GB RAM recommended

### Install Docker

**Ubuntu/Debian:**

```bash
# Install Docker
curl -fsSL https://get.docker.com | sh

# Add current user to docker group (logout/login required)
sudo usermod -aG docker $USER

# Verify installation
docker --version
docker compose version
```

## Cloud VPS Deployment

### Step 1: Provision a Server

**Recommended specs:**
- 4 vCPUs, 8 GB RAM minimum (for production workloads)
- 2 vCPUs, 4 GB RAM minimum (for light usage/testing)
- 40 GB SSD storage
- Ubuntu 22.04 LTS

**Providers:**
- [DigitalOcean](https://www.digitalocean.com/) - Droplets starting at $24/mo (4 vCPU, 8 GB)
- [Linode](https://www.linode.com/) - Dedicated CPU starting at $30/mo
- [Hetzner](https://www.hetzner.com/) - Cost-effective EU hosting
- [Vultr](https://www.vultr.com/) - Global locations

### Step 2: Configure DNS

Point your domain to the server's IP address:

1. Get your server's public IP: `curl -4 ifconfig.me`
2. Add an A record in your DNS provider:
   - Type: A
   - Name: `nmr` (or `@` for root domain)
   - Value: Your server's IP address
   - TTL: 300 (5 minutes)

**Example:** `nmr.example.com` -> `203.0.113.10`

### Step 3: Deploy qm-nmr-calc

```bash
# Clone repository
git clone https://github.com/steinbeck/qm-nmr-calc.git
cd qm-nmr-calc

# Configure environment
cp .env.example .env

# Edit .env with your domain
nano .env
# Set: DOMAIN=nmr.example.com
# Set: ACME_EMAIL=admin@example.com (for Let's Encrypt notifications)

# Start services
docker compose up -d
```

### Step 4: Verify Deployment

```bash
# Check all services are running
docker compose ps

# Expected output:
# NAME                  SERVICE   STATUS    PORTS
# qm-nmr-calc-api-1     api       healthy   8000/tcp
# qm-nmr-calc-worker-1  worker    healthy
# qm-nmr-calc-caddy-1   caddy     running   0.0.0.0:80->80/tcp, 0.0.0.0:443->443/tcp

# Check Caddy obtained certificate
docker compose logs caddy | grep -i certificate
```

Visit `https://nmr.example.com` - you should see the qm-nmr-calc interface with a valid HTTPS certificate.

## Configuration Reference

### Environment Variables (.env)

| Variable | Default | Description |
|----------|---------|-------------|
| `DOMAIN` | (empty) | Domain for HTTPS. Empty = HTTP localhost mode |
| `ACME_EMAIL` | (empty) | Email for Let's Encrypt certificate notifications |
| `NWCHEM_NPROC` | 4 | MPI processes for NWChem (set to CPU count) |
| `OMP_NUM_THREADS` | 4 | OpenMP threads for CREST/xTB |
| `WORKER_MEMORY_LIMIT` | 8g | Container memory limit |
| `API_PORT` | 8000 | Internal API port (rarely needs changing) |

### Resource Recommendations

| Use Case | CPUs | RAM | NWCHEM_NPROC | OMP_NUM_THREADS |
|----------|------|-----|--------------|-----------------|
| Light testing | 2 | 4 GB | 2 | 2 |
| Standard use | 4 | 8 GB | 4 | 4 |
| Production | 8+ | 16 GB | 8 | 8 |

**Memory formula:** Allow ~2 GB per MPI process. For `NWCHEM_NPROC=4`, set `WORKER_MEMORY_LIMIT=8g`.

## HTTPS Configuration

### Automatic HTTPS (Recommended)

Caddy automatically obtains and renews Let's Encrypt certificates when `DOMAIN` is set:

```bash
# In .env
DOMAIN=nmr.example.com
ACME_EMAIL=admin@example.com
```

Requirements:
- Domain must resolve to your server's IP
- Ports 80 and 443 must be open (firewall)
- No other service using ports 80/443

### HTTP-Only Mode (Development)

Leave `DOMAIN` empty for localhost HTTP mode:

```bash
# In .env (or omit DOMAIN entirely)
DOMAIN=
```

Access via `http://localhost` (port 80).

### Firewall Configuration

**UFW (Ubuntu):**

```bash
sudo ufw allow 80/tcp
sudo ufw allow 443/tcp
sudo ufw allow 443/udp  # HTTP/3 (QUIC)
sudo ufw enable
```

**firewalld (Fedora/RHEL):**

```bash
sudo firewall-cmd --permanent --add-service=http
sudo firewall-cmd --permanent --add-service=https
sudo firewall-cmd --reload
```

## Troubleshooting

### Caddy Can't Obtain Certificate

**Symptoms:** HTTPS doesn't work, Caddy logs show ACME errors

**Check DNS propagation:**
```bash
# Should return your server's IP
dig +short nmr.example.com

# Or use web tool: https://dnschecker.org
```

**Check ports are open:**
```bash
# From another machine, or use online port checker
nc -zv your-server-ip 80
nc -zv your-server-ip 443
```

**Check Caddy logs:**
```bash
docker compose logs caddy | grep -i error
```

**Common causes:**
- DNS not propagated yet (wait 5-15 minutes)
- Firewall blocking ports 80/443
- Another service using ports (nginx, apache)

### Worker Container Crashes

**Symptoms:** Jobs fail, worker restarts repeatedly

**Check worker logs:**
```bash
docker compose logs worker --tail 100
```

**Common issues:**

**Out of memory:**
```
MemoryError
# Or: Killed
```
Fix: Increase `WORKER_MEMORY_LIMIT` in `.env` and restart

**MPI issues:**
```
OMPI_MCA_plm_rsh_agent: command not found
```
This is expected - the container uses OMPI_ALLOW_RUN_AS_ROOT internally.

**shm_size too small:**
```
MPI_Init failed
```
The default docker-compose.yml sets `shm_size: 512m` which should be sufficient.

### NWChem Calculation Fails

**Check job-specific logs:**
```bash
# Find job directory
docker compose exec api ls -la /app/data/jobs/

# View NWChem output for a specific job
docker compose exec api cat /app/data/jobs/{job_id}/nwchem.out
```

**Common NWChem errors:**

**Basis set not found:**
```
basis set "6-311+g(2d,p)" not found
```
This indicates a corrupted NWChem installation - rebuild the worker image.

**SCF convergence failure:**
```
convergence NOT achieved
```
The molecule geometry may be problematic. Try a different input structure.

### API Health Check Fails

**Check API logs:**
```bash
docker compose logs api --tail 50
```

**Test health endpoint directly:**
```bash
docker compose exec api curl -f http://localhost:8000/health
```

**Common causes:**
- Init container didn't set permissions (check `docker compose logs init`)
- Volume mount issues

### X11/RDKit Drawing Errors

**Symptoms:** Molecule images fail to generate

**In worker logs:**
```
cannot open display
libGL error
```

The worker container includes X11 libraries (libxrender1, libxext6, libexpat1) for headless rendering. If you see these errors:

```bash
# Rebuild worker image
docker compose build worker
docker compose up -d
```

### Container Won't Start

**Check for port conflicts:**
```bash
sudo lsof -i :80
sudo lsof -i :443
```

**Check Docker daemon:**
```bash
sudo systemctl status docker
```

**Check disk space:**
```bash
df -h
```

Worker images are large (~2.1 GB). Ensure sufficient disk space.

## Backup and Restore

### What to Backup

| Data | Location | Description |
|------|----------|-------------|
| Job data | `qm-nmr-calc-data` volume | Calculation inputs, outputs, NMR results |
| Caddy certificates | `qm-nmr-calc-caddy-data` volume | Let's Encrypt certificates |
| Configuration | `.env` file | Your deployment settings |

### Backup Job Data

```bash
# Create backup directory
mkdir -p ~/backups

# Backup job data volume
docker run --rm \
  -v qm-nmr-calc-data:/data:ro \
  -v ~/backups:/backup \
  alpine tar czf /backup/qm-nmr-calc-data-$(date +%Y%m%d).tar.gz -C /data .

# Backup .env file
cp .env ~/backups/.env.backup
```

### Restore Job Data

```bash
# Stop services
docker compose down

# Remove existing volume (CAUTION: destroys current data)
docker volume rm qm-nmr-calc-data

# Restore from backup
docker run --rm \
  -v qm-nmr-calc-data:/data \
  -v ~/backups:/backup \
  alpine tar xzf /backup/qm-nmr-calc-data-20260203.tar.gz -C /data

# Fix permissions
docker run --rm \
  -v qm-nmr-calc-data:/data \
  alpine chown -R 999:999 /data

# Restart services
docker compose up -d
```

### Migrate to New Server

1. **On old server:** Create backups (see above)
2. **Transfer files:**
   ```bash
   scp ~/backups/qm-nmr-calc-data-*.tar.gz newserver:~/backups/
   scp .env newserver:~/qm-nmr-calc/
   ```
3. **On new server:** Restore data (see above)
4. **Update DNS** to point to new server IP

### Automated Backups (Optional)

Create a cron job for daily backups:

```bash
# Edit crontab
crontab -e

# Add daily backup at 2 AM
0 2 * * * cd /path/to/qm-nmr-calc && docker run --rm -v qm-nmr-calc-data:/data:ro -v ~/backups:/backup alpine tar czf /backup/qm-nmr-calc-data-$(date +\%Y\%m\%d).tar.gz -C /data . 2>/dev/null

# Keep only last 7 days of backups
0 3 * * * find ~/backups -name "qm-nmr-calc-data-*.tar.gz" -mtime +7 -delete
```

## Updating

### Update to New Version

```bash
cd qm-nmr-calc

# Pull latest changes
git pull

# Pull new images
docker compose pull

# Restart with new images
docker compose up -d

# Verify
docker compose ps
```

### Build from Source (Development)

```bash
# Rebuild images locally
docker compose up -d --build
```

## Architecture

### Service Overview

```
                    [Internet]
                        |
                   [Port 80/443]
                        |
                    [Caddy] -----> Auto-HTTPS via Let's Encrypt
                        |
                   [Port 8000]
                        |
                    [API] --------> FastAPI web server
                        |
                    [Volume]
                        |
                    [Worker] -----> NWChem/CREST/xTB calculations
```

### Data Flow

1. User submits molecule via web UI or API
2. API creates job record in `qm-nmr-calc-data` volume
3. Worker picks up job from Huey queue
4. Worker runs NWChem DFT and NMR calculations
5. Results stored in job directory
6. User retrieves results via API

### Volumes

| Volume | Purpose |
|--------|---------|
| `qm-nmr-calc-data` | Job data, calculation results |
| `qm-nmr-calc-caddy-data` | Let's Encrypt certificates |
| `qm-nmr-calc-caddy-config` | Caddy runtime config |

## Related Documentation

- [README](../README.md) - Quick start and feature overview
- [Usage Guide](usage.md) - Web UI and API reference
- [Installation Guide](installation.md) - Development setup from source
- [Architecture](architecture.md) - Technical system design
