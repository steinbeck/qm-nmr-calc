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

## Google Cloud Platform Deployment

Deploy qm-nmr-calc to a cost-effective GCP Spot VM with lifecycle management scripts. This is ideal for teams that want cloud hosting without managed complexity.

### Why GCP Spot VMs?

Spot VMs offer **60-91% cost savings** compared to on-demand pricing. For NMR calculations that run intermittently, this makes cloud hosting very affordable:

- **e2-standard-4** (4 vCPU, 16 GB): ~$30-50/month Spot vs ~$100-150/month on-demand
- Pay only for compute time - stop the VM when not running calculations
- Persistent disk preserves all data between sessions

### One-Command Deployment

For first-time users, use the interactive quick deploy wizard:

```bash
cd gcp
./quick-deploy.sh
```

This guides you through the entire process: prerequisites, infrastructure, DNS, and deployment.

### Prerequisites

Before starting (the quick-deploy script will check these for you):

1. **GCP Account with Billing** - [Create account](https://console.cloud.google.com/)
2. **gcloud CLI** - [Install](https://cloud.google.com/sdk/docs/install) and authenticate:
   ```bash
   gcloud auth login
   gcloud config set project YOUR_PROJECT_ID
   ```
3. **Domain Name** - Required for HTTPS (e.g., nmr.example.com)
4. **Git** - To clone the repository

### Cost Estimates

| Resource | Spot Pricing | On-Demand Pricing | Notes |
|----------|--------------|-------------------|-------|
| e2-standard-4 (4 vCPU, 16 GB) | ~$30-50/month | ~$100-150/month | Primary compute cost |
| e2-standard-2 (2 vCPU, 8 GB) | ~$15-25/month | ~$50-75/month | Light workloads |
| n2-standard-4 (4 vCPU, 16 GB) | ~$50-80/month | ~$150-200/month | Higher performance |
| Static IP (attached) | Free | Free | While VM is running |
| Static IP (detached) | ~$7/month | ~$7/month | When VM is stopped |
| Persistent Disk (100 GB SSD) | ~$17/month | ~$17/month | Always billed |

**Monthly estimate for typical usage:** $35-70/month (Spot VM + disk + IP during downtime)

### Preemption and Job Loss

> **WARNING: GCP can reclaim Spot VMs with only 30 seconds notice.**

This has important implications for NMR calculations:

- **Jobs in progress will be lost** - There is no checkpoint/restart mechanism
- The shutdown script attempts graceful container stop within the 30-second window
- **Persistent disk survives** - All completed job data is preserved
- **Static IP is retained** - DNS configuration remains valid

**Best practice:** Use `./stop-vm.sh` when not running calculations to avoid surprise preemption. This also eliminates compute charges during idle periods.

### Quick Start

**Step 1: Clone and configure**

```bash
git clone https://github.com/steinbeck/qm-nmr-calc.git
cd qm-nmr-calc/gcp

# Create configuration
cp config.sh.example config.sh

# Edit config.sh with your GCP project ID
nano config.sh
# Set: GCP_PROJECT_ID="your-actual-project-id"
```

**Step 2: Check prerequisites**

```bash
./check-prerequisites.sh
```

This verifies:
- gcloud CLI is installed and authenticated
- Project exists and is accessible
- Billing is enabled
- Compute Engine API is enabled

**Step 3: Create infrastructure**

```bash
./setup-infrastructure.sh
```

This creates:
- Static external IP (for DNS)
- Firewall rules (HTTP, HTTPS, SSH)
- Persistent disk (100 GB SSD for job data)

Note the static IP address displayed at the end.

**Step 4: Configure DNS**

Point your domain to the static IP (see [DNS Configuration](#dns-configuration) below).

**Step 5: Verify DNS propagation**

```bash
./check-dns.sh your-domain.com
```

Wait until DNS resolves to the correct IP before proceeding.

**Step 6: Deploy VM**

```bash
./deploy-vm.sh
```

You'll be prompted for:
- Region and zone (defaults from config.sh)
- Machine type (e2-standard-4 recommended)
- Domain name (for HTTPS certificate)

The script automatically validates DNS before deploying.

**Step 7: Verify deployment**

```bash
./verify-deployment.sh
```

This checks:
- VM is running
- Containers are healthy
- HTTPS endpoint responds
- Persistent disk is mounted

Access your deployment at `https://your-domain.com`

### DNS Configuration

After `./setup-infrastructure.sh` displays your static IP, configure DNS:

#### Cloudflare

1. Go to your domain's DNS settings in Cloudflare dashboard
2. Add an A record:
   - **Type:** A
   - **Name:** `nmr` (or `@` for root domain)
   - **IPv4 address:** Your static IP
   - **Proxy status:** DNS only (gray cloud) - required for initial setup
   - **TTL:** Auto
3. After Let's Encrypt certificate is obtained, you can enable Cloudflare proxy (orange cloud)

#### Namecheap

1. Go to Domain List > Manage > Advanced DNS
2. Add a new record:
   - **Type:** A Record
   - **Host:** `nmr` (or `@` for root domain)
   - **Value:** Your static IP
   - **TTL:** Automatic

#### General Instructions (Any Provider)

Create an A record pointing to your static IP:

```
nmr.example.com.  A  203.0.113.10  TTL=300
```

**Note:** DNS propagation can take 5-15 minutes. Use [dnschecker.org](https://dnschecker.org) to verify propagation before deploying.

### Lifecycle Management

Use these scripts from the `gcp/` directory to manage your deployment:

| Task | Command | Description |
|------|---------|-------------|
| Check Prerequisites | `./check-prerequisites.sh` | Verify GCP setup before first deployment |
| Check DNS | `./check-dns.sh domain` | Verify DNS resolves to static IP |
| Verify Deployment | `./verify-deployment.sh` | Check all components are healthy |
| Stop VM | `./stop-vm.sh` | Stop billing for compute, preserve data |
| Start VM | `./start-vm.sh` | Resume VM, services auto-start |
| Delete VM | `./delete-vm.sh` | Remove VM instance, keep disk and IP |
| Check Status | `./status-vm.sh` | Show VM state, IP, running containers |
| SSH Access | `./ssh-vm.sh` | Open interactive shell on VM |
| View Logs | `./logs-vm.sh` | Stream container logs |

**Cost-saving workflow:**

```bash
# Done with calculations - stop to save money
./stop-vm.sh

# Ready to run more calculations - resume
./start-vm.sh

# Wait 2-3 minutes for containers to start
./status-vm.sh
```

### Troubleshooting GCP Deployment

#### Prerequisites Failing

Run `./check-prerequisites.sh` to diagnose issues:

**Compute Engine API not enabled:**
```
Compute Engine API is not enabled
```
Fix: `gcloud services enable compute.googleapis.com`

**Billing not enabled:**
```
Billing is not enabled for project
```
Fix: Enable billing at https://console.cloud.google.com/billing/linkedaccount

**Not authenticated:**
```
Not authenticated with gcloud
```
Fix: `gcloud auth login && gcloud config set project YOUR_PROJECT`

#### VM Won't Start

**Spot capacity unavailable:**
```
ZONE_RESOURCE_POOL_EXHAUSTED
```
- Try a different zone in the same region
- Try during off-peak hours
- Consider a different machine type

**Quota exceeded:**
```
Quota 'CPUS' exceeded
```
- Request quota increase in GCP Console > IAM & Admin > Quotas
- Or use a smaller machine type

#### HTTPS Not Working

**Check DNS propagation:**
```bash
dig +short nmr.example.com
# Should return your static IP
```

**Check Caddy logs:**
```bash
./ssh-vm.sh
docker compose -f /opt/qm-nmr-calc/docker-compose.yml logs caddy | grep -i error
```

**Common causes:**
- DNS not propagated (wait 5-15 minutes)
- Cloudflare proxy enabled during initial setup (use DNS-only mode)
- Domain doesn't match what was entered during `./deploy-vm.sh`

#### Containers Not Running

**Check startup script completion:**
```bash
./ssh-vm.sh
tail -50 /var/log/startup-script.log
```

**Check container status:**
```bash
./status-vm.sh
# Or manually:
./ssh-vm.sh
docker compose -f /opt/qm-nmr-calc/docker-compose.yml ps
```

**Restart containers:**
```bash
./ssh-vm.sh
cd /opt/qm-nmr-calc
docker compose down
docker compose up -d
```

#### VM Preempted Unexpectedly

This is normal for Spot VMs. Simply restart:
```bash
./start-vm.sh
```

The persistent disk retains all completed job data. Only jobs that were running at preemption time are lost.

## ARM64 / Apple Silicon

qm-nmr-calc fully supports ARM64 architecture, including Apple Silicon Macs (M1/M2/M3) and ARM-based cloud instances (AWS Graviton, Ampere).

### How It Works

Docker images are published as multi-architecture manifests. When you run `docker compose up -d`, Docker automatically pulls the correct image for your platform:

- **x86_64 hosts**: Pull `linux/amd64` images
- **ARM64 hosts**: Pull `linux/arm64` images

No configuration changes needed - the same `docker-compose.yml` works on both architectures.

### Known Considerations

| Topic | Details |
|-------|---------|
| Memory | NWChem requires minimum 2 GB memory per MPI process for reliable DFT calculations |
| Threading | OpenBLAS configured with `OMP_NUM_THREADS=1` to avoid thread contention with MPI |
| CPU detection | `NWCHEM_NPROC` auto-detects available CPUs (capped at 40) if not explicitly set |
| Numerical results | ARM64 results match x86_64 within acceptable tolerances (0.5 ppm 1H, 2.0 ppm 13C) |

### CI/CD Note

ARM64 images are built using GitHub Actions' native `ubuntu-24.04-arm` runners. This is **free for public repositories only**. If you fork this repository as private, ARM64 builds will not run.

## Related Documentation

- [README](../README.md) - Quick start and feature overview
- [Usage Guide](usage.md) - Web UI and API reference
- [Installation Guide](installation.md) - Development setup from source
- [Architecture](architecture.md) - Technical system design
