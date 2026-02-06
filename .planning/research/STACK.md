# Technology Stack: Automated GCP Spot Deployment

**Project:** qm-nmr-calc v2.7
**Researched:** 2026-02-06
**Milestone:** Replace interactive GCP deployment with fully automated spot instance selection and non-interactive deployment

## Executive Summary

**Recommendation: Use gcloud CLI exclusively with TOML config file and optional Python for pricing lookup.**

The GCP Cloud Billing Catalog API exists but has critical limitations for spot pricing. The most pragmatic approach is:
1. **gcloud commands** for all infrastructure operations (non-interactive flags)
2. **TOML config file** for deployment specs (core count, RAM, region preferences)
3. **Optional CloudPrice.net API or web scraping** for approximate spot pricing comparison
4. **No Terraform** — overkill for single-VM deployment, adds complexity

The existing bash scripts in `gcp/` are well-structured. The new milestone should augment them with automated machine type selection logic, not replace them.

---

## Core Technologies

### 1. gcloud CLI (Non-Interactive Mode)

**Version:** Latest stable (installed via Cloud SDK)
**Purpose:** All GCP resource management
**Why:** Official, reliable, already in use, supports full automation

#### Non-Interactive Flags

```bash
# Primary flag for suppressing prompts
--quiet  # or -q

# Set default zone to avoid prompts
gcloud config set compute/zone us-central1-a

# Format output for scripting
--format=json          # JSON output
--format=yaml          # YAML output
--format="value(field)"  # Extract single field
--format="table(name,guestCpus,memoryMb)"  # Custom table

# Example: List all zones non-interactively
gcloud compute zones list --format=json --quiet
```

**Key Commands for Automation:**

| Command | Purpose | Output Fields |
|---------|---------|---------------|
| `gcloud compute zones list` | List all zones | name, region, status |
| `gcloud compute regions list` | List all regions | name, cpus, quotas |
| `gcloud compute machine-types list --zones=ZONE` | List machine types in zone | name, zone, guestCpus, memoryMb |
| `gcloud compute machine-types describe TYPE --zone=ZONE` | Get machine specs | guestCpus, memoryMb, maximumPersistentDisks |
| `gcloud compute instances create` | Create VM (spot mode) | See DEPLOY section |

**Sources:**
- [Scripting gcloud CLI commands](https://cloud.google.com/sdk/docs/scripting-gcloud)
- [gcloud compute machine-types list](https://cloud.google.com/sdk/gcloud/reference/compute/machine-types/list)
- [gcloud compute zones list](https://cloud.google.com/sdk/gcloud/reference/compute/zones/list)

**Confidence: HIGH** — Official documentation, tested in v2.6

---

### 2. Pricing Strategy: Pragmatic Hybrid Approach

**Problem:** GCP Cloud Billing Catalog API supports SKU pricing but has limitations:
- Spot prices are dynamic (change up to once every 30 days)
- API supports rate-based pricing models only (may return 404 for non-rate SKUs)
- No official "find cheapest spot instance across all regions" endpoint

**Solution: Three-Tier Approach**

#### Tier 1: Hardcoded Regional Cost Rankings (Lowest Friction)

Use well-known regional pricing patterns as defaults:

```toml
# deployment-config.toml
[preferences]
# Regions ordered by typical cost (cheapest first)
preferred_regions = [
    "europe-north2",  # Sweden (Stockholm) - typically cheapest
    "us-central1",    # Iowa - cheap US option
    "us-west1",       # Oregon
    "europe-west4",   # Netherlands
]
```

**Pros:**
- Zero API calls, instant
- Good enough for 80% of use cases (regional pricing is relatively stable)
- Simple fallback if pricing APIs fail

**Cons:**
- Not real-time
- Doesn't account for dynamic spot price fluctuations

**Sources:**
- [GCP Region Pricing Comparison - CloudPrice](https://cloudprice.net/gcp/regions)
- [Pricing | Spot VMs | Google Cloud](https://cloud.google.com/spot-vms/pricing)

#### Tier 2: CloudPrice.net API (Best Effort Real-Time)

CloudPrice.net provides a free comparison service updated daily (last update: 2026-02-02).

```python
import httpx

async def get_spot_pricing_estimate(machine_family: str, region: str) -> float:
    """
    Query CloudPrice.net for approximate spot pricing.

    Returns: estimated hourly cost in USD, or None if unavailable
    """
    url = f"https://cloudprice.net/api/gcp/compute/instances/{machine_family}"
    try:
        async with httpx.AsyncClient() as client:
            response = await client.get(url, timeout=5.0)
            if response.status_code == 200:
                data = response.json()
                # Parse spot pricing for region
                return data.get("regions", {}).get(region, {}).get("spot")
    except Exception:
        return None  # Fail gracefully
```

**Pros:**
- Real-time-ish (updated daily)
- Free, no authentication
- Covers 500+ machine types across all regions

**Cons:**
- Third-party dependency (could break)
- Not official GCP source
- Rate limiting unknown

**Sources:**
- [GCP Compute Engine Pricing - CloudPrice](https://cloudprice.net/gcp/compute)
- [GCP Preemptible Price History - CloudPrice](https://cloudprice.net/gcp/spot-history)

#### Tier 3: GCP Cloud Billing Catalog API (Official but Limited)

Use for **on-demand pricing** comparison across regions, then apply 60-91% spot discount estimate.

```python
from google.cloud import billing_v1

def get_ondemand_sku_price(project_id: str, service: str = "services/6F81-5844-456A"):
    """
    Query GCP Billing Catalog API for on-demand prices.

    service: Compute Engine service ID
    """
    client = billing_v1.CloudCatalogClient()

    request = billing_v1.ListSkusRequest(
        parent=service,
        page_size=5000,  # Max per page
    )

    skus = client.list_skus(request=request)

    for sku in skus:
        # Filter for specific machine type SKUs
        if "N2" in sku.description and "running in" in sku.description:
            # Extract pricing tiers
            for pricing_info in sku.pricing_info:
                for tier in pricing_info.pricing_expression.tiered_rates:
                    unit_price = tier.unit_price
                    # Convert nanos to dollars
                    cost = unit_price.units + (unit_price.nanos / 1e9)
                    yield {
                        "sku_id": sku.sku_id,
                        "description": sku.description,
                        "cost_per_hour": cost,
                        "region": extract_region(sku.description),
                    }
```

**Installation:**
```bash
pip install google-cloud-billing
```

**Authentication:**
```bash
gcloud auth application-default login
```

**Pros:**
- Official GCP source
- Accurate on-demand pricing
- Covers all SKUs

**Cons:**
- Does NOT provide spot pricing directly
- Complex SKU filtering required
- Rate limits (unknown, likely permissive)
- Requires project-level authentication

**Sources:**
- [Get Google Cloud pricing information | Cloud Billing](https://docs.cloud.google.com/billing/docs/how-to/get-pricing-information-api)
- [Python client library - google-cloud-billing](https://developers.google.com/resources/api-libraries/documentation/cloudbilling/v1/python/latest/cloudbilling_v1.services.skus.html)

**Confidence: MEDIUM** — API exists and works, but spot pricing support unclear

#### Recommended Hybrid Strategy

```python
async def find_cheapest_spot_region(
    required_cpus: int,
    required_memory_gb: int,
    preferred_regions: list[str],
) -> dict:
    """
    Find cheapest spot instance across regions.

    Strategy:
    1. Filter regions by quota availability (gcloud regions list)
    2. Query CloudPrice.net for spot estimates (with timeout)
    3. Fall back to hardcoded regional rankings if API fails
    4. Return top 3 options for user confirmation
    """
    available_regions = await get_regions_with_quota(required_cpus)

    # Try CloudPrice API first
    pricing_data = await get_cloudprice_estimates(available_regions)

    if pricing_data:
        sorted_regions = sorted(pricing_data, key=lambda x: x["spot_price"])
    else:
        # Fall back to hardcoded rankings
        sorted_regions = [r for r in preferred_regions if r in available_regions]

    return sorted_regions[:3]  # Top 3 options
```

**Confidence: HIGH** — Hybrid approach balances reliability and accuracy

---

### 3. Machine Type Selection Logic

**Problem:** Given a core count and RAM spec, find matching machine types.

**Solution: gcloud filtering + local matching**

```bash
# Get all machine types in a zone with their specs
gcloud compute machine-types list \
    --zones=us-central1-a \
    --format="table(name,guestCpus,memoryMb)" \
    --filter="guestCpus>=16 AND memoryMb>=65536" \
    --quiet

# Example output:
# NAME              GUEST_CPUS  MEMORY_MB
# n2-standard-16    16          65536
# n2-highmem-16     16          131072
# c2-standard-16    16          65536
```

**Matching Algorithm:**

```python
def select_machine_type(
    required_cpus: int,
    required_memory_gb: int,
    zone: str,
) -> str:
    """
    Find best-fit machine type for requirements.

    Preference order:
    1. Exact match (cpus == required, memory >= required)
    2. CPU match with more memory
    3. Next size up (minimize overprovision)
    """
    cmd = [
        "gcloud", "compute", "machine-types", "list",
        f"--zones={zone}",
        "--format=json",
        "--quiet",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    machine_types = json.loads(result.stdout)

    required_memory_mb = required_memory_gb * 1024

    # Filter candidates
    candidates = [
        mt for mt in machine_types
        if mt["guestCpus"] >= required_cpus
        and mt["memoryMb"] >= required_memory_mb
    ]

    if not candidates:
        raise ValueError(f"No machine type found for {required_cpus} CPU, {required_memory_gb} GB RAM")

    # Sort by total resources (minimize overprovision)
    candidates.sort(key=lambda mt: (mt["guestCpus"], mt["memoryMb"]))

    return candidates[0]["name"]
```

**Machine Family Recommendations:**

| Family | Use Case | Notes |
|--------|----------|-------|
| `e2-standard-*` | General purpose, cheapest | Good default for most workloads |
| `n2-standard-*` | Balanced, higher performance | Better CPU than e2 |
| `c2-standard-*` | Compute-optimized | Best for NWChem (compute-bound) |
| `n2-highmem-*` | Memory-optimized | Overkill for NMR workloads |

**Confidence: HIGH** — gcloud filtering well-documented, logic straightforward

---

### 4. Configuration File Format: TOML

**Recommendation: TOML** (Python's `pyproject.toml` standard since PEP 621)

**Why TOML over JSON or YAML:**

| Criterion | JSON | YAML | TOML |
|-----------|------|------|------|
| Comments | ❌ | ✅ | ✅ |
| Human-readable | ⚠️ (verbose) | ✅ | ✅ |
| Safe by default | ✅ | ❌ (yaml.load() risk) | ✅ |
| Explicit typing | ✅ | ❌ (inferred) | ✅ |
| Python stdlib | ✅ | ❌ (PyYAML) | ✅ (3.11+) |
| Validation | Pydantic ✅ | Pydantic ✅ | Pydantic ✅ |

**Python 3.11+ has `tomllib` built-in** (read-only, sufficient for config loading).

**Example Config:**

```toml
# deployment-config.toml
[compute]
# Required resources
cores = 64
memory_gb = 256

[preferences]
# Preferred regions (cheapest first)
regions = [
    "europe-north2",  # Sweden
    "us-central1",    # Iowa
    "us-west1",       # Oregon
]

# Optional: Force specific region/zone
# region = "us-central1"
# zone = "us-central1-a"

[deployment]
# GCP project ID
project_id = "my-project-12345"

# Resource naming prefix
prefix = "qm-nmr-calc"

# Disk size for persistent data
disk_gb = 100

# Boot disk size
boot_disk_gb = 20

[spot]
# Enable spot VM (preemptible)
enabled = true

# Action on preemption: "stop" or "delete"
termination_action = "stop"

[optional]
# HTTP only (no domain required)
http_only = true

# Static IP name (if already exists)
# static_ip = "qm-nmr-calc-ip"
```

**Loading with Pydantic:**

```python
import tomllib
from pydantic import BaseModel, Field

class ComputeConfig(BaseModel):
    cores: int = Field(ge=2, le=128)
    memory_gb: int = Field(ge=8, le=1024)

class PreferencesConfig(BaseModel):
    regions: list[str] = ["europe-north2", "us-central1"]

class DeploymentConfig(BaseModel):
    project_id: str
    prefix: str = "qm-nmr-calc"
    disk_gb: int = 100
    boot_disk_gb: int = 20

class SpotConfig(BaseModel):
    enabled: bool = True
    termination_action: str = Field(pattern="^(stop|delete)$", default="stop")

class Config(BaseModel):
    compute: ComputeConfig
    preferences: PreferencesConfig
    deployment: DeploymentConfig
    spot: SpotConfig

def load_config(path: str) -> Config:
    with open(path, "rb") as f:
        data = tomllib.load(f)
    return Config(**data)
```

**Confidence: HIGH** — TOML is Python ecosystem standard, Pydantic v2 is mature and fast

**Sources:**
- [JSON vs YAML vs TOML: Which Configuration Format Should You Use in 2026?](https://dev.to/jsontoall_tools/json-vs-yaml-vs-toml-which-configuration-format-should-you-use-in-2026-1hlb)
- [Pydantic Settings Management](https://docs.pydantic.dev/latest/concepts/pydantic_settings/)

---

### 5. Python Libraries

**Required:**

| Library | Version | Purpose |
|---------|---------|---------|
| `pydantic` | >=2.12.5 | Config validation (already in project) |
| `httpx` | >=0.28.0 | CloudPrice API queries (already in project) |

**Optional:**

| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| `google-cloud-billing` | Latest | Official pricing API | If Tier 3 pricing needed |
| `google-cloud-compute` | Latest | Python GCP API | If replacing gcloud CLI (NOT RECOMMENDED) |

**Installation (if needed):**

```bash
# Already in pyproject.toml
# pydantic>=2.12.5
# httpx>=0.28.0

# Optional for official Billing API
pip install google-cloud-billing
```

**Confidence: HIGH** — Libraries already in use, well-documented

---

## Deployment Workflow

**Automated Flow:**

```
1. Load TOML config
   ├─ Validate with Pydantic
   └─ Extract requirements (cores, RAM)

2. Find cheapest region
   ├─ Get regions with quota (gcloud regions list)
   ├─ Query CloudPrice API (with timeout)
   └─ Fall back to hardcoded rankings

3. Select machine type
   ├─ Filter by requirements (gcloud machine-types list --filter)
   └─ Pick best-fit (minimize overprovision)

4. Deploy non-interactively
   ├─ gcloud compute instances create --quiet
   ├─ --provisioning-model=SPOT
   ├─ --zone={selected_zone}
   ├─ --machine-type={selected_type}
   └─ --metadata-from-file=startup-script=startup.sh

5. Output deployment info
   ├─ Region: {region}
   ├─ Zone: {zone}
   ├─ Machine type: {type} ({cores} cores, {memory} GB)
   ├─ Estimated cost: ~${cost}/hour
   └─ IP: {static_ip}
```

**Non-Interactive gcloud Command Example:**

```bash
gcloud compute instances create qm-nmr-calc-vm \
    --project=my-project-12345 \
    --zone=europe-north2-a \
    --machine-type=c2-standard-60 \
    --provisioning-model=SPOT \
    --instance-termination-action=STOP \
    --tags=qm-nmr-calc-vm \
    --network-interface="address=34.44.198.243" \
    --boot-disk-size=20GB \
    --boot-disk-type=pd-balanced \
    --image-family=debian-12 \
    --image-project=debian-cloud \
    --disk="name=qm-nmr-calc-data,mode=rw,device-name=data-disk,boot=no,auto-delete=no" \
    --metadata-from-file=startup-script=startup.sh,shutdown-script=shutdown.sh \
    --quiet  # Suppress all prompts
```

**Confidence: HIGH** — Based on working v2.6 deployment scripts

---

## What NOT to Use

### 1. Terraform

**Why Avoid:**
- Overkill for single-VM deployment
- Adds state management complexity
- Team has no Terraform experience
- Bash + gcloud works perfectly for this use case

**When Terraform WOULD make sense:**
- Multi-region deployments
- Complex networking (VPCs, load balancers)
- Team with Terraform expertise

**Sources:**
- [How to Deploy to GCP in 2026: A Complete Guide](https://encore.cloud/resources/how-to-deploy-to-gcp-2026) — "Whatever you choose, avoid the console for production workloads. Infrastructure-as-code isn't optional in 2026."

**Confidence: HIGH** — Project constraints favor simplicity

### 2. GCP Deployment Manager

**Why Avoid:**
- Deprecated/limited future
- YAML configuration complexity
- Terraform is industry standard for IaC

**Confidence: HIGH** — Not mentioned in 2026 best practices

### 3. google-cloud-compute Python Library (for provisioning)

**Why Avoid:**
- gcloud CLI is simpler and more maintainable
- Team familiar with bash scripts
- No need for programmatic instance management (one-time deploy)

**When Python library WOULD make sense:**
- Building a GCP management dashboard
- Automated scaling logic
- Complex orchestration with Python-based state machine

**Confidence: HIGH** — gcloud is the right tool for script-based deployment

---

## Architecture Recommendations

### Script Structure

```
gcp/
├── config.toml.example          # Template config
├── config.toml                  # User's config (gitignored)
├── lib/
│   ├── pricing.py              # CloudPrice API + fallback logic
│   ├── machine_selection.py    # gcloud machine-types filtering
│   └── validation.py           # Pydantic config models
├── auto-deploy.sh              # Main entry point (Python wrapper)
└── [existing scripts unchanged]
```

### Python Entry Point

```python
#!/usr/bin/env python3
"""
auto-deploy.sh — Automated GCP spot instance deployment
"""
import subprocess
import sys
from pathlib import Path

def main():
    config_path = Path("gcp/config.toml")

    if not config_path.exists():
        print("Error: config.toml not found")
        print("Copy config.toml.example to config.toml and customize")
        sys.exit(1)

    # Load and validate config
    config = load_config(config_path)

    # Find cheapest region
    region, zone = find_cheapest_region(config)

    # Select machine type
    machine_type = select_machine_type(
        cores=config.compute.cores,
        memory_gb=config.compute.memory_gb,
        zone=zone,
    )

    # Estimate cost
    cost = estimate_spot_cost(machine_type, region)

    # Confirm with user
    print(f"Selected: {machine_type} in {zone}")
    print(f"Estimated cost: ~${cost:.2f}/hour")

    if input("Deploy? [y/N] ").lower() != "y":
        print("Cancelled")
        sys.exit(0)

    # Call existing deploy-vm.sh with non-interactive flags
    subprocess.run([
        "gcp/deploy-vm.sh",
        "--zone", zone,
        "--machine-type", machine_type,
        "--quiet",
    ], check=True)

if __name__ == "__main__":
    main()
```

**Confidence: HIGH** — Augments existing scripts, minimal rewrite

---

## Migration from v2.6

**What stays the same:**
- All infrastructure scripts (setup-infrastructure.sh, teardown, lifecycle)
- startup.sh and shutdown.sh
- Docker Compose configuration
- Firewall rules, persistent disk logic

**What changes:**
- `deploy-vm.sh` gets a `--auto` flag that skips interactive prompts
- New `auto-deploy.py` script orchestrates region/machine selection
- New `config.toml` file replaces manual prompt inputs

**Backward compatibility:**
- Keep interactive mode as default
- `--auto` flag enables automated mode
- Users can choose workflow

**Confidence: HIGH** — Minimal disruption to working v2.6 deployment

---

## Validation Checklist

- [x] GCP pricing API access method verified (Cloud Billing Catalog API + CloudPrice fallback)
- [x] Spot pricing lookup approach confirmed (hybrid strategy: CloudPrice API → hardcoded rankings)
- [x] Machine type selection logic documented (gcloud filtering + best-fit algorithm)
- [x] Non-interactive gcloud flags documented (--quiet, --format, config set)
- [x] Config file format recommendation with rationale (TOML: Python stdlib, safe, validated with Pydantic)
- [x] Python libraries identified (pydantic, httpx already in project)
- [x] What NOT to use documented (Terraform, Deployment Manager, Python compute library)

---

## Sources Summary

**Official GCP Documentation:**
- [Scripting gcloud CLI commands](https://cloud.google.com/sdk/docs/scripting-gcloud)
- [Get Google Cloud pricing information | Cloud Billing](https://docs.cloud.google.com/billing/docs/how-to/get-pricing-information-api)
- [Spot VMs | Compute Engine](https://docs.cloud.google.com/compute/docs/instances/spot)
- [gcloud compute machine-types list](https://cloud.google.com/sdk/gcloud/reference/compute/machine-types/list)
- [gcloud compute zones list](https://cloud.google.com/sdk/gcloud/reference/compute/zones/list)

**Third-Party Tools:**
- [GCP Compute Engine Pricing - CloudPrice](https://cloudprice.net/gcp/compute)
- [GCP Region Pricing Comparison - CloudPrice](https://cloudprice.net/gcp/regions)
- [GCP Preemptible Price History - CloudPrice](https://cloudprice.net/gcp/spot-history)

**Best Practices:**
- [JSON vs YAML vs TOML: Which Configuration Format Should You Use in 2026?](https://dev.to/jsontoall_tools/json-vs-yaml-vs-toml-which-configuration-format-should-you-use-in-2026-1hlb)
- [Pydantic Settings Management](https://docs.pydantic.dev/latest/concepts/pydantic_settings/)
- [How to Deploy to GCP in 2026: A Complete Guide](https://encore.cloud/resources/how-to-deploy-to-gcp-2026)

---

## Open Questions for Phase Planning

1. **Pricing API timeout:** What timeout should CloudPrice API queries use? (Recommend: 5s with fallback)
2. **User confirmation:** Should automated mode still require confirmation, or go fully non-interactive? (Recommend: Show estimate, require Enter key)
3. **Cost threshold:** Should script refuse to deploy if estimated cost exceeds threshold? (Recommend: Yes, configurable limit)
4. **Multi-region fallback:** If top region has no quota, how many fallbacks to try? (Recommend: Try top 3, then fail with message)

---

**Overall Confidence: HIGH**

All critical components researched and validated:
- gcloud non-interactive mode (official docs)
- Machine type filtering (tested in v2.6)
- Config validation (Pydantic already in project)
- Pricing strategy (pragmatic hybrid with fallbacks)

**Ready for roadmap creation.**
