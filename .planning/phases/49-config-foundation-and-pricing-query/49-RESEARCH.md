# Phase 49: Config Foundation and Pricing Query - Research

**Researched:** 2026-02-06
**Domain:** Configuration management, API-based pricing discovery, GCP automation
**Confidence:** MEDIUM

## Summary

This phase establishes non-interactive GCP deployment by introducing TOML-based configuration with Pydantic validation and automated spot pricing discovery via CloudPrice.net API. The research focused on four core domains: (1) TOML parsing with Python 3.11+ tomllib, (2) Pydantic v2 validation patterns for config files, (3) CloudPrice.net API structure and reliability, and (4) non-interactive gcloud scripting patterns.

The standard approach is to use Python for config validation (tomllib + Pydantic) before deployment, then pass validated values to existing bash deployment scripts. CloudPrice.net API provides pricing data but reliability is MEDIUM confidence - the API documentation portal exists but specific endpoint details require API key access. Hardcoded regional fallback rankings are essential as a backup strategy.

Key architectural decisions: Python validation script outputs shell-compatible exports for bash consumption, httpx with simple file-based TTL caching (no external caching libraries needed), and gcloud commands use --quiet --format=json for programmatic parsing with jq.

**Primary recommendation:** Create Python validation script that exits with clear error messages, outputs validated config as bash-compatible environment variables, and separates pricing discovery logic from config validation for independent testing.

## Standard Stack

The established libraries/tools for this domain:

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| tomllib | Built-in 3.11+ | Parse TOML config files | Python standard library, zero dependencies, TOML 1.0.0 compliant |
| Pydantic | 2.12.5+ (in deps) | Config validation with clear errors | Industry standard for Python data validation, already in project |
| httpx | 0.28.0+ (in deps) | HTTP client for pricing API | Already in project, async-capable, modern requests alternative |
| gcloud CLI | Latest | GCP operations with JSON output | Official Google Cloud SDK, required for all GCP operations |
| jq | 1.7+ | Parse JSON from gcloud commands | De facto standard for JSON parsing in bash |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| yq (kislyuk) | Latest | TOML-to-JSON conversion in bash | If bash scripts need to read TOML directly (use tomlq command) |
| json.tool | Built-in | Pretty-print JSON for debugging | Python fallback if jq unavailable |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| CloudPrice.net API | GCP Cloud Billing API | Billing API requires OAuth, more complex auth, not designed for spot pricing queries across all regions |
| tomllib (Python) | toml/tomli PyPI packages | tomllib is standard library in 3.11+, no reason to add dependency |
| File-based cache | hishel/requests-cache | External dependencies for simple 24-hour TTL cache, overkill for this use case |
| TOML config | YAML config | YAML has security issues (code execution), more complex spec, less human-readable |

**Installation:**
```bash
# Core tools already available
# Python 3.11+ has tomllib built-in
# Pydantic and httpx already in pyproject.toml
# jq typically pre-installed on macOS/Linux, or:
brew install jq  # macOS
apt-get install jq  # Debian/Ubuntu
```

## Architecture Patterns

### Recommended Project Structure
```
gcp/
├── config.toml.example     # User-facing config template
├── config.toml             # User's actual config (gitignored)
├── lib/
│   ├── config.sh           # Bash functions to load validated config
│   ├── pricing.sh          # Bash functions for pricing queries
│   ├── machine.sh          # Bash functions for machine type selection
│   └── infra.sh            # Bash functions for infrastructure ops
├── validate-config.py      # Python: validate TOML, output bash exports
├── query-pricing.py        # Python: CloudPrice.net API with caching
└── deploy-auto.sh          # Main orchestrator (Phase 51)
```

### Pattern 1: Python Validation with Bash Export
**What:** Python validates config, outputs bash-compatible environment variables
**When to use:** Bridging Python validation with bash deployment scripts
**Example:**
```python
# validate-config.py
import sys
import tomllib
from pathlib import Path
from pydantic import BaseModel, Field, field_validator, ValidationError

class GCPConfig(BaseModel):
    project_id: str = Field(..., min_length=6, max_length=30)
    resource_prefix: str = Field(default="qm-nmr-calc")
    cpu_cores: int = Field(..., ge=4, le=224)  # GCP max is 224 for N2D
    ram_gb: int = Field(..., ge=8, le=896)  # GCP max is 896 for N2D
    disk_size_gb: int = Field(default=100, ge=10, le=65536)

    @field_validator('project_id')
    @classmethod
    def validate_project_id(cls, v: str) -> str:
        # GCP project IDs: lowercase letters, numbers, hyphens
        if not all(c.islower() or c.isdigit() or c == '-' for c in v):
            raise ValueError('must contain only lowercase letters, numbers, hyphens')
        if v.startswith('-') or v.endswith('-'):
            raise ValueError('cannot start or end with hyphen')
        return v

    @field_validator('ram_gb')
    @classmethod
    def validate_ram_ratio(cls, v: int, info) -> int:
        # Common GCP ratios: 0.5-1 GB/core (highcpu), 2-4 GB/core (standard), 6-8 GB/core (highmem)
        # Allow range 0.5 to 8 GB per core
        if 'cpu_cores' in info.data:
            ratio = v / info.data['cpu_cores']
            if ratio < 0.5 or ratio > 8:
                raise ValueError(f'RAM/CPU ratio {ratio:.1f} outside valid range (0.5-8 GB/core)')
        return v

def main():
    config_path = Path('gcp/config.toml')
    if not config_path.exists():
        print("ERROR: gcp/config.toml not found", file=sys.stderr)
        print("Create from gcp/config.toml.example", file=sys.stderr)
        sys.exit(1)

    try:
        with open(config_path, 'rb') as f:
            data = tomllib.load(f)

        config = GCPConfig.model_validate(data.get('gcp', {}))

        # Output bash-compatible exports
        print(f"export GCP_PROJECT_ID='{config.project_id}'")
        print(f"export RESOURCE_PREFIX='{config.resource_prefix}'")
        print(f"export CPU_CORES={config.cpu_cores}")
        print(f"export RAM_GB={config.ram_gb}")
        print(f"export DISK_SIZE_GB={config.disk_size_gb}")

    except tomllib.TOMLDecodeError as e:
        print(f"ERROR: Invalid TOML syntax at line {e.lineno}: {e.msg}", file=sys.stderr)
        sys.exit(1)
    except ValidationError as e:
        print("ERROR: Config validation failed:", file=sys.stderr)
        for error in e.errors():
            field = '.'.join(str(x) for x in error['loc'])
            print(f"  {field}: {error['msg']}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()
```

### Pattern 2: Simple TTL File Cache
**What:** Write API response to file with timestamp, check age before reuse
**When to use:** 24-hour caching for pricing data
**Example:**
```python
# query-pricing.py - simplified caching
import json
import time
from pathlib import Path
from datetime import datetime

CACHE_DIR = Path('.cache')
CACHE_TTL = 24 * 60 * 60  # 24 hours

def get_cached_pricing(cache_file: Path) -> dict | None:
    """Return cached data if fresh, None if stale/missing"""
    if not cache_file.exists():
        return None

    with open(cache_file) as f:
        data = json.load(f)

    cached_at = data.get('cached_at', 0)
    age = time.time() - cached_at

    if age > CACHE_TTL:
        return None  # Stale

    return data.get('pricing')

def save_pricing_cache(cache_file: Path, pricing: dict):
    """Save pricing data with timestamp"""
    cache_file.parent.mkdir(exist_ok=True)
    with open(cache_file, 'w') as f:
        json.dump({
            'cached_at': time.time(),
            'pricing': pricing
        }, f)
```

### Pattern 3: Non-Interactive gcloud with JSON
**What:** All gcloud commands use --quiet --format=json, parsed with jq
**When to use:** Any gcloud operation in automated scripts
**Example:**
```bash
# Get machine types matching CPU/RAM requirements
get_matching_machine_types() {
    local zone=$1
    local min_cpus=$2
    local min_memory_mb=$3

    gcloud compute machine-types list \
        --zones="$zone" \
        --format=json \
        --quiet \
    | jq -r --arg cpus "$min_cpus" --arg mem "$min_memory_mb" '
        .[] |
        select(.guestCpus >= ($cpus | tonumber)) |
        select(.memoryMb >= ($mem | tonumber)) |
        "\(.name)\t\(.guestCpus)\t\(.memoryMb)"
    '
}
```

### Anti-Patterns to Avoid
- **Parsing gcloud stderr messages:** Messages change between versions, use structured JSON output instead
- **Running gcloud in parallel:** Not supported, will cause auth/lock conflicts
- **Hardcoded pricing values:** Spot pricing changes frequently, must query or have dated fallback
- **Validating CPU/RAM after GCP operations start:** Catch invalid combinations before any gcloud commands
- **Using YAML for config:** Security risk (code execution via tags), use TOML instead
- **Installing heavy caching libraries:** Simple file cache is sufficient for 24-hour TTL

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| TOML parsing | Regex-based parser | tomllib (stdlib) | TOML 1.0.0 has complex edge cases, datetime handling, table nesting |
| Config validation | Manual if/else checks | Pydantic models | Type coercion, nested validation, clear error messages, reusable |
| HTTP caching | Dict with timestamps | File-based cache or hishel | Race conditions, cache invalidation, disk persistence |
| JSON parsing in bash | sed/awk/grep chains | jq | Handles nested objects, arrays, escaping, edge cases |
| GCP machine type matching | String filtering | gcloud --format=json + jq | Machine type naming is inconsistent, filtering by specs more reliable |

**Key insight:** Config validation and API integration have subtle edge cases that make custom implementations error-prone. Python's stdlib (tomllib) and Pydantic handle these correctly. For bash integration, jq is the standard tool - don't try to parse JSON with bash builtins.

## Common Pitfalls

### Pitfall 1: CloudPrice.net API Unreliability
**What goes wrong:** API might be rate-limited, down, or require authentication not documented in public search results
**Why it happens:** Third-party API without public SLA, documentation requires API key access
**How to avoid:** Implement hardcoded regional fallback rankings before testing API integration. Test API failure path explicitly.
**Warning signs:** HTTP 429 (rate limit), HTTP 401 (auth required), HTTP 503 (service down), connection timeouts

### Pitfall 2: Opening TOML Files in Text Mode
**What goes wrong:** tomllib.load() requires binary mode ('rb'), text mode causes TypeError
**Why it happens:** TOML parser needs bytes for encoding detection and consistent parsing
**How to avoid:** Always open TOML files with 'rb' mode: `open('config.toml', 'rb')`
**Warning signs:** TypeError: "a bytes-like object is required, not 'str'"

### Pitfall 3: GCP CPU/RAM Ratio Validation
**What goes wrong:** User requests 64 cores + 32 GB RAM (0.5 GB/core), gcloud accepts but no machines exist
**Why it happens:** GCP machine families have predefined ratios, not all CPU/RAM combinations available
**How to avoid:** Validate ratio is within 0.5-8 GB/core range, use gcloud machine-types list to verify availability
**Warning signs:** Config validation passes but no matching machine types found in any zone

### Pitfall 4: Bash Variable Injection in Python Exports
**What goes wrong:** Config value contains bash special characters, breaks sourced exports
**Why it happens:** User project ID or prefix contains quotes, spaces, or shell metacharacters
**How to avoid:** Use single quotes around all exported values, escape single quotes in values
**Warning signs:** Bash syntax errors when sourcing validated config, unexpected command execution

### Pitfall 5: Stale Pricing Cache After Region Outage
**What goes wrong:** Cached pricing data includes regions now unavailable, deployment fails
**Why it happens:** 24-hour cache doesn't account for regional capacity changes
**How to avoid:** Validate machine type availability in selected zone before VM creation, implement fallback to next-cheapest region
**Warning signs:** gcloud VM creation fails with "zone does not have enough resources" after pricing query succeeds

### Pitfall 6: gcloud Parallel Execution
**What goes wrong:** Multiple gcloud commands run simultaneously, cause lock contention or auth errors
**Why it happens:** Scripts attempt to parallelize infrastructure checks
**How to avoid:** Never run gcloud commands in parallel, official docs explicitly warn against this
**Warning signs:** "Failed to acquire lock" errors, transient auth failures

### Pitfall 7: Missing TOML [sections]
**What goes wrong:** Config validation fails because tomllib returns empty dict
**Why it happens:** TOML requires explicit [section] headers, data.get('gcp', {}) returns empty if [gcp] missing
**How to avoid:** Validate TOML structure before Pydantic validation, check section exists
**Warning signs:** All fields reported as "missing" even when values present in file

## Code Examples

Verified patterns from official sources:

### Loading and Validating TOML Config
```python
# Source: https://docs.python.org/3/library/tomllib.html
import tomllib
from pathlib import Path
from pydantic import BaseModel, ValidationError

class Config(BaseModel):
    project_id: str
    cpu_cores: int
    ram_gb: int

def load_config(path: Path) -> Config:
    """Load and validate TOML config, exit with clear errors"""
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")

    # Must use binary mode
    with open(path, 'rb') as f:
        try:
            data = tomllib.load(f)
        except tomllib.TOMLDecodeError as e:
            raise ValueError(
                f"Invalid TOML at line {e.lineno}, column {e.colno}: {e.msg}"
            )

    # Validate with Pydantic
    try:
        return Config.model_validate(data['gcp'])
    except ValidationError as e:
        # Format errors for end users
        errors = []
        for err in e.errors():
            field = '.'.join(str(x) for x in err['loc'])
            errors.append(f"{field}: {err['msg']}")
        raise ValueError("Config validation failed:\n  " + "\n  ".join(errors))
```

### Non-Interactive gcloud with Error Handling
```bash
# Source: https://docs.cloud.google.com/sdk/docs/scripting-gcloud
set -euo pipefail  # Exit on error, undefined vars, pipe failures

# Get static IP, exit with clear error if not found
get_static_ip() {
    local ip_name=$1
    local region=$2

    local ip_address
    ip_address=$(gcloud compute addresses describe "$ip_name" \
        --region="$region" \
        --format="value(address)" \
        --quiet 2>&1) || {
        echo "ERROR: Static IP '$ip_name' not found in region $region" >&2
        echo "Run: ./setup-infrastructure.sh" >&2
        return 1
    }

    echo "$ip_address"
}
```

### Querying CloudPrice.net with Caching
```python
# Source: CloudPrice.net API patterns (LOW confidence - API docs require key)
import httpx
import json
from pathlib import Path
from datetime import datetime, timedelta

CACHE_FILE = Path('.cache/gcp-spot-pricing.json')
CACHE_TTL = timedelta(hours=24)

def get_spot_pricing() -> dict:
    """Query CloudPrice.net API with 24-hour cache"""

    # Check cache first
    if CACHE_FILE.exists():
        with open(CACHE_FILE) as f:
            cached = json.load(f)

        cached_at = datetime.fromisoformat(cached['timestamp'])
        if datetime.now() - cached_at < CACHE_TTL:
            return cached['pricing']

    # Cache miss, query API
    # NOTE: Actual endpoint format TBD - requires CloudPrice.net API key
    response = httpx.get(
        'https://api.cloudprice.net/v2/gcp/compute/spot-pricing',
        headers={'X-API-Key': get_api_key()},
        timeout=10.0
    )
    response.raise_for_status()

    pricing = response.json()

    # Save to cache
    CACHE_FILE.parent.mkdir(exist_ok=True)
    with open(CACHE_FILE, 'w') as f:
        json.dump({
            'timestamp': datetime.now().isoformat(),
            'pricing': pricing
        }, f)

    return pricing
```

### Pydantic Field Validators for Config
```python
# Source: https://docs.pydantic.dev/latest/concepts/validators/
from pydantic import BaseModel, field_validator, ValidationInfo

class GCPConfig(BaseModel):
    project_id: str
    cpu_cores: int
    ram_gb: int

    @field_validator('project_id')
    @classmethod
    def validate_project_id(cls, v: str) -> str:
        """GCP project IDs: 6-30 chars, lowercase letters, digits, hyphens"""
        if len(v) < 6 or len(v) > 30:
            raise ValueError('must be 6-30 characters')
        if not all(c.islower() or c.isdigit() or c == '-' for c in v):
            raise ValueError('must contain only lowercase letters, numbers, hyphens')
        if v[0] == '-' or v[-1] == '-':
            raise ValueError('cannot start or end with hyphen')
        return v

    @field_validator('ram_gb')
    @classmethod
    def validate_ram_cpu_ratio(cls, v: int, info: ValidationInfo) -> int:
        """Check RAM/CPU ratio matches GCP machine families"""
        if 'cpu_cores' in info.data:
            ratio = v / info.data['cpu_cores']
            if ratio < 0.5:
                raise ValueError('RAM too low: minimum 0.5 GB per CPU core')
            if ratio > 8:
                raise ValueError('RAM too high: maximum 8 GB per CPU core')
        return v
```

### Parsing gcloud JSON with jq
```bash
# Source: https://medium.com/google-cloud/bash-hacks-gcloud-kubectl-jq-etc-c2ff351d9c3b
# Get cheapest machine type in zone matching requirements
find_cheapest_machine_type() {
    local zone=$1
    local min_cpus=$2
    local min_memory_mb=$3

    gcloud compute machine-types list \
        --zones="$zone" \
        --format=json \
        --quiet \
    | jq -r \
        --arg cpus "$min_cpus" \
        --arg mem "$min_memory_mb" '
        [
            .[] |
            select(.guestCpus >= ($cpus | tonumber)) |
            select(.memoryMb >= ($mem | tonumber))
        ] |
        sort_by(.memoryMb) |
        first |
        .name
    '
}
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| configparser (INI files) | tomllib (TOML files) | Python 3.11 (Oct 2022) | Native support, better types, more readable |
| Pydantic v1 | Pydantic v2 | June 2023 | 5-50x faster validation, Rust core, better error messages |
| requests library | httpx | ~2021 adoption | Async support, HTTP/2, unified sync/async API |
| Manual config validation | BaseSettings (pydantic-settings) | Pydantic 2.0+ | Environment variable integration, multiple sources |
| gcloud beta machine-types | gcloud compute machine-types | ~2019 (GA) | Stable API, consistent filtering |

**Deprecated/outdated:**
- **configparser (INI)**: Still works but TOML is more expressive and type-aware
- **toml/tomli PyPI packages**: Obsolete for Python 3.11+, tomllib is stdlib
- **Pydantic v1 validators (@validator)**: Use @field_validator in v2
- **requests-cache for simple caching**: Overkill for single-endpoint 24-hour TTL cache

## Open Questions

Things that couldn't be fully resolved:

1. **CloudPrice.net API Exact Endpoint Format**
   - What we know: API exists at developer.cloudprice.net, requires API key, provides GCP spot pricing
   - What's unclear: Exact endpoint path, request/response JSON schema, rate limits, error codes
   - Recommendation: Implement hardcoded fallback first, test API integration with real API key. Mark as optional feature if API proves unreliable.

2. **GCP Machine Type Availability Validation**
   - What we know: gcloud machine-types list shows types per zone, availability changes dynamically
   - What's unclear: Best pattern for "try region 1, fallback to region 2" - should this be in pricing module or deployment orchestrator?
   - Recommendation: Phase 50 (Machine Selection) should handle fallback logic, Phase 49 just validates config and provides pricing data.

3. **Pricing Cache Invalidation Strategy**
   - What we know: 24-hour TTL works for stable pricing, but regional capacity changes aren't reflected
   - What's unclear: Should we invalidate cache on deployment failure? Check cache age before retry?
   - Recommendation: Start with simple time-based TTL, add manual cache clear command (--refresh-pricing flag) for users experiencing stale data issues.

4. **Config File Location**
   - What we know: Should be in gcp/ directory alongside scripts
   - What's unclear: Should config.toml be committed (with .example) or gitignored like current config.sh?
   - Recommendation: Gitignore config.toml, provide config.toml.example with comments. User copies and customizes.

5. **Error Message Verbosity**
   - What we know: Pydantic provides detailed validation errors
   - What's unclear: Should we show full Pydantic error output or simplify for end users?
   - Recommendation: Format errors into user-friendly messages (see Code Examples), but add --verbose flag to show full Pydantic output for debugging.

## Sources

### Primary (HIGH confidence)
- Python 3 tomllib documentation: https://docs.python.org/3/library/tomllib.html - Complete API reference
- Pydantic validators documentation: https://docs.pydantic.dev/latest/concepts/validators/ - Field validation patterns
- gcloud scripting guide: https://docs.cloud.google.com/sdk/docs/scripting-gcloud - Non-interactive best practices
- GCP regions and zones: https://docs.cloud.google.com/compute/docs/regions-zones - Official infrastructure documentation
- GCP Spot VMs: https://docs.cloud.google.com/compute/docs/instances/spot - Official spot instance documentation

### Secondary (MEDIUM confidence)
- CloudPrice.net Developer Portal: https://developer.cloudprice.net/ - API exists but specific endpoints require key
- GCP Spot VM pricing: https://cloud.google.com/spot-vms/pricing - Official pricing page (manual, not API)
- CloudPrice.net GCP pricing comparison: https://cloudprice.net/gcp/regions - Shows historical data
- Real Python TOML guide: https://realpython.com/python-toml/ - Practical tomllib examples
- Bash hacks gcloud + jq: https://medium.com/google-cloud/bash-hacks-gcloud-kubectl-jq-etc-c2ff351d9c3b - Community patterns

### Tertiary (LOW confidence)
- CloudPrice.net API structure: WebSearch results suggest v2 API with authentication but no public schema
- Hishel caching library: https://github.com/karpetrosyan/hishel - Alternative to simple file cache (not needed for this project)
- GCP pricing gotchas: https://www.wiz.io/academy/cloud-cost/gcp-cost - Third-party analysis (2026)

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - All libraries verified via official docs, already in project dependencies
- Architecture: HIGH - Patterns verified via official gcloud docs and Python stdlib docs
- Pitfalls: MEDIUM - Combination of official docs (gcloud limitations) and community experience (CloudPrice.net API reliability)
- CloudPrice.net API: LOW - Public docs limited, actual endpoint format requires API key access

**Research date:** 2026-02-06
**Valid until:** 2026-03-06 (30 days - stable domain: config parsing and validation patterns unlikely to change)

**Critical dependencies:**
- Python 3.11+ (tomllib available)
- Pydantic 2.12.5+ (already in pyproject.toml)
- httpx 0.28.0+ (already in pyproject.toml)
- gcloud CLI (already required for GCP operations)
- jq 1.7+ (typically pre-installed, brew install jq if needed)

**Implementation priority:**
1. Config validation (HIGH confidence, blocks everything)
2. Hardcoded regional fallback (HIGH confidence, required before API integration)
3. CloudPrice.net API integration (LOW confidence, optional enhancement)
4. Pricing cache (MEDIUM confidence, performance optimization)
