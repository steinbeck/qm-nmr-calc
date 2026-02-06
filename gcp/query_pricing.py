"""GCP spot pricing query with CloudPrice.net API, caching, and hardcoded fallback.

This module queries spot instance pricing across GCP regions to find the cheapest
deployment target. It implements a three-tier resolution strategy:

    1. File-based cache (24-hour TTL) -- fastest, avoids redundant API calls
    2. CloudPrice.net API query -- live pricing data sorted by cost
    3. Hardcoded regional fallback -- always works, based on historical patterns

The module is both importable and runnable as a CLI script:

    python3 gcp/query_pricing.py --cpu-cores 8 --ram-gb 32
    python3 gcp/query_pricing.py --refresh  # bypass cache

Output is always valid JSON parseable by jq.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import httpx

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

FALLBACK_REGIONS: list[dict] = [
    {"region": "us-central1", "zone": "us-central1-a", "rank": 1},
    {"region": "us-east4", "zone": "us-east4-a", "rank": 2},
    {"region": "us-west1", "zone": "us-west1-a", "rank": 3},
    {"region": "us-east1", "zone": "us-east1-b", "rank": 4},
    {"region": "us-west4", "zone": "us-west4-a", "rank": 5},
    {"region": "europe-west1", "zone": "europe-west1-b", "rank": 6},
    {"region": "europe-west4", "zone": "europe-west4-a", "rank": 7},
    {"region": "asia-east1", "zone": "asia-east1-a", "rank": 8},
    {"region": "asia-southeast1", "zone": "asia-southeast1-a", "rank": 9},
    {"region": "southamerica-east1", "zone": "southamerica-east1-a", "rank": 10},
]

CACHE_TTL: int = 86400  # 24 hours in seconds

DEFAULT_CACHE_DIR: Path = Path("gcp/.cache")

CLOUDPRICE_API_BASE = "https://api.cloudprice.net/v1/gcp/compute"
CLOUDPRICE_TIMEOUT = 10.0  # seconds


# ---------------------------------------------------------------------------
# Cache Functions
# ---------------------------------------------------------------------------


def load_cache(cache_file: Path) -> list | None:
    """Read cached pricing data if file exists and is fresh (< 24h old).

    Returns the cached regions list if fresh, None if stale/missing/corrupt.
    """
    if not cache_file.exists():
        return None

    try:
        with open(cache_file) as f:
            data = json.load(f)
    except (json.JSONDecodeError, OSError):
        return None

    cached_at = data.get("cached_at", 0)
    age = time.time() - cached_at

    if age > CACHE_TTL:
        return None  # Stale

    return data.get("regions")


def save_cache(regions: list, cache_file: Path) -> None:
    """Write pricing data to cache file with current timestamp.

    Creates parent directories if they don't exist.
    """
    cache_file.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_file, "w") as f:
        json.dump(
            {"cached_at": time.time(), "regions": regions},
            f,
            indent=2,
        )


# ---------------------------------------------------------------------------
# API Functions
# ---------------------------------------------------------------------------


def query_cloudprice_api(cpu_cores: int, ram_gb: int) -> list | None:
    """Query CloudPrice.net API for GCP spot pricing.

    Returns a list of ``{"region", "zone", "price_per_hour"}`` dicts sorted by
    price ascending, or None on any error (timeout, HTTP error, auth error,
    connection error).

    The API key is read from the ``CLOUDPRICE_API_KEY`` environment variable
    (optional -- the API may work without authentication for basic queries).
    """
    headers: dict[str, str] = {}
    api_key = os.environ.get("CLOUDPRICE_API_KEY")
    if api_key:
        headers["X-API-Key"] = api_key

    params = {
        "type": "spot",
        "provider": "gcp",
        "min_cpu": str(cpu_cores),
        "min_ram_gb": str(ram_gb),
    }

    try:
        response = httpx.get(
            CLOUDPRICE_API_BASE,
            headers=headers,
            params=params,
            timeout=CLOUDPRICE_TIMEOUT,
        )
        response.raise_for_status()

        raw = response.json()

        # Normalize API response into our standard format
        results = []
        if isinstance(raw, list):
            for entry in raw:
                region = entry.get("region", "")
                zone = entry.get("zone", region + "-a" if region else "")
                price = entry.get("price") or entry.get("price_per_hour", 0.0)
                if region:
                    results.append(
                        {
                            "region": region,
                            "zone": zone,
                            "price_per_hour": float(price),
                        }
                    )

        # Sort by price ascending (cheapest first)
        results.sort(key=lambda x: x["price_per_hour"])
        return results if results else None

    except httpx.HTTPStatusError as exc:
        print(
            f"WARNING: Pricing API HTTP error {exc.response.status_code}",
            file=sys.stderr,
        )
        return None
    except httpx.HTTPError as exc:
        print(
            f"WARNING: Pricing API request failed: {exc}",
            file=sys.stderr,
        )
        return None
    except (ValueError, KeyError, TypeError) as exc:
        print(
            f"WARNING: Pricing API response parse error: {exc}",
            file=sys.stderr,
        )
        return None


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------


def get_ranked_regions(
    cpu_cores: int = 8,
    ram_gb: int = 32,
    cache_dir: Path | None = None,
    refresh: bool = False,
) -> list:
    """Return ranked list of GCP regions for spot deployment.

    Resolution order:
        1. Cached data (if fresh and not --refresh)
        2. CloudPrice.net API (live query)
        3. Hardcoded FALLBACK_REGIONS (always succeeds)

    Never fails -- always returns a non-empty list.
    """
    if cache_dir is None:
        cache_dir = DEFAULT_CACHE_DIR
    cache_file = cache_dir / "spot-pricing.json"

    # 1. Try cache (unless refresh requested)
    if not refresh:
        cached = load_cache(cache_file)
        if cached is not None:
            return cached

    # 2. Try API query
    api_result = query_cloudprice_api(cpu_cores, ram_gb)
    if api_result is not None:
        save_cache(api_result, cache_file)
        return api_result

    # 3. Fall back to hardcoded regions (always works)
    print(
        "WARNING: Using hardcoded regional fallback (pricing API unavailable)",
        file=sys.stderr,
    )
    return FALLBACK_REGIONS


def format_pricing_output(regions: list) -> str:
    """Return JSON string of the regions list, parseable by jq."""
    return json.dumps(regions, indent=2)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main() -> None:
    """CLI entry point for pricing queries."""
    parser = argparse.ArgumentParser(
        description="Query GCP spot pricing across regions"
    )
    parser.add_argument(
        "--cpu-cores",
        type=int,
        default=8,
        help="Minimum CPU cores (default: 8)",
    )
    parser.add_argument(
        "--ram-gb",
        type=int,
        default=32,
        help="Minimum RAM in GB (default: 32)",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=DEFAULT_CACHE_DIR,
        help="Cache directory (default: gcp/.cache)",
    )
    parser.add_argument(
        "--refresh",
        action="store_true",
        help="Bypass cache and force fresh API query",
    )

    args = parser.parse_args()

    regions = get_ranked_regions(
        cpu_cores=args.cpu_cores,
        ram_gb=args.ram_gb,
        cache_dir=args.cache_dir,
        refresh=args.refresh,
    )

    print(format_pricing_output(regions))


if __name__ == "__main__":
    main()
