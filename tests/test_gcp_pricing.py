"""Tests for GCP spot pricing query module.

Tests cover:
- Hardcoded fallback regions (structure, coverage, ordering)
- File-based cache with 24-hour TTL
- CloudPrice.net API integration (mocked)
- Integration: cache -> API -> fallback cascade
- JSON output formatting for bash consumption
"""

import json
import time

import httpx
import pytest

from gcp.query_pricing import (
    CACHE_TTL,
    FALLBACK_REGIONS,
    format_pricing_output,
    get_ranked_regions,
    load_cache,
    query_cloudprice_api,
    save_cache,
)


# ---------------------------------------------------------------------------
# Hardcoded Fallback Tests
# ---------------------------------------------------------------------------


class TestFallbackRegions:
    """Tests for FALLBACK_REGIONS constant."""

    def test_fallback_regions_exist(self):
        """FALLBACK_REGIONS is a list of dicts with region, zone, rank keys."""
        assert isinstance(FALLBACK_REGIONS, list)
        assert len(FALLBACK_REGIONS) >= 10
        for entry in FALLBACK_REGIONS:
            assert isinstance(entry, dict)
            assert "region" in entry
            assert "zone" in entry
            assert "rank" in entry

    def test_fallback_regions_have_required_keys(self):
        """Each entry has exactly the expected keys."""
        required_keys = {"region", "zone", "rank"}
        for entry in FALLBACK_REGIONS:
            assert required_keys.issubset(entry.keys()), (
                f"Missing keys in {entry}"
            )

    def test_fallback_regions_sorted_by_rank(self):
        """List is sorted by rank ascending (cheapest first)."""
        ranks = [entry["rank"] for entry in FALLBACK_REGIONS]
        assert ranks == sorted(ranks), "FALLBACK_REGIONS must be sorted by rank"

    def test_fallback_covers_us_regions(self):
        """At least 3 US regions in fallback."""
        us_regions = [
            e for e in FALLBACK_REGIONS if e["region"].startswith("us-")
        ]
        assert len(us_regions) >= 3, (
            f"Expected >=3 US regions, got {len(us_regions)}"
        )

    def test_fallback_covers_europe(self):
        """At least 2 European regions in fallback."""
        eu_regions = [
            e for e in FALLBACK_REGIONS if e["region"].startswith("europe-")
        ]
        assert len(eu_regions) >= 2, (
            f"Expected >=2 EU regions, got {len(eu_regions)}"
        )

    def test_get_ranked_regions_returns_fallback(self, monkeypatch, tmp_path):
        """get_ranked_regions returns fallback list when API unavailable."""
        monkeypatch.setattr(
            httpx,
            "get",
            lambda *a, **kw: (_ for _ in ()).throw(
                httpx.ConnectError("Connection refused")
            ),
        )
        result = get_ranked_regions(
            cpu_cores=8, ram_gb=32, cache_dir=tmp_path
        )
        assert isinstance(result, list)
        assert len(result) > 0
        # Fallback entries have rank key instead of price_per_hour
        assert "region" in result[0]
        assert "zone" in result[0]


# ---------------------------------------------------------------------------
# Cache Tests
# ---------------------------------------------------------------------------


class TestCache:
    """Tests for file-based pricing cache with 24-hour TTL."""

    def test_save_and_load_cache(self, tmp_path):
        """save_cache writes JSON, load_cache reads it back."""
        cache_file = tmp_path / "spot-pricing.json"
        data = [
            {"region": "us-central1", "zone": "us-central1-a", "price_per_hour": 0.05}
        ]
        save_cache(data, cache_file)
        loaded = load_cache(cache_file)
        assert loaded is not None
        assert loaded == data

    def test_cache_fresh(self, tmp_path):
        """Cache written with current timestamp is considered fresh."""
        cache_file = tmp_path / "spot-pricing.json"
        data = [
            {"region": "us-central1", "zone": "us-central1-a", "price_per_hour": 0.05}
        ]
        save_cache(data, cache_file)
        result = load_cache(cache_file)
        assert result is not None, "Fresh cache should return data, not None"

    def test_cache_stale(self, tmp_path):
        """Cache written 25 hours ago is considered stale (returns None)."""
        cache_file = tmp_path / "spot-pricing.json"
        stale_data = {
            "cached_at": time.time() - (25 * 3600),  # 25 hours ago
            "regions": [
                {"region": "us-central1", "zone": "us-central1-a", "price_per_hour": 0.05}
            ],
        }
        cache_file.write_text(json.dumps(stale_data))
        result = load_cache(cache_file)
        assert result is None, "Stale cache (25h old) should return None"

    def test_cache_missing(self, tmp_path):
        """load_cache with nonexistent path returns None."""
        cache_file = tmp_path / "does-not-exist.json"
        result = load_cache(cache_file)
        assert result is None

    def test_cache_ttl_default(self):
        """Default TTL is 86400 seconds (24 hours)."""
        assert CACHE_TTL == 86400


# ---------------------------------------------------------------------------
# API Integration Tests (mocked)
# ---------------------------------------------------------------------------


class TestApiQuery:
    """Tests for CloudPrice.net API queries (all mocked)."""

    def test_query_api_success(self, monkeypatch):
        """Mock httpx.get to return pricing data, verify sorted by price."""

        class MockResponse:
            status_code = 200

            def json(self):
                return [
                    {
                        "region": "europe-west1",
                        "zone": "europe-west1-b",
                        "price": 0.08,
                    },
                    {
                        "region": "us-central1",
                        "zone": "us-central1-a",
                        "price": 0.05,
                    },
                    {
                        "region": "asia-east1",
                        "zone": "asia-east1-a",
                        "price": 0.12,
                    },
                ]

            def raise_for_status(self):
                pass

        monkeypatch.setattr(httpx, "get", lambda *a, **kw: MockResponse())
        result = query_cloudprice_api(cpu_cores=8, ram_gb=32)
        assert result is not None
        assert isinstance(result, list)
        assert len(result) > 0
        # Should be sorted by price ascending
        prices = [r["price_per_hour"] for r in result]
        assert prices == sorted(prices), "Results must be sorted by price ascending"

    def test_query_api_timeout(self, monkeypatch):
        """Mock httpx.get to raise TimeoutException. Returns None."""
        monkeypatch.setattr(
            httpx,
            "get",
            lambda *a, **kw: (_ for _ in ()).throw(
                httpx.TimeoutException("Timed out")
            ),
        )
        result = query_cloudprice_api(cpu_cores=8, ram_gb=32)
        assert result is None

    def test_query_api_http_error(self, monkeypatch):
        """Mock httpx.get to return 500 status. Returns None."""

        class MockResponse:
            status_code = 500

            def raise_for_status(self):
                raise httpx.HTTPStatusError(
                    "Server Error",
                    request=httpx.Request("GET", "https://example.com"),
                    response=self,
                )

        monkeypatch.setattr(httpx, "get", lambda *a, **kw: MockResponse())
        result = query_cloudprice_api(cpu_cores=8, ram_gb=32)
        assert result is None

    def test_query_api_auth_error(self, monkeypatch):
        """Mock httpx.get to return 401 status. Returns None."""

        class MockResponse:
            status_code = 401

            def raise_for_status(self):
                raise httpx.HTTPStatusError(
                    "Unauthorized",
                    request=httpx.Request("GET", "https://example.com"),
                    response=self,
                )

        monkeypatch.setattr(httpx, "get", lambda *a, **kw: MockResponse())
        result = query_cloudprice_api(cpu_cores=8, ram_gb=32)
        assert result is None


# ---------------------------------------------------------------------------
# Integration Tests
# ---------------------------------------------------------------------------


class TestIntegration:
    """Tests for get_ranked_regions orchestrator function."""

    def test_get_ranked_regions_uses_cache(self, monkeypatch, tmp_path):
        """Fresh cache is used; API is not called."""
        cache_file = tmp_path / "spot-pricing.json"
        cached_data = [
            {"region": "us-central1", "zone": "us-central1-a", "price_per_hour": 0.05}
        ]
        save_cache(cached_data, cache_file)

        # API should NOT be called -- raise if it is
        monkeypatch.setattr(
            httpx,
            "get",
            lambda *a, **kw: (_ for _ in ()).throw(
                AssertionError("API should not be called when cache is fresh")
            ),
        )

        result = get_ranked_regions(
            cpu_cores=8, ram_gb=32, cache_dir=tmp_path
        )
        assert result == cached_data

    def test_get_ranked_regions_refreshes_stale_cache(
        self, monkeypatch, tmp_path
    ):
        """Stale cache (25h old) triggers fresh API query."""
        cache_file = tmp_path / "spot-pricing.json"
        stale_data = {
            "cached_at": time.time() - (25 * 3600),
            "regions": [
                {"region": "old-region", "zone": "old-zone", "price_per_hour": 9.99}
            ],
        }
        cache_file.write_text(json.dumps(stale_data))

        # Mock API returns fresh data
        class MockResponse:
            status_code = 200

            def json(self):
                return [
                    {
                        "region": "us-central1",
                        "zone": "us-central1-a",
                        "price": 0.03,
                    },
                ]

            def raise_for_status(self):
                pass

        monkeypatch.setattr(httpx, "get", lambda *a, **kw: MockResponse())

        result = get_ranked_regions(
            cpu_cores=8, ram_gb=32, cache_dir=tmp_path
        )
        assert result is not None
        assert len(result) > 0
        assert result[0]["region"] == "us-central1"
        assert result[0]["price_per_hour"] == 0.03

    def test_get_ranked_regions_fallback_on_total_failure(
        self, monkeypatch, tmp_path
    ):
        """No cache + API error = fallback regions returned."""
        monkeypatch.setattr(
            httpx,
            "get",
            lambda *a, **kw: (_ for _ in ()).throw(
                httpx.ConnectError("Connection refused")
            ),
        )

        result = get_ranked_regions(
            cpu_cores=8, ram_gb=32, cache_dir=tmp_path
        )
        assert isinstance(result, list)
        assert len(result) >= 10
        # Should be fallback regions (have rank key)
        assert "rank" in result[0]

    def test_format_pricing_output(self):
        """format_pricing_output returns valid JSON array."""
        regions = [
            {
                "region": "us-central1",
                "zone": "us-central1-a",
                "price_per_hour": 0.05,
            },
            {
                "region": "europe-west1",
                "zone": "europe-west1-b",
                "price_per_hour": 0.08,
            },
        ]
        output = format_pricing_output(regions)
        # Must be valid JSON
        parsed = json.loads(output)
        assert isinstance(parsed, list)
        assert len(parsed) == 2
        # Each entry must have region, zone, price_per_hour
        for entry in parsed:
            assert "region" in entry
            assert "zone" in entry
            assert "price_per_hour" in entry
