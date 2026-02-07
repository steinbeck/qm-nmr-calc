"""Tests for GCP machine selection and resource calculation (Phase 50, Plan 01)."""

import json
from io import StringIO
from unittest.mock import patch

import pytest

from gcp.select_machine import (
    calculate_docker_resources,
    find_available_zone,
    generate_startup_script,
    select_machine_type,
    validate_machine_in_zone,
)


# --- select_machine_type tests ---


class TestSelectMachineType:
    """Tests for select_machine_type function."""

    def test_select_machine_type_returns_name(self):
        """Mock gcloud to return valid JSON, assert returns machine type name."""
        mock_response = json.dumps(
            [{"name": "n2-standard-8", "guestCpus": 8, "memoryMb": 32768}]
        )
        with patch("gcp.select_machine._run_gcloud", return_value=mock_response):
            result = select_machine_type(8, 32, "us-central1-a")
            assert result == "n2-standard-8"

    def test_select_machine_type_empty_result(self):
        """Mock gcloud to return empty list, assert raises ValueError."""
        mock_response = json.dumps([])
        with patch("gcp.select_machine._run_gcloud", return_value=mock_response):
            with pytest.raises(ValueError, match="No machine types found"):
                select_machine_type(8, 32, "us-central1-a")

    def test_select_machine_type_gcloud_failure(self):
        """Mock gcloud to raise RuntimeError, assert raises RuntimeError."""
        with patch(
            "gcp.select_machine._run_gcloud",
            side_effect=RuntimeError("gcloud failed"),
        ):
            with pytest.raises(RuntimeError, match="gcloud failed"):
                select_machine_type(8, 32, "us-central1-a")


# --- validate_machine_in_zone tests ---


class TestValidateMachineInZone:
    """Tests for validate_machine_in_zone function."""

    def test_validate_machine_exists(self):
        """Mock gcloud to return machine name, assert returns True."""
        with patch("gcp.select_machine._run_gcloud", return_value="n2-standard-8"):
            result = validate_machine_in_zone("n2-standard-8", "us-central1-a")
            assert result is True

    def test_validate_machine_missing(self):
        """Mock gcloud to return empty string, assert returns False."""
        with patch("gcp.select_machine._run_gcloud", return_value=""):
            result = validate_machine_in_zone("n2-standard-8", "us-central1-a")
            assert result is False


# --- find_available_zone tests ---


class TestFindAvailableZone:
    """Tests for find_available_zone function."""

    def test_find_available_zone_first_succeeds(self):
        """Mock gcloud to succeed on first region, assert returns zone info."""
        mock_machine_response = json.dumps(
            [{"name": "n2-standard-8", "guestCpus": 8, "memoryMb": 32768}]
        )
        mock_regions = [
            {"zone": "us-central1-a", "region": "us-central1"},
            {"zone": "us-east4-a", "region": "us-east4"},
        ]

        with patch("gcp.select_machine._run_gcloud", return_value=mock_machine_response):
            with patch(
                "gcp.select_machine.get_ranked_regions", return_value=mock_regions
            ):
                result = find_available_zone(8, 32)
                assert "zone" in result
                assert "region" in result
                assert "machine_type" in result
                assert result["zone"] == "us-central1-a"
                assert result["machine_type"] == "n2-standard-8"

    def test_find_available_zone_fallback_to_second(self):
        """Mock gcloud to fail first, succeed second, assert returns second zone."""
        mock_machine_response = json.dumps(
            [{"name": "n2-standard-8", "guestCpus": 8, "memoryMb": 32768}]
        )
        mock_regions = [
            {"zone": "us-central1-a", "region": "us-central1"},
            {"zone": "us-east4-a", "region": "us-east4"},
        ]

        # First call raises ValueError, second call succeeds
        with patch(
            "gcp.select_machine._run_gcloud",
            side_effect=[
                ValueError("No machine types found"),
                mock_machine_response,
            ],
        ):
            with patch(
                "gcp.select_machine.get_ranked_regions", return_value=mock_regions
            ):
                result = find_available_zone(8, 32)
                assert result["zone"] == "us-east4-a"

    def test_find_available_zone_all_exhausted(self):
        """Mock gcloud to always fail, assert raises RuntimeError."""
        mock_regions = [{"zone": "us-central1-a", "region": "us-central1"}]

        with patch(
            "gcp.select_machine._run_gcloud",
            side_effect=ValueError("No machine types found"),
        ):
            with patch(
                "gcp.select_machine.get_ranked_regions", return_value=mock_regions
            ):
                with pytest.raises(RuntimeError, match="All regions exhausted"):
                    find_available_zone(8, 32)

    def test_find_available_zone_calls_get_ranked_regions(self):
        """Mock both functions, assert get_ranked_regions called with correct params."""
        mock_machine_response = json.dumps(
            [{"name": "n2-standard-8", "guestCpus": 8, "memoryMb": 32768}]
        )
        mock_regions = [{"zone": "us-central1-a", "region": "us-central1"}]

        with patch("gcp.select_machine._run_gcloud", return_value=mock_machine_response):
            with patch(
                "gcp.select_machine.get_ranked_regions", return_value=mock_regions
            ) as mock_get_ranked:
                find_available_zone(8, 32)
                mock_get_ranked.assert_called_once_with(cpu_cores=8, ram_gb=32)


# --- calculate_docker_resources tests ---


class TestCalculateDockerResources:
    """Tests for calculate_docker_resources function."""

    def test_calculate_resources_standard(self):
        """Mock gcloud describe to return 8 CPU / 32GB RAM, assert correct calculation."""
        mock_response = json.dumps({"guestCpus": 8, "memoryMb": 32768})
        with patch("gcp.select_machine._run_gcloud", return_value=mock_response):
            result = calculate_docker_resources("n2-standard-8", "us-central1-a")
            assert result["worker_memory_limit"] == "24g"
            assert result["nwchem_nproc"] == 8
            assert result["total_ram_gb"] == 32
            assert result["total_cpus"] == 8

    def test_calculate_resources_highmem(self):
        """Mock gcloud describe to return 4 CPU / 32GB RAM, assert correct calculation."""
        mock_response = json.dumps({"guestCpus": 4, "memoryMb": 32768})
        with patch("gcp.select_machine._run_gcloud", return_value=mock_response):
            result = calculate_docker_resources("n2-highmem-4", "us-central1-a")
            assert result["worker_memory_limit"] == "24g"
            assert result["nwchem_nproc"] == 4

    def test_calculate_resources_insufficient_ram(self):
        """Mock gcloud to return 8GB total RAM, assert raises ValueError (insufficient after overhead)."""
        mock_response = json.dumps({"guestCpus": 2, "memoryMb": 8192})
        with patch("gcp.select_machine._run_gcloud", return_value=mock_response):
            with pytest.raises(ValueError, match="Insufficient RAM"):
                calculate_docker_resources("e2-small", "us-central1-a")


# --- generate_startup_script tests ---


class TestGenerateStartupScript:
    """Tests for generate_startup_script function."""

    def test_startup_script_contains_memory_limit(self):
        """Assert output contains WORKER_MEMORY_LIMIT=24g."""
        script = generate_startup_script("24g", "qm-nmr-calc", 100)
        assert "WORKER_MEMORY_LIMIT=24g" in script

    def test_startup_script_uses_nproc(self):
        """Assert output contains $(nproc) for runtime CPU detection."""
        script = generate_startup_script("24g", "qm-nmr-calc", 100)
        assert "$(nproc)" in script

    def test_startup_script_contains_docker_compose(self):
        """Assert output contains 'docker compose'."""
        script = generate_startup_script("24g", "qm-nmr-calc", 100)
        assert "docker compose" in script

    def test_startup_script_has_strict_mode(self):
        """Assert output contains 'set -euo pipefail'."""
        script = generate_startup_script("24g", "qm-nmr-calc", 100)
        assert "set -euo pipefail" in script

    def test_startup_script_disables_caddy(self):
        """Assert docker-compose override disables Caddy service."""
        script = generate_startup_script("24g", "qm-nmr-calc", 100)
        assert "caddy: null" in script


# --- CLI tests ---


class TestCLI:
    """Tests for CLI functionality."""

    def test_cli_help(self):
        """Run main() with --help, assert shows help text."""
        from gcp.select_machine import main

        with pytest.raises(SystemExit) as exc_info:
            with patch("sys.argv", ["select_machine.py", "--help"]):
                main()
        assert exc_info.value.code == 0

    def test_cli_json_output(self):
        """Mock gcloud and get_ranked_regions, capture stdout, assert valid JSON output."""
        from gcp.select_machine import main

        mock_machine_response = json.dumps(
            [{"name": "n2-standard-8", "guestCpus": 8, "memoryMb": 32768}]
        )
        mock_describe_response = json.dumps({"guestCpus": 8, "memoryMb": 32768})
        mock_regions = [{"zone": "us-central1-a", "region": "us-central1"}]

        with patch(
            "gcp.select_machine._run_gcloud",
            side_effect=[mock_machine_response, mock_describe_response],
        ):
            with patch(
                "gcp.select_machine.get_ranked_regions", return_value=mock_regions
            ):
                with patch("sys.argv", ["select_machine.py"]):
                    with patch("sys.stdout", new_callable=StringIO) as mock_stdout:
                        main()
                        output = mock_stdout.getvalue()
                        # Assert valid JSON
                        parsed = json.loads(output)
                        assert "zone" in parsed
                        assert "machine_type" in parsed
                        assert "resources" in parsed
