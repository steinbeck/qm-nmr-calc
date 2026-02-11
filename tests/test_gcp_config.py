"""Tests for GCP TOML config validation (Phase 49, Plan 01)."""

import sys
from pathlib import Path

import pytest
from pydantic import ValidationError

# Add gcp directory to path for imports
gcp_path = Path(__file__).parent.parent / "gcp"
sys.path.insert(0, str(gcp_path))

from validate_config import GCPConfig, format_exports, load_config


# --- GCPConfig model validation ---


class TestGCPConfigValid:
    """Tests for valid configurations."""

    def test_valid_config(self):
        """Valid config with all fields parses successfully."""
        config = GCPConfig(
            project_id="my-test-project",
            cpu_cores=8,
            ram_gb=32,
            disk_size_gb=100,
            resource_prefix="qm-nmr-calc",
        )
        assert config.project_id == "my-test-project"
        assert config.cpu_cores == 8
        assert config.ram_gb == 32
        assert config.disk_size_gb == 100
        assert config.resource_prefix == "qm-nmr-calc"

    def test_valid_config_minimal(self):
        """Only required fields work, defaults applied."""
        config = GCPConfig(
            project_id="my-test-project",
            cpu_cores=8,
            ram_gb=32,
        )
        assert config.resource_prefix == "qm-nmr-calc"
        assert config.disk_size_gb == 100


class TestGCPConfigProjectId:
    """Tests for project_id validation."""

    def test_missing_project_id(self):
        """Config without project_id raises ValidationError."""
        with pytest.raises(ValidationError):
            GCPConfig(cpu_cores=8, ram_gb=32)

    def test_project_id_too_short(self):
        """project_id less than 6 chars raises ValidationError."""
        with pytest.raises(ValidationError):
            GCPConfig(project_id="ab", cpu_cores=8, ram_gb=32)

    def test_project_id_too_long(self):
        """project_id with 31+ chars raises ValidationError."""
        with pytest.raises(ValidationError):
            GCPConfig(project_id="a" * 31, cpu_cores=8, ram_gb=32)

    def test_project_id_uppercase(self):
        """project_id with uppercase raises ValidationError."""
        with pytest.raises(ValidationError, match="(?i)lowercase"):
            GCPConfig(project_id="My-Project", cpu_cores=8, ram_gb=32)

    def test_project_id_starts_with_hyphen(self):
        """project_id starting with hyphen raises ValidationError."""
        with pytest.raises(ValidationError):
            GCPConfig(project_id="-my-project", cpu_cores=8, ram_gb=32)

    def test_project_id_ends_with_hyphen(self):
        """project_id ending with hyphen raises ValidationError."""
        with pytest.raises(ValidationError):
            GCPConfig(project_id="my-project-", cpu_cores=8, ram_gb=32)


class TestGCPConfigRamCpuRatio:
    """Tests for RAM/CPU ratio validation."""

    def test_ram_cpu_ratio_too_low(self):
        """CPU 64 + RAM 16 (ratio 0.25) raises ValidationError."""
        with pytest.raises(ValidationError):
            GCPConfig(project_id="my-test-project", cpu_cores=64, ram_gb=16)

    def test_ram_cpu_ratio_too_high(self):
        """CPU 4 + RAM 64 (ratio 16) raises ValidationError."""
        with pytest.raises(ValidationError):
            GCPConfig(project_id="my-test-project", cpu_cores=4, ram_gb=64)

    def test_ram_cpu_ratio_boundary_low(self):
        """CPU 8 + RAM 4 (ratio 0.5, boundary) passes."""
        config = GCPConfig(
            project_id="my-test-project", cpu_cores=8, ram_gb=4
        )
        assert config.ram_gb == 4

    def test_ram_cpu_ratio_boundary_high(self):
        """CPU 4 + RAM 32 (ratio 8, boundary) passes."""
        config = GCPConfig(
            project_id="my-test-project", cpu_cores=4, ram_gb=32
        )
        assert config.ram_gb == 32


class TestGCPConfigBounds:
    """Tests for field boundary validation."""

    def test_cpu_cores_too_low(self):
        """cpu_cores below min (4) raises ValidationError."""
        with pytest.raises(ValidationError):
            GCPConfig(project_id="my-test-project", cpu_cores=2, ram_gb=8)

    def test_disk_size_too_small(self):
        """disk_size_gb below min (10) raises ValidationError."""
        with pytest.raises(ValidationError):
            GCPConfig(
                project_id="my-test-project",
                cpu_cores=8,
                ram_gb=32,
                disk_size_gb=5,
            )


# --- TOML file loading ---


class TestLoadToml:
    """Tests for load_config() TOML file loading."""

    def test_load_toml_file(self, tmp_path):
        """Valid TOML with [gcp] section loads correctly."""
        config_file = tmp_path / "config.toml"
        config_file.write_text(
            '[gcp]\n'
            'project_id = "test-project-123"\n'
            'cpu_cores = 8\n'
            'ram_gb = 32\n'
            'disk_size_gb = 200\n'
        )
        config = load_config(config_file)
        assert isinstance(config, GCPConfig)
        assert config.project_id == "test-project-123"
        assert config.cpu_cores == 8
        assert config.ram_gb == 32
        assert config.disk_size_gb == 200

    def test_load_toml_missing_section(self, tmp_path):
        """TOML file without [gcp] section raises ValueError."""
        config_file = tmp_path / "config.toml"
        config_file.write_text('project_id = "test-project-123"\n')
        with pytest.raises(ValueError, match="\\[gcp\\] section"):
            load_config(config_file)

    def test_load_toml_invalid_syntax(self, tmp_path):
        """TOML file with broken syntax raises ValueError."""
        config_file = tmp_path / "config.toml"
        config_file.write_text("[gcp\nbroken = !!!\n")
        with pytest.raises(ValueError):
            load_config(config_file)

    def test_load_toml_file_not_found(self, tmp_path):
        """Non-existent path raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            load_config(tmp_path / "nonexistent.toml")


# --- Export formatting ---


class TestFormatExports:
    """Tests for format_exports() bash export generation."""

    def test_format_exports(self):
        """Exports contain all required bash export statements."""
        config = GCPConfig(
            project_id="my-test-project",
            cpu_cores=8,
            ram_gb=32,
            disk_size_gb=100,
            resource_prefix="qm-nmr-calc",
        )
        exports = format_exports(config)
        assert "export GCP_PROJECT_ID=" in exports
        assert "export CPU_CORES=" in exports
        assert "export RAM_GB=" in exports
        assert "export DISK_SIZE_GB=" in exports
        assert "export RESOURCE_PREFIX=" in exports
