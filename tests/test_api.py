"""Integration tests for the QM NMR Calculator API."""
from unittest import mock

import pytest
from fastapi.testclient import TestClient

from qm_nmr_calc.api.app import app

client = TestClient(app)


class TestHealth:
    """Health endpoint tests."""

    def test_liveness(self):
        """GET /health returns alive status."""
        response = client.get("/health")
        assert response.status_code == 200
        assert response.json()["status"] == "alive"

    def test_readiness(self):
        """GET /health/ready checks dependencies."""
        response = client.get("/health/ready")
        # May be 200 or 503 depending on environment
        assert response.status_code in [200, 503]
        assert "status" in response.json()

    def test_health_ready_includes_crest_available(self):
        """GET /health/ready includes crest_available boolean field."""
        # Mock CREST detection to ensure predictable result
        with mock.patch("qm_nmr_calc.conformers.crest_generator.detect_crest_available") as mock_detect:
            # Test with CREST available
            mock_detect.return_value = True
            response = client.get("/health/ready")
            data = response.json()
            assert "crest_available" in data
            assert data["crest_available"] is True
            assert "checks" in data
            assert data["checks"]["crest_available"] is True

            # Test with CREST not available
            mock_detect.return_value = False
            response = client.get("/health/ready")
            data = response.json()
            assert "crest_available" in data
            assert data["crest_available"] is False
            assert "checks" in data
            assert data["checks"]["crest_available"] is False


class TestJobSubmission:
    """Job submission endpoint tests."""

    def test_submit_valid_smiles(self):
        """POST /api/v1/jobs with valid SMILES returns 202."""
        response = client.post(
            "/api/v1/jobs",
            json={"smiles": "CCO", "name": "Ethanol", "solvent": "chcl3"}
        )
        assert response.status_code == 202
        data = response.json()
        assert "job_id" in data
        assert data["status"] == "queued"
        assert data["input_smiles"] == "CCO"
        assert data["input_name"] == "Ethanol"
        assert data["preset"] == "production"
        assert data["solvent"] == "chcl3"
        # Check headers
        assert "Location" in response.headers
        assert "Retry-After" in response.headers

    def test_submit_invalid_smiles(self):
        """POST /api/v1/jobs with invalid SMILES returns 422."""
        response = client.post(
            "/api/v1/jobs",
            json={"smiles": "not-a-valid-smiles-string", "solvent": "chcl3"}
        )
        assert response.status_code == 422
        data = response.json()["detail"]
        assert data["status"] == 422
        assert "title" in data  # RFC 7807

    def test_submit_smiles_without_name(self):
        """POST /api/v1/jobs works without optional name."""
        response = client.post(
            "/api/v1/jobs",
            json={"smiles": "c1ccccc1", "solvent": "dmso"}
        )
        assert response.status_code == 202
        data = response.json()
        assert data["input_name"] is None


class TestJobStatus:
    """Job status endpoint tests."""

    def test_get_existing_job(self):
        """GET /api/v1/jobs/{id} returns status for existing job."""
        # First create a job
        create_response = client.post(
            "/api/v1/jobs",
            json={"smiles": "CCO", "solvent": "chcl3"}
        )
        job_id = create_response.json()["job_id"]

        # Then get its status
        response = client.get(f"/api/v1/jobs/{job_id}")
        assert response.status_code == 200
        data = response.json()
        assert data["job_id"] == job_id
        assert data["status"] == "queued"

    def test_get_nonexistent_job(self):
        """GET /api/v1/jobs/{id} returns 404 for unknown job."""
        response = client.get("/api/v1/jobs/nonexistent123")
        assert response.status_code == 404
        data = response.json()["detail"]
        assert data["status"] == 404
        assert "title" in data  # RFC 7807


class TestFileUpload:
    """File upload endpoint tests."""

    def test_upload_mol_file(self):
        """POST /api/v1/jobs/upload with MOL file returns 202."""
        # Simple MOL file content for ethane
        mol_content = """
  ethane

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5400    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
"""
        response = client.post(
            "/api/v1/jobs/upload",
            files={"file": ("ethane.mol", mol_content, "chemical/x-mdl-molfile")},
            data={"name": "Ethane", "solvent": "chcl3"}
        )
        assert response.status_code == 202
        data = response.json()
        assert "job_id" in data
        assert data["status"] == "queued"

    def test_upload_invalid_file_type(self):
        """POST /api/v1/jobs/upload with wrong extension returns 422."""
        response = client.post(
            "/api/v1/jobs/upload",
            files={"file": ("test.txt", "not a molecule", "text/plain")}
        )
        assert response.status_code == 422


class TestOpenAPI:
    """OpenAPI documentation tests."""

    def test_docs_available(self):
        """GET /docs returns Swagger UI."""
        response = client.get("/docs")
        assert response.status_code == 200
        assert "swagger" in response.text.lower()

    def test_openapi_json(self):
        """GET /api/v1/openapi.json returns OpenAPI spec."""
        response = client.get("/api/v1/openapi.json")
        assert response.status_code == 200
        data = response.json()
        assert data["info"]["title"] == "QM NMR Calculator API"
        assert "/api/v1/jobs" in str(data["paths"])


class TestEnsembleSchemas:
    """Tests for ensemble API schema extensions."""

    def test_ensemble_metadata_response_model(self):
        """Verify EnsembleMetadataResponse model can be instantiated."""
        from qm_nmr_calc.api.schemas import EnsembleMetadataResponse

        metadata = EnsembleMetadataResponse(
            conformer_count=5,
            total_generated=20,
            method="rdkit_kdg",
            temperature_k=298.15,
            energy_range_kcal=2.5,
            top_populations=[
                {"id": "conf_001", "population": 0.45, "energy_kcal": 0.0},
                {"id": "conf_002", "population": 0.30, "energy_kcal": 0.5},
            ],
        )
        assert metadata.conformer_count == 5
        assert metadata.method == "rdkit_kdg"

    def test_conformer_progress_response_model(self):
        """Verify ConformerProgressResponse model can be instantiated."""
        from qm_nmr_calc.api.schemas import ConformerProgressResponse

        progress = ConformerProgressResponse(
            conformer_id="conf_001",
            status="nmr_complete",
            energy_kcal=0.0,
            population=0.45,
        )
        assert progress.conformer_id == "conf_001"
        assert progress.status == "nmr_complete"

    def test_job_status_response_ensemble_fields(self):
        """Verify JobStatusResponse has ensemble fields."""
        from qm_nmr_calc.api.schemas import JobStatusResponse

        fields = JobStatusResponse.model_fields
        assert "conformer_mode" in fields
        assert "conformer_method" in fields
        assert "conformer_count" in fields
        assert "conformer_progress" in fields
        assert "ensemble_metadata" in fields
        assert "conformer_method_warning" in fields

    def test_nmr_results_response_ensemble_metadata(self):
        """Verify NMRResultsResponse has ensemble_metadata field."""
        from qm_nmr_calc.api.schemas import NMRResultsResponse

        fields = NMRResultsResponse.model_fields
        assert "ensemble_metadata" in fields


class TestGeometryEndpoint:
    """Tests for geometry.json endpoint."""

    def test_geometry_endpoint_requires_existing_job(self):
        """GET /api/v1/jobs/{id}/geometry.json returns 404 for unknown job."""
        response = client.get("/api/v1/jobs/nonexistent123/geometry.json")
        assert response.status_code == 404
        data = response.json()["detail"]
        assert data["status"] == 404

    def test_geometry_endpoint_includes_conformer_mode(self):
        """GET /api/v1/jobs/{id}/geometry.json includes conformer_mode field."""
        # Create a single-conformer job
        create_response = client.post(
            "/api/v1/jobs",
            json={"smiles": "CCO", "solvent": "chcl3", "conformer_mode": "single"}
        )
        job_id = create_response.json()["job_id"]

        response = client.get(f"/api/v1/jobs/{job_id}/geometry.json")
        assert response.status_code == 200
        data = response.json()
        assert "conformer_mode" in data
        assert data["conformer_mode"] == "single"
        assert "conformers" in data
        assert data["conformers"] is None

    def test_geometry_endpoint_single_conformer_structure(self):
        """Verify single conformer geometry response has expected fields."""
        # Create job
        create_response = client.post(
            "/api/v1/jobs",
            json={"smiles": "C", "solvent": "chcl3", "conformer_mode": "single"}
        )
        job_id = create_response.json()["job_id"]

        response = client.get(f"/api/v1/jobs/{job_id}/geometry.json")
        assert response.status_code == 200
        data = response.json()

        # Check required fields
        assert "job_id" in data
        assert "status" in data
        assert "xyz" in data
        assert "sdf" in data
        assert "h1_assignments" in data
        assert "c13_assignments" in data
        assert "conformer_mode" in data
        assert "conformers" in data

    def test_geometry_endpoint_function_signature(self):
        """Verify get_geometry_data function exists with correct signature."""
        from qm_nmr_calc.api.routers.jobs import get_geometry_data
        import inspect

        sig = inspect.signature(get_geometry_data)
        assert "job_id" in sig.parameters

    def test_xyz_to_sdf_helper_exists(self):
        """Verify _xyz_to_sdf helper function exists."""
        from qm_nmr_calc.api.routers.jobs import _xyz_to_sdf

        # Test with valid input
        xyz_content = """5
Methane
C   0.0000   0.0000   0.0000
H   1.0890   0.0000   0.0000
H  -0.3630   1.0270   0.0000
H  -0.3630  -0.5135   0.8892
H  -0.3630  -0.5135  -0.8892
"""
        smiles = "C"
        result = _xyz_to_sdf(xyz_content, smiles)
        assert result is not None
        assert "M  END" in result  # SDF/MOL block terminator

    def test_xyz_to_sdf_helper_handles_invalid_smiles(self):
        """Verify _xyz_to_sdf returns None for invalid SMILES."""
        from qm_nmr_calc.api.routers.jobs import _xyz_to_sdf

        xyz_content = "2\ntest\nC 0 0 0\nH 1 0 0\n"
        result = _xyz_to_sdf(xyz_content, "not-valid-smiles")
        assert result is None

    def test_conformer_data_schema_fields(self):
        """Verify expected conformer data fields in ensemble response."""
        # This is a schema verification test
        # Full integration test would require running NWChem
        expected_conformer_fields = {"id", "xyz", "sdf", "energy_kcal", "population"}

        # Verify by examining the endpoint code
        from qm_nmr_calc.api.routers.jobs import get_geometry_data
        import inspect

        source = inspect.getsource(get_geometry_data)
        # Check that all expected fields are in the response construction
        for field in expected_conformer_fields:
            assert f'"{field}"' in source, f"Field '{field}' not found in conformer data"


class TestNMReDataEndpoint:
    """Tests for NMReData SDF download endpoint."""

    def test_nmredata_endpoint_returns_404_for_nonexistent_job(self):
        """GET /api/v1/jobs/{id}/nmredata.sdf returns 404 for unknown job."""
        response = client.get("/api/v1/jobs/nonexistent123/nmredata.sdf")
        assert response.status_code == 404
        data = response.json()["detail"]
        assert data["status"] == 404
        assert "title" in data  # RFC 7807

    def test_nmredata_endpoint_returns_409_for_incomplete_job(self):
        """GET /api/v1/jobs/{id}/nmredata.sdf returns 409 for incomplete job."""
        # Create a job (will be in "queued" state)
        create_response = client.post(
            "/api/v1/jobs",
            json={"smiles": "CCO", "solvent": "chcl3"}
        )
        job_id = create_response.json()["job_id"]

        # Immediately try to download NMReData (job still queued)
        response = client.get(f"/api/v1/jobs/{job_id}/nmredata.sdf")
        assert response.status_code == 409
        data = response.json()["detail"]
        assert data["status"] == 409
        assert "title" in data  # RFC 7807
