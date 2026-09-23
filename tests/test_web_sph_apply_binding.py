"""Web SPH results retain the numerical report's binding and artifact type."""

import json

import h5py
import pytest

from openmc2donjon.openmc_provenance import file_sha256
from openmc2donjon.physical_sph_contract import physical_sph_issues
from openmc2donjon.web.server import create_app
from tests.test_sph_apply import _write_mgxs, _write_openmc_native_mgxs, _write_sidecar

pytest.importorskip("httpx")
TestClient = pytest.importorskip("fastapi.testclient").TestClient


@pytest.mark.parametrize("mode,input_format,verified", [
    ("converter-final-exact-input", "converter", True),
    ("converter-unbound", "converter", False),
    ("openmc-mgxs-intermediate-unbound", "openmc-mgxs", False),
])
def test_live_apply_response_matches_hdf5_and_summary(tmp_path, mode, input_format, verified):
    source, native, sidecar, output, summary = (tmp_path / name for name in (
        "in.h5", "native.h5", "sph.h5", "out.h5", "summary.json",
    ))
    _write_mgxs(source)
    _write_openmc_native_mgxs(native)
    _write_sidecar(sidecar)
    if mode != "converter-unbound":
        with h5py.File(sidecar, "r+") as h5:
            # Native setN iteration bytes are intentionally not this bound input.
            h5.attrs["sph_input_h5_sha256"] = file_sha256(source)
    input_path = source if input_format == "converter" else native
    with TestClient(create_app(workspace_root=tmp_path)) as client:
        response = client.post("/api/execute/apply-sph", json={
            "input_h5": str(input_path), "sph_source": str(sidecar),
            "input_format": input_format, "output_path": str(output),
            "summary_json": str(summary), "force": False,
        })
    assert response.status_code == 200, response.text
    result = response.json()
    assert result["ok"] is True
    assert result["mock_mode"] is False
    assert result["input_format"] == input_format
    assert result["binding_mode"] == mode
    assert result["sidecar_input_hash_verified"] is verified
    saved_summary = json.loads(summary.read_text())
    for name in ("input_format", "binding_mode", "sidecar_input_hash_verified"):
        assert result[name] == saved_summary[name]
    with h5py.File(output) as h5:
        assert result["binding_mode"] == h5.attrs["sph_apply_binding_mode"]
        stored_verified = bool(h5.attrs["sph_apply_sidecar_input_hash_verified"])
        assert result["sidecar_input_hash_verified"] is stored_verified
    if input_format == "converter":
        # These arithmetic fixtures lack the full CE/MG provenance. Even a
        # matching input hash must not be presented as physical acceptance.
        assert physical_sph_issues(output)


@pytest.mark.parametrize("input_format", ["converter", "openmc-mgxs"])
def test_mock_apply_never_claims_a_verified_binding(input_format):
    with TestClient(create_app(mock_mode=True)) as client:
        response = client.post("/api/execute/apply-sph", json={
            "input_h5": "/mock/in.h5", "sph_source": "/mock/sph.h5",
            "input_format": input_format, "output_path": "/mock/out.h5", "force": False,
        })
    assert response.status_code == 200
    result = response.json()
    assert result["mock_mode"] is True
    assert result["input_format"] == input_format
    assert result["binding_mode"] == "mock-unverified"
    assert result["sidecar_input_hash_verified"] is False


@pytest.mark.parametrize("digest", ["malformed", "0" * 64])
def test_web_does_not_downgrade_invalid_binding_to_success(tmp_path, digest):
    source, sidecar, output = (tmp_path / name for name in ("in.h5", "sph.h5", "out.h5"))
    _write_mgxs(source)
    _write_sidecar(sidecar)
    with h5py.File(sidecar, "r+") as h5:
        h5.attrs["sph_input_h5_sha256"] = digest
    with TestClient(create_app(workspace_root=tmp_path)) as client:
        response = client.post("/api/execute/apply-sph", json={
            "input_h5": str(source), "sph_source": str(sidecar),
            "input_format": "converter", "output_path": str(output), "force": False,
        })
    assert response.status_code == 422
    assert "SHA-256" in response.json()["detail"]
    assert not output.exists()
