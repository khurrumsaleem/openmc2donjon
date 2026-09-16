from __future__ import annotations

import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
import pytest

from openmc2donjon.openmc_hex_current_export import (
    _check_paired_net_mean,
    _manifest_definition_sha256,
    _paired_statistics,
    _statistics,
    export_hex_currents,
)


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _fixture(tmp_path: Path, nearly_equal: bool = False) -> tuple[Path, list[Path], np.ndarray, np.ndarray]:
    signs = np.array([1, 1, 1, -1, -1, -1, -1, 1])
    faces = [
        {"index": i, "name": f"F{i}", "surface_id": 100 + i, "outward_sign": int(sign),
         "outward_normal": [float(sign), 0.0, 0.0], "area_cm2": float(i + 2),
         "boundary_type": "transmission"}
        for i, sign in enumerate(signs)
    ]
    manifest = {
        "schema": "openmc2donjon.hex-face-current.v1", "energy_bounds_eV": [0.0, 1.0, 10.0],
        "nodes": [{"node": [0, 0, 0], "cell_id": 10, "cell_name": "test", "faces": faces,
                   "tallies": {"outgoing": {"id": 1}, "incoming": {"id": 2}}}],
        "input_files": {},
    }
    for role in ("geometry", "materials", "settings", "tallies"):
        path = tmp_path / f"{role}.xml"
        content = ("<settings><run_mode>fixed source</run_mode><source strength='1' "
                   "type='independent'/></settings>") if role == "settings" else f"<{role}/>"
        path.write_text(content)
        manifest["input_files"][role] = {"path": path.name, "sha256": _digest(path)}
    out = np.broadcast_to(np.array([0.0, 1.0, 0.0, 1.0])[:, None, None, None], (4, 1, 8, 2)).copy()
    incoming = 1.0 - out
    if nearly_equal:
        out = np.broadcast_to(np.random.default_rng(3).random((4, 1, 1, 1)) * 1e-3, out.shape).copy()
        incoming = out[::-1].copy()
    out[..., 1] *= 2.0
    incoming[..., 1] *= 3.0
    paths = []
    for n in range(1, 5):
        path = tmp_path / f"statepoint.{n}.h5"
        with h5py.File(path, "w") as handle:
            handle.attrs["filetype"] = np.bytes_("statepoint")
            handle.attrs["photon_transport"] = False
            handle.attrs["openmc_version"] = [0, 15, 4]
            handle.attrs["version"] = [18, 2]
            for key, value in {"n_realizations": n, "current_batch": n, "run_mode": np.bytes_("fixed source"),
                               "energy_mode": np.bytes_("continuous-energy"), "seed": 12, "n_particles": 100,
                               "n_batches": 4, "stride": 152917, "generations_per_batch": 1}.items():
                handle.create_dataset(key, data=value)
            filters = handle.create_group("tallies/filters")
            for fid, kind, bins, count in ((1, "surface", list(range(100, 108)), 8),
                                           (2, "cellfrom", [10], 1), (3, "energy", [0.0, 1.0, 10.0], 2),
                                           (4, "cell", [10], 1)):
                group = filters.create_group(f"filter {fid}")
                group.create_dataset("type", data=np.bytes_(kind))
                group.create_dataset("bins", data=bins)
                group.create_dataset("n_bins", data=count)
            for tid, samples, orientation, filter_ids in ((1, out, signs, [1, 2, 3]),
                                                          (2, incoming, -signs, [1, 4, 3])):
                tally = handle.create_group(f"tallies/tally {tid}")
                tally.create_dataset("n_realizations", data=n)
                tally.create_dataset("score_bins", data=[np.bytes_("current")])
                tally.create_dataset("nuclides", data=[np.bytes_("total")])
                tally.create_dataset("estimator", data=np.bytes_("analog"))
                tally.create_dataset("filters", data=filter_ids)
                total = (samples[:n].sum(axis=0)[0] * orientation[:, None]).reshape(16, 1)
                squared = (samples[:n] ** 2).sum(axis=0)[0].reshape(16, 1)
                tally.create_dataset("results", data=np.stack([total, squared], axis=-1))
        paths.append(path)
    receipt = {"returncode": 0, "timed_out": False, "input_files": manifest["input_files"],
               "command": ["openmc", "--threads", "1"],
               "manifest_definition_sha256": _manifest_definition_sha256(manifest),
               "output_files": {path.name: _digest(path) for path in paths}}
    receipt_path = tmp_path / "receipt.json"
    receipt_path.write_text(json.dumps(receipt))
    manifest["execution_receipt"] = {"path": receipt_path.name, "sha256": _digest(receipt_path)}
    manifest_path = tmp_path / "currents.json"
    manifest_path.write_text(json.dumps(manifest))
    return manifest_path, paths, out, incoming


def _rebind(manifest_path: Path, paths: list[Path]) -> None:
    manifest = json.loads(manifest_path.read_text())
    for record in manifest["input_files"].values():
        record["sha256"] = _digest(manifest_path.parent / record["path"])
    receipt_path = manifest_path.parent / manifest["execution_receipt"]["path"]
    receipt = json.loads(receipt_path.read_text())
    receipt["input_files"] = manifest["input_files"]
    receipt["output_files"] = {path.name: _digest(path) for path in paths}
    receipt["manifest_definition_sha256"] = _manifest_definition_sha256(manifest)
    receipt_path.write_text(json.dumps(receipt))
    manifest["execution_receipt"]["sha256"] = _digest(receipt_path)
    manifest_path.write_text(json.dumps(manifest))


def test_paired_statistics_retains_negative_covariance() -> None:
    out = np.array([[0.0], [1.0], [0.0], [1.0]])
    incoming = 1.0 - out
    result = _paired_statistics(out.cumsum(0), incoming.cumsum(0), (out**2).sum(0), (incoming**2).sum(0))
    expected = (out - incoming).std(axis=0, ddof=1) / 2.0
    np.testing.assert_allclose(result["net_std_dev"], expected)
    _, out_std = _statistics(out.sum(0), (out**2).sum(0), 4)
    np.testing.assert_allclose(result["net_std_dev"], 2.0 * out_std)
    assert result["out_in_covariance_of_mean"][0] < 0


def test_zero_signed_mean_can_have_nonzero_standard_error() -> None:
    mean, error = _statistics(np.array([0.0]), np.array([4.0]), 4)
    assert mean[0] == 0
    assert error[0] > 0


def test_nearly_cancelled_net_retains_canonical_native_mean(tmp_path: Path) -> None:
    manifest, paths, out, incoming = _fixture(tmp_path, nearly_equal=True)
    output = tmp_path / "current.h5"
    export_hex_currents(paths, manifest, output)
    canonical = out.sum(axis=0) / len(out) - incoming.sum(axis=0) / len(incoming)
    paired = _paired_statistics(out.cumsum(axis=0), incoming.cumsum(axis=0),
                                (out**2).sum(axis=0), (incoming**2).sum(axis=0))
    assert np.any(paired["net_mean"] != canonical)
    with h5py.File(output) as handle:
        np.testing.assert_array_equal(handle["currents/net_mean"], canonical[..., ::-1])
        np.testing.assert_allclose(handle["currents/net_std_dev"], paired["net_std_dev"][..., ::-1])


def test_net_consistency_bound_rejects_non_roundoff_errors() -> None:
    outgoing, incoming = np.array([4e-4]), np.array([4e-4])
    canonical = outgoing - incoming
    _check_paired_net_mean(np.array([1e-20]), canonical, outgoing, incoming, 30)
    with pytest.raises(ValueError, match="does not match"):
        _check_paired_net_mean(np.array([1e-14]), canonical, outgoing, incoming, 30)


def test_complete_export_and_energy_density_mapping(tmp_path: Path) -> None:
    manifest, paths, out, incoming = _fixture(tmp_path)
    output = tmp_path / "current.h5"
    summary = export_hex_currents(paths[::-1], manifest, output, tmp_path / "summary.json")
    assert summary["net_std_available"] is True
    assert summary["precision_accepted"] is False
    with h5py.File(output) as handle:
        np.testing.assert_allclose(handle["currents/out_mean"], out.mean(axis=0)[..., ::-1])
        np.testing.assert_allclose(handle["currents/in_mean"], incoming.mean(axis=0)[..., ::-1])
        net = out - incoming
        np.testing.assert_allclose(handle["currents/net_std_dev"], net.std(axis=0, ddof=1)[..., ::-1] / 2.0)
        np.testing.assert_allclose(handle["current_density/net_mean"],
                                   net.mean(axis=0)[..., ::-1] / np.arange(2.0, 10.0)[None, :, None])
        np.testing.assert_array_equal(handle["group_bounds_eV"], [[1.0, 10.0], [0.0, 1.0]])
        assert handle.attrs["energy_order"] == "high-to-low"
    assert json.loads((tmp_path / "summary.json").read_text())["output_h5"]["sha256"] == _digest(output)


@pytest.mark.parametrize("alias", ["current.h5", "nested/../current.h5", "linked/current.h5"])
def test_output_and_summary_alias_rejected_without_writing(tmp_path: Path, alias: str) -> None:
    manifest, paths, _, _ = _fixture(tmp_path)
    (tmp_path / "nested").mkdir()
    (tmp_path / "linked").symlink_to(tmp_path, target_is_directory=True)
    output = tmp_path / "current.h5"
    with pytest.raises(ValueError, match="distinct output paths"):
        export_hex_currents(paths, manifest, output, tmp_path / alias)
    assert not output.exists()


def test_final_only_omits_net_error(tmp_path: Path) -> None:
    manifest, paths, _, _ = _fixture(tmp_path)
    output = tmp_path / "current.h5"
    summary = export_hex_currents(paths[-1:], manifest, output)
    assert summary["net_std_available"] is False
    with h5py.File(output) as handle:
        assert "currents/net_std_dev" not in handle
        assert "currents/out_std_dev" in handle


@pytest.mark.parametrize("selected", [[0, 2, 3], [1, 2, 3], [0, 1, 3]])
def test_incomplete_checkpoint_series_rejected(tmp_path: Path, selected: list[int]) -> None:
    manifest, paths, _, _ = _fixture(tmp_path)
    with pytest.raises(ValueError, match="every active realization"):
        export_hex_currents([paths[i] for i in selected], manifest, tmp_path / "current.h5")


def test_checkpoint_substitution_fails_marginal_second_moment(tmp_path: Path) -> None:
    manifest, paths, _, _ = _fixture(tmp_path)
    with h5py.File(paths[1], "r+") as handle:
        results = handle["tallies/tally 1/results"]
        values = results[()]
        values[:, 0, 0] *= 0.5
        results[...] = values
    _rebind(manifest, paths)
    with pytest.raises(ValueError, match="sum_sq"):
        export_hex_currents(paths, manifest, tmp_path / "current.h5")


def test_statepoint_receipt_hash_is_required(tmp_path: Path) -> None:
    manifest, paths, _, _ = _fixture(tmp_path)
    with h5py.File(paths[-1], "r+") as handle:
        handle["seed"][...] = 42
    with pytest.raises(ValueError, match="not bound"):
        export_hex_currents(paths, manifest, tmp_path / "current.h5")


def test_different_seed_series_is_rejected_even_when_rebound(tmp_path: Path) -> None:
    manifest, paths, _, _ = _fixture(tmp_path)
    with h5py.File(paths[1], "r+") as handle:
        handle["seed"][...] = 42
    _rebind(manifest, paths)
    with pytest.raises(ValueError, match="inconsistent run identities"):
        export_hex_currents(paths, manifest, tmp_path / "current.h5")


def test_native_energy_bins_must_match_exactly(tmp_path: Path) -> None:
    manifest, paths, _, _ = _fixture(tmp_path)
    with h5py.File(paths[-1], "r+") as handle:
        handle["tallies/filters/filter 3/bins"][...] = [0.0, 2.0, 10.0]
    _rebind(manifest, paths)
    with pytest.raises(ValueError, match="filter bins differ"):
        export_hex_currents(paths[-1:], manifest, tmp_path / "current.h5")


@pytest.mark.parametrize("setting", ["<no_reduce>true</no_reduce>", "<source strength='1'/>"])
def test_unsupported_normalization_is_rejected(tmp_path: Path, setting: str) -> None:
    manifest, paths, _, _ = _fixture(tmp_path)
    settings = tmp_path / "settings.xml"
    settings.write_text(settings.read_text().replace("</settings>", setting + "</settings>"))
    _rebind(manifest, paths)
    with pytest.raises(ValueError, match="no_reduce|unit-normalized"):
        export_hex_currents(paths, manifest, tmp_path / "current.h5")


def test_wrong_orientation_is_not_hidden_with_abs(tmp_path: Path) -> None:
    manifest, paths, _, _ = _fixture(tmp_path)
    payload = json.loads(manifest.read_text())
    payload["nodes"][0]["faces"][0]["outward_sign"] = -1
    manifest.write_text(json.dumps(payload))
    _rebind(manifest, paths)
    with pytest.raises(ValueError, match="partial current is negative"):
        export_hex_currents(paths, manifest, tmp_path / "current.h5")


def test_refuse_existing_output(tmp_path: Path) -> None:
    manifest, paths, _, _ = _fixture(tmp_path)
    output = tmp_path / "current.h5"
    output.write_bytes(b"user content")
    with pytest.raises(FileExistsError):
        export_hex_currents(paths, manifest, output)
    assert output.read_bytes() == b"user content"


@pytest.mark.parametrize("field,value", [("area_cm2", 100.0), ("outward_normal", [0.0, 1.0, 0.0]),
                                        ("surface_id", 999)])
def test_postrun_manifest_edits_are_rejected(tmp_path: Path, field: str, value: object) -> None:
    manifest, paths, _, _ = _fixture(tmp_path)
    payload = json.loads(manifest.read_text())
    payload["nodes"][0]["faces"][0][field] = value
    manifest.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="complete current manifest definition"):
        export_hex_currents(paths, manifest, tmp_path / "current.h5")


@pytest.mark.parametrize("setting", ["<photon_transport>true</photon_transport>",
                                    "<source strength='0' particle='photon'/>"])
def test_photon_settings_rejected_before_statepoint(tmp_path: Path, setting: str) -> None:
    manifest, paths, _, _ = _fixture(tmp_path)
    settings = tmp_path / "settings.xml"
    settings.write_text(settings.read_text().replace("</settings>", setting + "</settings>"))
    _rebind(manifest, paths)
    with pytest.raises(ValueError, match="photon transport settings|neutron-only sources"):
        export_hex_currents(paths, manifest, tmp_path / "current.h5")


@pytest.mark.parametrize("command", [["openmc", "--no-reduce"], "openmc --no-reduce"])
def test_command_line_no_reduce_rejected(tmp_path: Path, command: object) -> None:
    manifest, paths, _, _ = _fixture(tmp_path)
    receipt = tmp_path / "receipt.json"
    payload = json.loads(receipt.read_text())
    payload["command"] = command
    receipt.write_text(json.dumps(payload))
    _rebind(manifest, paths)
    with pytest.raises(ValueError, match="Command-line --no-reduce"):
        export_hex_currents(paths, manifest, tmp_path / "current.h5")


def test_statepoint_photon_dataset_rejected(tmp_path: Path) -> None:
    manifest, paths, _, _ = _fixture(tmp_path)
    with h5py.File(paths[-1], "r+") as handle:
        handle.create_dataset("photon_transport", data=True)
    _rebind(manifest, paths)
    with pytest.raises(ValueError, match="photon transport runs"):
        export_hex_currents(paths[-1:], manifest, tmp_path / "current.h5")
