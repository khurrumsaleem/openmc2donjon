"""Export native hex-face current tallies, retaining paired-batch covariance.

This is a current diagnostic sidecar, not a MGXS library, SPH correction,
surface-flux estimate, or accepted physical model.  A final statepoint alone
does not contain the out/in cross moment required for a net-current error bar.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shlex
from typing import Any, Sequence
import xml.etree.ElementTree as ET

import h5py
import numpy as np


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _text(value: Any) -> str:
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8")
    return str(value)


def _manifest_definition_sha256(manifest: dict[str, Any]) -> str:
    definition = {key: value for key, value in manifest.items() if key != "execution_receipt"}
    encoded = json.dumps(definition, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    return hashlib.sha256(encoded).hexdigest()


def _statistics(total: np.ndarray, squared: np.ndarray, count: int) -> tuple[np.ndarray, np.ndarray]:
    """Return mean and standard error using native normalized batch moments."""
    if count < 2:
        raise ValueError("At least two active realizations are required for current error bars")
    total = np.asarray(total, dtype=float)
    squared = np.asarray(squared, dtype=float)
    if total.shape != squared.shape or not np.all(np.isfinite(total)) or not np.all(np.isfinite(squared)):
        raise ValueError("Invalid native tally moments")
    mean = total / count
    numerator = squared / count - mean**2
    tolerance = 5e-12 * np.maximum(np.maximum(np.abs(squared / count), mean**2), 1e-30)
    if np.any(numerator < -tolerance):
        raise ValueError("Native tally moments imply negative variance")
    return mean, np.sqrt(np.maximum(numerator, 0.0) / (count - 1))


def _paired_statistics(
    cumulative_out: np.ndarray,
    cumulative_in: np.ndarray,
    final_squared_out: np.ndarray,
    final_squared_in: np.ndarray,
) -> dict[str, np.ndarray]:
    """Reconstruct paired normalized batch samples from complete running sums.

    The input has active realization as its first axis; out/in have already
    been converted to non-negative partial-current magnitudes.  Both signs
    must be transformed before differencing, never inferred with abs().
    """
    out = np.asarray(cumulative_out, dtype=float)
    incoming = np.asarray(cumulative_in, dtype=float)
    if out.shape != incoming.shape or out.ndim < 2 or len(out) < 2:
        raise ValueError("Paired current samples require matching complete batch arrays")
    out_batches = np.diff(out, axis=0, prepend=np.zeros_like(out[:1]))
    in_batches = np.diff(incoming, axis=0, prepend=np.zeros_like(incoming[:1]))
    for label, samples, squares in (
        ("outgoing", out_batches, final_squared_out),
        ("incoming", in_batches, final_squared_in),
    ):
        if not np.all(np.isfinite(samples)):
            raise ValueError(f"Non-finite {label} current batch values")
        tolerance = 5e-12 * np.maximum(np.max(np.abs(samples), axis=0), 1e-30)
        if np.any(samples < -tolerance):
            raise ValueError(f"Negative {label} partial current: invalid orientation or checkpoint series")
        reconstructed = np.sum(samples**2, axis=0)
        if not np.allclose(reconstructed, squares, rtol=2e-9, atol=1e-25):
            raise ValueError(f"{label} batch series does not reproduce final native sum_sq")
    net_batches = out_batches - in_batches
    return {
        "net_mean": net_batches.mean(axis=0),
        "net_std_dev": net_batches.std(axis=0, ddof=1) / np.sqrt(len(net_batches)),
        "out_in_covariance_of_mean": np.sum(
            (out_batches - out_batches.mean(axis=0)) * (in_batches - in_batches.mean(axis=0)),
            axis=0,
        ) / (len(out) * (len(out) - 1)),
    }


def _check_paired_net_mean(
    paired_mean: np.ndarray,
    canonical_mean: np.ndarray,
    outgoing_mean: np.ndarray,
    incoming_mean: np.ndarray,
    count: int,
) -> None:
    """Bound only floating-point reconstruction error, not a physics tolerance.

    Net currents may nearly cancel. Their own magnitude is therefore not a
    valid roundoff scale. Differencing N nonnegative cumulative partial sums,
    subtracting the paired partials, and averaging requires O(N) operations.
    The conservative gamma_(4N+8) bound is applied to the *partial* magnitude;
    the subnormal term covers underflow of the same number of operations.
    """
    operations = 4 * count + 8
    epsilon = np.finfo(float).eps
    if count < 2 or operations * epsilon >= 1:
        raise ValueError("Invalid realization count for floating-point consistency bound")
    gamma = operations * epsilon / (1.0 - operations * epsilon)
    scale = np.abs(outgoing_mean) + np.abs(incoming_mean)
    tolerance = gamma * scale + operations * np.nextafter(0.0, 1.0)
    difference = np.abs(paired_mean - canonical_mean)
    if not np.all(np.isfinite(difference)) or np.any(difference > tolerance):
        raise ValueError("Paired-batch net mean does not match final native partial-current means")


def _verified_artifact(record: dict[str, Any], base: Path) -> Path:
    if not isinstance(record, dict) or not record.get("path") or not record.get("sha256"):
        raise ValueError("Each bound artifact requires a path and SHA-256 digest")
    path = Path(record["path"])
    if not path.is_absolute():
        path = base / path
    if _sha256(path) != record["sha256"]:
        raise ValueError(f"Artifact checksum mismatch: {path}")
    return path.resolve()


def _verify_provenance(manifest: dict[str, Any], base: Path, statepoints: list[Path]) -> dict[str, Any]:
    inputs = manifest.get("input_files", {})
    resolved: dict[str, Path] = {}
    for role in ("geometry", "materials", "settings", "tallies"):
        record = inputs.get(role, inputs.get(f"{role}.xml"))
        if record is None:
            raise ValueError(f"Missing bound input artifact: {role}")
        resolved[role] = _verified_artifact(record, base)
    receipt_path = _verified_artifact(manifest.get("execution_receipt", {}), base)
    receipt = json.loads(receipt_path.read_text())
    if receipt.get("returncode") != 0 or receipt.get("timed_out", False):
        raise ValueError("Current export requires a successfully completed execution receipt")
    definition_hash = _manifest_definition_sha256(manifest)
    if receipt.get("manifest_definition_sha256") != definition_hash:
        raise ValueError("Execution receipt does not bind the complete current manifest definition")
    command = receipt.get("command", [])
    command_args = shlex.split(command) if isinstance(command, str) else command
    if any(str(argument).split("=", 1)[0] == "--no-reduce" for argument in command_args):
        raise ValueError("Command-line --no-reduce tally realizations are not supported")
    outputs = receipt.get("output_files", {})
    hashes = {}
    for path in statepoints:
        expected = outputs.get(path.name, outputs.get(str(path)))
        if isinstance(expected, dict):
            expected = expected.get("sha256")
        actual = _sha256(path)
        if actual != expected:
            raise ValueError(f"Statepoint is not bound to the execution receipt: {path}")
        hashes[str(path)] = actual
    receipt_inputs = receipt.get("input_files", {})
    for role, path in resolved.items():
        expected = receipt_inputs.get(role, receipt_inputs.get(f"{role}.xml"))
        if isinstance(expected, dict):
            expected = expected.get("sha256")
        if _sha256(path) != expected:
            raise ValueError(f"Execution receipt does not bind the declared {role} input")
    settings = ET.parse(resolved["settings"]).getroot()
    if settings.findtext("no_reduce", "false").strip().lower() in ("true", "1"):
        raise ValueError("Rank-local no_reduce tally realizations are not supported")
    if settings.findtext("photon_transport", "false").strip().lower() in ("true", "1"):
        raise ValueError("Neutron-only current export does not accept photon transport settings")
    sources = settings.findall("source")
    if any(source.get("particle", source.findtext("particle", "neutron")) != "neutron" for source in sources):
        raise ValueError("Current export requires neutron-only sources")
    run_mode = settings.findtext("run_mode", "eigenvalue").strip()
    if run_mode == "fixed source":
        strengths = np.array([float(source.get("strength", "1.0")) for source in sources])
        strength = strengths.sum()
        if not sources or not np.isclose(strength, 1.0, rtol=0.0, atol=1e-12):
            raise ValueError("Fixed-source currents require an explicitly unit-normalized source strength")
        if np.any(strengths < 0) or any(source.get("particle", "neutron") != "neutron" for source in sources):
            raise ValueError("Current export requires nonnegative neutron source strengths")
        if any(source.get("type", "independent") != "independent" for source in sources):
            raise ValueError("Only independent fixed sources have verified source normalization here")
    elif run_mode != "eigenvalue":
        raise ValueError(f"Unsupported current run mode: {run_mode}")
    return {
        "input_files": {
            role: {"path": str(path), "sha256": _sha256(path)} for role, path in resolved.items()
        },
        "execution_receipt": {"path": str(receipt_path), "sha256": _sha256(receipt_path)},
        "manifest_definition_sha256": definition_hash,
        "statepoints": hashes,
        "run_mode": run_mode,
        "source_normalization": "per unit starting-source neutron; no power normalization",
    }


def _validate_manifest(manifest: dict[str, Any]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if manifest.get("schema") != "openmc2donjon.hex-face-current.v1":
        raise ValueError("Unsupported hex-face current manifest schema")
    bounds = np.asarray(manifest.get("energy_bounds_eV", []), dtype=float)
    if bounds.ndim != 1 or bounds.size < 2 or not np.all(np.isfinite(bounds)):
        raise ValueError("Invalid current energy boundaries")
    if bounds[0] < 0 or np.any(np.diff(bounds) <= 0):
        raise ValueError("Current energy boundaries must be nonnegative and strictly increasing")
    nodes = manifest.get("nodes", [])
    if not nodes or len({int(node["cell_id"]) for node in nodes}) != len(nodes):
        raise ValueError("Current domains must have distinct cell IDs")
    signs, areas = [], []
    for node in nodes:
        faces = node["faces"]
        if len(faces) != 8 or [face["index"] for face in faces] != list(range(8)):
            raise ValueError("Each hexagonal node requires eight ordered faces")
        if len({face["surface_id"] for face in faces}) != 8:
            raise ValueError("Each hexagonal node requires eight distinct bounding surfaces")
        normals = np.asarray([face["outward_normal"] for face in faces], dtype=float)
        if normals.shape != (8, 3) or not np.allclose(np.linalg.norm(normals, axis=1), 1.0):
            raise ValueError("Face outward normals must be finite unit vectors")
        if any(face.get("boundary_type") != "transmission" for face in faces):
            raise ValueError("Only transmission current faces are supported")
        signs.append([face["outward_sign"] for face in faces])
        areas.append([face["area_cm2"] for face in faces])
    signs_array, areas_array = np.asarray(signs, dtype=float), np.asarray(areas, dtype=float)
    if not np.all(np.isin(signs_array, (-1, 1))):
        raise ValueError("Each face requires outward_sign +1 or -1")
    if not np.all(np.isfinite(areas_array)) or np.any(areas_array <= 0):
        raise ValueError("Each face requires a positive finite area")
    return bounds, signs_array[..., None], areas_array[..., None]


def _read_checkpoint(path: Path, manifest: dict[str, Any], bounds: np.ndarray) -> dict[str, Any]:
    shape = (len(manifest["nodes"]), 8, len(bounds) - 1)
    sums = {direction: np.zeros(shape) for direction in ("outgoing", "incoming")}
    squares = {direction: np.zeros(shape) for direction in sums}
    with h5py.File(path, "r") as handle:
        if _text(handle.attrs.get("filetype", "")) != "statepoint":
            raise ValueError(f"Not an OpenMC statepoint: {path}")
        photon_dataset = handle["photon_transport"][()] if "photon_transport" in handle else False
        if bool(handle.attrs.get("photon_transport", False)) or bool(photon_dataset):
            raise ValueError("Neutron-only current export does not accept photon transport runs")
        count = int(handle["n_realizations"][()])
        header: dict[str, Any] = {}
        for key in ("run_mode", "energy_mode", "seed", "n_particles", "n_batches", "stride",
                    "generations_per_batch", "n_inactive"):
            if key in handle:
                value = handle[key][()]
                header[key] = _text(value) if isinstance(value, (bytes, np.bytes_)) else int(value)
        header["openmc_version"] = np.asarray(handle.attrs["openmc_version"]).tolist()
        header["statepoint_version"] = np.asarray(handle.attrs["version"]).tolist()
        current_batch = int(handle["current_batch"][()])
        if count != current_batch - int(header.get("n_inactive", 0)):
            raise ValueError("Tally realizations do not correspond to consecutive active batches")
        for node_index, node in enumerate(manifest["nodes"]):
            surface_ids = [int(face["surface_id"]) for face in node["faces"]]
            for direction in sums:
                definition = node["tallies"][direction]
                tally = handle[f"tallies/tally {int(definition['id'])}"]
                if int(tally["n_realizations"][()]) != count:
                    raise ValueError("Current tallies have different realization counts")
                if [_text(value) for value in tally["score_bins"][()]] != ["current"]:
                    raise ValueError("Expected current-only tally")
                if [_text(value) for value in tally["nuclides"][()]] != ["total"]:
                    raise ValueError("Current tally must score total particles, not nuclide bins")
                if _text(tally["estimator"][()]) != "analog":
                    raise ValueError("Expected analog current tally")
                filter_ids = np.asarray(tally["filters"][()]).tolist()
                expected_types = ["surface", "cellfrom" if direction == "outgoing" else "cell", "energy"]
                if len(filter_ids) != 3:
                    raise ValueError("Expected exactly surface, cell and energy current filters")
                filters = [handle[f"tallies/filters/filter {int(value)}"] for value in filter_ids]
                if [_text(item["type"][()]) for item in filters] != expected_types:
                    raise ValueError("Native current tally filter types or order differ from the manifest")
                expected_bins = (
                    np.asarray(surface_ids, dtype=int),
                    np.asarray([int(node["cell_id"])], dtype=int),
                    bounds,
                )
                for item, expected in zip(filters, expected_bins, strict=True):
                    if not np.array_equal(item["bins"][()], expected):
                        raise ValueError("Native current filter bins differ from the manifest")
                expected_counts = [8, 1, len(bounds) - 1]
                if [int(item["n_bins"][()]) for item in filters] != expected_counts:
                    raise ValueError("Native current filter dimensions differ from the manifest")
                results = np.asarray(tally["results"][()], dtype=float)
                if results.shape != (8 * (len(bounds) - 1), 1, 2) or not np.all(np.isfinite(results)):
                    raise ValueError("Invalid native current tally result shape or values")
                sums[direction][node_index] = results[:, 0, 0].reshape(8, len(bounds) - 1)
                squares[direction][node_index] = results[:, 0, 1].reshape(8, len(bounds) - 1)
    return {"path": str(path), "count": count, "header": header, "sum": sums, "sum_sq": squares}


def export_hex_currents(
    statepoints: Sequence[str | Path],
    manifest_path: str | Path,
    output_h5: str | Path,
    summary_json: str | Path | None = None,
) -> dict[str, Any]:
    """Write a separate, provenance-bound native hex-current diagnostic HDF5.

    Supply one final statepoint for partial means/errors and net means only.
    Supply every active-batch statepoint (realizations 1 through N) to obtain
    covariance-aware net errors. An incomplete multi-file series is rejected.
    All files must be hashed in a successful execution receipt referenced by
    the tally manifest. Fixed sources must have total source strength one.
    """
    manifest_path, output_h5 = Path(manifest_path).resolve(), Path(output_h5).resolve()
    summary_path = Path(summary_json).resolve() if summary_json is not None else None
    if summary_path == output_h5:
        raise ValueError("Current sidecar and summary must use distinct output paths")
    if output_h5.exists() or (summary_path is not None and summary_path.exists()):
        raise FileExistsError("Refusing to overwrite an existing current sidecar or summary")
    paths = [Path(path).resolve() for path in statepoints]
    if not paths or len(set(paths)) != len(paths):
        raise ValueError("Specify a final statepoint or distinct complete active-batch statepoints")
    manifest = json.loads(manifest_path.read_text())
    bounds, signs, areas = _validate_manifest(manifest)
    provenance = _verify_provenance(manifest, manifest_path.parent, paths)
    checkpoints = sorted(
        (_read_checkpoint(path, manifest, bounds) for path in paths), key=lambda x: x["count"]
    )
    final = checkpoints[-1]
    if final["header"]["run_mode"] != provenance["run_mode"]:
        raise ValueError("Statepoint run mode differs from bound settings")
    if any(item["header"] != final["header"] for item in checkpoints):
        raise ValueError("Statepoints have inconsistent run identities")
    n = final["count"]
    out_mean, out_std = _statistics(signs * final["sum"]["outgoing"], final["sum_sq"]["outgoing"], n)
    in_mean, in_std = _statistics(-signs * final["sum"]["incoming"], final["sum_sq"]["incoming"], n)
    if np.any(out_mean < 0) or np.any(in_mean < 0):
        raise ValueError("Converted partial current is negative; check normals and tally selection")
    arrays = {"out_mean": out_mean, "out_std_dev": out_std, "in_mean": in_mean,
              "in_std_dev": in_std, "net_mean": out_mean - in_mean}
    complete = len(checkpoints) > 1
    if complete:
        if [item["count"] for item in checkpoints] != list(range(1, n + 1)):
            raise ValueError("Net-current uncertainty requires every active realization from 1 through N")
        paired = _paired_statistics(
            np.stack([signs * item["sum"]["outgoing"] for item in checkpoints]),
            np.stack([-signs * item["sum"]["incoming"] for item in checkpoints]),
            final["sum_sq"]["outgoing"], final["sum_sq"]["incoming"],
        )
        _check_paired_net_mean(paired.pop("net_mean"), arrays["net_mean"], out_mean, in_mean, n)
        arrays.update(paired)
    summary: dict[str, Any] = {
        "schema": "openmc2donjon.hex-face-current-results.v1",
        "status": "diagnostic_not_physics_acceptance",
        "shape": list(out_mean.shape), "array_axes": ["node", "face", "energy_group"],
        "energy_order": "high-to-low", "active_realizations": n,
        "net_std_available": complete,
        "net_std_method": "paired normalized active-batch samples" if complete else None,
        "net_std_unavailable_reason": None if complete else "final statepoint lacks out/in covariance",
        "units": "weighted neutron crossings per starting-source neutron",
        "current_density_units": "weighted neutron crossings per starting-source neutron per cm^2",
        "zero_out_bins": int(np.count_nonzero(out_mean == 0)),
        "zero_in_bins": int(np.count_nonzero(in_mean == 0)),
        "precision_accepted": False,
        "caveat": ("Zero bins are unscored, not validated zero physics. "
                   "Batch errors do not prove convergence."),
        "run": final["header"], "provenance": provenance,
        "manifest": {"path": str(manifest_path), "sha256": _sha256(manifest_path)},
    }
    output_h5.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output_h5, "x") as handle:
        handle.attrs["schema"] = summary["schema"]
        handle.attrs["status"] = summary["status"]
        handle.attrs["net_std_available"] = complete
        handle.attrs["energy_order"] = "high-to-low"
        handle.attrs["array_axes"] = "node,face,energy_group"
        handle.attrs["summary_json"] = json.dumps(summary, sort_keys=True)
        handle.attrs["manifest_json"] = json.dumps(manifest, sort_keys=True)
        handle.create_dataset("energy_bounds_eV", data=bounds)
        handle.create_dataset("group_bounds_eV", data=np.column_stack((bounds[:-1], bounds[1:]))[::-1])
        nodes = manifest["nodes"]
        handle.create_dataset("nodes/indices", data=[node["node"] for node in nodes])
        handle.create_dataset("nodes/cell_ids", data=[node["cell_id"] for node in nodes])
        strings = h5py.string_dtype("utf-8")
        handle.create_dataset("nodes/names", data=[node["cell_name"] for node in nodes], dtype=strings)
        handle.create_dataset("faces/names", data=[[face["name"] for face in node["faces"]]
                                                   for node in nodes], dtype=strings)
        handle.create_dataset("faces/surface_ids", data=[[face["surface_id"] for face in node["faces"]]
                                                         for node in nodes])
        handle.create_dataset("faces/outward_signs", data=signs[..., 0])
        handle.create_dataset("faces/outward_normals", data=[
            [face["outward_normal"] for face in node["faces"]] for node in nodes
        ])
        handle.create_dataset("faces/areas_cm2", data=areas[..., 0])
        for name, values in arrays.items():
            dataset = handle.create_dataset(f"currents/{name}", data=values[..., ::-1])
            dataset.attrs["units"] = summary["units"] + (" squared" if "covariance" in name else "")
            divisor = areas**2 if "covariance" in name else areas
            dataset = handle.create_dataset(f"current_density/{name}", data=(values / divisor)[..., ::-1])
            dataset.attrs["units"] = (
                summary["current_density_units"] + (" squared" if "covariance" in name else "")
            )
    summary["output_h5"] = {"path": str(output_h5), "sha256": _sha256(output_h5)}
    if summary_path is not None:
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        with summary_path.open("x") as stream:
            stream.write(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary
