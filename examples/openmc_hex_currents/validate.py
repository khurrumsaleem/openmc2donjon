"""Native hex-face current checks: two adjacent hexagons and two axial layers.

All current faces are interior transmission interfaces. A surrounding void box
has a separate vacuum boundary. No collaborator geometry or nuclear data is
included. Run with a matching OpenMC Python module/executable and an installed
openmc2donjon checkout. Outputs are diagnostic currents, not MGXS/SPH objects.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import time

import h5py
import numpy as np
import openmc

from openmc2donjon.openmc_hex_currents import build_hex_current_surfaces, bind_hex_current_domains


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def face_center(node, face):
    cx, cy = node["center_xy_cm"]
    z0, z1 = node["z_bounds_cm"]
    normal = np.asarray(face["outward_normal"], dtype=float)
    center = np.array([cx, cy, (z0 + z1) / 2])
    return center + normal * (1.0 if face["index"] < 6 else (z1 - z0) / 2)


def analytic_events(nodes, origin, direction):
    """Independent ray/convex-polyhedron calculation, not OpenMC geometry.find."""
    events = []
    for n, node in enumerate(nodes):
        entry, exit_ = -math.inf, math.inf
        entry_face = exit_face = None
        missed = False
        for f, face in enumerate(node["faces"]):
            normal = np.asarray(face["outward_normal"])
            numerator = float(np.dot(normal, face_center(node, face) - origin))
            denominator = float(np.dot(normal, direction))
            if abs(denominator) < 1e-12:
                if numerator < -1e-10:
                    missed = True
                    break
                continue
            distance = numerator / denominator
            if denominator > 0 and distance < exit_:
                exit_, exit_face = distance, f
            elif denominator < 0 and distance > entry:
                entry, entry_face = distance, f
        if missed or entry > exit_ + 1e-10 or exit_ < 0:
            continue
        if entry > 1e-8:
            events.append([n, entry_face, "incoming"])
        if exit_ > 1e-8:
            events.append([n, exit_face, "outgoing"])
    return events


def make_case(case, mode, *, batches=12, particles=8000):
    case.mkdir(parents=True, exist_ok=False)
    registry = build_hex_current_surfaces([(0, 0), (0, 1)], [0, 1, 2], 2.0)
    # Shared nested fine universe tests ancestor cell filters, not just leaf cells.
    split = openmc.XPlane(x0=0.123)
    fine = openmc.Universe(cells=[openmc.Cell(region=-split), openmc.Cell(region=+split)])
    domains = {
        (i, j, k): openmc.Cell(name=f"HEX_{i}_{j}_Z{k}", region=registry.node_region((i, j), k), fill=fine)
        for i, j in registry.sites
        for k in range(2)
    }
    box = openmc.model.RectangularParallelepiped(-10, 10, -10, 10, -5, 5, boundary_type="vacuum")
    union = openmc.Union([cell.region for cell in domains.values()])
    remainder = openmc.Cell(name="void_outside_current_domains", region=-box & ~union)
    geometry = openmc.Geometry(openmc.Universe(cells=[*domains.values(), remainder]))
    # 64 distinct energies identify all 4 nodes × 8 faces × 2 crossing directions.
    edges = np.arange(0.5, 65.0, 1.0) * 1e5 if mode == "directed" else np.array([0.0, 1e6, 2e6])
    plan = bind_hex_current_domains(registry, domains, energy_bounds=edges, geometry=geometry)
    nodes = plan.manifest["nodes"]
    rays = []
    if mode == "directed":
        for n, node in enumerate(nodes):
            for f, face in enumerate(node["faces"]):
                normal = np.asarray(face["outward_normal"], dtype=float)
                point = face_center(node, face)
                for sense in (1, -1):
                    origin = point - sense * normal * 1e-4
                    direction = sense * normal
                    rays.append(
                        {
                            "origin": origin.tolist(),
                            "direction": direction.tolist(),
                            "energy_eV": (len(rays) + 1) * 1e5,
                            "target": [n, f, "outgoing" if sense == 1 else "incoming"],
                            "events": analytic_events(nodes, origin, direction),
                        }
                    )
    else:
        for origin, direction in (((0, 0, 0.5), (0, 1, 0)), ((0, 2, 0.5), (0, -1, 0))):
            rays.append(
                {
                    "origin": list(origin),
                    "direction": list(direction),
                    "energy_eV": 5e5,
                    "events": analytic_events(nodes, np.array(origin), np.array(direction)),
                }
            )
    sources = [
        openmc.IndependentSource(
            space=openmc.stats.Point(ray["origin"]),
            angle=openmc.stats.Monodirectional(ray["direction"]),
            energy=openmc.stats.Discrete([ray["energy_eV"]], [1.0]),
            strength=1 / len(rays),
        )
        for ray in rays
    ]
    settings = openmc.Settings()
    settings.run_mode = "fixed source"
    settings.particles, settings.batches = particles, batches
    settings.seed = 20260916
    settings.source = sources
    settings.statepoint = {"batches": list(range(1, batches + 1))}
    settings.output = {"tallies": False}
    # Unit-strength source and standard reduced per-batch tallies are required
    # for the paired-batch covariance validation below.
    # Default reduced tallies; this runner never passes OpenMC's --no-reduce.
    geometry.export_to_xml(case / "geometry.xml")
    openmc.Materials().export_to_xml(case / "materials.xml")
    settings.export_to_xml(case / "settings.xml")
    plan.tallies.export_to_xml(case / "tallies.xml")
    manifest = plan.manifest
    manifest["input_files"] = {
        role: {"path": f"{role}.xml", "sha256": sha(case / f"{role}.xml")}
        for role in ("geometry", "materials", "settings", "tallies")
    }
    manifest["validation_case"] = mode
    manifest["source_strength_sum"] = sum(s.strength for s in sources)
    manifest["production_physics_accepted"] = False
    dump(case / "current_manifest.json", manifest)
    dump(case / "analytic_rays.json", rays)
    return manifest, rays


def run_case(case, exe, *, threads=2, timeout=180):
    log = case / "openmc.log"
    if log.exists() or list(case.glob("statepoint.*.h5")):
        raise ValueError("Refusing to overwrite a previous run")
    manifest = json.loads((case / "current_manifest.json").read_text())
    command = [exe, "-s", str(threads)]
    started = time.time()
    env = dict(os.environ, OMPI_MCA_btl="self", OMPI_MCA_pml="ob1")
    with log.open("w") as stream:
        try:
            result = subprocess.run(
                command, cwd=case, env=env, stdout=stream, stderr=subprocess.STDOUT, timeout=timeout
            )
            code, timed_out = result.returncode, False
        except subprocess.TimeoutExpired:
            code, timed_out = 124, True
    outputs = {
        path.name: sha(path)
        for path in [*case.glob("statepoint.*.h5"), case / "summary.h5"]
        if path.is_file()
    }
    receipt = {
        "command": command,
        "threads": threads,
        "mpi_ranks": 1,
        "returncode": code,
        "timed_out": timed_out,
        "elapsed_seconds": time.time() - started,
        "input_files": manifest["input_files"],
        "output_files": outputs,
        "manifest_definition_sha256": hashlib.sha256(
            json.dumps(
                {key: value for key, value in manifest.items() if key != "execution_receipt"},
                sort_keys=True,
                separators=(",", ":"),
                allow_nan=False,
            ).encode()
        ).hexdigest(),
    }
    dump(case / "execution_receipt.json", receipt)
    manifest["execution_receipt"] = {
        "path": "execution_receipt.json",
        "sha256": sha(case / "execution_receipt.json"),
    }
    dump(case / "current_manifest.json", manifest)
    if code:
        raise RuntimeError(f"OpenMC failed ({code}); inspect {log}")


def validate_export(case, manifest, rays, mode):
    from openmc2donjon.openmc_hex_current_export import export_hex_currents

    statepoints = sorted(case.glob("statepoint.*.h5"), key=lambda p: int(p.name.split(".")[1]))
    output = case / "hex_currents.h5"
    export_hex_currents(
        statepoints, case / "current_manifest.json", output, summary_json=case / "export_summary.json"
    )
    with h5py.File(output, "r") as h5:
        out = h5["currents/out_mean"][:]
        inc = h5["currents/in_mean"][:]
        out_std = h5["currents/out_std_dev"][:]
        in_std = h5["currents/in_std_dev"][:]
        net = h5["currents/net_mean"][:]
        net_std = h5["currents/net_std_dev"][:]
    checks = {
        "net_mean_equals_out_minus_in": bool(np.allclose(net, out - inc, atol=1e-12, rtol=1e-12)),
        "all_partial_means_nonnegative": bool(np.all(out >= 0) and np.all(inc >= 0)),
    }
    with openmc.StatePoint(statepoints[-1]) as statepoint:
        leakage = next(
            row for row in statepoint.global_tallies
            if (row["name"].decode() if isinstance(row["name"], bytes) else row["name"]) == "leakage"
        )
        checks["void_global_leakage_equals_source"] = bool(np.isclose(leakage["mean"], 1.0, atol=1e-12))
    nodes = manifest["nodes"]
    index_by_cell = {node["cell_id"]: n for n, node in enumerate(nodes)}
    for n, node in enumerate(nodes):
        for f, face in enumerate(node["faces"]):
            if face["neighbor_cell_id"] is not None:
                m, opposite = index_by_cell[face["neighbor_cell_id"]], face["neighbor_face_index"]
                checks[f"shared_face_{n}_{f}_out_in"] = bool(
                    np.allclose(out[n, f], inc[m, opposite], atol=1e-12, rtol=1e-12)
                )
                checks[f"shared_face_{n}_{f}_net_opposite"] = bool(
                    np.allclose(net[n, f], -net[m, opposite], atol=1e-12, rtol=1e-12)
                )
    if mode == "directed":
        anchors = []
        for beam, ray in enumerate(rays):
            # HDF5 exports high-to-low group order, while ray energies increase.
            group = len(rays) - 1 - beam
            expected = {tuple(event) for event in ray["events"]}
            anchor_node, anchor_face, anchor_direction = ray["target"]
            anchor = (
                out[anchor_node, anchor_face, group]
                if anchor_direction == "outgoing"
                else inc[anchor_node, anchor_face, group]
            )
            anchors.append(float(anchor))
            checks[f"beam_{beam}_observed"] = bool(anchor > 0)
            for direction, values in (("outgoing", out), ("incoming", inc)):
                target = np.zeros((len(nodes), 8))
                for n, f, d in expected:
                    if d == direction:
                        target[n, f] = anchor
                checks[f"beam_{beam}_{direction}_exact_pattern"] = bool(
                    np.allclose(values[:, :, group], target, atol=1e-12, rtol=1e-12)
                )
        checks["directed_source_weights_sum_one"] = bool(np.isclose(sum(anchors), 1.0, atol=1e-12))
    else:
        n = next(n for n, node in enumerate(nodes) if node["node"] == [0, 0, 0])
        f, group = 0, 1  # north shared face, lower-energy bin in high-to-low order
        checks["opposing_beams_out_plus_in_one"] = bool(
            np.isclose(out[n, f, group] + inc[n, f, group], 1, atol=1e-12)
        )
        paired = net_std[n, f, group]
        independent = math.hypot(out_std[n, f, group], in_std[n, f, group])
        checks["net_uncertainty_preserves_anticorrelation"] = bool(
            paired > 0 and np.isclose(paired / independent, math.sqrt(2), rtol=1e-8)
        )
        checks["upper_layer_not_scored"] = bool(
            all(
                np.all(out[n] == 0) and np.all(inc[n] == 0)
                for n, node in enumerate(nodes)
                if node["node"][2] == 1
            )
        )
        checks["void_nodal_balance_matches_source"] = bool(np.isclose(np.sum(net), 1.0, atol=1e-12))
    report = {
        "case": mode,
        "checks": checks,
        "passed": all(checks.values()),
        "scope": "Finite-face binning, direction and batch covariance; not reactor physics acceptance",
    }
    dump(case / "validation.json", report)
    if not report["passed"]:
        raise AssertionError([name for name, passed in checks.items() if not passed])
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--openmc", default="openmc")
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--build-only", action="store_true")
    args = parser.parse_args()
    root = args.work_dir.resolve()
    if root.exists() and any(root.iterdir()):
        parser.error("Choose a new output directory")
    if args.threads <= 0:
        parser.error("Threads must be positive")
    executable = shutil.which(args.openmc)
    if not args.build_only and executable is None:
        parser.error("OpenMC executable not found")
    root.mkdir(parents=True, exist_ok=True)
    reports = []
    for mode in ("directed", "opposing"):
        case = root / mode
        manifest, rays = make_case(case, mode)
        print(f"Built {mode}: {case}", flush=True)
        if not args.build_only:
            run_case(case, executable, threads=args.threads)
            reports.append(validate_export(case, manifest, rays, mode))
            print(f"Passed {mode}: {len(reports[-1]['checks'])} checks", flush=True)
    if not args.build_only:
        dump(
            root / "validation_summary.json",
            {
                "passed": all(r["passed"] for r in reports),
                "cases": reports,
                "production_physics_accepted": False,
            },
        )


if __name__ == "__main__":
    main()
