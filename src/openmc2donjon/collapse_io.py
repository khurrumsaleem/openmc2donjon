"""Failure-safe publication and source binding for pre-equivalence collapse."""

from __future__ import annotations

from contextlib import ExitStack
from dataclasses import dataclass
import json
import os
from pathlib import Path
import shutil
import tempfile
from typing import Callable

import h5py
import numpy as np

from .hdf5_names import read_mixture_names
from .mgxs_fields import COLLAPSE_VECTOR_FIELDS
from .mgxs_input_contract import validate_openmc_volume_flux
from .mgxs_input_report import InputReport
from .mgxs_input_uncertainty import UncertaintyConfig
from .collapse_substitution import read_substitution_records
from .openmc_provenance import (
    OPENMC_PROVENANCE_GROUP,
    file_sha256,
    provenance_before_hdf5_mutation,
    refresh_openmc_provenance_after_hdf5_mutation,
)


HISTORY_PATH = "provenance/collapse/history_json"
HISTORY_SCHEMA = "openmc2donjon.collapse-history.v1"


def run_collapse(
    input_h5: str | Path,
    output_h5: str | Path,
    *,
    writer: Callable[[Path, Path], dict[str, object]],
    force: bool,
    summary_json: str | Path | None,
) -> dict[str, object]:
    """Build and bind a sibling file before atomically publishing the HDF5.

    The input is never overwritten. An existing forced destination survives
    calculation, write, and provenance-refresh failures byte for byte.
    """

    source, destination = Path(input_h5), Path(output_h5)
    summary = None if summary_json is None else Path(summary_json)
    if not source.is_file():
        raise FileNotFoundError(f"input HDF5 does not exist: {source}")
    if _same_file(source, destination):
        raise ValueError("collapse output must be different from input HDF5")
    if destination.exists() and (not force or not destination.is_file()):
        raise FileExistsError(f"output already exists; use --force for a file: {destination}")
    if summary is not None:
        if any(_same_file(summary, path) for path in (source, destination)):
            raise ValueError("summary JSON must be different from input/output HDF5")
        if summary.exists() and not summary.is_file():
            raise ValueError("summary JSON must name a file, not a directory")

    parent_sha256 = file_sha256(source)
    record = provenance_before_hdf5_mutation(source)
    with h5py.File(source, "r") as h5:
        if record is None and (
            OPENMC_PROVENANCE_GROUP in h5
            or any(str(key).startswith("openmc_provenance_") for key in h5.attrs)
        ):
            raise ValueError("cannot collapse malformed embedded OpenMC provenance")
        _require_pre_equivalence(h5)
        _validate_collapse_metadata(h5)
        history = _read_history(h5)

    destination.parent.mkdir(parents=True, exist_ok=True)
    with ExitStack() as stack:
        staging = Path(stack.enter_context(tempfile.TemporaryDirectory(
            prefix=".openmc2donjon-collapse-", dir=destination.parent,
        )))
        staged_h5 = staging / destination.name
        report = writer(source, staged_h5)
        report["output_h5"] = str(destination)
        # Do not bind an output to a source that changed during processing.
        if file_sha256(source) != parent_sha256:
            raise ValueError("collapse input changed during processing; output not published")
        history["steps"].append({
            "operation": report["schema"],
            "input_h5": str(source.resolve()),
            "input_h5_sha256": parent_sha256,
            "mapping": report["groups"],
            "weight": report["weight"],
            "parameters": {
                key: report[key] for key in (
                    "energy_group_structure", "source_energy_groups", "energy_groups",
                    "source_mixture_count", "component_count",
                ) if key in report
            },
        })
        with h5py.File(staged_h5, "r+") as h5:
            _validate_collapse_metadata(h5)
            if HISTORY_PATH in h5:
                del h5[HISTORY_PATH]
            group = h5.require_group("provenance/collapse")
            group.create_dataset(
                "history_json", data=json.dumps(history, sort_keys=True),
                dtype=h5py.string_dtype(encoding="utf-8"),
            )
        # The history is inside the payload covered by this refreshed binding.
        # Reuse only the verified input record; never repair damaged provenance.
        refresh_openmc_provenance_after_hdf5_mutation(staged_h5, record)
        report["input_h5_sha256"] = parent_sha256
        report["output_h5_sha256"] = file_sha256(staged_h5)
        report["collapse_history_steps"] = len(history["steps"])

        staged_summary = None
        if summary is not None:
            summary.parent.mkdir(parents=True, exist_ok=True)
            summary_staging = Path(stack.enter_context(tempfile.TemporaryDirectory(
                prefix=".openmc2donjon-collapse-summary-", dir=summary.parent,
            )))
            staged_summary = summary_staging / "new.json"
            staged_summary.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        summary_backup = None
        if staged_summary is not None:
            if summary.exists() or summary.is_symlink():
                summary_backup = staged_summary.with_name("previous.json")
                shutil.copy2(summary, summary_backup, follow_symlinks=False)
            os.replace(staged_summary, summary)
        try:
            if force:
                os.replace(staged_h5, destination)
            else:
                # Sibling staging is on the same filesystem. link() also protects
                # a destination created after the initial check.
                os.link(staged_h5, destination)
        except BaseException:
            # Publish the HDF5 last. If publication fails, restore the old
            # summary as well, so it cannot announce an output never delivered.
            if staged_summary is not None:
                if summary_backup is None:
                    summary.unlink()
                else:
                    os.replace(summary_backup, summary)
            raise
    return report


def _same_file(left: Path, right: Path) -> bool:
    return left.resolve() == right.resolve() or (
        left.exists() and right.exists() and left.samefile(right)
    )


def _require_pre_equivalence(h5) -> None:
    def check(group):
        datasets = {"sph", "SPH", "NSPH", "adf", "ADF", "discontinuity_factors"}.intersection(group)
        metadata = []
        for key, value in group.attrs.items():
            if key == "sph_applied":
                text = value.decode() if isinstance(value, bytes) else str(value)
                if text.strip().lower() in {"false", "0"}:
                    continue
            if str(key).startswith("sph_"):
                metadata.append(key)
        if datasets or metadata:
            raise ValueError(
                f"{group.name}: collapse requires pre-equivalence data; SPH/ADF data "
                "or applied SPH cannot be carried across a changed model. Collapse "
                "the uncorrected reference first, then recompute equivalence; "
                "for an unchanged model, use the original HDF5 directly."
            )
        for name in group:
            if str(name).endswith("_std_dev"):
                mean_name = str(name)[:-len("_std_dev")]
                if mean_name not in group or not isinstance(group[mean_name], h5py.Dataset):
                    raise ValueError(f"{group.name}/{name} has no matching mean dataset")
                read_standard_deviation(group, mean_name, group[mean_name].shape)

    check(h5)
    if "mixtures" in h5:
        check(h5["mixtures"])
        h5["mixtures"].visititems(lambda _name, obj: check(obj) if isinstance(obj, h5py.Group) else None)


def _read_history(h5) -> dict:
    if HISTORY_PATH not in h5:
        return {"schema": HISTORY_SCHEMA, "steps": []}
    try:
        result = json.loads(h5[HISTORY_PATH].asstr()[()])
    except (TypeError, ValueError, AttributeError) as exc:
        raise ValueError("invalid collapse processing history") from exc
    if not isinstance(result, dict) or result.get("schema") != HISTORY_SCHEMA:
        raise ValueError("invalid collapse processing history schema")
    if not isinstance(result.get("steps"), list) or not all(
        isinstance(step, dict) for step in result["steps"]
    ):
        raise ValueError("invalid collapse processing history steps")
    return result


def read_standard_deviation(group, mean_name: str, shape: tuple[int, ...]) -> np.ndarray | None:
    name = f"{mean_name}_std_dev"
    if name not in group:
        return None
    values = np.asarray(group[name][()], dtype=float)
    if values.shape != shape or not np.all(np.isfinite(values)) or np.any(values < 0.0):
        raise ValueError(f"{group.name}/{name} must match {shape} and contain finite non-negative values")
    return values


def _validate_collapse_metadata(h5) -> None:
    raw = h5.attrs.get("energy_groups")
    if isinstance(raw, (bool, np.bool_)) or not isinstance(raw, (int, np.integer)) or raw <= 0:
        raise ValueError("collapse requires a positive integer energy_groups declaration")
    ngroups = int(raw)
    names = read_mixture_names(h5)
    report = InputReport(path=str(h5.filename))
    validate_openmc_volume_flux(
        h5, names, ngroups, report, require_openmc_volume_flux=True,
        uncertainty=UncertaintyConfig(warn_threshold=None),
    )
    if not report.ok:
        raise ValueError("collapse reference flux contract: " + "; ".join(report.issues))
    for name in names:
        read_substitution_records(h5["mixtures"][name], ngroups)


@dataclass(frozen=True)
class CollapseVector:
    values: np.ndarray
    std: np.ndarray | None
    source_name: str
    std_source_name: str | None


def read_collapse_vector(group, canonical: str, ngroups: int) -> CollapseVector | None:
    """Resolve Converter aliases without discarding conflicting data/errors."""
    names = [name for name in COLLAPSE_VECTOR_FIELDS[canonical] if name in group]
    if not names:
        return None
    arrays = []
    deviations = []
    for name in names:
        if not isinstance(group[name], h5py.Dataset):
            raise ValueError(f"{group.name}/{name} must be a vector dataset")
        values = np.asarray(group[name][()], dtype=float)
        if values.shape != (ngroups,) or not np.all(np.isfinite(values)):
            raise ValueError(f"{group.name}/{name} must be a finite group-wise vector")
        arrays.append(values)
        std = read_standard_deviation(group, name, values.shape)
        if std is not None:
            deviations.append((name, std))
    if any(not np.array_equal(arrays[0], values) for values in arrays[1:]):
        raise ValueError(f"{group.name}: conflicting aliases for {canonical}: {names}")
    if any(not np.array_equal(deviations[0][1], std) for _, std in deviations[1:]):
        raise ValueError(f"{group.name}: conflicting alias uncertainties for {canonical}")
    return CollapseVector(
        arrays[0], deviations[0][1] if deviations else None, names[0],
        deviations[0][0] if deviations else None,
    )
