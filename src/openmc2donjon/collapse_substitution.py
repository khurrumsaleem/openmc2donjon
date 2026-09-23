"""Keep zero-flux donor provenance through spatial and energy collapse.

Indices are zero-based DONJON groups. A derived index means that the output
contains substituted contributions, not that its entire XS was donor-filled.
"""

from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np


RECORD_DATASET = "zero_flux_substitution_records"
RECORD_SCHEMA = "openmc2donjon.collapsed-substitution.v1"
SEMANTICS = "contains-substituted-contributions"


def _indices(value, ngroups: int, label: str) -> list[int]:
    values = np.asarray(value)
    if values.ndim != 1 or (values.size and values.dtype.kind not in "iu"):
        raise ValueError(f"{label}: zero-flux group indices must be one-dimensional integers")
    result = [int(index) for index in values]
    if len(set(result)) != len(result) or any(index < 0 or index >= ngroups for index in result):
        raise ValueError(f"{label}: zero-flux group indices must be unique and in [0, {ngroups})")
    return sorted(result)


def _json_value(value):
    if isinstance(value, np.ndarray):
        return _json_value(value.tolist())
    if isinstance(value, np.generic):
        return _json_value(value.item())
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    return value


def read_substitution_records(group, ngroups: int) -> list[dict]:
    metadata = {key: _json_value(value) for key, value in group.attrs.items()
                if key.startswith("zero_flux_")}
    if not metadata and RECORD_DATASET not in group:
        return []
    if "zero_flux_filled_groups" not in metadata:
        raise ValueError(f"{group.name}: zero-flux metadata requires zero_flux_filled_groups")
    affected = _indices(metadata["zero_flux_filled_groups"], ngroups, group.name)
    if RECORD_DATASET not in group:
        if "zero_flux_fill_semantics" in metadata:
            raise ValueError(f"{group.name}: derived zero-flux metadata is missing its origin records")
        # Preserve all donor metadata, including unknown future fields, without
        # asserting missing donor properties. The source payload is bound by
        # the transformation history's file SHA256.
        try:
            json.dumps(metadata, allow_nan=False)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{group.name}: invalid zero-flux donor metadata") from exc
        return [{
            "source_h5": str(Path(group.file.filename).resolve()),
            "source_mixture": group.name,
            "source_energy_groups": ngroups,
            "filled_groups": affected,
            "donor_metadata": metadata,
            "current_groups": affected,
        }] if affected else []

    if metadata.get("zero_flux_fill_semantics") != SEMANTICS:
        raise ValueError(f"{group.name}: derived zero-flux records require contribution semantics")
    try:
        record = json.loads(group[RECORD_DATASET].asstr()[()])
        if not isinstance(record, dict) or record.get("schema") != RECORD_SCHEMA:
            raise ValueError("unknown schema")
        if record.get("group_order") != "mgxs_donjon" or record.get("index_base") != 0:
            raise ValueError("invalid group order/index base")
        origins = record["origins"]
        if not isinstance(origins, list) or not origins:
            raise ValueError("missing origin records")
        covered = set()
        for origin in origins:
            if not isinstance(origin, dict):
                raise ValueError("invalid origin")
            for key in ("source_h5", "source_mixture"):
                if not isinstance(origin[key], str) or not origin[key]:
                    raise ValueError(f"missing {key}")
            count = origin["source_energy_groups"]
            if type(count) is not int or count <= 0:
                raise ValueError("invalid origin group count")
            filled = _indices(origin["filled_groups"], count, group.name)
            current = _indices(origin["current_groups"], ngroups, group.name)
            donor = origin["donor_metadata"]
            if not filled or not current or not isinstance(donor, dict):
                raise ValueError("invalid origin contribution")
            if _indices(donor["zero_flux_filled_groups"], count, group.name) != filled:
                raise ValueError("origin donor groups disagree")
            json.dumps(donor, allow_nan=False)
            covered.update(current)
        if sorted(covered) != affected:
            raise ValueError("affected groups disagree with origin records")
    except (KeyError, TypeError, ValueError, AttributeError) as exc:
        raise ValueError(f"{group.name}: invalid derived zero-flux substitution records: {exc}") from exc
    return origins


def write_collapsed_substitutions(target, sources, source_groups: int, mapping) -> None:
    """Map source group indices; keep all contributing mixtures' origins."""
    fine_to_coarse = {int(fine): coarse for coarse, members in enumerate(mapping) for fine in members}
    origins = []
    for source in sources:
        for origin in read_substitution_records(source, source_groups):
            origins.append({
                **origin,
                "current_groups": sorted({fine_to_coarse[index] for index in origin["current_groups"]}),
            })
    if not origins:
        return
    target.attrs["zero_flux_filled_groups"] = sorted({
        index for origin in origins for index in origin["current_groups"]
    })
    target.attrs["zero_flux_fill_semantics"] = SEMANTICS
    # A scalar string dataset avoids HDF5's small-attribute size limit for
    # components combining hundreds of source regions.
    target.create_dataset(RECORD_DATASET, data=json.dumps({
        "schema": RECORD_SCHEMA, "group_order": "mgxs_donjon", "index_base": 0,
        "origins": origins,
    }, sort_keys=True, allow_nan=False), dtype=h5py.string_dtype(encoding="utf-8"))
