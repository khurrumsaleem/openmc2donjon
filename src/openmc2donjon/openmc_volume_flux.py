"""Helpers for writing OpenMC volume-flux reference maps."""

from __future__ import annotations

from dataclasses import dataclass
import json
import operator
from os import PathLike
from pathlib import Path
from typing import Iterable

import numpy as np

from . import __version__
from .constants import MGXS_DONJON_GROUP_ORDER
from .hdf5_names import read_mixture_names


DATASET_NAME = "openmc_volume_flux"
STD_DEV_DATASET_NAME = "openmc_volume_flux_std_dev"
SCHEMA = "openmc2donjon.openmc-volume-flux.v1"
PASS_DECISION = "openmc2donjon_volume_flux_export_passed"
DEFAULT_SOURCE_GROUP_ORDER = "openmc_energy_filter_reversed"


@dataclass(frozen=True)
class OpenMCVolumeFluxReport:
    output_h5: Path
    dataset: str
    mixture_names: tuple[str, ...]
    energy_groups: int
    source_group_order: str
    minimum: float
    maximum: float
    std_dev_dataset: str | None = None
    max_relative_std_dev: float | None = None
    statepoint: Path | None = None
    tally_name: str | None = None
    allow_zero_flux: bool = False
    energy_bounds_verified: bool = False
    spatial_domain_order_verified: bool = False
    source_filter_order: tuple[str, ...] = ()
    source_domain_ids: tuple[int, ...] = ()


@dataclass(frozen=True)
class _VolumeFluxMetadata:
    mixture_names: tuple[str, ...]
    energy_groups: int
    energy_bounds: np.ndarray | None = None
    source_domain_ids: tuple[int, ...] | None = None


def export_openmc_volume_flux(
    statepoint: str | PathLike[str],
    output_h5: str | PathLike[str],
    *,
    mgxs_h5: str | PathLike[str] | None = None,
    tally_name: str = DATASET_NAME,
    dataset_name: str = DATASET_NAME,
    std_dev_dataset_name: str | None = None,
    mixture_names: Iterable[str] | None = None,
    energy_groups: int | None = None,
    source_domain_ids: Iterable[int] | None = None,
    source_group_order: str = DEFAULT_SOURCE_GROUP_ORDER,
    allow_zero: bool = False,
    force: bool = False,
    summary_json: str | PathLike[str] | None = None,
) -> OpenMCVolumeFluxReport:
    """Export an OpenMC statepoint volume-flux tally as canonical HDF5.

    The written dataset has shape ``(mixture, group)`` and carries the
    ``group_order='mgxs_donjon'`` and ``mixture_names`` attributes required by
    ``make-openmc-sph-sidecar``.  OpenMC EnergyFilter tally bins are
    low-to-high; values are reversed along the energy axis by default.
    """

    try:
        import openmc
    except ImportError as exc:  # pragma: no cover - depends on user environment
        raise RuntimeError(
            "OpenMC is required to export a volume-flux tally from a statepoint"
        ) from exc

    statepoint_path = Path(statepoint)
    output_path = Path(output_h5)
    if not statepoint_path.exists():
        raise FileNotFoundError(f"statepoint does not exist: {statepoint_path}")
    if output_path.exists() and not force:
        raise FileExistsError(f"output already exists; use --force to overwrite: {output_path}")

    metadata = _resolve_metadata(
        mgxs_h5=None if mgxs_h5 is None else Path(mgxs_h5),
        mixture_names=tuple(mixture_names) if mixture_names is not None else None,
        energy_groups=energy_groups,
        source_domain_ids=(
            tuple(source_domain_ids) if source_domain_ids is not None else None
        ),
    )
    with openmc.StatePoint(str(statepoint_path)) as sp:
        tally = sp.get_tally(name=tally_name)
        mean = tally.get_values(scores=["flux"], value="mean")
        std_dev = tally.get_values(scores=["flux"], value="std_dev")
        filter_order = _tally_filter_order(tally)
        if mgxs_h5 is not None:
            flux = _validated_mgxs_tally_flux(
                tally,
                mean,
                metadata=metadata,
                value_name="mean",
            )
            flux_std_dev = _validated_mgxs_tally_flux(
                tally,
                std_dev,
                metadata=metadata,
                value_name="std_dev",
            )
        else:
            # Explicit, no-MGXS exports remain available as a diagnostic
            # compatibility path.  Without MGXS bounds and source-domain IDs,
            # their energy and spatial ordering cannot be certified.
            flux = reverse_openmc_energy_filter_flux(
                mean,
                mixture_count=len(metadata.mixture_names),
                energy_groups=metadata.energy_groups,
            )
            flux_std_dev = reverse_openmc_energy_filter_flux(
                std_dev,
                mixture_count=len(metadata.mixture_names),
                energy_groups=metadata.energy_groups,
            )

    report = _write_openmc_flux_hdf5(
        output_path,
        flux,
        mixture_names=metadata.mixture_names,
        std_dev=flux_std_dev,
        dataset_name=dataset_name,
        std_dev_dataset_name=(
            std_dev_dataset_name
            if std_dev_dataset_name is not None
            else f"{dataset_name}_std_dev"
        ),
        source_group_order=source_group_order,
        allow_zero=allow_zero,
        replace=True,
        statepoint=statepoint_path,
        tally_name=tally_name,
        energy_bounds_verified=mgxs_h5 is not None,
        spatial_domain_order_verified=mgxs_h5 is not None,
        source_filter_order=filter_order,
        source_domain_ids=metadata.source_domain_ids or (),
    )
    print_report(report)
    if summary_json is not None:
        write_summary(Path(summary_json), report)
    return report


def reverse_openmc_energy_filter_flux(
    values: np.ndarray | Iterable[float],
    *,
    mixture_count: int,
    energy_groups: int,
) -> np.ndarray:
    """Return a volume-flux tally in MGXS/DONJON group order."""

    if mixture_count <= 0:
        raise ValueError("mixture_count must be positive")
    if energy_groups <= 0:
        raise ValueError("energy_groups must be positive")
    array = np.asarray(values, dtype=float).squeeze()
    expected = (int(mixture_count), int(energy_groups))
    try:
        reshaped = array.reshape(expected)
    except ValueError as exc:
        raise ValueError(
            f"OpenMC volume-flux tally cannot be reshaped to {expected}; "
            f"got {array.shape}"
        ) from exc
    return reshaped[:, ::-1].copy()


def write_openmc_volume_flux_hdf5(
    output_h5: str | PathLike[str],
    values: np.ndarray | Iterable[Iterable[float]],
    *,
    mixture_names: Iterable[str],
    std_dev: np.ndarray | Iterable[Iterable[float]] | None = None,
    dataset_name: str = DATASET_NAME,
    std_dev_dataset_name: str = STD_DEV_DATASET_NAME,
    source_group_order: str = DEFAULT_SOURCE_GROUP_ORDER,
    allow_zero: bool = False,
    replace: bool = True,
) -> OpenMCVolumeFluxReport:
    """Append a canonical ``/openmc_volume_flux`` dataset to an MGXS HDF5 file."""

    return write_openmc_flux_hdf5(
        output_h5,
        values,
        mixture_names=mixture_names,
        std_dev=std_dev,
        dataset_name=dataset_name,
        std_dev_dataset_name=std_dev_dataset_name,
        source_group_order=source_group_order,
        allow_zero=allow_zero,
        replace=replace,
    )


def write_openmc_flux_hdf5(
    output_h5: str | PathLike[str],
    values: np.ndarray | Iterable[Iterable[float]],
    *,
    mixture_names: Iterable[str],
    std_dev: np.ndarray | Iterable[Iterable[float]] | None = None,
    dataset_name: str = DATASET_NAME,
    std_dev_dataset_name: str = STD_DEV_DATASET_NAME,
    source_group_order: str = DEFAULT_SOURCE_GROUP_ORDER,
    allow_zero: bool = False,
    replace: bool = True,
    statepoint: str | PathLike[str] | None = None,
    tally_name: str | None = None,
) -> OpenMCVolumeFluxReport:
    """Append an unverified, caller-supplied OpenMC flux matrix to HDF5.

    Direct/manual writers cannot prove the tally filter layout, so the two
    verification attributes are always false.  Statepoint exports backed by
    ``--mgxs`` use the private verified writer only after filter validation.
    """

    return _write_openmc_flux_hdf5(
        output_h5,
        values,
        mixture_names=mixture_names,
        std_dev=std_dev,
        dataset_name=dataset_name,
        std_dev_dataset_name=std_dev_dataset_name,
        source_group_order=source_group_order,
        allow_zero=allow_zero,
        replace=replace,
        statepoint=statepoint,
        tally_name=tally_name,
        energy_bounds_verified=False,
        spatial_domain_order_verified=False,
        source_filter_order=(),
        source_domain_ids=(),
    )


def _write_openmc_flux_hdf5(
    output_h5: str | PathLike[str],
    values: np.ndarray | Iterable[Iterable[float]],
    *,
    mixture_names: Iterable[str],
    std_dev: np.ndarray | Iterable[Iterable[float]] | None = None,
    dataset_name: str = DATASET_NAME,
    std_dev_dataset_name: str = STD_DEV_DATASET_NAME,
    source_group_order: str = DEFAULT_SOURCE_GROUP_ORDER,
    allow_zero: bool = False,
    replace: bool = True,
    statepoint: str | PathLike[str] | None = None,
    tally_name: str | None = None,
    energy_bounds_verified: bool = False,
    spatial_domain_order_verified: bool = False,
    source_filter_order: Iterable[str] = (),
    source_domain_ids: Iterable[int] = (),
) -> OpenMCVolumeFluxReport:
    """Append a canonical OpenMC flux matrix and its verification evidence.

    ``dataset_name`` lets the same writer produce both CE reference flux
    (typically ``openmc_volume_flux``) and MG macro flux
    (for example ``openmc_mg_flux``) for OpenMC-side SPH equivalence.
    """

    import h5py

    path = Path(output_h5)
    names = _as_mixture_names(mixture_names)
    filter_order = tuple(str(value) for value in source_filter_order)
    domain_ids = tuple(int(value) for value in source_domain_ids)
    if bool(energy_bounds_verified) != bool(spatial_domain_order_verified):
        raise ValueError("energy and spatial layout verification must succeed together")
    if energy_bounds_verified and (not filter_order or not domain_ids):
        raise ValueError(
            "verified flux datasets require source_filter_order and source_domain_ids"
        )
    if domain_ids and len(domain_ids) != len(names):
        raise ValueError(
            "source_domain_ids must match mixture_names: "
            f"{len(domain_ids)} != {len(names)}"
        )
    if len(set(domain_ids)) != len(domain_ids):
        raise ValueError("source_domain_ids must be unique")
    flux = _as_flux_array(values, names=names, dataset_name=dataset_name, allow_zero=allow_zero)
    flux_std_dev = _as_std_dev_array(
        std_dev,
        expected_shape=flux.shape,
        dataset_name=std_dev_dataset_name,
    )

    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "a") as h5:
        if dataset_name in h5:
            if not replace:
                raise FileExistsError(f"{path}: /{dataset_name} already exists")
            del h5[dataset_name]
        if std_dev_dataset_name in h5 and replace:
            del h5[std_dev_dataset_name]
        if std_dev_dataset_name in h5 and flux_std_dev is not None:
            raise FileExistsError(f"{path}: /{std_dev_dataset_name} already exists")
        _write_flux_dataset(
            h5,
            dataset_name,
            flux,
            names=names,
            source_group_order=source_group_order,
            energy_bounds_verified=energy_bounds_verified,
            spatial_domain_order_verified=spatial_domain_order_verified,
            source_filter_order=filter_order,
            source_domain_ids=domain_ids,
        )
        if flux_std_dev is not None:
            _write_flux_dataset(
                h5,
                std_dev_dataset_name,
                flux_std_dev,
                names=names,
                source_group_order=source_group_order,
                energy_bounds_verified=energy_bounds_verified,
                spatial_domain_order_verified=spatial_domain_order_verified,
                source_filter_order=filter_order,
                source_domain_ids=domain_ids,
            )

    return OpenMCVolumeFluxReport(
        output_h5=path,
        dataset=dataset_name,
        mixture_names=names,
        energy_groups=int(flux.shape[1]),
        source_group_order=str(source_group_order),
        minimum=float(np.min(flux)),
        maximum=float(np.max(flux)),
        std_dev_dataset=std_dev_dataset_name if flux_std_dev is not None else None,
        max_relative_std_dev=(
            None
            if flux_std_dev is None
            else _max_relative_std_dev(flux, flux_std_dev)
        ),
        statepoint=None if statepoint is None else Path(statepoint),
        tally_name=tally_name,
        allow_zero_flux=allow_zero,
        energy_bounds_verified=bool(energy_bounds_verified),
        spatial_domain_order_verified=bool(spatial_domain_order_verified),
        source_filter_order=filter_order,
        source_domain_ids=domain_ids,
    )


def print_report(report: OpenMCVolumeFluxReport) -> None:
    print("OpenMC-to-DONJON volume flux export")
    print(f"  schema: {SCHEMA}")
    if report.statepoint is not None:
        print(f"  statepoint: {report.statepoint}")
    if report.tally_name is not None:
        print(f"  tally: {report.tally_name}")
    print(f"  output: {report.output_h5}::{report.dataset}")
    print(
        f"  mixtures={len(report.mixture_names)} groups={report.energy_groups} "
        f"group_order={MGXS_DONJON_GROUP_ORDER} "
        f"source_group_order={report.source_group_order}"
    )
    print(
        "  verification: "
        f"energy_bounds={str(report.energy_bounds_verified).lower()} "
        f"spatial_domain_order={str(report.spatial_domain_order_verified).lower()} "
        "source_filter_order="
        f"{','.join(report.source_filter_order) if report.source_filter_order else 'unavailable'}"
    )
    if report.source_domain_ids:
        print(
            "  canonical domains: "
            + ", ".join(
                f"{name}={domain_id}"
                for name, domain_id in zip(
                    report.mixture_names,
                    report.source_domain_ids,
                    strict=True,
                )
            )
        )
    print(
        "  flux range: "
        f"min={report.minimum:.6g} max={report.maximum:.6g}"
    )
    if report.std_dev_dataset is not None:
        print(
            f"  std_dev: {report.std_dev_dataset} "
            f"max_rel={report.max_relative_std_dev:.6g}"
        )
    print()
    print("Volume flux export decision")
    print(f"  {PASS_DECISION}")


def write_summary(path: Path, report: OpenMCVolumeFluxReport) -> None:
    payload = {
        "schema": SCHEMA,
        "package_version": __version__,
        "decision": PASS_DECISION,
        "statepoint": None if report.statepoint is None else str(report.statepoint),
        "tally_name": report.tally_name,
        "output_h5": str(report.output_h5),
        "dataset": report.dataset,
        "std_dev_dataset": report.std_dev_dataset,
        "mixture_count": len(report.mixture_names),
        "mixture_names": list(report.mixture_names),
        "energy_groups": report.energy_groups,
        "group_order": MGXS_DONJON_GROUP_ORDER,
        "source_group_order": report.source_group_order,
        "min": report.minimum,
        "max": report.maximum,
        "max_relative_std_dev": report.max_relative_std_dev,
        "allow_zero_flux": report.allow_zero_flux,
        "energy_bounds_verified": report.energy_bounds_verified,
        "spatial_domain_order_verified": report.spatial_domain_order_verified,
        "source_filter_order": list(report.source_filter_order),
        "source_domain_ids": list(report.source_domain_ids),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _as_mixture_names(values: Iterable[str]) -> tuple[str, ...]:
    names = tuple(str(value) for value in values)
    if not names:
        raise ValueError("mixture_names must not be empty")
    if any(not name for name in names):
        raise ValueError("mixture_names must not contain empty names")
    if len(set(names)) != len(names):
        raise ValueError("mixture_names must be unique")
    return names


def _resolve_metadata(
    *,
    mgxs_h5: Path | None,
    mixture_names: tuple[str, ...] | None,
    energy_groups: int | None,
    source_domain_ids: tuple[int, ...] | None,
) -> _VolumeFluxMetadata:
    mgxs_names: tuple[str, ...] | None = None
    mgxs_groups: int | None = None
    mgxs_bounds: np.ndarray | None = None
    if mgxs_h5 is not None:
        import h5py

        with h5py.File(mgxs_h5, "r") as h5:
            if "mixtures" not in h5:
                raise ValueError(f"{mgxs_h5}: missing /mixtures group")
            mgxs_names = read_mixture_names(h5)
            if "energy_bounds" not in h5:
                raise ValueError(
                    f"{mgxs_h5}: missing /energy_bounds dataset required to verify "
                    "the OpenMC EnergyFilter"
                )
            mgxs_bounds = _as_energy_bounds(
                np.asarray(h5["energy_bounds"][:], dtype=float),
                source=f"{mgxs_h5}:/energy_bounds",
            )
            mgxs_groups = int(mgxs_bounds.size - 1)
            if "energy_groups" in h5.attrs:
                declared_groups = int(h5.attrs["energy_groups"])
                if declared_groups != mgxs_groups:
                    raise ValueError(
                        f"{mgxs_h5}: energy_groups={declared_groups} disagrees with "
                        f"/energy_bounds ({mgxs_groups} groups)"
                    )
            if source_domain_ids is None:
                source_domain_ids = _read_source_domain_ids(
                    h5,
                    mgxs_h5=mgxs_h5,
                    mixture_names=mgxs_names,
                )

    if mixture_names is not None:
        names = _as_mixture_names(mixture_names)
        if mgxs_names is not None and names != mgxs_names:
            raise ValueError(
                "--mixture-names disagrees with --mgxs mixture order: "
                f"{names!r} != {mgxs_names!r}"
            )
    else:
        names = mgxs_names
    if names is None:
        raise ValueError("mixture names must be supplied, either via --mgxs or --mixture-names")

    if energy_groups is not None:
        groups = int(energy_groups)
        if mgxs_groups is not None and groups != mgxs_groups:
            raise ValueError(
                "--energy-groups disagrees with --mgxs /energy_bounds: "
                f"{groups} != {mgxs_groups}"
            )
    else:
        groups = int(mgxs_groups or 0)
    if groups <= 0:
        raise ValueError("energy group count must be supplied, either via --mgxs or --energy-groups")
    if source_domain_ids is not None:
        source_domain_ids = _as_source_domain_ids(
            source_domain_ids,
            mixture_names=_as_mixture_names(names),
            option="--source-domain-ids",
        )
    return _VolumeFluxMetadata(
        mixture_names=_as_mixture_names(names),
        energy_groups=groups,
        energy_bounds=mgxs_bounds,
        source_domain_ids=source_domain_ids,
    )


def _read_source_domain_ids(
    h5,
    *,
    mgxs_h5: Path,
    mixture_names: tuple[str, ...],
) -> tuple[int, ...]:
    ids: list[int] = []
    for name in mixture_names:
        group = h5["mixtures"][name]
        if "source_domain_id" not in group.attrs:
            raise ValueError(
                f"{mgxs_h5}: mixture {name!r} has no scalar source_domain_id; "
                "spatial tally order cannot be verified"
            )
        raw_value = group.attrs["source_domain_id"]
        array = np.asarray(raw_value)
        if array.ndim != 0:
            raise ValueError(
                f"{mgxs_h5}: mixture {name!r} source_domain_id is not scalar; "
                "spatial tally order cannot be verified"
            )
        try:
            value = operator.index(array.item())
        except TypeError as exc:
            raise ValueError(
                f"{mgxs_h5}: mixture {name!r} source_domain_id is not an integer; "
                "spatial tally order cannot be verified"
            ) from exc
        ids.append(int(value))
    if len(set(ids)) != len(ids):
        raise ValueError(
            f"{mgxs_h5}: mixture source_domain_id values are not unique; "
            "spatial tally order cannot be verified"
        )
    return tuple(ids)


def _as_source_domain_ids(
    values: Iterable[int],
    *,
    mixture_names: tuple[str, ...],
    option: str,
) -> tuple[int, ...]:
    ids: list[int] = []
    for raw_value in values:
        try:
            value = operator.index(raw_value)
        except TypeError as exc:
            raise ValueError(f"{option} must contain integer domain IDs") from exc
        ids.append(int(value))
    if len(ids) != len(mixture_names):
        raise ValueError(
            f"{option} must provide one ID per canonical mixture name: "
            f"{len(ids)} != {len(mixture_names)}"
        )
    if len(set(ids)) != len(ids):
        raise ValueError(f"{option} domain IDs must be unique")
    return tuple(ids)


def _as_energy_bounds(values: np.ndarray, *, source: str) -> np.ndarray:
    bounds = np.asarray(values, dtype=float).reshape(-1)
    if bounds.size < 2:
        raise ValueError(f"{source} must contain at least two energy bounds")
    if not np.all(np.isfinite(bounds)):
        raise ValueError(f"{source} must contain only finite values")
    if not np.all(np.diff(bounds) > 0.0):
        raise ValueError(f"{source} must be strictly ascending")
    return bounds


def _tally_filter_order(tally) -> tuple[str, ...]:
    filters = getattr(tally, "filters", None)
    if filters is None:
        return ()
    return tuple(type(filter_).__name__ for filter_ in filters)


def _validated_mgxs_tally_flux(
    tally,
    values: np.ndarray,
    *,
    metadata: _VolumeFluxMetadata,
    value_name: str,
) -> np.ndarray:
    """Return a verified tally matrix in ``[mixture, DONJON group]`` order."""

    if metadata.energy_bounds is None or metadata.source_domain_ids is None:
        raise ValueError("internal error: MGXS verification metadata is incomplete")

    filters = getattr(tally, "filters", None)
    if filters is None:
        raise ValueError(
            "selected tally does not expose filters; MGXS spatial and energy "
            "ordering cannot be verified"
        )
    filters = tuple(filters)
    if len(filters) != 2:
        raise ValueError(
            "selected tally must have exactly one EnergyFilter and one scalar "
            f"spatial/domain filter; found {_describe_filters(filters)}"
        )

    energy_indices = [
        index for index, filter_ in enumerate(filters) if _is_energy_filter(filter_)
    ]
    spatial_indices = [
        index for index, filter_ in enumerate(filters) if _is_spatial_domain_filter(filter_)
    ]
    if len(energy_indices) != 1 or len(spatial_indices) != 1:
        raise ValueError(
            "selected tally must have exactly one EnergyFilter and one supported "
            "spatial/domain filter (CellFilter, CellFromFilter, MaterialFilter, "
            f"or UniverseFilter); found {_describe_filters(filters)}"
        )
    energy_index = energy_indices[0]
    spatial_index = spatial_indices[0]
    if energy_index == spatial_index:
        raise ValueError("selected tally energy and spatial filters are ambiguous")

    energy_filter = filters[energy_index]
    spatial_filter = filters[spatial_index]
    _verify_energy_filter_bounds(energy_filter, metadata.energy_bounds)
    spatial_ids = _scalar_filter_ids(spatial_filter)
    expected_ids = metadata.source_domain_ids
    if len(spatial_ids) != len(expected_ids):
        raise ValueError(
            f"{type(spatial_filter).__name__} has {len(spatial_ids)} bins but "
            f"--mgxs declares {len(expected_ids)} mixtures"
        )
    if spatial_ids != expected_ids:
        raise ValueError(
            f"{type(spatial_filter).__name__} bin order {spatial_ids!r} does not "
            f"match --mgxs source_domain_id order {expected_ids!r}"
        )

    filter_counts = tuple(_filter_bin_count(filter_) for filter_ in filters)
    expected_count = int(np.prod(filter_counts, dtype=int))
    array = np.asarray(values, dtype=float)
    if array.size != expected_count:
        raise ValueError(
            f"selected tally {value_name} has {array.size} values, but verified "
            f"filter layout {filter_counts!r} requires {expected_count}"
        )
    filter_ordered = array.reshape(filter_counts)
    if spatial_index == 0:
        spatial_energy = filter_ordered
    else:
        spatial_energy = np.transpose(filter_ordered, (1, 0))
    expected_shape = (len(metadata.mixture_names), metadata.energy_groups)
    if spatial_energy.shape != expected_shape:
        raise ValueError(
            f"verified tally layout produced {spatial_energy.shape}, expected "
            f"{expected_shape}"
        )
    # OpenMC EnergyFilter bins are low-to-high.  The handoff contract stores
    # cross sections and flux in high-to-low MGXS/DONJON group-index order.
    return spatial_energy[:, ::-1].copy()


def _is_energy_filter(filter_) -> bool:
    return type(filter_).__name__ == "EnergyFilter"


def _is_spatial_domain_filter(filter_) -> bool:
    return type(filter_).__name__ in {
        "CellFilter",
        "CellFromFilter",
        "MaterialFilter",
        "UniverseFilter",
    }


def _describe_filters(filters: tuple[object, ...]) -> str:
    if not filters:
        return "no filters"
    return ", ".join(type(filter_).__name__ for filter_ in filters)


def _filter_bin_count(filter_) -> int:
    bins = getattr(filter_, "bins", None)
    if bins is None:
        raise ValueError(f"{type(filter_).__name__} does not expose bins")
    try:
        count = len(bins)
    except TypeError as exc:
        raise ValueError(f"{type(filter_).__name__} bins are not inspectable") from exc
    declared = getattr(filter_, "num_bins", count)
    try:
        declared_count = int(declared)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"{type(filter_).__name__} num_bins is not an integer"
        ) from exc
    if declared_count != count:
        raise ValueError(
            f"{type(filter_).__name__} num_bins={declared_count} disagrees with "
            f"its {count} exposed bins"
        )
    if count <= 0:
        raise ValueError(f"{type(filter_).__name__} has no bins")
    return count


def _verify_energy_filter_bounds(energy_filter, expected_bounds: np.ndarray) -> None:
    bins = getattr(energy_filter, "bins", None)
    if bins is None:
        raise ValueError("EnergyFilter does not expose bins")
    try:
        pairs = np.asarray(bins, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("EnergyFilter bins are not numeric edge pairs") from exc
    if pairs.ndim == 1 and pairs.size == 2:
        pairs = pairs.reshape((1, 2))
    if pairs.ndim != 2 or pairs.shape[1] != 2:
        raise ValueError(
            f"EnergyFilter bins must be edge pairs, got shape {pairs.shape}"
        )
    _filter_bin_count(energy_filter)
    if not np.all(np.isfinite(pairs)):
        raise ValueError("EnergyFilter bin edges must be finite")
    if np.any(pairs[:, 1] <= pairs[:, 0]):
        raise ValueError("EnergyFilter bins must have ascending lower/upper edges")
    if pairs.shape[0] > 1 and not np.allclose(
        pairs[:-1, 1],
        pairs[1:, 0],
        rtol=1.0e-12,
        atol=0.0,
    ):
        raise ValueError("EnergyFilter bins must be contiguous")
    actual_bounds = np.concatenate((pairs[:1, 0], pairs[:, 1]))
    if actual_bounds.shape != expected_bounds.shape or not np.allclose(
        actual_bounds,
        expected_bounds,
        rtol=1.0e-12,
        atol=0.0,
    ):
        raise ValueError(
            "EnergyFilter bin edges do not match --mgxs /energy_bounds: "
            f"{actual_bounds.tolist()} != {expected_bounds.tolist()}"
        )


def _scalar_filter_ids(spatial_filter) -> tuple[int, ...]:
    bins = getattr(spatial_filter, "bins", None)
    if bins is None:
        raise ValueError(f"{type(spatial_filter).__name__} does not expose bins")
    ids: list[int] = []
    for index, raw_value in enumerate(bins):
        array = np.asarray(raw_value)
        if array.ndim != 0:
            raise ValueError(
                f"{type(spatial_filter).__name__} bin {index} is not a scalar "
                "domain ID; spatial order cannot be verified"
            )
        try:
            value = operator.index(array.item())
        except TypeError as exc:
            raise ValueError(
                f"{type(spatial_filter).__name__} bin {index} is not an integer "
                "domain ID; spatial order cannot be verified"
            ) from exc
        ids.append(int(value))
    if len(set(ids)) != len(ids):
        raise ValueError(
            f"{type(spatial_filter).__name__} domain IDs are not unique; "
            "spatial order cannot be verified"
        )
    return tuple(ids)


def _as_flux_array(
    values: np.ndarray | Iterable[Iterable[float]],
    *,
    names: tuple[str, ...],
    dataset_name: str,
    allow_zero: bool = False,
) -> np.ndarray:
    flux = np.asarray(values, dtype=float)
    if flux.ndim != 2:
        raise ValueError(f"{dataset_name} values must have shape (mixture, group)")
    if flux.shape[0] != len(names):
        raise ValueError(
            f"{dataset_name} mixture axis must match mixture_names: "
            f"{flux.shape[0]} != {len(names)}"
        )
    if flux.shape[1] <= 0:
        raise ValueError(f"{dataset_name} must contain at least one energy group")
    if not np.all(np.isfinite(flux)):
        raise ValueError(f"{dataset_name} values must be finite")
    if allow_zero:
        if np.any(flux < 0.0):
            raise ValueError(f"{dataset_name} values must be non-negative")
    elif np.any(flux <= 0.0):
        raise ValueError(f"{dataset_name} values must be positive")
    return flux


def _max_relative_std_dev(flux: np.ndarray, std_dev: np.ndarray) -> float:
    # Zero-flux bins (allow_zero) carry zero Monte Carlo std_dev; exclude them
    # from the ratio so the report stays finite.
    rel = np.zeros_like(std_dev, dtype=float)
    np.divide(std_dev, np.abs(flux), out=rel, where=flux != 0.0)
    return float(np.max(rel))


def _as_std_dev_array(
    values: np.ndarray | Iterable[Iterable[float]] | None,
    *,
    expected_shape: tuple[int, int],
    dataset_name: str,
) -> np.ndarray | None:
    if values is None:
        return None
    std_dev = np.asarray(values, dtype=float)
    if std_dev.shape != expected_shape:
        raise ValueError(
            f"{dataset_name} shape must match openmc_volume_flux shape: "
            f"{std_dev.shape} != {expected_shape}"
        )
    if not np.all(np.isfinite(std_dev)):
        raise ValueError(f"{dataset_name} values must be finite")
    if np.any(std_dev < 0.0):
        raise ValueError(f"{dataset_name} values must be non-negative")
    return std_dev


def _write_flux_dataset(
    h5,
    name: str,
    values: np.ndarray,
    *,
    names: tuple[str, ...],
    source_group_order: str,
    energy_bounds_verified: bool,
    spatial_domain_order_verified: bool,
    source_filter_order: tuple[str, ...],
    source_domain_ids: tuple[int, ...],
) -> None:
    dataset = h5.create_dataset(name, data=values)
    dataset.attrs["schema"] = SCHEMA
    dataset.attrs["package_version"] = __version__
    dataset.attrs["layout"] = "[mixture, group]"
    dataset.attrs["group_order"] = MGXS_DONJON_GROUP_ORDER
    dataset.attrs["source_group_order"] = str(source_group_order)
    dataset.attrs["mixture_names"] = np.asarray(names, dtype="S")
    dataset.attrs["energy_bounds_verified"] = bool(energy_bounds_verified)
    dataset.attrs["spatial_domain_order_verified"] = bool(
        spatial_domain_order_verified
    )
    dataset.attrs["source_filter_order"] = np.asarray(source_filter_order, dtype="S")
    dataset.attrs["source_domain_ids"] = np.asarray(source_domain_ids, dtype=np.int64)
