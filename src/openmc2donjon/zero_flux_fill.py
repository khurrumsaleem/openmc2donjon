"""Fill zero-flux-group XS in a converter MGXS handoff from an OpenMC MG macrolib.

Fast-spectrum cores carry literally zero Monte Carlo flux in the thermal
groups of fine fast meshes (e.g. ECCO-33), so flux-weighted MGXS tallies are
0/0 -> 0 there at any statistics. For those groups this module substitutes
the exact material cross sections from the OpenMC MG macrolib that the
transport run consumed (track-length tallies already reproduce the macrolib
to machine precision wherever flux is nonzero, so the substitution is
consistent by construction).

For each mixture and each group g with total == 0, or transport_total <= 0
when that dataset is present, or one of the opt-in noise criteria is met:

- total_std_dev > threshold * total
- sum(P0 scatter row) > (1 + threshold) * total

- total, absorption, fission, nu_fission        <- macrolib material data
- scatter_matrix[:, g, :] (all shared orders)   <- matching macrolib scatter rows
- reduced_absorption                            <- total - sum(P0 nu-scatter)
  for an explicitly declared multiplicity-weighted handoff
- transport_total = total - sum_g' P1(g -> g')   (material-macrolib
  Legendre out-scatter correction, not an OpenMC CE TransportXS reconstruction)
- matching *_std_dev entries                     <- 0 (exact library value)

Filled group indices are recorded per mixture in the
``zero_flux_filled_groups`` attribute (converter order, index 0 = highest
energy group) together with the ``zero_flux_fill_source`` provenance path.
When ``transport_total`` is touched, ``zero_flux_transport_method``
records ``macrolib_p1_outscatter``, ``macrolib_nu_p1_outscatter``, or
``macrolib_p0_total``.

Mixtures are matched to macrolib materials through a label attribute on the
mixture group (default ``irena_mixture_label``); pass ``label_attr`` to use
another attribute name.

An ordinary-scatter handoff retains the legacy interpretation of an
unannotated OpenMC MG macrolib. A multiplicity-weighted (nu-scatter) handoff
is stricter: substitution is allowed only when the source macrolib explicitly
declares the same canonical ``openmc_scatter_*`` contract and carries a valid
OpenMC ``multiplicity_matrix``. The latter proves that the source retained
ordinary absorption alongside its stored nu-scatter matrix. Ordinary scatter
is never substituted for nu-scatter.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
import os
from pathlib import Path
import shutil
import tempfile
from typing import Any

import numpy as np

from .mgxs_physics_checks import (
    _ResolvedScatterContract,
    _scatter_balance_contract,
)
from .openmc_provenance import (
    provenance_before_hdf5_mutation,
    refresh_openmc_provenance_after_hdf5_mutation,
)


SCHEMA = "openmc2donjon.zero-flux-fill.v2"
DEFAULT_LABEL_ATTR = "irena_mixture_label"
_SCATTER_CONTRACT_ATTRS = (
    "openmc_scatter_mgxs_type",
    "openmc_scatter_multiplicity_weighted",
    "openmc_scatter_balance_dataset",
)


@dataclass(frozen=True)
class _FillPlan:
    fill: np.ndarray
    total: np.ndarray
    absorption: np.ndarray
    balance: np.ndarray
    fission: np.ndarray | None
    nu_fission: np.ndarray | None
    scatter: np.ndarray
    transport: np.ndarray | None
    transport_method: str | None
    scatter_format: str
    scatter_order: int


@dataclass(frozen=True)
class ZeroFluxFillReport:
    input_h5: Path
    macrolib: Path
    output_h5: Path
    label_attr: str
    max_total_rel_std_dev: float | None
    max_scatter_row_overshoot_rel: float | None
    mixture_count: int
    filled_per_mixture: tuple[tuple[str, int], ...]
    total_filled_bins: int


def print_report(report: ZeroFluxFillReport) -> None:
    """Print the user-facing zero-flux fill report."""

    print("OpenMC-to-DONJON zero-flux fill")
    print(f"  schema: {SCHEMA}")
    print(f"  input: {report.input_h5}")
    print(f"  macrolib: {report.macrolib}")
    print(f"  output: {report.output_h5}")
    print(f"  label attribute: {report.label_attr}")
    print(f"  max total relative std_dev: {report.max_total_rel_std_dev}")
    print(
        "  max P0 scatter-row overshoot: "
        f"{report.max_scatter_row_overshoot_rel}"
    )
    print(f"  mixtures: {report.mixture_count}")
    print(f"  mixtures filled: {len(report.filled_per_mixture)}")
    print(f"  filled (mixture, group) bins: {report.total_filled_bins}")
    for name, count in report.filled_per_mixture:
        print(f"    {name}: {count}")
    print("")
    print("Zero-flux fill decision")
    print("  openmc2donjon_zero_flux_fill_passed")


def write_summary(path: Path, report: ZeroFluxFillReport) -> None:
    """Write a machine-readable zero-flux fill summary."""

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(summary_payload(report), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def summary_payload(report: ZeroFluxFillReport) -> dict[str, Any]:
    """Return the JSON-serializable zero-flux fill payload."""

    return {
        "schema": SCHEMA,
        "decision": "openmc2donjon_zero_flux_fill_passed",
        "input_h5": str(report.input_h5),
        "macrolib": str(report.macrolib),
        "output_h5": str(report.output_h5),
        "label_attr": report.label_attr,
        "max_total_rel_std_dev": report.max_total_rel_std_dev,
        "max_scatter_row_overshoot_rel": report.max_scatter_row_overshoot_rel,
        "mixtures": report.mixture_count,
        "filled_per_mixture": {name: count for name, count in report.filled_per_mixture},
        "total_filled_bins": report.total_filled_bins,
    }


def fill_zero_flux_groups(
    input_h5: Path,
    *,
    macrolib: Path,
    output_h5: Path | None = None,
    in_place: bool = False,
    label_attr: str = DEFAULT_LABEL_ATTR,
    force: bool = False,
    max_total_rel_std_dev: float | None = None,
    max_scatter_row_overshoot_rel: float | None = None,
) -> ZeroFluxFillReport:
    """Substitute macrolib XS into zero-flux (mixture, group) bins.

    By default the input file is copied to ``output_h5`` and the copy is
    edited; pass ``in_place=True`` (and no ``output_h5``) to edit the input
    file directly.

    ``max_total_rel_std_dev`` opts in a noise criterion: bins whose total XS
    carries a relative standard deviation above the threshold are filled as
    well. Micro-flux bins can collect a handful of scores (flux > 0, so the
    zero-flux criterion misses them) whose rate/flux ratio explodes into
    unphysical cross sections (degenerate single-score bins carry
    rel std = sqrt(2)); the macrolib value is strictly better there.

    ``max_scatter_row_overshoot_rel`` opts in a solver-stability criterion
    for ordinary scattering matrices: a P0 out-scatter row that exceeds the
    total cross section by more than the requested relative tolerance is
    unphysical and can make a deterministic source iteration diverge before
    statistical uncertainty alone identifies the bad bin. This criterion is
    forbidden for multiplicity-weighted scattering, where a P0 row above the
    total cross section can be physical.
    """

    input_h5 = Path(input_h5)
    macrolib = Path(macrolib)
    _validate_threshold("max_total_rel_std_dev", max_total_rel_std_dev)
    _validate_threshold(
        "max_scatter_row_overshoot_rel", max_scatter_row_overshoot_rel
    )
    if not input_h5.exists():
        raise FileNotFoundError(f"input HDF5 does not exist: {input_h5}")
    if not macrolib.exists():
        raise FileNotFoundError(f"macrolib does not exist: {macrolib}")
    if in_place:
        if output_h5 is not None:
            raise ValueError("--in-place cannot be combined with an output path")
        target_h5 = input_h5
    else:
        if output_h5 is None:
            raise ValueError("an output path is required unless in_place is set")
        target_h5 = Path(output_h5)
        if target_h5.resolve() == input_h5.resolve():
            raise ValueError("output HDF5 must be different from input HDF5; use in_place instead")
        if target_h5.exists() and not force:
            raise FileExistsError(f"output already exists; use --force to overwrite: {target_h5}")

    library = _load_macrolib(macrolib)
    by_name = {xsdata.name: xsdata for xsdata in library.xsdatas}
    target_contracts, filled_contracts = _preflight_scatter_contracts(
        input_h5,
        max_total_rel_std_dev=max_total_rel_std_dev,
        max_scatter_row_overshoot_rel=max_scatter_row_overshoot_rel,
    )
    if len(filled_contracts) > 1:
        raise ValueError(
            "zero-flux fill cannot use one OpenMC MG macrolib for a mixture "
            "of different scatter contracts"
        )
    if filled_contracts:
        target_contract = next(iter(filled_contracts))
        source_contract = _macrolib_scatter_contract(macrolib)
        if target_contract.multiplicity_weighted:
            if not source_contract.declared:
                raise ValueError(
                    "multiplicity-weighted (nu-scatter) zero-flux fill requires "
                    "the OpenMC MG macrolib to declare all three root attributes: "
                    + ", ".join(_SCATTER_CONTRACT_ATTRS)
                    + "; ordinary or undeclared scatter is not a valid substitute"
                )
            if not source_contract.multiplicity_weighted:
                raise ValueError(
                    "multiplicity-weighted (nu-scatter) handoff cannot be filled "
                    "from an ordinary-scatter OpenMC MG macrolib"
                )
            if (
                target_contract.mgxs_type is not None
                and source_contract.mgxs_type != target_contract.mgxs_type
            ):
                raise ValueError(
                    "nu-scatter zero-flux fill requires the same OpenMC MGXS "
                    f"estimator contract: handoff={target_contract.mgxs_type!r}, "
                    f"macrolib={source_contract.mgxs_type!r}"
                )
        elif source_contract.declared and source_contract.multiplicity_weighted:
            raise ValueError(
                "ordinary-scatter handoff cannot be filled from a "
                "multiplicity-weighted (nu-scatter) OpenMC MG macrolib"
            )
    fill_plans = _preflight_fill_plans(
        input_h5,
        macrolib=macrolib,
        by_name=by_name,
        target_contracts=target_contracts,
        label_attr=label_attr,
        max_total_rel_std_dev=max_total_rel_std_dev,
        max_scatter_row_overshoot_rel=max_scatter_row_overshoot_rel,
    )
    openmc_provenance = provenance_before_hdf5_mutation(input_h5)

    # Complete all changes (including the provenance refresh) in a sibling
    # file. A validation or HDF5 write failure must leave an in-place source
    # and any pre-existing --force destination byte-for-byte intact.
    target_h5.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=".openmc2donjon-zero-flux-", dir=target_h5.parent
    ) as staging_dir:
        staged_h5 = Path(staging_dir) / target_h5.name
        shutil.copy2(input_h5, staged_h5)
        mixture_count, filled_per_mixture = _apply_fill_plans(
            staged_h5,
            macrolib=macrolib,
            target_contracts=target_contracts,
            fill_plans=fill_plans,
        )
        refresh_openmc_provenance_after_hdf5_mutation(
            staged_h5,
            openmc_provenance,
        )
        os.replace(staged_h5, target_h5)

    return ZeroFluxFillReport(
        input_h5=input_h5,
        macrolib=macrolib,
        output_h5=target_h5,
        label_attr=label_attr,
        max_total_rel_std_dev=max_total_rel_std_dev,
        max_scatter_row_overshoot_rel=max_scatter_row_overshoot_rel,
        mixture_count=mixture_count,
        filled_per_mixture=tuple(filled_per_mixture),
        total_filled_bins=sum(count for _name, count in filled_per_mixture),
    )


def _apply_fill_plans(
    target_h5: Path,
    *,
    macrolib: Path,
    target_contracts: dict[str, _ResolvedScatterContract],
    fill_plans: dict[str, _FillPlan],
) -> tuple[int, list[tuple[str, int]]]:
    import h5py

    mixture_count = 0
    filled_per_mixture: list[tuple[str, int]] = []
    with h5py.File(target_h5, "r+") as h5:
        for name, group in h5["mixtures"].items():
            mixture_count += 1
            scatter_contract = target_contracts[str(name)]
            plan = fill_plans.get(str(name))
            if plan is None:
                continue
            fill = plan.fill
            _fill_dataset(group, fill, "total", plan.total)
            _fill_dataset(group, fill, "absorption", plan.absorption)
            if scatter_contract.multiplicity_weighted:
                _fill_dataset(
                    group,
                    fill,
                    scatter_contract.balance_dataset,
                    plan.balance,
                )
            if plan.fission is not None:
                _fill_dataset(group, fill, "fission", plan.fission)
                assert plan.nu_fission is not None
                _fill_dataset(group, fill, "nu_fission", plan.nu_fission)

            matrix = group["scatter_matrix"][:]
            # Do not retain noisy higher moments in a substituted material row
            # when the source macrolib carries fewer Legendre orders.
            matrix[:, fill, :] = 0.0
            for order in range(min(matrix.shape[0], plan.scatter.shape[0])):
                matrix[order][fill, :] = plan.scatter[order][fill, :]
            group["scatter_matrix"][...] = matrix
            if "scatter_matrix_std_dev" in group:
                std = group["scatter_matrix_std_dev"][:]
                std[:, fill, :] = 0.0
                group["scatter_matrix_std_dev"][...] = std

            if "transport_total" in group:
                assert plan.transport is not None
                assert plan.transport_method is not None
                _fill_dataset(group, fill, "transport_total", plan.transport)
                group.attrs["zero_flux_transport_method"] = plan.transport_method
                group.attrs["zero_flux_scatter_format"] = plan.scatter_format
                group.attrs["zero_flux_scatter_order"] = plan.scatter_order

            # Preserve provenance when a file is filled in more than one pass
            # (for example zero-flux first and an opt-in noise criterion
            # later).  Replacing the attribute would make earlier substituted
            # groups indistinguishable from untouched Monte Carlo data.
            previous_fill = np.asarray(
                group.attrs.get("zero_flux_filled_groups", ()), dtype=np.int64
            )
            group.attrs["zero_flux_filled_groups"] = np.union1d(
                previous_fill, fill
            ).astype(np.int64)
            group.attrs["zero_flux_fill_source"] = str(macrolib)
            filled_per_mixture.append((name, len(fill)))

    return mixture_count, filled_per_mixture


def _validate_threshold(name: str, value: float | None) -> None:
    if value is None:
        return
    if not np.isfinite(value) or value < 0.0:
        raise ValueError(f"{name} must be finite and non-negative")


def _preflight_fill_plans(
    input_h5: Path,
    *,
    macrolib: Path,
    by_name: dict[str, Any],
    target_contracts: dict[str, _ResolvedScatterContract],
    label_attr: str,
    max_total_rel_std_dev: float | None,
    max_scatter_row_overshoot_rel: float | None,
) -> dict[str, _FillPlan]:
    """Validate every selected source row before copying or editing output."""

    import h5py

    plans: dict[str, _FillPlan] = {}
    with h5py.File(input_h5, "r") as h5:
        for name, group in h5["mixtures"].items():
            fill = np.where(
                _fill_mask(
                    group,
                    max_total_rel_std_dev=max_total_rel_std_dev,
                    max_scatter_row_overshoot_rel=max_scatter_row_overshoot_rel,
                )
            )[0]
            if not len(fill):
                continue
            label = _mixture_label(group, str(name), label_attr)
            if label not in by_name:
                raise ValueError(
                    f"{name}: macrolib {macrolib} has no material named {label!r} "
                    f"(from mixture attribute {label_attr!r})"
                )
            plans[str(name)] = _validated_fill_plan(
                by_name[label],
                label=label,
                expected_groups=int(np.asarray(group["total"][:]).size),
                fill=fill,
                contract=target_contracts[str(name)],
                include_transport="transport_total" in group,
                include_fission=bool(
                    getattr(by_name[label], "fissionable", False)
                    and "fission" in group
                ),
            )
    return plans


def _validated_fill_plan(
    xsdata: Any,
    *,
    label: str,
    expected_groups: int,
    fill: np.ndarray,
    contract: _ResolvedScatterContract,
    include_transport: bool,
    include_fission: bool,
) -> _FillPlan:
    if len(xsdata.temperatures) != 1:
        raise ValueError(f"{label}: expected a single-temperature macrolib")
    temp_idx = 0
    scatter_format = str(getattr(xsdata, "scatter_format", "")).strip().lower()
    if scatter_format != "legendre":
        raise ValueError(
            f"{label}: zero-flux fill requires a Legendre macrolib; "
            f"scatter_format={scatter_format or '<missing>'!r}"
        )

    mac_total = _group_vector(xsdata.total[temp_idx])
    mac_absorption = _group_vector(xsdata.absorption[temp_idx])
    # An OpenMC MGXS Library generated with nu-scatter stores the nu-scatter
    # transfer matrix directly in XSData.scatter_matrix; multiplicity_matrix
    # is the separate nu-scatter/scatter ratio used by OpenMC transport. The
    # explicit source contract above selects the stored representation, so it
    # must not be multiplied by that ratio a second time.
    scatter = _dense_scatter(xsdata, temp_idx)
    if (
        mac_total.shape != (expected_groups,)
        or mac_absorption.shape != (expected_groups,)
        or scatter.ndim != 3
        or scatter.shape[1:] != (expected_groups, expected_groups)
    ):
        raise ValueError(
            f"{label}: macrolib group dimensions do not match the "
            f"{expected_groups}-group handoff"
        )
    if (
        not np.all(np.isfinite(mac_total))
        or np.any(mac_total <= 0.0)
        or not np.all(np.isfinite(mac_absorption))
        or np.any(mac_absorption < 0.0)
        or not np.all(np.isfinite(scatter))
        or np.any(scatter[0] < 0.0)
    ):
        raise ValueError(
            f"{label}: macrolib total/absorption/P0 scatter data must "
            "be finite and physically signed"
        )

    p0_outscatter = scatter[0].sum(axis=1)
    if contract.multiplicity_weighted:
        multiplicity_tables = getattr(xsdata, "multiplicity_matrix", None)
        multiplicity = (
            None
            if multiplicity_tables is None
            else multiplicity_tables[temp_idx]
        )
        if multiplicity is None:
            raise ValueError(
                f"{label}: nu-scatter zero-flux fill requires the source "
                "OpenMC XSData multiplicity_matrix; without it, source "
                "absorption may already be reduced rather than ordinary"
            )
        multiplicity = np.asarray(multiplicity, dtype=float)
        if (
            multiplicity.shape != (expected_groups, expected_groups)
            or not np.all(np.isfinite(multiplicity))
            or np.any(multiplicity < 0.0)
        ):
            raise ValueError(
                f"{label}: source multiplicity_matrix must match the group "
                "structure and be finite and non-negative"
            )
        mac_balance = mac_total - p0_outscatter
        if not np.all(np.isfinite(mac_balance)):
            raise ValueError(
                f"{label}: macrolib total - P0 nu-scatter is non-finite"
            )
    else:
        mac_balance = mac_absorption
        p0_residual = mac_total - mac_balance - p0_outscatter
        p0_tolerance = 1.0e-12 * np.maximum(np.abs(mac_total), 1.0)
        if np.any(p0_residual < -p0_tolerance):
            worst = int(np.argmin(p0_residual))
            raise ValueError(
                f"{label}: macrolib P0 scatter plus absorption exceeds "
                f"total in group {worst + 1}"
            )

    transport_method: str | None = None
    mac_transport: np.ndarray | None = None
    if include_transport:
        if scatter.shape[0] > 1:
            correction = scatter[1].sum(axis=1)
            transport_method = (
                "macrolib_nu_p1_outscatter"
                if contract.multiplicity_weighted
                else "macrolib_p1_outscatter"
            )
        else:
            correction = np.zeros_like(mac_total)
            transport_method = "macrolib_p0_total"
        mac_transport = mac_total - correction
        selected_transport = mac_transport[fill]
        if (
            not np.all(np.isfinite(selected_transport))
            or np.any(selected_transport <= 0.0)
        ):
            raise ValueError(
                f"{label}: {transport_method} produced non-positive or "
                "non-finite transport_total in a group selected for fill"
            )

    mac_fission: np.ndarray | None = None
    mac_nu_fission: np.ndarray | None = None
    if include_fission:
        if xsdata.fission is None or xsdata.nu_fission is None:
            raise ValueError(
                f"{label}: fissionable macrolib material is missing fission data"
            )
        mac_fission = _group_vector(xsdata.fission[temp_idx])
        mac_nu_fission = _group_vector(xsdata.nu_fission[temp_idx])
        if (
            mac_fission.shape != (expected_groups,)
            or mac_nu_fission.shape != (expected_groups,)
            or not np.all(np.isfinite(mac_fission))
            or not np.all(np.isfinite(mac_nu_fission))
            or np.any(mac_fission < 0.0)
            or np.any(mac_nu_fission < 0.0)
        ):
            raise ValueError(
                f"{label}: macrolib fission data must match the group structure "
                "and be finite and non-negative"
            )

    return _FillPlan(
        fill=fill,
        total=mac_total,
        absorption=mac_absorption,
        balance=mac_balance,
        fission=mac_fission,
        nu_fission=mac_nu_fission,
        scatter=scatter,
        transport=mac_transport,
        transport_method=transport_method,
        scatter_format=scatter_format,
        scatter_order=int(getattr(xsdata, "order", scatter.shape[0] - 1)),
    )


def _preflight_scatter_contracts(
    input_h5: Path,
    *,
    max_total_rel_std_dev: float | None,
    max_scatter_row_overshoot_rel: float | None,
) -> tuple[
    dict[str, _ResolvedScatterContract], tuple[_ResolvedScatterContract, ...]
]:
    """Resolve contracts and identify contracts whose rows will be replaced."""

    import h5py

    contracts: dict[str, _ResolvedScatterContract] = {}
    filled_contracts: dict[
        tuple[str | None, bool, str, bool], _ResolvedScatterContract
    ] = {}
    with h5py.File(input_h5, "r") as h5:
        mixtures = h5["mixtures"]
        for name, group in mixtures.items():
            _validate_target_fill_arrays(group, h5)
            contract = _scatter_balance_contract(
                group,
                h5,
                None,
                transport_total_present="transport_total" in group,
            )
            if contract.issue:
                raise ValueError(f"{name}: {contract.issue}")
            contracts[str(name)] = contract
            if contract.multiplicity_weighted:
                if max_scatter_row_overshoot_rel is not None:
                    raise ValueError(
                        "max_scatter_row_overshoot_rel is defined only for "
                        "ordinary scatter; a nu-scatter row may physically "
                        "exceed total"
                    )
                if contract.balance_dataset not in group:
                    raise ValueError(
                        f"{name}: multiplicity-weighted scatter requires "
                        f"dataset {contract.balance_dataset!r} before zero-flux fill"
                    )
            if np.any(
                _fill_mask(
                    group,
                    max_total_rel_std_dev=max_total_rel_std_dev,
                    max_scatter_row_overshoot_rel=max_scatter_row_overshoot_rel,
                )
            ):
                key = (
                    contract.mgxs_type,
                    contract.multiplicity_weighted,
                    contract.balance_dataset,
                    contract.declared,
                )
                filled_contracts.setdefault(key, contract)
    return contracts, tuple(filled_contracts.values())


def _validate_target_fill_arrays(group: Any, h5: Any) -> None:
    """Check every array that row substitution reads or writes before staging."""

    if "total" not in group:
        raise ValueError(f"{group.name}: zero-flux fill requires total")
    total_shape = getattr(group["total"], "shape", None)
    if total_shape is None or len(total_shape) != 1 or total_shape[0] < 1:
        raise ValueError(f"{group.name}: total must be a non-empty group-wise vector")
    energy_groups = int(h5.attrs.get("energy_groups", total_shape[0]))
    expected = (energy_groups,)
    vector_names = (
        "total", "absorption", "reduced_absorption", "fission",
        "nu_fission", "transport_total",
    )
    for name in vector_names:
        if name not in group:
            if name in {"total", "absorption"}:
                raise ValueError(f"{group.name}: zero-flux fill requires {name}")
            continue
        _require_target_array(group[name], expected)
        std_name = f"{name}_std_dev"
        if std_name in group:
            _require_target_array(group[std_name], expected)

    if "scatter_matrix" not in group:
        raise ValueError(f"{group.name}: zero-flux fill requires scatter_matrix")
    scatter_shape = getattr(group["scatter_matrix"], "shape", None)
    if (
        scatter_shape is None
        or len(scatter_shape) != 3
        or scatter_shape[0] < 1
        or scatter_shape[1:] != (energy_groups, energy_groups)
    ):
        raise ValueError(
            f"{group.name}: scatter_matrix must have shape "
            f"(moment, {energy_groups}, {energy_groups})"
        )
    _require_target_array(group["scatter_matrix"], scatter_shape)
    if "scatter_matrix_std_dev" in group:
        _require_target_array(group["scatter_matrix_std_dev"], scatter_shape)
    for attrs in (h5.attrs, group.attrs):
        if "scatter_axes" in attrs:
            axes = attrs["scatter_axes"]
            if isinstance(axes, bytes):
                axes = axes.decode()
            if str(axes).replace(" ", "") != "moment,from,to":
                raise ValueError(
                    f"{group.name}: zero-flux fill requires "
                    "scatter_axes='moment,from,to'"
                )


def _require_target_array(dataset: Any, expected: tuple[int, ...]) -> None:
    if getattr(dataset, "shape", None) != expected:
        raise ValueError(f"{dataset.name}: expected shape {expected}")
    if dataset.dtype.kind != "f":
        raise ValueError(
            f"{dataset.name}: zero-flux fill requires floating-point storage "
            "to preserve cross-section values"
        )


def _fill_mask(
    group: Any,
    *,
    max_total_rel_std_dev: float | None,
    max_scatter_row_overshoot_rel: float | None,
) -> np.ndarray:
    total = np.asarray(group["total"][:], dtype=float)
    fill_mask = total == 0.0
    if "transport_total" in group:
        fill_mask |= np.asarray(group["transport_total"][:], dtype=float) <= 0.0
    if max_total_rel_std_dev is not None and "total_std_dev" in group:
        std = np.asarray(group["total_std_dev"][:], dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            rel = np.where(total > 0.0, std / total, 0.0)
        fill_mask |= rel > max_total_rel_std_dev
    if max_scatter_row_overshoot_rel is not None:
        matrix = np.asarray(group["scatter_matrix"][:], dtype=float)
        p0_outscatter = matrix[0].sum(axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            overshoot_rel = np.where(
                total > 0.0,
                (p0_outscatter - total) / total,
                0.0,
            )
        fill_mask |= overshoot_rel > max_scatter_row_overshoot_rel
    return fill_mask


def _macrolib_scatter_contract(macrolib: Path) -> _ResolvedScatterContract:
    """Read an optional explicit contract from an OpenMC MG macrolib root."""

    import h5py

    try:
        with h5py.File(macrolib, "r") as h5:
            present = [name for name in _SCATTER_CONTRACT_ATTRS if name in h5.attrs]
            if not present:
                contract = _scatter_balance_contract(
                    h5,
                    h5,
                    None,
                    transport_total_present=(
                        "openmc_transport_mgxs_type" in h5.attrs
                    ),
                )
            else:
                missing = [
                    name for name in _SCATTER_CONTRACT_ATTRS if name not in h5.attrs
                ]
                if missing:
                    raise ValueError(
                        f"{macrolib}: incomplete OpenMC MG macrolib scatter "
                        "contract; missing root attributes: "
                        + ", ".join(missing)
                    )
                contract = _scatter_balance_contract(
                    h5,
                    h5,
                    None,
                    transport_total_present=(
                        "openmc_transport_mgxs_type" in h5.attrs
                    ),
                )
    except OSError as exc:
        raise ValueError(f"cannot inspect OpenMC MG macrolib scatter contract: {macrolib}: {exc}") from exc
    if contract.issue:
        raise ValueError(f"{macrolib}: {contract.issue}")
    return contract


def _load_macrolib(macrolib: Path) -> Any:
    import openmc  # type: ignore[import-not-found]

    return openmc.MGXSLibrary.from_hdf5(str(macrolib))


def _mixture_label(group: Any, name: str, label_attr: str) -> str:
    if label_attr not in group.attrs:
        raise ValueError(f"{name}: mixture is missing the label attribute {label_attr!r}")
    label = group.attrs[label_attr]
    return label.decode() if isinstance(label, bytes) else str(label)


def _dense_scatter(xsdata: Any, temp_idx: int) -> np.ndarray:
    """Return dense (order, g_in, g_out) scatter matrix in converter order
    (index 0 = highest energy; macrolib storage is ascending energy)."""
    matrix = np.transpose(np.asarray(xsdata.scatter_matrix[temp_idx]), (2, 0, 1))
    return matrix[:, ::-1, ::-1]


def _group_vector(values: Any) -> np.ndarray:
    """Reverse an ascending-energy macrolib vector into converter order."""
    return np.asarray(values, dtype=float)[::-1]


def _fill_dataset(group: Any, fill: np.ndarray, key: str, values: np.ndarray) -> None:
    data = group[key][:]
    data[fill] = values[fill]
    group[key][...] = data
    std_key = f"{key}_std_dev"
    if std_key in group:
        std = group[std_key][:]
        std[fill] = 0.0
        group[std_key][...] = std
