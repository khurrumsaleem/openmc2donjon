"""Shared physics-consistency checks for converter-facing MGXS HDF5 files."""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from .constants import MGXS_DONJON_GROUP_ORDER
from .energy_groups import validate_energy_bounds_internal
from .hdf5_names import decode_hdf5_names
from .mgxs_input_scatter import (
    MOMENT_FIRST_SCATTER_AXES,
    MOMENT_LAST_SCATTER_AXES,
    normalize_axes,
    p0_scatter_matrix,
)


DEFAULT_SCATTER_ROW_BALANCE_REL = 5.0e-2
DEFAULT_CHI_SUM_TOLERANCE = 1.0e-6
# There is no universal physical interval for nu*Sigma_f / Sigma_f.  The
# effective neutron yield depends on nuclide mix and incident energy.  These
# names remain import-compatible, but range warnings are opt-in only.
DEFAULT_NU_RATIO_MINIMUM: float | None = None
DEFAULT_NU_RATIO_MAXIMUM: float | None = None
DEFAULT_TRANSPORT_P1_REL = 5.0e-2
LOCAL_ENERGY_BOUNDS_RTOL = 1.0e-10
LOCAL_ENERGY_BOUNDS_ATOL = 0.0
# Ordinary scattering balances against absorption. Multiplicity-weighted
# nu-scattering balances against OpenMC's reduced absorption.
ORDINARY_OPENMC_SCATTER_MGXS_TYPES = frozenset(
    {
        "scatter matrix",
        "consistent scatter matrix",
    }
)
NU_WEIGHTED_OPENMC_SCATTER_MGXS_TYPES = frozenset(
    {
        "nu scatter matrix",
        "consistent nu scatter matrix",
    }
)
OPENMC_SCATTER_BALANCE_DATASETS = frozenset(
    {
        "absorption",
        "reduced_absorption",
    }
)
OPENMC_TRANSPORT_MGXS_TYPES = frozenset({"transport", "nu transport"})


@dataclass(frozen=True)
class _ResolvedScatterContract:
    mgxs_type: str | None
    multiplicity_weighted: bool
    balance_dataset: str
    declared: bool
    transport_mgxs_type: str | None
    transport_declared: bool
    issue: str | None = None


@dataclass(frozen=True)
class MgxsPhysicsCheckReport:
    energy_bounds_local_count: int = 0
    energy_bounds_consistency_errors: tuple[str, ...] = ()
    scatter_row_balance_checked: int = 0
    scatter_row_balance_max_rel: float | None = None
    scatter_row_balance_max_abs: float | None = None
    scatter_row_balance_worst: str | None = None
    scatter_row_balance_warnings: tuple[str, ...] = ()
    scatter_row_balance_errors: tuple[str, ...] = ()
    openmc_scatter_mgxs_type: str | None = None
    openmc_scatter_multiplicity_weighted: bool | None = None
    openmc_scatter_balance_dataset: str | None = None
    openmc_scatter_contract_declared: bool | None = None
    openmc_scatter_contract_valid: bool | None = None
    openmc_transport_mgxs_type: str | None = None
    openmc_transport_contract_declared: bool | None = None
    chi_checked: int = 0
    chi_sum_max_abs_error: float | None = None
    chi_sum_worst: str | None = None
    chi_errors: tuple[str, ...] = ()
    nu_ratio_checked_bins: int = 0
    nu_ratio_min: float | None = None
    nu_ratio_max: float | None = None
    nu_ratio_worst: str | None = None
    nu_ratio_warning_count: int = 0
    nu_ratio_warnings: tuple[str, ...] = ()
    nu_ratio_support_mismatch_count: int = 0
    nu_ratio_errors: tuple[str, ...] = ()
    adf_calculations: int = 0
    adf_faces: tuple[str, ...] = ()
    adf_face_errors: tuple[str, ...] = ()
    transport_p1_checked: int = 0
    transport_p1_skipped: int = 0
    transport_p1_max_rel: float | None = None
    transport_p1_max_abs: float | None = None
    transport_p1_worst: str | None = None
    transport_p1_errors: tuple[str, ...] = ()


@dataclass
class _MutablePhysicsReport:
    energy_bounds_local_count: int = 0
    energy_bounds_consistency_errors: list[str] = field(default_factory=list)
    scatter_row_balance_checked: int = 0
    scatter_row_balance_max_rel: float | None = None
    scatter_row_balance_max_abs: float | None = None
    scatter_row_balance_worst: str | None = None
    scatter_row_balance_warnings: list[str] = field(default_factory=list)
    scatter_row_balance_errors: list[str] = field(default_factory=list)
    scatter_contracts: list[tuple[str, _ResolvedScatterContract]] = field(
        default_factory=list
    )
    scatter_contract_errors: list[str] = field(default_factory=list)
    chi_checked: int = 0
    chi_sum_max_abs_error: float | None = None
    chi_sum_worst: str | None = None
    chi_errors: list[str] = field(default_factory=list)
    nu_ratio_checked_bins: int = 0
    nu_ratio_min: float | None = None
    nu_ratio_max: float | None = None
    nu_ratio_worst: str | None = None
    nu_ratio_warnings: list[str] = field(default_factory=list)
    nu_ratio_support_mismatch_count: int = 0
    nu_ratio_errors: list[str] = field(default_factory=list)
    adf_calculations: int = 0
    adf_faces: tuple[str, ...] = ()
    adf_face_errors: list[str] = field(default_factory=list)
    transport_p1_checked: int = 0
    transport_p1_skipped: int = 0
    transport_p1_max_rel: float | None = None
    transport_p1_max_abs: float | None = None
    transport_p1_worst: str | None = None
    transport_p1_errors: list[str] = field(default_factory=list)

    def freeze(self) -> MgxsPhysicsCheckReport:
        scatter_contracts = {
            (
                contract.mgxs_type,
                contract.multiplicity_weighted,
                contract.balance_dataset,
                contract.declared,
            )
            for _label, contract in self.scatter_contracts
        }
        common_scatter = (
            next(iter(scatter_contracts)) if len(scatter_contracts) == 1 else None
        )
        resolved_transport = [
            contract.transport_mgxs_type
            for _label, contract in self.scatter_contracts
        ]
        present_transport_types = {
            value for value in resolved_transport if value is not None
        }
        if not present_transport_types:
            common_transport: tuple[str | None, bool] | None = (None, False)
        elif len(present_transport_types) == 1 and all(
            value is not None for value in resolved_transport
        ):
            common_transport = (
                next(iter(present_transport_types)),
                all(
                    contract.transport_declared
                    for _label, contract in self.scatter_contracts
                ),
            )
        elif len(present_transport_types) == 1:
            # A P0 handoff may legitimately carry transport_total for only a
            # subset of calculations. There is no file-wide transport source
            # contract in that case, but the scatter/removal contract remains
            # valid and each present dataset is still checked locally.
            common_transport = (None, False)
        else:
            common_transport = None
        contract_present = bool(self.scatter_contracts)
        return MgxsPhysicsCheckReport(
            energy_bounds_local_count=self.energy_bounds_local_count,
            energy_bounds_consistency_errors=tuple(self.energy_bounds_consistency_errors),
            scatter_row_balance_checked=self.scatter_row_balance_checked,
            scatter_row_balance_max_rel=self.scatter_row_balance_max_rel,
            scatter_row_balance_max_abs=self.scatter_row_balance_max_abs,
            scatter_row_balance_worst=self.scatter_row_balance_worst,
            scatter_row_balance_warnings=tuple(self.scatter_row_balance_warnings),
            scatter_row_balance_errors=tuple(self.scatter_row_balance_errors),
            openmc_scatter_mgxs_type=(
                None if common_scatter is None else common_scatter[0]
            ),
            openmc_scatter_multiplicity_weighted=(
                None if common_scatter is None else common_scatter[1]
            ),
            openmc_scatter_balance_dataset=(
                None if common_scatter is None else common_scatter[2]
            ),
            openmc_scatter_contract_declared=(
                None
                if not contract_present
                else (None if common_scatter is None else common_scatter[3])
            ),
            openmc_scatter_contract_valid=(
                None
                if not contract_present
                else not self.scatter_contract_errors
            ),
            openmc_transport_mgxs_type=(
                None if common_transport is None else common_transport[0]
            ),
            openmc_transport_contract_declared=(
                None
                if not contract_present
                else (None if common_transport is None else common_transport[1])
            ),
            chi_checked=self.chi_checked,
            chi_sum_max_abs_error=self.chi_sum_max_abs_error,
            chi_sum_worst=self.chi_sum_worst,
            chi_errors=tuple(self.chi_errors),
            nu_ratio_checked_bins=self.nu_ratio_checked_bins,
            nu_ratio_min=self.nu_ratio_min,
            nu_ratio_max=self.nu_ratio_max,
            nu_ratio_worst=self.nu_ratio_worst,
            nu_ratio_warning_count=len(self.nu_ratio_warnings),
            nu_ratio_warnings=tuple(self.nu_ratio_warnings),
            nu_ratio_support_mismatch_count=self.nu_ratio_support_mismatch_count,
            nu_ratio_errors=tuple(self.nu_ratio_errors),
            adf_calculations=self.adf_calculations,
            adf_faces=self.adf_faces,
            adf_face_errors=tuple(self.adf_face_errors),
            transport_p1_checked=self.transport_p1_checked,
            transport_p1_skipped=self.transport_p1_skipped,
            transport_p1_max_rel=self.transport_p1_max_rel,
            transport_p1_max_abs=self.transport_p1_max_abs,
            transport_p1_worst=self.transport_p1_worst,
            transport_p1_errors=tuple(self.transport_p1_errors),
        )


def evaluate_mgxs_physics(
    h5: Any,
    *,
    mixture_names: tuple[str, ...],
    energy_groups: int,
    legendre_order: int,
    root_energy_bounds: np.ndarray | None,
    energy_bounds_consistency: bool = False,
    scatter_row_balance_rel: float | None = None,
    scatter_row_balance_warn_rel: float | None = None,
    chi_sum_tolerance: float | None = None,
    require_adf_face_consistency: bool = False,
    transport_p1_rel: float | None = None,
    nu_ratio_minimum: float | None = DEFAULT_NU_RATIO_MINIMUM,
    nu_ratio_maximum: float | None = DEFAULT_NU_RATIO_MAXIMUM,
) -> MgxsPhysicsCheckReport:
    """Evaluate production physics guardrails without mutating the HDF5 file."""

    report = _MutablePhysicsReport()
    if energy_bounds_consistency and root_energy_bounds is not None:
        _check_local_energy_bounds(
            report,
            root_energy_bounds=root_energy_bounds,
            h5=h5,
            mixture_names=mixture_names,
            energy_groups=energy_groups,
        )

    reference_flux = _root_reference_flux_by_mixture(
        h5,
        mixture_names=mixture_names,
        energy_groups=energy_groups,
    )
    adf_names_by_calc: list[tuple[str, ...]] = []
    for mixture_name, label, group, parent_group in _iter_calculations(
        h5, mixture_names
    ):
        parent_attrs = None if parent_group is None else parent_group.attrs
        fissionable = bool(
            _attr_with_parent(group.attrs, parent_attrs, "fissionable", False)
        )
        axes = _scatter_axes(group, h5, parent_group)
        total = _vector_or_none(group, "total", energy_groups)
        scatter_contract = _scatter_balance_contract(
            group,
            h5,
            parent_group,
            transport_total_present="transport_total" in group,
        )
        report.scatter_contracts.append((label, scatter_contract))
        if scatter_contract.issue is not None:
            issue = f"{label}: {scatter_contract.issue}"
            report.scatter_row_balance_errors.append(issue)
            report.scatter_contract_errors.append(issue)
        scatter_balance_name = scatter_contract.balance_dataset
        scatter_balance = _vector_or_none(
            group,
            scatter_balance_name,
            energy_groups,
        )
        scatter_balance_issue = _scatter_balance_vector_issue(
            label=label,
            contract=scatter_contract,
            balance=scatter_balance,
            dataset_present=scatter_balance_name in group,
            energy_groups=energy_groups,
        )
        if scatter_balance_issue is not None:
            report.scatter_row_balance_errors.append(scatter_balance_issue)
            report.scatter_contract_errors.append(scatter_balance_issue)
        fission = _vector_or_none(group, "fission", energy_groups)
        nu_fission = _vector_or_none(group, "nu_fission", energy_groups)
        chi = _vector_or_none(group, "chi", energy_groups)
        scatter = _scatter_or_none(group, "scatter_matrix")

        if (
            scatter_row_balance_rel is not None
            or scatter_row_balance_warn_rel is not None
        ) and scatter_contract.issue is None and scatter_balance_issue is None:
            _check_scatter_row_balance(
                report,
                label=label,
                total=total,
                balance=scatter_balance,
                balance_name=scatter_balance_name,
                scatter=scatter,
                axes=axes,
                energy_groups=energy_groups,
                legendre_order=legendre_order,
                fail_threshold=scatter_row_balance_rel,
                warn_threshold=scatter_row_balance_warn_rel,
            )
        if chi_sum_tolerance is not None:
            _check_chi(
                report,
                label=label,
                chi=chi,
                fissionable=fissionable,
                tolerance=chi_sum_tolerance,
            )
        _check_nu_ratio(
            report,
            label=label,
            fission=fission,
            nu_fission=nu_fission,
            fissionable=fissionable,
            minimum=nu_ratio_minimum,
            maximum=nu_ratio_maximum,
        )
        if require_adf_face_consistency:
            adf_names_by_calc.append(_adf_names(group))
        if transport_p1_rel is not None:
            _check_transport_p1(
                report,
                label=label,
                total=total,
                transport_total=_vector_or_none(group, "transport_total", energy_groups),
                scatter=scatter,
                axes=axes,
                energy_groups=energy_groups,
                legendre_order=legendre_order,
                reference_flux=(
                    reference_flux.get(mixture_name)
                    if parent_group is None
                    else None
                ),
                sph_applied=h5.attrs.get("sph_applied", False),
                sph_apply_operator=h5.attrs.get("sph_apply_operator"),
                applied_sph=(
                    _vector_or_none(group, "applied_sph", energy_groups)
                    if getattr(group.get("applied_sph"), "shape", None) == (energy_groups,)
                    else None
                ),
                threshold=transport_p1_rel,
            )

    if require_adf_face_consistency:
        _finalize_adf_faces(report, adf_names_by_calc)
    _finalize_scatter_contract_consistency(report)
    return report.freeze()


def _check_local_energy_bounds(
    report: _MutablePhysicsReport,
    *,
    root_energy_bounds: np.ndarray,
    h5: Any,
    mixture_names: tuple[str, ...],
    energy_groups: int,
) -> None:
    mixtures = h5["mixtures"]
    for mixture_name in mixture_names:
        mixture = mixtures[mixture_name]
        if "energy_bounds" in mixture:
            _check_one_local_energy_bounds(
                report,
                label=f"{mixture_name}/energy_bounds",
                obj=mixture["energy_bounds"],
                root_energy_bounds=root_energy_bounds,
                energy_groups=energy_groups,
            )
        if "states" not in mixture:
            continue
        states = mixture["states"]
        for state_name in _sorted_state_names(states):
            state_group = states[state_name]
            if "energy_bounds" in state_group:
                _check_one_local_energy_bounds(
                    report,
                    label=f"{mixture_name}/states/{state_name}/energy_bounds",
                    obj=state_group["energy_bounds"],
                    root_energy_bounds=root_energy_bounds,
                    energy_groups=energy_groups,
                )


def _check_one_local_energy_bounds(
    report: _MutablePhysicsReport,
    *,
    label: str,
    obj: Any,
    root_energy_bounds: np.ndarray,
    energy_groups: int,
) -> None:
    report.energy_bounds_local_count += 1
    try:
        values = np.asarray(obj[:], dtype=float)
    except (TypeError, ValueError, OSError):
        report.energy_bounds_consistency_errors.append(f"{label} must be numeric")
        return
    issues = validate_energy_bounds_internal(
        values,
        expected_groups=energy_groups,
        expected_order="ascending",
    )
    if issues:
        report.energy_bounds_consistency_errors.extend(
            f"{label}: {issue}" for issue in issues
        )
        return
    if not np.allclose(
        values,
        root_energy_bounds,
        rtol=LOCAL_ENERGY_BOUNDS_RTOL,
        atol=LOCAL_ENERGY_BOUNDS_ATOL,
    ):
        index = int(np.argmax(np.abs(values - root_energy_bounds)))
        report.energy_bounds_consistency_errors.append(
            f"{label} differs from /energy_bounds at index {index}: "
            f"actual={values[index]:.12e} expected={root_energy_bounds[index]:.12e}"
        )


def _scatter_balance_vector_issue(
    *,
    label: str,
    contract: _ResolvedScatterContract,
    balance: np.ndarray | None,
    dataset_present: bool,
    energy_groups: int,
) -> str | None:
    if contract.issue is not None or contract.balance_dataset != "reduced_absorption":
        return None
    if balance is not None and np.all(np.isfinite(balance)):
        return None
    qualifier = "missing" if not dataset_present else "invalid"
    mgxs_label = contract.mgxs_type or "nu-weighted scatter matrix"
    return (
        f"{label}: explicit nu-weighted OpenMC scattering {mgxs_label!r} "
        "requires a finite reduced_absorption vector "
        f"with {energy_groups} values; dataset is {qualifier}"
    )


def _check_scatter_row_balance(
    report: _MutablePhysicsReport,
    *,
    label: str,
    total: np.ndarray | None,
    balance: np.ndarray | None,
    balance_name: str,
    scatter: np.ndarray | None,
    axes: str | None,
    energy_groups: int,
    legendre_order: int,
    fail_threshold: float | None,
    warn_threshold: float | None,
) -> None:
    if balance is None or not np.all(np.isfinite(balance)):
        return
    if total is None or scatter is None:
        return
    p0 = p0_scatter_matrix(scatter, axes, energy_groups, legendre_order)
    if p0 is None:
        return
    if not (np.all(np.isfinite(total)) and np.all(np.isfinite(p0))):
        return
    report.scatter_row_balance_checked += 1
    residual = total - balance - p0.sum(axis=1)
    rel, max_abs, max_rel, index = _relative_worst(residual, total)
    del rel
    _update_worst(
        current=report,
        attr_rel="scatter_row_balance_max_rel",
        attr_abs="scatter_row_balance_max_abs",
        attr_worst="scatter_row_balance_worst",
        max_rel=max_rel,
        max_abs=max_abs,
        worst=f"{label}: group={index + 1} residual={residual[index]:.6e}",
    )
    detail = (
        "scatter row-balance max relative residual "
        f"{max_rel:.6e} (abs {max_abs:.6e}) at {label}: group={index + 1} "
        f"using {balance_name}"
    )
    if fail_threshold is not None and max_rel > fail_threshold:
        report.scatter_row_balance_errors.append(
            f"{detail} exceeds fail threshold {fail_threshold:.6e}"
        )
    elif warn_threshold is not None and max_rel > warn_threshold:
        report.scatter_row_balance_warnings.append(
            f"{detail} exceeds warn threshold {warn_threshold:.6e}"
        )


def _finalize_scatter_contract_consistency(
    report: _MutablePhysicsReport,
) -> None:
    contracts: dict[tuple[str | None, bool, str, bool], list[str]] = {}
    for label, contract in report.scatter_contracts:
        key = (
            contract.mgxs_type,
            contract.multiplicity_weighted,
            contract.balance_dataset,
            contract.declared,
        )
        contracts.setdefault(key, []).append(label)
    if len(contracts) <= 1:
        return
    detail = "; ".join(
        f"type={mgxs_type!r}, multiplicity_weighted={weighted}, "
        f"balance={balance!r}, declared={declared} at {','.join(labels)}"
        for (mgxs_type, weighted, balance, declared), labels in contracts.items()
    )
    report.scatter_row_balance_errors.append(
        "inconsistent resolved OpenMC scatter contract across calculations: "
        f"{detail}"
    )
    report.scatter_contract_errors.append(
        "inconsistent resolved OpenMC scatter contract across calculations: "
        f"{detail}"
    )


def _check_chi(
    report: _MutablePhysicsReport,
    *,
    label: str,
    chi: np.ndarray | None,
    fissionable: bool,
    tolerance: float,
) -> None:
    if not fissionable:
        return
    report.chi_checked += 1
    if chi is None:
        report.chi_errors.append(f"mixture {label}: chi is required for fissionable data")
        return
    if not np.all(np.isfinite(chi)):
        report.chi_errors.append(f"mixture {label}: chi contains non-finite values")
        return
    if np.any(chi < 0.0):
        report.chi_errors.append(f"mixture {label}: chi must be non-negative")
        return
    error = abs(float(np.sum(chi)) - 1.0)
    if report.chi_sum_max_abs_error is None or error > report.chi_sum_max_abs_error:
        report.chi_sum_max_abs_error = error
        report.chi_sum_worst = f"{label}: sum(chi)={float(np.sum(chi)):.12e}"
    if error > tolerance:
        report.chi_errors.append(
            f"mixture {label}: chi sum error {error:.6e} exceeds "
            f"tolerance {tolerance:.6e}"
        )


def _check_nu_ratio(
    report: _MutablePhysicsReport,
    *,
    label: str,
    fission: np.ndarray | None,
    nu_fission: np.ndarray | None,
    fissionable: bool,
    minimum: float | None,
    maximum: float | None,
) -> None:
    if not fissionable or fission is None or nu_fission is None:
        return
    if not (np.all(np.isfinite(fission)) and np.all(np.isfinite(nu_fission))):
        return
    fission_positive = fission > 0.0
    nu_fission_positive = nu_fission > 0.0
    support_mismatch = fission_positive != nu_fission_positive
    if np.any(support_mismatch):
        mismatch_groups = np.flatnonzero(support_mismatch)
        report.nu_ratio_support_mismatch_count += int(mismatch_groups.size)
        rendered_groups = ", ".join(str(int(index) + 1) for index in mismatch_groups)
        report.nu_ratio_errors.append(
            f"mixture {label}: fission and nu_fission must have identical "
            f"positive group support; mismatch in group(s) {rendered_groups}"
        )

    # The ratio is the group-wise effective neutron yield.  Record it only
    # where both reaction vectors are positive; its magnitude is not subject
    # to a universal isotope- and energy-independent acceptance interval.
    mask = fission_positive & nu_fission_positive
    if not np.any(mask):
        return
    ratio = nu_fission[mask] / fission[mask]
    groups = np.nonzero(mask)[0]
    finite = np.isfinite(ratio) & (ratio > 0.0)
    if not np.all(finite):
        invalid_groups = groups[~finite]
        rendered_groups = ", ".join(str(int(index) + 1) for index in invalid_groups)
        report.nu_ratio_errors.append(
            f"mixture {label}: nu_fission/fission must be finite and positive "
            f"on common positive support; invalid group(s) {rendered_groups}"
        )
    ratio = ratio[finite]
    groups = groups[finite]
    if ratio.size == 0:
        return
    report.nu_ratio_checked_bins += int(ratio.size)
    local_min = float(np.min(ratio))
    local_max = float(np.max(ratio))
    report.nu_ratio_min = (
        local_min if report.nu_ratio_min is None else min(report.nu_ratio_min, local_min)
    )
    report.nu_ratio_max = (
        local_max if report.nu_ratio_max is None else max(report.nu_ratio_max, local_max)
    )
    low = np.zeros(ratio.shape, dtype=bool)
    high = np.zeros(ratio.shape, dtype=bool)
    if minimum is not None:
        low = ratio < minimum
    if maximum is not None:
        high = ratio > maximum
    if not np.any(low | high):
        return
    distance = np.zeros(ratio.shape, dtype=float)
    if minimum is not None:
        distance = np.maximum(distance, minimum - ratio)
    if maximum is not None:
        distance = np.maximum(distance, ratio - maximum)
    index = int(np.argmax(np.where(low | high, distance, -np.inf)))
    group_index = int(groups[index])
    value = float(ratio[index])
    report.nu_ratio_worst = f"{label}: group={group_index + 1} nu={value:.6e}"
    lower = "-inf" if minimum is None else f"{minimum:.6e}"
    upper = "+inf" if maximum is None else f"{maximum:.6e}"
    report.nu_ratio_warnings.append(
        f"mixture {label}: nu_fission/fission={value:.6e} in group "
        f"{group_index + 1} is outside explicitly configured "
        f"[{lower}, {upper}]"
    )


def _check_transport_p1(
    report: _MutablePhysicsReport,
    *,
    label: str,
    total: np.ndarray | None,
    transport_total: np.ndarray | None,
    scatter: np.ndarray | None,
    axes: str | None,
    energy_groups: int,
    legendre_order: int,
    reference_flux: np.ndarray | None,
    sph_applied: Any,
    sph_apply_operator: Any,
    applied_sph: np.ndarray | None,
    threshold: float,
) -> None:
    if total is None or transport_total is None or scatter is None:
        return
    p1 = scatter_moment_matrix(scatter, axes, energy_groups, legendre_order, moment=1)
    if p1 is None:
        return
    try:
        has_applied_sph = _scatter_multiplicity_bool(sph_applied)
    except (TypeError, ValueError, OverflowError):
        report.transport_p1_errors.append(
            f"{label}: transport/P1 validation requires an explicit boolean sph_applied"
        )
        return
    if has_applied_sph:
        if _attr_text(sph_apply_operator) != "divide-xs-by-nsph":
            report.transport_p1_errors.append(
                f"{label}: transport/P1 validation cannot interpret SPH-applied "
                "cross sections without sph_apply_operator='divide-xs-by-nsph'"
            )
            return
        if (
            applied_sph is None
            or not np.all(np.isfinite(applied_sph))
            or np.any(applied_sph <= 0.0)
        ):
            report.transport_p1_errors.append(
                f"{label}: transport/P1 validation of SPH-applied cross sections "
                f"requires a finite positive applied_sph vector with {energy_groups} values"
            )
            return
        if reference_flux is not None:
            # The supported SPH operator divides total, TransportXS, and
            # incoming P1 rows by NSPH. Its TransportXS identity therefore
            # uses phi_check = NSPH * phi_CE: the incoming factors cancel
            # while the outgoing denominator supplies the same NSPH divisor
            # as the stored transport vector. This algebraic audit leaves
            # the frozen CE reference flux and all XS datasets unchanged.
            with np.errstate(over="ignore", invalid="ignore"):
                reference_flux = reference_flux * applied_sph
            if not np.all(np.isfinite(reference_flux)) or np.any(reference_flux <= 0.0):
                report.transport_p1_errors.append(
                    f"{label}: applied_sph times reference flux is non-positive or non-finite "
                    "during transport/P1 validation"
                )
                return
    # OpenMC's transport correction is an outgoing-group quantity.  It needs
    # the incoming-group flux used to tally P1; a bare row sum is not the
    # TransportXS definition and can produce large false failures.  A
    # root-level flux is deliberately not attached to stateful calculations,
    # because it is file-global rather than state-bound.
    if reference_flux is None:
        report.transport_p1_skipped += 1
        return
    if not (
        np.all(np.isfinite(total))
        and np.all(np.isfinite(transport_total))
        and np.all(np.isfinite(p1))
        and np.all(np.isfinite(reference_flux))
        and np.all(reference_flux > 0.0)
    ):
        report.transport_p1_skipped += 1
        return
    report.transport_p1_checked += 1
    correction = np.sum(reference_flux[:, np.newaxis] * p1, axis=0) / reference_flux
    derived = total - correction
    residual = transport_total - derived
    _, max_abs, max_rel, index = _relative_worst(residual, transport_total)
    _update_worst(
        current=report,
        attr_rel="transport_p1_max_rel",
        attr_abs="transport_p1_max_abs",
        attr_worst="transport_p1_worst",
        max_rel=max_rel,
        max_abs=max_abs,
        worst=(
            f"{label}: group={index + 1} transport_total={transport_total[index]:.6e} "
            f"p1_derived={derived[index]:.6e}"
        ),
    )
    if max_rel > threshold:
        report.transport_p1_errors.append(
            "flux-weighted transport_total/P1 max relative residual "
            f"{max_rel:.6e} (abs {max_abs:.6e}) at {label}: group={index + 1} "
            f"exceeds fail threshold {threshold:.6e}"
        )


def scatter_moment_matrix(
    values: np.ndarray,
    axes: str | None,
    ngroups: int,
    legendre_order: int,
    *,
    moment: int,
) -> np.ndarray | None:
    expected_moments = legendre_order + 1
    if moment >= expected_moments or values.ndim != 3:
        return None
    normalized = normalize_axes(axes)
    moment_first = values.shape == (expected_moments, ngroups, ngroups)
    moment_last = values.shape == (ngroups, ngroups, expected_moments)
    if normalized in MOMENT_FIRST_SCATTER_AXES and moment_first:
        return np.asarray(values[moment])
    if normalized in MOMENT_LAST_SCATTER_AXES and moment_last:
        return np.asarray(values[:, :, moment])
    if axes is not None:
        return None
    if moment_first and not moment_last:
        return np.asarray(values[moment])
    if moment_last and not moment_first:
        return np.asarray(values[:, :, moment])
    return None


def _finalize_adf_faces(
    report: _MutablePhysicsReport,
    adf_names_by_calc: list[tuple[str, ...]],
) -> None:
    if not adf_names_by_calc:
        return
    present = [names for names in adf_names_by_calc if names]
    report.adf_calculations = len(present)
    if not present:
        return
    first = present[0]
    report.adf_faces = first
    for index, names in enumerate(adf_names_by_calc, start=1):
        if names != first:
            report.adf_face_errors.append(
                f"calculation index {index}: ADF faces {names!r} do not match {first!r}"
            )
            return


def _iter_calculations(
    h5: Any, mixture_names: tuple[str, ...]
) -> Iterator[tuple[str, str, Any, Any]]:
    mixtures = h5["mixtures"]
    for mixture_name in mixture_names:
        mixture = mixtures[mixture_name]
        if "states" in mixture:
            states = mixture["states"]
            for state_name in _sorted_state_names(states):
                yield (
                    mixture_name,
                    f"{mixture_name}/states/{state_name}",
                    states[state_name],
                    mixture,
                )
        else:
            yield mixture_name, mixture_name, mixture, None


def _root_reference_flux_by_mixture(
    h5: Any,
    *,
    mixture_names: tuple[str, ...],
    energy_groups: int,
) -> dict[str, np.ndarray]:
    """Return canonical direct-mixture OpenMC flux rows, or no rows.

    The transport/P1 identity is meaningful only when the flux row is bound to
    the same mixture ordering and the same high-to-low group convention as the
    MGXS arrays.  Fail closed on incomplete metadata instead of guessing.
    """

    if "openmc_volume_flux" not in h5:
        return {}
    dataset = h5["openmc_volume_flux"]
    try:
        values = np.asarray(dataset[:], dtype=float)
    except (TypeError, ValueError, OSError):
        return {}
    if values.shape != (len(mixture_names), energy_groups):
        return {}
    if "mixture_names" not in dataset.attrs:
        return {}
    try:
        declared_names = decode_hdf5_names(dataset.attrs["mixture_names"])
    except (TypeError, ValueError, UnicodeDecodeError):
        return {}
    if declared_names != mixture_names:
        return {}
    if _attr_text(dataset.attrs.get("group_order")) != MGXS_DONJON_GROUP_ORDER:
        return {}
    if not np.all(np.isfinite(values)) or np.any(values <= 0.0):
        return {}
    return {
        mixture_name: np.asarray(values[index], dtype=float)
        for index, mixture_name in enumerate(mixture_names)
    }


def _sorted_state_names(states_group: Any) -> list[str]:
    def key(name: str) -> tuple[int, int | str]:
        try:
            return (0, int(name))
        except ValueError:
            return (1, name)

    return sorted(states_group.keys(), key=key)


def _vector_or_none(group: Any, name: str, ngroups: int) -> np.ndarray | None:
    if name not in group:
        return None
    try:
        values = np.asarray(group[name][:], dtype=float).reshape(-1)
    except (TypeError, ValueError, OSError):
        return None
    if values.shape != (ngroups,):
        return None
    return values


def _scatter_or_none(group: Any, name: str) -> np.ndarray | None:
    if name not in group:
        return None
    try:
        return np.asarray(group[name][:], dtype=float)
    except (TypeError, ValueError, OSError):
        return None


def _scatter_axes(group: Any, h5: Any, parent_group: Any | None) -> str | None:
    sources = [group.attrs]
    if parent_group is not None:
        sources.append(parent_group.attrs)
    sources.append(h5.attrs)
    for source in sources:
        for key in ("scatter_axes", "axes"):
            if key in source:
                return _attr_text(source[key])
    return None


def _scatter_balance_contract(
    group: Any,
    h5: Any,
    parent_group: Any | None,
    *,
    transport_total_present: bool,
) -> _ResolvedScatterContract:
    mgxs_type, issues = _consistent_inherited_scatter_attr(
        group,
        h5,
        parent_group,
        "openmc_scatter_mgxs_type",
        _canonical_scatter_mgxs_type,
    )
    weighted, weighted_issues = _consistent_inherited_scatter_attr(
        group,
        h5,
        parent_group,
        "openmc_scatter_multiplicity_weighted",
        _scatter_multiplicity_bool,
    )
    balance, balance_issues = _consistent_inherited_scatter_attr(
        group,
        h5,
        parent_group,
        "openmc_scatter_balance_dataset",
        _scatter_balance_dataset,
    )
    issues.extend(weighted_issues)
    issues.extend(balance_issues)
    transport, transport_issues = _consistent_inherited_scatter_attr(
        group,
        h5,
        parent_group,
        "openmc_transport_mgxs_type",
        _canonical_transport_mgxs_type,
    )
    issues.extend(transport_issues)
    transport_declared = transport is not None
    declared = mgxs_type is not None or weighted is not None or balance is not None

    if mgxs_type is not None:
        expected_weighted = (
            _normalize_openmc_mgxs_type(mgxs_type)
            in NU_WEIGHTED_OPENMC_SCATTER_MGXS_TYPES
        )
        expected_balance = (
            "reduced_absorption" if expected_weighted else "absorption"
        )
        if weighted is not None and weighted != expected_weighted:
            issues.append(
                "openmc_scatter_multiplicity_weighted contradicts "
                f"openmc_scatter_mgxs_type {mgxs_type!r}"
            )
        if balance is not None and balance != expected_balance:
            issues.append(
                "openmc_scatter_balance_dataset contradicts "
                f"openmc_scatter_mgxs_type {mgxs_type!r}; expected "
                f"{expected_balance!r}"
            )
        weighted = expected_weighted
        balance = expected_balance
    else:
        # Files written before estimator metadata existed used ordinary,
        # non-multiplicity-weighted scattering. Preserve that legacy meaning.
        if weighted is None and balance is None:
            weighted = False
            balance = "absorption"
        elif weighted is None:
            weighted = balance == "reduced_absorption"
        elif balance is None:
            balance = "reduced_absorption" if weighted else "absorption"
        elif weighted != (balance == "reduced_absorption"):
            issues.append(
                "openmc_scatter_multiplicity_weighted contradicts "
                f"openmc_scatter_balance_dataset {balance!r}"
            )

    expected_transport = "nu-transport" if weighted else "transport"
    if transport is not None:
        if transport != expected_transport:
            issues.append(
                f"openmc_transport_mgxs_type {transport!r} contradicts the "
                f"resolved scattering convention; expected {expected_transport!r}"
            )
        if not transport_total_present:
            issues.append(
                f"openmc_transport_mgxs_type {transport!r} is declared but "
                "transport_total is missing"
            )
    elif transport_total_present:
        if weighted:
            issues.append(
                "nu-weighted scattering with transport_total requires an explicit "
                "openmc_transport_mgxs_type='nu-transport'; an older ordinary "
                "TransportXS cannot be inferred as multiplicity-weighted"
            )
        else:
            # Legacy ordinary-scatter handoffs predate transport estimator
            # metadata. Their transport_total can only have come from the
            # ordinary OpenMC TransportXS supported at that time.
            transport = "transport"

    return _ResolvedScatterContract(
        mgxs_type=mgxs_type,
        multiplicity_weighted=bool(weighted),
        balance_dataset=str(balance),
        declared=declared,
        transport_mgxs_type=transport,
        transport_declared=transport_declared,
        issue="; ".join(issues) or None,
    )


def _consistent_inherited_scatter_attr(
    group: Any,
    h5: Any,
    parent_group: Any | None,
    name: str,
    parser: Any,
) -> tuple[Any | None, list[str]]:
    sources = [("calculation", group.attrs)]
    if parent_group is not None:
        sources.append(("mixture", parent_group.attrs))
    sources.append(("root", h5.attrs))
    declarations: list[tuple[str, Any]] = []
    issues: list[str] = []
    for scope, attrs in sources:
        if name not in attrs:
            continue
        try:
            declarations.append((scope, parser(attrs[name])))
        except (TypeError, ValueError, OverflowError) as exc:
            issues.append(f"invalid {scope} {name}: {exc}")
    if not declarations:
        return None, issues
    first = declarations[0][1]
    if any(value != first for _scope, value in declarations[1:]):
        detail = ", ".join(
            f"{scope}={value!r}" for scope, value in declarations
        )
        issues.append(f"contradictory inherited {name} declarations: {detail}")
    return first, issues


def _canonical_scatter_mgxs_type(value: Any) -> str:
    normalized = _normalize_openmc_mgxs_type(_attr_text(value))
    if normalized in ORDINARY_OPENMC_SCATTER_MGXS_TYPES:
        return normalized
    if normalized in NU_WEIGHTED_OPENMC_SCATTER_MGXS_TYPES:
        return normalized.replace("nu scatter", "nu-scatter")
    raise ValueError(
        f"unsupported value {_attr_text(value)!r}; expected ordinary "
        "'scatter matrix' / 'consistent scatter matrix' or nu-weighted "
        "'nu-scatter matrix' / 'consistent nu-scatter matrix'"
    )


def _canonical_transport_mgxs_type(value: Any) -> str:
    normalized = _normalize_openmc_mgxs_type(_attr_text(value))
    if normalized not in OPENMC_TRANSPORT_MGXS_TYPES:
        raise ValueError(
            f"unsupported value {_attr_text(value)!r}; expected "
            "'transport' or 'nu-transport'"
        )
    return "nu-transport" if normalized == "nu transport" else "transport"


def _scatter_multiplicity_bool(value: Any) -> bool:
    array = np.asarray(value)
    if array.size != 1:
        raise ValueError("must be a scalar boolean")
    scalar = array.reshape(-1)[0]
    if isinstance(scalar, (bool, np.bool_)):
        return bool(scalar)
    if isinstance(scalar, (int, np.integer)) and scalar in (0, 1):
        return bool(scalar)
    raise ValueError("must be boolean true/false (or integer 1/0)")


def _scatter_balance_dataset(value: Any) -> str:
    balance = _attr_text(value).strip()
    if balance not in OPENMC_SCATTER_BALANCE_DATASETS:
        raise ValueError("must be 'absorption' or 'reduced_absorption'")
    return balance


def _normalize_openmc_mgxs_type(value: str) -> str:
    return " ".join(value.strip().lower().replace("_", " ").replace("-", " ").split())


def _adf_names(group: Any) -> tuple[str, ...]:
    for name in ("adf", "ADF", "discontinuity_factors"):
        if name not in group:
            continue
        obj = group[name]
        if hasattr(obj, "keys"):
            return tuple(str(face_name) for face_name in obj)
        values = np.asarray(obj[:], dtype=float)
        for key in ("names", "face_names", "adf_names"):
            if key in obj.attrs:
                return tuple(
                    _attr_text(value)
                    for value in np.asarray(obj.attrs[key]).reshape(-1)
                )
        if values.ndim == 1:
            return ("FD_B",)
        if values.ndim == 2:
            return tuple(f"FD_{index + 1:05d}" for index in range(values.shape[0]))
    return ()


def _attr_with_parent(
    attrs: Any,
    parent_attrs: Any | None,
    name: str,
    default: object,
) -> object:
    value = attrs.get(name)
    if value is None and parent_attrs is not None:
        value = parent_attrs.get(name)
    return default if value is None else value


def _attr_text(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, np.bytes_):
        return value.decode("utf-8")
    return str(value)


def _relative_worst(
    residual: np.ndarray,
    denominator: np.ndarray,
) -> tuple[np.ndarray, float, float, int]:
    abs_residual = np.abs(residual)
    relative = abs_residual / np.maximum(np.abs(denominator), 1.0e-30)
    index = int(np.argmax(relative))
    return relative, float(abs_residual[index]), float(relative[index]), index


def _update_worst(
    *,
    current: _MutablePhysicsReport,
    attr_rel: str,
    attr_abs: str,
    attr_worst: str,
    max_rel: float,
    max_abs: float,
    worst: str,
) -> None:
    current_rel = getattr(current, attr_rel)
    if current_rel is None or max_rel > current_rel:
        setattr(current, attr_rel, max_rel)
        setattr(current, attr_abs, max_abs)
        setattr(current, attr_worst, worst)
