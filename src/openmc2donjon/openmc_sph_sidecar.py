"""Build OpenMC-side SPH sidecars from CE/MG flux comparisons."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any

from . import __version__
from .openmc_provenance import file_sha256
from .sph_augment import SphSidecarReport, create_table_sph_sidecar
from .sph_iteration import SphUpdateTableReport, create_sph_update_table


SCHEMA = "openmc2donjon.openmc-sph-sidecar.v1"
PASS_DECISION = "openmc2donjon_openmc_sph_sidecar_passed"
SOURCE_BINDING_SCHEMA = "openmc2donjon.openmc-sph-source-bindings.v1"


@dataclass(frozen=True)
class OpenmcSphSidecarReport:
    input_h5: Path
    output_h5: Path
    output_table: Path
    reference_flux: str | Path
    mg_flux: str | Path
    update: SphUpdateTableReport
    sidecar: SphSidecarReport
    source_bindings: dict[str, Any]


def create_openmc_sph_sidecar(
    input_h5: Path,
    output_h5: Path,
    *,
    reference_flux: str | Path,
    mg_flux: str | Path,
    table_output: Path | None = None,
    previous_sph: str | Path | None = None,
    damping: float = 1.0,
    clip_min: float | None = None,
    clip_max: float | None = None,
    flux_normalization: str = "auto",
    sph_target: str = "rate",
    zero_flux_policy: str = "reject",
    flux_floor_rel: float | None = None,
    freeze_groups: tuple[int, ...] | None = None,
    tie_mixture_groups: tuple[tuple[str, ...], ...] | None = None,
    require_reference_flux_std_dev: bool = False,
    max_reference_flux_std_dev_rel: float | None = None,
    require_mg_flux_std_dev: bool = False,
    max_mg_flux_std_dev_rel: float | None = None,
    source_label: str = "openmc-ce-mg-sph",
    sph_kind: str = "openmc-ce-mg",
    sph_real: bool = True,
    sph_applied: bool = False,
    force: bool = False,
    summary_json: Path | None = None,
) -> OpenmcSphSidecarReport:
    """Compute OpenMC CE/MG SPH factors and write a sidecar HDF5.

    ``reference_flux`` is the OpenMC continuous-energy reference flux and
    ``mg_flux`` is the OpenMC multi-group macro flux.  The CE calculation uses
    the detailed fine geometry and the MG calculation uses the homogenized
    coarse geometry.  Their comparison-domain ordering, energy-group
    structure, physical state, and boundary conditions must be aligned with
    the declared ``input_h5`` handoff.
    The default target preserves reaction rates and the default normalization
    resolves group-wise H-FACTOR/kappa-fission data to power normalization.
    ``sph_target="flux"`` and non-power normalizations are retained for
    explicit diagnostic studies.
    """

    input_h5 = Path(input_h5)
    output_h5 = Path(output_h5)
    output_table = (
        Path(table_output)
        if table_output is not None
        else output_h5.with_suffix(".sph.csv")
    )
    if output_table == output_h5:
        raise ValueError("--table-output must be different from --output")
    _reject_source_output_aliases(
        input_h5=input_h5,
        reference_flux=reference_flux,
        mg_flux=mg_flux,
        previous_sph=previous_sph,
        output_h5=output_h5,
        output_table=output_table,
    )
    if output_h5.exists() and not force:
        raise FileExistsError(f"output already exists; use --force to overwrite: {output_h5}")

    update = create_sph_update_table(
        input_h5,
        output_table,
        reference_flux=reference_flux,
        low_order_flux=mg_flux,
        previous_sph=previous_sph,
        damping=damping,
        clip_min=clip_min,
        clip_max=clip_max,
        flux_normalization=flux_normalization,
        sph_target=sph_target,
        zero_flux_policy=zero_flux_policy,
        flux_floor_rel=flux_floor_rel,
        freeze_groups=freeze_groups,
        tie_mixture_groups=tie_mixture_groups,
        require_reference_flux_std_dev=require_reference_flux_std_dev,
        max_reference_flux_std_dev_rel=max_reference_flux_std_dev_rel,
        require_low_order_flux_std_dev=require_mg_flux_std_dev,
        max_low_order_flux_std_dev_rel=max_mg_flux_std_dev_rel,
        source_label=source_label,
        force=force,
        summary_json=None,
    )
    source_bindings = _source_bindings(
        update,
        require_reference_flux_std_dev=require_reference_flux_std_dev,
        max_reference_flux_std_dev_rel=max_reference_flux_std_dev_rel,
        require_mg_flux_std_dev=require_mg_flux_std_dev,
        max_mg_flux_std_dev_rel=max_mg_flux_std_dev_rel,
    )
    sidecar = create_table_sph_sidecar(
        input_h5,
        output_h5,
        table=output_table,
        force=force,
        sph_kind=sph_kind,
        sph_real=sph_real,
        sph_applied=sph_applied,
        summary_json=None,
    )
    source_bindings = _write_physics_provenance(
        output_h5,
        update,
        source_bindings=source_bindings,
    )
    report = OpenmcSphSidecarReport(
        input_h5=input_h5,
        output_h5=output_h5,
        output_table=output_table,
        reference_flux=reference_flux,
        mg_flux=mg_flux,
        update=update,
        sidecar=sidecar,
        source_bindings=source_bindings,
    )
    print_report(report)
    if summary_json is not None:
        write_summary(summary_json, report)
    return report


def _write_physics_provenance(
    path: Path,
    update: SphUpdateTableReport,
    *,
    source_bindings: dict[str, Any],
) -> dict[str, Any]:
    """Persist the derivation and convergence evidence beside the factors."""

    import h5py

    raw_minimum = float(update.raw_update_minimum)
    raw_maximum = float(update.raw_update_maximum)
    max_residual = max(abs(raw_minimum - 1.0), abs(raw_maximum - 1.0))
    derivation = (
        "rate-preserving-ce-mg-fixed-point"
        if update.sph_target == "rate"
        else "ce-mg-flux-fixed-point"
    )
    with h5py.File(path, "r+") as h5:
        h5.attrs["sph_derivation"] = derivation
        h5.attrs["sph_target"] = update.sph_target
        h5.attrs["sph_raw_update_minimum"] = raw_minimum
        h5.attrs["sph_raw_update_maximum"] = raw_maximum
        h5.attrs["sph_max_update_residual"] = max_residual
        h5.attrs["sph_flux_normalization"] = update.flux_normalization
        h5.attrs["sph_zero_flux_policy"] = update.zero_flux_policy
        h5.attrs["sph_identity_bin_count"] = update.identity_bin_count
        h5.attrs["sph_floored_bin_count"] = update.floored_bin_count
        h5.attrs["sph_frozen_group_bin_count"] = update.frozen_group_bin_count
        h5.attrs["sph_tie_mixture_groups"] = json.dumps(
            [list(group) for group in update.tie_mixture_groups],
            separators=(",", ":"),
        )
        h5.attrs["sph_tied_bin_count"] = update.tied_bin_count
        h5.attrs["sph_clipped_count"] = update.clipped_count
        for name, value in source_bindings.items():
            if value is not None:
                h5.attrs[name] = value
    return source_bindings


def _reject_source_output_aliases(
    *,
    input_h5: Path,
    reference_flux: str | Path,
    mg_flux: str | Path,
    previous_sph: str | Path | None,
    output_h5: Path,
    output_table: Path,
) -> None:
    sources = {
        "input HDF5": input_h5,
        "reference flux": _source_file(reference_flux),
        "MG flux": _source_file(mg_flux),
    }
    if previous_sph is not None:
        sources["previous SPH"] = _source_file(previous_sph)
    for label, source in sources.items():
        for output in (output_h5, output_table):
            if source.expanduser().resolve() == output.expanduser().resolve():
                raise ValueError(
                    f"{label} source must be different from output path {output}"
                )


def _source_file(source: str | Path) -> Path:
    text = str(source)
    path_text, separator, _dataset = text.partition("::")
    return Path(path_text if separator else text)


def _source_bindings(
    update: SphUpdateTableReport,
    *,
    require_reference_flux_std_dev: bool,
    max_reference_flux_std_dev_rel: float | None,
    require_mg_flux_std_dev: bool,
    max_mg_flux_std_dev_rel: float | None,
) -> dict[str, Any]:
    input_h5 = update.input_h5.expanduser().resolve()
    reference_flux = update.reference_flux_source.expanduser().resolve()
    mg_flux = update.low_order_flux_source.expanduser().resolve()
    previous_sph = (
        None
        if update.previous_sph_source is None
        else update.previous_sph_source.expanduser().resolve()
    )
    bindings = {
        "sph_source_binding_schema": SOURCE_BINDING_SCHEMA,
        "sph_input_h5_path": str(input_h5),
        "sph_input_h5_sha256": file_sha256(input_h5),
        "sph_reference_flux_path": str(reference_flux),
        "sph_reference_flux_sha256": file_sha256(reference_flux),
        "sph_reference_flux_dataset": update.reference_flux_dataset,
        "sph_reference_flux_layout_verified": _flux_layout_verified(
            reference_flux,
            update.reference_flux_dataset,
        ),
        "sph_mg_flux_path": str(mg_flux),
        "sph_mg_flux_sha256": file_sha256(mg_flux),
        "sph_mg_flux_dataset": update.low_order_flux_dataset,
        "sph_mg_flux_layout_verified": _flux_layout_verified(
            mg_flux,
            update.low_order_flux_dataset,
        ),
        "sph_previous_sph_used": previous_sph is not None,
        "sph_previous_sph_path": (
            None if previous_sph is None else str(previous_sph)
        ),
        "sph_previous_sph_sha256": (
            None if previous_sph is None else file_sha256(previous_sph)
        ),
        "sph_previous_sph_dataset": update.previous_sph_dataset,
    }
    bindings.update(
        _uncertainty_bindings(
            prefix="sph_reference_flux",
            require_coverage=require_reference_flux_std_dev,
            limit=max_reference_flux_std_dev_rel,
            observed=update.reference_flux_max_relative_std_dev,
            std_dev_dataset=update.reference_flux_std_dev_dataset,
        )
    )
    bindings.update(
        _uncertainty_bindings(
            prefix="sph_mg_flux",
            require_coverage=require_mg_flux_std_dev,
            limit=max_mg_flux_std_dev_rel,
            observed=update.low_order_flux_max_relative_std_dev,
            std_dev_dataset=update.low_order_flux_std_dev_dataset,
        )
    )
    return bindings


def _uncertainty_bindings(
    *,
    prefix: str,
    require_coverage: bool,
    limit: float | None,
    observed: float | None,
    std_dev_dataset: str | None,
) -> dict[str, Any]:
    coverage = bool(std_dev_dataset) and observed is not None
    passed = (
        bool(require_coverage)
        and coverage
        and limit is not None
        and observed is not None
        and _finite_nonnegative(limit)
        and _finite_nonnegative(observed)
        and observed <= limit
    )
    return {
        f"{prefix}_uncertainty_require_coverage": bool(require_coverage),
        f"{prefix}_uncertainty_coverage": coverage,
        f"{prefix}_uncertainty_limit": limit,
        f"{prefix}_uncertainty_observed_max_rel": observed,
        f"{prefix}_uncertainty_pass": passed,
        f"{prefix}_std_dev_dataset": std_dev_dataset,
    }


def _finite_nonnegative(value: float) -> bool:
    import math

    return math.isfinite(float(value)) and float(value) >= 0.0


def _flux_layout_verified(path: Path, dataset_path: str | None) -> bool:
    """Return true only for an explicit HDF5 tally-layout attestation."""

    if dataset_path is None:
        return False
    try:
        import h5py

        if not h5py.is_hdf5(path):
            return False
        with h5py.File(path, "r") as h5:
            if dataset_path not in h5:
                return False
            dataset = h5[dataset_path]
            if not hasattr(dataset, "attrs"):
                return False
            return (
                _boolean_value(dataset.attrs.get("energy_bounds_verified")) is True
                and _boolean_value(
                    dataset.attrs.get("spatial_domain_order_verified")
                )
                is True
            )
    except OSError:
        return False


def _boolean_value(value: Any) -> bool | None:
    if isinstance(value, bool):
        return value
    if isinstance(value, int) and value in (0, 1):
        return bool(value)
    if isinstance(value, bytes):
        value = value.decode("utf-8", errors="replace")
    text = "" if value is None else str(value).strip().lower()
    if text in {"true", "1"}:
        return True
    if text in {"false", "0"}:
        return False
    return None


def print_report(report: OpenmcSphSidecarReport) -> None:
    print("OpenMC CE/MG SPH sidecar")
    print(f"  schema: {SCHEMA}")
    print(f"  input: {report.input_h5}")
    print(f"  output: {report.output_h5}")
    print(f"  table: {report.output_table}")
    print(f"  reference_flux: {report.reference_flux}")
    print(f"  mg_flux: {report.mg_flux}")
    print(
        f"  mixtures={len(report.sidecar.mixture_names)} "
        f"groups={report.sidecar.energy_groups} "
        f"sph_range={report.sidecar.sph_min:g}..{report.sidecar.sph_max:g} "
        f"clipped={report.update.clipped_count}"
    )
    print()
    print("OpenMC SPH sidecar decision")
    print(f"  {PASS_DECISION}")


def write_summary(path: Path, report: OpenmcSphSidecarReport) -> None:
    payload = {
        "schema": SCHEMA,
        "package_version": __version__,
        "decision": PASS_DECISION,
        "input_h5": str(report.input_h5),
        "output_h5": str(report.output_h5),
        "output_table": str(report.output_table),
        "reference_flux": str(report.reference_flux),
        "reference_flux_dataset": report.update.reference_flux_dataset,
        "reference_flux_std_dev_dataset": report.update.reference_flux_std_dev_dataset,
        "reference_flux_max_relative_std_dev": (
            report.update.reference_flux_max_relative_std_dev
        ),
        "mg_flux": str(report.mg_flux),
        "mg_flux_dataset": report.update.low_order_flux_dataset,
        "mg_flux_std_dev_dataset": report.update.low_order_flux_std_dev_dataset,
        "mg_flux_max_relative_std_dev": report.update.low_order_flux_max_relative_std_dev,
        "previous_sph": None
        if report.update.previous_sph_source is None
        else str(report.update.previous_sph_source),
        "previous_sph_dataset": report.update.previous_sph_dataset,
        "mixture_count": len(report.sidecar.mixture_names),
        "mixture_names": list(report.sidecar.mixture_names),
        "energy_groups": report.sidecar.energy_groups,
        "damping": report.update.damping,
        "clip_min": report.update.clip_min,
        "clip_max": report.update.clip_max,
        "flux_normalization": report.update.flux_normalization,
        "sph_target": report.update.sph_target,
        "zero_flux_policy": report.update.zero_flux_policy,
        "identity_bin_count": report.update.identity_bin_count,
        "flux_floor_rel": report.update.flux_floor_rel,
        "floored_bin_count": report.update.floored_bin_count,
        "freeze_groups": None
        if report.update.freeze_groups is None
        else list(report.update.freeze_groups),
        "frozen_group_bin_count": report.update.frozen_group_bin_count,
        "tie_mixture_groups": [list(group) for group in report.update.tie_mixture_groups],
        "tied_bin_count": report.update.tied_bin_count,
        "normalization_factor": report.update.normalization_factor,
        "sph_min": report.sidecar.sph_min,
        "sph_max": report.sidecar.sph_max,
        "sph_kind": report.sidecar.sph_kind,
        "sph_real": report.sidecar.sph_real,
        "sph_applied": report.sidecar.sph_applied,
        "raw_update_minimum": report.update.raw_update_minimum,
        "raw_update_maximum": report.update.raw_update_maximum,
        "clipped_count": report.update.clipped_count,
        "source_label": report.update.source_label,
        "source_bindings": report.source_bindings,
        "formula": (
            "sph = previous_sph * "
            "(normalized_openmc_mg_flux / (previous_sph * openmc_ce_reference_flux)) "
            "** damping"
            if report.update.sph_target == "rate"
            else "sph = previous_sph * "
            "(openmc_ce_reference_flux / normalized_openmc_mg_flux) ** damping"
        ),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
