"""Strict physical provenance gates for production SPH handoffs.

The base contract is geometry-agnostic: one assembly, an assembly containing
many homogenization domains, an arbitrary colorset, and a full-core coarse
mesh are all valid shapes.  Benchmark/template topology belongs in a thin
specialization, never in the physical SPH contract itself.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np

from .hdf5_names import read_mixture_names
from .openmc_provenance import file_sha256


PHYSICAL_SPH_MAX_UPDATE_RESIDUAL = 0.02
EXPECTED_COLORSET_DOMAINS = 7
SPH_SOURCE_BINDING_SCHEMA = "openmc2donjon.openmc-sph-source-bindings.v1"
SPH_APPLY_BINDING_SCHEMA = "openmc2donjon.sph-apply-bindings.v1"
_ACTIVE_SPH_NAMES = frozenset({"sph", "nsph"})
_ADF_NAMES = frozenset({"adf", "discontinuity_factors"})


def physical_sph_issues(
    path: str | Path,
    *,
    max_update_residual: float = PHYSICAL_SPH_MAX_UPDATE_RESIDUAL,
) -> list[str]:
    """Return every reason an applied SPH handoff is not production-ready.

    This contract deliberately rejects empirical/global calibration provenance.
    The accepted route pairs a detailed OpenMC CE reference geometry with an
    OpenMC MG homogenized coarse geometry through an explicit comparison-domain
    mapping.  Its rate-preserving fixed point has already been applied to the
    declared coarse domains.  It intentionally imposes no component count,
    domain count, lattice, or core topology.
    """

    import h5py

    issues: list[str] = []
    with h5py.File(Path(path), "r") as h5:
        attrs = {str(name): value for name, value in h5.attrs.items()}
        try:
            mixture_names = read_mixture_names(h5)
        except ValueError as exc:
            issues.append(str(exc))
            mixture_names = ()

        if "mixture_names" not in h5:
            issues.append("/mixture_names must declare the SPH domain order")
        if not mixture_names:
            issues.append("physical SPH requires at least one declared domain")
        energy_groups = _energy_groups(h5, issues)
        _validate_applied_sph_payload(
            h5,
            mixture_names=mixture_names,
            energy_groups=energy_groups,
            issues=issues,
        )
        _reject_active_equivalence_payloads(h5, issues)
        filled_bins = sum(
            _attribute_item_count(
                h5["mixtures"][name].attrs.get("zero_flux_filled_groups")
            )
            for name in mixture_names
            if "mixtures" in h5 and name in h5["mixtures"]
        )

        applied = _boolean_attr(attrs.get("sph_applied"))
        source = _text_attr(attrs.get("sph_applied_source"))
        operator = _text_attr(attrs.get("sph_apply_operator"))
        kind = _text_attr(attrs.get("sph_kind")).lower()
        derivation = _text_attr(attrs.get("sph_derivation"))
        target = _text_attr(attrs.get("sph_target"))
        normalization = _text_attr(attrs.get("sph_flux_normalization"))
        zero_flux_policy = _text_attr(attrs.get("sph_zero_flux_policy"))
        identity_bin_count = _integer_attr(attrs.get("sph_identity_bin_count"))
        floored_bin_count = _integer_attr(attrs.get("sph_floored_bin_count"))
        frozen_bin_count = _integer_attr(attrs.get("sph_frozen_group_bin_count"))
        clipped_count = _integer_attr(attrs.get("sph_clipped_count"))
        is_real = _boolean_attr(attrs.get("sph_real"))
        residual = attrs.get("sph_max_update_residual")

    if applied is not True:
        issues.append("sph_applied=true is required; run apply-sph before Converter")
    if not source:
        issues.append("sph_applied_source must record the physical SPH sidecar")
    if operator != "divide-xs-by-nsph":
        issues.append("sph_apply_operator must be divide-xs-by-nsph")
    if is_real is not True or not kind.startswith("openmc-ce-mg"):
        issues.append(
            "SPH must be derived from a real OpenMC CE fine-to-MG coarse "
            "comparison-domain calculation"
        )
    if any(token in kind for token in ("global", "optical", "calibrat", "empirical")):
        issues.append("empirical, global, optical, or calibrated SPH provenance is forbidden")
    if derivation != "rate-preserving-ce-mg-fixed-point" or target != "rate":
        issues.append("SPH must use the rate-preserving CE/MG fixed-point derivation")
    if normalization != "power":
        issues.append("SPH must use H-FACTOR/kappa-fission power normalization")
    if zero_flux_policy != "reject":
        issues.append("SPH zero-flux policy must be reject; identity is forbidden")
    if identity_bin_count != 0:
        issues.append("SPH provenance must report zero identity-substituted bins")
    if floored_bin_count != 0:
        issues.append("SPH provenance must report zero flux-floored bins")
    if frozen_bin_count != 0:
        issues.append("SPH provenance must report zero frozen-group bins")
    if clipped_count != 0:
        issues.append("SPH provenance must report zero clipped update bins")
    if filled_bins != 0:
        issues.append(
            "physical SPH handoff must contain zero macrolib-filled XS bins; "
            f"found {filled_bins}"
        )
    try:
        numeric_residual = float(residual)
    except (TypeError, ValueError):
        numeric_residual = math.inf
    if not math.isfinite(numeric_residual) or numeric_residual < 0.0:
        issues.append(
            "SPH max update residual must be finite and non-negative; "
            f"got {numeric_residual!r}"
        )
    elif numeric_residual > max_update_residual:
        issues.append(
            "SPH fixed point is not converged: max update residual "
            f"{numeric_residual:.6g} exceeds {max_update_residual:.6g}"
        )
    _validate_provenance_bindings(attrs, source=source, issues=issues)
    return issues


def physical_colorset_sph_issues(
    path: str | Path,
    *,
    expected_domains: int = EXPECTED_COLORSET_DOMAINS,
    max_update_residual: float = PHYSICAL_SPH_MAX_UPDATE_RESIDUAL,
) -> list[str]:
    """Add the IRENA-style center-plus-six-neighbors colorset topology gate.

    This is a template contract, not the general definition of physical SPH.
    It remains available for existing IRENA manifests and other projects that
    explicitly choose the same seven-domain topology.
    """

    import h5py

    issues = physical_sph_issues(
        path,
        max_update_residual=max_update_residual,
    )
    with h5py.File(Path(path), "r") as h5:
        try:
            mixture_names = read_mixture_names(h5)
        except ValueError:
            mixture_names = ()
        if len(mixture_names) != expected_domains:
            issues.append(
                f"colorset must contain exactly {expected_domains} domains "
                f"(center target + six neighbors), found {len(mixture_names)}"
            )
        if mixture_names:
            first = h5["mixtures"][mixture_names[0]]
            try:
                first_index = int(first.attrs.get("source_domain_index", 0))
            except (TypeError, ValueError):
                first_index = 0
            if first_index != 1:
                issues.append(
                    "the center target must be the first declared domain "
                    "with source_domain_index=1"
                )
    return issues


def _text_attr(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace").strip()
    if isinstance(value, np.bytes_):
        return value.decode("utf-8", errors="replace").strip()
    return "" if value is None else str(value).strip()


def _boolean_attr(value: Any) -> bool | None:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)) and int(value) in (0, 1):
        return bool(value)
    text = _text_attr(value).lower()
    if text in {"true", "1"}:
        return True
    if text in {"false", "0"}:
        return False
    return None


def _integer_attr(value: Any) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError, OverflowError):
        return None


def _attribute_item_count(value: Any) -> int:
    if value is None:
        return 0
    try:
        return int(value.size)
    except (AttributeError, TypeError, ValueError):
        pass
    if isinstance(value, (str, bytes)):
        return 1 if value else 0
    try:
        return len(value)
    except TypeError:
        return 1


def _energy_groups(h5: Any, issues: list[str]) -> int | None:
    attr_groups: int | None = None
    bounds_groups: int | None = None
    if "energy_groups" in h5.attrs:
        try:
            attr_groups = int(h5.attrs["energy_groups"])
        except (TypeError, ValueError, OverflowError):
            issues.append("energy_groups must be a positive integer")
    if "energy_bounds" in h5:
        bounds = h5["energy_bounds"]
        if not hasattr(bounds, "shape") or len(bounds.shape) != 1:
            issues.append("/energy_bounds must be a one-dimensional dataset")
        else:
            bounds_groups = int(bounds.shape[0]) - 1
    groups = attr_groups if attr_groups is not None else bounds_groups
    if groups is None or groups <= 0:
        issues.append("physical SPH input must define a positive energy-group count")
        return None
    if (
        attr_groups is not None
        and bounds_groups is not None
        and attr_groups != bounds_groups
    ):
        issues.append(
            "energy_groups does not match /energy_bounds: "
            f"{attr_groups} != {bounds_groups}"
        )
        return None
    return groups


def _validate_applied_sph_payload(
    h5: Any,
    *,
    mixture_names: tuple[str, ...],
    energy_groups: int | None,
    issues: list[str],
) -> None:
    mixtures = h5.get("mixtures")
    if mixtures is None or not hasattr(mixtures, "keys"):
        issues.append("/mixtures must be an HDF5 group")
        return
    for mixture_name in mixture_names:
        if mixture_name not in mixtures:
            issues.append(f"declared mixture {mixture_name!r} is missing from /mixtures")
            continue
        mixture = mixtures[mixture_name]
        states = mixture.get("states") if hasattr(mixture, "get") else None
        if states is None:
            _validate_applied_sph_vector(
                mixture,
                label=f"/mixtures/{mixture_name}/applied_sph",
                energy_groups=energy_groups,
                issues=issues,
            )
            continue
        if not hasattr(states, "keys"):
            issues.append(f"/mixtures/{mixture_name}/states must be an HDF5 group")
            continue
        state_names = tuple(str(name) for name in states.keys())
        if len(state_names) != 1:
            issues.append(
                "multi-state physical SPH is not supported until state-specific "
                f"SPH factors are derived; mixture {mixture_name!r} declares "
                f"{len(state_names)} states"
            )
        _validate_applied_sph_vector(
            mixture,
            label=f"/mixtures/{mixture_name}/applied_sph",
            energy_groups=energy_groups,
            issues=issues,
        )
        for state_name in state_names:
            _validate_applied_sph_vector(
                states[state_name],
                label=(
                    f"/mixtures/{mixture_name}/states/{state_name}/applied_sph"
                ),
                energy_groups=energy_groups,
                issues=issues,
            )


def _validate_applied_sph_vector(
    group: Any,
    *,
    label: str,
    energy_groups: int | None,
    issues: list[str],
) -> None:
    if "applied_sph" not in group:
        issues.append(f"{label} is required")
        return
    dataset = group["applied_sph"]
    if not hasattr(dataset, "shape"):
        issues.append(f"{label} must be a dataset")
        return
    if energy_groups is None:
        return
    if tuple(dataset.shape) != (energy_groups,):
        issues.append(
            f"{label} must have shape ({energy_groups},), got {tuple(dataset.shape)}"
        )
        return
    try:
        values = np.asarray(dataset[:], dtype=float)
    except (TypeError, ValueError, OSError) as exc:
        issues.append(f"{label} is not a readable numeric vector: {exc}")
        return
    if not np.all(np.isfinite(values)):
        issues.append(f"{label} must contain only finite values")
    if np.any(values <= 0.0):
        issues.append(f"{label} must contain only positive values")


def _reject_active_equivalence_payloads(h5: Any, issues: list[str]) -> None:
    active_sph: list[str] = []
    adf_payloads: list[str] = []
    adf_attributes: list[str] = []

    for attribute in h5.attrs:
        normalized = str(attribute).strip().lower()
        if normalized.startswith("adf") or normalized == "discontinuity_factors":
            adf_attributes.append(f"/@{attribute}")

    def collect(name: str, obj: Any) -> None:
        basename = name.rsplit("/", 1)[-1].strip().lower()
        rendered = "/" + name.strip("/")
        if basename in _ACTIVE_SPH_NAMES:
            active_sph.append(rendered)
        if basename in _ADF_NAMES:
            adf_payloads.append(rendered)
        for attribute in obj.attrs:
            normalized = str(attribute).strip().lower()
            if normalized.startswith("adf") or normalized == "discontinuity_factors":
                adf_attributes.append(f"{rendered}@{attribute}")

    h5.visititems(collect)
    if active_sph:
        issues.append(
            "physical SPH input contains active sph/SPH/NSPH payload(s); "
            "apply-sph must remove them before Converter: "
            + ", ".join(sorted(set(active_sph)))
        )
    if adf_payloads or adf_attributes:
        rendered = sorted(set(adf_payloads + adf_attributes))
        issues.append(
            "physical SPH input must not contain any ADF payload or metadata: "
            + ", ".join(rendered)
        )


def _validate_provenance_bindings(
    attrs: dict[str, Any],
    *,
    source: str,
    issues: list[str],
) -> None:
    if _text_attr(attrs.get("sph_source_binding_schema")) != SPH_SOURCE_BINDING_SCHEMA:
        issues.append(
            f"sph_source_binding_schema must be {SPH_SOURCE_BINDING_SCHEMA}"
        )
    if _text_attr(attrs.get("sph_apply_binding_schema")) != SPH_APPLY_BINDING_SCHEMA:
        issues.append(
            f"sph_apply_binding_schema must be {SPH_APPLY_BINDING_SCHEMA}"
        )
    if _text_attr(attrs.get("sph_apply_binding_mode")) != (
        "converter-final-exact-input"
    ):
        issues.append(
            "sph_apply_binding_mode must be converter-final-exact-input; "
            "OpenMC-native setN applications are intermediate artifacts"
        )
    if _boolean_attr(
        attrs.get("sph_apply_sidecar_input_hash_verified")
    ) is not True:
        issues.append(
            "sph_apply_sidecar_input_hash_verified=true is required for the "
            "final Converter handoff"
        )

    digests: dict[str, str | None] = {}
    for prefix, label in (
        ("sph_input_h5", "sidecar input HDF5"),
        ("sph_reference_flux", "OpenMC CE reference flux"),
        ("sph_mg_flux", "OpenMC MG flux"),
        ("sph_apply_input_h5", "apply-sph input HDF5"),
        ("sph_apply_sidecar", "apply-sph sidecar"),
    ):
        digests[prefix] = _validate_file_binding(
            attrs,
            prefix=prefix,
            label=label,
            issues=issues,
        )

    if (
        digests["sph_input_h5"] is not None
        and digests["sph_apply_input_h5"] is not None
        and digests["sph_input_h5"] != digests["sph_apply_input_h5"]
    ):
        issues.append(
            "apply-sph input HDF5 SHA-256 does not match the HDF5 bound by "
            "the SPH sidecar"
        )

    previous_used = _boolean_attr(attrs.get("sph_previous_sph_used"))
    if previous_used is None:
        issues.append("sph_previous_sph_used must be an explicit boolean")
    elif previous_used:
        _validate_file_binding(
            attrs,
            prefix="sph_previous_sph",
            label="previous SPH",
            issues=issues,
        )

    for name, label in (
        ("sph_reference_flux_layout_verified", "OpenMC CE reference-flux layout"),
        ("sph_mg_flux_layout_verified", "OpenMC MG flux layout"),
    ):
        if _boolean_attr(attrs.get(name)) is not True:
            issues.append(
                f"{name}=true is required; {label} must attest verified "
                "energy bounds and spatial-domain order"
            )

    _validate_uncertainty_binding(
        attrs,
        prefix="sph_reference_flux",
        label="OpenMC CE reference flux",
        issues=issues,
    )
    _validate_uncertainty_binding(
        attrs,
        prefix="sph_mg_flux",
        label="OpenMC MG flux",
        issues=issues,
    )

    sidecar_path = _text_attr(attrs.get("sph_apply_sidecar_path"))
    if source and sidecar_path and not _same_file_reference(source, sidecar_path):
        issues.append(
            "sph_applied_source does not match the SHA-256-bound apply-sph sidecar"
        )


def _validate_uncertainty_binding(
    attrs: dict[str, Any],
    *,
    prefix: str,
    label: str,
    issues: list[str],
) -> None:
    require_coverage = _boolean_attr(
        attrs.get(f"{prefix}_uncertainty_require_coverage")
    )
    coverage = _boolean_attr(attrs.get(f"{prefix}_uncertainty_coverage"))
    passed = _boolean_attr(attrs.get(f"{prefix}_uncertainty_pass"))
    std_dev_dataset = _text_attr(attrs.get(f"{prefix}_std_dev_dataset"))
    limit = _finite_nonnegative_attr(attrs.get(f"{prefix}_uncertainty_limit"))
    observed = _finite_nonnegative_attr(
        attrs.get(f"{prefix}_uncertainty_observed_max_rel")
    )

    if require_coverage is not True:
        issues.append(
            f"{prefix}_uncertainty_require_coverage=true is required for "
            f"{label}"
        )
    if coverage is not True or not std_dev_dataset:
        issues.append(
            f"{prefix}_uncertainty_coverage=true and a std-dev dataset are "
            f"required for {label}"
        )
    if limit is None:
        issues.append(
            f"{prefix}_uncertainty_limit must be an explicit finite "
            "non-negative value"
        )
    if observed is None:
        issues.append(
            f"{prefix}_uncertainty_observed_max_rel must be finite and "
            "non-negative"
        )
    if limit is not None and observed is not None and observed > limit:
        issues.append(
            f"{label} relative uncertainty {observed:.6g} exceeds its "
            f"explicit limit {limit:.6g}"
        )
    if passed is not True:
        issues.append(f"{prefix}_uncertainty_pass=true is required")


def _validate_file_binding(
    attrs: dict[str, Any],
    *,
    prefix: str,
    label: str,
    issues: list[str],
) -> str | None:
    source = _text_attr(attrs.get(f"{prefix}_path"))
    digest = _text_attr(attrs.get(f"{prefix}_sha256")).lower()
    if not source:
        issues.append(f"{prefix}_path must record the {label} source")
    if not _is_sha256(digest):
        issues.append(f"{prefix}_sha256 must be a well-formed SHA-256 digest")
        return None
    if not source:
        return digest
    candidate = Path(source).expanduser()
    if not candidate.exists():
        return digest
    if not candidate.is_file():
        issues.append(f"{prefix}_path is not a regular file: {candidate}")
        return digest
    try:
        actual = file_sha256(candidate)
    except OSError as exc:
        issues.append(f"could not verify {label} SHA-256: {exc}")
    else:
        if actual.lower() != digest:
            issues.append(
                f"{prefix}_sha256 does not match the existing {label} source"
            )
    return digest


def _is_sha256(value: str) -> bool:
    if len(value) != 64:
        return False
    try:
        int(value, 16)
    except ValueError:
        return False
    return True


def _finite_nonnegative_attr(value: Any) -> float | None:
    try:
        numeric = float(value)
    except (TypeError, ValueError, OverflowError):
        return None
    if not math.isfinite(numeric) or numeric < 0.0:
        return None
    return numeric


def _same_file_reference(left: str, right: str) -> bool:
    try:
        return Path(left).expanduser().resolve() == Path(right).expanduser().resolve()
    except OSError:
        return left == right
