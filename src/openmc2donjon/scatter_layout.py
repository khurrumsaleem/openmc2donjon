"""Shared HDF5 scattering layout for conversion, collapse, and SPH application."""

from __future__ import annotations

import numpy as np

from .mgxs_input_scatter import (
    MOMENT_FIRST_SCATTER_AXES,
    MOMENT_LAST_SCATTER_AXES,
    normalize_axes,
)
from .scatter_order import validate_order_metadata


CANONICAL_SCATTER_AXES = "moment,from,to"


def scatter_axes_from_attrs(group_attrs, parent_attrs, root_attrs) -> str | None:
    """Use calculation → mixture → root precedence, including the axes alias."""

    for attrs in (group_attrs, parent_attrs, root_attrs):
        if attrs is None:
            continue
        for name in ("scatter_axes", "axes"):
            value = attrs.get(name)
            if value is not None:
                value = value.decode() if isinstance(value, bytes) else str(value)
                if value:
                    return value
    return None


def normalise_scatter_axes(
    raw: np.ndarray,
    ngroups: int,
    mix_name: str,
    *,
    expected_moments: int | None,
    axes: str | None,
) -> np.ndarray:
    """Normalize a 3D matrix using the Converter's explicit/inferred layout rules."""

    resolved = resolve_scatter_axes(
        raw, ngroups, mix_name, expected_moments=expected_moments, axes=axes,
    )
    return raw if resolved == CANONICAL_SCATTER_AXES else np.moveaxis(raw, -1, 0)


def resolve_scatter_axes(
    raw: np.ndarray,
    ngroups: int,
    mix_name: str,
    *,
    expected_moments: int | None,
    axes: str | None,
) -> str:
    """Resolve a 3D layout without changing its storage or guessing a cubic one."""

    if axes is not None:
        normalized = normalize_axes(axes)
        if normalized in MOMENT_FIRST_SCATTER_AXES:
            return CANONICAL_SCATTER_AXES
        if normalized in MOMENT_LAST_SCATTER_AXES:
            return "from,to,moment"
        raise ValueError(
            f"mixture {mix_name}: unsupported scatter_axes={axes!r}; "
            "expected 'moment,from,to' or 'from,to,moment'"
        )

    moment_first_shape = raw.shape[1:] == (ngroups, ngroups)
    moment_last_shape = raw.shape[:2] == (ngroups, ngroups)

    if expected_moments is not None:
        first_matches = moment_first_shape and raw.shape[0] == expected_moments
        last_matches = moment_last_shape and raw.shape[2] == expected_moments
        if first_matches and not last_matches:
            return CANONICAL_SCATTER_AXES
        if last_matches and not first_matches:
            return "from,to,moment"

    if moment_first_shape and moment_last_shape:
        raise ValueError(
            f"mixture {mix_name}: ambiguous scatter_matrix shape {raw.shape}; "
            "set scatter_axes='moment,from,to' or 'from,to,moment'"
        )

    if moment_first_shape and not moment_last_shape:
        return CANONICAL_SCATTER_AXES
    if moment_last_shape and not moment_first_shape:
        return "from,to,moment"

    raise ValueError(
        f"mixture {mix_name}: scatter_matrix shape {raw.shape} is not compatible "
        f"with {ngroups} groups"
    )


def read_collapse_scatter(group, ngroups: int) -> tuple[np.ndarray, np.ndarray | None]:
    """Read and validate means and matching uncertainties before any collapse.

    Collapse accepts direct (single-state) mixture groups. The returned arrays
    are always (moment, from, to), including 2D P0 input. Missing uncertainty is
    kept missing; it is never replaced with a zero-error estimate.
    """

    stored_order = int(group.file.attrs.get("legendre_order", 0))
    expected = (stored_order + 1, ngroups, ngroups)
    axes = scatter_axes_from_attrs(group.attrs, None, group.file.attrs)
    raw = np.asarray(group["scatter_matrix"][:], dtype=float)

    def canonical(values):
        if values.ndim == 2:
            result = values[np.newaxis]
        elif values.ndim == 3:
            result = normalise_scatter_axes(
                values, ngroups, group.name,
                expected_moments=stored_order + 1, axes=axes,
            )
        else:
            raise ValueError(f"mixture {group.name}: scatter_matrix must be 2D or 3D")
        if result.shape != expected:
            raise ValueError(
                f"mixture {group.name}: normalized scatter_matrix shape {result.shape} "
                f"expected {expected} from energy_groups/legendre_order"
            )
        return result

    scatter = canonical(raw)
    std = None
    if "scatter_matrix_std_dev" in group:
        raw_std = np.asarray(group["scatter_matrix_std_dev"][:], dtype=float)
        if raw_std.shape != raw.shape:
            raise ValueError(f"mixture {group.name}: scatter standard-deviation shape mismatch")
        std = canonical(raw_std)
    validate_order_metadata(
        dict(group.attrs), scatter, axes=CANONICAL_SCATTER_AXES,
        stored_order=stored_order, std_dev=std,
    )
    if not np.all(np.isfinite(scatter)) or np.any(scatter[0] < 0.0):
        raise ValueError(f"mixture {group.name}: scatter_matrix must be finite with non-negative P0")
    if std is not None and (not np.all(np.isfinite(std)) or np.any(std < 0.0)):
        raise ValueError(f"mixture {group.name}: scatter standard deviations must be finite and non-negative")
    return scatter, std


def write_canonical_scatter_axes(attrs) -> None:
    """Prevent copied input metadata from relabeling a normalized output."""

    attrs["scatter_axes"] = CANONICAL_SCATTER_AXES
    if "axes" in attrs:
        attrs["axes"] = CANONICAL_SCATTER_AXES
