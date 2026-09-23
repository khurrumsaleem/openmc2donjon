"""Explicit local Legendre truncation inside a common storage layout."""
from __future__ import annotations

from numbers import Integral
from typing import Any, Mapping

import numpy as np

ORDER_ATTRS = frozenset({"source_legendre_order", "scatter_padding"})


def validate_order_metadata(
    attrs: Mapping[str, Any], scatter: np.ndarray, *, axes: str | None,
    stored_order: int, std_dev: np.ndarray | None = None,
) -> int | None:
    """Check declared truncation, not infer unmeasured moments from zeros.

    Legacy inputs without either attribute remain supported. Zero padding is
    a deterministic truncation choice, never evidence of zero physical moments.
    """
    present = ORDER_ATTRS.intersection(attrs)
    if not present:
        return None
    if present != ORDER_ATTRS:
        raise ValueError("source_legendre_order and scatter_padding must be declared together")
    order = attrs["source_legendre_order"]
    if isinstance(order, (bool, np.bool_)) or not isinstance(order, Integral):
        raise ValueError("source_legendre_order must be an integer")
    order = int(order)
    if not 0 <= order <= stored_order:
        raise ValueError("source_legendre_order must lie between zero and stored legendre_order")
    policy = attrs["scatter_padding"]
    if isinstance(policy, bytes):
        policy = policy.decode()
    expected = "zero-truncation" if order < stored_order else "none"
    if policy != expected:
        raise ValueError(f"scatter_padding must be {expected!r} for this source/stored order")
    if order == stored_order:
        return order
    if isinstance(axes, bytes):
        axes = axes.decode()
    normalized = "" if axes is None else axes.lower().replace(" ", "").replace("_", "")
    first = normalized in {"moment,from,to", "moment,in,out", "moment,gin,gout",
                           "legendre,from,to", "legendre,gin,gout"}
    last = normalized in {"from,to,moment", "in,out,moment", "gin,gout,moment",
                          "from,to,legendre", "gin,gout,legendre"}
    if normalized and not first and not last:
        raise ValueError(f"unsupported scatter_axes={axes!r}")
    if not normalized:
        first_match = scatter.ndim == 3 and scatter.shape[0] == stored_order + 1
        last_match = scatter.ndim == 3 and scatter.shape[-1] == stored_order + 1
        if first_match == last_match:
            raise ValueError("explicit scatter_axes required to validate padded moments")
        last = last_match
    for label, array in (("scatter_matrix", scatter), ("scatter_matrix_std_dev", std_dev)):
        if array is None:
            continue
        array = np.asarray(array)
        if array.shape != scatter.shape or array.ndim != 3:
            raise ValueError(f"{label} shape does not match the padded scatter layout")
        if array.shape[-1 if last else 0] != stored_order + 1:
            raise ValueError(f"{label} moment dimension does not match stored legendre_order")
        tail = array[..., order + 1:] if last else array[order + 1:]
        if np.any(tail != 0.0):
            raise ValueError(f"{label} contains nonzero/nonfinite declared truncated moments")
    return order
