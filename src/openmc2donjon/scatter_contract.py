"""Serialize the scattering/removal contract carried by MGXS preflight reports."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any


SCATTER_CONTRACT_KEYS = (
    "openmc_scatter_mgxs_type",
    "openmc_scatter_multiplicity_weighted",
    "openmc_scatter_balance_dataset",
    "openmc_scatter_contract_declared",
    "openmc_scatter_contract_valid",
    "openmc_transport_mgxs_type",
    "openmc_transport_contract_declared",
)


def scatter_contract_from_preflight(
    preflight: Mapping[str, Any] | None,
) -> dict[str, Any] | None:
    """Return the common explicit scatter contract from a preflight payload.

    Converter currently accepts one HDF5 input, but the MGXS contract report
    format permits more than one.  Return a contract only when every reported
    input carries the same resolved convention; a receipt must not collapse
    conflicting source conventions into one apparently authoritative value.
    """

    if not isinstance(preflight, Mapping):
        return None
    inputs = preflight.get("inputs")
    if not isinstance(inputs, list) or not inputs:
        return None

    contracts: list[dict[str, Any]] = []
    for item in inputs:
        if not isinstance(item, Mapping):
            return None
        if any(key not in item for key in SCATTER_CONTRACT_KEYS):
            return None
        contract = {key: item.get(key) for key in SCATTER_CONTRACT_KEYS}
        # A null MGXS type may be the validator's explicit legacy ordinary
        # default or a contract declared through the weighting/balance fields.
        # The status fields preserve that distinction.  Never write a receipt
        # that turns a contradictory or incomplete declaration into an
        # apparently authoritative contract.
        if (
            contract["openmc_scatter_multiplicity_weighted"] is None
            or contract["openmc_scatter_balance_dataset"] is None
            or contract["openmc_scatter_contract_declared"] is None
            or contract["openmc_scatter_contract_valid"] is not True
            or contract["openmc_transport_contract_declared"] is None
        ):
            return None
        contracts.append(contract)

    first = contracts[0]
    if any(contract != first for contract in contracts[1:]):
        return None
    return first
