"""Opt-in, exact cell-domain Pn preparation using native OpenMC MGXS.

OpenMC is imported lazily. Configure a freshly built Library before accessing
its tallies. This helper does not change geometry or estimator definitions.
"""
from __future__ import annotations

import copy
import json
from numbers import Integral
from typing import Any, Mapping

import numpy as np


def set_domain_scatter_orders(library: Any, orders: Mapping[int, int]) -> None:
    """Assign explicit scattering orders to every unique cell domain, atomically."""
    import openmc

    if library.domain_type != "cell" or library.scatter_format != "legendre":
        raise ValueError("Mixed-order helper requires cell domains and Legendre scattering")
    domains = list(library.domains)
    ids = [domain.id for domain in domains]
    if not ids or len(set(ids)) != len(ids) or set(orders) != set(ids):
        raise ValueError("Declare exactly one order for every unique domain cell ID")
    if any(isinstance(key, bool) or not isinstance(key, Integral) for key in orders):
        raise ValueError("Domain IDs must be integers")
    if any(isinstance(v, bool) or not isinstance(v, Integral) or not 0 <= v <= 10
           for v in orders.values()):
        raise ValueError("Domain Legendre orders must be integers in [0, 10]")
    changes = []
    for domain in domains:
        objects = [library.get_mgxs(domain, kind) for kind in library.mgxs_types]
        matrices = [obj for obj in objects if isinstance(obj, openmc.mgxs.ScatterMatrixXS)]
        if not matrices:
            raise ValueError(f"No native scattering matrix declared for cell {domain.id}")
        for obj in objects:
            if any(getattr(obj, name, None) is not None
                   for name in ("_tallies", "_xs_tally", "_rxn_rate_tally")):
                raise ValueError("Set domain orders before creating or loading ANY MGXS tallies")
        for obj in matrices:
            if obj.scatter_format != "legendre" or obj.correction is not None:
                raise ValueError("Set correction=None and use native Legendre scattering")
            changes.append((obj, int(orders[domain.id])))
    library.legendre_order = max(orders.values())
    for obj, order in changes:
        obj.legendre_order = order
    library.domain_legendre_orders = dict(orders)


def exact_cell_mgxs_tallies(library: Any) -> Any:
    """Group identical native templates with exact, never approximate, cell IDs.

    Scalar templates shared by orders are combined. Legendre filters of
    different orders remain separate, so P1 cells never score P2/P3 here.
    """
    import openmc

    if library.domain_type != "cell" or library.by_nuclide:
        raise ValueError("Exact grouped tallies currently require macroscopic cell-domain MGXS")
    definitions: dict[str, tuple[Any, list[int]]] = {}
    for domain in library.domains:
        for kind in library.mgxs_types:
            for tally in library.get_mgxs(domain, kind).tallies.values():
                if tally.derivative is not None or tally.triggers or not tally.multiply_density:
                    raise ValueError("Unsupported derivative, trigger or density setting in native template")
                cells = [f for f in tally.filters if isinstance(f, openmc.CellFilter)]
                if len(cells) != 1 or list(cells[0].bins) != [domain.id]:
                    raise ValueError("Native template must select exactly its declared domain cell")
                filters = [(type(f).__name__, "CELL" if isinstance(f, openmc.CellFilter)
                            else np.asarray(f.bins).tolist()) for f in tally.filters]
                key = json.dumps([tally.estimator, list(tally.scores), list(tally.nuclides), filters])
                if key not in definitions:
                    definitions[key] = (tally, [])
                if domain.id not in definitions[key][1]:
                    definitions[key][1].append(domain.id)
    result = openmc.Tallies()
    filter_pool: dict[str, Any] = {}
    for index, (template, ids) in enumerate(definitions.values(), 1):
        tally = openmc.Tally(name=f"mgxs_{index:02d}")
        tally.estimator, tally.scores, tally.nuclides = (
            template.estimator, list(template.scores), list(template.nuclides))
        filters = []
        for f in template.filters:
            bins = ids if isinstance(f, openmc.CellFilter) else np.asarray(f.bins).tolist()
            key = json.dumps([type(f).__name__, bins])
            if key not in filter_pool:
                if isinstance(f, openmc.CellFilter):
                    clone = openmc.CellFilter(ids)
                else:
                    clone = copy.deepcopy(f)
                    clone.id = None
                filter_pool[key] = clone
            filters.append(filter_pool[key])
        tally.filters = filters
        result.append(tally, merge=False)
    return result
