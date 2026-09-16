"""Opt-in native OpenMC currents on explicit, shared hexagonal node faces.

This module does not implement a HexMesh and does not modify an existing model.
The returned surfaces must be used to construct the actual geometry. Domains are
unique, uncut hexagonal prisms in global coordinates; their fine-geometry fills
may be translated. Only transmission faces are supported. OpenMC is imported
only when constructing geometry or tallies.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import importlib
import math
import operator
from typing import Any, Iterable, Mapping

import numpy as np


SCHEMA = "openmc2donjon.hex-face-current.v1"
FACE_NAMES = ("F0_N", "F1_NE", "F2_SE", "F3_S", "F4_SW", "F5_NW", "BOTTOM", "TOP")
_NEIGHBORS = ((0, 1), (1, 0), (1, -1), (0, -1), (-1, 0), (-1, 1))
_NORMALS = ((0.0, 1.0), (math.sqrt(3) / 2, 0.5), (math.sqrt(3) / 2, -0.5))
NodeKey = tuple[int, int, int]
Site = tuple[int, int]


def _openmc() -> Any:
    try:
        return importlib.import_module("openmc")
    except ImportError as exc:  # pragma: no cover - depends on user environment
        raise RuntimeError("OpenMC is required to construct native hexagonal current tallies") from exc


def _integer(value: Any, name: str) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be an integer, not a boolean")
    try:
        return operator.index(value)
    except TypeError as exc:
        raise ValueError(f"{name} must be an integer") from exc


def _site(value: Iterable[int]) -> Site:
    values = tuple(value)
    if len(values) != 2:
        raise ValueError("each hex site must contain exactly two integer lattice indices")
    return _integer(values[0], "site index"), _integer(values[1], "site index")


def _node(value: Iterable[int]) -> NodeKey:
    values = tuple(value)
    if len(values) != 3:
        raise ValueError("each node must contain two lattice indices and an axial layer index")
    return tuple(_integer(item, "node index") for item in values)  # type: ignore[return-value]


@dataclass(frozen=True)
class HexCurrentSurfaces:
    """Shared global planes for a y-oriented OpenMC HexLattice (x-oriented tiles).

    Site ``(i, j)`` has center ``(cx + sqrt(3)*pitch*i/2,
    cy + pitch*(j+i/2))``. Six lateral faces are numbered clockwise from north.
    A plane is shared globally, even by collinear nonadjacent faces. Cell filters
    restrict its tally to the actual finite node face.
    """

    sites: tuple[Site, ...]
    z_edges: tuple[float, ...]
    pitch: float
    center: tuple[float, float]
    planes: Mapping[tuple[int, int], Any]
    z_planes: tuple[Any, ...]

    @property
    def num_layers(self) -> int:
        return len(self.z_edges) - 1

    def site_center(self, site: Iterable[int]) -> tuple[float, float]:
        i, j = self._check_site(site)
        return self.center[0] + math.sqrt(3) * self.pitch * i / 2, self.center[1] + self.pitch * (j + i / 2)

    def _check_site(self, site: Iterable[int]) -> Site:
        result = _site(site)
        if result not in self.sites:
            raise ValueError(f"hex site {result} is not registered")
        return result

    def _check_layer(self, layer: int) -> int:
        layer = _integer(layer, "axial layer")
        if not 0 <= layer < self.num_layers:
            raise ValueError(f"axial layer must be in [0, {self.num_layers})")
        return layer

    def side_surfaces(self, site: Iterable[int]) -> tuple[Any, ...]:
        i, j = self._check_site(site)
        offsets = (2 * j + i, 2 * i + j, i - j)
        return tuple(
            self.planes[(face % 3, offsets[face % 3] + (1 if face < 3 else -1))] for face in range(6)
        )

    def tile_region(self, site: Iterable[int]) -> Any:
        planes = self.side_surfaces(site)
        region = -planes[0]
        for face in range(1, 6):
            region &= -planes[face] if face < 3 else +planes[face]
        return region

    def layer_region(self, layer: int) -> Any:
        layer = self._check_layer(layer)
        return +self.z_planes[layer] & -self.z_planes[layer + 1]

    def node_region(self, site: Iterable[int], layer: int) -> Any:
        return self.tile_region(site) & self.layer_region(layer)

    def node_surfaces(self, site: Iterable[int], layer: int) -> tuple[Any, ...]:
        layer = self._check_layer(layer)
        return self.side_surfaces(site) + self.z_planes[layer : layer + 2]

    def node_vertices(self, site: Iterable[int], layer: int) -> np.ndarray:
        """Return the twelve prism vertices, used to validate ancestor clipping."""
        layer = self._check_layer(layer)
        x, y = self.site_center(site)
        radius = self.pitch / math.sqrt(3)
        return np.asarray(
            [
                (x + radius * math.cos(n * math.pi / 3), y + radius * math.sin(n * math.pi / 3), z)
                for z in self.z_edges[layer : layer + 2]
                for n in range(6)
            ]
        )


def build_hex_current_surfaces(
    sites: Iterable[Iterable[int]],
    z_edges: Iterable[float],
    pitch: float,
    center: Iterable[float] = (0.0, 0.0),
) -> HexCurrentSurfaces:
    """Build transmission planes without changing a user's existing geometry.

    Axial bounds and pitch are in cm. Axial bounds must be finite and strictly
    increasing. Only one uniform, y-oriented hexagonal lattice is supported.
    """
    site_values = tuple(_site(site) for site in sites)
    if not site_values or len(set(site_values)) != len(site_values):
        raise ValueError("sites must be nonempty and unique")
    pitch = float(pitch)
    if not math.isfinite(pitch) or pitch <= 0:
        raise ValueError("pitch must be finite and positive")
    center_values = tuple(float(value) for value in center)
    if len(center_values) != 2 or not all(math.isfinite(value) for value in center_values):
        raise ValueError("center must contain two finite coordinates")
    z_values = tuple(float(value) for value in z_edges)
    if (
        len(z_values) < 2
        or not all(math.isfinite(value) for value in z_values)
        or any(b <= a for a, b in zip(z_values, z_values[1:], strict=False))
    ):
        raise ValueError("z_edges must contain at least two finite, strictly increasing values")
    openmc = _openmc()
    plane_keys: set[tuple[int, int]] = set()
    for i, j in site_values:
        for family, offset in enumerate((2 * j + i, 2 * i + j, i - j)):
            plane_keys.update(((family, offset - 1), (family, offset + 1)))
    planes = {}
    for family, offset in sorted(plane_keys):
        nx, ny = _NORMALS[family]
        planes[(family, offset)] = openmc.Plane(
            a=nx,
            b=ny,
            c=0,
            d=nx * center_values[0] + ny * center_values[1] + pitch * offset / 2,
            name=f"hex-current-plane-{family}-{offset}",
            boundary_type="transmission",
        )
    z_planes = tuple(
        openmc.ZPlane(z0=z, name=f"hex-current-z-{k}", boundary_type="transmission")
        for k, z in enumerate(z_values)
    )
    return HexCurrentSurfaces(site_values, z_values, pitch, center_values, planes, z_planes)


@dataclass(frozen=True)
class HexCurrentPlan:
    """Native tallies and a JSON-serializable description of their exact bins."""

    tallies: Any
    manifest: dict[str, Any]


def _halfspaces(region: Any) -> list[Any]:
    if region is None:
        return []
    if hasattr(region, "surface") and hasattr(region, "side"):
        return [region]
    if type(region).__name__ != "Intersection":
        raise ValueError("current nodes and their ancestors require intersections of planar halfspaces")
    return [halfspace for child in region for halfspace in _halfspaces(child)]


def _plane_coefficients(surface: Any) -> np.ndarray:
    values = surface.coefficients
    kind = type(surface).__name__
    if kind == "Plane":
        return np.asarray([values["a"], values["b"], values["c"], values["d"]])
    if kind in ("XPlane", "YPlane", "ZPlane"):
        axis = "XYZ".index(kind[0])
        coeff = np.zeros(4)
        coeff[axis] = 1
        coeff[3] = values[kind[0].lower() + "0"]
        return coeff
    raise ValueError("current-node ancestors must use planar boundaries; curved clipping is unsupported")


def _flat_universes(value: Any) -> Iterable[Any]:
    if hasattr(value, "cells"):
        yield value
    else:
        for item in value:
            yield from _flat_universes(item)


def _validate_placement(geometry: Any, domains: Mapping[NodeKey, Any], registry: HexCurrentSurfaces) -> None:
    """Reject reused or transformed nodes without enumerating fine-cell paths."""
    target_ids = {cell.id for cell in domains.values()}
    reachable: dict[int, bool] = {}
    active: set[int] = set()

    def children(cell: Any) -> list[tuple[Any, bool]]:
        if cell.fill_type == "universe":
            return [(cell.fill, False)]
        if cell.fill_type == "lattice":
            result = [(universe, True) for universe in _flat_universes(cell.fill.universes)]
            if cell.fill.outer is not None:
                result.append((cell.fill.outer, True))
            return result
        return []

    def has_target(universe: Any) -> bool:
        if universe.id in reachable:
            return reachable[universe.id]
        if universe.id in active:
            raise ValueError("cyclic universe fills are unsupported")
        active.add(universe.id)
        found = False
        for cell in universe.cells.values():
            found |= cell.id in target_ids
            for child, _ in children(cell):
                found |= has_target(child)
        active.remove(universe.id)
        reachable[universe.id] = found
        return found

    has_target(geometry.root_universe)
    counts: Counter[int] = Counter()
    lookup = {cell.id: key for key, cell in domains.items()}

    def walk(universe: Any, ancestors: list[Any]) -> None:
        for cell in universe.cells.values():
            if cell.id in target_ids:
                counts[cell.id] += 1
                if counts[cell.id] != 1:
                    raise ValueError(
                        f"current domain cell {cell.id} is reused; unique spatial cells are required"
                    )
                key = lookup[cell.id]
                if cell is not domains[key]:
                    raise ValueError(f"geometry cell {cell.id} is not the supplied domain object")
                vertices = registry.node_vertices(key[:2], key[2])
                node_surfaces = registry.node_surfaces(key[:2], key[2])
                for ancestor in ancestors:
                    for halfspace in _halfspaces(ancestor.region):
                        coefficients = _plane_coefficients(halfspace.surface)
                        signed = vertices @ coefficients[:3] - coefficients[3]
                        if halfspace.side == "-":
                            signed *= -1
                        tolerance = 1e-9 * max(1.0, float(np.linalg.norm(coefficients[:3])))
                        if np.any(signed < -tolerance):
                            raise ValueError(f"ancestor cell {ancestor.id} clips current domain {cell.id}")
                        if halfspace.surface.boundary_type != "transmission" and np.any(signed <= tolerance):
                            raise ValueError("a non-transmission ancestor boundary touches a current domain")
                        if np.count_nonzero(np.abs(signed) <= tolerance) >= 4:
                            normal = coefficients / np.linalg.norm(coefficients[:3])
                            for node_surface in node_surfaces:
                                expected = _plane_coefficients(node_surface)
                                expected /= np.linalg.norm(expected[:3])
                                if np.allclose(normal, expected, rtol=0, atol=1e-10) or np.allclose(
                                    normal, -expected, rtol=0, atol=1e-10
                                ):
                                    if halfspace.surface is not node_surface:
                                        raise ValueError(
                                            "coincident ancestor faces must share registry surface objects"
                                        )
            for child, through_lattice in children(cell):
                if not reachable[child.id]:
                    continue
                if through_lattice:
                    raise ValueError(
                        "current-domain ancestors may not be lattices; use explicit global node cells"
                    )
                if (cell.translation is not None and np.any(np.asarray(cell.translation) != 0)) or (
                    cell.rotation is not None and np.any(np.asarray(cell.rotation) != 0)
                ):
                    raise ValueError("current-domain ancestors may not be translated or rotated")
                walk(child, ancestors + [cell])

    walk(geometry.root_universe, [])
    if set(counts) != target_ids:
        raise ValueError(f"current domain cells absent from geometry: {sorted(target_ids - set(counts))}")


def bind_hex_current_domains(
    registry: HexCurrentSurfaces,
    domains: Mapping[NodeKey, Any],
    *,
    energy_bounds: Iterable[float],
    geometry: Any,
    name_prefix: str = "hex-current",
) -> HexCurrentPlan:
    """Create outgoing/incoming native ``current`` tallies for unique nodes.

    Each domain's region must be exactly ``registry.node_region(site, layer)``.
    Its cell may translate/rotate its *fill*, but no ancestor may transform the
    global face planes. Ancestors must be planar convex regions that fully
    contain the prism. Nodes may not lie within a lattice or use non-transmission
    faces. These restrictions make the area and crossing-direction contract
    checkable, rather than silently accepting a clipped or repeated domain.

    Native signed surface current is converted to positive partial-current
    magnitudes using the manifest's outward sign. Do not merge these tallies.
    ``net_std_dev`` is intentionally unavailable without correlated batch data.
    """
    if not domains:
        raise ValueError("at least one current domain is required")
    if not isinstance(name_prefix, str) or not name_prefix.strip():
        raise ValueError("name_prefix must be nonempty")
    domain_values = {_node(key): cell for key, cell in domains.items()}
    if len(domain_values) != len(domains):
        raise ValueError("duplicate current node keys")
    cell_ids = [_integer(cell.id, "domain cell ID") for cell in domain_values.values()]
    if min(cell_ids) <= 0 or len(set(cell_ids)) != len(cell_ids):
        raise ValueError("current domains require unique, positive cell IDs")
    bounds = np.asarray(tuple(energy_bounds), dtype=float)
    if bounds.ndim != 1 or len(bounds) < 2 or not np.all(np.isfinite(bounds)):
        raise ValueError("energy_bounds must contain at least two finite values")
    if bounds[0] < 0 or np.any(np.diff(bounds) <= 0):
        raise ValueError("energy_bounds must be nonnegative and strictly increasing")
    registered_surfaces = list(registry.planes.values()) + list(registry.z_planes)
    if len({surface.id for surface in registered_surfaces}) != len(registered_surfaces):
        raise ValueError("registry surfaces must have unique IDs")
    for (family, offset), surface in registry.planes.items():
        nx, ny = _NORMALS[family]
        expected = [
            nx,
            ny,
            0,
            nx * registry.center[0] + ny * registry.center[1] + registry.pitch * offset / 2,
        ]
        if not np.array_equal(_plane_coefficients(surface), expected):
            raise ValueError("registry plane coefficients were modified after construction")
    for surface, z in zip(registry.z_planes, registry.z_edges, strict=True):
        if not np.array_equal(_plane_coefficients(surface), [0, 0, 1, z]):
            raise ValueError("registry axial plane coefficients were modified after construction")
    for (i, j, k), cell in domain_values.items():
        expected_surfaces = registry.node_surfaces((i, j), k)
        halfspaces = _halfspaces(cell.region)
        expected_sides = ("-", "-", "-", "+", "+", "+", "+", "-")
        if Counter((id(h.surface), h.side) for h in halfspaces) != Counter(
            (id(surface), side) for surface, side in zip(expected_surfaces, expected_sides, strict=True)
        ):
            raise ValueError(f"domain cell {cell.id} must use the registry's exact eight node halfspaces")
        if any(surface.boundary_type != "transmission" for surface in expected_surfaces):
            raise ValueError("only transmission current faces are supported")
    _validate_placement(geometry, domain_values, registry)
    openmc = _openmc()
    tallies = openmc.Tallies()
    nodes = []
    energy_filter = openmc.EnergyFilter(bounds)
    for (i, j, k), cell in domain_values.items():
        surfaces = registry.node_surfaces((i, j), k)
        surface_filter = openmc.SurfaceFilter(surfaces)
        node_tallies = {}
        for direction, cell_filter in (
            ("outgoing", openmc.CellFromFilter([cell])),
            ("incoming", openmc.CellFilter([cell])),
        ):
            tally = openmc.Tally(name=f"{name_prefix}-{i}-{j}-{k}-{direction}")
            tally.filters = [surface_filter, cell_filter, energy_filter]
            tally.scores = ["current"]
            # OpenMC chooses analog for explicit surface-current tallies.
            tally.estimator = "analog"
            tallies.append(tally, merge=False)
            node_tallies[direction] = {
                "id": tally.id,
                "name": tally.name,
                "filter_order": ["surface", "cellfrom" if direction == "outgoing" else "cell", "energy"],
            }
        side_area = registry.pitch / math.sqrt(3) * (registry.z_edges[k + 1] - registry.z_edges[k])
        cap_area = math.sqrt(3) / 2 * registry.pitch**2
        faces = []
        for face, surface in enumerate(surfaces):
            sign = 1 if face < 3 or face == 7 else -1
            if face < 6:
                nx, ny = _NORMALS[face % 3]
                normal = [sign * nx, sign * ny, 0.0]
                di, dj = _NEIGHBORS[face]
                neighbor = (i + di, j + dj, k)
                opposite = (face + 3) % 6
            else:
                normal = [0.0, 0.0, float(sign)]
                neighbor = (i, j, k + sign)
                opposite = 7 if face == 6 else 6
            faces.append(
                {
                    "index": face,
                    "name": FACE_NAMES[face],
                    "surface_id": surface.id,
                    "outward_sign": sign,
                    "outward_normal": normal,
                    "area_cm2": side_area if face < 6 else cap_area,
                    "neighbor_node": list(neighbor) if neighbor in domain_values else None,
                    "neighbor_cell_id": domain_values[neighbor].id if neighbor in domain_values else None,
                    "neighbor_face_index": opposite if neighbor in domain_values else None,
                    "boundary_type": "transmission",
                }
            )
        nodes.append(
            {
                "node": [i, j, k],
                "cell_id": cell.id,
                "cell_name": cell.name,
                "center_xy_cm": list(registry.site_center((i, j))),
                "z_bounds_cm": list(registry.z_edges[k : k + 2]),
                "faces": faces,
                "tallies": node_tallies,
            }
        )
    manifest = {
        "schema": SCHEMA,
        "status": "native_tally_definition_not_physics_acceptance",
        "lattice_orientation": "y",
        "tile_orientation": "x",
        "pitch_cm": registry.pitch,
        "center_xy_cm": list(registry.center),
        "z_edges_cm": list(registry.z_edges),
        "energy_bounds_eV": bounds.tolist(),
        "energy_order": "low-to-high (native OpenMC EnergyFilter)",
        "score": "current",
        "estimator": "analog",
        "particle": "neutron (caller must run neutron-only model)",
        "normalization": "weighted crossings per starting source particle; not area-normalized",
        "current_density_normalization": "divide each face current by its area_cm2",
        "outgoing_formula": "outward_sign * native outgoing current",
        "incoming_formula": "-outward_sign * native incoming current",
        "net_formula": "outgoing - incoming; positive outward",
        "net_std_dev": None,
        "net_std_dev_reason": "partial-current estimates are correlated; batch covariance is not supplied",
        "scope": "explicit unique uncut global hex-prism cells; transmission faces only",
        "nodes": nodes,
    }
    plan = HexCurrentPlan(tallies, manifest)
    validate_hex_current_tally_xml(plan, tallies.to_xml_element())
    return plan


def validate_hex_current_tally_xml(plan: HexCurrentPlan, root: Any) -> None:
    """Check serialized bins exactly (not OpenMC's approximate Filter equality)."""
    filters = {int(element.get("id")): element for element in root.findall("filter")}
    tallies = {int(element.get("id")): element for element in root.findall("tally")}
    expected_energy = np.asarray(plan.manifest["energy_bounds_eV"])
    for node in plan.manifest["nodes"]:
        surfaces = [face["surface_id"] for face in node["faces"]]
        for direction, definition in node["tallies"].items():
            tally = tallies.get(definition["id"])
            if tally is None or tally.get("name") != definition["name"]:
                raise ValueError("current tally was omitted or renamed during serialization")
            if tally.findtext("scores", "").split() != ["current"]:
                raise ValueError("current tally has incorrect scores")
            if tally.findtext("estimator") != "analog":
                raise ValueError("current tally requires the analog estimator")
            references = [int(value) for value in tally.findtext("filters", "").split()]
            if len(references) != 3 or any(value not in filters for value in references):
                raise ValueError("current tally requires exact surface, cell and energy filters")
            selected = [filters[value] for value in references]
            expected_types = ["surface", "cellfrom" if direction == "outgoing" else "cell", "energy"]
            if [value.get("type") for value in selected] != expected_types:
                raise ValueError("current tally filter order or types changed")
            actual_surface = [int(value) for value in selected[0].findtext("bins", "").split()]
            actual_cell = [int(value) for value in selected[1].findtext("bins", "").split()]
            actual_energy = np.asarray([float(value) for value in selected[2].findtext("bins", "").split()])
            if (
                actual_surface != surfaces
                or actual_cell != [node["cell_id"]]
                or not np.array_equal(actual_energy, expected_energy)
            ):
                raise ValueError("current tally bins changed during serialization")


def signed_partial_currents(
    outgoing_mean: Any,
    incoming_mean: Any,
    outgoing_std_dev: Any,
    incoming_std_dev: Any,
    *,
    outward_signs: Iterable[int],
    face_areas: Iterable[float] | None = None,
) -> dict[str, Any]:
    """Convert native signed arrays ``[face, group]`` to outward/inward magnitudes.

    Arrays retain native low-to-high energy order. A zero tally remains zero;
    this function neither fills unscored bins nor declares statistical adequacy.
    Negative converted partial means are rejected. Net standard deviation is
    not reconstructed from marginal standard deviations because covariance is
    unknown. Optional area normalization yields current density, not scalar flux.
    """
    arrays = [
        np.asarray(value, dtype=float)
        for value in (outgoing_mean, incoming_mean, outgoing_std_dev, incoming_std_dev)
    ]
    if arrays[0].ndim != 2 or arrays[0].size == 0 or any(value.shape != arrays[0].shape for value in arrays):
        raise ValueError("all current mean/std_dev arrays must have the same nonempty [face, group] shape")
    if any(not np.all(np.isfinite(value)) for value in arrays):
        raise ValueError("current mean/std_dev values must be finite")
    signs = np.asarray(tuple(outward_signs), dtype=float)
    if signs.shape != (arrays[0].shape[0],) or not np.all(np.isin(signs, (-1, 1))):
        raise ValueError("outward_signs must contain one +1/-1 per face")
    if any(np.any(value < 0) for value in arrays[2:]):
        raise ValueError("current standard deviations must be nonnegative")
    outgoing, incoming = signs[:, None] * arrays[0], -signs[:, None] * arrays[1]
    if np.any(outgoing < 0) or np.any(incoming < 0):
        raise ValueError("converted partial current is negative; check normal directions and tally selection")
    scale = np.ones((len(signs), 1))
    if face_areas is not None:
        areas = np.asarray(tuple(face_areas), dtype=float)
        if areas.shape != signs.shape or not np.all(np.isfinite(areas)) or np.any(areas <= 0):
            raise ValueError("face_areas must contain one finite positive area per face")
        scale = areas[:, None]
    return {
        "outgoing_mean": outgoing / scale,
        "incoming_mean": incoming / scale,
        "outgoing_std_dev": arrays[2] / scale,
        "incoming_std_dev": arrays[3] / scale,
        "net_mean": (outgoing - incoming) / scale,
        "net_std_dev": None,
        "net_std_dev_reason": "requires covariance from correlated batch estimates",
        "area_normalized": face_areas is not None,
        "energy_order": "low-to-high (native OpenMC EnergyFilter)",
    }
