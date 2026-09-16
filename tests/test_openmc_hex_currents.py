from __future__ import annotations

import json
import math
import unittest

import numpy as np

from openmc2donjon.openmc_hex_currents import (
    FACE_NAMES,
    bind_hex_current_domains,
    build_hex_current_surfaces,
    signed_partial_currents,
    validate_hex_current_tally_xml,
)

try:
    import openmc
except ImportError:
    openmc = None


class PartialCurrentTests(unittest.TestCase):
    def test_sign_area_and_energy_order(self) -> None:
        result = signed_partial_currents(
            [[2, 4], [-6, -8]],
            [[-1, -3], [2, 5]],
            [[0.2, 0.4], [0.6, 0.8]],
            [[0.1, 0.3], [0.2, 0.5]],
            outward_signs=(1, -1),
            face_areas=(2, 4),
        )
        np.testing.assert_array_equal(result["outgoing_mean"], [[1, 2], [1.5, 2]])
        np.testing.assert_array_equal(result["incoming_mean"], [[0.5, 1.5], [0.5, 1.25]])
        np.testing.assert_array_equal(result["net_mean"], [[0.5, 0.5], [1, 0.75]])
        np.testing.assert_allclose(result["outgoing_std_dev"], [[0.1, 0.2], [0.15, 0.2]])
        self.assertIsNone(result["net_std_dev"])
        self.assertIn("covariance", result["net_std_dev_reason"])
        self.assertIn("low-to-high", result["energy_order"])

    def test_zero_currents_remain_zero_without_acceptance_claim(self) -> None:
        zeros = np.zeros((8, 3))
        result = signed_partial_currents(
            zeros, zeros, zeros, zeros, outward_signs=(1, 1, 1, -1, -1, -1, -1, 1)
        )
        np.testing.assert_array_equal(result["net_mean"], zeros)
        self.assertFalse(result["area_normalized"])
        self.assertIsNone(result["net_std_dev"])

    def test_invalid_arrays_rejected(self) -> None:
        valid = np.ones((2, 1))
        cases = [
            (valid[:, 0], -valid, valid, valid, (1, 1), None),
            (valid, -valid, np.zeros((1, 1)), valid, (1, 1), None),
            (valid * np.nan, -valid, valid, valid, (1, 1), None),
            (valid, -valid, -valid, valid, (1, 1), None),
            (valid, -valid, valid, valid, (0, 1), None),
            (valid, -valid, valid, valid, (1,), None),
            (valid, -valid, valid, valid, (1, 1), (0, 2)),
            (valid, -valid, valid, valid, (1, 1), (1,)),
            (-valid, -valid, valid, valid, (1, 1), None),
            (valid, valid, valid, valid, (1, 1), None),
        ]
        for out_mean, in_mean, out_std, in_std, signs, areas in cases:
            with self.subTest(signs=signs, areas=areas, out=out_mean.tolist()):
                with self.assertRaises(ValueError):
                    signed_partial_currents(
                        out_mean, in_mean, out_std, in_std, outward_signs=signs, face_areas=areas
                    )

    def test_invalid_geometry_parameters_fail_before_openmc_import(self) -> None:
        for sites, z, pitch, center in [
            ([], [0, 1], 2, (0, 0)),
            ([(0, 0), (0, 0)], [0, 1], 2, (0, 0)),
            ([(0.5, 0)], [0, 1], 2, (0, 0)),
            ([(True, 0)], [0, 1], 2, (0, 0)),
            ([(0, 0, 0)], [0, 1], 2, (0, 0)),
            ([(0, 0)], [0, 0], 2, (0, 0)),
            ([(0, 0)], [0, float("inf")], 2, (0, 0)),
            ([(0, 0)], [0], 2, (0, 0)),
            ([(0, 0)], [0, 1], -2, (0, 0)),
            ([(0, 0)], [0, 1], 2, (0,)),
        ]:
            with self.subTest(sites=sites, z=z, pitch=pitch, center=center):
                with self.assertRaises(ValueError):
                    build_hex_current_surfaces(sites, z, pitch, center)


@unittest.skipIf(openmc is None, "native OpenMC Python package is optional")
class NativeHexCurrentTests(unittest.TestCase):
    def setUp(self) -> None:
        openmc.reset_auto_ids()
        self.registry = build_hex_current_surfaces([(0, 0), (0, 1)], [0, 1, 3], 2)
        self.domains = {
            (i, j, k): openmc.Cell(name=f"node-{i}-{j}-{k}", region=self.registry.node_region((i, j), k))
            for i, j in self.registry.sites
            for k in range(2)
        }
        self.geometry = openmc.Geometry(list(self.domains.values()))

    def bind(self, **kwargs):
        return bind_hex_current_domains(
            self.registry, self.domains, energy_bounds=(0, 1, 1e7), geometry=self.geometry, **kwargs
        )

    def test_planes_shared_across_side_and_axial_neighbors(self) -> None:
        first = self.registry.node_surfaces((0, 0), 0)
        lateral = self.registry.node_surfaces((0, 1), 0)
        axial = self.registry.node_surfaces((0, 0), 1)
        self.assertIs(first[0], lateral[3])
        self.assertIs(first[7], axial[6])
        self.assertEqual(len(set(surface.id for surface in first)), 8)
        self.assertEqual(self.registry.site_center((0, 1)), (0, 2))
        self.assertIn((0, 0, 0.5), self.domains[(0, 0, 0)].region)
        self.assertNotIn((0, 1.01, 0.5), self.domains[(0, 0, 0)].region)

    def test_all_six_neighbor_maps_and_global_offsets(self) -> None:
        sites = [(0, 0), (0, 1), (1, 0), (1, -1), (0, -1), (-1, 0), (-1, 1)]
        registry = build_hex_current_surfaces(sites, [-2, 5], 5.1, center=(3, -7))
        central = registry.side_surfaces((0, 0))
        for side, site in enumerate(sites[1:]):
            self.assertIs(central[side], registry.side_surfaces(site)[(side + 3) % 6])
        self.assertEqual(registry.site_center((0, 0)), (3, -7))
        for site in sites:
            x, y = registry.site_center(site)
            self.assertIn((x, y, 0), registry.node_region(site, 0))

    def test_tallies_and_manifest_are_exact_and_serializable(self) -> None:
        plan = self.bind()
        self.assertEqual(len(plan.tallies), 8)
        json.dumps(plan.manifest, allow_nan=False)
        node = plan.manifest["nodes"][0]
        self.assertEqual([face["name"] for face in node["faces"]], list(FACE_NAMES))
        self.assertEqual([face["outward_sign"] for face in node["faces"]], [1, 1, 1, -1, -1, -1, -1, 1])
        self.assertEqual(node["faces"][0]["neighbor_node"], [0, 1, 0])
        self.assertEqual(node["faces"][0]["neighbor_face_index"], 3)
        self.assertEqual(node["faces"][7]["neighbor_node"], [0, 0, 1])
        self.assertIsNone(node["faces"][6]["neighbor_node"])
        self.assertAlmostEqual(node["faces"][0]["area_cm2"], 2 / math.sqrt(3))
        self.assertAlmostEqual(node["faces"][6]["area_cm2"], 2 * math.sqrt(3))
        self.assertIsNone(plan.manifest["net_std_dev"])
        for tally in plan.tallies:
            self.assertEqual(tally.scores, ["current"])
            self.assertEqual(tally.estimator, "analog")
            self.assertEqual([f.num_bins for f in tally.filters], [8, 1, 2])
        validate_hex_current_tally_xml(plan, plan.tallies.to_xml_element())

    def test_domain_may_translate_its_fine_fill(self) -> None:
        domain = self.domains[(0, 1, 0)]
        domain.fill = openmc.Universe(cells=[openmc.Cell()])
        domain.translation = (0, 2, 0)
        self.bind()

    def test_ancestors_may_not_transform_the_global_geometry(self) -> None:
        parent = openmc.Cell(fill=self.geometry.root_universe)
        parent.translation = (2, 0, 0)
        self.geometry = openmc.Geometry([parent])
        with self.assertRaisesRegex(ValueError, "ancestors may not be translated"):
            self.bind()

    def test_reused_domains_are_rejected(self) -> None:
        inner = self.geometry.root_universe
        self.geometry = openmc.Geometry([openmc.Cell(fill=inner), openmc.Cell(fill=inner)])
        with self.assertRaisesRegex(ValueError, "is reused"):
            self.bind()

    def test_clipped_domains_are_rejected(self) -> None:
        clip = openmc.XPlane(x0=0.5)
        self.geometry = openmc.Geometry([openmc.Cell(region=-clip, fill=self.geometry.root_universe)])
        with self.assertRaisesRegex(ValueError, "clips current domain"):
            self.bind()

    def test_convex_containing_ancestor_is_supported(self) -> None:
        bottom, top = openmc.ZPlane(z0=-2), openmc.ZPlane(z0=6)
        self.geometry = openmc.Geometry(
            [openmc.Cell(region=+bottom & -top, fill=self.geometry.root_universe)]
        )
        self.bind()

    def test_ancestor_vacuum_touching_node_is_rejected(self) -> None:
        bottom = openmc.ZPlane(z0=0, boundary_type="vacuum")
        self.geometry = openmc.Geometry([openmc.Cell(region=+bottom, fill=self.geometry.root_universe)])
        with self.assertRaisesRegex(ValueError, "non-transmission ancestor boundary"):
            self.bind()

    def test_ancestor_duplicate_coincident_face_is_rejected(self) -> None:
        bottom = openmc.ZPlane(z0=0)
        self.geometry = openmc.Geometry([openmc.Cell(region=+bottom, fill=self.geometry.root_universe)])
        with self.assertRaisesRegex(ValueError, "coincident ancestor faces"):
            self.bind()

    def test_ancestor_shared_face_is_supported(self) -> None:
        self.geometry = openmc.Geometry(
            [openmc.Cell(region=+self.registry.z_planes[0], fill=self.geometry.root_universe)]
        )
        self.bind()

    def test_detached_equivalent_surfaces_are_rejected(self) -> None:
        other = build_hex_current_surfaces([(0, 0)], [0, 1], 2)
        self.domains[(0, 0, 0)].region = other.node_region((0, 0), 0)
        with self.assertRaisesRegex(ValueError, "exact eight node halfspaces"):
            self.bind()

    def test_absent_domains_and_non_transmission_are_rejected(self) -> None:
        self.geometry = openmc.Geometry([openmc.Cell()])
        with self.assertRaisesRegex(ValueError, "absent from geometry"):
            self.bind()
        self.geometry = openmc.Geometry(list(self.domains.values()))
        self.registry.z_planes[0].boundary_type = "reflective"
        with self.assertRaisesRegex(ValueError, "transmission"):
            self.bind()

    def test_wrong_surface_sense_is_rejected(self) -> None:
        self.domains[(0, 0, 0)].region = ~self.domains[(0, 0, 0)].region
        with self.assertRaises(ValueError):
            self.bind()

    def test_unknown_nodes_are_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "not registered"):
            self.registry.node_region((3, 4), 0)
        with self.assertRaisesRegex(ValueError, "axial layer"):
            self.registry.node_region((0, 0), 2)

    def test_tampered_xml_is_rejected(self) -> None:
        plan = self.bind()
        root = plan.tallies.to_xml_element()
        root.find("filter[@type='cellfrom']/bins").text = "900001"
        with self.assertRaisesRegex(ValueError, "bins changed"):
            validate_hex_current_tally_xml(plan, root)

    def test_modified_plane_coefficients_are_rejected(self) -> None:
        self.registry.z_planes[0].z0 = -0.5
        with self.assertRaisesRegex(ValueError, "coefficients were modified"):
            self.bind()

    def test_high_cell_ids_remain_exact_without_merging(self) -> None:
        for cell, new_id in zip(self.domains.values(), range(900001, 900005), strict=True):
            cell.id = new_id
        self.geometry = openmc.Geometry(list(self.domains.values()))
        plan = self.bind()
        self.assertEqual([node["cell_id"] for node in plan.manifest["nodes"]], list(range(900001, 900005)))
        validate_hex_current_tally_xml(plan, plan.tallies.to_xml_element())


if __name__ == "__main__":
    unittest.main()
