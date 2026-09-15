from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

import h5py
import numpy as np

from openmc2donjon.mixture_collapse import collapse_components


class MixtureCollapseTests(unittest.TestCase):
    def test_preserves_vector_and_scatter_reaction_rates(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "source.h5"
            output = Path(tmp) / "components.h5"
            _write_source(source)

            collapse_components(
                source,
                output,
                groups=(("CENTER", ("A",)), ("RING", ("B", "C"))),
            )

            with h5py.File(source, "r") as original, h5py.File(output, "r") as collapsed:
                flux = np.asarray(original["openmc_volume_flux"])[1:]
                total = np.stack(
                    [np.asarray(original["mixtures"][name]["total"]) for name in ("B", "C")]
                )
                reduced_absorption = np.stack(
                    [
                        np.asarray(original["mixtures"][name]["reduced_absorption"])
                        for name in ("B", "C")
                    ]
                )
                scatter = np.stack(
                    [
                        np.asarray(original["mixtures"][name]["scatter_matrix"])
                        for name in ("B", "C")
                    ]
                )
                ring_flux = np.asarray(collapsed["openmc_volume_flux"])[1]
                ring = collapsed["mixtures"]["RING"]

                np.testing.assert_allclose(ring_flux, flux.sum(axis=0))
                np.testing.assert_allclose(
                    np.asarray(ring["total"]) * ring_flux,
                    np.sum(total * flux, axis=0),
                )
                np.testing.assert_allclose(
                    np.asarray(ring["reduced_absorption"]) * ring_flux,
                    np.sum(reduced_absorption * flux, axis=0),
                )
                np.testing.assert_allclose(
                    np.asarray(ring["scatter_matrix"])[0] * ring_flux[:, np.newaxis],
                    np.sum(scatter[:, 0] * flux[:, :, np.newaxis], axis=0),
                )
                expected_total_std = np.sum(0.01 * total * flux, axis=0) / ring_flux
                np.testing.assert_allclose(
                    np.asarray(ring["total_std_dev"]),
                    expected_total_std,
                )
                self.assertIn("scatter_matrix_std_dev", ring)
                self.assertIn("chi_std_dev", ring)
                self.assertEqual(
                    ring["total_std_dev"].attrs["component_uncertainty_method"],
                    "conservative-l1-source-xs-bound-no-covariance",
                )
                self.assertEqual(float(ring.attrs["volume"]), 5.0)
                self.assertEqual(tuple(collapsed["mixture_names"].asstr()[:]), ("CENTER", "RING"))
                self.assertFalse(
                    bool(collapsed.attrs["openmc_scatter_multiplicity_weighted"])
                )
                self.assertEqual(
                    collapsed.attrs["openmc_scatter_balance_dataset"],
                    "absorption",
                )
                self.assertFalse(
                    bool(ring.attrs["openmc_scatter_multiplicity_weighted"])
                )
                self.assertEqual(
                    ring.attrs["openmc_scatter_balance_dataset"],
                    "absorption",
                )

    def test_requires_exact_source_coverage(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "source.h5"
            _write_source(source)
            with self.assertRaisesRegex(ValueError, "exactly once"):
                collapse_components(
                    source,
                    Path(tmp) / "bad.h5",
                    groups=(("ONLY", ("A", "B")),),
                )

    def test_rejects_mixed_scatter_contracts_before_writing(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "source.h5"
            output = Path(tmp) / "components.h5"
            _write_source(source)
            with h5py.File(source, "a") as h5:
                ordinary = h5["mixtures/A"].attrs
                ordinary["openmc_scatter_mgxs_type"] = "scatter matrix"
                ordinary["openmc_scatter_multiplicity_weighted"] = False
                ordinary["openmc_scatter_balance_dataset"] = "absorption"
                ordinary["openmc_transport_mgxs_type"] = "transport"
                for name in ("B", "C"):
                    weighted = h5[f"mixtures/{name}"].attrs
                    weighted["openmc_scatter_mgxs_type"] = (
                        "consistent nu-scatter matrix"
                    )
                    weighted["openmc_scatter_multiplicity_weighted"] = True
                    weighted["openmc_scatter_balance_dataset"] = (
                        "reduced_absorption"
                    )
                    weighted["openmc_transport_mgxs_type"] = "nu-transport"

            with self.assertRaisesRegex(
                ValueError,
                "one coherent OpenMC scatter/removal contract",
            ):
                collapse_components(
                    source,
                    output,
                    groups=(("CENTER", ("A",)), ("RING", ("B", "C"))),
                )
            self.assertFalse(output.exists())

    def test_rejects_nu_scatter_missing_reduced_absorption_before_writing(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "source.h5"
            output = Path(tmp) / "components.h5"
            _write_source(source)
            with h5py.File(source, "a") as h5:
                h5.attrs["openmc_scatter_mgxs_type"] = (
                    "consistent nu-scatter matrix"
                )
                h5.attrs["openmc_scatter_multiplicity_weighted"] = True
                h5.attrs["openmc_scatter_balance_dataset"] = "reduced_absorption"
                h5.attrs["openmc_transport_mgxs_type"] = "nu-transport"
                del h5["mixtures/B/reduced_absorption"]

            with self.assertRaisesRegex(
                ValueError,
                "requires a finite reduced_absorption vector",
            ):
                collapse_components(
                    source,
                    output,
                    groups=(("CENTER", ("A",)), ("RING", ("B", "C"))),
                )
            self.assertFalse(output.exists())

    def test_writes_one_canonical_nu_scatter_contract(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "source.h5"
            output = Path(tmp) / "components.h5"
            _write_source(source)
            with h5py.File(source, "a") as h5:
                for attrs in (
                    h5.attrs,
                    *(group.attrs for group in h5["mixtures"].values()),
                ):
                    attrs["openmc_scatter_mgxs_type"] = (
                        "consistent_nu_scatter_matrix"
                    )
                    attrs["openmc_scatter_multiplicity_weighted"] = True
                    attrs["openmc_scatter_balance_dataset"] = (
                        "reduced_absorption"
                    )
                    attrs["openmc_transport_mgxs_type"] = "nu-transport"

            collapse_components(
                source,
                output,
                groups=(("CENTER", ("A",)), ("RING", ("B", "C"))),
            )

            with h5py.File(output, "r") as collapsed:
                expected = "consistent nu-scatter matrix"
                self.assertEqual(
                    collapsed.attrs["openmc_scatter_mgxs_type"], expected
                )
                self.assertTrue(
                    bool(collapsed.attrs["openmc_scatter_multiplicity_weighted"])
                )
                self.assertEqual(
                    collapsed.attrs["openmc_scatter_balance_dataset"],
                    "reduced_absorption",
                )
                self.assertEqual(
                    collapsed.attrs["openmc_transport_mgxs_type"],
                    "nu-transport",
                )
                for group in collapsed["mixtures"].values():
                    self.assertEqual(
                        group.attrs["openmc_scatter_mgxs_type"], expected
                    )
                    self.assertTrue(
                        bool(
                            group.attrs[
                                "openmc_scatter_multiplicity_weighted"
                            ]
                        )
                    )
                    self.assertEqual(
                        group.attrs["openmc_scatter_balance_dataset"],
                        "reduced_absorption",
                    )
                    self.assertEqual(
                        group.attrs["openmc_transport_mgxs_type"],
                        "nu-transport",
                    )

    def test_rejects_mixed_declared_and_inferred_transport_contracts(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "source.h5"
            output = Path(tmp) / "components.h5"
            _write_source(source)
            with h5py.File(source, "a") as h5:
                h5["mixtures/A"].attrs["openmc_transport_mgxs_type"] = (
                    "transport"
                )

            with self.assertRaisesRegex(
                ValueError,
                "consistently declared or legacy-inferred",
            ):
                collapse_components(
                    source,
                    output,
                    groups=(("CENTER", ("A",)), ("RING", ("B", "C"))),
                )
            self.assertFalse(output.exists())


def _write_source(path: Path) -> None:
    names = ("A", "B", "C")
    flux = np.asarray([[1.0, 2.0], [2.0, 4.0], [3.0, 5.0]])
    with h5py.File(path, "w") as h5:
        h5.attrs["energy_groups"] = 2
        h5.attrs["legendre_order"] = 0
        h5.create_dataset("energy_bounds", data=[1.0e-5, 1.0, 1.0e7])
        h5.create_dataset("mixture_names", data=np.asarray(names, dtype="S"))
        flux_ds = h5.create_dataset("openmc_volume_flux", data=flux)
        flux_ds.attrs["mixture_names"] = np.asarray(names, dtype="S")
        flux_std = h5.create_dataset("openmc_volume_flux_std_dev", data=0.1 * flux)
        flux_std.attrs["mixture_names"] = np.asarray(names, dtype="S")
        mixtures = h5.create_group("mixtures")
        for index, name in enumerate(names):
            group = mixtures.create_group(name)
            group.attrs["volume"] = float(index + 1)
            group.attrs["source_domain_index"] = index + 1
            group.attrs["fissionable"] = True
            total = np.asarray([1.0 + index, 2.0 + index])
            for dataset in (
                "total",
                "transport_total",
                "absorption",
                "reduced_absorption",
                "fission",
                "nu_fission",
                "kappa_fission",
            ):
                group.create_dataset(dataset, data=total)
                group.create_dataset(f"{dataset}_std_dev", data=0.01 * total)
            group.create_dataset("chi", data=[1.0, 0.0])
            group.create_dataset("chi_std_dev", data=[0.01, 0.0])
            group.create_dataset(
                "scatter_matrix",
                data=[[[0.1 + index, 0.2], [0.0, 0.3 + index]]],
            )
            group.create_dataset(
                "scatter_matrix_std_dev",
                data=[[[0.01, 0.01], [0.0, 0.01]]],
            )
