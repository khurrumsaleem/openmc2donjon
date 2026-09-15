from __future__ import annotations

import contextlib
import io
import json
from pathlib import Path
import sys
import tempfile
import types
import unittest

import h5py
import numpy as np

from openmc2donjon.cli import main as cli_main
from openmc2donjon.constants import MGXS_DONJON_GROUP_ORDER
from openmc2donjon.openmc_volume_flux import (
    DATASET_NAME,
    DEFAULT_SOURCE_GROUP_ORDER,
    SCHEMA,
    STD_DEV_DATASET_NAME,
    export_openmc_volume_flux,
    reverse_openmc_energy_filter_flux,
    write_openmc_flux_hdf5,
    write_openmc_volume_flux_hdf5,
)


class OpenMCVolumeFluxTests(unittest.TestCase):
    def test_reverses_openmc_energy_filter_order(self) -> None:
        raw = np.array(
            [
                [[1.0], [2.0], [3.0]],
                [[4.0], [5.0], [6.0]],
            ]
        )

        values = reverse_openmc_energy_filter_flux(
            raw,
            mixture_count=2,
            energy_groups=3,
        )

        np.testing.assert_allclose(values, [[3.0, 2.0, 1.0], [6.0, 5.0, 4.0]])

    def test_writes_canonical_openmc_volume_flux_dataset(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "mgxs.h5"
            flux = np.array([[10.0, 20.0], [30.0, 40.0]])
            std_dev = np.array([[0.1, 0.2], [0.3, 0.4]])

            report = write_openmc_volume_flux_hdf5(
                path,
                flux,
                mixture_names=("ASM_A", "ASM_B"),
                std_dev=std_dev,
            )

            with h5py.File(path, "r") as h5:
                dataset = h5[DATASET_NAME]
                values = dataset[:]
                attrs = dict(dataset.attrs)
                std_dataset = h5[STD_DEV_DATASET_NAME]
                std_values = std_dataset[:]
                std_attrs = dict(std_dataset.attrs)

        self.assertEqual(report.dataset, DATASET_NAME)
        self.assertEqual(report.std_dev_dataset, STD_DEV_DATASET_NAME)
        self.assertEqual(report.mixture_names, ("ASM_A", "ASM_B"))
        self.assertEqual(report.energy_groups, 2)
        self.assertEqual(report.minimum, 10.0)
        self.assertEqual(report.maximum, 40.0)
        self.assertAlmostEqual(report.max_relative_std_dev or 0.0, 0.01)
        np.testing.assert_allclose(values, flux)
        np.testing.assert_allclose(std_values, std_dev)
        self.assertEqual(attrs["schema"], SCHEMA)
        self.assertEqual(attrs["group_order"], MGXS_DONJON_GROUP_ORDER)
        self.assertEqual(attrs["source_group_order"], DEFAULT_SOURCE_GROUP_ORDER)
        self.assertEqual(attrs["layout"], "[mixture, group]")
        self.assertFalse(attrs["energy_bounds_verified"])
        self.assertFalse(attrs["spatial_domain_order_verified"])
        self.assertEqual(
            tuple(_decode(value) for value in attrs["mixture_names"]),
            ("ASM_A", "ASM_B"),
        )
        self.assertEqual(std_attrs["schema"], SCHEMA)
        self.assertEqual(std_attrs["group_order"], MGXS_DONJON_GROUP_ORDER)
        self.assertEqual(std_attrs["source_group_order"], DEFAULT_SOURCE_GROUP_ORDER)
        self.assertFalse(std_attrs["energy_bounds_verified"])
        self.assertFalse(std_attrs["spatial_domain_order_verified"])

    def test_writes_custom_openmc_flux_dataset_for_mg_macro_flux(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "mg_flux.h5"
            flux = np.array([[1.0, 2.0], [3.0, 4.0]])

            report = write_openmc_flux_hdf5(
                path,
                flux,
                mixture_names=("FUEL", "MOD"),
                dataset_name="openmc_mg_flux",
                std_dev_dataset_name="openmc_mg_flux_std_dev",
            )

            with h5py.File(path, "r") as h5:
                dataset = h5["openmc_mg_flux"]
                values = dataset[:]
                attrs = dict(dataset.attrs)

        self.assertEqual(report.dataset, "openmc_mg_flux")
        self.assertIsNone(report.std_dev_dataset)
        self.assertFalse(report.energy_bounds_verified)
        self.assertFalse(report.spatial_domain_order_verified)
        np.testing.assert_allclose(values, flux)
        self.assertEqual(attrs["group_order"], MGXS_DONJON_GROUP_ORDER)
        self.assertEqual(
            tuple(_decode(value) for value in attrs["mixture_names"]),
            ("FUEL", "MOD"),
        )

    def test_exports_statepoint_tally_to_openmc_flux_contract(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            statepoint = tmp / "statepoint.10.h5"
            statepoint.write_bytes(b"fake")
            mgxs = tmp / "mgxs.h5"
            output = tmp / "ce_flux.h5"
            summary = tmp / "ce_flux_summary.json"
            _write_mgxs_metadata(mgxs)
            fake_openmc = _fake_openmc_module(
                mean=np.array([[1.0], [2.0], [3.0], [4.0]]),
                std_dev=np.array([[0.1], [0.2], [0.3], [0.4]]),
            )

            stream = io.StringIO()
            with _patched_openmc(fake_openmc), contextlib.redirect_stdout(stream):
                report = export_openmc_volume_flux(
                    statepoint,
                    output,
                    mgxs_h5=mgxs,
                    tally_name="ce_flux",
                    dataset_name="openmc_ce_flux",
                    summary_json=summary,
                )

            payload = json.loads(summary.read_text(encoding="utf-8"))
            with h5py.File(output, "r") as h5:
                values = h5["openmc_ce_flux"][:]
                std_values = h5["openmc_ce_flux_std_dev"][:]
                attrs = dict(h5["openmc_ce_flux"].attrs)

        self.assertEqual(report.dataset, "openmc_ce_flux")
        self.assertEqual(report.std_dev_dataset, "openmc_ce_flux_std_dev")
        self.assertEqual(report.statepoint, statepoint)
        self.assertEqual(report.tally_name, "ce_flux")
        self.assertIn("openmc2donjon_volume_flux_export_passed", stream.getvalue())
        np.testing.assert_allclose(values, [[2.0, 1.0], [4.0, 3.0]])
        np.testing.assert_allclose(std_values, [[0.2, 0.1], [0.4, 0.3]])
        self.assertEqual(attrs["group_order"], MGXS_DONJON_GROUP_ORDER)
        self.assertEqual(payload["decision"], "openmc2donjon_volume_flux_export_passed")
        self.assertEqual(payload["dataset"], "openmc_ce_flux")
        self.assertEqual(payload["tally_name"], "ce_flux")
        self.assertTrue(payload["energy_bounds_verified"])
        self.assertTrue(payload["spatial_domain_order_verified"])
        self.assertEqual(payload["source_filter_order"], ["CellFilter", "EnergyFilter"])
        self.assertEqual(payload["source_domain_ids"], [101, 202])
        self.assertTrue(attrs["energy_bounds_verified"])
        self.assertTrue(attrs["spatial_domain_order_verified"])
        self.assertEqual(
            tuple(_decode(value) for value in attrs["source_filter_order"]),
            ("CellFilter", "EnergyFilter"),
        )
        self.assertEqual(tuple(attrs["source_domain_ids"]), (101, 202))

    def test_canonicalizes_reversed_energy_and_spatial_filter_axes(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            statepoint = tmp / "statepoint.10.h5"
            statepoint.write_bytes(b"fake")
            mgxs = tmp / "mgxs.h5"
            output = tmp / "ce_flux.h5"
            _write_mgxs_metadata(mgxs)
            # Filter order is [energy, cell].  In that order, the low-to-high
            # matrix is [[FUEL-g0, MOD-g0], [FUEL-g1, MOD-g1]].
            fake_openmc = _fake_openmc_module(
                mean=np.array([1.0, 3.0, 2.0, 4.0]),
                std_dev=np.array([0.1, 0.3, 0.2, 0.4]),
                filters=(
                    EnergyFilter(((1.0e-5, 1.0), (1.0, 1.0e7))),
                    CellFilter((101, 202)),
                ),
            )

            with _patched_openmc(fake_openmc):
                report = export_openmc_volume_flux(
                    statepoint,
                    output,
                    mgxs_h5=mgxs,
                    tally_name="ce_flux",
                )

            with h5py.File(output, "r") as h5:
                values = h5[DATASET_NAME][:]

        np.testing.assert_allclose(values, [[2.0, 1.0], [4.0, 3.0]])
        self.assertEqual(report.source_filter_order, ("EnergyFilter", "CellFilter"))
        self.assertTrue(report.energy_bounds_verified)
        self.assertTrue(report.spatial_domain_order_verified)

    def test_rejects_energy_filter_bounds_that_disagree_with_mgxs(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            statepoint = tmp / "statepoint.10.h5"
            statepoint.write_bytes(b"fake")
            mgxs = tmp / "mgxs.h5"
            _write_mgxs_metadata(mgxs)
            fake_openmc = _fake_openmc_module(
                mean=np.arange(1.0, 5.0),
                std_dev=np.full(4, 0.1),
                filters=(
                    CellFilter((101, 202)),
                    EnergyFilter(((1.0e-5, 2.0), (2.0, 1.0e7))),
                ),
            )

            with _patched_openmc(fake_openmc), self.assertRaisesRegex(
                ValueError,
                "EnergyFilter bin edges do not match",
            ):
                export_openmc_volume_flux(
                    statepoint,
                    tmp / "out.h5",
                    mgxs_h5=mgxs,
                    tally_name="ce_flux",
                )

    def test_rejects_swapped_spatial_domain_ids(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            statepoint = tmp / "statepoint.10.h5"
            statepoint.write_bytes(b"fake")
            mgxs = tmp / "mgxs.h5"
            _write_mgxs_metadata(mgxs)
            fake_openmc = _fake_openmc_module(
                mean=np.arange(1.0, 5.0),
                std_dev=np.full(4, 0.1),
                filters=(
                    CellFilter((202, 101)),
                    EnergyFilter(((1.0e-5, 1.0), (1.0, 1.0e7))),
                ),
            )

            with _patched_openmc(fake_openmc), self.assertRaisesRegex(
                ValueError,
                "does not match --mgxs source_domain_id order",
            ):
                export_openmc_volume_flux(
                    statepoint,
                    tmp / "out.h5",
                    mgxs_h5=mgxs,
                    tally_name="ce_flux",
                )

    def test_explicit_source_domain_ids_support_a_different_tally_geometry(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            statepoint = tmp / "statepoint.10.h5"
            statepoint.write_bytes(b"fake")
            mgxs = tmp / "mgxs.h5"
            output = tmp / "ce_flux.h5"
            _write_mgxs_metadata(mgxs)
            fake_openmc = _fake_openmc_module(
                mean=np.arange(1.0, 5.0),
                std_dev=np.full(4, 0.1),
                filters=(
                    CellFilter((7001, 7002)),
                    EnergyFilter(((1.0e-5, 1.0), (1.0, 1.0e7))),
                ),
            )

            with _patched_openmc(fake_openmc):
                report = export_openmc_volume_flux(
                    statepoint,
                    output,
                    mgxs_h5=mgxs,
                    tally_name="ce_flux",
                    source_domain_ids=(7001, 7002),
                )

            with h5py.File(output, "r") as h5:
                attrs = dict(h5[DATASET_NAME].attrs)

        self.assertEqual(report.source_domain_ids, (7001, 7002))
        self.assertEqual(tuple(attrs["source_domain_ids"]), (7001, 7002))
        self.assertEqual(
            tuple(_decode(value) for value in attrs["mixture_names"]),
            ("FUEL", "MOD"),
        )

    def test_rejects_missing_or_unverifiable_tally_filters_with_mgxs(self) -> None:
        cases = (
            ((), "exactly one EnergyFilter"),
            (
                (
                    CellFilter((101, 202)),
                    MaterialFilter((101, 202)),
                ),
                "exactly one EnergyFilter",
            ),
            (
                (
                    MeshFilter(((1, 1), (2, 1))),
                    EnergyFilter(((1.0e-5, 1.0), (1.0, 1.0e7))),
                ),
                "supported spatial/domain filter",
            ),
        )
        for filters, message in cases:
            with self.subTest(filters=tuple(type(value).__name__ for value in filters)):
                with tempfile.TemporaryDirectory() as tmpdir:
                    tmp = Path(tmpdir)
                    statepoint = tmp / "statepoint.10.h5"
                    statepoint.write_bytes(b"fake")
                    mgxs = tmp / "mgxs.h5"
                    _write_mgxs_metadata(mgxs)
                    fake_openmc = _fake_openmc_module(
                        mean=np.arange(1.0, 5.0),
                        std_dev=np.full(4, 0.1),
                        filters=filters,
                    )

                    with _patched_openmc(fake_openmc), self.assertRaisesRegex(
                        ValueError,
                        message,
                    ):
                        export_openmc_volume_flux(
                            statepoint,
                            tmp / "out.h5",
                            mgxs_h5=mgxs,
                            tally_name="ce_flux",
                        )

    def test_rejects_mgxs_metadata_overrides_that_disagree(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            statepoint = tmp / "statepoint.10.h5"
            statepoint.write_bytes(b"fake")
            mgxs = tmp / "mgxs.h5"
            _write_mgxs_metadata(mgxs)
            fake_openmc = _fake_openmc_module(
                mean=np.arange(1.0, 5.0),
                std_dev=np.full(4, 0.1),
            )

            with _patched_openmc(fake_openmc), self.assertRaisesRegex(
                ValueError,
                "--mixture-names disagrees",
            ):
                export_openmc_volume_flux(
                    statepoint,
                    tmp / "names.h5",
                    mgxs_h5=mgxs,
                    tally_name="ce_flux",
                    mixture_names=("MOD", "FUEL"),
                )
            with _patched_openmc(fake_openmc), self.assertRaisesRegex(
                ValueError,
                "--energy-groups disagrees",
            ):
                export_openmc_volume_flux(
                    statepoint,
                    tmp / "groups.h5",
                    mgxs_h5=mgxs,
                    tally_name="ce_flux",
                    energy_groups=3,
                )

    def test_rejects_nonpositive_or_mismatched_values(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "mgxs.h5"

            with self.assertRaisesRegex(ValueError, "must be positive"):
                write_openmc_volume_flux_hdf5(
                    path,
                    [[1.0, 0.0]],
                    mixture_names=("ASM_A",),
                )

            with self.assertRaisesRegex(ValueError, "mixture axis"):
                write_openmc_volume_flux_hdf5(
                    path,
                    [[1.0, 2.0]],
                    mixture_names=("ASM_A", "ASM_B"),
                )

            with self.assertRaisesRegex(ValueError, "non-negative"):
                write_openmc_volume_flux_hdf5(
                    path,
                    [[1.0, 2.0]],
                    mixture_names=("ASM_A",),
                    std_dev=[[0.0, -0.1]],
                )

    def test_allow_zero_accepts_exactly_zero_flux_bins(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "mgxs.h5"
            flux = np.array([[10.0, 0.0], [30.0, 40.0]])
            std_dev = np.array([[0.1, 0.0], [0.3, 0.4]])

            report = write_openmc_volume_flux_hdf5(
                path,
                flux,
                mixture_names=("ASM_A", "ASM_B"),
                std_dev=std_dev,
                allow_zero=True,
            )

            with h5py.File(path, "r") as h5:
                values = h5[DATASET_NAME][:]

        self.assertEqual(report.minimum, 0.0)
        self.assertEqual(report.maximum, 40.0)
        self.assertTrue(report.allow_zero_flux)
        self.assertAlmostEqual(report.max_relative_std_dev or 0.0, 0.01)
        np.testing.assert_allclose(values, flux)

    def test_allow_zero_still_rejects_negative_or_nonfinite_values(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "mgxs.h5"

            with self.assertRaisesRegex(ValueError, "must be non-negative"):
                write_openmc_volume_flux_hdf5(
                    path,
                    [[1.0, -2.0]],
                    mixture_names=("ASM_A",),
                    allow_zero=True,
                )

            with self.assertRaisesRegex(ValueError, "must be finite"):
                write_openmc_volume_flux_hdf5(
                    path,
                    [[1.0, np.nan]],
                    mixture_names=("ASM_A",),
                    allow_zero=True,
                )

    def test_cli_allow_zero_flux_flag_threads_through_export(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            statepoint = tmp / "statepoint.10.h5"
            statepoint.write_bytes(b"fake")
            mgxs = tmp / "mgxs.h5"
            output = tmp / "ce_flux.h5"
            summary = tmp / "ce_flux_summary.json"
            _write_mgxs_metadata(mgxs)
            fake_openmc = _fake_openmc_module(
                mean=np.array([[0.0], [1.0], [0.0], [2.0]]),
                std_dev=np.array([[0.0], [0.1], [0.0], [0.2]]),
                filters=(
                    CellFilter((7001, 7002)),
                    EnergyFilter(((1.0e-5, 1.0), (1.0, 1.0e7))),
                ),
            )

            stream = io.StringIO()
            with _patched_openmc(fake_openmc), contextlib.redirect_stdout(stream):
                rc = cli_main(
                    [
                        "export-volume-flux",
                        str(statepoint),
                        "-o",
                        str(output),
                        "--mgxs",
                        str(mgxs),
                        "--tally-name",
                        "ce_flux",
                        "--source-domain-ids",
                        "7001,7002",
                        "--allow-zero-flux",
                        "--summary-json",
                        str(summary),
                    ]
                )

            payload = json.loads(summary.read_text(encoding="utf-8"))
            with h5py.File(output, "r") as h5:
                values = h5[DATASET_NAME][:]

        self.assertEqual(rc, 0)
        self.assertTrue(payload["allow_zero_flux"])
        self.assertEqual(payload["min"], 0.0)
        self.assertEqual(payload["source_domain_ids"], [7001, 7002])
        np.testing.assert_allclose(values, [[1.0, 0.0], [2.0, 0.0]])


def _decode(value: object) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def _write_mgxs_metadata(path: Path) -> None:
    with h5py.File(path, "w") as h5:
        h5.attrs["energy_groups"] = 2
        h5.create_dataset("energy_bounds", data=np.array([1.0e-5, 1.0, 1.0e7]))
        h5.create_dataset("mixture_names", data=np.asarray(("FUEL", "MOD"), dtype="S"))
        mixtures = h5.create_group("mixtures")
        fuel = mixtures.create_group("FUEL")
        fuel.attrs["source_domain_id"] = 101
        moderator = mixtures.create_group("MOD")
        moderator.attrs["source_domain_id"] = 202


class EnergyFilter:
    def __init__(self, bins: tuple[tuple[float, float], ...]) -> None:
        self.bins = bins
        self.num_bins = len(bins)


class CellFilter:
    def __init__(self, bins: tuple[int, ...]) -> None:
        self.bins = bins
        self.num_bins = len(bins)


class MaterialFilter(CellFilter):
    pass


class MeshFilter:
    def __init__(self, bins: tuple[tuple[int, int], ...]) -> None:
        self.bins = bins
        self.num_bins = len(bins)


def _fake_openmc_module(
    *,
    mean: np.ndarray,
    std_dev: np.ndarray,
    filters: tuple[object, ...] | None = None,
):
    selected_filters = (
        filters
        if filters is not None
        else (
            CellFilter((101, 202)),
            EnergyFilter(((1.0e-5, 1.0), (1.0, 1.0e7))),
        )
    )

    class FakeTally:
        def __init__(self) -> None:
            self.filters = selected_filters

        def get_values(self, *, scores=None, value: str = "mean"):
            if scores != ["flux"]:
                raise AssertionError(f"unexpected scores: {scores!r}")
            if value == "mean":
                return mean
            if value == "std_dev":
                return std_dev
            raise AssertionError(f"unexpected value selector: {value!r}")

    class FakeStatePoint:
        def __init__(self, path: str) -> None:
            self.path = path

        def __enter__(self):
            return self

        def __exit__(self, *_exc: object) -> None:
            return None

        def get_tally(self, *, name: str):
            if name != "ce_flux":
                raise AssertionError(f"unexpected tally name: {name!r}")
            return FakeTally()

    return types.SimpleNamespace(StatePoint=FakeStatePoint)


@contextlib.contextmanager
def _patched_openmc(fake_openmc):
    previous = sys.modules.get("openmc")
    sys.modules["openmc"] = fake_openmc
    try:
        yield
    finally:
        if previous is None:
            sys.modules.pop("openmc", None)
        else:
            sys.modules["openmc"] = previous


if __name__ == "__main__":
    unittest.main()
