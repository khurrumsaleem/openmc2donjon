from __future__ import annotations

import inspect
import json
from pathlib import Path
import tempfile
import unittest

import h5py
import numpy as np

from openmc2donjon.cli import main as cli_main
from openmc2donjon.commands.sph import (
    build_make_openmc_sph_sidecar_parser,
    build_make_sph_update_table_parser,
)
from openmc2donjon.openmc_provenance import file_sha256
from openmc2donjon.openmc_sph_sidecar import create_openmc_sph_sidecar
from openmc2donjon.sph_iteration import create_sph_update_table


class SphIterationTests(unittest.TestCase):
    def test_openmc_sph_entrypoints_share_production_defaults(self) -> None:
        sidecar_args = build_make_openmc_sph_sidecar_parser().parse_args(
            [
                "mgxs.h5",
                "-o",
                "sph.h5",
                "--reference-flux",
                "ce.h5",
                "--mg-flux",
                "mg.h5",
            ]
        )
        table_args = build_make_sph_update_table_parser().parse_args(
            [
                "mgxs.h5",
                "-o",
                "sph.csv",
                "--reference-flux",
                "ce.h5",
                "--low-order-flux",
                "mg.h5",
            ]
        )
        self.assertEqual(sidecar_args.sph_target, "rate")
        self.assertEqual(sidecar_args.flux_normalization, "auto")
        self.assertEqual(table_args.sph_target, "rate")
        self.assertEqual(table_args.flux_normalization, "auto")

        for function in (create_openmc_sph_sidecar, create_sph_update_table):
            parameters = inspect.signature(function).parameters
            self.assertEqual(parameters["sph_target"].default, "rate")
            self.assertEqual(parameters["flux_normalization"].default, "auto")

    def test_rate_sph_can_pool_a_declared_mixture_symmetry_class(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "mgxs.h5"
            reference_flux = root / "reference_flux.csv"
            mg_flux = root / "mg_flux.csv"
            sidecar = root / "sph.h5"
            summary = root / "summary.json"
            write_mgxs(mgxs)
            reference_flux.write_text(
                "mixture,group,flux\n"
                "fuel,1,2.0\nfuel,2,1.0\n"
                "moderator,1,4.0\nmoderator,2,3.0\n",
                encoding="utf-8",
            )
            mg_flux.write_text(
                "mixture,group,flux\n"
                "fuel,1,3.0\nfuel,2,1.0\n"
                "moderator,1,3.0\nmoderator,2,7.0\n",
                encoding="utf-8",
            )

            self.assertEqual(
                cli_main(
                    [
                        "make-openmc-sph-sidecar",
                        str(mgxs),
                        "-o",
                        str(sidecar),
                        "--reference-flux",
                        str(reference_flux),
                        "--mg-flux",
                        str(mg_flux),
                        "--sph-target",
                        "rate",
                        "--flux-normalization",
                        "none",
                        "--damping",
                        "0.5",
                        "--tie-mixtures",
                        "fuel,moderator",
                        "--summary-json",
                        str(summary),
                    ]
                ),
                0,
            )

            with h5py.File(sidecar, "r") as h5:
                np.testing.assert_allclose(
                    h5["sph"][:],
                    np.asarray([[1.0, np.sqrt(2.0)], [1.0, np.sqrt(2.0)]]),
                )
                self.assertEqual(int(h5.attrs["sph_tied_bin_count"]), 4)
            payload = json.loads(summary.read_text(encoding="utf-8"))
            self.assertEqual(payload["tie_mixture_groups"], [["fuel", "moderator"]])
            self.assertEqual(payload["tied_bin_count"], 4)

    def test_cli_builds_openmc_sph_sidecar_in_one_step(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "mgxs.h5"
            reference_flux = root / "openmc_ce_flux.csv"
            mg_flux = root / "openmc_mg_flux.h5"
            sidecar = root / "openmc_sph.h5"
            table = root / "openmc_sph.csv"
            summary = root / "openmc_sph_summary.json"
            write_mgxs(
                mgxs,
                h_factor={"fuel": np.asarray([19.0, 21.0])},
            )
            reference_flux.write_text(
                "mixture,group,reference_flux\n"
                "fuel,1,1.21\nfuel,2,0.81\n"
                "moderator,1,0.64\nmoderator,2,1.44\n",
                encoding="utf-8",
            )
            with h5py.File(mg_flux, "w") as h5:
                dataset = h5.create_dataset("openmc_mg_flux", data=np.ones((2, 2)))
                dataset.attrs["mixture_names"] = np.asarray(("fuel", "moderator"), dtype="S")
                dataset.attrs["group_order"] = "mgxs_donjon"

            self.assertEqual(
                cli_main(
                    [
                        "make-openmc-sph-sidecar",
                        str(mgxs),
                        "-o",
                        str(sidecar),
                        "--reference-flux",
                        str(reference_flux),
                        "--mg-flux",
                        f"{mg_flux}::openmc_mg_flux",
                        "--table-output",
                        str(table),
                        "--damping",
                        "0.5",
                        "--summary-json",
                        str(summary),
                    ]
                ),
                0,
            )

            expected = np.array(
                [
                    [1.0 / np.sqrt(1.21), 1.0 / np.sqrt(0.81)],
                    [1.0 / np.sqrt(0.64), 1.0 / np.sqrt(1.44)],
                ]
            )
            self.assertTrue(table.exists())
            with h5py.File(sidecar, "r") as h5:
                self.assertEqual(h5.attrs["sph_kind"], "openmc-ce-mg")
                self.assertEqual(bool(h5.attrs["sph_real"]), True)
                self.assertEqual(bool(h5.attrs["sph_applied"]), False)
                self.assertEqual(h5.attrs["source_table"], str(table))
                self.assertEqual(h5.attrs["sph_target"], "rate")
                self.assertEqual(h5.attrs["sph_flux_normalization"], "power")
                self.assertEqual(
                    h5.attrs["sph_derivation"],
                    "rate-preserving-ce-mg-fixed-point",
                )
                self.assertEqual(h5.attrs["sph_input_h5_sha256"], file_sha256(mgxs))
                self.assertEqual(
                    h5.attrs["sph_reference_flux_sha256"],
                    file_sha256(reference_flux),
                )
                self.assertEqual(h5.attrs["sph_mg_flux_sha256"], file_sha256(mg_flux))
                self.assertFalse(bool(h5.attrs["sph_previous_sph_used"]))
                np.testing.assert_allclose(h5["sph"][:], expected, rtol=1.0e-11)
            payload = json.loads(summary.read_text(encoding="utf-8"))
            self.assertEqual(payload["schema"], "openmc2donjon.openmc-sph-sidecar.v1")
            self.assertEqual(payload["decision"], "openmc2donjon_openmc_sph_sidecar_passed")
            self.assertEqual(payload["output_h5"], str(sidecar))
            self.assertEqual(payload["output_table"], str(table))
            self.assertEqual(payload["mg_flux_dataset"], "openmc_mg_flux")
            self.assertEqual(payload["source_label"], "openmc-ce-mg-sph")
            self.assertEqual(payload["sph_target"], "rate")
            self.assertEqual(payload["flux_normalization"], "power")
            self.assertIn("openmc_ce_reference_flux", payload["formula"])

    def test_openmc_sph_sidecar_records_flux_uncertainty(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "mgxs.h5"
            reference_flux = root / "openmc_ce_flux.h5"
            mg_flux = root / "openmc_mg_flux.h5"
            sidecar = root / "openmc_sph.h5"
            table = root / "openmc_sph.csv"
            summary = root / "openmc_sph_summary.json"
            write_mgxs(mgxs)
            _write_flux_source(
                reference_flux,
                "openmc_volume_flux",
                values=np.array([[2.0, 1.0], [4.0, 2.0]]),
                std_dev=np.array([[0.02, 0.02], [0.08, 0.02]]),
            )
            _write_flux_source(
                mg_flux,
                "openmc_mg_flux",
                values=np.array([[1.0, 1.0], [1.0, 1.0]]),
                std_dev=np.array([[0.01, 0.02], [0.03, 0.04]]),
            )

            self.assertEqual(
                cli_main(
                    [
                        "make-openmc-sph-sidecar",
                        str(mgxs),
                        "-o",
                        str(sidecar),
                        "--reference-flux",
                        f"{reference_flux}::openmc_volume_flux",
                        "--mg-flux",
                        f"{mg_flux}::openmc_mg_flux",
                        "--table-output",
                        str(table),
                        "--flux-normalization",
                        "none",
                        "--require-reference-flux-std-dev",
                        "--max-reference-flux-std-dev-rel",
                        "0.03",
                        "--require-mg-flux-std-dev",
                        "--max-mg-flux-std-dev-rel",
                        "0.05",
                        "--summary-json",
                        str(summary),
                    ]
                ),
                0,
            )

            payload = json.loads(summary.read_text(encoding="utf-8"))
            self.assertEqual(
                payload["reference_flux_std_dev_dataset"],
                "openmc_volume_flux_std_dev",
            )
            self.assertAlmostEqual(payload["reference_flux_max_relative_std_dev"], 0.02)
            self.assertEqual(payload["mg_flux_std_dev_dataset"], "openmc_mg_flux_std_dev")
            self.assertAlmostEqual(payload["mg_flux_max_relative_std_dev"], 0.04)
            with h5py.File(sidecar, "r") as h5:
                for prefix, expected_limit, expected_observed in (
                    ("sph_reference_flux", 0.03, 0.02),
                    ("sph_mg_flux", 0.05, 0.04),
                ):
                    self.assertTrue(
                        bool(h5.attrs[f"{prefix}_uncertainty_require_coverage"])
                    )
                    self.assertTrue(bool(h5.attrs[f"{prefix}_uncertainty_coverage"]))
                    self.assertAlmostEqual(
                        float(h5.attrs[f"{prefix}_uncertainty_limit"]),
                        expected_limit,
                    )
                    self.assertAlmostEqual(
                        float(h5.attrs[f"{prefix}_uncertainty_observed_max_rel"]),
                        expected_observed,
                    )
                    self.assertTrue(bool(h5.attrs[f"{prefix}_uncertainty_pass"]))
                    self.assertTrue(bool(h5.attrs[f"{prefix}_layout_verified"]))

    def test_openmc_sph_sidecar_binds_previous_sph_file(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "mgxs.h5"
            reference_flux = root / "reference_flux.h5"
            mg_flux = root / "mg_flux.h5"
            previous_sph = root / "previous_sph.csv"
            sidecar = root / "sph.h5"
            write_mgxs(mgxs)
            _write_flux_source(
                reference_flux,
                "reference_flux",
                values=np.ones((2, 2)),
            )
            _write_flux_source(
                mg_flux,
                "mg_flux",
                values=np.ones((2, 2)),
            )
            previous_sph.write_text(
                "mixture,g1,g2\nfuel,1.0,1.0\nmoderator,1.0,1.0\n",
                encoding="utf-8",
            )

            create_openmc_sph_sidecar(
                mgxs,
                sidecar,
                reference_flux=f"{reference_flux}::reference_flux",
                mg_flux=f"{mg_flux}::mg_flux",
                previous_sph=previous_sph,
                flux_normalization="none",
            )

            with h5py.File(sidecar, "r") as h5:
                self.assertTrue(bool(h5.attrs["sph_previous_sph_used"]))
                self.assertEqual(
                    h5.attrs["sph_previous_sph_sha256"],
                    file_sha256(previous_sph),
                )
                self.assertEqual(
                    Path(str(h5.attrs["sph_previous_sph_path"])),
                    previous_sph.resolve(),
                )

    def test_openmc_sph_sidecar_can_require_flux_uncertainty(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "mgxs.h5"
            reference_flux = root / "openmc_ce_flux.h5"
            mg_flux = root / "openmc_mg_flux.h5"
            sidecar = root / "openmc_sph.h5"
            write_mgxs(mgxs)
            _write_flux_source(reference_flux, "openmc_volume_flux", values=np.ones((2, 2)))
            _write_flux_source(
                mg_flux,
                "openmc_mg_flux",
                values=np.ones((2, 2)),
                std_dev=np.full((2, 2), 0.01),
            )

            with self.assertRaises(SystemExit) as ctx:
                cli_main(
                    [
                        "make-openmc-sph-sidecar",
                        str(mgxs),
                        "-o",
                        str(sidecar),
                        "--reference-flux",
                        f"{reference_flux}::openmc_volume_flux",
                        "--mg-flux",
                        f"{mg_flux}::openmc_mg_flux",
                        "--flux-normalization",
                        "none",
                        "--require-reference-flux-std-dev",
                    ]
                )
            self.assertEqual(ctx.exception.code, 1)

    def test_openmc_sph_sidecar_rejects_noisy_flux_uncertainty(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "mgxs.h5"
            reference_flux = root / "openmc_ce_flux.h5"
            mg_flux = root / "openmc_mg_flux.h5"
            sidecar = root / "openmc_sph.h5"
            write_mgxs(mgxs)
            _write_flux_source(
                reference_flux,
                "openmc_volume_flux",
                values=np.ones((2, 2)),
                std_dev=np.full((2, 2), 0.20),
            )
            _write_flux_source(
                mg_flux,
                "openmc_mg_flux",
                values=np.ones((2, 2)),
                std_dev=np.full((2, 2), 0.01),
            )

            with self.assertRaises(SystemExit) as ctx:
                cli_main(
                    [
                        "make-openmc-sph-sidecar",
                        str(mgxs),
                        "-o",
                        str(sidecar),
                        "--reference-flux",
                        f"{reference_flux}::openmc_volume_flux",
                        "--mg-flux",
                        f"{mg_flux}::openmc_mg_flux",
                        "--flux-normalization",
                        "none",
                        "--max-reference-flux-std-dev-rel",
                        "0.05",
                    ]
                )
            self.assertEqual(ctx.exception.code, 1)

    def test_cli_builds_damped_sph_update_table(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "mgxs.h5"
            reference_flux = root / "reference_flux.csv"
            low_order_flux = root / "low_order_flux.h5"
            previous_sph = root / "previous_sph.csv"
            table = root / "next_sph.csv"
            sidecar = root / "next_sph.h5"
            summary = root / "summary.json"
            write_mgxs(
                mgxs,
                h_factor={"fuel": np.asarray([19.0, 21.0])},
            )
            reference_flux.write_text(
                "\n".join(
                    [
                        "mixture,group,reference_flux",
                        "moderator,2,1.44",
                        "fuel,1,1.21",
                        "moderator,1,0.64",
                        "fuel,2,0.81",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            previous_sph.write_text(
                "mixture,g1,g2\nfuel,1.0,1.1\nmoderator,0.9,1.0\n",
                encoding="utf-8",
            )
            with h5py.File(low_order_flux, "w") as h5:
                data = np.array([[1.0, 1.0], [1.0, 1.0]])
                dataset = h5.create_dataset("volume_flux", data=data)
                dataset.attrs["mixture_names"] = np.asarray(("fuel", "moderator"), dtype="S")
                dataset.attrs["group_order"] = "mgxs_donjon"

            self.assertEqual(
                cli_main(
                    [
                        "make-sph-update-table",
                        str(mgxs),
                        "-o",
                        str(table),
                        "--reference-flux",
                        str(reference_flux),
                        "--low-order-flux",
                        f"{low_order_flux}::volume_flux",
                        "--previous-sph",
                        str(previous_sph),
                        "--damping",
                        "0.5",
                        "--summary-json",
                        str(summary),
                    ]
                ),
                0,
            )
            self.assertEqual(
                cli_main(
                    [
                        "make-sph-sidecar",
                        str(mgxs),
                        "-o",
                        str(sidecar),
                        "--mode",
                        "table",
                        "--table",
                        str(table),
                    ]
                ),
                0,
            )

            normalization_factor = 40.0 / (19.0 + 21.0 / 1.1)
            expected = np.array(
                [
                    [
                        1.0 * np.sqrt(normalization_factor / (1.0 * 1.21)),
                        1.1 * np.sqrt(normalization_factor / (1.1 * 0.81)),
                    ],
                    [
                        0.9 * np.sqrt(normalization_factor / (0.9 * 0.64)),
                        1.0 * np.sqrt(normalization_factor / (1.0 * 1.44)),
                    ],
                ]
            )
            with h5py.File(sidecar, "r") as h5:
                np.testing.assert_allclose(h5["sph"][:], expected)
            payload = json.loads(summary.read_text(encoding="utf-8"))
            self.assertEqual(payload["decision"], "openmc2donjon_sph_iteration_table_passed")
            self.assertEqual(
                payload["formula"],
                "next_sph = previous_sph * "
                "(normalized_low_order_flux / (previous_sph * reference_flux)) ** damping",
            )
            self.assertEqual(payload["sph_target"], "rate")
            self.assertEqual(payload["flux_normalization"], "power")
            self.assertAlmostEqual(
                payload["normalization_factor"],
                normalization_factor,
            )
            self.assertEqual(payload["energy_groups"], 2)
            self.assertEqual(payload["clipped_count"], 0)
            self.assertEqual(payload["clipped_bins"], [])
            self.assertEqual(payload["diagnostic_bin_limit"], 10)
            worst = payload["worst_residual_bins"][0]
            self.assertEqual(worst["mixture"], "moderator")
            self.assertEqual(worst["group"], 1)
            self.assertAlmostEqual(
                worst["raw_update"],
                normalization_factor / (0.9 * 0.64),
            )
            self.assertAlmostEqual(
                worst["residual"],
                normalization_factor / (0.9 * 0.64) - 1.0,
            )

    def test_rate_power_normalization_preserves_nonuniform_sph_fixed_point(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "mgxs.h5"
            reference_flux = root / "reference_flux.csv"
            low_order_flux = root / "low_order_flux.csv"
            previous_sph = root / "previous_sph.csv"
            table = root / "next_sph.csv"
            write_mgxs(
                mgxs,
                h_factor={
                    "fuel": np.asarray([2.0, 3.0]),
                    "moderator": np.asarray([5.0, 7.0]),
                },
            )
            previous = np.asarray([[0.5, 2.0], [1.5, 0.8]])
            reference_flux.write_text(
                "mixture,group,flux\n"
                "fuel,1,1.0\nfuel,2,2.0\n"
                "moderator,1,3.0\nmoderator,2,4.0\n",
                encoding="utf-8",
            )
            low_order_flux.write_text(
                "mixture,group,flux\n"
                "fuel,1,0.5\nfuel,2,4.0\n"
                "moderator,1,4.5\nmoderator,2,3.2\n",
                encoding="utf-8",
            )
            previous_sph.write_text(
                "mixture,g1,g2\nfuel,0.5,2.0\nmoderator,1.5,0.8\n",
                encoding="utf-8",
            )

            report = create_sph_update_table(
                mgxs,
                table,
                reference_flux=reference_flux,
                low_order_flux=low_order_flux,
                previous_sph=previous_sph,
                sph_target="rate",
                flux_normalization="power",
            )

            np.testing.assert_allclose(_read_sph_table(table), previous, rtol=1.0e-12)
            self.assertAlmostEqual(report.normalization_factor, 1.0)
            self.assertAlmostEqual(
                report.reference_normalization_integral,
                report.low_order_normalization_integral,
            )
            self.assertAlmostEqual(report.raw_update_minimum, 1.0)
            self.assertAlmostEqual(report.raw_update_maximum, 1.0)

    def test_power_normalization_scales_low_order_flux_with_h_factor(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "mgxs.h5"
            reference_flux = root / "reference_flux.csv"
            low_order_flux = root / "low_order_flux.csv"
            table = root / "next_sph.csv"
            summary = root / "summary.json"
            write_mgxs(
                mgxs,
                h_factor={
                    "fuel": np.asarray([10.0, 100.0]),
                    "moderator": np.asarray([1.0, 1.0]),
                },
            )
            reference_flux.write_text(
                "mixture,group,flux\n"
                "fuel,1,10.0\nfuel,2,20.0\n"
                "moderator,1,30.0\nmoderator,2,40.0\n",
                encoding="utf-8",
            )
            low_order_flux.write_text(
                "mixture,group,flux\n"
                "fuel,1,1.0\nfuel,2,1.0\n"
                "moderator,1,1.0\nmoderator,2,1.0\n",
                encoding="utf-8",
            )

            report = create_sph_update_table(
                mgxs,
                table,
                reference_flux=reference_flux,
                low_order_flux=low_order_flux,
                flux_normalization="power",
                sph_target="flux",
                summary_json=summary,
            )

            factor = (10.0 * 10.0 + 20.0 * 100.0 + 30.0 + 40.0) / (10.0 + 100.0 + 1.0 + 1.0)
            expected = np.asarray([[10.0 / factor, 20.0 / factor], [30.0 / factor, 40.0 / factor]])
            self.assertAlmostEqual(report.normalization_factor, factor)
            self.assertEqual(report.flux_normalization, "power")
            rows = table.read_text(encoding="utf-8").strip().splitlines()[1:]
            actual = np.asarray([float(row.split(",")[2]) for row in rows]).reshape(2, 2)
            np.testing.assert_allclose(actual, expected, rtol=1.0e-11)
            payload = json.loads(summary.read_text(encoding="utf-8"))
            self.assertEqual(payload["flux_normalization"], "power")
            self.assertAlmostEqual(payload["normalization_factor"], factor)
            self.assertEqual(payload["normalization_weight_source"], "H-FACTOR/kappa_fission")
            self.assertAlmostEqual(payload["reference_normalization_integral"], 2170.0)
            self.assertAlmostEqual(payload["low_order_normalization_integral"], 112.0)

    def test_power_normalization_requires_h_factor(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "mgxs.h5"
            reference_flux = root / "reference_flux.csv"
            low_order_flux = root / "low_order_flux.csv"
            table = root / "next_sph.csv"
            write_mgxs(mgxs)
            reference_flux.write_text(
                "mixture,group,flux\nfuel,1,1.0\nfuel,2,1.0\nmoderator,1,1.0\nmoderator,2,1.0\n",
                encoding="utf-8",
            )
            low_order_flux.write_text(
                "mixture,group,flux\nfuel,1,1.0\nfuel,2,1.0\nmoderator,1,1.0\nmoderator,2,1.0\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "requires group-wise H-FACTOR"):
                create_sph_update_table(
                    mgxs,
                    table,
                    reference_flux=reference_flux,
                    low_order_flux=low_order_flux,
                    flux_normalization="power",
                )

    def test_power_normalization_allows_missing_nonfissionable_h_factor(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "mgxs.h5"
            reference_flux = root / "reference_flux.csv"
            low_order_flux = root / "low_order_flux.csv"
            table = root / "next_sph.csv"
            write_mgxs(
                mgxs,
                h_factor={
                    "fuel": np.asarray([10.0, 100.0]),
                },
            )
            reference_flux.write_text(
                "mixture,group,flux\n"
                "fuel,1,10.0\nfuel,2,20.0\n"
                "moderator,1,30.0\nmoderator,2,40.0\n",
                encoding="utf-8",
            )
            low_order_flux.write_text(
                "mixture,group,flux\n"
                "fuel,1,1.0\nfuel,2,1.0\n"
                "moderator,1,1.0\nmoderator,2,1.0\n",
                encoding="utf-8",
            )

            report = create_sph_update_table(
                mgxs,
                table,
                reference_flux=reference_flux,
                low_order_flux=low_order_flux,
                flux_normalization="power",
                sph_target="flux",
            )

            factor = (10.0 * 10.0 + 20.0 * 100.0) / (10.0 + 100.0)
            expected = np.asarray(
                [
                    [10.0 / factor, 20.0 / factor],
                    [30.0 / factor, 40.0 / factor],
                ]
            )
            self.assertAlmostEqual(report.normalization_factor, factor)
            self.assertAlmostEqual(report.reference_normalization_integral, 2100.0)
            self.assertAlmostEqual(report.low_order_normalization_integral, 110.0)
            rows = table.read_text(encoding="utf-8").strip().splitlines()[1:]
            actual = np.asarray([float(row.split(",")[2]) for row in rows]).reshape(2, 2)
            np.testing.assert_allclose(actual, expected, rtol=1.0e-11)

    def test_default_normalization_resolves_to_power_when_h_factor_exists(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "mgxs.h5"
            reference_flux = root / "reference_flux.csv"
            low_order_flux = root / "low_order_flux.csv"
            table = root / "next_sph.csv"
            summary = root / "summary.json"
            write_mgxs(
                mgxs,
                h_factor={
                    "fuel": np.asarray([10.0, 100.0]),
                    "moderator": np.asarray([1.0, 1.0]),
                },
            )
            reference_flux.write_text(
                "mixture,group,flux\n"
                "fuel,1,10.0\nfuel,2,20.0\n"
                "moderator,1,30.0\nmoderator,2,40.0\n",
                encoding="utf-8",
            )
            low_order_flux.write_text(
                "mixture,group,flux\n"
                "fuel,1,1.0\nfuel,2,1.0\n"
                "moderator,1,1.0\nmoderator,2,1.0\n",
                encoding="utf-8",
            )

            report = create_sph_update_table(
                mgxs,
                table,
                reference_flux=reference_flux,
                low_order_flux=low_order_flux,
                summary_json=summary,
            )

            self.assertEqual(report.flux_normalization, "power")
            self.assertEqual(report.sph_target, "rate")
            self.assertEqual(
                report.normalization_weight_source,
                "H-FACTOR/kappa_fission (auto)",
            )
            payload = json.loads(summary.read_text(encoding="utf-8"))
            self.assertEqual(payload["flux_normalization"], "power")
            self.assertEqual(
                payload["normalization_weight_source"],
                "H-FACTOR/kappa_fission (auto)",
            )

    def test_default_normalization_requires_h_factor(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "mgxs.h5"
            reference_flux = root / "reference_flux.csv"
            low_order_flux = root / "low_order_flux.csv"
            table = root / "next_sph.csv"
            write_mgxs(mgxs)
            reference_flux.write_text(
                "mixture,group,flux\nfuel,1,1.0\nfuel,2,1.0\nmoderator,1,1.0\nmoderator,2,1.0\n",
                encoding="utf-8",
            )
            low_order_flux.write_text(
                "mixture,group,flux\nfuel,1,1.0\nfuel,2,1.0\nmoderator,1,1.0\nmoderator,2,1.0\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "auto flux normalization requires"):
                create_sph_update_table(
                    mgxs,
                    table,
                    reference_flux=reference_flux,
                    low_order_flux=low_order_flux,
                )

    def test_rejects_nonpositive_low_order_flux(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "mgxs.h5"
            reference_flux = root / "reference_flux.csv"
            low_order_flux = root / "low_order_flux.csv"
            table = root / "next_sph.csv"
            write_mgxs(mgxs)
            reference_flux.write_text(
                "mixture,group,flux\nfuel,1,1.0\nfuel,2,1.0\nmoderator,1,1.0\nmoderator,2,1.0\n",
                encoding="utf-8",
            )
            low_order_flux.write_text(
                "mixture,group,flux\nfuel,1,1.0\nfuel,2,0.0\nmoderator,1,1.0\nmoderator,2,1.0\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "low-order flux values must be positive"):
                create_sph_update_table(
                    mgxs,
                    table,
                    reference_flux=reference_flux,
                    low_order_flux=low_order_flux,
                )

    def test_rejects_hdf5_flux_with_wrong_group_order(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "mgxs.h5"
            reference_flux = root / "reference_flux.csv"
            low_order_flux = root / "low_order_flux.h5"
            table = root / "next_sph.csv"
            write_mgxs(mgxs)
            reference_flux.write_text(
                "mixture,group,flux\n"
                "fuel,1,1.0\nfuel,2,1.0\n"
                "moderator,1,1.0\nmoderator,2,1.0\n",
                encoding="utf-8",
            )
            with h5py.File(low_order_flux, "w") as h5:
                dataset = h5.create_dataset(
                    "low_order_flux",
                    data=np.asarray([[1.0, 1.0], [1.0, 1.0]]),
                )
                dataset.attrs["mixture_names"] = np.asarray(
                    ("fuel", "moderator"),
                    dtype="S",
                )
                dataset.attrs["group_order"] = "ascending_energy"

            with self.assertRaisesRegex(
                ValueError,
                "low-order flux: group_order must be 'mgxs_donjon'",
            ):
                create_sph_update_table(
                    mgxs,
                    table,
                    reference_flux=reference_flux,
                    low_order_flux=f"{low_order_flux}::low_order_flux",
                )

    def test_cli_accepts_mesh_shaped_hdf5_flux(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "mgxs.h5"
            flux = root / "mesh_flux.h5"
            table = root / "next_sph.csv"
            sidecar = root / "next_sph.h5"
            write_mgxs(mgxs)
            with h5py.File(flux, "w") as h5:
                names = np.asarray([["moderator", "fuel"]], dtype="S")
                h5.create_dataset("mixture_names", data=names)
                reference = h5.create_dataset(
                    "reference_flux",
                    data=np.asarray([[[4.0, 9.0], [16.0, 25.0]]]),
                )
                reference.attrs["group_order"] = "mgxs_donjon"
                low_order = h5.create_dataset(
                    "low_order_flux",
                    data=np.asarray([[[1.0, 1.0], [1.0, 1.0]]]),
                )
                low_order.attrs["group_order"] = "mgxs_donjon"

            self.assertEqual(
                cli_main(
                    [
                        "make-sph-update-table",
                        str(mgxs),
                        "-o",
                        str(table),
                        "--reference-flux",
                        f"{flux}::reference_flux",
                        "--low-order-flux",
                        f"{flux}::low_order_flux",
                        "--sph-target",
                        "flux",
                        "--flux-normalization",
                        "none",
                        "--damping",
                        "0.5",
                    ]
                ),
                0,
            )
            self.assertEqual(
                cli_main(
                    [
                        "make-sph-sidecar",
                        str(mgxs),
                        "-o",
                        str(sidecar),
                        "--mode",
                        "table",
                        "--table",
                        str(table),
                    ]
                ),
                0,
            )

            expected = np.asarray([[4.0, 5.0], [2.0, 3.0]])
            with h5py.File(sidecar, "r") as h5:
                np.testing.assert_allclose(h5["sph"][:], expected)

    def test_summary_records_clipped_bins(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "mgxs.h5"
            reference_flux = root / "reference_flux.csv"
            low_order_flux = root / "low_order_flux.csv"
            table = root / "next_sph.csv"
            summary = root / "summary.json"
            write_mgxs(mgxs)
            reference_flux.write_text(
                "mixture,group,flux\n"
                "fuel,1,3.0\nfuel,2,1.0\n"
                "moderator,1,1.0\nmoderator,2,2.0\n",
                encoding="utf-8",
            )
            low_order_flux.write_text(
                "mixture,group,flux\n"
                "fuel,1,1.0\nfuel,2,1.0\n"
                "moderator,1,1.0\nmoderator,2,1.0\n",
                encoding="utf-8",
            )

            report = create_sph_update_table(
                mgxs,
                table,
                reference_flux=reference_flux,
                low_order_flux=low_order_flux,
                sph_target="flux",
                flux_normalization="none",
                clip_max=1.5,
                summary_json=summary,
            )

            self.assertEqual(report.clipped_count, 2)
            payload = json.loads(summary.read_text(encoding="utf-8"))
            self.assertEqual(payload["clipped_count"], 2)
            clipped = payload["clipped_bins"][0]
            self.assertEqual(clipped["mixture"], "fuel")
            self.assertEqual(clipped["group"], 1)
            self.assertAlmostEqual(clipped["unclipped_sph"], 3.0)
            self.assertAlmostEqual(clipped["sph"], 1.5)
            self.assertTrue(clipped["clipped"])

    def test_cli_accepts_previous_sph_sidecar(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "mgxs.h5"
            previous = root / "previous_sph.h5"
            reference_flux = root / "reference_flux.csv"
            low_order_flux = root / "low_order_flux.csv"
            table = root / "next_sph.csv"
            write_mgxs(mgxs)
            with h5py.File(previous, "w") as h5:
                dataset = h5.create_dataset(
                    "sph",
                    data=np.asarray([[1.1, 1.2], [0.9, 1.0]]),
                )
                dataset.attrs["mixture_names"] = np.asarray(("fuel", "moderator"), dtype="S")
            reference_flux.write_text(
                "mixture,group,flux\nfuel,1,4.0\nfuel,2,9.0\nmoderator,1,16.0\nmoderator,2,25.0\n",
                encoding="utf-8",
            )
            low_order_flux.write_text(
                "mixture,group,flux\nfuel,1,1.0\nfuel,2,1.0\nmoderator,1,1.0\nmoderator,2,1.0\n",
                encoding="utf-8",
            )

            self.assertEqual(
                cli_main(
                    [
                        "make-sph-update-table",
                        str(mgxs),
                        "-o",
                        str(table),
                        "--reference-flux",
                        str(reference_flux),
                        "--low-order-flux",
                        str(low_order_flux),
                        "--previous-sph",
                        str(previous),
                        "--sph-target",
                        "flux",
                        "--flux-normalization",
                        "none",
                        "--damping",
                        "0.5",
                    ]
                ),
                0,
            )

            rows = table.read_text(encoding="utf-8").strip().splitlines()
            self.assertEqual(rows[0], "mixture,group,sph")
            self.assertIn("fuel,1,2.2", rows)
            self.assertIn("fuel,2,3.6", rows)
            self.assertIn("moderator,1,3.6", rows)
            self.assertIn("moderator,2,5", rows)

    def test_rejects_hdf5_flux_without_required_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "mgxs.h5"
            reference_flux = root / "reference_flux.csv"
            low_order_flux = root / "low_order_flux.h5"
            table = root / "next_sph.csv"
            write_mgxs(mgxs)
            reference_flux.write_text(
                "mixture,group,flux\n"
                "fuel,1,1.0\nfuel,2,1.0\n"
                "moderator,1,1.0\nmoderator,2,1.0\n",
                encoding="utf-8",
            )
            with h5py.File(low_order_flux, "w") as h5:
                h5.create_dataset(
                    "low_order_flux",
                    data=np.asarray([[1.0, 1.0], [1.0, 1.0]]),
                )

            with self.assertRaisesRegex(
                ValueError,
                "low-order flux: group_order must be 'mgxs_donjon'",
            ):
                create_sph_update_table(
                    mgxs,
                    table,
                    reference_flux=reference_flux,
                    low_order_flux=f"{low_order_flux}::low_order_flux",
                )

            with h5py.File(low_order_flux, "a") as h5:
                h5["low_order_flux"].attrs["group_order"] = "mgxs_donjon"

            with self.assertRaisesRegex(
                ValueError,
                "low-order flux: HDF5 flux sources must declare mixture_names",
            ):
                create_sph_update_table(
                    mgxs,
                    table,
                    reference_flux=reference_flux,
                    low_order_flux=f"{low_order_flux}::low_order_flux",
                )

    def test_identity_zero_flux_policy_passes_through_matched_zero_bins(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "mgxs.h5"
            reference_flux = root / "reference_flux.csv"
            low_order_flux = root / "low_order_flux.csv"
            control_reference_flux = root / "control_reference_flux.csv"
            control_low_order_flux = root / "control_low_order_flux.csv"
            table = root / "next_sph.csv"
            control_table = root / "control_sph.csv"
            summary = root / "summary.json"
            write_mgxs(mgxs)
            reference_flux.write_text(
                "mixture,group,flux\n"
                "fuel,1,4.0\nfuel,2,0.0\n"
                "moderator,1,16.0\nmoderator,2,25.0\n",
                encoding="utf-8",
            )
            low_order_flux.write_text(
                "mixture,group,flux\n"
                "fuel,1,1.0\nfuel,2,0.0\n"
                "moderator,1,1.0\nmoderator,2,1.0\n",
                encoding="utf-8",
            )
            control_reference_flux.write_text(
                "mixture,group,flux\n"
                "fuel,1,4.0\nfuel,2,9.0\n"
                "moderator,1,16.0\nmoderator,2,25.0\n",
                encoding="utf-8",
            )
            control_low_order_flux.write_text(
                "mixture,group,flux\n"
                "fuel,1,1.0\nfuel,2,1.0\n"
                "moderator,1,1.0\nmoderator,2,1.0\n",
                encoding="utf-8",
            )

            report = create_sph_update_table(
                mgxs,
                table,
                reference_flux=reference_flux,
                low_order_flux=low_order_flux,
                sph_target="flux",
                flux_normalization="none",
                zero_flux_policy="identity",
                summary_json=summary,
            )
            control_report = create_sph_update_table(
                mgxs,
                control_table,
                reference_flux=control_reference_flux,
                low_order_flux=control_low_order_flux,
                sph_target="flux",
                flux_normalization="none",
            )

            self.assertEqual(report.zero_flux_policy, "identity")
            self.assertEqual(report.identity_bin_count, 1)
            self.assertEqual(control_report.zero_flux_policy, "reject")
            self.assertEqual(control_report.identity_bin_count, 0)
            actual = _read_sph_table(table)
            control = _read_sph_table(control_table)
            # The identity bin keeps the previous SPH value (unity here).
            self.assertEqual(actual[0, 1], 1.0)
            np.testing.assert_allclose(actual[0, 0], control[0, 0], rtol=1.0e-11)
            np.testing.assert_allclose(actual[1, :], control[1, :], rtol=1.0e-11)
            payload = json.loads(summary.read_text(encoding="utf-8"))
            self.assertEqual(payload["zero_flux_policy"], "identity")
            self.assertEqual(payload["identity_bin_count"], 1)

    def test_identity_zero_flux_policy_keeps_previous_sph_value(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "mgxs.h5"
            reference_flux = root / "reference_flux.csv"
            low_order_flux = root / "low_order_flux.csv"
            previous_sph = root / "previous_sph.csv"
            table = root / "next_sph.csv"
            write_mgxs(mgxs)
            reference_flux.write_text(
                "mixture,group,flux\n"
                "fuel,1,4.0\nfuel,2,0.0\n"
                "moderator,1,16.0\nmoderator,2,25.0\n",
                encoding="utf-8",
            )
            low_order_flux.write_text(
                "mixture,group,flux\n"
                "fuel,1,1.0\nfuel,2,0.0\n"
                "moderator,1,1.0\nmoderator,2,1.0\n",
                encoding="utf-8",
            )
            previous_sph.write_text(
                "mixture,g1,g2\nfuel,1.0,1.1\nmoderator,0.9,1.0\n",
                encoding="utf-8",
            )

            create_sph_update_table(
                mgxs,
                table,
                reference_flux=reference_flux,
                low_order_flux=low_order_flux,
                previous_sph=previous_sph,
                sph_target="flux",
                flux_normalization="none",
                zero_flux_policy="identity",
            )

            actual = _read_sph_table(table)
            self.assertEqual(actual[0, 1], 1.1)
            np.testing.assert_allclose(actual[0, 0], 1.0 * 4.0, rtol=1.0e-11)

    def test_identity_zero_flux_policy_rejects_one_sided_zero(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "mgxs.h5"
            reference_flux = root / "reference_flux.csv"
            low_order_flux = root / "low_order_flux.csv"
            table = root / "next_sph.csv"
            write_mgxs(mgxs)
            reference_flux.write_text(
                "mixture,group,flux\n"
                "fuel,1,4.0\nfuel,2,0.0\n"
                "moderator,1,16.0\nmoderator,2,25.0\n",
                encoding="utf-8",
            )
            low_order_flux.write_text(
                "mixture,group,flux\n"
                "fuel,1,1.0\nfuel,2,1.0\n"
                "moderator,1,1.0\nmoderator,2,1.0\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(
                ValueError,
                r"CE/MG inconsistency: fuel g2 \(zero reference flux\)",
            ):
                create_sph_update_table(
                    mgxs,
                    table,
                    reference_flux=reference_flux,
                    low_order_flux=low_order_flux,
                    zero_flux_policy="identity",
                )

    def test_default_policy_rejects_matched_zero_bins(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "mgxs.h5"
            reference_flux = root / "reference_flux.csv"
            low_order_flux = root / "low_order_flux.csv"
            table = root / "next_sph.csv"
            write_mgxs(mgxs)
            reference_flux.write_text(
                "mixture,group,flux\n"
                "fuel,1,4.0\nfuel,2,0.0\n"
                "moderator,1,16.0\nmoderator,2,25.0\n",
                encoding="utf-8",
            )
            low_order_flux.write_text(
                "mixture,group,flux\n"
                "fuel,1,1.0\nfuel,2,0.0\n"
                "moderator,1,1.0\nmoderator,2,1.0\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "reference flux values must be positive"):
                create_sph_update_table(
                    mgxs,
                    table,
                    reference_flux=reference_flux,
                    low_order_flux=low_order_flux,
                )

            with self.assertRaisesRegex(ValueError, "--zero-flux-policy must be one of"):
                create_sph_update_table(
                    mgxs,
                    table,
                    reference_flux=reference_flux,
                    low_order_flux=low_order_flux,
                    zero_flux_policy="bogus",
                )

    def test_identity_zero_flux_policy_excludes_identity_bins_from_std_dev_gate(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "mgxs.h5"
            reference_flux = root / "openmc_ce_flux.h5"
            mg_flux = root / "openmc_mg_flux.h5"
            table = root / "next_sph.csv"
            write_mgxs(mgxs)
            _write_flux_source(
                reference_flux,
                "openmc_volume_flux",
                values=np.array([[2.0, 0.0], [4.0, 2.0]]),
                std_dev=np.array([[0.02, 0.0], [0.08, 0.02]]),
            )
            _write_flux_source(
                mg_flux,
                "openmc_mg_flux",
                values=np.array([[1.0, 0.0], [1.0, 1.0]]),
                std_dev=np.array([[0.01, 0.0], [0.03, 0.04]]),
            )

            report = create_sph_update_table(
                mgxs,
                table,
                reference_flux=f"{reference_flux}::openmc_volume_flux",
                low_order_flux=f"{mg_flux}::openmc_mg_flux",
                sph_target="flux",
                flux_normalization="none",
                zero_flux_policy="identity",
                require_reference_flux_std_dev=True,
                max_reference_flux_std_dev_rel=0.03,
                require_low_order_flux_std_dev=True,
                max_low_order_flux_std_dev_rel=0.05,
            )

            self.assertEqual(report.identity_bin_count, 1)
            self.assertAlmostEqual(report.reference_flux_max_relative_std_dev or 0.0, 0.02)
            self.assertAlmostEqual(report.low_order_flux_max_relative_std_dev or 0.0, 0.04)

    def test_cli_openmc_sph_sidecar_accepts_zero_flux_policy(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "mgxs.h5"
            reference_flux = root / "openmc_ce_flux.h5"
            mg_flux = root / "openmc_mg_flux.h5"
            sidecar = root / "openmc_sph.h5"
            table = root / "openmc_sph.csv"
            summary = root / "openmc_sph_summary.json"
            write_mgxs(mgxs)
            _write_flux_source(
                reference_flux,
                "openmc_volume_flux",
                values=np.array([[4.0, 0.0], [16.0, 25.0]]),
            )
            _write_flux_source(
                mg_flux,
                "openmc_mg_flux",
                values=np.array([[1.0, 0.0], [1.0, 1.0]]),
            )

            self.assertEqual(
                cli_main(
                    [
                        "make-openmc-sph-sidecar",
                        str(mgxs),
                        "-o",
                        str(sidecar),
                        "--reference-flux",
                        f"{reference_flux}::openmc_volume_flux",
                        "--mg-flux",
                        f"{mg_flux}::openmc_mg_flux",
                        "--table-output",
                        str(table),
                        "--sph-target",
                        "flux",
                        "--flux-normalization",
                        "none",
                        "--zero-flux-policy",
                        "identity",
                        "--summary-json",
                        str(summary),
                    ]
                ),
                0,
            )

            expected = np.array([[4.0, 1.0], [16.0, 25.0]])
            with h5py.File(sidecar, "r") as h5:
                np.testing.assert_allclose(h5["sph"][:], expected, rtol=1.0e-11)
            payload = json.loads(summary.read_text(encoding="utf-8"))
            self.assertEqual(payload["zero_flux_policy"], "identity")
            self.assertEqual(payload["identity_bin_count"], 1)

    def test_flux_floor_freezes_low_flux_bins_and_keeps_previous_sph(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "mgxs.h5"
            reference_flux = root / "reference_flux.csv"
            low_order_flux = root / "low_order_flux.csv"
            control_reference_flux = root / "control_reference_flux.csv"
            control_low_order_flux = root / "control_low_order_flux.csv"
            previous_sph = root / "previous_sph.csv"
            table = root / "next_sph.csv"
            control_table = root / "control_sph.csv"
            summary = root / "summary.json"
            write_mgxs(mgxs)
            reference_flux.write_text(
                "mixture,group,flux\n"
                "fuel,1,4.0\nfuel,2,1.0e-5\n"
                "moderator,1,16.0\nmoderator,2,25.0\n",
                encoding="utf-8",
            )
            low_order_flux.write_text(
                "mixture,group,flux\n"
                "fuel,1,1.0\nfuel,2,0.5\n"
                "moderator,1,1.0\nmoderator,2,1.0\n",
                encoding="utf-8",
            )
            control_reference_flux.write_text(
                "mixture,group,flux\n"
                "fuel,1,4.0\nfuel,2,9.0\n"
                "moderator,1,16.0\nmoderator,2,25.0\n",
                encoding="utf-8",
            )
            control_low_order_flux.write_text(
                "mixture,group,flux\n"
                "fuel,1,1.0\nfuel,2,1.0\n"
                "moderator,1,1.0\nmoderator,2,1.0\n",
                encoding="utf-8",
            )
            previous_sph.write_text(
                "mixture,g1,g2\nfuel,1.0,1.1\nmoderator,0.9,1.0\n",
                encoding="utf-8",
            )

            report = create_sph_update_table(
                mgxs,
                table,
                reference_flux=reference_flux,
                low_order_flux=low_order_flux,
                previous_sph=previous_sph,
                sph_target="flux",
                flux_normalization="none",
                flux_floor_rel=1.0e-3,
                summary_json=summary,
            )
            control_report = create_sph_update_table(
                mgxs,
                control_table,
                reference_flux=control_reference_flux,
                low_order_flux=control_low_order_flux,
                previous_sph=previous_sph,
                sph_target="flux",
                flux_normalization="none",
            )

            self.assertEqual(report.flux_floor_rel, 1.0e-3)
            self.assertEqual(report.floored_bin_count, 1)
            self.assertEqual(report.identity_bin_count, 0)
            self.assertIsNone(control_report.flux_floor_rel)
            self.assertEqual(control_report.floored_bin_count, 0)
            actual = _read_sph_table(table)
            control = _read_sph_table(control_table)
            # The floored bin keeps the previous SPH value.
            self.assertEqual(actual[0, 1], 1.1)
            np.testing.assert_allclose(actual[0, 0], control[0, 0], rtol=1.0e-11)
            np.testing.assert_allclose(actual[1, :], control[1, :], rtol=1.0e-11)
            payload = json.loads(summary.read_text(encoding="utf-8"))
            self.assertEqual(payload["flux_floor_rel"], 1.0e-3)
            self.assertEqual(payload["floored_bin_count"], 1)

    def test_flux_floor_exempts_one_sided_zero_below_floor(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "mgxs.h5"
            reference_flux = root / "reference_flux.csv"
            low_order_flux = root / "low_order_flux.csv"
            table = root / "next_sph.csv"
            identity_table = root / "identity_sph.csv"
            write_mgxs(mgxs)
            reference_flux.write_text(
                "mixture,group,flux\n"
                "fuel,1,4.0\nfuel,2,1.0e-6\n"
                "moderator,1,16.0\nmoderator,2,25.0\n",
                encoding="utf-8",
            )
            low_order_flux.write_text(
                "mixture,group,flux\n"
                "fuel,1,1.0\nfuel,2,0.0\n"
                "moderator,1,1.0\nmoderator,2,1.0\n",
                encoding="utf-8",
            )

            report = create_sph_update_table(
                mgxs,
                table,
                reference_flux=reference_flux,
                low_order_flux=low_order_flux,
                sph_target="flux",
                flux_normalization="none",
                flux_floor_rel=1.0e-3,
            )
            identity_report = create_sph_update_table(
                mgxs,
                identity_table,
                reference_flux=reference_flux,
                low_order_flux=low_order_flux,
                sph_target="flux",
                flux_normalization="none",
                zero_flux_policy="identity",
                flux_floor_rel=1.0e-3,
            )

            self.assertEqual(report.floored_bin_count, 1)
            self.assertEqual(report.identity_bin_count, 0)
            self.assertEqual(identity_report.floored_bin_count, 1)
            self.assertEqual(identity_report.identity_bin_count, 0)
            actual = _read_sph_table(table)
            self.assertEqual(actual[0, 1], 1.0)
            np.testing.assert_allclose(actual, _read_sph_table(identity_table))

    def test_flux_floor_one_sided_zero_above_floor_still_raises(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "mgxs.h5"
            reference_flux = root / "reference_flux.csv"
            low_order_flux = root / "low_order_flux.csv"
            table = root / "next_sph.csv"
            write_mgxs(mgxs)
            reference_flux.write_text(
                "mixture,group,flux\n"
                "fuel,1,4.0\nfuel,2,3.0\n"
                "moderator,1,16.0\nmoderator,2,25.0\n",
                encoding="utf-8",
            )
            low_order_flux.write_text(
                "mixture,group,flux\n"
                "fuel,1,1.0\nfuel,2,0.0\n"
                "moderator,1,1.0\nmoderator,2,1.0\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(
                ValueError,
                r"CE/MG inconsistency: fuel g2 \(zero low-order flux\)",
            ):
                create_sph_update_table(
                    mgxs,
                    table,
                    reference_flux=reference_flux,
                    low_order_flux=low_order_flux,
                    zero_flux_policy="identity",
                    flux_floor_rel=1.0e-3,
                )

    def test_rejects_invalid_flux_floor_rel(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "mgxs.h5"
            reference_flux = root / "reference_flux.csv"
            low_order_flux = root / "low_order_flux.csv"
            table = root / "next_sph.csv"
            write_mgxs(mgxs)
            reference_flux.write_text(
                "mixture,group,flux\nfuel,1,1.0\nfuel,2,1.0\nmoderator,1,1.0\nmoderator,2,1.0\n",
                encoding="utf-8",
            )
            low_order_flux.write_text(
                "mixture,group,flux\nfuel,1,1.0\nfuel,2,1.0\nmoderator,1,1.0\nmoderator,2,1.0\n",
                encoding="utf-8",
            )

            for invalid in (0.0, 1.0, -0.5):
                with self.assertRaisesRegex(
                    ValueError,
                    "--flux-floor-rel must be strictly between 0 and 1",
                ):
                    create_sph_update_table(
                        mgxs,
                        table,
                        reference_flux=reference_flux,
                        low_order_flux=low_order_flux,
                        flux_floor_rel=invalid,
                    )

    def test_cli_openmc_sph_sidecar_accepts_flux_floor_rel(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "mgxs.h5"
            reference_flux = root / "openmc_ce_flux.h5"
            mg_flux = root / "openmc_mg_flux.h5"
            sidecar = root / "openmc_sph.h5"
            table = root / "openmc_sph.csv"
            summary = root / "openmc_sph_summary.json"
            write_mgxs(mgxs)
            _write_flux_source(
                reference_flux,
                "openmc_volume_flux",
                values=np.array([[4.0, 1.0e-5], [16.0, 25.0]]),
            )
            _write_flux_source(
                mg_flux,
                "openmc_mg_flux",
                values=np.array([[1.0, 0.0], [1.0, 1.0]]),
            )

            self.assertEqual(
                cli_main(
                    [
                        "make-openmc-sph-sidecar",
                        str(mgxs),
                        "-o",
                        str(sidecar),
                        "--reference-flux",
                        f"{reference_flux}::openmc_volume_flux",
                        "--mg-flux",
                        f"{mg_flux}::openmc_mg_flux",
                        "--table-output",
                        str(table),
                        "--sph-target",
                        "flux",
                        "--flux-normalization",
                        "none",
                        "--flux-floor-rel",
                        "1.0e-3",
                        "--summary-json",
                        str(summary),
                    ]
                ),
                0,
            )

            expected = np.array([[4.0, 1.0], [16.0, 25.0]])
            with h5py.File(sidecar, "r") as h5:
                np.testing.assert_allclose(h5["sph"][:], expected, rtol=1.0e-11)
            payload = json.loads(summary.read_text(encoding="utf-8"))
            self.assertEqual(payload["flux_floor_rel"], 1.0e-3)
            self.assertEqual(payload["floored_bin_count"], 1)
            self.assertEqual(payload["zero_flux_policy"], "reject")

    def test_freeze_groups_holds_group_at_previous_value_for_all_mixtures(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "mgxs.h5"
            reference_flux = root / "reference_flux.csv"
            low_order_flux = root / "low_order_flux.csv"
            previous_sph = root / "previous_sph.csv"
            table = root / "next_sph.csv"
            summary = root / "summary.json"
            write_mgxs(mgxs)
            reference_flux.write_text(
                "mixture,group,flux\n"
                "fuel,1,4.0\nfuel,2,9.0\n"
                "moderator,1,16.0\nmoderator,2,25.0\n",
                encoding="utf-8",
            )
            low_order_flux.write_text(
                "mixture,group,flux\n"
                "fuel,1,1.0\nfuel,2,1.0\n"
                "moderator,1,1.0\nmoderator,2,1.0\n",
                encoding="utf-8",
            )
            previous_sph.write_text(
                "mixture,g1,g2\nfuel,1.0,1.1\nmoderator,0.9,1.2\n",
                encoding="utf-8",
            )

            report = create_sph_update_table(
                mgxs,
                table,
                reference_flux=reference_flux,
                low_order_flux=low_order_flux,
                previous_sph=previous_sph,
                sph_target="flux",
                flux_normalization="none",
                freeze_groups=(2,),
                summary_json=summary,
            )

            self.assertEqual(report.freeze_groups, (2,))
            self.assertEqual(report.frozen_group_bin_count, 2)
            self.assertEqual(report.floored_bin_count, 0)
            self.assertEqual(report.identity_bin_count, 0)
            # The frozen group is matched by the table's own group labels.
            labeled = _read_sph_table_by_label(table)
            self.assertEqual(labeled[("fuel", 2)], 1.1)
            self.assertEqual(labeled[("moderator", 2)], 1.2)
            np.testing.assert_allclose(labeled[("fuel", 1)], 1.0 * 4.0, rtol=1.0e-11)
            np.testing.assert_allclose(labeled[("moderator", 1)], 0.9 * 16.0, rtol=1.0e-11)
            payload = json.loads(summary.read_text(encoding="utf-8"))
            self.assertEqual(payload["freeze_groups"], [2])
            self.assertEqual(payload["frozen_group_bin_count"], 2)

    def test_rejects_invalid_freeze_groups(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "mgxs.h5"
            reference_flux = root / "reference_flux.csv"
            low_order_flux = root / "low_order_flux.csv"
            table = root / "next_sph.csv"
            write_mgxs(mgxs)
            reference_flux.write_text(
                "mixture,group,flux\nfuel,1,1.0\nfuel,2,1.0\nmoderator,1,1.0\nmoderator,2,1.0\n",
                encoding="utf-8",
            )
            low_order_flux.write_text(
                "mixture,group,flux\nfuel,1,1.0\nfuel,2,1.0\nmoderator,1,1.0\nmoderator,2,1.0\n",
                encoding="utf-8",
            )

            for invalid in ((0,), (3,)):
                with self.assertRaisesRegex(
                    ValueError,
                    r"--freeze-groups group -?\d+ outside 1\.\.2",
                ):
                    create_sph_update_table(
                        mgxs,
                        table,
                        reference_flux=reference_flux,
                        low_order_flux=low_order_flux,
                        freeze_groups=invalid,
                    )

            with self.assertRaisesRegex(
                ValueError,
                "--freeze-groups must not contain duplicate groups",
            ):
                create_sph_update_table(
                    mgxs,
                    table,
                    reference_flux=reference_flux,
                    low_order_flux=low_order_flux,
                    freeze_groups=(1, 1),
                )

    def test_cli_openmc_sph_sidecar_accepts_freeze_groups(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "mgxs.h5"
            reference_flux = root / "openmc_ce_flux.h5"
            mg_flux = root / "openmc_mg_flux.h5"
            sidecar = root / "openmc_sph.h5"
            table = root / "openmc_sph.csv"
            summary = root / "openmc_sph_summary.json"
            write_mgxs(mgxs)
            _write_flux_source(
                reference_flux,
                "openmc_volume_flux",
                values=np.array([[4.0, 9.0], [16.0, 25.0]]),
            )
            _write_flux_source(
                mg_flux,
                "openmc_mg_flux",
                values=np.array([[1.0, 1.0], [1.0, 1.0]]),
            )

            self.assertEqual(
                cli_main(
                    [
                        "make-openmc-sph-sidecar",
                        str(mgxs),
                        "-o",
                        str(sidecar),
                        "--reference-flux",
                        f"{reference_flux}::openmc_volume_flux",
                        "--mg-flux",
                        f"{mg_flux}::openmc_mg_flux",
                        "--table-output",
                        str(table),
                        "--sph-target",
                        "flux",
                        "--flux-normalization",
                        "none",
                        "--freeze-groups",
                        "2",
                        "--summary-json",
                        str(summary),
                    ]
                ),
                0,
            )

            expected = np.array([[4.0, 1.0], [16.0, 1.0]])
            with h5py.File(sidecar, "r") as h5:
                np.testing.assert_allclose(h5["sph"][:], expected, rtol=1.0e-11)
            payload = json.loads(summary.read_text(encoding="utf-8"))
            self.assertEqual(payload["freeze_groups"], [2])
            self.assertEqual(payload["frozen_group_bin_count"], 2)

    def test_rate_target_update_targets_reaction_rates(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "mgxs.h5"
            reference_flux = root / "reference_flux.csv"
            low_order_flux = root / "low_order_flux.csv"
            previous_sph = root / "previous_sph.csv"
            table = root / "next_sph.csv"
            summary = root / "summary.json"
            write_mgxs(mgxs)
            reference_flux.write_text(
                "mixture,group,flux\n"
                "fuel,1,4.0\nfuel,2,9.0\n"
                "moderator,1,16.0\nmoderator,2,25.0\n",
                encoding="utf-8",
            )
            low_order_flux.write_text(
                "mixture,group,flux\n"
                "fuel,1,2.0\nfuel,2,3.0\n"
                "moderator,1,4.0\nmoderator,2,5.0\n",
                encoding="utf-8",
            )
            previous_sph.write_text(
                "mixture,g1,g2\nfuel,1.0,0.5\nmoderator,2.0,0.2\n",
                encoding="utf-8",
            )

            report = create_sph_update_table(
                mgxs,
                table,
                reference_flux=reference_flux,
                low_order_flux=low_order_flux,
                previous_sph=previous_sph,
                damping=0.5,
                sph_target="rate",
                flux_normalization="none",
                summary_json=summary,
            )

            # raw = phi_mg / (s_prev * phi_ref):
            #   fuel g1 2/(1*4)=0.5, fuel g2 3/(0.5*9)=2/3,
            #   moderator g1 4/(2*16)=1/8, moderator g2 5/(0.2*25)=1 (fixed point).
            expected = np.array(
                [
                    [1.0 * np.sqrt(0.5), 0.5 * np.sqrt(2.0 / 3.0)],
                    [2.0 * np.sqrt(0.125), 0.2 * 1.0],
                ]
            )
            self.assertEqual(report.sph_target, "rate")
            actual = _read_sph_table(table)
            np.testing.assert_allclose(actual, expected, rtol=1.0e-11)
            # A bin already at the rate-preserving fixed point stays put.
            self.assertEqual(actual[1, 1], 0.2)
            payload = json.loads(summary.read_text(encoding="utf-8"))
            self.assertEqual(payload["sph_target"], "rate")
            self.assertEqual(
                payload["formula"],
                "next_sph = previous_sph * "
                "(normalized_low_order_flux / (previous_sph * reference_flux)) ** damping",
            )

    def test_default_rate_target_matches_explicit_rate_target(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "mgxs.h5"
            reference_flux = root / "reference_flux.csv"
            low_order_flux = root / "low_order_flux.csv"
            previous_sph = root / "previous_sph.csv"
            default_table = root / "default_sph.csv"
            explicit_table = root / "explicit_sph.csv"
            diagnostic_table = root / "diagnostic_flux_sph.csv"
            write_mgxs(
                mgxs,
                h_factor={"fuel": np.asarray([1.0, 1.0])},
            )
            reference_flux.write_text(
                "mixture,group,flux\n"
                "fuel,1,4.0\nfuel,2,9.0\n"
                "moderator,1,16.0\nmoderator,2,25.0\n",
                encoding="utf-8",
            )
            low_order_flux.write_text(
                "mixture,group,flux\n"
                "fuel,1,1.0\nfuel,2,1.0\n"
                "moderator,1,1.0\nmoderator,2,1.0\n",
                encoding="utf-8",
            )
            previous_sph.write_text(
                "mixture,g1,g2\nfuel,1.1,1.2\nmoderator,0.9,1.0\n",
                encoding="utf-8",
            )

            default_report = create_sph_update_table(
                mgxs,
                default_table,
                reference_flux=reference_flux,
                low_order_flux=low_order_flux,
                previous_sph=previous_sph,
                damping=0.5,
            )
            create_sph_update_table(
                mgxs,
                explicit_table,
                reference_flux=reference_flux,
                low_order_flux=low_order_flux,
                previous_sph=previous_sph,
                damping=0.5,
                sph_target="rate",
                flux_normalization="auto",
            )
            create_sph_update_table(
                mgxs,
                diagnostic_table,
                reference_flux=reference_flux,
                low_order_flux=low_order_flux,
                previous_sph=previous_sph,
                damping=0.5,
                sph_target="flux",
                flux_normalization="auto",
            )

            self.assertEqual(default_report.sph_target, "rate")
            self.assertEqual(default_report.flux_normalization, "power")
            self.assertEqual(
                default_table.read_text(encoding="utf-8"),
                explicit_table.read_text(encoding="utf-8"),
            )
            self.assertNotEqual(
                default_table.read_text(encoding="utf-8"),
                diagnostic_table.read_text(encoding="utf-8"),
            )

            with self.assertRaisesRegex(ValueError, "--sph-target must be one of"):
                create_sph_update_table(
                    mgxs,
                    root / "bogus_sph.csv",
                    reference_flux=reference_flux,
                    low_order_flux=low_order_flux,
                    sph_target="bogus",
                )

    def test_rate_target_frozen_bins_keep_previous_value(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "mgxs.h5"
            reference_flux = root / "reference_flux.csv"
            low_order_flux = root / "low_order_flux.csv"
            previous_sph = root / "previous_sph.csv"
            table = root / "next_sph.csv"
            write_mgxs(mgxs)
            reference_flux.write_text(
                "mixture,group,flux\n"
                "fuel,1,4.0\nfuel,2,9.0\n"
                "moderator,1,16.0\nmoderator,2,25.0\n",
                encoding="utf-8",
            )
            low_order_flux.write_text(
                "mixture,group,flux\n"
                "fuel,1,2.0\nfuel,2,3.0\n"
                "moderator,1,4.0\nmoderator,2,5.0\n",
                encoding="utf-8",
            )
            previous_sph.write_text(
                "mixture,g1,g2\nfuel,1.0,1.1\nmoderator,2.0,1.2\n",
                encoding="utf-8",
            )

            report = create_sph_update_table(
                mgxs,
                table,
                reference_flux=reference_flux,
                low_order_flux=low_order_flux,
                previous_sph=previous_sph,
                sph_target="rate",
                flux_normalization="none",
                freeze_groups=(2,),
            )

            self.assertEqual(report.sph_target, "rate")
            self.assertEqual(report.frozen_group_bin_count, 2)
            actual = _read_sph_table(table)
            # Frozen group 2 keeps the previous values for all mixtures.
            self.assertEqual(actual[0, 1], 1.1)
            self.assertEqual(actual[1, 1], 1.2)
            # Active bins follow the rate update with damping 1.
            np.testing.assert_allclose(actual[0, 0], 1.0 * (2.0 / 4.0), rtol=1.0e-11)
            np.testing.assert_allclose(actual[1, 0], 2.0 * (4.0 / 32.0), rtol=1.0e-11)

    def test_cli_openmc_sph_sidecar_accepts_explicit_flux_diagnostic_target(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "mgxs.h5"
            reference_flux = root / "openmc_ce_flux.h5"
            mg_flux = root / "openmc_mg_flux.h5"
            sidecar = root / "openmc_sph.h5"
            table = root / "openmc_sph.csv"
            summary = root / "openmc_sph_summary.json"
            write_mgxs(mgxs)
            _write_flux_source(
                reference_flux,
                "openmc_volume_flux",
                values=np.array([[4.0, 9.0], [16.0, 25.0]]),
            )
            _write_flux_source(
                mg_flux,
                "openmc_mg_flux",
                values=np.array([[2.0, 3.0], [4.0, 5.0]]),
            )

            self.assertEqual(
                cli_main(
                    [
                        "make-openmc-sph-sidecar",
                        str(mgxs),
                        "-o",
                        str(sidecar),
                        "--reference-flux",
                        f"{reference_flux}::openmc_volume_flux",
                        "--mg-flux",
                        f"{mg_flux}::openmc_mg_flux",
                        "--table-output",
                        str(table),
                        "--sph-target",
                        "flux",
                        "--flux-normalization",
                        "none",
                        "--summary-json",
                        str(summary),
                    ]
                ),
                0,
            )

            expected = np.array([[2.0, 3.0], [4.0, 5.0]])
            with h5py.File(sidecar, "r") as h5:
                np.testing.assert_allclose(h5["sph"][:], expected, rtol=1.0e-11)
                self.assertEqual(
                    h5.attrs["sph_derivation"],
                    "ce-mg-flux-fixed-point",
                )
                self.assertEqual(h5.attrs["sph_target"], "flux")
                self.assertAlmostEqual(
                    float(h5.attrs["sph_max_update_residual"]),
                    4.0,
                )
            payload = json.loads(summary.read_text(encoding="utf-8"))
            self.assertEqual(payload["sph_target"], "flux")
            self.assertEqual(
                payload["formula"],
                "sph = previous_sph * "
                "(openmc_ce_reference_flux / normalized_openmc_mg_flux) ** damping",
            )


def _read_sph_table_by_label(path: Path) -> dict[tuple[str, int], float]:
    rows = path.read_text(encoding="utf-8").strip().splitlines()[1:]
    labeled: dict[tuple[str, int], float] = {}
    for row in rows:
        mixture, group, value = row.split(",")
        labeled[(mixture, int(group))] = float(value)
    return labeled


def _read_sph_table(path: Path) -> np.ndarray:
    rows = path.read_text(encoding="utf-8").strip().splitlines()[1:]
    return np.asarray([float(row.split(",")[2]) for row in rows]).reshape(2, 2)


def write_mgxs(path: Path, *, h_factor: dict[str, np.ndarray] | None = None) -> None:
    with h5py.File(path, "w") as h5:
        h5.attrs["energy_groups"] = 2
        h5.create_dataset("energy_bounds", data=np.array([1.0e-5, 1.0, 1.0e7]))
        mixtures = h5.create_group("mixtures")
        for name in ("fuel", "moderator"):
            group = mixtures.create_group(name)
            group.attrs["fissionable"] = name == "fuel"
            group.attrs["volume"] = 1.0
            group.create_dataset("total", data=np.ones(2))
            group.create_dataset("absorption", data=np.full(2, 0.1))
            group.create_dataset("fission", data=np.zeros(2))
            group.create_dataset("nu_fission", data=np.zeros(2))
            group.create_dataset("chi", data=np.zeros(2))
            group.create_dataset("scatter_matrix", data=np.zeros((1, 2, 2)))
            if h_factor and name in h_factor:
                group.create_dataset("kappa_fission", data=np.asarray(h_factor[name], dtype=float))


def _write_flux_source(
    path: Path,
    dataset_name: str,
    *,
    values: np.ndarray,
    std_dev: np.ndarray | None = None,
) -> None:
    with h5py.File(path, "w") as h5:
        dataset = h5.create_dataset(dataset_name, data=np.asarray(values, dtype=float))
        dataset.attrs["group_order"] = "mgxs_donjon"
        dataset.attrs["mixture_names"] = np.asarray(("fuel", "moderator"), dtype="S")
        dataset.attrs["energy_bounds_verified"] = True
        dataset.attrs["spatial_domain_order_verified"] = True
        if std_dev is not None:
            std_dataset = h5.create_dataset(
                f"{dataset_name}_std_dev",
                data=np.asarray(std_dev, dtype=float),
            )
            std_dataset.attrs["group_order"] = "mgxs_donjon"
            std_dataset.attrs["mixture_names"] = np.asarray(("fuel", "moderator"), dtype="S")
            std_dataset.attrs["std_dev_of"] = dataset_name


if __name__ == "__main__":
    unittest.main()
