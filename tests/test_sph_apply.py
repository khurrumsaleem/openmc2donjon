from __future__ import annotations

from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import tempfile
import unittest

import h5py
import numpy as np

from openmc2donjon.cli import main as cli_main
from openmc2donjon.mgxs_input_contract import validate_input
from openmc2donjon.openmc_provenance import file_sha256
from openmc2donjon.openmc_sph_sidecar import create_openmc_sph_sidecar
from openmc2donjon.physical_sph_contract import physical_sph_issues
from openmc2donjon.sph_apply import (
    apply_sph_to_openmc_mgxs_hdf5,
    apply_sph_to_hdf5,
    apply_sph_to_mixture_arrays,
    apply_sph_to_scatter_matrix,
)


class SphApplyTests(unittest.TestCase):
    def test_transport_identity_survives_groupwise_sph_for_ordinary_and_nu_scatter(self) -> None:
        for nu_weighted in (False, True):
            with self.subTest(nu_weighted=nu_weighted), tempfile.TemporaryDirectory() as tmpdir:
                root = Path(tmpdir)
                mgxs, sidecar, output = root / "mgxs.h5", root / "sph.h5", root / "applied.h5"
                _write_transport_mgxs(mgxs, nu_weighted=nu_weighted)
                _write_sidecar(sidecar, mixture_names=("fuel",))
                with h5py.File(sidecar, "r+") as h5:
                    h5.attrs["sph_input_h5_sha256"] = file_sha256(mgxs)
                baseline = validate_input(mgxs, transport_p1_fail=1.0e-12)
                self.assertTrue(baseline.ok, baseline.issues)
                self.assertEqual(baseline.transport_p1_checked, 1)

                apply_sph_to_hdf5(mgxs, sph_source=sidecar, output_h5=output)
                applied = validate_input(output, transport_p1_fail=1.0e-12)
                self.assertTrue(applied.ok, applied.issues)
                self.assertEqual(applied.transport_p1_checked, 1)
                self.assertLess(float(applied.transport_p1_max_rel), 1.0e-12)
                with h5py.File(output, "r+") as h5:
                    np.testing.assert_array_equal(h5["openmc_volume_flux"][:], [[2.0, 1.0]])
                    np.testing.assert_allclose(
                        h5["mixtures/fuel/transport_total"][:], [0.9375, 2.2]
                    )
                    h5["mixtures/fuel/transport_total"][:] *= 1.2
                corrupted = validate_input(output, transport_p1_fail=5.0e-2)
                self.assertFalse(corrupted.ok)
                self.assertEqual(corrupted.transport_p1_checked, 1)
                self.assertTrue(any("transport_total/P1" in issue for issue in corrupted.issues))

    def test_sph_transport_audit_rejects_missing_or_invalid_applied_operator(self) -> None:
        cases = (
            ("missing factors", None, "divide-xs-by-nsph", "applied_sph vector"),
            ("wrong shape", [2.0], "divide-xs-by-nsph", "applied_sph vector"),
            ("matrix shape", [[2.0, 0.5]], "divide-xs-by-nsph", "applied_sph vector"),
            ("nonfinite", [np.nan, 0.5], "divide-xs-by-nsph", "applied_sph vector"),
            ("zero", [0.0, 0.5], "divide-xs-by-nsph", "applied_sph vector"),
            ("negative", [-1.0, 0.5], "divide-xs-by-nsph", "applied_sph vector"),
            ("missing operator", [2.0, 0.5], None, "sph_apply_operator"),
            ("unsupported", [2.0, 0.5], "multiply-xs", "sph_apply_operator"),
        )
        for label, factors, operator, issue in cases:
            with self.subTest(case=label), tempfile.TemporaryDirectory() as tmpdir:
                path = Path(tmpdir) / "applied.h5"
                _write_transport_mgxs(path, nu_weighted=True)
                with h5py.File(path, "r+") as h5:
                    h5.attrs["sph_applied"] = True
                    if operator is not None:
                        h5.attrs["sph_apply_operator"] = operator
                    if factors is not None:
                        h5["mixtures/fuel"].create_dataset("applied_sph", data=factors)
                report = validate_input(path, transport_p1_fail=5.0e-2)
                self.assertFalse(report.ok)
                self.assertEqual(report.transport_p1_skipped, 0)
                self.assertTrue(any(issue in detail for detail in report.issues), report.issues)

    def test_apply_sph_to_mixture_arrays_uses_nsph_divisor(self) -> None:
        sph = np.array([2.0, 0.5])
        scatter = np.array(
            [
                [[4.0, 6.0], [8.0, 10.0]],
                [[0.4, 0.6], [0.8, 1.0]],
            ]
        )
        datasets = {
            "total": np.array([10.0, 20.0]),
            "absorption": np.array([2.0, 3.0]),
            "reduced_absorption": np.array([1.0, -0.5]),
            "nu_fission": np.array([1.0, 0.0]),
            "total_std_dev": np.array([0.2, 0.4]),
            "scatter_matrix": scatter,
            "chi": np.array([0.25, 0.75]),
        }

        applied = apply_sph_to_mixture_arrays(datasets, sph)

        np.testing.assert_allclose(applied.datasets["total"], [5.0, 40.0])
        np.testing.assert_allclose(applied.datasets["absorption"], [1.0, 6.0])
        np.testing.assert_allclose(applied.datasets["reduced_absorption"], [0.5, -1.0])
        np.testing.assert_allclose(applied.datasets["nu_fission"], [0.5, 0.0])
        np.testing.assert_allclose(applied.datasets["total_std_dev"], [0.1, 0.8])
        np.testing.assert_allclose(applied.datasets["chi"], [0.25, 0.75])
        np.testing.assert_allclose(applied.datasets["scatter_matrix"][0], [[2.0, 3.0], [16.0, 20.0]])
        np.testing.assert_allclose(applied.datasets["scatter_matrix"][1], [[0.2, 0.3], [1.6, 2.0]])
        self.assertEqual(
            applied.scaled_names,
            (
                "total",
                "absorption",
                "reduced_absorption",
                "nu_fission",
                "total_std_dev",
                "scatter_matrix",
            ),
        )

    def test_apply_sph_to_scatter_matrix_supports_openmc_axis_order(self) -> None:
        values = np.array(
            [
                [[2.0, 20.0], [4.0, 40.0]],
                [[6.0, 60.0], [8.0, 80.0]],
            ]
        )

        corrected = apply_sph_to_scatter_matrix(
            values,
            np.array([2.0, 4.0]),
            scatter_axes="from,to,moment",
        )

        np.testing.assert_allclose(
            corrected,
            [
                [[1.0, 10.0], [2.0, 20.0]],
                [[1.5, 15.0], [2.0, 20.0]],
            ],
        )

    def test_apply_sph_to_hdf5_writes_corrected_copy_without_active_sph(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "mgxs.h5"
            sidecar = root / "sph.h5"
            output = root / "corrected.h5"
            _write_mgxs(mgxs)
            _write_sidecar(sidecar)

            report = apply_sph_to_hdf5(mgxs, sph_source=sidecar, output_h5=output)

            self.assertEqual(report.energy_groups, 2)
            self.assertEqual(report.mixture_names, ("fuel", "moderator"))
            self.assertEqual(report.operator, "divide-xs-by-nsph")
            self.assertEqual(report.scaled_dataset_count, 16)
            self.assertAlmostEqual(report.sph_min, 0.25)
            self.assertAlmostEqual(report.sph_max, 2.0)
            with h5py.File(output, "r") as h5:
                self.assertTrue(bool(h5.attrs["sph_applied"]))
                self.assertEqual(h5.attrs["sph_apply_operator"], "divide-xs-by-nsph")
                self.assertEqual(
                    h5.attrs["sph_derivation"],
                    "rate-preserving-ce-mg-fixed-point",
                )
                self.assertAlmostEqual(float(h5.attrs["sph_max_update_residual"]), 0.01)
                self.assertEqual(h5.attrs["sph_zero_flux_policy"], "reject")
                self.assertEqual(int(h5.attrs["sph_clipped_count"]), 0)
                self.assertEqual(
                    h5.attrs["sph_apply_input_h5_sha256"],
                    file_sha256(mgxs),
                )
                self.assertEqual(
                    h5.attrs["sph_apply_sidecar_sha256"],
                    file_sha256(sidecar),
                )
                fuel = h5["mixtures/fuel"]
                moderator = h5["mixtures/moderator"]
                np.testing.assert_allclose(fuel["total"][:], [5.0, 40.0])
                np.testing.assert_allclose(fuel["absorption"][:], [1.0, 6.0])
                np.testing.assert_allclose(fuel["nu_fission"][:], [0.5, 0.0])
                np.testing.assert_allclose(fuel["H-FACTOR"][:], [50.0, 400.0])
                np.testing.assert_allclose(fuel["total_std_dev"][:], [0.05, 0.4])
                np.testing.assert_allclose(fuel["chi"][:], [0.4, 0.6])
                np.testing.assert_allclose(fuel["scatter_matrix"][0], [[2.0, 3.0], [16.0, 20.0]])
                np.testing.assert_allclose(fuel["applied_sph"][:], [2.0, 0.5])
                self.assertNotIn("sph", fuel)
                np.testing.assert_allclose(moderator["total"][:], [5.0, 24.0])
                np.testing.assert_allclose(moderator["scatter_matrix"][0], [[4.0, 6.0], [32.0, 40.0]])
                np.testing.assert_allclose(moderator["applied_sph"][:], [1.0, 0.25])

    def test_apply_sph_to_hdf5_handles_one_state_per_mixture(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "states.h5"
            sidecar = root / "sph.h5"
            output = root / "corrected.h5"
            _write_state_mgxs(mgxs, state_count=1)
            _write_sidecar(sidecar, mixture_names=("fuel", "moderator"))

            apply_sph_to_hdf5(mgxs, sph_source=sidecar, output_h5=output)

            with h5py.File(output, "r") as h5:
                fuel = h5["mixtures/fuel"]
                np.testing.assert_allclose(fuel["applied_sph"][:], [2.0, 0.5])
                self.assertNotIn("sph", fuel)
                state = fuel["states/00000001"]
                np.testing.assert_allclose(state["total"][:], [5.0, 40.0])
                np.testing.assert_allclose(state["applied_sph"][:], [2.0, 0.5])
                self.assertNotIn("sph", state)

    def test_apply_sph_rejects_multi_state_without_state_specific_factors(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "states.h5"
            sidecar = root / "sph.h5"
            output = root / "corrected.h5"
            _write_state_mgxs(mgxs, state_count=2)
            _write_sidecar(sidecar, mixture_names=("fuel", "moderator"))

            with self.assertRaisesRegex(ValueError, "only one state per mixture"):
                apply_sph_to_hdf5(mgxs, sph_source=sidecar, output_h5=output)

    def test_apply_sph_rejects_nonpositive_factor(self) -> None:
        with self.assertRaisesRegex(ValueError, "positive"):
            apply_sph_to_mixture_arrays({"total": np.array([1.0, 2.0])}, np.array([1.0, 0.0]))

    def test_apply_sph_preserves_false_string_as_false(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "mgxs.h5"
            sidecar = root / "sph.h5"
            output = root / "corrected.h5"
            _write_mgxs(mgxs)
            _write_sidecar(sidecar)
            with h5py.File(sidecar, "r+") as h5:
                h5.attrs["sph_real"] = "false"

            apply_sph_to_hdf5(mgxs, sph_source=sidecar, output_h5=output)

            with h5py.File(output, "r") as h5:
                self.assertFalse(bool(h5.attrs["sph_real"]))

    def test_apply_sph_rejects_sidecar_bound_to_another_input(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "mgxs.h5"
            sidecar = root / "sph.h5"
            output = root / "corrected.h5"
            _write_mgxs(mgxs)
            _write_sidecar(sidecar)
            with h5py.File(sidecar, "r+") as h5:
                h5.attrs["sph_input_h5_sha256"] = "a" * 64

            with self.assertRaisesRegex(ValueError, "different input HDF5"):
                apply_sph_to_hdf5(mgxs, sph_source=sidecar, output_h5=output)

    def test_apply_sph_cli_writes_corrected_hdf5_and_summary(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "mgxs.h5"
            sidecar = root / "sph.h5"
            output = root / "corrected.h5"
            summary = root / "apply_summary.json"
            _write_mgxs(mgxs)
            _write_sidecar(sidecar)

            stdout = io.StringIO()
            with redirect_stdout(stdout):
                rc = cli_main(
                    [
                        "apply-sph",
                        str(mgxs),
                        "--sph-source",
                        str(sidecar),
                        "-o",
                        str(output),
                        "--summary-json",
                        str(summary),
                    ]
                )

            self.assertEqual(rc, 0)
            self.assertIn("SPH-applied MGXS", stdout.getvalue())
            self.assertIn("openmc2donjon_sph_apply_passed", stdout.getvalue())
            self.assertTrue(output.exists())
            payload = json.loads(summary.read_text(encoding="utf-8"))
            self.assertEqual(payload["schema"], "openmc2donjon.sph-apply.v1")
            self.assertEqual(payload["decision"], "openmc2donjon_sph_apply_passed")
            self.assertEqual(payload["operator"], "divide-xs-by-nsph")
            self.assertEqual(payload["mixtures"], ["fuel", "moderator"])
            with h5py.File(output, "r") as h5:
                np.testing.assert_allclose(h5["mixtures/fuel/total"][:], [5.0, 40.0])
                self.assertNotIn("sph", h5["mixtures/fuel"])

    def test_apply_sph_to_openmc_native_mgxs_uses_set_order(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "openmc_mgxs.h5"
            sidecar = root / "sph.h5"
            output = root / "openmc_mgxs_corrected.h5"
            _write_openmc_native_mgxs(mgxs)
            _write_sidecar(sidecar)
            with h5py.File(sidecar, "r+") as h5:
                # The sidecar was derived from a Converter-layout MGXS, not
                # this native setN iteration file.
                h5.attrs["sph_input_h5_sha256"] = "a" * 64

            report = apply_sph_to_openmc_mgxs_hdf5(mgxs, sph_source=sidecar, output_h5=output)

            self.assertEqual(report.input_format, "openmc-mgxs")
            self.assertEqual(report.binding_mode, "openmc-mgxs-intermediate-unbound")
            self.assertFalse(report.sidecar_input_hash_verified)
            self.assertEqual(report.scaled_dataset_count, 12)
            with h5py.File(output, "r") as h5:
                self.assertTrue(bool(h5.attrs["sph_applied"]))
                self.assertEqual(h5.attrs["sph_apply_input_format"], "openmc-mgxs")
                self.assertEqual(
                    h5.attrs["sph_apply_binding_mode"],
                    "openmc-mgxs-intermediate-unbound",
                )
                self.assertFalse(
                    bool(h5.attrs["sph_apply_sidecar_input_hash_verified"])
                )
                np.testing.assert_allclose(h5["set1/294K/total"][:], [5.0, 40.0])
                np.testing.assert_allclose(h5["set1/294K/absorption"][:], [1.0, 6.0])
                np.testing.assert_allclose(h5["set1/294K/nu-fission"][:], [0.5, 0.0])
                np.testing.assert_allclose(h5["set1/294K/kappa-fission"][:], [50.0, 400.0])
                np.testing.assert_allclose(h5["set1/294K/chi"][:], [0.4, 0.6])
                np.testing.assert_allclose(
                    h5["set1/294K/scatter_data/scatter_matrix"][:],
                    [2.0, 3.0, 16.0, 20.0],
                )
                np.testing.assert_allclose(
                    h5["set1/294K/scatter_data/multiplicity_matrix"][:],
                    [1.0, 1.0, 1.0, 1.0],
                )
                np.testing.assert_allclose(h5["set2/294K/total"][:], [5.0, 24.0])
                np.testing.assert_allclose(
                    h5["set2/294K/scatter_data/scatter_matrix"][:],
                    [4.0, 6.0, 32.0, 40.0],
                )
                np.testing.assert_allclose(
                    h5["set1/294K/applied_sph"][:],
                    [2.0, 0.5],
                )

    def test_apply_sph_cli_supports_openmc_native_mgxs(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "openmc_mgxs.h5"
            sidecar = root / "sph.h5"
            output = root / "openmc_mgxs_corrected.h5"
            _write_openmc_native_mgxs(mgxs)
            _write_sidecar(sidecar)

            with redirect_stdout(io.StringIO()):
                rc = cli_main(
                    [
                        "apply-sph",
                        str(mgxs),
                        "--input-format",
                        "openmc-mgxs",
                        "--sph-source",
                        str(sidecar),
                        "-o",
                        str(output),
                    ]
                )

            self.assertEqual(rc, 0)
            with h5py.File(output, "r") as h5:
                np.testing.assert_allclose(h5["set2/294K/total"][:], [5.0, 24.0])

    def test_verified_sidecar_apply_produces_strict_physical_handoff(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            mgxs = root / "mgxs.h5"
            reference_flux = root / "reference_flux.h5"
            mg_flux = root / "mg_flux.h5"
            sidecar = root / "sph.h5"
            corrected = root / "corrected.h5"
            _write_mgxs(mgxs)
            _write_verified_flux(reference_flux, "reference_flux")
            _write_verified_flux(mg_flux, "mg_flux")

            create_openmc_sph_sidecar(
                mgxs,
                sidecar,
                reference_flux=f"{reference_flux}::reference_flux",
                mg_flux=f"{mg_flux}::mg_flux",
                require_reference_flux_std_dev=True,
                max_reference_flux_std_dev_rel=0.02,
                require_mg_flux_std_dev=True,
                max_mg_flux_std_dev_rel=0.02,
            )
            apply_sph_to_hdf5(
                mgxs,
                sph_source=sidecar,
                output_h5=corrected,
            )

            self.assertEqual(physical_sph_issues(corrected), [])


def _write_transport_mgxs(path: Path, *, nu_weighted: bool) -> None:
    with h5py.File(path, "w") as h5:
        h5.attrs["energy_groups"] = 2
        h5.attrs["legendre_order"] = 1
        h5.attrs["openmc_scatter_mgxs_type"] = (
            "consistent nu-scatter matrix" if nu_weighted else "scatter matrix"
        )
        h5.attrs["openmc_transport_mgxs_type"] = "nu-transport" if nu_weighted else "transport"
        h5.create_dataset("energy_bounds", data=[1.0e-5, 1.0, 2.0])
        h5.create_dataset("mixture_names", data=np.asarray(["fuel"], dtype="S"))
        phi = np.array([2.0, 1.0])
        flux = h5.create_dataset("openmc_volume_flux", data=phi[None, :])
        flux.attrs["mixture_names"] = np.asarray(["fuel"], dtype="S")
        flux.attrs["group_order"] = "mgxs_donjon"
        group = h5.create_group("mixtures/fuel")
        group.attrs["volume"] = 1.0
        group.attrs["fissionable"] = False
        group.attrs["scatter_axes"] = "moment,from,to"
        p0 = np.array([[1.0, 0.4], [0.2, 1.1]]) if nu_weighted else np.array([[0.6, 0.4], [0.2, 0.8]])
        p1 = np.array([[0.1, 0.4], [0.05, 0.1]])
        total = np.array([2.0, 2.0])
        group.create_dataset("total", data=total)
        group.create_dataset("absorption", data=[1.0, 1.0])
        if nu_weighted:
            group.create_dataset("reduced_absorption", data=total - p0.sum(axis=1))
        group.create_dataset("scatter_matrix", data=np.stack([p0, p1]))
        group.create_dataset(
            "transport_total", data=total - np.sum(phi[:, None] * p1, axis=0) / phi
        )
        for name in ("fission", "nu_fission", "chi"):
            group.create_dataset(name, data=[0.0, 0.0])


def _write_mgxs(path: Path) -> None:
    with h5py.File(path, "w") as h5:
        h5.attrs["energy_groups"] = 2
        h5.create_dataset("energy_bounds", data=[0.0, 1.0, 2.0])
        h5.create_dataset(
            "mixture_names",
            data=np.asarray(("fuel", "moderator"), dtype="S"),
        )
        mixtures = h5.create_group("mixtures")
        _write_mix(mixtures.create_group("fuel"), fissionable=True)
        _write_mix(mixtures.create_group("moderator"), fissionable=False, total=(5.0, 6.0))


def _write_state_mgxs(path: Path, *, state_count: int) -> None:
    with h5py.File(path, "w") as h5:
        h5.attrs["energy_groups"] = 2
        h5.create_dataset("energy_bounds", data=[0.0, 1.0, 2.0])
        mixtures = h5.create_group("mixtures")
        for mixture_name, fissionable in (("fuel", True), ("moderator", False)):
            states = mixtures.create_group(mixture_name).create_group("states")
            for state_index in range(state_count):
                _write_mix(
                    states.create_group(f"{state_index + 1:08d}"),
                    fissionable=fissionable,
                )


def _write_mix(group, *, fissionable: bool, total: tuple[float, float] = (10.0, 20.0)) -> None:
    group.attrs["fissionable"] = fissionable
    group.attrs["scatter_axes"] = "moment,from,to"
    group.create_dataset("total", data=np.array(total))
    group.create_dataset("total_std_dev", data=np.array([0.1, 0.2]))
    group.create_dataset("absorption", data=np.array([2.0, 3.0]))
    group.create_dataset("fission", data=np.array([0.2, 0.0]))
    group.create_dataset("nu_fission", data=np.array([1.0, 0.0]))
    group.create_dataset("H-FACTOR", data=np.array([100.0, 200.0]))
    group.create_dataset("chi", data=np.array([0.4, 0.6]))
    group.create_dataset(
        "scatter_matrix",
        data=np.array([[[4.0, 6.0], [8.0, 10.0]]]),
    )
    group.create_dataset(
        "scatter_matrix_std_dev",
        data=np.array([[[0.04, 0.06], [0.08, 0.10]]]),
    )
    group.create_dataset("sph", data=np.array([9.0, 9.0]))


def _write_sidecar(path: Path, *, mixture_names: tuple[str, ...] = ("fuel", "moderator")) -> None:
    values = np.array([[2.0, 0.5], [1.0, 0.25]])
    with h5py.File(path, "w") as h5:
        h5.attrs["sph_kind"] = "openmc-ce-mg"
        h5.attrs["sph_real"] = True
        h5.attrs["sph_applied"] = False
        h5.attrs["sph_derivation"] = "rate-preserving-ce-mg-fixed-point"
        h5.attrs["sph_target"] = "rate"
        h5.attrs["sph_flux_normalization"] = "power"
        h5.attrs["sph_raw_update_minimum"] = 0.99
        h5.attrs["sph_raw_update_maximum"] = 1.01
        h5.attrs["sph_max_update_residual"] = 0.01
        h5.attrs["sph_zero_flux_policy"] = "reject"
        h5.attrs["sph_identity_bin_count"] = 0
        h5.attrs["sph_floored_bin_count"] = 0
        h5.attrs["sph_frozen_group_bin_count"] = 0
        h5.attrs["sph_clipped_count"] = 0
        dataset = h5.create_dataset("sph", data=values[: len(mixture_names)])
        dataset.attrs["mixture_names"] = np.asarray(mixture_names, dtype="S")
        dataset.attrs["group_order"] = "mgxs_donjon"


def _write_openmc_native_mgxs(path: Path) -> None:
    with h5py.File(path, "w") as h5:
        h5.attrs["filetype"] = np.bytes_("mgxs")
        h5.attrs["energy_groups"] = 2
        h5.attrs["group structure"] = np.array([0.0, 1.0, 2.0])
        h5.create_group("settings")
        _write_openmc_set(h5.create_group("set1"), total=(10.0, 20.0))
        _write_openmc_set(h5.create_group("set2"), total=(5.0, 6.0))


def _write_verified_flux(path: Path, dataset_name: str) -> None:
    with h5py.File(path, "w") as h5:
        dataset = h5.create_dataset(dataset_name, data=np.ones((2, 2)))
        dataset.attrs["mixture_names"] = np.asarray(("fuel", "moderator"), dtype="S")
        dataset.attrs["group_order"] = "mgxs_donjon"
        dataset.attrs["energy_bounds_verified"] = True
        dataset.attrs["spatial_domain_order_verified"] = True
        std_dev = h5.create_dataset(f"{dataset_name}_std_dev", data=np.full((2, 2), 0.01))
        std_dev.attrs["mixture_names"] = np.asarray(("fuel", "moderator"), dtype="S")
        std_dev.attrs["group_order"] = "mgxs_donjon"
        std_dev.attrs["std_dev_of"] = dataset_name


def _write_openmc_set(group, *, total: tuple[float, float]) -> None:
    group.attrs["scatter_format"] = np.bytes_("histogram")
    group.attrs["scatter_shape"] = np.bytes_("[G][G'][Order]")
    group.create_group("kTs").create_dataset("294K", data=294.0)
    temperature = group.create_group("294K")
    temperature.create_dataset("total", data=np.array(total))
    temperature.create_dataset("absorption", data=np.array([2.0, 3.0]))
    temperature.create_dataset("fission", data=np.array([0.2, 0.0]))
    temperature.create_dataset("nu-fission", data=np.array([1.0, 0.0]))
    temperature.create_dataset("kappa-fission", data=np.array([100.0, 200.0]))
    temperature.create_dataset("chi", data=np.array([0.4, 0.6]))
    scatter = temperature.create_group("scatter_data")
    scatter.create_dataset("g_min", data=np.array([1, 1]))
    scatter.create_dataset("g_max", data=np.array([2, 2]))
    scatter.create_dataset("multiplicity_matrix", data=np.array([1.0, 1.0, 1.0, 1.0]))
    scatter.create_dataset("scatter_matrix", data=np.array([4.0, 6.0, 8.0, 10.0]))
