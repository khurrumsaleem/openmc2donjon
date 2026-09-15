from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np

from openmc2donjon.openmc_provenance import file_sha256
from openmc2donjon.physical_sph_contract import (
    SPH_APPLY_BINDING_SCHEMA,
    SPH_SOURCE_BINDING_SCHEMA,
    physical_colorset_sph_issues,
    physical_sph_issues,
)


class PhysicalColorsetSphContractTests(unittest.TestCase):
    def test_generic_contract_accepts_one_or_many_declared_domains(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for domains in (1, 17, 91):
                path = root / f"physical_{domains}.h5"
                _write_colorset(path, domains=domains)
                self.assertEqual(physical_sph_issues(path), [])

    def test_accepts_seven_domain_applied_rate_sph_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "colorset.h5"
            _write_colorset(path)
            self.assertEqual(physical_colorset_sph_issues(path), [])

    def test_rejects_global_or_unconverged_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "colorset.h5"
            _write_colorset(path, kind="openmc-ce-mg-global", residual=0.08)
            issues = physical_colorset_sph_issues(path)
            self.assertTrue(any("global" in issue for issue in issues))
            self.assertTrue(any("not converged" in issue for issue in issues))

    def test_rejects_a_non_colorset_domain_count(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "fullcore.h5"
            _write_colorset(path, domains=91)
            issues = physical_colorset_sph_issues(path)
            self.assertTrue(any("exactly 7 domains" in issue for issue in issues))
            self.assertEqual(physical_sph_issues(path), [])

    def test_rejects_false_string_booleans(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "colorset.h5"
            _write_colorset(path)
            with h5py.File(path, "r+") as h5:
                h5.attrs["sph_applied"] = "false"
                h5.attrs["sph_real"] = "false"
            issues = physical_sph_issues(path)
            self.assertTrue(any("sph_applied=true" in issue for issue in issues))
            self.assertTrue(any("real OpenMC CE" in issue for issue in issues))

    def test_rejects_intermediate_openmc_mgxs_binding_mode(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "colorset.h5"
            _write_colorset(path)
            with h5py.File(path, "r+") as h5:
                h5.attrs["sph_apply_binding_mode"] = (
                    "openmc-mgxs-intermediate-unbound"
                )
                h5.attrs["sph_apply_sidecar_input_hash_verified"] = False
            issues = physical_sph_issues(path)
            self.assertTrue(
                any("converter-final-exact-input" in issue for issue in issues)
            )
            self.assertTrue(
                any("sidecar_input_hash_verified=true" in issue for issue in issues)
            )

    def test_rejects_nonfinite_or_negative_residual(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for label, residual in (("nan", np.nan), ("inf", np.inf), ("negative", -0.1)):
                path = root / f"{label}.h5"
                _write_colorset(path, residual=residual)
                issues = physical_sph_issues(path)
                self.assertTrue(
                    any("finite and non-negative" in issue for issue in issues),
                    (label, issues),
                )

    def test_rejects_missing_malformed_or_nonpositive_applied_sph(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            missing = root / "missing.h5"
            _write_colorset(missing)
            with h5py.File(missing, "r+") as h5:
                del h5["mixtures/domain_1/applied_sph"]
            self.assertTrue(
                any("applied_sph is required" in issue for issue in physical_sph_issues(missing))
            )

            malformed = root / "malformed.h5"
            _write_colorset(malformed)
            with h5py.File(malformed, "r+") as h5:
                del h5["mixtures/domain_1/applied_sph"]
                h5["mixtures/domain_1"].create_dataset("applied_sph", data=[1.0])
            self.assertTrue(
                any("must have shape (2,)" in issue for issue in physical_sph_issues(malformed))
            )

            invalid = root / "invalid.h5"
            _write_colorset(invalid)
            with h5py.File(invalid, "r+") as h5:
                h5["mixtures/domain_1/applied_sph"][:] = [np.nan, -1.0]
            issues = physical_sph_issues(invalid)
            self.assertTrue(any("finite values" in issue for issue in issues))
            self.assertTrue(any("positive values" in issue for issue in issues))

    def test_rejects_active_sph_and_any_adf_payload(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "colorset.h5"
            _write_colorset(path)
            with h5py.File(path, "r+") as h5:
                h5["mixtures/domain_1"].create_dataset("NSPH", data=np.ones(2))
                h5["mixtures/domain_2"].create_group("adf").create_dataset(
                    "FD_XMIN", data=np.ones(2)
                )
            issues = physical_sph_issues(path)
            self.assertTrue(any("active sph/SPH/NSPH" in issue for issue in issues))
            self.assertTrue(any("ADF payload" in issue for issue in issues))

    def test_single_state_requires_state_specific_vector_and_multistate_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            single = Path(tmp) / "single.h5"
            _write_colorset(single, domains=1, states=1)
            self.assertEqual(physical_sph_issues(single), [])
            with h5py.File(single, "r+") as h5:
                del h5["mixtures/domain_1/states/00000001/applied_sph"]
            self.assertTrue(
                any(
                    "states/00000001/applied_sph is required" in issue
                    for issue in physical_sph_issues(single)
                )
            )

            multiple = Path(tmp) / "multiple.h5"
            _write_colorset(multiple, domains=1, states=2)
            self.assertTrue(
                any(
                    "multi-state physical SPH is not supported" in issue
                    for issue in physical_sph_issues(multiple)
                )
            )

    def test_rejects_missing_malformed_or_mismatched_sha256_binding(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            missing = root / "missing.h5"
            _write_colorset(missing)
            with h5py.File(missing, "r+") as h5:
                del h5.attrs["sph_reference_flux_sha256"]
            self.assertTrue(
                any(
                    "sph_reference_flux_sha256 must be" in issue
                    for issue in physical_sph_issues(missing)
                )
            )

            malformed = root / "malformed.h5"
            _write_colorset(malformed)
            with h5py.File(malformed, "r+") as h5:
                h5.attrs["sph_mg_flux_sha256"] = "not-a-digest"
            self.assertTrue(
                any(
                    "sph_mg_flux_sha256 must be" in issue
                    for issue in physical_sph_issues(malformed)
                )
            )

            previous = root / "previous.h5"
            _write_colorset(previous)
            with h5py.File(previous, "r+") as h5:
                h5.attrs["sph_previous_sph_used"] = True
            self.assertTrue(
                any(
                    "sph_previous_sph_sha256 must be" in issue
                    for issue in physical_sph_issues(previous)
                )
            )

            mismatch = root / "mismatch.h5"
            _write_colorset(mismatch)
            with h5py.File(mismatch, "r") as h5:
                input_source = Path(str(h5.attrs["sph_apply_input_h5_path"]))
            input_source.write_bytes(b"tampered")
            self.assertTrue(
                any(
                    "does not match the existing apply-sph input HDF5" in issue
                    for issue in physical_sph_issues(mismatch)
                )
            )

    def test_uncertainty_gate_requires_coverage_limit_and_passing_observation(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            missing = root / "missing.h5"
            _write_colorset(missing)
            with h5py.File(missing, "r+") as h5:
                h5.attrs["sph_reference_flux_uncertainty_coverage"] = False
                del h5.attrs["sph_reference_flux_uncertainty_limit"]
            issues = physical_sph_issues(missing)
            self.assertTrue(any("uncertainty_coverage=true" in issue for issue in issues))
            self.assertTrue(any("uncertainty_limit must be" in issue for issue in issues))

            over = root / "over.h5"
            _write_colorset(over)
            with h5py.File(over, "r+") as h5:
                h5.attrs["sph_mg_flux_uncertainty_limit"] = 0.02
                h5.attrs["sph_mg_flux_uncertainty_observed_max_rel"] = 0.03
                h5.attrs["sph_mg_flux_uncertainty_pass"] = False
            issues = physical_sph_issues(over)
            self.assertTrue(any("exceeds its explicit limit" in issue for issue in issues))
            self.assertTrue(any("uncertainty_pass=true" in issue for issue in issues))

    def test_rejects_numerical_exemptions(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "colorset.h5"
            _write_colorset(path)
            with h5py.File(path, "r+") as h5:
                h5.attrs["sph_zero_flux_policy"] = "identity"
                h5.attrs["sph_identity_bin_count"] = 1
                h5.attrs["sph_floored_bin_count"] = 2
                h5.attrs["sph_frozen_group_bin_count"] = 3
                h5.attrs["sph_clipped_count"] = 4
            issues = physical_colorset_sph_issues(path)
            self.assertTrue(any("identity is forbidden" in issue for issue in issues))
            self.assertTrue(any("identity-substituted" in issue for issue in issues))
            self.assertTrue(any("flux-floored" in issue for issue in issues))
            self.assertTrue(any("frozen-group" in issue for issue in issues))
            self.assertTrue(any("clipped" in issue for issue in issues))

    def test_rejects_macrolib_filled_cross_section_bins(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "colorset.h5"
            _write_colorset(path)
            with h5py.File(path, "r+") as h5:
                h5["mixtures"]["domain_1"].attrs["zero_flux_filled_groups"] = [1, 2]
            issues = physical_colorset_sph_issues(path)
            self.assertTrue(any("macrolib-filled XS bins" in issue for issue in issues))


def _write_colorset(
    path: Path,
    *,
    domains: int = 7,
    kind: str = "openmc-ce-mg-rate",
    residual: float = 0.01,
    states: int | None = None,
) -> None:
    names = [f"domain_{index + 1}" for index in range(domains)]
    sources = _write_binding_sources(path)
    with h5py.File(path, "w") as h5:
        h5.attrs["energy_groups"] = 2
        h5.create_dataset("energy_bounds", data=np.asarray([0.0, 1.0, 2.0]))
        h5.create_dataset("mixture_names", data=np.asarray(names, dtype="S"))
        mixtures = h5.create_group("mixtures")
        for index, name in enumerate(names, start=1):
            group = mixtures.create_group(name)
            group.attrs["source_domain_index"] = index
            _write_minimal_xs(group)
            group.create_dataset("applied_sph", data=np.ones(2))
            if states is not None:
                state_group = group.create_group("states")
                for state_index in range(states):
                    state = state_group.create_group(f"{state_index + 1:08d}")
                    _write_minimal_xs(state)
                    state.create_dataset("applied_sph", data=np.ones(2))
        h5.attrs["sph_applied"] = True
        h5.attrs["sph_applied_source"] = str(sources["sidecar"])
        h5.attrs["sph_apply_operator"] = "divide-xs-by-nsph"
        h5.attrs["sph_kind"] = kind
        h5.attrs["sph_real"] = True
        h5.attrs["sph_derivation"] = "rate-preserving-ce-mg-fixed-point"
        h5.attrs["sph_target"] = "rate"
        h5.attrs["sph_flux_normalization"] = "power"
        h5.attrs["sph_zero_flux_policy"] = "reject"
        h5.attrs["sph_identity_bin_count"] = 0
        h5.attrs["sph_floored_bin_count"] = 0
        h5.attrs["sph_frozen_group_bin_count"] = 0
        h5.attrs["sph_clipped_count"] = 0
        h5.attrs["sph_max_update_residual"] = residual
        h5.attrs["sph_source_binding_schema"] = SPH_SOURCE_BINDING_SCHEMA
        h5.attrs["sph_apply_binding_schema"] = SPH_APPLY_BINDING_SCHEMA
        h5.attrs["sph_apply_binding_mode"] = "converter-final-exact-input"
        h5.attrs["sph_apply_sidecar_input_hash_verified"] = True
        _write_binding_attrs(h5, sources)
        for prefix, dataset in (
            ("sph_reference_flux", "/openmc_ce_flux"),
            ("sph_mg_flux", "/openmc_mg_flux"),
        ):
            h5.attrs[f"{prefix}_layout_verified"] = True
            h5.attrs[f"{prefix}_uncertainty_require_coverage"] = True
            h5.attrs[f"{prefix}_uncertainty_coverage"] = True
            h5.attrs[f"{prefix}_uncertainty_limit"] = 0.05
            h5.attrs[f"{prefix}_uncertainty_observed_max_rel"] = 0.01
            h5.attrs[f"{prefix}_uncertainty_pass"] = True
            h5.attrs[f"{prefix}_std_dev_dataset"] = f"{dataset}_std_dev"
        h5.attrs["sph_previous_sph_used"] = False


def _write_minimal_xs(group: h5py.Group) -> None:
    group.create_dataset("total", data=np.asarray([1.0, 0.8]))
    group.create_dataset("absorption", data=np.asarray([0.1, 0.2]))
    group.create_dataset("nu_fission", data=np.asarray([0.05, 0.0]))
    group.create_dataset("chi", data=np.asarray([1.0, 0.0]))
    group.create_dataset("scatter_matrix", data=np.zeros((1, 2, 2)))


def _write_binding_sources(path: Path) -> dict[str, Path]:
    sources = {
        "input_h5": path.with_name(f"{path.stem}_input.bin").resolve(),
        "reference_flux": path.with_name(f"{path.stem}_ce_flux.bin").resolve(),
        "mg_flux": path.with_name(f"{path.stem}_mg_flux.bin").resolve(),
        "sidecar": path.with_name(f"{path.stem}_sidecar.bin").resolve(),
    }
    for label, source in sources.items():
        source.write_bytes(f"{path.name}:{label}\n".encode())
    return sources


def _write_binding_attrs(h5: h5py.File, sources: dict[str, Path]) -> None:
    h5.attrs["sph_input_h5_path"] = str(sources["input_h5"])
    h5.attrs["sph_input_h5_sha256"] = file_sha256(sources["input_h5"])
    h5.attrs["sph_reference_flux_path"] = str(sources["reference_flux"])
    h5.attrs["sph_reference_flux_sha256"] = file_sha256(sources["reference_flux"])
    h5.attrs["sph_mg_flux_path"] = str(sources["mg_flux"])
    h5.attrs["sph_mg_flux_sha256"] = file_sha256(sources["mg_flux"])
    h5.attrs["sph_apply_input_h5_path"] = str(sources["input_h5"])
    h5.attrs["sph_apply_input_h5_sha256"] = file_sha256(sources["input_h5"])
    h5.attrs["sph_apply_sidecar_path"] = str(sources["sidecar"])
    h5.attrs["sph_apply_sidecar_sha256"] = file_sha256(sources["sidecar"])


if __name__ == "__main__":
    unittest.main()
