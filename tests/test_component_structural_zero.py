"""Structural nonfission zeros are not missing Monte Carlo statistics."""

import shutil

import h5py
import numpy as np
import pytest

from openmc2donjon.mgxs_input_contract import validate_input
from openmc2donjon.mgxs_input_uncertainty import UncertaintyConfig
from openmc2donjon.mixture_collapse import collapse_components
from openmc2donjon.energy_collapse import collapse_energy_groups
from openmc2donjon.sph_apply import apply_sph_to_hdf5
from tests.test_collapse_safety import _add_uncertainty
from tests.test_collapse_scatter_layout import _source
from tests.test_sph_apply import _write_sidecar


FAMILY = ("fission", "nu_fission", "chi")


def _case(path):
    _source(path)
    _add_uncertainty(path)
    with h5py.File(path, "r+") as h5:
        low = h5["mixtures/LOW"]
        low.attrs["fissionable"] = False
        for name in FAMILY:
            low[name][:] = 0
            del low[name + "_std_dev"]
        low["kappa_fission"][:] = 0
        low["kappa_fission_std_dev"][:] = 0


@pytest.mark.parametrize("members", [("HIGH", "LOW"), ("LOW", "HIGH")])
@pytest.mark.parametrize("declaration", [False, 0])
def test_zero_placeholder_exemption_retains_fuel_errors_and_chi(tmp_path, members, declaration):
    source, explicit, output, reference = [tmp_path / name for name in (
        "source.h5", "explicit.h5", "output.h5", "reference.h5",
    )]
    _case(source)
    with h5py.File(source, "r+") as h5:
        h5["mixtures/LOW"].attrs["fissionable"] = declaration
    cfg = UncertaintyConfig(require_coverage=True)
    report = validate_input(source, uncertainty=cfg)
    assert report.ok, report.issues
    shutil.copy2(source, explicit)
    with h5py.File(explicit, "r+") as h5:
        for name in FAMILY:
            h5["mixtures/LOW"].create_dataset(name + "_std_dev", data=np.zeros(4))
    before = source.read_bytes()
    collapse_components(source, output, groups=[("BOTH", members)])
    collapse_components(explicit, reference, groups=[("BOTH", members)])
    assert source.read_bytes() == before
    report = validate_input(output, uncertainty=cfg)
    assert report.ok, report.issues
    with h5py.File(source) as src, h5py.File(output) as out, h5py.File(reference) as ref:
        assert out["mixtures/BOTH"].attrs["fissionable"]
        for name in FAMILY:
            for field in (name, name + "_std_dev"):
                np.testing.assert_allclose(out[f"mixtures/BOTH/{field}"][:], ref[f"mixtures/BOTH/{field}"][:])
            std = out[f"mixtures/BOTH/{name}_std_dev"]
            assert tuple(std.attrs["component_structural_zero_sources"]) == (b"/mixtures/LOW",)
            assert std.attrs["component_structural_zero_basis"] == (
                "declared-nonfissionable-and-zero-fission-family"
            )
        phi = src["openmc_volume_flux"][:]
        for name in ("fission", "nu_fission"):
            expected = src[f"mixtures/HIGH/{name}_std_dev"][:] * phi[0] / phi.sum(axis=0)
            np.testing.assert_allclose(out[f"mixtures/BOTH/{name}_std_dev"][:], expected)
        assert np.any(out["mixtures/BOTH/chi_std_dev"][:] > 0)


@pytest.mark.parametrize("fault", ["no_flag", "fissile", "string", "fission", "nu_fission", "chi", "missing"])
def test_zero_tallies_or_incomplete_declarations_are_not_structural(tmp_path, fault):
    source, output = tmp_path / "in.h5", tmp_path / "out.h5"
    _case(source)
    with h5py.File(source, "r+") as h5:
        low = h5["mixtures/LOW"]
        if fault == "no_flag":
            del low.attrs["fissionable"]
        elif fault in {"fissile", "string"}:
            low.attrs["fissionable"] = True if fault == "fissile" else "false"
        elif fault == "missing":
            del low["chi"]
        else:
            low[fault][0] = 1e-30
    before = source.read_bytes()
    output.write_bytes(b"old output")
    with pytest.raises(ValueError, match="partial fission_std_dev"):
        collapse_components(source, output, groups=[("BOTH", ("HIGH", "LOW"))], force=True)
    assert source.read_bytes() == before
    assert output.read_bytes() == b"old output"


@pytest.mark.parametrize("field", ["total", "absorption", "kappa_fission", "inverse_velocity"])
def test_structural_exemption_does_not_extend_to_other_responses(tmp_path, field):
    source, output = tmp_path / "in.h5", tmp_path / "out.h5"
    _case(source)
    with h5py.File(source, "r+") as h5:
        if field == "inverse_velocity":
            for group in h5["mixtures"].values():
                group.create_dataset(field, data=np.zeros(4))
            h5["mixtures/HIGH"].create_dataset(field + "_std_dev", data=np.ones(4) * .01)
        else:
            del h5[f"mixtures/LOW/{field}_std_dev"]
    with pytest.raises(ValueError, match=f"partial {field}_std_dev"):
        collapse_components(source, output, groups=[("BOTH", ("HIGH", "LOW"))])
    assert not output.exists()


def test_entirely_absent_uncertainty_remains_absent(tmp_path):
    source, output = tmp_path / "in.h5", tmp_path / "out.h5"
    _case(source)
    with h5py.File(source, "r+") as h5:
        for field in FAMILY:
            del h5[f"mixtures/HIGH/{field}_std_dev"]
    collapse_components(source, output, groups=[("BOTH", ("HIGH", "LOW"))])
    with h5py.File(output) as h5:
        assert all(field + "_std_dev" not in h5["mixtures/BOTH"] for field in FAMILY)


def test_unknown_chi_error_is_not_silently_dropped(tmp_path):
    source, output = tmp_path / "in.h5", tmp_path / "out.h5"
    _source(source)
    _add_uncertainty(source)
    with h5py.File(source, "r+") as h5:
        del h5["mixtures/LOW/chi_std_dev"]
    with pytest.raises(ValueError, match="partial chi_std_dev"):
        collapse_components(source, output, groups=[("BOTH", ("HIGH", "LOW"))])
    assert not output.exists()


def test_structural_zero_audit_survives_energy_collapse_and_sph(tmp_path):
    source, spatial, energy, sidecar, applied = [tmp_path / name for name in (
        "source.h5", "spatial.h5", "energy.h5", "sph.h5", "applied.h5",
    )]
    _case(source)
    collapse_components(source, spatial, groups=[("BOTH", ("HIGH", "LOW"))])
    collapse_energy_groups(spatial, energy, groups=((1, 2), (3, 4)))
    _write_sidecar(sidecar, mixture_names=("BOTH",))
    apply_sph_to_hdf5(energy, sph_source=sidecar, output_h5=applied)
    report = validate_input(applied, uncertainty=UncertaintyConfig(require_coverage=True))
    assert report.ok, report.issues
    with h5py.File(energy) as before, h5py.File(applied) as after:
        for name in FAMILY:
            original = before[f"mixtures/BOTH/{name}_std_dev"]
            corrected = after[f"mixtures/BOTH/{name}_std_dev"]
            for key in original.attrs:
                np.testing.assert_array_equal(original.attrs[key], corrected.attrs[key])
            factor = np.array([2., .5]) if name != "chi" else np.ones(2)
            np.testing.assert_allclose(corrected[:], original[:] / factor)
