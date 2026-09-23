"""SPH application preserves metadata and cannot publish partially scaled XS."""

import h5py
import numpy as np
import pytest

from openmc2donjon import sph_apply
from openmc2donjon.collapse_substitution import RECORD_DATASET
from openmc2donjon.energy_collapse import collapse_energy_groups
from openmc2donjon.mgxs_input_contract import validate_input
from openmc2donjon.mixture_collapse import collapse_components
from openmc2donjon.openmc_provenance import file_sha256, read_openmc_provenance
from openmc2donjon.physical_sph_contract import physical_sph_issues
from tests.test_collapse_contract import _mark_filled
from tests.test_collapse_safety import _bind
from tests.test_collapse_scatter_layout import _source
from tests.test_sph_apply import _write_mgxs, _write_openmc_native_mgxs, _write_sidecar


@pytest.mark.parametrize("operation", ["component", "energy"])
@pytest.mark.parametrize("bound", [False, True])
def test_collapsed_donor_metadata_survives_sph_without_physical_acceptance(tmp_path, operation, bound):
    source, collapsed, sidecar, output = [tmp_path / name for name in (
        "source.h5", "collapsed.h5", "sph.h5", "output.h5",
    )]
    ngroups = 2 if operation == "component" else 4
    _source(source, ngroups=ngroups, order=1)
    with h5py.File(source, "r+") as h5:
        _mark_filled(h5["mixtures/LOW"], [ngroups - 1], "diagnostic-donor.h5")
    if bound:
        _bind(source, tmp_path)
    if operation == "component":
        collapse_components(source, collapsed, groups=[("BOTH", ("HIGH", "LOW"))])
        names = ("BOTH",)
    else:
        collapse_energy_groups(source, collapsed, groups=((1, 2), (3, 4)))
        names = ("HIGH", "LOW")
    _write_sidecar(sidecar, mixture_names=names)
    with h5py.File(sidecar, "r+") as h5:
        h5.attrs["sph_input_h5_sha256"] = file_sha256(collapsed)
        factors = h5["sph"][:]
    before, sidecar_before = collapsed.read_bytes(), sidecar.read_bytes()
    output.write_bytes(b"previous output")
    sph_apply.apply_sph_to_hdf5(collapsed, sph_source=sidecar, output_h5=output, force=True)
    assert collapsed.read_bytes() == before
    assert sidecar.read_bytes() == sidecar_before
    assert validate_input(output).ok
    with h5py.File(collapsed) as src, h5py.File(output) as dst:
        assert dst.attrs["sph_applied"]
        for i, name in enumerate(names):
            original, corrected = src[f"mixtures/{name}"], dst[f"mixtures/{name}"]
            np.testing.assert_allclose(corrected["total"][:], original["total"][:] / factors[i])
            for field in ("scatter_matrix", "scatter_matrix_std_dev"):
                np.testing.assert_allclose(
                    corrected[field][:], original[field][:] / factors[i][None, :, None],
                )
                assert dict(corrected[field].attrs) == dict(original[field].attrs)
            if RECORD_DATASET in original:
                assert original[RECORD_DATASET][()] == corrected[RECORD_DATASET][()]
                np.testing.assert_array_equal(original.attrs["zero_flux_filled_groups"],
                                              corrected.attrs["zero_flux_filled_groups"])
                assert (
                    original.attrs["zero_flux_fill_semantics"] == corrected.attrs["zero_flux_fill_semantics"]
                )
    # Successful arithmetic application is NOT physical acceptance of donor-filled XS.
    assert any("zero macrolib-filled XS bins" in issue for issue in physical_sph_issues(output))
    if bound:
        assert read_openmc_provenance(output)["integrity"]["ok"]


def test_sph_leaves_unscaled_metadata_and_observables_untouched(tmp_path):
    source, sidecar, output = (tmp_path / name for name in ("in.h5", "sph.h5", "out.h5"))
    _write_mgxs(source)
    _write_sidecar(sidecar)
    with h5py.File(source, "r+") as h5:
        group = h5["mixtures/fuel"]
        group.create_dataset("note", data="not a cross section", dtype=h5py.string_dtype())
        group.create_dataset("scalar_metadata", data=42)
        group.create_dataset("inverse_velocity", data=[.1, .2])
        group["total_std_dev"].attrs["uncertainty_method"] = "conditional-test-bound"
        group["total"].attrs["units"] = "cm^-1"
        group["chi"].attrs["role"] = "outgoing-spectrum"
    sph_apply.apply_sph_to_hdf5(source, sph_source=sidecar, output_h5=output)
    with h5py.File(source) as src, h5py.File(output) as dst:
        for name in ("note", "scalar_metadata", "inverse_velocity", "chi"):
            left, right = src[f"mixtures/fuel/{name}"], dst[f"mixtures/fuel/{name}"]
            np.testing.assert_array_equal(left[()], right[()])
            assert dict(left.attrs) == dict(right.attrs)
        for name in ("total", "total_std_dev"):
            assert dict(src[f"mixtures/fuel/{name}"].attrs) == dict(dst[f"mixtures/fuel/{name}"].attrs)


def _case(tmp_path, native):
    source, sidecar, output = (tmp_path / name for name in ("in.h5", "sph.h5", "out.h5"))
    (_write_openmc_native_mgxs if native else _write_mgxs)(source)
    _write_sidecar(sidecar)
    apply = sph_apply.apply_sph_to_openmc_mgxs_hdf5 if native else sph_apply.apply_sph_to_hdf5
    return source, sidecar, output, apply


@pytest.mark.parametrize("native", [False, True])
@pytest.mark.parametrize("existing", [False, True])
@pytest.mark.parametrize("failure", ["second_region", "write", "refresh", "publish"])
def test_late_sph_failures_preserve_inputs_and_old_output(tmp_path, monkeypatch, native, existing, failure):
    source, sidecar, output, apply = _case(tmp_path, native)
    if failure == "second_region":
        path = "set2/294K/total" if native else "mixtures/moderator/total"
        with h5py.File(source, "r+") as h5:
            del h5[path]
            h5.create_dataset(path, data=[1.0])
    before, sidecar_before = source.read_bytes(), sidecar.read_bytes()
    if existing:
        output.write_bytes(b"previous output")

    def fail(*args, **kwargs):
        raise RuntimeError("injected failure")

    if failure == "write":
        original = sph_apply._replace_dataset

        def write_then_fail(*args, **kwargs):
            original(*args, **kwargs)
            fail()

        monkeypatch.setattr(sph_apply, "_replace_dataset", write_then_fail)
    elif failure == "refresh":
        monkeypatch.setattr(sph_apply, "refresh_openmc_provenance_after_hdf5_mutation", fail)
    elif failure == "publish":
        monkeypatch.setattr(sph_apply.os, "replace" if existing else "link", fail)
    with pytest.raises((RuntimeError, ValueError)):
        apply(source, sph_source=sidecar, output_h5=output, force=existing)
    assert source.read_bytes() == before
    assert sidecar.read_bytes() == sidecar_before
    if existing:
        assert output.read_bytes() == b"previous output"
    else:
        assert not output.exists()
    assert not list(tmp_path.glob(".openmc2donjon-sph-*"))


@pytest.mark.parametrize("native", [False, True])
@pytest.mark.parametrize("target", ["input", "sidecar"])
@pytest.mark.parametrize("alias", ["same", "hardlink", "symlink"])
def test_sph_output_cannot_overwrite_input_or_sidecar(tmp_path, native, target, alias):
    source, sidecar, output, apply = _case(tmp_path, native)
    protected = source if target == "input" else sidecar
    before, sidecar_before = source.read_bytes(), sidecar.read_bytes()
    if alias == "same":
        output = protected
    elif alias == "hardlink":
        output.hardlink_to(protected)
    else:
        output.symlink_to(protected)
    with pytest.raises(ValueError, match="different from"):
        apply(source, sph_source=sidecar, output_h5=output, force=True)
    assert source.read_bytes() == before
    assert sidecar.read_bytes() == sidecar_before


@pytest.mark.parametrize("native", [False, True])
@pytest.mark.parametrize("change", ["input", "sidecar", "destination"])
def test_concurrent_file_change_is_not_silently_published(tmp_path, monkeypatch, native, change):
    source, sidecar, output, apply = _case(tmp_path, native)
    function = "_apply_to_openmc_macro_group" if native else "_apply_to_mixture_group"
    original = getattr(sph_apply, function)
    changed = False

    def apply_and_change(*args, **kwargs):
        nonlocal changed
        result = original(*args, **kwargs)
        if not changed:
            if change == "destination":
                output.write_bytes(b"concurrent output")
            else:
                path = source if change == "input" else sidecar
                with h5py.File(path, "r+") as h5:
                    h5.attrs["concurrent_edit"] = True
            changed = True
        return result

    monkeypatch.setattr(sph_apply, function, apply_and_change)
    with pytest.raises((ValueError, FileExistsError)):
        apply(source, sph_source=sidecar, output_h5=output)
    if change == "destination":
        assert output.read_bytes() == b"concurrent output"
    else:
        assert not output.exists()
    assert not list(tmp_path.glob(".openmc2donjon-sph-*"))
