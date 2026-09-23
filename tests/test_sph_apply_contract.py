"""SPH application shares Converter layouts and reports only verified bindings."""

import h5py
import numpy as np
import pytest

from openmc2donjon.mgxs_input_contract import validate_input
from openmc2donjon.multicompo import read_mgxs_hdf5
from openmc2donjon.openmc_provenance import file_sha256
from openmc2donjon.openmc_sph_sidecar import create_openmc_sph_sidecar
from openmc2donjon.physical_sph_contract import physical_sph_issues
from openmc2donjon.sph_apply import (
    apply_sph_to_hdf5,
    print_report,
    summary_payload,
)
from tests.test_collapse_scatter_layout import _source
from tests.test_sph_apply import _write_mgxs, _write_sidecar, _write_verified_flux


def _move_datasets_to_states(path):
    with h5py.File(path, "r+") as h5:
        # A root reference-flux field is not supported in the states/ contract.
        del h5["openmc_volume_flux"]
        for group in h5["mixtures"].values():
            names = list(group)
            state = group.create_group("states/00000001")
            for name in names:
                h5.move(f"{group.name}/{name}", f"{state.name}/{name}")


def _bound_sidecar(path, source, ngroups):
    factors = np.stack([2.0 ** np.arange(1, ngroups + 1),
                        2.0 ** np.arange(ngroups, 0, -1)])
    with h5py.File(path, "w") as h5:
        h5.attrs["sph_input_h5_sha256"] = file_sha256(source)
        dataset = h5.create_dataset("sph", data=factors)
        dataset.attrs["mixture_names"] = np.asarray(["HIGH", "LOW"], dtype="S")
        dataset.attrs["group_order"] = "mgxs_donjon"
    return factors


@pytest.mark.parametrize("states", [False, True])
@pytest.mark.parametrize("ngroups,order,layout", [
    (4, 3, "first"), (4, 3, "last"), (4, 3, "root_last"),
    (2, 1, "root_last"), (4, 3, "mixed_alias"),
    (3, 3, "inferred_first"), (3, 3, "inferred_last"), (3, 0, "2d"),
])
def test_sph_preserves_layout_and_scales_only_incoming_groups(tmp_path, states, ngroups, order, layout):
    source, sidecar, output = (tmp_path / name for name in ("in.h5", "sph.h5", "out.h5"))
    means, deviations, _ = _source(source, ngroups=ngroups, order=order, layout=layout)
    if states:
        _move_datasets_to_states(source)
    assert validate_input(source).ok
    factors = _bound_sidecar(sidecar, source, ngroups)
    source_before, sidecar_before = source.read_bytes(), sidecar.read_bytes()

    apply_sph_to_hdf5(source, sph_source=sidecar, output_h5=output)

    assert source.read_bytes() == source_before
    assert sidecar.read_bytes() == sidecar_before
    assert validate_input(output).ok
    mixtures, _ = read_mgxs_hdf5(output)
    with h5py.File(source) as src, h5py.File(output) as dst:
        for i, name in enumerate(("HIGH", "LOW")):
            expected = means[i] / factors[i][None, :, None]
            np.testing.assert_allclose(mixtures[i].scatter_matrix, expected)
            path = f"mixtures/{name}" + ("/states/00000001" if states else "")
            np.testing.assert_allclose(
                dst[f"{path}/total"][:],
                dst[f"{path}/absorption"][:] + mixtures[i].scatter_matrix[0].sum(axis=1),
            )
            for field, original in (("scatter_matrix", means), ("scatter_matrix_std_dev", deviations)):
                left, right = src[f"{path}/{field}"], dst[f"{path}/{field}"]
                assert left.shape == right.shape
                for key, value in left.attrs.items():
                    np.testing.assert_array_equal(right.attrs[key], value)
                # The fixture annotates each dataset independently of inherited
                # group axes, so the reference does not reuse the implementation.
                stored = right[:]
                canonical = stored[None] if stored.ndim == 2 else (
                    np.moveaxis(stored, -1, 0)
                    if left.attrs["axes"] == "from,to,moment" else stored
                )
                np.testing.assert_allclose(canonical, original[i] / factors[i][None, :, None])


@pytest.mark.parametrize("attribute", ["scatter_axes", "axes"])
def test_state_axes_override_mixture_and_root_declarations(tmp_path, attribute):
    source, sidecar, output = (tmp_path / name for name in ("in.h5", "sph.h5", "out.h5"))
    means, _, _ = _source(source, layout="last")
    _move_datasets_to_states(source)
    with h5py.File(source, "r+") as h5:
        h5.attrs["scatter_axes"] = "moment,from,to"
        for group in h5["mixtures"].values():
            group.attrs["scatter_axes"] = "moment,from,to"
            group["states/00000001"].attrs[attribute] = np.bytes_("G_in, G_out, legendre")
    factors = _bound_sidecar(sidecar, source, 4)
    assert validate_input(source).ok
    apply_sph_to_hdf5(source, sph_source=sidecar, output_h5=output)
    mixtures, _ = read_mgxs_hdf5(output)
    for i, mixture in enumerate(mixtures):
        np.testing.assert_allclose(mixture.scatter_matrix, means[i] / factors[i][None, :, None])


@pytest.mark.parametrize("existing", [False, True])
@pytest.mark.parametrize("invalid", ["ambiguous", "unsupported", "shape"])
def test_unresolved_scatter_layout_cannot_publish_output(tmp_path, existing, invalid):
    source, sidecar, output = (tmp_path / name for name in ("in.h5", "sph.h5", "out.h5"))
    _source(source)
    with h5py.File(source, "r+") as h5:
        group = h5["mixtures/LOW"]
        if invalid == "ambiguous":
            del group.attrs["scatter_axes"]
        elif invalid == "unsupported":
            group.attrs["scatter_axes"] = "from,moment,to"
        else:
            del group["scatter_matrix"]
            group.create_dataset("scatter_matrix", data=np.ones((4, 3, 4)))
    _bound_sidecar(sidecar, source, 4)
    before = source.read_bytes(), sidecar.read_bytes()
    if existing:
        output.write_bytes(b"previous output")
    message = "ambiguous scatter_matrix" if invalid == "ambiguous" else "scatter"
    with pytest.raises(ValueError, match=message):
        apply_sph_to_hdf5(source, sph_source=sidecar, output_h5=output, force=existing)
    assert (source.read_bytes(), sidecar.read_bytes()) == before
    if existing:
        assert output.read_bytes() == b"previous output"
    else:
        assert not output.exists()
    assert not list(tmp_path.glob(".openmc2donjon-sph-*"))


@pytest.mark.parametrize("binding", ["missing", "empty", "whitespace", "matching", "uppercase"])
def test_binding_status_is_truthful_in_hdf5_report_and_summary(tmp_path, capsys, binding):
    source, sidecar, output = (tmp_path / name for name in ("in.h5", "sph.h5", "out.h5"))
    _write_mgxs(source)
    _write_sidecar(sidecar)
    if binding != "missing":
        value = {"empty": "", "whitespace": "  ", "matching": file_sha256(source),
                 "uppercase": np.bytes_(file_sha256(source).upper())}[binding]
        with h5py.File(sidecar, "r+") as h5:
            h5.attrs["sph_input_h5_sha256"] = value
    before = source.read_bytes(), sidecar.read_bytes()
    report = apply_sph_to_hdf5(source, sph_source=sidecar, output_h5=output)
    verified = binding in {"matching", "uppercase"}
    mode = "converter-final-exact-input" if verified else "converter-unbound"
    assert report.sidecar_input_hash_verified is verified
    assert report.binding_mode == mode
    summary = summary_payload(report)
    assert summary["sidecar_input_hash_verified"] is verified
    assert summary["binding_mode"] == mode
    with h5py.File(output) as h5:
        assert bool(h5.attrs["sph_apply_sidecar_input_hash_verified"]) is verified
        assert h5.attrs["sph_apply_binding_mode"] == mode
        np.testing.assert_allclose(h5["mixtures/fuel/total"][:], [5.0, 40.0])
    print_report(report)
    text = capsys.readouterr().out
    assert f"sidecar input hash verified: {str(verified).lower()}" in text
    if not verified:
        assert "not a verified physical-SPH handoff" in text
    assert (source.read_bytes(), sidecar.read_bytes()) == before


@pytest.mark.parametrize("existing", [False, True])
@pytest.mark.parametrize("digest", ["not-a-digest", "a" * 64])
def test_invalid_or_mismatched_hash_is_not_downgraded_to_unbound(tmp_path, existing, digest):
    source, sidecar, output = (tmp_path / name for name in ("in.h5", "sph.h5", "out.h5"))
    _write_mgxs(source)
    _write_sidecar(sidecar)
    with h5py.File(sidecar, "r+") as h5:
        h5.attrs["sph_input_h5_sha256"] = digest
    if existing:
        output.write_bytes(b"previous output")
    with pytest.raises(ValueError, match="SHA-256"):
        apply_sph_to_hdf5(source, sph_source=sidecar, output_h5=output, force=existing)
    if existing:
        assert output.read_bytes() == b"previous output"
    else:
        assert not output.exists()


@pytest.mark.parametrize("remove_hash", [False, True])
def test_physical_gate_still_requires_verified_final_binding(tmp_path, remove_hash):
    source, sidecar, output, ce, mg = (tmp_path / name for name in (
        "in.h5", "sph.h5", "out.h5", "ce.h5", "mg.h5",
    ))
    _write_mgxs(source)
    _write_verified_flux(ce, "reference_flux")
    _write_verified_flux(mg, "mg_flux")
    create_openmc_sph_sidecar(
        source, sidecar, reference_flux=f"{ce}::reference_flux", mg_flux=f"{mg}::mg_flux",
        require_reference_flux_std_dev=True, max_reference_flux_std_dev_rel=.02,
        require_mg_flux_std_dev=True, max_mg_flux_std_dev_rel=.02,
    )
    if remove_hash:
        with h5py.File(sidecar, "r+") as h5:
            del h5.attrs["sph_input_h5_sha256"]
    report = apply_sph_to_hdf5(source, sph_source=sidecar, output_h5=output)
    issues = physical_sph_issues(output)
    if remove_hash:
        assert report.sidecar_input_hash_verified is False
        assert any("converter-final-exact-input" in issue for issue in issues)
        assert any("sidecar_input_hash_verified=true" in issue for issue in issues)
        assert any("sph_input_h5_sha256" in issue for issue in issues)
    else:
        assert report.sidecar_input_hash_verified is True
        assert issues == []
