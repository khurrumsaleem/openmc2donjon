"""Layout-invariance regressions, including cubic 4-group/P3 matrices."""

import h5py
import numpy as np
import pytest

from openmc2donjon.energy_collapse import collapse_energy_groups
from openmc2donjon.mgxs_input_contract import validate_input
from openmc2donjon.mixture_collapse import collapse_components
from openmc2donjon.multicompo import convert_mgxs_hdf5, read_mgxs_hdf5


def _source(path, ngroups=4, order=3, layout="first"):
    """Use off-diagonal transfers and unequal fluxes to expose wrong-axis weights."""
    flux = np.stack([2.0 ** np.arange(ngroups), 2.0 ** np.arange(ngroups)[::-1]])
    means = np.zeros((2, order + 1, ngroups, ngroups))
    p0 = np.roll(np.eye(ngroups), 1, axis=1) * 0.4
    means[:, 0] = p0
    if order >= 1:
        means[0, 1], means[1, 1] = p0 * 0.1, p0 * 0.05
    if order == 3:
        means[0, 2], means[0, 3] = p0 * 0.05, p0 * -0.005
    std = np.abs(means) * 0.03
    with h5py.File(path, "w") as h5:
        h5.attrs["energy_groups"] = ngroups
        h5.attrs["legendre_order"] = order
        h5.create_dataset("energy_bounds", data=np.geomspace(1e-5, 1e7, ngroups + 1))
        h5.create_dataset("mixture_names", data=np.asarray(["HIGH", "LOW"], dtype="S"))
        weights = h5.create_dataset("openmc_volume_flux", data=flux)
        weights.attrs["group_order"] = "mgxs_donjon"
        weights.attrs["mixture_names"] = h5["mixture_names"][:]
        if layout in {"root_last", "mixed_alias"}:
            h5.attrs["scatter_axes"] = "from,to,moment"
            h5.attrs["axes"] = "from,to,moment"
        for i, name in enumerate(("HIGH", "LOW")):
            group = h5.create_group(f"mixtures/{name}")
            group.attrs["volume"] = 1.0
            group.attrs["fissionable"] = False
            local_order = min(order, 1) if name == "LOW" else order
            group.attrs["source_legendre_order"] = local_order
            group.attrs["scatter_padding"] = "none" if local_order == order else "zero-truncation"
            for field, values in (
                ("total", np.ones(ngroups)),
                ("absorption", np.full(ngroups, 0.6)),
                ("transport_total", np.ones(ngroups) - means[i, 1].sum(axis=1)
                 if order else np.ones(ngroups)),
                ("fission", np.zeros(ngroups)),
                ("nu_fission", np.zeros(ngroups)),
                ("chi", np.zeros(ngroups)),
            ):
                group.create_dataset(field, data=values)
            last = layout in {"last", "root_last", "inferred_last"} or (
                layout in {"mixed", "mixed_alias"} and name == "LOW"
            )
            axes = "from,to,moment" if last else "moment,from,to"
            if layout == "mixed_alias" and name == "HIGH":
                group.attrs["axes"] = np.bytes_("moment,G_in,G_out")
            elif layout not in {"root_last", "mixed_alias", "inferred_first", "inferred_last", "2d"}:
                group.attrs["scatter_axes"] = axes
            for field, array in (("scatter_matrix", means[i]), ("scatter_matrix_std_dev", std[i])):
                stored = array[0] if layout == "2d" else np.moveaxis(array, 0, -1) if last else array
                dataset = group.create_dataset(field, data=stored)
                # Old dataset annotations must not survive with a stale layout.
                dataset.attrs["axes"] = axes
                dataset.attrs["scatter_axes"] = axes
    return means, std, flux


def _collapse(operation, source, output, ngroups, *, force=False):
    if operation == "component":
        return collapse_components(source, output, groups=[("BOTH", ("HIGH", "LOW"))], force=force)
    split = max(1, ngroups // 2)
    mapping = (tuple(range(1, split + 1)),)
    if split < ngroups:
        mapping += (tuple(range(split + 1, ngroups + 1)),)
    return collapse_energy_groups(source, output, groups=mapping, force=force)


def _expected(operation, values, flux):
    if operation == "component":
        return [np.sum(values * flux[:, None, :, None], axis=0) / flux.sum(axis=0)[None, :, None]]
    ngroups = flux.shape[1]
    split = max(1, ngroups // 2)
    mapping = [list(range(split))]
    if split < ngroups:
        mapping.append(list(range(split, ngroups)))
    result = np.zeros((2, values.shape[1], len(mapping), len(mapping)))
    # An independent scalar reference checks every moment's transfer rate.
    for region in range(2):
        for moment in range(values.shape[1]):
            for a, incoming in enumerate(mapping):
                for b, outgoing in enumerate(mapping):
                    result[region, moment, a, b] = sum(
                        values[region, moment, g, h] * flux[region, g]
                        for g in incoming for h in outgoing
                    ) / sum(flux[region, g] for g in incoming)
    return result


@pytest.mark.parametrize("operation", ["component", "energy"])
@pytest.mark.parametrize("ngroups,order", [(4, 3), (3, 3), (2, 1), (1, 0)])
@pytest.mark.parametrize("layout", ["first", "last", "root_last", "mixed", "mixed_alias"])
def test_collapse_is_layout_invariant(tmp_path, operation, ngroups, order, layout):
    _check_layout(tmp_path, operation, ngroups, order, layout)


@pytest.mark.parametrize("operation", ["component", "energy"])
@pytest.mark.parametrize("layout", ["inferred_first", "inferred_last", "2d"])
def test_inferred_layout_and_2d_p0(tmp_path, operation, layout):
    _check_layout(tmp_path, operation, 3, 0 if layout == "2d" else 3, layout)


def _check_layout(tmp_path, operation, ngroups, order, layout):
    source, output = tmp_path / "source.h5", tmp_path / "collapsed.h5"
    mean, std, flux = _source(source, ngroups, order, layout)
    source_bytes = source.read_bytes()
    report = validate_input(source)
    assert report.ok, report.issues
    _collapse(operation, source, output, ngroups)
    report = validate_input(output)
    assert report.ok, report.issues
    expected_mean = _expected(operation, mean, flux)
    expected_std = _expected(operation, std, flux)
    mixtures, _ = read_mgxs_hdf5(output)
    names = ("BOTH",) if operation == "component" else ("HIGH", "LOW")
    with h5py.File(output) as h5:
        assert h5.attrs["scatter_axes"] == "moment,from,to"
        for i, name in enumerate(names):
            group = h5[f"mixtures/{name}"]
            assert group.attrs["source_legendre_order"] == (min(order, 1) if name == "LOW" else order)
            expected_padding = "zero-truncation" if name == "LOW" and order == 3 else "none"
            assert group.attrs["scatter_padding"] == expected_padding
            for attrs in (h5.attrs, group.attrs, group["scatter_matrix"].attrs,
                          group["scatter_matrix_std_dev"].attrs):
                assert attrs["scatter_axes"] == "moment,from,to"
                if "axes" in attrs:
                    assert attrs["axes"] == "moment,from,to"
            np.testing.assert_allclose(group["scatter_matrix"][:], expected_mean[i], rtol=1e-14, atol=1e-16)
            np.testing.assert_allclose(
                group["scatter_matrix_std_dev"][:], expected_std[i], rtol=1e-14, atol=1e-16,
            )
            np.testing.assert_allclose(mixtures[i].scatter_matrix, expected_mean[i], rtol=1e-14, atol=1e-16)
            if operation == "energy":
                assert group["scatter_matrix_std_dev"].attrs["energy_uncertainty_method"] == (
                    "conservative-l1-source-xs-bound-no-covariance"
                )
    convert_mgxs_hdf5(output, tmp_path / "collapsed.mcompo.txt")
    assert source.read_bytes() == source_bytes


@pytest.mark.parametrize("operation", ["component", "energy"])
@pytest.mark.parametrize("existing_output", [False, True])
@pytest.mark.parametrize("corruption", ["ambiguous", "bad_axes", "shape", "std_shape", "std_tail",
                                        "negative_std", "nonfinite_std", "tail"])
def test_invalid_scatter_rejected_before_writing(tmp_path, operation, existing_output, corruption):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    _source(source, layout="last")
    with h5py.File(source, "r+") as h5:
        group = h5["mixtures/LOW"]
        if corruption == "ambiguous":
            del group.attrs["scatter_axes"]
        elif corruption == "bad_axes":
            group.attrs["scatter_axes"] = "from,moment,to"
        elif corruption in {"shape", "std_shape"}:
            name = "scatter_matrix" if corruption == "shape" else "scatter_matrix_std_dev"
            del group[name]
            group.create_dataset(name, data=np.zeros((3, 4, 4)))
        elif corruption == "std_tail":
            group["scatter_matrix_std_dev"][0, 1, 2] = 0.001
        elif corruption == "tail":
            group["scatter_matrix"][0, 1, 2] = 0.001
        else:
            group["scatter_matrix_std_dev"][0, 1, 0] = -0.001 if corruption == "negative_std" else np.nan
    before = source.read_bytes()
    if existing_output:
        output.write_bytes(b"existing destination must survive")
    with pytest.raises(ValueError):
        _collapse(operation, source, output, 4, force=existing_output)
    assert source.read_bytes() == before
    if existing_output:
        assert output.read_bytes() == b"existing destination must survive"
    else:
        assert not output.exists()


@pytest.mark.parametrize("operation", ["component", "energy"])
def test_absent_uncertainty_is_not_invented(tmp_path, operation):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    _source(source, layout="last")
    with h5py.File(source, "r+") as h5:
        del h5["mixtures/LOW/scatter_matrix_std_dev"]
    _collapse(operation, source, output, 4)
    with h5py.File(output) as h5:
        name = "BOTH" if operation == "component" else "LOW"
        assert "scatter_matrix_std_dev" not in h5[f"mixtures/{name}"]


def test_component_and_energy_collapse_can_be_chained(tmp_path):
    source, spatial, energy = (tmp_path / name for name in ("source.h5", "spatial.h5", "energy.h5"))
    _source(source, layout="mixed_alias")
    collapse_components(source, spatial, groups=[("BOTH", ("HIGH", "LOW"))])
    collapse_energy_groups(spatial, energy, groups=((1, 2), (3, 4)))
    assert validate_input(energy).ok
    mixtures, _ = read_mgxs_hdf5(energy)
    assert mixtures[0].scatter_matrix.shape == (4, 2, 2)
