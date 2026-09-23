from pathlib import Path

import h5py
import numpy as np
import pytest

from openmc2donjon.export_openmc_mgxs import DomainExportSpec, export_openmc_mgxs_library
from openmc2donjon.mgxs_input_contract import validate_input
from openmc2donjon.multicompo import read_mgxs_hdf5
from openmc2donjon.scatter_order import validate_order_metadata
from tests.test_export_openmc_mgxs import KeywordOnlyFakeLibrary


def mixed_library():
    lib = KeywordOnlyFakeLibrary()
    a, b = np.zeros((3, 3, 4)), np.zeros((3, 3, 2))
    a[:, :, 0], b[:, :, 0] = np.diag([.45, .54, .63]), np.diag([.19, .28, .37])
    a[:, :, 1], b[:, :, 1] = np.eye(3) * .01, np.eye(3) * .01
    a[:, :, 2], a[:, :, 3] = np.eye(3) * .005, np.eye(3) * -.002
    lib._data[(101, "scatter matrix")], lib._data[(102, "scatter matrix")] = a, b
    lib._data[(101, "transport")] = np.array([.49, .59, .69])
    lib._data[(102, "transport")] = np.array([.19, .29, .39])
    # Provide complete optional fields for the component-collapse tests.
    lib._data[(102, "kappa-fission")] = np.zeros(3)
    lib._data[(102, "inverse-velocity")] = np.array([1.1e-8, 2.2e-7, 3.3e-6])
    return lib


def test_mixed_export_and_preflight(tmp_path):
    lib = mixed_library()
    path = tmp_path / "mixed.h5"
    export_openmc_mgxs_library(lib, path)
    report = validate_input(path)
    assert report.ok, report.issues
    assert report.source_legendre_orders == {"ASM_1": 3, "MOD": 1}
    mixtures, _ = read_mgxs_hdf5(path)
    np.testing.assert_array_equal(mixtures[0].scatter_matrix,
                                  np.moveaxis(lib._data[(101, "scatter matrix")], -1, 0))
    np.testing.assert_array_equal(mixtures[1].scatter_matrix[:2],
                                  np.moveaxis(lib._data[(102, "scatter matrix")], -1, 0))
    assert np.all(mixtures[1].scatter_matrix[2:] == 0)


@pytest.mark.parametrize("corruption", ["order", "policy", "missing", "tail", "std_tail"])
def test_mixed_contract_rejects_bad_metadata_or_padding(tmp_path, corruption):
    path = tmp_path / "mixed.h5"
    export_openmc_mgxs_library(mixed_library(), path)
    with h5py.File(path, "r+") as handle:
        group = handle["mixtures/MOD"]
        if corruption == "order":
            group.attrs["source_legendre_order"] = 4
        elif corruption == "policy":
            group.attrs["scatter_padding"] = "none"
        elif corruption == "missing":
            del group.attrs["source_legendre_order"]
        elif corruption == "tail":
            group["scatter_matrix"][3, 0, 0] = 1e-40
        else:
            values = np.zeros((4, 3, 3))
            values[2, 0, 0] = 1e-40
            group.create_dataset("scatter_matrix_std_dev", data=values)
    assert not validate_input(path).ok
    with pytest.raises(ValueError):
        read_mgxs_hdf5(path)


@pytest.mark.parametrize("order", [True, 1.5, "1", -1, 4])
def test_invalid_source_orders(order):
    with pytest.raises(ValueError):
        validate_order_metadata({"source_legendre_order": order, "scatter_padding": "zero-truncation"},
                                np.zeros((4, 2, 2)), axes="moment,from,to", stored_order=3)


def test_moment_last_and_legacy():
    attrs = {"source_legendre_order": 1, "scatter_padding": "zero-truncation"}
    assert validate_order_metadata(attrs, np.zeros((2, 2, 4)),
                                   axes="from,to,moment", stored_order=3) == 1
    assert validate_order_metadata({}, np.zeros((4, 2, 2)), axes=None, stored_order=3) is None


def test_order_attrs_are_exporter_owned(tmp_path):
    lib = mixed_library()
    with pytest.raises(ValueError, match="exporter-owned"):
        export_openmc_mgxs_library(lib, tmp_path / "bad.h5", domain_specs=[
            DomainExportSpec(d, attrs={"source_legendre_order": 3}) for d in lib.domains])


@pytest.mark.parametrize("bad", ["selected_moment", "missing_moment", "std_shape"])
def test_export_rejects_incomplete_moments(tmp_path, monkeypatch, bad):
    lib = mixed_library()
    original = lib.get_mgxs

    def get_mgxs(*, domain, mgxs_type):
        obj = original(domain=domain, mgxs_type=mgxs_type)
        if mgxs_type == "scatter matrix" and domain.id == 101:
            if bad == "missing_moment":
                obj.legendre_order = 4
            if bad == "std_shape":
                obj.std_dev = np.zeros((3, 3, 2))
        return obj

    monkeypatch.setattr(lib, "get_mgxs", get_mgxs)
    specs = [DomainExportSpec(d, xs_kwargs={"moment": 0}) for d in lib.domains]
    output = tmp_path / "bad.h5"
    with pytest.raises(ValueError, match="moment|Moment"):
        export_openmc_mgxs_library(lib, output,
                                  domain_specs=specs if bad == "selected_moment" else None)
    assert not output.exists()


@pytest.mark.parametrize("name", [
    "source_legendre_order", "scatter_padding", "legendre_order", "energy_groups",
])
def test_root_cannot_override_computed_orders(tmp_path, name):
    with pytest.raises(ValueError, match="exporter-owned"):
        export_openmc_mgxs_library(mixed_library(), tmp_path / "bad.h5", root_attrs={name: 1})


def test_collapse_keeps_high_order_when_first_member_is_p1(tmp_path):
    from openmc2donjon.mixture_collapse import collapse_components
    source, output = tmp_path / "mixed.h5", tmp_path / "collapsed.h5"
    export_openmc_mgxs_library(mixed_library(), source)
    with h5py.File(source, "r+") as handle:
        flux = handle.create_dataset("openmc_volume_flux", data=np.ones((2, 3)))
        flux.attrs["group_order"] = "mgxs_donjon"
        flux.attrs["mixture_names"] = handle["mixture_names"][:]
    collapse_components(source, output, groups=[("COMBINED", ("MOD", "ASM_1"))])
    with h5py.File(output) as handle:
        group = handle["mixtures/COMBINED"]
        assert group.attrs["source_legendre_order"] == 3
        assert group.attrs["scatter_padding"] == "none"
        np.testing.assert_array_equal(group.attrs["collapsed_source_legendre_orders"], [1, 3])
        assert np.any(group["scatter_matrix"][3] != 0)
    read_mgxs_hdf5(output)


def test_collapse_rejects_corrupt_truncation_before_writing(tmp_path):
    from openmc2donjon.mixture_collapse import collapse_components
    source, output = tmp_path / "mixed.h5", tmp_path / "collapsed.h5"
    export_openmc_mgxs_library(mixed_library(), source)
    with h5py.File(source, "r+") as handle:
        flux = handle.create_dataset("openmc_volume_flux", data=np.ones((2, 3)))
        flux.attrs["group_order"] = "mgxs_donjon"
        flux.attrs["mixture_names"] = handle["mixture_names"][:]
        handle["mixtures/MOD/scatter_matrix"][2, 0, 0] = 0.01
    with pytest.raises(ValueError, match="truncated moments"):
        collapse_components(source, output, groups=[("COMBINED", ("MOD", "ASM_1"))])
    assert not output.exists()


@pytest.mark.parametrize("legacy_names", [("ASM_1",), ("MOD",), ("MOD", "ASM_1")])
def test_collapse_preserves_unknown_orders_including_repeat(tmp_path, legacy_names, capsys):
    from openmc2donjon.mixture_collapse import collapse_components
    from openmc2donjon.mgxs_input_report import input_report_payload, print_report
    source, output, repeat = (tmp_path / name for name in ("source.h5", "collapsed.h5", "repeat.h5"))
    export_openmc_mgxs_library(mixed_library(), source)
    with h5py.File(source, "r+") as handle:
        for name in legacy_names:
            group = handle[f"mixtures/{name}"]
            del group.attrs["source_legendre_order"]
            del group.attrs["scatter_padding"]
            group["scatter_matrix"][2:] = 0
        flux = handle.create_dataset("openmc_volume_flux", data=np.ones((2, 3)))
        flux.attrs["group_order"] = "mgxs_donjon"
        flux.attrs["mixture_names"] = handle["mixture_names"][:]
    assert validate_input(source).ok
    collapse_components(source, output, groups=[("COMBINED", ("MOD", "ASM_1"))])
    with h5py.File(output) as handle:
        group = handle["mixtures/COMBINED"]
        assert "source_legendre_order" not in group.attrs
        assert "scatter_padding" not in group.attrs
        np.testing.assert_array_equal(group.attrs["collapsed_source_legendre_orders"],
                                      [-1 if n in legacy_names else order
                                       for n, order in (("MOD", 1), ("ASM_1", 3))])
    report = validate_input(output)
    assert report.ok, report.issues
    assert report.source_legendre_orders == {}
    assert input_report_payload(report)["source_legendre_orders_unknown"] == ["COMBINED"]
    print_report(report)
    assert "source scattering orders: unknown: 1" in capsys.readouterr().out
    read_mgxs_hdf5(output)
    collapse_components(output, repeat, groups=[("AGAIN", ("COMBINED",))])
    with h5py.File(repeat) as handle:
        group = handle["mixtures/AGAIN"]
        assert "source_legendre_order" not in group.attrs
        np.testing.assert_array_equal(group.attrs["collapsed_source_legendre_orders"], [-1])
    assert validate_input(repeat).ok


def mixed_fill_case(tmp_path, monkeypatch):
    from tests.test_zero_flux_fill import _FakeXSData, _FakeMGXSLibrary, _touch_macrolib
    source = tmp_path / "mixed.h5"
    export_openmc_mgxs_library(mixed_library(), source)
    donors = []
    with h5py.File(source, "r+") as handle:
        for name, group in handle["mixtures"].items():
            group.attrs["fill_label"] = name
            scatter = group["scatter_matrix"][:]
            group.create_dataset("scatter_matrix_std_dev", data=np.abs(scatter) * .01)
            scatter[2], scatter[3] = np.eye(3) * .005, np.eye(3) * -.002
            donors.append(_FakeXSData(name, total=group["total"][:][::-1],
                absorption=group["absorption"][:][::-1],
                scatter=np.moveaxis(scatter[:, ::-1, ::-1], 0, -1),
                fission=group["fission"][:][::-1], nu_fission=group["nu_fission"][:][::-1],
                fissionable=bool(group.attrs["fissionable"])))
            group["total"][1], group["transport_total"][1] = 0, 0
    monkeypatch.setattr("openmc2donjon.zero_flux_fill._load_macrolib", lambda _: _FakeMGXSLibrary(donors))
    return source, _touch_macrolib(tmp_path)


@pytest.mark.parametrize("in_place", [False, True])
def test_fill_preserves_local_p1_truncation_and_p3_values(tmp_path, monkeypatch, in_place):
    from openmc2donjon.zero_flux_fill import fill_zero_flux_groups
    source, donor = mixed_fill_case(tmp_path, monkeypatch)
    output = source if in_place else tmp_path / "filled.h5"
    source_bytes = source.read_bytes()
    report = fill_zero_flux_groups(source, macrolib=donor, label_attr="fill_label",
                                  in_place=in_place, output_h5=None if in_place else output)
    assert report.total_filled_bins == 2
    if not in_place:
        assert source.read_bytes() == source_bytes
    with h5py.File(output) as handle:
        low, high = handle["mixtures/MOD"], handle["mixtures/ASM_1"]
        assert low.attrs["source_legendre_order"] == 1
        assert low.attrs["scatter_padding"] == "zero-truncation"
        assert np.all(low["scatter_matrix"][2:] == 0)
        assert np.all(low["scatter_matrix_std_dev"][2:] == 0)
        assert low.attrs["zero_flux_scatter_order"] == 3
        assert low.attrs["zero_flux_applied_scatter_order"] == 1
        assert high.attrs["zero_flux_applied_scatter_order"] == 3
        np.testing.assert_array_equal(high["scatter_matrix"][2:, 1, 1], [.005, -.002])
    check = validate_input(output)
    assert check.ok, check.issues
    read_mgxs_hdf5(output)


@pytest.mark.parametrize("in_place", [False, True])
@pytest.mark.parametrize("corrupt", ["tail", "std_tail", "missing", "bad_order"])
def test_fill_rejects_bad_local_order_without_changing_files(tmp_path, monkeypatch, in_place, corrupt):
    from openmc2donjon.zero_flux_fill import fill_zero_flux_groups
    source, donor = mixed_fill_case(tmp_path, monkeypatch)
    with h5py.File(source, "r+") as handle:
        group = handle["mixtures/MOD"]
        if corrupt == "tail":
            group["scatter_matrix"][2, 1, 1] = .01
        elif corrupt == "std_tail":
            group["scatter_matrix_std_dev"][3, 0, 0] = .01
        elif corrupt == "missing":
            del group.attrs["scatter_padding"]
        else:
            group.attrs["source_legendre_order"] = 4
    output = tmp_path / "existing.h5"
    export_openmc_mgxs_library(mixed_library(), output)
    source_bytes, output_bytes = source.read_bytes(), output.read_bytes()
    with pytest.raises(ValueError):
        fill_zero_flux_groups(source, macrolib=donor, label_attr="fill_label",
                              in_place=in_place, output_h5=None if in_place else output, force=True)
    assert source.read_bytes() == source_bytes
    assert output.read_bytes() == output_bytes


def test_noop_fill_still_rejects_invalid_order_metadata(tmp_path, monkeypatch):
    from openmc2donjon.zero_flux_fill import fill_zero_flux_groups
    source, donor = mixed_fill_case(tmp_path, monkeypatch)
    with h5py.File(source, "r+") as handle:
        for group in handle["mixtures"].values():
            group["total"][1], group["transport_total"][1] = 1.0, .9
        handle["mixtures/MOD"].attrs["source_legendre_order"] = 4
    before = source.read_bytes()
    with pytest.raises(ValueError, match="source_legendre_order"):
        fill_zero_flux_groups(source, macrolib=donor, label_attr="fill_label", in_place=True)
    assert source.read_bytes() == before


def test_fill_without_transport_retains_order_audit(tmp_path, monkeypatch):
    from openmc2donjon.zero_flux_fill import fill_zero_flux_groups
    source, donor = mixed_fill_case(tmp_path, monkeypatch)
    with h5py.File(source, "r+") as handle:
        del handle["mixtures/MOD/transport_total"]
        del handle["mixtures/MOD"].attrs["openmc_transport_mgxs_type"]
        del handle.attrs["openmc_transport_mgxs_type"]
    fill_zero_flux_groups(source, macrolib=donor, label_attr="fill_label", in_place=True)
    with h5py.File(source) as handle:
        group = handle["mixtures/MOD"]
        assert group.attrs["zero_flux_scatter_order"] == 3
        assert group.attrs["zero_flux_applied_scatter_order"] == 1
        assert np.all(group["scatter_matrix"][2:] == 0)


def test_fill_does_not_invent_original_order_for_legacy_region(tmp_path, monkeypatch):
    from openmc2donjon.zero_flux_fill import fill_zero_flux_groups
    source, donor = mixed_fill_case(tmp_path, monkeypatch)
    with h5py.File(source, "r+") as handle:
        group = handle["mixtures/MOD"]
        del group.attrs["source_legendre_order"]
        del group.attrs["scatter_padding"]
    fill_zero_flux_groups(source, macrolib=donor, label_attr="fill_label", in_place=True)
    with h5py.File(source) as handle:
        group = handle["mixtures/MOD"]
        assert "source_legendre_order" not in group.attrs
        assert "scatter_padding" not in group.attrs
        assert group.attrs["zero_flux_applied_scatter_order"] == 3
    report = validate_input(source)
    assert report.ok, report.issues
    assert report.source_legendre_orders_unknown == ["MOD"]


def native_library():
    openmc = pytest.importorskip("openmc")
    pytest.importorskip("openmc.mgxs")
    a, b = openmc.Cell(), openmc.Cell()
    lib = openmc.mgxs.Library(openmc.Geometry([a, b]))
    lib.domain_type, lib.domains = "cell", [a, b]
    lib.energy_groups = openmc.mgxs.EnergyGroups([1e-5, 1e3, 2e7])
    lib.mgxs_types = ["total", "absorption", "nu-transport", "consistent nu-scatter matrix"]
    lib.correction, lib.legendre_order = None, 3
    lib.build_library()
    return openmc, lib


def test_native_grouped_tallies_keep_high_moments_off_p1_cells(tmp_path):
    from openmc2donjon.openmc_mixed_order import exact_cell_mgxs_tallies, set_domain_scatter_orders
    openmc, lib = native_library()
    a, b = lib.domains
    set_domain_scatter_orders(lib, {a.id: 1, b.id: 3})
    assert lib.legendre_order == 3
    tallies = exact_cell_mgxs_tallies(lib)
    path = tmp_path / "tallies.xml"
    tallies.export_to_xml(path)
    # Intentional reload of the same IDs in this process, not ID collisions
    # in the exported tally definition.
    with pytest.warns(openmc.IDWarning):
        reread = openmc.Tallies.from_xml(Path(path))
    by_order = {}
    for tally in reread:
        for f in tally.filters:
            if isinstance(f, openmc.LegendreFilter) and f.order > 1:
                assert list(tally.find_filter(openmc.CellFilter).bins) == [b.id]
            if isinstance(f, openmc.LegendreFilter):
                by_order.setdefault(f.order, set()).update(tally.find_filter(openmc.CellFilter).bins)
    assert by_order[3] == {b.id}
    assert a.id in by_order[1]
    assert lib.get_mgxs(a, "consistent nu-scatter matrix").legendre_order == 1


@pytest.mark.parametrize("bad", ["missing", "extra", "bool", "negative", "cached"])
def test_native_order_configuration_fails_closed(bad):
    from openmc2donjon.openmc_mixed_order import set_domain_scatter_orders
    _, lib = native_library()
    a, b = lib.domains
    orders = {a.id: 1, b.id: 3}
    if bad == "missing":
        orders.pop(b.id)
    elif bad == "extra":
        orders[999999] = 1
    elif bad == "bool":
        orders[a.id] = True
    elif bad == "negative":
        orders[a.id] = -1
    else:
        _ = lib.get_mgxs(a, "total").tallies
    with pytest.raises(ValueError):
        set_domain_scatter_orders(lib, orders)
