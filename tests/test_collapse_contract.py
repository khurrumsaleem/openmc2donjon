"""Collapse must preserve field meaning, ordering, and donor provenance."""

import json

import h5py
import numpy as np
import pytest

from openmc2donjon.collapse_substitution import RECORD_DATASET, SEMANTICS, read_substitution_records
from openmc2donjon.energy_collapse import collapse_energy_groups
from openmc2donjon.mgxs_fields import H_FACTOR_DATASETS, INVERSE_VELOCITY_DATASETS
from openmc2donjon.mgxs_input_contract import validate_input
from openmc2donjon.mixture_collapse import collapse_components
from openmc2donjon.multicompo import convert_mgxs_hdf5, read_mgxs_hdf5
from openmc2donjon.openmc_provenance import read_openmc_provenance
from openmc2donjon.zero_flux_fill import fill_zero_flux_groups
from tests.test_collapse_safety import _bind, _run
from tests.test_collapse_scatter_layout import _source
from tests.test_zero_flux_fill import _fake_library, _fake_openmc, _touch_macrolib


@pytest.mark.parametrize("operation", ["component", "energy"])
@pytest.mark.parametrize("dataset", ["openmc_volume_flux", "openmc_volume_flux_std_dev"])
@pytest.mark.parametrize("fault", ["missing_order", "reversed_order", "missing_names", "reversed_names"])
def test_invalid_flux_labels_fail_without_repair(tmp_path, operation, dataset, fault):
    source, output, summary = (tmp_path / name for name in ("in.h5", "out.h5", "summary.json"))
    _source(source)
    with h5py.File(source, "r+") as h5:
        std = h5.create_dataset("openmc_volume_flux_std_dev", data=h5["openmc_volume_flux"][:] * .01)
        for key, value in h5["openmc_volume_flux"].attrs.items():
            std.attrs[key] = value
        field = "group_order" if fault.endswith("order") else "mixture_names"
        if fault.startswith("missing"):
            del h5[dataset].attrs[field]
        else:
            h5[dataset].attrs[field] = (
                "openmc" if field == "group_order" else np.array(["LOW", "HIGH"], dtype="S")
            )
            if field == "mixture_names":
                h5[dataset][:] = h5[dataset][:][::-1]
    _assert_rejected_unchanged(operation, source, output, summary, "reference flux contract")


def _assert_rejected_unchanged(operation, source, output, summary, match):
    before = source.read_bytes()
    output.write_bytes(b"old output")
    summary.write_bytes(b"old summary")
    with pytest.raises(ValueError, match=match):
        _run(operation, source, output, force=True, summary_json=summary)
    assert source.read_bytes() == before
    assert output.read_bytes() == b"old output"
    assert summary.read_bytes() == b"old summary"


ALIASES = [(name, "kappa_fission", "h_factor") for name in H_FACTOR_DATASETS] + [
    (name, "inverse_velocity", "inverse_velocity") for name in INVERSE_VELOCITY_DATASETS
]


@pytest.mark.parametrize("operation", ["component", "energy"])
@pytest.mark.parametrize("alias,canonical,converter_field", ALIASES)
def test_all_converter_aliases_keep_values_and_uncertainty(
    tmp_path, operation, alias, canonical, converter_field,
):
    source, output = tmp_path / "in.h5", tmp_path / "out.h5"
    _, _, flux = _source(source)
    values = np.array([[1., 2., 3., 4.], [8., 4., 2., 1.]])
    with h5py.File(source, "r+") as h5:
        for i, name in enumerate(("HIGH", "LOW")):
            group = h5[f"mixtures/{name}"]
            # Different valid spellings in the same output component.
            field = alias if i == 0 else canonical
            group.create_dataset(field, data=values[i])
            group.create_dataset(field + "_std_dev", data=values[i] * .02)
    _run(operation, source, output)
    assert validate_input(output).ok
    mixtures, _ = read_mgxs_hdf5(output)
    if operation == "component":
        expected = [(values * flux).sum(axis=0) / flux.sum(axis=0)]
    else:
        expected = [[np.dot(v[g:g+2], f[g:g+2]) / f[g:g+2].sum() for g in (0, 2)]
                    for v, f in zip(values, flux, strict=True)]
    with h5py.File(output) as h5:
        for i, mix in enumerate(mixtures):
            group = h5[f"mixtures/{mix.name}"]
            np.testing.assert_allclose(group[canonical][:], expected[i])
            np.testing.assert_allclose(group[canonical + "_std_dev"][:], np.asarray(expected[i]) * .02)
            np.testing.assert_allclose(getattr(mix, converter_field), expected[i])
            if alias != canonical:
                assert alias not in group
    convert_mgxs_hdf5(output, tmp_path / "converted.mcompo.txt")


@pytest.mark.parametrize("operation", ["component", "energy"])
@pytest.mark.parametrize("conflict", ["mean", "std", "none"])
def test_duplicate_aliases_cannot_hide_conflicts_or_uncertainty(tmp_path, operation, conflict):
    source, output, summary = (tmp_path / name for name in ("in.h5", "out.h5", "summary.json"))
    _source(source)
    with h5py.File(source, "r+") as h5:
        for group in h5["mixtures"].values():
            group.create_dataset("h_factor", data=np.ones(4))
            group.create_dataset("kappa_fission", data=np.ones(4) * (2 if conflict == "mean" else 1))
            group.create_dataset("kappa_fission_std_dev", data=np.ones(4) * .02)
            if conflict == "std":
                group.create_dataset("h_factor_std_dev", data=np.ones(4) * .03)
    if conflict != "none":
        _assert_rejected_unchanged(operation, source, output, summary, "conflicting alias")
    else:
        _run(operation, source, output)
        with h5py.File(output) as h5:
            for group in h5["mixtures"].values():
                np.testing.assert_allclose(group["kappa_fission_std_dev"][:], .02)


@pytest.mark.parametrize("canonical", ["kappa_fission", "inverse_velocity"])
@pytest.mark.parametrize("partial", ["mean", "std"])
def test_partial_optional_coverage_is_not_silently_dropped(tmp_path, canonical, partial):
    source, output, summary = (tmp_path / name for name in ("in.h5", "out.h5", "summary.json"))
    _source(source)
    with h5py.File(source, "r+") as h5:
        h5["mixtures/HIGH"].create_dataset(canonical, data=np.ones(4))
        if partial == "std":
            h5["mixtures/LOW"].create_dataset(canonical, data=np.ones(4))
            h5["mixtures/HIGH"].create_dataset(canonical + "_std_dev", data=np.ones(4) * .02)
    _assert_rejected_unchanged("component", source, output, summary, "partial " + canonical)


def _mark_filled(group, groups, donor):
    group.attrs["zero_flux_filled_groups"] = groups
    group.attrs["zero_flux_fill_source"] = donor
    group.attrs["zero_flux_scatter_order"] = 3
    group.attrs["zero_flux_applied_scatter_order"] = int(group.attrs["source_legendre_order"])
    group.attrs["zero_flux_transport_method"] = "macrolib_p1_outscatter"


@pytest.mark.parametrize("first", ["component", "energy"])
@pytest.mark.parametrize("first_member_filled", [False, True])
def test_substitutions_follow_repeated_mappings_and_remain_bound(tmp_path, first, first_member_filled):
    source, stage, output, final = (tmp_path / name for name in ("in.h5", "stage.h5", "out.h5", "final.h5"))
    _source(source)
    with h5py.File(source, "r+") as h5:
        _mark_filled(h5["mixtures/LOW"], [3], "low-donor.h5")
        if first_member_filled:
            _mark_filled(h5["mixtures/HIGH"], [0, 2], "high-donor.h5")
    _bind(source, tmp_path)
    before = source.read_bytes()
    _run(first, source, stage)
    if first == "component":
        collapse_energy_groups(stage, output, groups=((1, 2), (3, 4)))
    else:
        collapse_components(stage, output, groups=[("BOTH", ("HIGH", "LOW"))])
    with h5py.File(output) as h5:
        group = h5["mixtures/BOTH"]
        assert group.attrs["zero_flux_fill_semantics"] == SEMANTICS
        assert "zero_flux_fill_source" not in group.attrs
        np.testing.assert_array_equal(group.attrs["zero_flux_filled_groups"],
                                      [0, 1] if first_member_filled else [1])
        origins = read_substitution_records(group, 2)
        assert len(origins) == (2 if first_member_filled else 1)
        low = next(item for item in origins if item["source_mixture"].endswith("/LOW"))
        assert low["filled_groups"] == [3]
        assert low["current_groups"] == [1]
        assert low["source_h5"] == str(source.resolve())
        assert low["donor_metadata"]["zero_flux_fill_source"] == "low-donor.h5"
        assert low["donor_metadata"]["zero_flux_scatter_order"] == 3
        assert low["donor_metadata"]["zero_flux_applied_scatter_order"] == 1
    collapse_energy_groups(output, final, groups=((1, 2),))
    assert validate_input(final).ok
    assert read_openmc_provenance(final)["integrity"]["ok"]
    with h5py.File(final, "r+") as h5:
        origins = read_substitution_records(h5["mixtures/BOTH"], 1)
        assert all(item["current_groups"] == [0] for item in origins)
        assert all(item["source_h5"] == str(source.resolve()) for item in origins)
        record = json.loads(h5[f"mixtures/BOTH/{RECORD_DATASET}"].asstr()[()])
        record["origins"][0]["donor_metadata"]["zero_flux_fill_source"] = "tampered.h5"
        h5[f"mixtures/BOTH/{RECORD_DATASET}"][()] = json.dumps(record)
    assert not read_openmc_provenance(final)["integrity"]["ok"]
    assert source.read_bytes() == before


@pytest.mark.parametrize("operation", ["component", "energy"])
@pytest.mark.parametrize("bad", [[4], [-1], [1.5], [True], [1, 1], [[1]]])
def test_invalid_substitution_indices_rejected(tmp_path, operation, bad):
    source, output, summary = (tmp_path / name for name in ("in.h5", "out.h5", "summary.json"))
    _source(source)
    with h5py.File(source, "r+") as h5:
        _mark_filled(h5["mixtures/LOW"], bad, "donor.h5")
    _assert_rejected_unchanged(operation, source, output, summary, "zero-flux group indices")


@pytest.mark.parametrize("fault", ["missing_records", "summary", "current", "origin", "donor", "schema"])
def test_invalid_derived_records_cannot_be_carried_forward(tmp_path, fault):
    source, stage, output = (tmp_path / name for name in ("in.h5", "stage.h5", "out.h5"))
    _source(source)
    with h5py.File(source, "r+") as h5:
        _mark_filled(h5["mixtures/LOW"], [3], "donor.h5")
    _run("energy", source, stage)
    with h5py.File(stage, "r+") as h5:
        group = h5["mixtures/LOW"]
        if fault == "missing_records":
            del group[RECORD_DATASET]
        elif fault == "summary":
            group.attrs["zero_flux_filled_groups"] = [0]
        else:
            record = json.loads(group[RECORD_DATASET].asstr()[()])
            if fault == "current":
                record["origins"][0]["current_groups"] = [3]
            elif fault == "origin":
                record["origins"][0]["filled_groups"] = [4]
            elif fault == "donor":
                record["origins"][0]["donor_metadata"]["zero_flux_filled_groups"] = [2]
            else:
                record["schema"] = "unknown"
            group[RECORD_DATASET][()] = json.dumps(record)
    before = stage.read_bytes()
    with pytest.raises(ValueError, match="zero-flux"):
        collapse_components(stage, output, groups=[("BOTH", ("HIGH", "LOW"))])
    assert stage.read_bytes() == before
    assert not output.exists()


@pytest.mark.parametrize("needs_fill", [False, True])
def test_refill_cannot_overwrite_collapsed_donor_records(tmp_path, needs_fill):
    source, output = tmp_path / "in.h5", tmp_path / "out.h5"
    _source(source)
    with h5py.File(source, "r+") as h5:
        _mark_filled(h5["mixtures/LOW"], [3], "donor.h5")
    _run("component", source, output)
    with h5py.File(output, "r+") as h5:
        if needs_fill:
            h5["mixtures/BOTH/total"][0] = 0
        records = h5[f"mixtures/BOTH/{RECORD_DATASET}"][()]
    before = output.read_bytes()
    macrolib = _touch_macrolib(tmp_path)
    with _fake_openmc(_fake_library()):
        if needs_fill:
            with pytest.raises(ValueError, match="perform substitution before collapse"):
                fill_zero_flux_groups(output, macrolib=macrolib, in_place=True)
            assert output.read_bytes() == before
        else:
            report = fill_zero_flux_groups(output, macrolib=macrolib, in_place=True)
            assert report.total_filled_bins == 0
            with h5py.File(output) as h5:
                assert h5[f"mixtures/BOTH/{RECORD_DATASET}"][()] == records
