"""Collapse is pre-equivalence, provenance-bound, and failure-safe."""

import json
from pathlib import Path

import h5py
import numpy as np
import pytest

from openmc2donjon import collapse_io, energy_collapse, mixture_collapse
from openmc2donjon.mgxs_input_contract import validate_input
from openmc2donjon.mgxs_input_uncertainty import UncertaintyConfig
from openmc2donjon.openmc_provenance import (
    collect_openmc_provenance,
    file_sha256,
    read_openmc_provenance,
    write_openmc_provenance,
)
from tests.test_collapse_scatter_layout import _source
from tests.test_openmc_provenance import OpenmcProvenanceTests


def _run(operation, source, output, **kwargs):
    if operation == "component":
        return mixture_collapse.collapse_components(
            source, output, groups=[("BOTH", ("HIGH", "LOW"))], **kwargs,
        )
    return energy_collapse.collapse_energy_groups(source, output, groups=((1, 2), (3, 4)), **kwargs)


@pytest.mark.parametrize("operation", ["component", "energy"])
@pytest.mark.parametrize("equivalence", ["sph", "SPH", "NSPH", "applied", "applied_text", "binding", "adf"])
def test_equivalence_rejected_without_mutation(tmp_path, operation, equivalence):
    source, output, summary = (tmp_path / name for name in ("source.h5", "output.h5", "summary.json"))
    _source(source)
    with h5py.File(source, "r+") as h5:
        if equivalence.startswith("applied"):
            h5.attrs["sph_applied"] = True if equivalence == "applied" else "true"
        elif equivalence == "binding":
            h5.attrs["sph_apply_operator"] = "divide-xs-by-nsph"
        else:
            for group in h5["mixtures"].values():
                group.create_dataset(equivalence, data=[.9, 1.1, .8, 1.2])
    before = source.read_bytes()
    output.write_bytes(b"previous output")
    summary.write_bytes(b"previous summary")
    with pytest.raises(ValueError, match="pre-equivalence"):
        _run(operation, source, output, force=True, summary_json=summary)
    assert source.read_bytes() == before
    assert output.read_bytes() == b"previous output"
    assert summary.read_bytes() == b"previous summary"


@pytest.mark.parametrize("operation", ["component", "energy"])
def test_identity_mapping_does_not_silently_discard_sph(tmp_path, operation):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    _source(source)
    with h5py.File(source, "r+") as h5:
        for group in h5["mixtures"].values():
            group.create_dataset("sph", data=[.9, 1.1, .8, 1.2])
    with pytest.raises(ValueError, match="original HDF5 directly"):
        if operation == "component":
            mixture_collapse.collapse_components(
                source, output, groups=[("HIGH", ("HIGH",)), ("LOW", ("LOW",))],
            )
        else:
            energy_collapse.collapse_energy_groups(source, output, groups=((1,), (2,), (3,), (4,)))
    assert not output.exists()


def _bind(source, tmp_path):
    files = OpenmcProvenanceTests._write_complete_case(tmp_path)
    record = collect_openmc_provenance(
        recipe_path=files["recipe"], statepoint_path=files["statepoint"], statepoint_loaded=True,
        declared_files={key: files[key] for key in ("geometry", "materials", "settings")},
        declared_metadata={"input_closure_complete": True},
    )
    return write_openmc_provenance(source, record)


@pytest.mark.parametrize("operation", ["component", "energy"])
def test_bound_input_remains_bound_and_repeated_history_is_retained(tmp_path, operation):
    source, output, repeat = (tmp_path / name for name in ("source.h5", "output.h5", "repeat.h5"))
    _source(source, layout="mixed_alias")
    original_record = _bind(source, tmp_path)
    source_bytes = source.read_bytes()
    assert validate_input(source).ok
    report = _run(operation, source, output)
    assert report["input_h5_sha256"] == file_sha256(source)
    assert report["output_h5_sha256"] == file_sha256(output)
    assert validate_input(output).ok
    record = read_openmc_provenance(output)
    assert record["integrity"]["ok"]
    assert record["status"] == original_record["status"]
    assert record["artifacts"] == original_record["artifacts"]
    assert record["handoff"]["payload_sha256"] != original_record["handoff"]["payload_sha256"]
    if operation == "component":
        energy_collapse.collapse_energy_groups(output, repeat, groups=((1, 2), (3, 4)))
    else:
        mixture_collapse.collapse_components(output, repeat, groups=[("BOTH", ("HIGH", "LOW"))])
    assert validate_input(repeat).ok
    assert read_openmc_provenance(repeat)["integrity"]["ok"]
    with h5py.File(repeat) as h5:
        history = json.loads(h5[collapse_io.HISTORY_PATH].asstr()[()])
        assert history["schema"] == collapse_io.HISTORY_SCHEMA
        assert len(history["steps"]) == 2
        assert history["steps"][0]["input_h5_sha256"] == file_sha256(source)
        assert history["steps"][1]["input_h5_sha256"] == file_sha256(output)
        assert history["steps"][0]["mapping"] == report["groups"]
    assert source.read_bytes() == source_bytes
    # Processing history is itself covered by the final payload binding.
    with h5py.File(repeat, "r+") as h5:
        h5[collapse_io.HISTORY_PATH][()] = '{"schema":"tampered","steps":[]}'
    assert not read_openmc_provenance(repeat)["integrity"]["ok"]


@pytest.mark.parametrize("operation", ["component", "energy"])
@pytest.mark.parametrize("corruption", ["payload", "missing_record", "malformed_record"])
def test_damaged_provenance_is_not_repaired(tmp_path, operation, corruption):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    _source(source)
    _bind(source, tmp_path)
    with h5py.File(source, "r+") as h5:
        if corruption == "payload":
            h5["mixtures/HIGH/total"][0] += .01
        elif corruption == "missing_record":
            del h5["provenance/openmc"]
        else:
            h5["provenance/openmc/record_json"][()] = "{}"
    output.write_bytes(b"preserve this")
    before = source.read_bytes()
    with pytest.raises(ValueError, match="provenance"):
        _run(operation, source, output, force=True)
    assert output.read_bytes() == b"preserve this"
    assert source.read_bytes() == before


@pytest.mark.parametrize("operation", ["component", "energy"])
@pytest.mark.parametrize("existing", [False, True])
@pytest.mark.parametrize("failure", ["write", "refresh", "summary", "publish"])
def test_late_failures_keep_old_files(tmp_path, monkeypatch, operation, existing, failure):
    source, output, summary = (tmp_path / name for name in ("source.h5", "output.h5", "summary.json"))
    _source(source)
    before = source.read_bytes()
    if existing:
        output.write_bytes(b"previous output")
        summary.write_bytes(b"previous summary")

    def fail(*args, **kwargs):
        raise RuntimeError("injected failure")

    if failure == "write":
        module = mixture_collapse if operation == "component" else energy_collapse
        name = "_write_component_datasets" if operation == "component" else "_write_collapsed_mixture"
        original = getattr(module, name)

        def fail_after_write(*args, **kwargs):
            original(*args, **kwargs)
            fail()

        monkeypatch.setattr(module, name, fail_after_write)
    elif failure == "refresh":
        monkeypatch.setattr(collapse_io, "refresh_openmc_provenance_after_hdf5_mutation", fail)
    else:
        original_replace, original_link = collapse_io.os.replace, collapse_io.os.link

        def replace(src, dest):
            if Path(dest) == (summary if failure == "summary" else output):
                fail()
            return original_replace(src, dest)

        def link(src, dest):
            if Path(dest) == output and failure == "publish":
                fail()
            return original_link(src, dest)

        monkeypatch.setattr(collapse_io.os, "replace", replace)
        monkeypatch.setattr(collapse_io.os, "link", link)
    with pytest.raises(RuntimeError, match="injected failure"):
        _run(operation, source, output, force=existing, summary_json=summary)
    assert source.read_bytes() == before
    if existing:
        assert output.read_bytes() == b"previous output"
        assert summary.read_bytes() == b"previous summary"
    else:
        assert not output.exists()
        assert not summary.exists()
    assert not list(tmp_path.glob(".openmc2donjon-collapse*"))


@pytest.mark.parametrize("operation", ["component", "energy"])
@pytest.mark.parametrize("alias", ["output", "summary_input", "summary_output", "hardlink"])
def test_input_and_output_aliases_are_rejected(tmp_path, operation, alias):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    _source(source)
    before = source.read_bytes()
    if alias == "hardlink":
        output.hardlink_to(source)
    summary = source if alias == "summary_input" else output if alias == "summary_output" else None
    with pytest.raises(ValueError, match="different"):
        _run(operation, source, source if alias == "output" else output, force=True,
             summary_json=summary)
    assert source.read_bytes() == before


def _add_uncertainty(source):
    with h5py.File(source, "r+") as h5:
        for group in h5["mixtures"].values():
            group.attrs["fissionable"] = True
            group["fission"][:] = .1
            group["nu_fission"][:] = .25
            group["chi"][:] = [.4, .3, .2, .1]
            group.create_dataset("kappa_fission", data=np.arange(1, 5, dtype=float))
            for name in tuple(group):
                if name == "scatter_matrix" or name.endswith("_std_dev"):
                    continue
                group.create_dataset(f"{name}_std_dev", data=np.abs(group[name][:]) * .01)


def test_energy_collapse_keeps_full_xs_uncertainty_coverage(tmp_path):
    source, output, summary = (tmp_path / name for name in ("source.h5", "output.h5", "summary.json"))
    _source(source)
    _add_uncertainty(source)
    cfg = UncertaintyConfig(require_coverage=True)
    before = validate_input(source, uncertainty=cfg)
    assert before.ok, before.issues
    report = _run("energy", source, output, summary_json=summary)
    after = validate_input(output, uncertainty=cfg)
    assert after.ok, after.issues
    assert before.uncertainty_datasets == after.uncertainty_datasets
    assert after.uncertainty_datasets == after.uncertainty_expected_datasets
    assert json.loads(summary.read_text()) == report
    with h5py.File(source) as src, h5py.File(output) as out:
        fine_flux = src["openmc_volume_flux"][0]
        for name in energy_collapse.VECTOR_XS:
            if name not in src["mixtures/HIGH"]:
                continue
            std = src[f"mixtures/HIGH/{name}_std_dev"][:]
            expected = [sum(std[g] * fine_flux[g] for g in pair) / sum(fine_flux[list(pair)])
                        for pair in ((0, 1), (2, 3))]
            np.testing.assert_allclose(out[f"mixtures/HIGH/{name}_std_dev"][:], expected)
        # Chi=(.7,.3); first-order L1 ratio bounds retain denominator correlation.
        np.testing.assert_allclose(out["mixtures/HIGH/chi_std_dev"][:], [.0042, .0042])
        assert out["mixtures/HIGH/chi_std_dev"].attrs["energy_uncertainty_method"] == (
            "first-order-l1-normalized-chi-no-covariance"
        )


@pytest.mark.parametrize("operation", ["component", "energy"])
@pytest.mark.parametrize("field", ["total", "chi", "openmc_volume_flux"])
@pytest.mark.parametrize("bad", ["negative", "nan", "shape"])
def test_invalid_uncertainty_is_rejected_without_mutation(tmp_path, operation, field, bad):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    _source(source)
    _add_uncertainty(source)
    with h5py.File(source, "r+") as h5:
        group = h5 if field == "openmc_volume_flux" else h5["mixtures/HIGH"]
        name = f"{field}_std_dev"
        if name in group:
            del group[name]
        values = (np.zeros((1,)) if bad == "shape"
                  else np.full(group[field].shape, -.1 if bad == "negative" else np.nan))
        group.create_dataset(name, data=values)
    before = source.read_bytes()
    output.write_bytes(b"previous output")
    with pytest.raises(ValueError, match="finite non-negative"):
        _run(operation, source, output, force=True)
    assert source.read_bytes() == before
    assert output.read_bytes() == b"previous output"


@pytest.mark.parametrize("operation", ["component", "energy"])
def test_false_applied_flag_does_not_block_uncorrected_data(tmp_path, operation):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    _source(source)
    with h5py.File(source, "r+") as h5:
        h5.attrs["sph_applied"] = False
    _run(operation, source, output)
    assert validate_input(output).ok


@pytest.mark.parametrize("operation", ["component", "energy"])
def test_orphan_uncertainty_is_rejected(tmp_path, operation):
    source, output = tmp_path / "source.h5", tmp_path / "output.h5"
    _source(source)
    with h5py.File(source, "r+") as h5:
        h5["mixtures/LOW"].create_dataset("missing_std_dev", data=np.ones(4))
    with pytest.raises(ValueError, match="no matching mean"):
        _run(operation, source, output)
    assert not output.exists()
