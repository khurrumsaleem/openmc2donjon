"""Small native mixed-P1/P3 transfer test; not a reactor physics benchmark."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import time

import h5py
import numpy as np
import openmc
import openmc.mgxs

from openmc2donjon.export_openmc_mgxs import DomainExportSpec, export_openmc_mgxs_library
from openmc2donjon.macrolib import convert_mgxs_hdf5_to_macrolib, read_macrolib_ascii
from openmc2donjon.mgxs_input_contract import validate_input
from openmc2donjon.mgxs_input_report import print_report
from openmc2donjon.multicompo import convert_mgxs_hdf5, read_mgxs_hdf5
from openmc2donjon.openmc_mixed_order import exact_cell_mgxs_tallies, set_domain_scatter_orders


SCATTER = "consistent nu-scatter matrix"


def build():
    carbon, lead = openmc.Material(name="carbon"), openmc.Material(name="lead")
    carbon.add_nuclide("C12", 1.0)
    carbon.set_density("g/cm3", 1.7)
    lead.add_nuclide("Pb208", 1.0)
    lead.set_density("g/cm3", 11.34)
    box = openmc.model.RectangularParallelepiped(-10, 10, -10, 10, -10, 10,
                                                boundary_type="vacuum")
    mid = openmc.XPlane(x0=0)
    cells = [openmc.Cell(name="carbon_P1", fill=carbon, region=-box & -mid),
             openmc.Cell(name="lead_P3", fill=lead, region=-box & +mid)]
    geometry = openmc.Geometry(cells)
    library = openmc.mgxs.Library(geometry)
    library.domain_type, library.domains = "cell", cells
    library.energy_groups = openmc.mgxs.EnergyGroups([1e-5, 1e5, 2e7])
    library.mgxs_types = ["total", "absorption", "reduced absorption", "nu-transport", SCATTER]
    library.correction, library.legendre_order = None, 3
    library.build_library()
    set_domain_scatter_orders(library, {cells[0].id: 1, cells[1].id: 3})
    tallies = exact_cell_mgxs_tallies(library)
    # Check scoring scope, not only the final padded array shape.
    high = [t for t in tallies if any(isinstance(f, openmc.LegendreFilter) and f.order > 1
                                    for f in t.filters)]
    assert high
    assert all(list(t.find_filter(openmc.CellFilter).bins) == [cells[1].id] for t in high)
    settings = openmc.Settings()
    settings.run_mode, settings.batches, settings.particles = "fixed source", 20, 4000
    settings.seed = 271828
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Box((-10, -10, -10), (10, 10, 10)),
        energy=openmc.stats.Discrete([1e3, 2e6], [0.5, 0.5]))
    return openmc.Model(geometry, openmc.Materials([carbon, lead]), settings, tallies), library


def check_macrolib(path, mixtures):
    macro = read_macrolib_ascii(path)
    assert macro.nmixtures == 2 and macro.ngroups == 2
    assert set(macro.scatter) == {0, 1, 2, 3}, "P2/P3 lost during transfer"
    maximum = 0.0
    for moment in range(4):
        expected = np.stack([mix.scatter_matrix[moment] for mix in mixtures])
        np.testing.assert_allclose(macro.scatter[moment], expected, rtol=3e-6, atol=1e-12)
        maximum = max(maximum, float(np.max(np.abs(macro.scatter[moment] - expected))))
    return maximum


def donjon_ingest(root, work, mixtures):
    root = root.resolve()
    if not (root / "rdonjon").is_file():
        raise ValueError("--donjon-root must name the Donjon directory containing rdonjon")
    parent = root / "data" / "openmc2donjon"
    parent.mkdir(parents=True, exist_ok=True)
    run = Path(tempfile.mkdtemp(prefix="mixed_pn_", dir=parent))
    stem = run.name
    deck = run / f"{stem}.x2m"
    # CLE-2000 limits both line and character-token lengths. Keep file names
    # short even when the user's work directory has a long absolute path.
    with tempfile.TemporaryDirectory(prefix="o2dpn_", dir="/tmp") as io_name:
        io = Path(io_name)
        shutil.copy2(work / "mixed.mcompo.txt", io / "input.txt")
        output = io / "output.txt"
        deck.write_text(f"""MODULE NCR: END: ;
LINKED_LIST CPO MACRO ;
SEQ_ASCII CPO_ASC ::
 FILE '{io / 'input.txt'}' ;
SEQ_ASCII MAC_ASC ::
 FILE '{output}' ;
CPO := CPO_ASC ;
MACRO := NCR: CPO :: EDIT 1 MACRO NMIX 2
 COMPO CPO CPO
 MIX 1 USE ENDMIX
 MIX 2 USE ENDMIX ;
MAC_ASC := MACRO ;
END: ;
""")
        shutil.copy2(deck, work / "donjon_ingest.x2m")
        with (work / "donjon.log").open("w") as log:
            result = subprocess.run(["./rdonjon", "-q", str(deck.relative_to(root / "data"))],
                                    cwd=root, stdout=log, stderr=subprocess.STDOUT, timeout=120)
        # The wrapper's architecture directory varies between installations.
        listings = list(root.glob(f"*/{stem}.result"))
        if len(listings) != 1:
            raise RuntimeError("Expected exactly one DONJON listing; inspect donjon.log")
        listing = listings[0]
        shutil.copy2(listing, work / "donjon.result")
        if result.returncode != 0 or "normal end" not in listing.read_text().lower():
            raise RuntimeError("DONJON did not finish normally; inspect donjon.log/result")
        shutil.copy2(output, work / "donjon.macrolib.txt")
    return check_macrolib(work / "donjon.macrolib.txt", mixtures)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--openmc", default="openmc")
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--donjon-root", type=Path)
    args = parser.parse_args()
    work = args.work_dir.resolve()
    if work.exists() and any(work.iterdir()):
        parser.error("Use a new or empty work directory; existing results are not overwritten")
    if args.threads < 1:
        parser.error("--threads must be positive")
    work.mkdir(parents=True, exist_ok=True)
    receipt = {"schema": "openmc2donjon.mixed-order-diagnostic.v1", "status": "incomplete",
               "production_physics_accepted": False, "donjon_solve_completed": False,
               "openmc_version": openmc.__version__, "particle_histories": 80000,
               "threads": args.threads, "source_orders": {"carbon_P1": 1, "lead_P3": 3},
               "purpose": "Native scoring, export and optional DONJON NCR ingestion; no eigenvalue test"}
    start = time.monotonic()
    try:
        model, library = build()
        model.export_to_xml(directory=work)
        command = [args.openmc, "-s", str(args.threads)]
        receipt["openmc_command"] = command
        with (work / "openmc.log").open("w") as log:
            subprocess.run(command, cwd=work, stdout=log, stderr=subprocess.STDOUT,
                           timeout=300, check=True)
        with openmc.StatePoint(work / "statepoint.20.h5") as sp:
            library.load_from_statepoint(sp)
        hdf = work / "mgxs.h5"
        export_openmc_mgxs_library(library, hdf, scatter_mgxs_type=SCATTER,
            domain_specs=[DomainExportSpec(cell, name=cell.name, volume=4000.0)
                          for cell in library.domains], root_attrs={"domain_mode": "component"})
        with h5py.File(hdf) as handle:
            for cell in library.domains:
                order = library.domain_legendre_orders[cell.id]
                group = handle[f"mixtures/{cell.name}"]
                assert group.attrs["source_legendre_order"] == order
                for value, suffix in (("mean", ""), ("std_dev", "_std_dev")):
                    native = library.get_mgxs(cell, SCATTER).get_xs(
                        xs_type="macro", moment="all", value=value)
                    expected = np.zeros((4, 2, 2))
                    expected[:order+1] = np.moveaxis(native, -1, 0)
                    np.testing.assert_array_equal(group["scatter_matrix"+suffix][:], expected)
            assert np.any(handle["mixtures/lead_P3/scatter_matrix"][2:] != 0)
        report = validate_input(hdf, require_volume=True, require_transport_dataset=True)
        print_report(report)
        if not report.ok:
            raise RuntimeError(f"Preflight failed: {report.issues}")
        receipt["preflight_warnings"] = report.warnings
        mixtures, _ = read_mgxs_hdf5(hdf)
        convert_mgxs_hdf5(hdf, work / "mixed.mcompo.txt")
        convert_mgxs_hdf5_to_macrolib(hdf, work / "mixed.macrolib.txt")
        receipt["macrolib_max_abs_scatter_error"] = check_macrolib(work / "mixed.macrolib.txt", mixtures)
        receipt["donjon_ingest_checked"] = args.donjon_root is not None
        if args.donjon_root:
            receipt["donjon_max_abs_scatter_error"] = donjon_ingest(args.donjon_root, work, mixtures)
        receipt["status"] = "passed"
    except Exception as exc:
        receipt["status"], receipt["error"] = "failed", str(exc)
        raise
    finally:
        receipt["elapsed_seconds"] = time.monotonic() - start
        receipt["cross_sections_xml"] = os.environ.get("OPENMC_CROSS_SECTIONS")
        receipt["sha256"] = {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                             for path in work.iterdir() if path.is_file() and path.name != "receipt.json"}
        (work / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(f"PASS: mixed native orders retained; diagnostic only. Receipt: {work / 'receipt.json'}")


if __name__ == "__main__":
    main()
