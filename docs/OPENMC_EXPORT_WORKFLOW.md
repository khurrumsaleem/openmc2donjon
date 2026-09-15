# OpenMC Export Workflow

This is the production-facing path for creating the openmc2donjon HDF5 handoff
from a real OpenMC MGXS run.

## Recipe-Based Statepoint Export

Use `openmc2donjon-export` with a small Python recipe:

```sh
openmc2donjon-export \
  --recipe export_recipe.py \
  --statepoint statepoint.120.h5 \
  -o mgxs_library.h5
```

The recipe owns OpenMC-specific modeling details: geometry loading, group
structure, MGXS domain type, spatial domain partition, and stable mixture names.
The CLI owns the package handoff: loading the statepoint into the recipe's
library and writing the documented HDF5 contract.

For a ready-to-edit starting point, copy
[`examples/openmc_recipe_template/`](../examples/openmc_recipe_template/) into a
case directory and edit `export_recipe.py`.

Check the runtime environment and recipe import path:

```sh
openmc2donjon doctor --recipe export_recipe.py
```

Before running a long OpenMC job, dry-run the recipe:

```sh
openmc2donjon-export --recipe export_recipe.py --no-load-statepoint --dry-run
openmc2donjon-export --recipe export_recipe.py --no-load-statepoint --dry-run --strict-dry-run
```

The dry-run builds the recipe library, reports the group count, Legendre order,
domain type, MGXS types, root attributes, and the first mixture names. It also
prints a production checklist for required MGXS types, transport availability,
domain mapping, volume provenance, and `domain_mode`. It does not read MGXS
values or write an HDF5 file. Use `--strict-dry-run` in automation when recipe
warnings should fail early, before an expensive OpenMC run is launched.

Use the same recipe to write the OpenMC MGXS tallies before running OpenMC:

```sh
openmc2donjon-export \
  --recipe export_recipe.py \
  --write-tallies tallies.xml
```

The recipe can define `extra_tallies(...)` to append case-specific non-MGXS
tallies, such as surface-current tallies used later for ADF/DF generation.

The recipe is a separate case adapter, not a patch that must be pasted into the
OpenMC driver, and it may live anywhere addressable by `--recipe`. An existing
`main.py` needs no openmc2donjon-specific lines when it:

1. exports the same materials/geometry/settings used by the recipe;
2. leaves the generated `tallies.xml` in place; and
3. runs OpenMC in that XML directory.

If `main.py` instead constructs an in-memory `openmc.Model`, assigns
`model.tallies`, or exports XML after `--write-tallies`, it may ignore or
overwrite the generated MGXS tallies. In that case, add the recipe library's
MGXS tallies to the model or change the order to `export model XML -> generate
tallies.xml -> run OpenMC`. The recipe library's geometry and domain IDs must
match the statepoint. Adding a recipe after the transport run cannot recover
tallies that the statepoint never recorded.

To export and immediately write DONJON ASCII in one command:

```sh
openmc2donjon-from-openmc \
  --recipe export_recipe.py \
  --dry-run \
  --run-dir runs/case1 \
  --check \
  --strict-dry-run
```

The one-step dry-run reports the same recipe/domain metadata plus the planned
DONJON format, ASCII output, HDF5 handoff path, summary paths, and preflight
requirements. It does not write HDF5, summary JSON, or DONJON ASCII files.

```sh
openmc2donjon-from-openmc \
  --recipe export_recipe.py \
  --statepoint statepoint.120.h5 \
  --run-dir runs/case1 \
  --check
```

With `--run-dir`, the command writes `mgxs_library.h5`, `out.mcompo.txt`,
`run_summary.json`, `openmc_provenance.json`, optional `check_summary.json`, a
recipe/model-source copy, and `manifest.json`. The summary JSON records recipe,
statepoint, HDF5, output,
group count, Legendre order, and mixture names. The summary schema is documented
in [From-OpenMC summary JSON](FROM_OPENMC_SUMMARY_SCHEMA.md). Existing managed
run-directory files are refused unless `--force-run-dir` is set.

Before claiming academic replay, declare every model dependency through the
recipe `provenance_files()` hook and set `input_closure_complete=true` only
after that list is complete. The final handoff binds both those source hashes
and a recomputed digest of the actual MGXS HDF5 numerical payload.

Additional production side artifacts can be copied into the same manifest
during the one-step run:

```sh
openmc2donjon-from-openmc \
  --recipe export_recipe.py \
  --statepoint statepoint.120.h5 \
  --run-dir runs/case1 \
  --extra-artifact surface-flux=runs/case1/openmc_surface_flux.h5 \
  --extra-artifact low-order-driver=runs/case1/low_order_driver.h5 \
  --extra-artifact homogeneous-face-flux=runs/case1/homogeneous_face_flux.h5 \
  --check
```

To add extra files to an existing handoff manifest:

```sh
openmc2donjon bundle \
  --output-dir runs/case1 \
  --mgxs runs/case1/mgxs_library.h5 \
  --mcompo runs/case1/out.mcompo.txt \
  --run-summary runs/case1/run_summary.json \
  --extra notes=notes.txt \
  --force
```

If ADF/DF values are produced by a separate OpenMC or nodal post-processing
step, inject them into the HDF5 handoff before conversion:

```sh
openmc2donjon make-adf-sidecar runs/case1/mgxs_library.h5 \
  -o runs/case1/adf_sidecar.h5 \
  --mode unity

openmc2donjon augment-adf runs/case1/mgxs_library.h5 \
  --adf-source runs/case1/adf_sidecar.h5 \
  -o runs/case1/mgxs_with_adf.h5 \
  --faces FD_XMIN,FD_XMAX,FD_YMIN,FD_YMAX \
  --summary-json runs/case1/adf_summary.json
```

The generated unity sidecar is marked `adf_real=false`; it verifies the
interface and should be replaced by case-specific physics ADF/DF values for
production neutronics.

When heterogeneous and homogeneous face fluxes are available, generate the
sidecar directly from their ratio:

```sh
openmc2donjon export-surface-flux statepoint.120.h5 \
  --mgxs runs/case1/mgxs_library.h5 \
  -o runs/case1/openmc_surface_flux.h5 \
  --tally-name openmc2donjon_surface_current_mu \
  --mesh-shape 1,2 \
  --mu-edges 0.0,0.25,0.5,0.75,1.0 \
  --face-area 4.0

openmc2donjon make-low-order-driver runs/case1/mgxs_library.h5 \
  -o runs/case1/low_order_driver.h5 \
  --raw-driver runs/case1/raw_low_order_driver.h5 \
  --net-current-sign-convention auto \
  --faces FD_XMIN,FD_XMAX,FD_YMIN,FD_YMAX

openmc2donjon check-low-order-driver \
  runs/case1/mgxs_library.h5 runs/case1/low_order_driver.h5 \
  --faces FD_XMIN,FD_XMAX,FD_YMIN,FD_YMAX \
  --face-widths 4.0

openmc2donjon make-homogeneous-face-flux runs/case1/mgxs_library.h5 \
  -o runs/case1/homogeneous_face_flux.h5 \
  --volume-flux runs/case1/low_order_driver.h5 \
  --net-current runs/case1/low_order_driver.h5 \
  --faces FD_XMIN,FD_XMAX,FD_YMIN,FD_YMAX \
  --face-widths 4.0

openmc2donjon check-face-flux runs/case1/mgxs_library.h5 \
  --surface-flux runs/case1/openmc_surface_flux.h5 \
  --homogeneous-face-flux runs/case1/homogeneous_face_flux.h5 \
  --faces FD_XMIN,FD_XMAX,FD_YMIN,FD_YMAX

openmc2donjon make-adf-sidecar runs/case1/mgxs_library.h5 \
  -o runs/case1/adf_sidecar.h5 \
  --mode flux-ratio \
  --surface-flux runs/case1/openmc_surface_flux.h5 \
  --homogeneous-face-flux runs/case1/homogeneous_face_flux.h5 \
  --faces FD_XMIN,FD_XMAX,FD_YMIN,FD_YMAX
```

For the one-step production path, let `openmc2donjon-from-openmc` build and
inject the sidecar inside the managed run directory. It also writes and bundles
the surface-flux, low-order driver, low-order contract check, homogeneous
face-flux, face-flux contract check, and ADF-sidecar summaries:

```sh
openmc2donjon-from-openmc \
  --recipe export_recipe.py \
  --statepoint statepoint.120.h5 \
  --run-dir runs/case1 \
  --build-flux-ratio-adf \
  --export-surface-flux \
  --surface-flux-tally-name openmc2donjon_surface_current_mu \
  --surface-flux-mesh-shape 1,2 \
  --surface-flux-mu-edges 0.0,0.25,0.5,0.75,1.0 \
  --surface-flux-face-area 4.0 \
  --low-order-raw-driver raw_low_order_driver.h5 \
  --low-order-net-current-sign-convention auto \
  --adf-faces FD_XMIN,FD_XMAX,FD_YMIN,FD_YMAX \
  --adf-face-widths 4.0 \
  --adf-invalid-fill 1.0 \
  --require-volume \
  --require-transport-dataset
```

If a trusted low-order or nodal solve has already written the homogeneous
face-flux denominator, the same one-step path can consume it directly. In that
mode the run directory bundles the external denominator and skips the
low-order-driver and homogeneous-face-flux reconstruction artifacts:

```sh
openmc2donjon-from-openmc \
  --recipe export_recipe.py \
  --statepoint statepoint.120.h5 \
  --run-dir runs/case1 \
  --build-flux-ratio-adf \
  --export-surface-flux \
  --surface-flux-tally-name openmc2donjon_surface_current_mu \
  --surface-flux-mesh-shape 1,2 \
  --surface-flux-mu-edges 0.0,0.25,0.5,0.75,1.0 \
  --surface-flux-face-area 4.0 \
  --homogeneous-face-flux runs/case1/homogeneous_face_flux.h5 \
  --adf-faces FD_XMIN,FD_XMAX,FD_YMIN,FD_YMAX \
  --adf-invalid-fill 1.0 \
  --require-volume \
  --require-transport-dataset
```

`--adf-surface-flux` and `--homogeneous-face-flux` both accept
`FILE::/dataset/path` references when the required datasets live under
case-specific names.

For a compact, deterministic example of an external low-order handoff with
case-specific dataset paths, mixture/face reordering, and positive-inward
current conversion, run:

```sh
bash examples/external_low_order_handoff/run_smoke.sh
```

If the external solver already writes homogeneous face fluxes, use the direct
face-flux adapter contract instead:

```sh
bash examples/external_face_flux_adapter/run_smoke.sh
```

The formal denominator contract is in `docs/EXTERNAL_FACE_FLUX_CONTRACT.md`.

If the ADF sidecar was produced separately, pass it directly:

```sh
openmc2donjon-from-openmc \
  --recipe export_recipe.py \
  --statepoint statepoint.120.h5 \
  --run-dir runs/case1 \
  --adf-source runs/case1/adf_sidecar.h5 \
  --adf-faces FD_XMIN,FD_XMAX,FD_YMIN,FD_YMAX \
  --check --require-adf
```

## OpenMC-Side SPH From CE/MG Fluxes

The standard production SPH route uses OpenMC MG as the equivalence operator
upstream of Converter and DONJON. Hold the heterogeneous OpenMC CE model fixed
as the fine reference, build a homogenized OpenMC MG coarse model, and iterate
the rate-preserving update against the MG solution.

The CE and MG geometries are deliberately different. The CE geometry resolves
the physical heterogeneity; the MG geometry replaces each declared comparison
domain with its homogenized region. The CE tallies must use the MG transport
group boundaries, and both calculations must represent the same physical state
and boundary conditions. A conservative mapping assigns each fine volume and
its integrated flux/rates to exactly one coarse comparison domain. A shared
`(mixture, group)` output shape proves array alignment, not geometric identity.

After convergence, apply the factors to the Converter-facing HDF5 and pass that
corrected file through Converter. DONJON verifies the converted object; it is
not the standard SPH feedback operator.

Export the region/group flux tally from each OpenMC statepoint with the same
MGXS handoff metadata. The fine and coarse models may use different native
OpenMC IDs. Set `CE_DOMAIN_IDS` and `MG_DOMAIN_IDS` to comma-separated IDs in
the canonical `mixture_names` order; the exporter verifies each model's own
mapping before writing that common order:

```sh
openmc2donjon export-volume-flux ce_statepoint.h5 \
  --mgxs runs/case1/mgxs_library.h5 \
  --tally-name openmc_ce_volume_flux \
  --dataset-name openmc_volume_flux \
  --source-domain-ids "${CE_DOMAIN_IDS}" \
  -o runs/case1/openmc_ce_flux.h5 \
  --summary-json runs/case1/openmc_ce_flux_summary.json

openmc2donjon export-volume-flux mg_statepoint.h5 \
  --mgxs runs/case1/mgxs_library.h5 \
  --tally-name openmc_mg_volume_flux \
  --dataset-name openmc_mg_flux \
  --source-domain-ids "${MG_DOMAIN_IDS}" \
  -o runs/case1/openmc_mg_flux.h5 \
  --summary-json runs/case1/openmc_mg_flux_summary.json
```

The exporter writes datasets in `(mixture, group)` order and tags them with
`group_order=mgxs_donjon` plus the canonical `mixture_names` list. OpenMC
EnergyFilter tally bins are reversed to the high-to-low group order used by
the MGXS handoff and DONJON.

Then build the OpenMC-side SPH factors. Set
`CE_FLUX_REL_SIGMA_LIMIT` and `MG_FLUX_REL_SIGMA_LIMIT` from the project's
predeclared statistical uncertainty budget; neither is a universal product
constant:

```sh
openmc2donjon make-openmc-sph-sidecar runs/case1/mgxs_library.h5 \
  -o runs/case1/openmc_sph_sidecar.h5 \
  --reference-flux runs/case1/openmc_ce_flux.h5::openmc_volume_flux \
  --mg-flux runs/case1/openmc_mg_flux.h5::openmc_mg_flux \
  --table-output runs/case1/openmc_sph.csv \
  --damping 0.5 \
  --sph-target rate \
  --flux-normalization power \
  --require-reference-flux-std-dev \
  --max-reference-flux-std-dev-rel "${CE_FLUX_REL_SIGMA_LIMIT}" \
  --require-mg-flux-std-dev \
  --max-mg-flux-std-dev-rel "${MG_FLUX_REL_SIGMA_LIMIT}" \
  --summary-json runs/case1/openmc_sph_summary.json
```

For each next OpenMC MG iteration, apply `XS / NSPH` to the OpenMC-native
`setN` library, rerun the homogenized MG model, and recompute the sidecar with
`--previous-sph` until the declared update residual passes:

```sh
openmc2donjon apply-sph runs/case1/mg_case/mgxs_unapplied.h5 \
  --input-format openmc-mgxs \
  --sph-source runs/case1/openmc_sph_sidecar.h5 \
  -o runs/case1/mg_case/mgxs.h5
```

The OpenMC-native `setN` file is a different serialization from the
Converter-layout sidecar input, so this intermediate application is explicitly
marked `openmc-mgxs-intermediate-unbound`. It is valid only for the next MG
iteration and cannot satisfy the strict Converter physical-SPH gate. The final
application must target the exact Converter-layout HDF5 bound by the sidecar
SHA-256.

Once converged, apply the final factors to the Converter-facing HDF5. This
rewrites the macroscopic cross sections; `augment-sph` only attaches factor
records and is not a substitute for this production step:

```sh
openmc2donjon apply-sph runs/case1/mgxs_library.h5 \
  --sph-source runs/case1/openmc_sph_sidecar.h5 \
  -o runs/case1/mgxs_sph_applied.h5 \
  --summary-json runs/case1/sph_apply_summary.json

openmc2donjon runs/case1/mgxs_sph_applied.h5 \
  -o runs/case1/out.mcompo.txt --production --require-physical-sph \
  --summary-json runs/case1/out.mcompo.txt.convert.json
```

For a portable fixture-backed check of this route:

```sh
bash examples/openmc_sph_sidecar_minicase/run_smoke.sh
```

For a small workflow check before using a real OpenMC model:

```sh
bash scripts/run_recipe_export_smoke.sh
```

That smoke uses `examples/recipe_export_smoke/minimal_recipe.py`, a tiny
MGXS-like recipe that exercises the same CLI hooks without being a physics
benchmark.

Minimal recipe shape:

```python
from pathlib import Path

import openmc
import openmc.mgxs as mgxs

from openmc2donjon import DomainExportSpec

CELL_VOLUMES_CM3 = {
    101: 1.0,
}


def build_library():
    materials = openmc.Materials.from_xml("materials.xml")
    geometry = openmc.Geometry.from_xml("geometry.xml", materials=materials)

    library = mgxs.Library(geometry)
    library.energy_groups = mgxs.EnergyGroups([1.0e-5, 1.0, 1.0e3, 1.0e7])
    library.mgxs_types = [
        "total",
        "absorption",
        "fission",
        "nu-fission",
        "chi",
        "scatter matrix",
        "transport",
    ]
    library.domain_type = "cell"
    library.domains = list(geometry.get_all_cells().values())
    library.by_nuclide = False
    library.legendre_order = 1
    library.correction = None
    library.build_library()
    return library


def domain_names(library):
    return {cell.id: f"CELL_{cell.id}" for cell in library.domains}


def domain_specs(library):
    return [
        DomainExportSpec(
            domain=cell,
            name=f"CELL_{cell.id}",
            volume=CELL_VOLUMES_CM3[cell.id],
            attrs={"source_domain_id": int(cell.id)},
        )
        for cell in library.domains
    ]


def root_attrs():
    return {"domain_mode": "cell"}
```

The general static policy uses `"absorption"` paired with ordinary
`"scatter matrix"` (or the ordinary `"consistent scatter matrix"` estimator).
Set `library.correction = None` before `library.build_library()` for either
policy. OpenMC otherwise defaults to a diagonal transport correction when the
Legendre order is zero. Converter writes raw total XS, so it rejects active
P0-corrected scattering rather than mixing corrected and uncorrected terms.
This setting does not disable a separately tallied `transport` or `nu-transport`
used for the diffusion coefficient.

For a fast-spectrum static calculation in which
multiplicative `(n,xn)` scattering must contribute to neutron balance, include
the complete matched set below in `library.mgxs_types`:

```python
"absorption",
"reduced absorption",
"consistent nu-scatter matrix",
"nu-transport",
```

and select the nu-weighted matrix explicitly with a recipe hook:

```python
def scatter_mgxs_type():
    return "consistent nu-scatter matrix"
```

The ordinary and nu-weighted policies are distinct. The exporter does not infer
one from the energy range, silently substitute one matrix for another, or
construct reduced absorption from independent reaction-rate arithmetic. For
the fast-spectrum policy, OpenMC supplies `reduced absorption`, the
`consistent nu-scatter matrix`, and the matching `nu-transport`; the exporter
checks and preserves that complete set. Ordinary scattering instead pairs with
ordinary `transport`.
Do not manually edit the absorption vector or scattering matrix.

Converter writes `NTOT0` and the selected `SCAT` records. DONJON therefore sees
the static net absorption implicitly as total minus the P0 scattering row sum.
The reduced-absorption/consistent-nu-scatter pair retains neutron multiplication
from `(n,xn)` in this balance. It does not create separate `N2N` or `N3N`
depletion-reaction records; generating those records is outside this static
macroscopic handoff.

For mesh or other subdomain exports, return explicit `DomainExportSpec` objects:

```python
def domain_specs(library):
    mesh = library.domains[0]
    return [
        DomainExportSpec(
            domain=mesh,
            name="ASM_Y01_X01",
            xs_kwargs={"subdomains": [(1, 1, 1)]},
            volume=assembly_volume,
            attrs={"mesh_index": [1, 1, 1]},
        ),
    ]
```

The CLI performs the default OpenMC load:

```python
with openmc.StatePoint(str(statepoint_path)) as sp:
    library.load_from_statepoint(sp)
```

If a case needs custom loading, define this in the recipe:

```python
def load_statepoint(library, statepoint_path):
    with openmc.StatePoint(str(statepoint_path)) as sp:
        library.load_from_statepoint(sp)
        print(f"keff = {sp.keff}")
```

Optional recipe hooks:

| Function | Purpose |
| --- | --- |
| `build_library()` | Required. Return an OpenMC `mgxs.Library`-like object. |
| `domain_specs(library)` | Optional. Return explicit spatial subdomain specs. |
| `domain_names(library)` | Optional. Return stable names keyed by domain object, id, or name. |
| `root_attrs(library)` | Optional. Return root HDF5 attributes such as `domain_mode`. |
| `scatter_mgxs_type(library)` | Optional. Explicitly select a non-default scattering MGXS type. |
| `load_statepoint(library, statepoint_path)` | Optional. Override default OpenMC statepoint loading. |
| `provenance_files(...)` | Optional but strongly recommended. Map stable roles to every Python/XML/CAD/mesh file needed to rebuild the fine model. Small files are copied into a managed run bundle. |
| `provenance_metadata(...)` | Optional. Record launcher-only facts such as the actual thread and MPI-rank counts. Values must come from the original run receipt, not the current export shell. |
| `postprocess_hdf5(output_path, library)` | Optional. Add case-specific payloads such as ADF. |

Optional hooks may declare only the arguments they need. Supported argument names
are `library`, `recipe_path`, `statepoint_path`, `output_path`, and `summary`.

## OpenMC Reference Flux for SPH Loops

The standard fixed-reference OpenMC CE/MG workflow needs a CE volume-flux
matrix projected conservatively from the heterogeneous fine geometry onto the
declared comparison domains. Its shape is `(mixture, group)`, and its
comparison-domain/group order matches the homogenized MG flux and MGXS handoff.
A recipe can add that matrix in `postprocess_hdf5`:

```python
def postprocess_hdf5(output_path, summary):
    import h5py

    names = [domain.name for domain in summary.domains]
    with h5py.File(output_path, "a") as h5:
        flux = h5.create_dataset("openmc_volume_flux", data=openmc_flux)
        flux.attrs["mixture_names"] = names
        flux.attrs["group_order"] = "mgxs_donjon"

        std_dev = h5.create_dataset(
            "openmc_volume_flux_std_dev",
            data=openmc_flux_std_dev,
        )
        std_dev.attrs["mixture_names"] = names
        std_dev.attrs["group_order"] = "mgxs_donjon"
        std_dev.attrs["std_dev_of"] = "openmc_volume_flux"
```

OpenMC-side SPH tooling can use this payload as the CE reference flux. When the
sibling `openmc_volume_flux_std_dev` dataset is present, preserve it in the SPH
derivation record so the corrected HDF5 or sidecar can prove reference-flux
uncertainty coverage. Enable strict production checks when the case policy
requires both MGXS and reference-flux uncertainty to be present:

```sh
openmc2donjon make-openmc-sph-sidecar mgxs_library.h5 \
  -o sph_sidecar.h5 \
  --reference-flux openmc_ce_flux.h5::openmc_volume_flux \
  --mg-flux openmc_mg_flux.h5::openmc_volume_flux \
  --table-output sph_openmc_ce_mg.csv \
  --sph-target rate \
  --flux-normalization power
```

These are two different uncertainty paths. MGXS `*_std_dev` datasets audit the
cross sections exported to DONJON; `openmc_volume_flux_std_dev` audits the
OpenMC CE reference flux used to compute SPH factors.

The command above writes both the auditable CSV table and the HDF5 sidecar.
After convergence and review, apply it to the Converter-facing HDF5:

```sh
openmc2donjon apply-sph mgxs_library.h5 \
  --sph-source sph_sidecar.h5 \
  -o mgxs_sph_applied.h5
```

`augment-sph` remains available when a downstream `L_MACROLIB` workflow
explicitly consumes `GROUP/*/NSPH`; it attaches records without rewriting cross
sections and does not satisfy the standard applied-SPH production contract.

When `external_reference_flux.h5` also contains
`/openmc_volume_flux_std_dev`, the scaffold copies and audits it the same way.

## Existing Lower-Level API

If a script already has a loaded OpenMC `mgxs.Library` object, it can call the
package API directly:

```python
from openmc2donjon import export_openmc_mgxs_library

export_openmc_mgxs_library(library, "mgxs_library.h5")
```

The pickle mode remains available for small local debugging fixtures:

```sh
openmc2donjon-export library.pkl -o mgxs_library.h5
```
