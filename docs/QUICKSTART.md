# Quickstart

This page is the shortest path from a fresh checkout to a working
OpenMC-to-DONJON conversion.

## Install

From the repository root:

```sh
python -m pip install -e .
```

Or run directly from source:

```sh
export PYTHONPATH=src
```

## Check The User Entry Point

Check the local Python/package environment:

```sh
openmc2donjon doctor
```

Run the tiny recipe/statepoint smoke:

```sh
bash scripts/run_recipe_export_smoke.sh
```

This checks:

- environment doctor;
- recipe dry-run metadata preflight;
- one-step dry-run conversion-plan preflight;
- recipe/statepoint export to the HDF5 handoff contract;
- HDF5 inventory inspect;
- HDF5 preflight;
- HDF5 baseline diff;
- managed run directory and bundle manifest;
- `L_MULTICOMPO` write/readback;
- root `L_MACROLIB` write/readback;
- one-command `openmc2donjon-from-openmc` conversion.

On a machine with OpenMC and continuous-energy data configured, run the minimal
production-style case:

```sh
bash scripts/run_production_minicase_smoke.sh
```

## Open The Local Web UI

The CLI remains the authoritative execution path, but the localhost Web UI is
the fastest way to inspect a handoff, dry-run conversion gates, preview command
lines, and build OpenMC-side ADF/SPH sidecar commands:

```sh
python -m pip install -e ".[web]"
openmc2donjon serve

# In another shell:
cd web
npm install
npm run dev
```

Open <http://localhost:3000>. The main surfaces are:

- `Commands` - workflow map plus every CLI command linked to a Web surface;
- `Convert` - direct HDF5 -> MULTICOMPO/MACROLIB dry-run and conversion;
- `Inspect` - HDF5 mixture, energy mesh, scatter, ADF/SPH, and std_dev hints;
- `Equivalence` - copyable ADF/SPH sidecar and augmentation commands;
- `Builder` - copyable command builders for diagnostics and bundles;
- `DONJON` - starter handoff guidance for downstream cards.

For fixture-backed UI work without local files, start the backend with
`openmc2donjon serve --mock`.

## One-Step OpenMC To DONJON

For a real OpenMC case, write a small recipe that builds the case's
`openmc.mgxs.Library`. You can start from
[`examples/openmc_recipe_template/`](../examples/openmc_recipe_template/) or inspect
[`examples/production_minicase/`](../examples/production_minicase/). Then run:

```sh
openmc2donjon-from-openmc \
  --recipe export_recipe.py \
  --dry-run \
  --run-dir runs/case1 \
  --check \
  --strict-dry-run
```

The dry-run output includes a production checklist for MGXS coverage,
transport/STRD readiness, domain-to-mixture mapping, volumes, and `domain_mode`.
Strict dry-run returns non-zero if that checklist still has warnings or
failures.

Before running OpenMC, write the exact MGXS tallies declared by the recipe:

```sh
openmc2donjon-export --recipe export_recipe.py --write-tallies tallies.xml
```

The recipe need not be beside `main.py`, and `main.py` needs no
openmc2donjon-specific code if it uses this `tallies.xml`. If it replaces
`model.tallies` or exports XML again, integrate the same library tallies or
change the order to `export model XML -> generate tallies.xml -> run OpenMC`.
The recipe geometry/domain IDs must match the statepoint.

Keep the static absorption/scattering/transport estimators paired. The ordinary
policy uses `absorption` + `scatter matrix` + `transport`. For a fast-spectrum
case where `(n,xn)` neutron multiplicity matters, also tally
`reduced absorption`, explicitly select `consistent nu-scatter matrix`, and use
`nu-transport` instead of ordinary `transport`, as documented in
[OpenMC Export Workflow](OPENMC_EXPORT_WORKFLOW.md). No manual absorption or
scatter-matrix edit is needed; the static handoff does not generate separate
`N2N`/`N3N` depletion records.

```sh
openmc2donjon-from-openmc \
  --recipe export_recipe.py \
  --statepoint statepoint.120.h5 \
  --run-dir runs/case1 \
  --check
```

The run directory contains the HDF5 handoff, DONJON ASCII output, summary JSON,
check summary when `--check` is enabled, and `manifest.json`. Existing managed
files are refused unless `--force-run-dir` is set. Add side artifacts to the
same manifest with repeatable `--extra-artifact LABEL=PATH` options:

```sh
openmc2donjon-from-openmc \
  --recipe export_recipe.py \
  --statepoint statepoint.120.h5 \
  --run-dir runs/case1 \
  --extra-artifact surface-flux=openmc_surface_flux.h5 \
  --extra-artifact low-order-driver=low_order_driver.h5 \
  --extra-artifact homogeneous-face-flux=homogeneous_face_flux.h5 \
  --check
```

For production ADF/DF handoff, the same entry point can build the flux-ratio
ADF sidecar from OpenMC surface flux plus a low-order driver and inject it
before conversion:

```sh
openmc2donjon-from-openmc \
  --recipe export_recipe.py \
  --statepoint statepoint.120.h5 \
  --run-dir runs/case1 \
  --build-flux-ratio-adf \
  --export-surface-flux \
  --surface-flux-mu-edges 0.0,0.25,0.5,0.75,1.0 \
  --low-order-raw-driver raw_low_order_driver.h5 \
  --low-order-net-current-sign-convention auto \
  --adf-faces FD_XMIN,FD_XMAX,FD_YMIN,FD_YMAX \
  --adf-face-widths 4.0 \
  --require-volume \
  --require-transport-dataset
```

When the homogeneous face flux has already been computed, pass it directly as
the denominator and skip low-order reconstruction:

```sh
openmc2donjon-from-openmc \
  --recipe export_recipe.py \
  --statepoint statepoint.120.h5 \
  --run-dir runs/case1 \
  --build-flux-ratio-adf \
  --export-surface-flux \
  --surface-flux-mu-edges 0.0,0.25,0.5,0.75,1.0 \
  --homogeneous-face-flux homogeneous_face_flux.h5 \
  --adf-faces FD_XMIN,FD_XMAX,FD_YMIN,FD_YMAX \
  --require-volume \
  --require-transport-dataset
```

Use `FILE::/dataset/path` for either external face-flux argument when a source
file uses case-specific dataset names.
The managed `--build-flux-ratio-adf` workflow writes
`face_flux_check_summary.json` and bundles it in `manifest.json`; standalone
users can run `openmc2donjon check-face-flux` before `make-adf-sidecar` for the
same numerator/denominator contract check.

A complete deterministic external low-order handoff example is available here:

```sh
bash examples/external_low_order_handoff/run_smoke.sh
```

When the external solver already computes the homogeneous face flux directly,
start from the adapter template instead of reconstructing the denominator from
volume flux and net current:

```sh
bash examples/external_face_flux_adapter/run_smoke.sh
```

For standard physical SPH, keep a heterogeneous OpenMC CE model as the fixed
fine reference and solve a homogenized OpenMC MG coarse model. The geometries
are different, but the CE tally bins and MG transport groups, state, boundary
conditions, and conservative fine-to-coarse comparison-domain map must agree.
Their native OpenMC domain IDs need not be equal. Set `CE_DOMAIN_IDS` and
`MG_DOMAIN_IDS` to the comma-separated native tally-bin IDs that correspond,
in order, to the canonical `mixture_names` in `mgxs_library.h5`. Export aligned
volume fluxes from both statepoints; the complete, non-overlapping map must
preserve physical volume, integrated reference flux, and reaction rates.
Also set `CE_FLUX_REL_SIGMA_LIMIT` and `MG_FLUX_REL_SIGMA_LIMIT` from the
project's predeclared uncertainty budget; the product does not invent universal
values. Compute the rate-preserving sidecar and iterate the OpenMC MG solve:

```sh
openmc2donjon export-volume-flux ce_statepoint.h5 \
  --mgxs mgxs_library.h5 \
  --tally-name openmc_ce_volume_flux \
  --dataset-name openmc_volume_flux \
  --source-domain-ids "${CE_DOMAIN_IDS}" \
  -o openmc_ce_flux.h5

openmc2donjon export-volume-flux mg_statepoint.h5 \
  --mgxs mgxs_library.h5 \
  --tally-name openmc_mg_volume_flux \
  --dataset-name openmc_mg_flux \
  --source-domain-ids "${MG_DOMAIN_IDS}" \
  -o openmc_mg_flux.h5

openmc2donjon make-openmc-sph-sidecar mgxs_library.h5 \
  -o openmc_sph_sidecar.h5 \
  --reference-flux openmc_ce_flux.h5::openmc_volume_flux \
  --mg-flux openmc_mg_flux.h5::openmc_mg_flux \
  --table-output openmc_sph.csv \
  --sph-target rate \
  --flux-normalization power \
  --require-reference-flux-std-dev \
  --max-reference-flux-std-dev-rel "${CE_FLUX_REL_SIGMA_LIMIT}" \
  --require-mg-flux-std-dev \
  --max-mg-flux-std-dev-rel "${MG_FLUX_REL_SIGMA_LIMIT}"
```

Apply `XS / NSPH` to the OpenMC-native `setN` library, rerun OpenMC MG, and
recompute the sidecar with `--previous-sph` until converged:

```sh
openmc2donjon apply-sph mg_case/mgxs_unapplied.h5 \
  --input-format openmc-mgxs \
  --sph-source openmc_sph_sidecar.h5 \
  -o mg_case/mgxs.h5
```

This OpenMC-native `setN` application is an iteration artifact, so it is
recorded as `openmc-mgxs-intermediate-unbound`; it is not a Converter handoff
and cannot pass `--require-physical-sph`. The final application below uses the
exact Converter-layout HDF5 that the sidecar hash binds.

Apply the converged factors to the Converter-facing HDF5, then use Converter as
the mandatory formal handoff boundary. `apply-sph` records
`sph_applied=true`:

```sh
openmc2donjon apply-sph mgxs_library.h5 \
  --sph-source openmc_sph_sidecar.h5 \
  -o mgxs_sph_applied.h5

openmc2donjon mgxs_sph_applied.h5 -o out.mcompo.txt \
  --production --require-physical-sph
```

The small portable smoke for this route is:

```sh
bash examples/openmc_sph_sidecar_minicase/run_smoke.sh
```

To keep explicit paths instead of using a managed run directory:

```sh
openmc2donjon-from-openmc \
  --recipe export_recipe.py \
  --statepoint statepoint.120.h5 \
  --keep-hdf5 mgxs_library.h5 \
  -o out.mcompo.txt \
  --summary-json run_summary.json
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

For root `L_MACROLIB` output:

```sh
openmc2donjon-from-openmc \
  --recipe export_recipe.py \
  --statepoint statepoint.120.h5 \
  --format macrolib \
  -o out.macrolib.txt \
  --summary-json run_summary.json
```

## Two-Step Workflow

If you want to inspect or archive the HDF5 before conversion:

```sh
openmc2donjon-export \
  --recipe export_recipe.py \
  --statepoint statepoint.120.h5 \
  -o mgxs_library.h5

openmc2donjon inspect mgxs_library.h5
openmc2donjon mgxs_library.h5 -o out.mcompo.txt --check
```

To confirm a regenerated handoff matches a baseline:

```sh
openmc2donjon diff accepted_mgxs.h5 mgxs_library.h5
```

To inject computed ADF/DF values before conversion:

```sh
openmc2donjon make-adf-sidecar mgxs_library.h5 \
  -o adf_sidecar.h5 \
  --mode unity

openmc2donjon augment-adf mgxs_library.h5 \
  --adf-source adf_sidecar.h5 \
  -o mgxs_with_adf.h5 \
  --faces FD_XMIN,FD_XMAX,FD_YMIN,FD_YMAX
```

`make-adf-sidecar --mode unity` writes identity ADF values marked
`adf_real=false`; use it to verify plumbing before replacing the sidecar with
case-specific physics ADF/DF values.

For a physics sidecar, provide heterogeneous and homogeneous face-flux HDF5
datasets:

```sh
openmc2donjon export-surface-flux statepoint.120.h5 \
  --mgxs mgxs_library.h5 \
  -o openmc_surface_flux.h5 \
  --tally-name openmc2donjon_surface_current_mu \
  --mesh-shape 1,2 \
  --mu-edges 0.0,0.25,0.5,0.75,1.0 \
  --face-area 4.0

openmc2donjon make-low-order-driver mgxs_library.h5 \
  -o low_order_driver.h5 \
  --raw-driver raw_low_order_driver.h5 \
  --net-current-sign-convention auto \
  --faces FD_XMIN,FD_XMAX,FD_YMIN,FD_YMAX

openmc2donjon check-low-order-driver mgxs_library.h5 low_order_driver.h5 \
  --faces FD_XMIN,FD_XMAX,FD_YMIN,FD_YMAX \
  --face-widths 4.0

openmc2donjon make-homogeneous-face-flux mgxs_library.h5 \
  -o homogeneous_face_flux.h5 \
  --volume-flux low_order_driver.h5 \
  --net-current low_order_driver.h5 \
  --faces FD_XMIN,FD_XMAX,FD_YMIN,FD_YMAX \
  --face-widths 4.0

openmc2donjon make-adf-sidecar mgxs_library.h5 \
  -o adf_sidecar.h5 \
  --mode flux-ratio \
  --surface-flux openmc_surface_flux.h5 \
  --homogeneous-face-flux homogeneous_face_flux.h5 \
  --faces FD_XMIN,FD_XMAX,FD_YMIN,FD_YMAX
```

Or inject the sidecar during the one-step OpenMC export:

```sh
openmc2donjon-from-openmc \
  --recipe export_recipe.py \
  --statepoint statepoint.120.h5 \
  --run-dir runs/case1 \
  --adf-source adf_sidecar.h5 \
  --adf-faces FD_XMIN,FD_XMAX,FD_YMIN,FD_YMAX \
  --check --require-adf
```

Run preflight on the HDF5:

```sh
openmc2donjon check mgxs_library.h5
```

## C5G7 Accepted Check

Run the portable converter-side acceptance check:

```sh
bash scripts/release_check.sh --skip-tests
```

On a machine with the local DRAGON/DONJON checkout and staged data, run the full
DONJON-side check:

```sh
bash scripts/release_check.sh --run-donjon
```

## Next Documents

- [OpenMC export workflow](OPENMC_EXPORT_WORKFLOW.md)
- [Fast-spectrum workflow](FAST_SPECTRUM_WORKFLOW.md)
- [HDF5 input contract](HDF5_INPUT_CONTRACT.md)
- [From-OpenMC summary JSON](FROM_OPENMC_SUMMARY_SCHEMA.md)
- [Validation summary](VALIDATION.md)
