# Converter Product Model

## Product boundary

`openmc2donjon` is a general OpenMC MGXS to DRAGON/DONJON Converter. It does
not define a universal reactor, component count, colorset geometry,
equivalence method, or full-core loading map.

The stable operation is:

```text
one validated MGXS HDF5
  -> Converter contract and preflight
  -> L_MULTICOMPO or L_MACROLIB
  -> output-specific, hash-linked Converter receipt
```

The built-in ASCII writer is the default. PyGan/LCM is an optional writer and
semantic cross-validation backend for the same conversion.

## Product architecture

- **Converter is the mandatory core.** Every formal OpenMC-to-DRAGON/DONJON
  handoff passes through the same preflight, writer, and hash-linked receipt.
- **OpenMC MGXS is optional input preparation.** Use it when the converter-ready
  MGXS HDF5 does not already exist.
- **SPH is optional but first-class physical equivalence.** The standard route
  holds a heterogeneous OpenMC CE calculation fixed as the fine reference and
  iterates a homogenized OpenMC MG calculation with the rate-preserving SPH
  update. The geometries are intentionally different, but the energy groups,
  state, boundary conditions, and conservative comparison-domain mapping must
  agree. That map is complete and non-overlapping and preserves physical volume,
  integrated reference flux, and reaction rates. The corrected HDF5 records
  `sph_applied=true` and then passes through Converter like every other formal
  handoff. SPH never fits a global multiplier to k-effective or another
  requested answer.
- **Projects are optional coordination.** A manifest groups any number of
  repeated component/state handoffs without changing Converter's one-input
  contract.
- **Inspect is an independent read-only tool.** It visualizes compatible MGXS
  data and still exposes the structure and root metadata of other HDF5 files;
  using it does not require Converter.
- **DONJON is the downstream consumer and validation surface.** Its geometry,
  mixture map, solver family (`SN` or `SPN`, for example), boundary conditions,
  and acceptance criteria belong to the user's model or an explicit template.

ASCII and PyGan/LCM are two writer paths inside Converter, not separate product
pipelines.

Generic Converter inputs may preserve explicitly supplied ADF/DF records for a
separately declared workflow. Physical-SPH acceptance permits no ADF
substitution, empirical/global coefficient, fitted k-effective correction, or
numerical exemption such as identity substitution, floors, frozen groups, or
clipping.

An advanced project may explicitly select an external native DRAGON `SPH:`
solver. That route starts from an uncorrected Converter `L_MACROLIB` and receipt;
the project-owned DRAGON deck performs the iteration, while
`validate-native-sph` audits the resulting artifacts. It is a supported
external-solver contract, not the standard built-in SPH operator.

## Manifest-driven projects

A durable project may add `openmc2donjon.project.json` at its root. The
manifest, rather than frontend code, declares:

- any positive number of components;
- each component's label, input HDF5, output object, and required/optional role;
- the input physics contract (`converter-hdf5`, the standard OpenMC CE/MG
  `physical-sph`, or the advanced `native-sph` contract for a Converter
  reference plus externally solved DRAGON SPH);
- optional evidence paths and output identity;
- the downstream consumer and any result artifacts the project wants tracked.
- an optional project-owned acceptance decision whose model-specific criteria
  link to evidence files and hashes.

For a `native-sph` component, `receipt` is always the production Converter
receipt for the uncorrected reference MACROLIB. `physics_summary` is a separate
field for the native-SPH validation summary. The validator receives and checks
both paths independently; they may never name the same file. Legacy manifests
that used `receipt` for a physics summary remain readable for migration, but
their component stays on HOLD until both fields are declared distinctly.

The one-component starting manifest is
`examples/project_templates/minimal/openmc2donjon.project.json`.

## Product invariants

- Converter never invents missing project components or physics evidence.
- Every accepted output is cryptographically linked to its exact input.
- Changing the input or output after conversion invalidates the receipt.
- The manifest can strengthen an input contract but cannot bypass Converter
  preflight.
- Solver results are validation observables; Converter never fits a correction
  coefficient to a requested answer.
- A completed downstream solve is not automatically a physics acceptance.
  Closure criteria belong to the project and its independent reference.
- A project is accepted only when its explicit decision says accepted, every
  declared criterion passes, required Converter outputs are ready, and all
  referenced evidence exists with matching hashes.

## IRENA-30 example templates

`examples/project_templates/irena30_fullcore/openmc2donjon.project.json` is an
advanced native-DRAGON research template and deliberately starts on HOLD. It
represents one position-resolved full-core handoff, not five reusable material
components, and is not an accepted IRENA result or the standard SPH route.

A standard physical full-core candidate must model all 91 heterogeneous
assemblies in OpenMC CE and may pool tallies during transport on 21 exact
global D3 symmetry orbits (or keep 91 independent domains). Those fine
comparison domains map conservatively onto the homogenized OpenMC MG coarse
geometry; the rate-preserving update must converge before the corrected HDF5
records `sph_applied=true` and passes through Converter. No IRENA full-core SPH
result is currently accepted.
Promotion would additionally require k-effective, leakage, 91-position power,
statistical-quality, and every numerical-convergence gate to pass together.

`examples/project_templates/irena30/openmc2donjon.project.json` is the withdrawn
five-colorset diagnostic. Its five material labels (INT, EXT, CSD, DSDF, PNL)
and 91-position map are benchmark data, never frontend or Converter defaults,
and its acceptance decision is permanently rejected.

That older component-library path remains available only for diagnostic
inspection.
`assemble-component-library` selects a record from each qualified component
result and `expand-component-library` copies records onto an explicit map.
Copying changes no cross section or NSPH value, but it also cannot create the
missing global leakage environment. It must not be treated as IRENA full-core
physics acceptance merely because a downstream solve terminates.

Those five components and 91 positions are facts about that example only. A
different user can declare one component, a different colorset family, many
state-dependent libraries, no SPH, a different solver, or an external
consumer without changing Converter itself.
