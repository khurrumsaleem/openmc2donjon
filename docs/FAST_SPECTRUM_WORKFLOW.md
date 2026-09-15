# IRENA Fast-Spectrum Example Workflow

This document defines the strict physics contract selected by the IRENA-30
full-core example manifest. It is not a universal Converter workflow.
Historical scripts that use local colorsets, post-transport record reuse,
zero-flux filling, identity substitutions, flux floors, frozen groups, clipping,
eigenvalue fitting, or a global multiplier remain useful as research records,
but they do not produce an accepted IRENA full-core result.

The current qualification unit is one position-resolved full core. The fine
transport either keeps 91 independent domains or pools tallies on 21 exact
global D3 symmetry orbits while particles are being transported:

```text
heterogeneous OpenMC CE 91-position fine full core
  -> 91 independent domains or 21 exact D3 transport-time tally pools
  -> conservative fine-to-coarse comparison-domain map
  -> homogenized OpenMC MG full-core coarse solve
  -> rate-preserving SPH update, XS / NSPH, and MG rerun to convergence
  -> corrected MGXS HDF5
  -> Converter + hash-linked receipt
  -> DONJON k-effective, leakage, and 91-position power verification
```

Historical IRENA local studies declare colorsets such as `int_ext`, `ext_int`,
`csd_int`, `dsdf_int`, and `pnl_ext`. These labels and counts are template
facts, not Converter defaults. They are not the current full-core production
mapping because identical material labels can occupy different global leakage
environments.

## 1. Strict full-core SPH contract

An accepted IRENA full-core handoff must satisfy all of the following:

- The fine reference contains all 91 physical OpenMC CE positions with exact
  integrated flux, reaction rates, energy coverage, uncertainty, and boundary
  leakage evidence.
- The homogenized OpenMC MG geometry is a different coarse representation. It
  uses the CE tally group boundaries and the same state and boundary
  conditions as the CE problem,
  and its comparison domains are connected to the fine model by a complete,
  non-overlapping, volume/rate-conservative map.
- The rate-preserving OpenMC CE/MG update converges on those comparison domains
  before any formal handoff is accepted.
- Converter validates the corrected HDF5, writes the declared
  `L_MULTICOMPO` or `L_MACROLIB`, and records the hash-linked receipt.
- Zero or unusable bins are rejected. Identity substitution, floors, frozen
  groups, clipping, and ADF are absent.
- Converter rate balance and final DONJON eigenvalue both pass the predeclared
  OpenMC-uncertainty gate.
- Full-core k-effective, leakage, and power are validation observables only; no
  empirical scalar or fitted global coefficient participates in factor
  generation.

If statistics are insufficient to evaluate a bin, the remedy is better tally
design, a physically justified energy structure, or more histories. A numerical
exception is not a physical SPH result.

The local workspace currently has no accepted IRENA full-core physical closure.
Earlier PNL/EXT and INT/EXT summaries had converged SPH fixed points but
unconverged final transport solves, so they are retained only as negative
evidence.

## 2. Static `(n,xn)` cross-section contract

IRENA uses the multiplicity-weighted fast-spectrum policy. Its OpenMC recipe
must tally the following matched set; ordinary absorption is retained for its
reaction-rate observable, while reduced absorption is the balance partner of
the consistent nu-scatter matrix:

```python
MGXS_TYPES = [
    "total",
    "absorption",
    "reduced absorption",
    "fission",
    "kappa-fission",
    "nu-fission",
    "chi",
    "consistent nu-scatter matrix",
    "nu-transport",
]


def scatter_mgxs_type():
    return "consistent nu-scatter matrix"
```

Do not manually subtract `(n,2n)` or `(n,3n)` rates from absorption and do not
edit scattering matrices in `main.py`, the recipe, or an HDF5 postprocessor.
OpenMC constructs the paired `reduced absorption`,
`consistent nu-scatter matrix`, and `nu-transport` MGXS estimators. The
exporter preserves the set, and Converter serializes total plus the selected
multiplicity-weighted scattering records.

In the resulting static DONJON object, net absorption is implicit in `NTOT0`
minus the outgoing P0 `SCAT` row. Consequently `(n,xn)` neutron multiplication
is retained by the reduced-absorption/nu-scatter balance; Converter does not
need a fitted correction or separate matrix rewrite. This is a static
macroscopic transport handoff: it does **not** emit separate `N2N` or `N3N`
depletion-reaction records. A depletion workflow requiring those reaction
channels needs a separate, explicitly validated depletion-data contract.

## 3. Converter is the controlled boundary

The standard SPH iteration happens between the fixed OpenMC CE reference and
the homogenized OpenMC MG coarse solve. After convergence, `apply-sph` folds
the factors into the Converter-facing HDF5. Converter remains mandatory: it
runs the production contract, writes the downstream object, and receipts its
exact input and output. Neither a sidecar nor a terminating transport solve can
bypass that boundary or establish physics acceptance.

An advanced project may instead choose an external native DRAGON `SPH:` solve.
That project-specific route must first use Converter to write and receipt an
uncorrected reference `L_MACROLIB`; the external solver performs the iteration,
and `validate-native-sph` audits the deck and artifacts. It is not the standard
IRENA or product operator, and no native-DRAGON IRENA full-core result is
accepted.

## 4. Full-core DONJON model

The IRENA full-core solve is a coarse transport model over 91 physical core
positions (52 are fuel), not a 91-mixture fit and not 91 fuel assemblies. The
fine reference keeps all 91 heterogeneous assemblies. It may retain 91
independent homogenized domains or pool tallies during OpenMC transport on the
21 exact global D3 symmetry orbits. The older five-material and 13-local-
signature maps are diagnostic only because they overmerge distinct global
environments. The downstream solver may be SN or SPN; a special `SN8` choice
is not part of the general physical contract.

Full-core k-effective, leakage, and power shape are validation observables only.
They may reveal that the coarse model is inadequate, but they may never be
used to tune a global SPH coefficient.

## 5. Historical evidence

- `examples/irena30_sph_stage2_csd` records the earlier seven-assembly
  colorset experiments. Its identity/floor/freeze/clip prescriptions are
  archived and production-rejected; its reproduction runners are guarded and
  their summaries are permanently marked withdrawn diagnostics.
- `examples/irena30_sph_stage3_fullcore` records a rejected OpenMC-MG full-core
  SPH research line. Its sparse/fill/floor/frozen-group artifacts are not an
  accepted full-core input.
- `examples/irena30_native_fullcore` records the advanced external-DRAGON
  candidate and validator work. It remains useful implementation evidence, but
  it is not the standard route and has no accepted IRENA full-core result.
- `examples/irena30_zrefl_hex` remains useful converter/DONJON benchmark
  evidence. Its older 91-position representation does not prove an accepted
  five-component physical SPH model.
- C5G7 ADF validation remains a separate capability. ADF is not part of the
  IRENA product route described here.
