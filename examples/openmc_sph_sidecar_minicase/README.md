# OpenMC CE/MG SPH Sidecar Minicase

This portable smoke exercises the sidecar and conversion mechanics without
requiring a real OpenMC run:

```text
OpenMC CE reference flux
  + OpenMC MG macro flux on the same regions
  -> make-openmc-sph-sidecar
  -> augment-sph
  -> openmc2donjon conversion to L_MULTICOMPO / L_MACROLIB ASCII
```

The fixture is deliberately tiny: two output regions and two energy groups.
`make_inputs.py` writes:

- `mgxs_library.h5` - converter-facing MGXS handoff
- `openmc_ce_flux.h5` - stand-in for the OpenMC continuous-energy reference
  region/group flux
- `openmc_mg_flux.h5` - stand-in for the OpenMC multi-group macro flux already
  mapped to the same comparison-domain and group order
- `reference_expected.h5` - expected SPH factors for the smoke validator

Run it from the repository root:

```sh
bash examples/openmc_sph_sidecar_minicase/run_smoke.sh
```

The example intentionally does **not** use a DONJON feedback loop. It also
explicitly selects the diagnostic `flux` target and `none` normalization, so it proves
the handoff mechanics for OpenMC-side SPH factors: compute factors from CE/MG
flux comparison, augment the HDF5 with them as `NSPH`, and verify that both ASCII writer
formats carry the factors through. It is not physical-equivalence acceptance;
the standard route uses a fine heterogeneous CE geometry, a homogenized MG
coarse geometry, and a rate-preserving fixed point.
