# OpenMC-Side SPH Update-Table Example

This deterministic example exercises the `make-sph-update-table` path:

```text
OpenMC CE reference flux
OpenMC MG macro flux
previous SPH table
  -> make-sph-update-table
  -> make-sph-sidecar
  -> augment-sph
  -> convert to L_MACROLIB
```

It is a deterministic mechanics fixture: the input matrices are already
mapped to a shared `(comparison domain, energy group)` order, and the script
explicitly selects the diagnostic `flux` target and `none` normalization to
keep its reference numbers stable. The standard physical workflow instead uses
a detailed heterogeneous CE geometry, a homogenized MG coarse geometry, and a
rate-preserving, power-normalized update before Converter writes the
DRAGON/DONJON object.

This is not a DONJON feedback loop. DONJON is not run by this example, and
no DONJON flux is fed back into the SPH update.

Run it from the repository root:

```sh
examples/openmc_sph_update_table_example/run_smoke.sh
```

The default run directory is
`/private/tmp/openmc2donjon_openmc_sph_update_table_example`; set `RUN_DIR`
to override it.
