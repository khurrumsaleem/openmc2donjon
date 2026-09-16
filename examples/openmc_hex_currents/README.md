# Native hexagonal face-current diagnostic

Experimental, opt-in Python support for explicit hexagonal current domains in
standard OpenMC. It does **not** require a custom HexMesh extension, change the
Converter's input requirements, or implement SPH/ADF. It is not yet a Web UI
workflow or a general automatic adapter for arbitrary OpenMC models.

Each unique position/axial-layer cell has six side faces plus bottom and top.
Shared global planes, explicit neighbor identities, outward normals and face
areas define the mapping. Surface filters combined with the cell before/after
crossing restrict statistics to each finite face; a bare infinite-plane tally
would not be sufficient.

## Run the public verification

Use matching OpenMC Python and executable installations, NumPy, h5py, and an
installed checkout of openmc2donjon:

```sh
python examples/openmc_hex_currents/validate.py \
  --work-dir /path/to/new/hex-current-check \
  --openmc /path/to/openmc --threads 2
```

The directory must be new. `--build-only` writes inputs without running.
OpenMC may require a configured `OPENMC_CROSS_SECTIONS` XML even though this
verification uses void geometry and no nuclear interaction data.

The recorded native runs were tested with OpenMC
`0.15.4.dev160+ge542b2f03`. Compatibility with older releases, especially their
lattice-crossing current behavior, has not been established by these runs.

Two adjacent hexagonal columns, each split into two axial layers, exercise:

1. **Directed crossings:** 64 individually energy-labeled sources cover all
   four nodes, eight faces, and both directions. An independent analytic ray
   intersection calculation verifies the entire scored pattern, including
   zero entries outside the correct face/layer. Nested fine cells test parent
   cell filters. Shared-face currents and unit global leakage are checked.
2. **Opposing random sources:** two equally weighted beams cross the same face
   in opposite directions. Their partial currents are anticorrelated. The
   net-current standard error must be `sqrt(2)` times the incorrect independent
   quadrature estimate for this specific test. Upper-layer bins remain zero.

Each test uses 12 batches × 8,000 particles, one MPI rank and the requested
thread count. It saves every active-batch statepoint. These are geometry and
statistics-method tests, **not a reactor-physics benchmark**.

## Python integration

```python
from openmc2donjon.openmc_hex_currents import (
    build_hex_current_surfaces,
    bind_hex_current_domains,
)

registry = build_hex_current_surfaces(
    sites=[(0, 0), (0, 1)], z_edges=[0.0, 1.0, 2.0], pitch=2.0,
)
# Build actual unique OpenMC cells using registry.node_region(site, layer).
# Retain the intended fine universe as each cell's fill.
plan = bind_hex_current_domains(
    registry, domains, energy_bounds=energy_bounds, geometry=geometry,
)
# Add plan.tallies to the calculation without approximate tally merging.
# Save plan.manifest with the audited inputs and execution receipt.
```

See `validate.py` for a complete runnable example, including provenance.
The registry's surfaces must be used in the **actual geometry**, not simply
created as detached scoring planes. Adding the tally calls alone to an existing
HexLattice does not expose its implicit tile faces as Surface IDs.

Supported first-release scope:

- One uniform y-oriented hexagonal grid (matching x-oriented hexagonal tiles).
- Unique, uncut, explicit global-coordinate cell domains, optionally containing
  translated/rotated fine fills. The current domains themselves must not be
  repeated through a lattice or hidden under transformed ancestors.
- Transmission faces only. An exterior void collar in the example keeps the
  physical vacuum boundary separate. Reflective, periodic and coincident vacuum
  current faces need dedicated handling and are deliberately not supported.
- Neutron-only analog current statistics. Fixed-source strength must sum to one;
  eigenvalue currents are per starting-source neutron, not power-normalized.

Ancestor clipping, duplicated/coincident face definitions and changed serialized
filter bins are rejected. Model adaptation must independently verify that the
fine material geometry and physical boundary conditions were preserved.

## Export and uncertainty

```python
from openmc2donjon.openmc_hex_current_export import export_hex_currents

export_hex_currents(
    active_batch_statepoints,
    "current_manifest.json",
    "hex_currents.h5",
    summary_json="current_export_summary.json",
)
```

The manifest must bind geometry/materials/settings/tallies XML by SHA-256 and
reference a successful execution receipt. The receipt binds the full manifest
definition and each statepoint. The example runner implements this contract.
Existing output files are not overwritten.

Output arrays use `[node, face, group]` with **high-to-low energy group order**:

- `currents/out_mean`, `out_std_dev`: positive outward partial current.
- `currents/in_mean`, `in_std_dev`: positive inward partial current.
- `currents/net_mean`: outward minus inward.
- `currents/net_std_dev` and `out_in_covariance_of_mean`: only when all active
  batch checkpoints, starting with realization 1, are supplied.
- `current_density/…`: the same quantities divided by face area (area squared
  for covariance). This is current density, **not scalar surface flux**.
- Face IDs, names, areas, normals, node identities and energy bounds accompany
  the arrays. Bounds are ascending; `group_bounds_eV` explicitly maps each
  high-to-low array row.

Consecutive cumulative tally-sum differences reconstruct paired batch samples.
Their covariance is retained in the net standard error and checked against the
native final marginal moments. With only the final statepoint, partial-current
errors and the net mean can be exported, but the net error is explicitly absent.
It is never invented by treating incoming/outgoing estimates as independent.

This is a separate diagnostic sidecar, **not a Converter-ready MGXS HDF5, SPH
coefficient file, ADF, or DONJON object**. Zero-scored bins remain visible; they
do not establish precision. Batch standard errors do not establish fission-source
convergence or eliminate inter-batch correlations in eigenvalue problems.
