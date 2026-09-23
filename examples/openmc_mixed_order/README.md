# Native mixed P1/P3 scattering orders

Different homogenization regions can score different Legendre orders. For
example, explicitly selected lead regions may score P3 while the remaining
regions score P1. Selection is by **domain ID**, not an automatic material-name
heuristic; a homogenized region can contain several materials.

## Integrate into an OpenMC recipe

Build a fresh native cell-domain MGXS library with `correction=None`, then set
orders **before accessing any MGXS tallies**:

```python
from openmc2donjon.openmc_mixed_order import (
    set_domain_scatter_orders,
    exact_cell_mgxs_tallies,
)

library.build_library()
orders = {cell.id: 1 for cell in library.domains}
orders[lead_cell.id] = 3
set_domain_scatter_orders(library, orders)
model.tallies = exact_cell_mgxs_tallies(library)
```

The helper requires macroscopic cell-domain MGXS, unique domain IDs, and an
explicit order for every domain. It groups identical native tally definitions
with exact cell sets, while keeping different Legendre filters separate. The
export recipe must reconstruct **the same per-domain orders** before loading
the statepoint. Existing P1 statepoints cannot provide previously unscored P2/P3;
rerun OpenMC with the new tallies. Keep any independent model tallies as well.

The scattering/removal convention is independent of this choice. Ordinary
scattering still pairs with absorption and transport; the example uses
`consistent nu-scatter matrix`, `reduced absorption` and `nu-transport`.

## What Converter stores

- The HDF5 and DRAGON/DONJON object use a common **maximum** moment dimension.
  A mixed P1/P3 input therefore reports root `legendre_order=3`.
- Each mixture records `source_legendre_order` and `scatter_padding`.
  P1 regions carry explicit zero truncation for P2/P3; P3 regions retain the
  measured moments, including their signs and statistical uncertainties.
- Padded zeros, including their zero standard deviations, are **not measured
  physical zeros or evidence of zero statistical uncertainty**.
- Preflight logs/JSON distinguish the source orders from the storage maximum.
  Inconsistent declarations and nonzero padded tails are rejected. Older HDF5
  files without these attributes remain readable; their local scoring order
  cannot be inferred from zero values.
- Component collapse retains the maximum contributing order and records the
  contributing source orders, rather than copying the first region's order.
  If any contributor has an unknown original order, the combined original order
  remains unknown. `collapsed_source_legendre_orders` records `-1` for unknown
  immediate contributors; the order/padding declaration is omitted. Preflight
  lists unknown regions explicitly, including after repeated collapse.
- `fill-zero-flux` respects each declared local truncation: a P3 donor cannot
  insert P2/P3 into a P1 region. The donor order and applied order are recorded
  separately; invalid declarations or truncated tails fail before mutation.
  This operation is only for the matching OpenMC MG material library that the
  transport run consumed; it cannot repair unscored CE data or create P3 scores.

This reduces native high-order tally storage compared with scoring P3 everywhere.
It does not promise a proportional runtime gain or a smaller final DONJON object.
The downstream solver's angular approximation and scattering-order settings are
separate choices; verify them in your own calculation.

## Run the transfer diagnostic

Use matching OpenMC Python/executable installations, an installed source checkout
of openmc2donjon, and a CE library containing C12 and Pb208:

```sh
export OPENMC_CROSS_SECTIONS=/path/to/cross_sections.xml
python examples/openmc_mixed_order/validate.py \
  --work-dir /tmp/mixed-pn-check \
  --openmc /path/to/openmc --threads 2 \
  --donjon-root /path/to/dragon-5.1/Donjon
```

Omit `--donjon-root` to check native OpenMC, HDF5 and direct ASCII conversion only.
Use a new/empty work directory. The optional DONJON step uses the local `rdonjon`
wrapper and creates a unique deck under its `data/openmc2donjon/` directory.

The two-cell, two-group fixed-source case uses 80,000 histories. It checks that
the P1 cell has no P3 tally, compares **means and standard deviations** with native
OpenMC MGXS, and compares P0–P3 after native DONJON `NCR:` ingestion of the
`L_MULTICOMPO`. Outputs include the HDF5, both ASCII objects, logs, the DONJON
readback and a hash-indexed `receipt.json`.

This is a transfer/implementation check, **not** a k-effective, source-convergence,
SPH, full-core or production-accuracy validation. Sparse terms may have large
relative uncertainty; warnings remain visible. A successful receipt does not
approve a user's local P1/P3 truncation choices or production statistics.
