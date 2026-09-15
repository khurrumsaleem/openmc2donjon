# Handoff Note

For a current repo-wide status snapshot, start with
[Handoff Snapshot](HANDOFF_SNAPSHOT.md). This note keeps the accepted C5G7
handoff baseline details.

## Stable Tag

```text
v0.1.2-openmc-workflow
```

This tag is the current internal handoff point for `openmc2donjon`.

## What Is Accepted

The accepted validation line is C5G7 assembly-wise homogenization:

```text
OpenMC statepoint
  -> OpenMC mgxs.Library exporter
  -> HDF5 input contract
  -> L_MULTICOMPO / L_MACROLIB
  -> DONJON assembly-wise solves
```

The accepted HDF5 snapshot is:

```text
examples/donjon_openmc2donjon/c5g7_assembly_p1_adf_production.h5
```

It contains:

- 9 assembly-wise mixtures;
- 7 energy groups;
- P1 scattering;
- production ADF payload copied from the OpenMC/DONJON face-flux workflow.

## Locked Results

| Case | k-effective |
| --- | ---: |
| OpenMC reference | `1.18798` |
| DONJON diffusion | `1.1896194220` |
| DONJON SPN3 | `1.1912802458` |
| DONJON SPN3 with first-order scattering retained | `1.1912822723` |
| DONJON 2-group ADF smoke | `1.18533289` |
| DONJON 2-group NODF smoke | `1.20179343` |

The accepted snapshot now declares OpenMC's ordinary
`consistent scatter matrix` paired with `absorption`. An audit of the saved
source statepoint found its ordinary and nu-weighted P0 scattering
reaction-rate tallies identical in all 441 bins, so this thermal C5G7 fixture
does not exercise an `(n,xn)` multiplicity handoff. Its numerical HDF5 datasets
remain unchanged; only the explicit scattering-contract metadata was corrected.

The saved source statepoint predates the current transport-complete C5G7 recipe.
It does not contain the complete tally set now required, including an explicit
OpenMC `TransportXS`, and therefore cannot directly reproduce the accepted HDF5
through the current exporter. Generate `tallies.xml` with the current recipe,
rerun OpenMC, and use that new statepoint before claiming exporter parity. No P1
row-sum reconstruction is an acceptable substitute.

## Reproduction

Portable release check:

```sh
bash scripts/release_check.sh
```

Full local acceptance with DONJON decks:

```sh
bash scripts/release_check.sh --run-donjon
```

The full check covers package tests, CLI smoke, C5G7 converter readback,
accepted baseline validation, C5G7 DONJON face-flux regeneration when local
dumps are available, and DONJON locked decks. Direct parity regeneration from
the saved C5G7 statepoint is intentionally not claimed until a new statepoint
has been produced with the current transport-complete recipe.

## Scope

- One state point by default.
- Experimental one-dimensional `BURN` history serialization exists, but no
  burnup/history/branch parameter axis is accepted as a physics validation line
  yet. A tiny DONJON `NCR:` smoke verifies the serializer's `BURN` selection
  plumbing.
- General multi-axis branch libraries are rejected by preflight and converter
  code rather than being silently ignored.
- Spatial mapping is one OpenMC MGXS domain or subdomain to one DONJON mixture.
- Hex-domain conversion/modeling support exists, but no accepted hex benchmark
  is included in this release.
- OpenMC remains the source of homogenized MGXS data; this package is the
  HDF5-to-DONJON handoff layer.
