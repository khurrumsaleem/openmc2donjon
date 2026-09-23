# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

For narrative release notes (rationale, validation status, migration tips),
see [`RELEASE_NOTES.md`](RELEASE_NOTES.md). This file is a chronological
machine-readable index; cross-reference the release notes for context.

## [Unreleased]

- Fixed SPH scattering-axis resolution to share Converter's inherited/alias
  and unambiguous-shape rules while preserving the stored layout. Legacy
  sidecars without an input hash are now reported as unbound, not verified;
  malformed or mismatched hashes still fail, and physical acceptance is unchanged.
  The web response and result card now retain binding and input-format evidence;
  unbound, native MG iteration, incomplete, and mock results cannot offer the
  final Converter handoff link.
- Made SPH application metadata-aware and failure-safe in Converter and native
  OpenMC-MG layouts, retaining substitution records and uncertainty annotations
  without exposing partially scaled outputs. Component collapse now handles
  explicitly declared nonfission structural zeros without discarding the
  fuel's fission/nu-fission/chi uncertainty; genuine missing statistics still fail.
- Validate flux group/domain labels before both collapse paths, share Converter's
  heat-factor and inverse-velocity aliases, and reject conflicting aliases or
  partially supplied component vectors instead of losing fields. Collapse now
  maps and unions zero-flux substitution indices and preserves per-origin donor
  records through repeated transformations, with explicit contribution semantics.
- Made both collapse paths failure-safe: stage the HDF5, verify and refresh
  existing OpenMC provenance, retain a hash-bound processing history, then
  publish. SPH-bearing/applied and ADF-bearing inputs are explicitly rejected
  rather than silently losing model-specific equivalence. Energy collapse now
  retains vector-XS and chi uncertainties with declared propagation limits.
- Fixed component and energy collapse for moment-last scattering matrices,
  including ambiguous-sized 4-group/P3 inputs with explicit axes. Both paths
  now share Converter's layout rules, normalize means and matching standard
  deviations, and write canonical axis metadata. Invalid layouts or padded
  moments are rejected before opening the destination. Energy collapse retains
  available scattering uncertainty as a fixed-flux-weight conservative L1 bound.
- Added opt-in per-cell Legendre orders for native OpenMC MGXS (for example,
  P3 in selected lead regions and P1 elsewhere), exact grouped tally preparation,
  explicit truncation metadata, and a mixed-order OpenMC/DONJON transfer example.
  Preflight and readers reject inconsistent source-order declarations; component
  collapse preserves the contributing orders and combined fissionability.
  Mixed-order zero-flux filling now preserves the target's local truncation;
  component collapse retains unknown legacy orders instead of inferring them
  from the stored maximum.
- Added an experimental native OpenMC hexagonal face-current Python helper,
  provenance-checked HDF5 export, and a runnable verification example. It scores
  incoming/outgoing/net currents on explicit per-position/per-layer transmission
  faces and retains paired-batch covariance for net-current uncertainty. This is
  an optional diagnostic, not a new MGXS, SPH, or DONJON acceptance result.
- Extended `fill-zero-flux` with opt-in total-XS uncertainty and P0
  scatter-row overshoot criteria for solver-destabilizing micro-flux bins;
  repeated fills now retain the union of substituted-group provenance, and
  the standalone summary schema is now v2.
- Added the IRENA-30 full-core Stage 3 fine CE versus assembly-homogenized
  MG+DONJON full-core workflow, resumable SPH iterations, exact final-sidecar MG evaluation,
  120-degree local-factor regularization, and normalized power-shape closure.
  The empirical eigenvalue-fitted global scalar has been withdrawn from the
  frontend, web execution API, and runnable example. Stage 3 now permits only
  the rate-preserving CE/MG SPH fixed-point iteration with uncertainty and
  convergence gates. The existing physical full-core result is documented as
  **not accepted** pending convergence and independent closure.

## [0.1.4] - 2026-07-10

SPH physics release: update-direction fix (the pre-existing loop was
structurally divergent and even one-shot corrections carried the wrong
sign), rate-preserving `--sph-target rate` (k-pinning, DRAGON
semantics), fast-spectrum zero-flux/floor/freeze policies and the
`fill-zero-flux` command, IRENA Stage 1/2 colorset examples with the
PNL prescription validated to core level, regenerated minicase
evidence, verified DONJON NSPH consumption routes, and web-frontend
wiring for the new options.

See [`RELEASE_NOTES.md`](RELEASE_NOTES.md) for the full narrative.

## [0.1.3] - 2026-07-07

Accepted hex benchmark release. Adds the IRENA-30 ZREFL 91-hex benchmark
(the first accepted hex validation line, SN8 k-eff within Monte Carlo
statistics and fission-source shape within 1.3 % worst / 0.5 % RMS),
the SPH sidecar/augmentation and external face-flux workflows, and the
C5G7 accepted-artifact parity refresh for the uncertainty-preserving
exporter.

See [`RELEASE_NOTES.md`](RELEASE_NOTES.md) for the full narrative.

## [0.1.2] - 2026-05-22

OpenMC end-to-end workflow release. Adds the recipe-based statepoint
exporter, the one-step `openmc2donjon-from-openmc` entry point, ADF and
SPH equivalence carry-through, and
the accepted C5G7 production handoff with flux-ratio ADF and SPH factors.

See [`RELEASE_NOTES.md`](RELEASE_NOTES.md) for the full narrative.

## [0.1.1] - C5G7 handoff snapshot

Tagged as `v0.1.1-c5g7-handoff`. Mid-stream snapshot used to publish the
accepted C5G7 ADF + SPH handoff artifact before the OpenMC workflow
release.

## [0.1.0] - C5G7 accepted baseline

Tagged as `v0.1.0-c5g7-accepted`. First accepted C5G7 assembly-wise
DONJON/DRAGON validation: diffusion and SPN3 k-effective against an
OpenMC reference, with documented bias margins.

[Unreleased]: https://github.com/arbreCui/openmc2donjon/compare/v0.1.4-sph-rate...HEAD
[0.1.4]: https://github.com/arbreCui/openmc2donjon/releases/tag/v0.1.4-sph-rate
[0.1.3]: https://github.com/arbreCui/openmc2donjon/releases/tag/v0.1.3-hex-accepted
[0.1.2]: https://github.com/arbreCui/openmc2donjon/releases/tag/v0.1.2-openmc-workflow
[0.1.1]: https://github.com/arbreCui/openmc2donjon/releases/tag/v0.1.1-c5g7-handoff
[0.1.0]: https://github.com/arbreCui/openmc2donjon/releases/tag/v0.1.0-c5g7-accepted
