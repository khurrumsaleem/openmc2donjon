import type { OpenmcEquivalenceMode, OpenmcWorkflowKind } from "./api";

export type OpenmcEntryPointId = "direct-mgxs" | "openmc-sph";

export interface OpenmcEntryPoint {
  id: OpenmcEntryPointId;
  eyebrow: string;
  title: string;
  body: string;
  primaryLabel: string;
  secondaryHref: string;
  secondaryLabel: string;
  workflow: OpenmcWorkflowKind;
  equivalence: OpenmcEquivalenceMode;
  production: boolean;
  check: boolean;
}

export const OPENMC_ENTRY_POINTS: readonly OpenmcEntryPoint[] = [
  {
    id: "direct-mgxs",
    eyebrow: "Need HDF5 first",
    title: "Prepare OpenMC MGXS HDF5",
    body:
      "Start here when your input is an OpenMC recipe/statepoint and you still need the MGXS HDF5. Skip to Converter only when the existing HDF5 is already ready, including any required upstream correction.",
    primaryLabel: "Plan HDF5 export",
    secondaryHref: "/convert?intent=direct-convert&format=multicompo&check=1&production=1",
    secondaryLabel: "Already have HDF5? Open Converter",
    workflow: "two-step",
    equivalence: "direct",
    production: true,
    check: true,
  },
  {
    id: "openmc-sph",
    eyebrow: "Recommended SPH equivalence",
    title: "Preserve rates from CE fine to MG coarse",
    body:
      "Compare a heterogeneous fine-geometry CE reference with a homogenized coarse-geometry MG model. Match the CE tally bins to the MG transport group boundaries, then align physical state, boundary conditions, and the declared fine-to-coarse domain mapping; iterate rate-preserving NSPH, apply the converged factors, then send the corrected HDF5 to Converter.",
    primaryLabel: "Plan CE/MG SPH route",
    secondaryHref:
      "/openmc?workflow=two-step&equivalence=sph&format=multicompo&production=1#openmc-sph-summary",
    secondaryLabel: "Open SPH summary",
    workflow: "two-step",
    equivalence: "sph",
    production: true,
    check: true,
  },
] as const;

export function openmcEntryPoint(id: OpenmcEntryPointId): OpenmcEntryPoint {
  return OPENMC_ENTRY_POINTS.find((item) => item.id === id) ?? OPENMC_ENTRY_POINTS[0];
}

export function activeOpenmcEntryPoint(
  workflow: OpenmcWorkflowKind,
  equivalence: OpenmcEquivalenceMode,
): OpenmcEntryPointId {
  if (workflow === "two-step" && equivalence === "sph") return "openmc-sph";
  return "direct-mgxs";
}
