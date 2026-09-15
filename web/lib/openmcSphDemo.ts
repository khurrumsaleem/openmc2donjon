import type { ConvertFormat, OpenmcEquivalenceMode, OpenmcWorkflowKind } from "./api";

export interface OpenmcSphDemoPreset {
  id: "mock-openmc-sph" | "live-openmc-sph";
  label: string;
  description: string;
  runRoot: string;
  recipe: string;
  ceStatepoint: string;
  sphSidecar: string;
  // Raw MGXS export target (pre-SPH). The workflow planner derives the
  // SPH-augmented handoff name by appending `_sph` to this stem, so this
  // is named such that the augmentation yields exactly `augmentedH5`.
  exportH5: string;
  augmentedH5: string;
  ascii: string;
  physicsSummary: string;
}

export interface OpenmcSphPlannerPrefill {
  workflow: OpenmcWorkflowKind;
  equivalence: OpenmcEquivalenceMode;
  format: ConvertFormat;
  production: boolean;
  check: boolean;
  runDir: string;
  keepHdf5Path: string;
  outputPath: string;
  sphSource: string;
  recipePath: string;
  statepointPath: string;
  loadStatepoint: boolean;
}

export const MOCK_OPENMC_SPH_DEMO: OpenmcSphDemoPreset = {
  id: "mock-openmc-sph",
  label: "Mock CE/MG SPH diagnostic",
  description:
    "Prefill the CE/MG SPH inputs with bundled mock paths. Apply the converged factors before entering Converter.",
  runRoot: "/mock/home/openmc-runs/openmc-sph-minicase",
  recipe: "/mock/home/openmc-runs/openmc-sph-minicase/export_recipe.py",
  ceStatepoint: "/mock/home/openmc-runs/openmc-sph-minicase/ce_statepoint.h5",
  sphSidecar: "/mock/home/openmc-runs/openmc-sph-minicase/openmc_sph_sidecar.h5",
  exportH5: "/mock/home/openmc-runs/openmc-sph-minicase/mgxs_with_openmc.h5",
  augmentedH5: "/mock/home/openmc-runs/openmc-sph-minicase/mgxs_with_openmc_sph.h5",
  ascii: "/mock/home/openmc-runs/openmc-sph-minicase/out.macrolib.txt",
  physicsSummary: "/mock/home/openmc-runs/openmc-sph-minicase/physics_summary.json",
};

export const LIVE_OPENMC_SPH_DEMO: OpenmcSphDemoPreset = {
  id: "live-openmc-sph",
  label: "Two-region OpenMC-side SPH diagnostic",
  description:
    "Review the high-statistics same-partition flux-target example: two output regions produce two SPH factors per group and exercise HDF5-to-DONJON data carriage. This is not physical-SPH acceptance.",
  runRoot: "/private/tmp/openmc2donjon_two_region_production_20260709",
  recipe: "examples/openmc_ce_mg_33g_sph_minicase/export_recipe.py",
  ceStatepoint:
    "/private/tmp/openmc2donjon_two_region_production_20260709/ce_case/statepoint.80.h5",
  sphSidecar:
    "/private/tmp/openmc2donjon_two_region_production_20260709/handoff/openmc_sph_sidecar.h5",
  exportH5:
    "/private/tmp/openmc2donjon_two_region_production_20260709/handoff/mgxs_with_openmc.h5",
  augmentedH5:
    "/private/tmp/openmc2donjon_two_region_production_20260709/handoff/mgxs_with_openmc_sph.h5",
  ascii:
    "/private/tmp/openmc2donjon_two_region_production_20260709/handoff/out_with_openmc_sph.macrolib.txt",
  physicsSummary:
    "/private/tmp/openmc2donjon_two_region_production_20260709/handoff/physics_summary.json",
};

export function openmcSphPlannerPrefill(
  preset: OpenmcSphDemoPreset,
): OpenmcSphPlannerPrefill {
  return {
    workflow: "two-step",
    equivalence: "sph",
    format: "macrolib",
    production: true,
    check: true,
    runDir: preset.runRoot,
    keepHdf5Path: preset.exportH5,
    outputPath: preset.ascii,
    sphSource: preset.sphSidecar,
    recipePath: preset.recipe,
    statepointPath: preset.ceStatepoint,
    loadStatepoint: true,
  };
}

export function openmcSphEvidenceHref(preset: OpenmcSphDemoPreset): string {
  const params = new URLSearchParams({
    workflow: "two-step",
    equivalence: "sph",
    format: "macrolib",
    production: "1",
    summary: preset.physicsSummary,
  });
  return `/openmc?${params.toString()}#openmc-sph-summary`;
}
