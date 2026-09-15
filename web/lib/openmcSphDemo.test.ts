import { describe, expect, it } from "vitest";
import {
  LIVE_OPENMC_SPH_DEMO,
  MOCK_OPENMC_SPH_DEMO,
  openmcSphEvidenceHref,
  openmcSphPlannerPrefill,
} from "./openmcSphDemo";

describe("openmcSphDemo", () => {
  it("prefills the OpenMC planner for the SPH route", () => {
    const prefill = openmcSphPlannerPrefill(MOCK_OPENMC_SPH_DEMO);

    expect(prefill.workflow).toBe("two-step");
    expect(prefill.equivalence).toBe("sph");
    expect(prefill.format).toBe("macrolib");
    expect(prefill.production).toBe(true);
    expect(prefill.keepHdf5Path).toBe(MOCK_OPENMC_SPH_DEMO.exportH5);
    expect(prefill.outputPath).toBe(MOCK_OPENMC_SPH_DEMO.ascii);
    expect(prefill.outputPath).toContain(".macrolib.txt");
    expect(prefill.sphSource).toBe(MOCK_OPENMC_SPH_DEMO.sphSidecar);
  });

  it("keeps the intermediate HDF5 raw so diagnostic augmentation happens once", () => {
    for (const preset of [MOCK_OPENMC_SPH_DEMO, LIVE_OPENMC_SPH_DEMO]) {
      const prefill = openmcSphPlannerPrefill(preset);

      // The planned export must not overwrite the factor-bearing diagnostic
      // artifact or attach the same sidecar twice.
      expect(prefill.keepHdf5Path).not.toBe(preset.augmentedH5);
      // The backend derives the augmented handoff by appending `_sph` to
      // the intermediate stem; the result must be the diagnostic artifact
      // the mock tree and fixture name.
      expect(prefill.keepHdf5Path.replace(/\.h5$/, "_sph.h5")).toBe(
        preset.augmentedH5,
      );
    }
  });

  it("fills the required recipe and a consistent statepoint", () => {
    const prefill = openmcSphPlannerPrefill(MOCK_OPENMC_SPH_DEMO);

    expect(prefill.recipePath).toBe(
      "/mock/home/openmc-runs/openmc-sph-minicase/export_recipe.py",
    );
    // A filled statepoint field must not be silently ignored.
    expect(prefill.loadStatepoint).toBe(prefill.statepointPath.trim().length > 0);
    expect(prefill.statepointPath).toBe(MOCK_OPENMC_SPH_DEMO.ceStatepoint);
  });

  it("links demo evidence without bypassing apply-sph into Converter", () => {
    const evidence = new URL(
      openmcSphEvidenceHref(MOCK_OPENMC_SPH_DEMO),
      "http://localhost:3000",
    );
    expect(evidence.pathname).toBe("/openmc");
    expect(evidence.hash).toBe("#openmc-sph-summary");
    expect(evidence.searchParams.get("equivalence")).toBe("sph");
    expect(evidence.searchParams.get("summary")).toBe(
      MOCK_OPENMC_SPH_DEMO.physicsSummary,
    );
    expect(MOCK_OPENMC_SPH_DEMO.description).toContain(
      "before entering Converter",
    );
  });

  it("calls the augmented artifacts SPH-augmented, not corrected", () => {
    for (const preset of [MOCK_OPENMC_SPH_DEMO, LIVE_OPENMC_SPH_DEMO]) {
      expect(preset.description).not.toMatch(/corrected/i);
    }
  });

  it("points the live minicase to the repository smoke output directory", () => {
    expect(LIVE_OPENMC_SPH_DEMO.runRoot).toBe(
      "/private/tmp/openmc2donjon_two_region_production_20260709",
    );
    expect(LIVE_OPENMC_SPH_DEMO.ceStatepoint).toContain("ce_case/statepoint.80.h5");
    expect(LIVE_OPENMC_SPH_DEMO.physicsSummary).toContain("physics_summary.json");
  });
});
