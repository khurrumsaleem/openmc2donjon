import { describe, expect, it } from "vitest";
import {
  OPENMC_ENTRY_POINTS,
  activeOpenmcEntryPoint,
  openmcEntryPoint,
} from "./openmcEntryPoints";

describe("openmcEntryPoints", () => {
  it("offers direct MGXS and OpenMC-side SPH as the two primary entries", () => {
    expect(OPENMC_ENTRY_POINTS.map((entry) => entry.id)).toEqual([
      "direct-mgxs",
      "openmc-sph",
    ]);
  });

  it("maps the SPH entry to two-step OpenMC-side SPH planning", () => {
    const sph = openmcEntryPoint("openmc-sph");

    expect(sph.workflow).toBe("two-step");
    expect(sph.equivalence).toBe("sph");
    expect(sph.production).toBe(true);
    expect(sph.secondaryHref).toContain("equivalence=sph");
    expect(sph.secondaryHref).toContain("format=multicompo");
    expect(sph.secondaryHref).not.toContain("colorset=");
    // The label promises the summary card, so the href must land on it.
    expect(sph.secondaryHref).toContain("#openmc-sph-summary");
  });

  it("uses production checks for the direct Converter shortcut", () => {
    const direct = openmcEntryPoint("direct-mgxs");
    expect(direct.production).toBe(true);
    expect(direct.secondaryHref).toContain("production=1");
  });

  it("describes record-attachment as attach, never inject", () => {
    for (const entry of OPENMC_ENTRY_POINTS) {
      expect(entry.body).not.toMatch(/inject/i);
    }
  });

  it("keeps the SPH entry generic about geometry and component count", () => {
    const sph = openmcEntryPoint("openmc-sph");
    expect(sph.body).toContain("fine-geometry");
    expect(sph.body).toContain("coarse-geometry");
    expect(sph.body).toContain("CE tally bins");
    expect(sph.body).toContain("MG transport group boundaries");
    expect(sph.body).toContain("physical state");
    expect(sph.body).toContain("boundary conditions");
    expect(sph.body).toContain("fine-to-coarse domain mapping");
    expect(sph.body).toContain("corrected HDF5 to Converter");
    expect(sph.body).not.toMatch(/seven|colorset|five|91-position/i);
  });

  it("identifies the active entry from the planner state", () => {
    expect(activeOpenmcEntryPoint("two-step", "sph")).toBe("openmc-sph");
    expect(activeOpenmcEntryPoint("one-step", "sph")).toBe("direct-mgxs");
    expect(activeOpenmcEntryPoint("two-step", "direct")).toBe("direct-mgxs");
  });
});
