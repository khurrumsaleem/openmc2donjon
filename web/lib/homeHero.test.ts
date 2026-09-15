import { describe, expect, it } from "vitest";
import {
  HOME_EQUIVALENCE_FLOW,
  HOME_HANDOFF_FLOW,
  HOME_HERO,
} from "./homeHero";

describe("home hero", () => {
  it("leads with the concrete Converter product", () => {
    expect(HOME_HERO.kicker).toBe("OpenMC → DRAGON / DONJON handoff");
    expect(HOME_HERO.heading).toBe(
      "Convert OpenMC MGXS into a traceable DRAGON/DONJON object.",
    );
    expect(HOME_HERO.paragraph).toContain("L_MULTICOMPO or L_MACROLIB");
    expect(HOME_HERO.paragraph).toContain("hash-linked receipt");
  });

  it("shows the recommended SPH route before the formal Converter boundary", () => {
    expect(HOME_EQUIVALENCE_FLOW.map((stage) => stage.label)).toEqual([
      "OpenMC CE fine + MG coarse",
      "Rate-preserving SPH",
      "Corrected / Converter-ready HDF5",
    ]);
    expect(HOME_EQUIVALENCE_FLOW[0].qualifier).toContain("different geometries");
    expect(HOME_EQUIVALENCE_FLOW[1].qualifier).toContain("recommended");
    expect(HOME_EQUIVALENCE_FLOW[2].qualifier).toContain("Converter-ready");

    expect(HOME_HANDOFF_FLOW.map((stage) => stage.label)).toEqual([
      "Converter-ready HDF5",
      "Converter",
      "L_MULTICOMPO / L_MACROLIB",
      "DONJON",
    ]);
    expect(HOME_HANDOFF_FLOW[1].qualifier).toContain("required");
    expect(HOME_HANDOFF_FLOW[2].qualifier).toContain("hash-linked receipt");
  });

  it("places Converter at the center without erasing SPH or projects", () => {
    expect(HOME_HERO.paragraph).toContain("formal handoff boundary");
    expect(HOME_HERO.supporting).toContain("heterogeneous OpenMC CE");
    expect(HOME_HERO.supporting).toContain("homogenized OpenMC MG");
    expect(HOME_HERO.supporting).toContain("only then enter Converter");
    expect(HOME_HERO.supporting).toContain("advanced, project-specific");
    expect(HOME_HERO.supporting).toContain("PyGan/LCM is optional");
  });
});
