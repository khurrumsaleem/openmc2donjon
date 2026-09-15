import { describe, expect, it } from "vitest";
import {
  OPENMC_SPH_FIXED_POLICY,
  OPENMC_SPH_UPDATE_GATE,
  openmcSphUncertaintyGate,
  parseOpenmcSphDamping,
} from "./openmcSphExecution";

describe("OpenMC-side SPH execution policy", () => {
  it("exposes the fixed 2% update gate without calling it final acceptance", () => {
    expect(OPENMC_SPH_UPDATE_GATE).toBe(0.02);
    expect(OPENMC_SPH_FIXED_POLICY).toContainEqual([
      "Numerical exemptions",
      "None: no clipping, floors, or frozen groups",
    ]);
    expect(OPENMC_SPH_FIXED_POLICY).toContainEqual([
      "Uncertainty limits",
      "Explicit CE and MG project thresholds; no default",
    ]);
  });

  it("fails closed on missing, non-finite, or out-of-range damping", () => {
    for (const value of ["", "abc", "NaN", "Infinity", "-0.1", "1.01"]) {
      expect(parseOpenmcSphDamping(value).ok, value).toBe(false);
    }
    expect(parseOpenmcSphDamping("0")).toEqual({ ok: true, value: 0 });
    expect(parseOpenmcSphDamping("0.6")).toEqual({ ok: true, value: 0.6 });
    expect(parseOpenmcSphDamping("1")).toEqual({ ok: true, value: 1 });
  });

  it("holds production SPH until both explicit project thresholds are valid", () => {
    const missing = openmcSphUncertaintyGate("rate", "", "");
    expect(missing.ok).toBe(false);
    expect(missing.issues).toEqual([
      "CE max relative std dev is required for rate-preserving physical SPH.",
      "MG max relative std dev is required for rate-preserving physical SPH.",
    ]);

    const invalid = openmcSphUncertaintyGate("rate", "-0.01", "not-a-number");
    expect(invalid.ok).toBe(false);
    expect(invalid.issues).toEqual([
      "CE max relative std dev must be a finite non-negative number.",
      "MG max relative std dev must be a finite non-negative number.",
    ]);

    expect(openmcSphUncertaintyGate("rate", "0.025", "0.04")).toEqual({
      required: true,
      ok: true,
      ceMaxRelativeStdDev: 0.025,
      mgMaxRelativeStdDev: 0.04,
      issues: [],
    });
  });

  it("allows diagnostic flux targeting to omit uncertainty thresholds", () => {
    expect(openmcSphUncertaintyGate("flux", "", "")).toEqual({
      required: false,
      ok: true,
      ceMaxRelativeStdDev: null,
      mgMaxRelativeStdDev: null,
      issues: [],
    });
    expect(openmcSphUncertaintyGate("flux", "bad", "").ok).toBe(false);
  });
});
