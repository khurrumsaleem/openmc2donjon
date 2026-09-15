export const OPENMC_SPH_UPDATE_GATE = 0.02;

export const OPENMC_SPH_FIXED_POLICY = [
  ["Target", "Reaction-rate preserving"],
  ["Normalization", "H-factor / kappa-fission power (auto)"],
  ["Zero-flux bins", "Reject"],
  ["Uncertainty limits", "Explicit CE and MG project thresholds; no default"],
  ["Numerical exemptions", "None: no clipping, floors, or frozen groups"],
] as const;

export type DampingParseResult =
  | { ok: true; value: number }
  | { ok: false; message: string };

export function parseOpenmcSphDamping(value: string): DampingParseResult {
  const trimmed = value.trim();
  if (!trimmed) {
    return { ok: false, message: "Damping is required." };
  }
  const parsed = Number(trimmed);
  if (!Number.isFinite(parsed) || parsed < 0 || parsed > 1) {
    return {
      ok: false,
      message: "Damping must be a finite number within 0..1.",
    };
  }
  return { ok: true, value: parsed };
}

export interface OpenmcSphUncertaintyGate {
  required: boolean;
  ok: boolean;
  ceMaxRelativeStdDev: number | null;
  mgMaxRelativeStdDev: number | null;
  issues: string[];
}

/**
 * Production rate-preserving SPH fails closed until the project declares
 * independent CE and MG tally-uncertainty limits. Diagnostic flux targeting
 * may omit them, but any supplied limit must still be finite and non-negative.
 */
export function openmcSphUncertaintyGate(
  sphTarget: "flux" | "rate",
  ceValue: string,
  mgValue: string,
): OpenmcSphUncertaintyGate {
  const required = sphTarget === "rate";
  const ce = parseUncertaintyThreshold(
    ceValue,
    "CE max relative std dev",
    required,
  );
  const mg = parseUncertaintyThreshold(
    mgValue,
    "MG max relative std dev",
    required,
  );
  const issues = [...ce.issues, ...mg.issues];
  return {
    required,
    ok: issues.length === 0,
    ceMaxRelativeStdDev: ce.value,
    mgMaxRelativeStdDev: mg.value,
    issues,
  };
}

function parseUncertaintyThreshold(
  raw: string,
  label: string,
  required: boolean,
): { value: number | null; issues: string[] } {
  const value = raw.trim();
  if (value === "") {
    return {
      value: null,
      issues: required
        ? [`${label} is required for rate-preserving physical SPH.`]
        : [],
    };
  }
  const parsed = Number(value);
  if (!Number.isFinite(parsed) || parsed < 0) {
    return {
      value: null,
      issues: [`${label} must be a finite non-negative number.`],
    };
  }
  return { value: parsed, issues: [] };
}
