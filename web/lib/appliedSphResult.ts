import type { SphExecutionResponse } from "./api";

export interface AppliedSphStatus {
  kind: "converter-bound" | "converter-unbound" | "mg-iteration" | "unverified" | "mock";
  title: string;
  detail: string;
  next: string;
  canCheckConverter: boolean;
}

/** Numerical success is not physical acceptance; use the returned artifact's evidence. */
export function appliedSphStatus(result: SphExecutionResponse): AppliedSphStatus {
  if (result.mock_mode === true) {
    return {
      kind: "mock",
      title: "Demo result — no physical evidence",
      detail: "This response is simulated. It does not verify an input binding or a physical SPH result.",
      next: "Run apply-sph with real files in the live backend before checking a handoff.",
      canCheckConverter: false,
    };
  }
  if (result.ok === true && result.operation === "apply-sph" && result.output_path.trim()) {
    if (
      result.input_format === "converter" &&
      result.binding_mode === "converter-final-exact-input" &&
      result.sidecar_input_hash_verified === true
    ) {
      return {
        kind: "converter-bound",
        title: "Input binding verified — physical checks still required",
        detail: "The sidecar matches the exact Converter input. This confirms the binding, not physical acceptance.",
        next: "Run Converter's physical-SPH checks before writing the DRAGON/DONJON object.",
        canCheckConverter: true,
      };
    }
    if (
      result.input_format === "openmc-mgxs" &&
      result.binding_mode === "openmc-mgxs-intermediate-unbound" &&
      result.sidecar_input_hash_verified === false
    ) {
      return {
        kind: "mg-iteration",
        title: "OpenMC MG iteration file — not a Converter handoff",
        detail: "This native setN library is an intermediate input for the homogenized OpenMC MG model.",
        next: "Rerun that MG model with this file, export the new MG flux, then compute the next SPH update. After convergence and independent validation, apply the final factors to the exact bound Converter-layout input.",
        canCheckConverter: false,
      };
    }
    if (
      result.input_format === "converter" &&
      result.binding_mode === "converter-unbound" &&
      result.sidecar_input_hash_verified === false
    ) {
      return {
        kind: "converter-unbound",
        title: "Numerical result only — input binding not verified",
        detail: "Cross sections were scaled, but the sidecar does not carry a verified hash of this input. This is not a verified physical-SPH handoff.",
        next: "Recompute the sidecar against the exact uncorrected Converter-layout input, then apply it again. Do not use this corrected output as the new reference input.",
        canCheckConverter: false,
      };
    }
  }
  return {
    kind: "unverified",
    title: "Handoff status unverified",
    detail: "The response does not contain consistent input-format and binding evidence. An output path alone is not enough.",
    next: "Check that the backend is up to date and rerun apply-sph to obtain an explicit binding report.",
    canCheckConverter: false,
  };
}
