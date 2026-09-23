import { createElement } from "react";
import { renderToStaticMarkup } from "react-dom/server";
import { describe, expect, it } from "vitest";
import type { SphExecutionResponse } from "../lib/api";
import { appliedSphStatus } from "../lib/appliedSphResult";
import AppliedSphResult from "./AppliedSphResult";

function response(overrides: Partial<SphExecutionResponse> = {}): SphExecutionResponse {
  return {
    schema: "openmc2donjon.web-sph-execution.v1",
    ok: true,
    operation: "apply-sph",
    output_path: "/runs/model/corrected.h5",
    summary_path: "/runs/model/apply.json",
    mixtures: 2,
    energy_groups: 4,
    sph_min: 0.99,
    sph_max: 1.01,
    input_format: "converter",
    binding_mode: "converter-final-exact-input",
    sidecar_input_hash_verified: true,
    mock_mode: false,
    ...overrides,
  };
}

function render(result: SphExecutionResponse) {
  return renderToStaticMarkup(createElement(AppliedSphResult, {
    result,
    projectRoot: "/runs/model with spaces",
    componentId: "assembly-a",
  }));
}

describe("SPH application result", () => {
  it("offers physical checks, not acceptance, only for a verified Converter binding", () => {
    const result = response();
    expect(appliedSphStatus(result).canCheckConverter).toBe(true);
    const html = render(result);
    expect(html).toContain('data-sph-status="converter-bound"');
    expect(html).toContain("physical checks still required");
    expect(html).toContain("not physical acceptance");
    expect(html).toContain("Check physical SPH in Converter");
    expect(html).toContain('href="/convert?contract=physical-sph&amp;check=1&amp;production=1');
    expect(html).toContain("input=%2Fruns%2Fmodel%2Fcorrected.h5");
    expect(html).toContain("project=%2Fruns%2Fmodel+with+spaces&amp;component=assembly-a");
    expect(html).toContain("/runs/model/apply.json");
  });

  it("warns on unbound numerical output and keeps it out of the final handoff route", () => {
    const result = response({ binding_mode: "converter-unbound", sidecar_input_hash_verified: false });
    const html = render(result);
    expect(html).toContain('data-sph-status="converter-unbound"');
    expect(html).toContain("input binding not verified");
    expect(html).toContain("exact uncorrected Converter-layout input");
    expect(html).toContain("Do not use this corrected output as the new reference input");
    expect(html).toContain("border-amber");
    expect(html).not.toContain('href="/convert');
  });

  it("directs native MG output to the next MG run, not to Converter", () => {
    const result = response({
      input_format: "openmc-mgxs",
      binding_mode: "openmc-mgxs-intermediate-unbound",
      sidecar_input_hash_verified: false,
    });
    const html = render(result);
    expect(html).toContain('data-sph-status="mg-iteration"');
    expect(html).toContain("not a Converter handoff");
    expect(html).toContain("Rerun that MG model");
    expect(html).toContain("export the new MG flux");
    expect(html).not.toContain('href="/convert');
  });

  it.each([
    { input_format: undefined },
    { binding_mode: undefined },
    { sidecar_input_hash_verified: undefined },
    { input_format: "openmc-mgxs" as const },
    { binding_mode: "converter-unbound" },
    { binding_mode: "future-unknown-mode" },
    { sidecar_input_hash_verified: false },
    { sidecar_input_hash_verified: "true" as unknown as boolean },
    { operation: "sph-sidecar" as const },
    { output_path: "  " },
    { ok: false },
  ])("fails closed on incomplete or contradictory evidence: %j", (overrides) => {
    const result = response(overrides);
    expect(appliedSphStatus(result)).toMatchObject({ kind: "unverified", canCheckConverter: false });
    const html = render(result);
    expect(html).toContain("Handoff status unverified");
    expect(html).not.toContain('href="/convert');
    expect(html).toContain("border-amber");
  });

  it("does not infer evidence from a legacy success response", () => {
    const result = response({
      input_format: undefined,
      binding_mode: undefined,
      sidecar_input_hash_verified: undefined,
      mock_mode: undefined,
    });
    expect(appliedSphStatus(result).canCheckConverter).toBe(false);
    expect(render(result)).not.toContain('href="/convert');
  });

  it("never treats mock success as real evidence, even with contradictory verified fields", () => {
    const result = response({ mock_mode: true });
    expect(appliedSphStatus(result)).toMatchObject({ kind: "mock", canCheckConverter: false });
    const html = render(result);
    expect(html).toContain("Demo result");
    expect(html).toContain("Demo output");
    expect(html).not.toContain('href="/convert');
  });
});
