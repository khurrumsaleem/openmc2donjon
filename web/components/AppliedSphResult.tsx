import React from "react";
import Link from "next/link";
import type { SphExecutionResponse } from "../lib/api";
import { appliedSphStatus } from "../lib/appliedSphResult";
import { equivalenceAppliedHandoffHref } from "../lib/equivalenceRoutes";

export default function AppliedSphResult({
  result,
  projectRoot,
  componentId,
}: {
  result: SphExecutionResponse;
  projectRoot?: string;
  componentId?: string | null;
}) {
  const status = appliedSphStatus(result);
  return (
    <div
      role="status"
      aria-live="polite"
      data-sph-status={status.kind}
      className={
        "mt-3 rounded-lg border p-3 text-[12px] leading-5 " +
        (status.canCheckConverter
          ? "border-emerald-300/20 bg-emerald-300/[0.06] text-emerald-100"
          : "border-amber-300/25 bg-amber-300/[0.06] text-amber-100")
      }
    >
      <p className="font-semibold">{status.title}</p>
      <p className="mt-1 break-all font-mono text-[11px]">
        {status.kind === "mock" ? "Demo output: " : "Output: "}{result.output_path}
      </p>
      <p className="mt-1 text-[11px] opacity-80">
        {result.mixtures} mixtures · {result.energy_groups} groups · SPH {result.sph_min.toFixed(6)}…{result.sph_max.toFixed(6)}
      </p>
      <p className="mt-2">{status.detail}</p>
      <p className="mt-2"><span className="font-semibold">Next: </span>{status.next}</p>
      {result.summary_path ? (
        <p className="mt-2 break-all font-mono text-[11px] opacity-80">Summary: {result.summary_path}</p>
      ) : null}
      {status.canCheckConverter ? (
        <Link
          href={equivalenceAppliedHandoffHref({
            inputH5: result.output_path,
            projectRoot,
            componentId,
          })}
          className="mt-3 inline-flex min-h-9 items-center font-semibold underline-offset-4 hover:underline"
        >
          Check physical SPH in Converter →
        </Link>
      ) : null}
    </div>
  );
}
