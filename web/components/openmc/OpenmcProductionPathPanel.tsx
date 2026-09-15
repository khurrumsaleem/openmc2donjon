"use client";

import Link from "next/link";
import { CopyCliButton } from "@/components/commands/CopyCliButton";
import type {
  ConvertFormat,
  OpenmcEquivalenceMode,
  OpenmcExportExecutionResponse,
  OpenmcWorkflowPlan,
} from "@/lib/api";
import type { OpenmcSphDemoPreset } from "@/lib/openmcSphDemo";
import { openmcSphEvidenceHref } from "@/lib/openmcSphDemo";
import {
  OPENMC_SPH_APPLY_FORM_HREF,
  openmcDirectConvertHref,
  OPENMC_SPH_SIDECAR_FORM_HREF,
  openmcWalkthroughStatuses,
  type OpenmcWalkthroughRun,
  type OpenmcWalkthroughStatus,
} from "@/lib/openmcWorkflowWalkthrough";
import type { ProjectComponentRouteContext } from "@/lib/projectWorkspace";
import { openmcExportProvenanceVerified } from "@/lib/openmcExportGate";

export type OpenmcPlannerState =
  | { kind: "idle" }
  | { kind: "loading" }
  | { kind: "ok"; data: OpenmcWorkflowPlan }
  | { kind: "error"; message: string; status?: number };

export interface OpenmcSphDemoActions {
  preset: OpenmcSphDemoPreset;
  mode: "mock" | "live";
  onFill: () => void;
  onReview: () => void;
}

export type OpenmcExportState =
  | { kind: "idle" }
  | { kind: "loading" }
  | { kind: "ok"; data: OpenmcExportExecutionResponse }
  | { kind: "error"; message: string };

export default function OpenmcProductionPathPanel({
  state,
  exportState,
  equivalence,
  format,
  production,
  recipePath,
  statepointPath,
  loadStatepoint,
  runDir,
  projectContext,
  demo = null,
}: {
  state: OpenmcPlannerState;
  exportState: OpenmcExportState;
  equivalence: OpenmcEquivalenceMode;
  format: ConvertFormat;
  production: boolean;
  recipePath: string;
  statepointPath: string;
  loadStatepoint: boolean;
  runDir: string;
  projectContext?: ProjectComponentRouteContext;
  demo?: OpenmcSphDemoActions | null;
}) {
  const planned = state.kind === "ok" && state.data.ok;
  const plan = state.kind === "ok" ? state.data : null;
  const plannedStatuses = openmcWalkthroughStatuses({
    hasRecipe: recipePath.trim().length > 0,
    hasStatepoint: statepointPath.trim().length > 0,
    loadStatepoint,
    hasRunDir: runDir.trim().length > 0,
    run: openmcWalkthroughRunFromState(state),
  });
  const object = format === "macrolib" ? "L_MACROLIB" : "L_MULTICOMPO";
  const objectShort = format === "macrolib" ? "MACROLIB" : "MULTICOMPO";
  const isOpenmcSph = equivalence === "sph";
  const isSphExport = isOpenmcSph && plan?.plan_scope === "export";
  const exported =
    exportState.kind === "ok" &&
    exportState.data.ok &&
    !exportState.data.mock_mode;
  const provenanceVerified =
    exported && openmcExportProvenanceVerified(exportState.data, production);
  const actualHdf5Path = exported ? exportState.data.output_path : null;
  const inspectHref = actualHdf5Path
    ? `/inspect?path=${encodeURIComponent(actualHdf5Path)}`
    : null;
  const plannedAsciiPath =
    plan?.artifacts.find((artifact) => artifact.kind === "ascii")?.path ?? "";
  const convertHref =
    planned && plan && provenanceVerified && !isOpenmcSph && actualHdf5Path
      ? openmcDirectConvertHref(
          actualHdf5Path,
          plannedAsciiPath,
          format,
          production,
          projectContext,
        )
      : null;
  const exportStatus: OpenmcWalkthroughStatus =
    exportState.kind === "loading"
      ? "planning"
      : exportState.kind === "error"
        ? "blocked"
        : exported
          ? "written"
          : exportState.kind === "ok"
            ? "planned"
            : planned
              ? "ready"
              : plannedStatuses.run;
  const verificationStatus: OpenmcWalkthroughStatus = provenanceVerified
    ? "verified"
    : exported || exportState.kind === "error"
      ? "blocked"
      : "needed";
  const statuses = {
    ...plannedStatuses,
    run: exportStatus,
    review: verificationStatus,
    bundle: "needed" as OpenmcWalkthroughStatus,
  };
  const sphDemo = isOpenmcSph && demo ? demo : null;
  const items = isSphExport
    ? [
        {
          id: "export",
          label: "01",
          eyebrow: "OpenMC export",
          title: "Write the MGXS HDF5",
          body:
            "Run the generated export command with the fine-reference CE statepoint and the project recipe that declares the equivalence domains.",
          status: statuses.run,
          href: undefined,
          hrefLabel: undefined,
          onClick: undefined,
        },
        {
          id: "inspect",
          label: "02",
          eyebrow: "Confirm domains",
          title: "Inspect the exported HDF5",
          body:
            "Confirm the declared domains and ordering, energy groups, scattering order, uncertainty data, and expected mesh before solving SPH.",
          status: statuses.review,
          href: inspectHref ?? undefined,
          hrefLabel: inspectHref ? "Inspect exported HDF5" : undefined,
          onClick: undefined,
        },
        {
          id: "sph",
          label: "03",
          eyebrow: "Next page",
          title: "Build and apply SPH",
          body:
            "Use a heterogeneous fine CE reference and a homogenized coarse MG model. Score CE tallies on the MG group boundaries and align state, BC, and domain mapping. Converter comes after apply-sph, not directly after export.",
          status: statuses.review,
          href: OPENMC_SPH_SIDECAR_FORM_HREF,
          hrefLabel: "Continue to SPH",
          onClick: undefined,
        },
      ]
    : isOpenmcSph
    ? [
        {
          id: "run-openmc",
          label: "01",
          eyebrow: "Run physics",
          title: "Run OpenMC CE/MG SPH",
          body:
            "Run the heterogeneous fine-reference and homogenized coarse-MG OpenMC models on intentionally different geometries. Score CE tallies on the MG group boundaries and align state, boundary conditions, and fine-to-coarse domain mapping; converge SPH(domain, group), then apply NSPH to write the corrected HDF5.",
          status: statuses.run,
          href: undefined,
          hrefLabel: undefined,
          onClick: undefined,
        },
        {
          id: "summary",
          label: "02",
          eyebrow: "Review evidence",
          title: "Review SPH evidence",
          body:
            "Load physics_summary.json and check CE/MG flux uncertainty, SPH factor range, reaction-rate preservation, and NSPH handoff status.",
          status: statuses.run,
          href: sphDemo ? openmcSphEvidenceHref(sphDemo.preset) : undefined,
          hrefLabel: sphDemo ? "Load demo evidence" : undefined,
          onClick: sphDemo?.onReview,
        },
        {
          id: "apply",
          label: "03",
          eyebrow: "Corrected handoff",
          title: "Apply SPH before Converter",
          body:
            `Run apply-sph to write the corrected HDF5. Only then enter Converter to write ${object}; Converter is the formal handoff boundary, not part of the SPH iteration.`,
          status: statuses.review,
          href: OPENMC_SPH_APPLY_FORM_HREF,
          hrefLabel: "Open apply-sph",
          onClick: undefined,
        },
      ]
    : [
        {
          id: "source",
          label: "01",
          eyebrow: "Export source",
          title: loadStatepoint ? "Export OpenMC MGXS" : "Prepare export command",
          body: loadStatepoint
            ? "Select the export recipe and statepoint that define the OpenMC MGXS handoff."
            : "Use the recipe without loading a statepoint when you only need the generated command scaffold.",
          status: statuses.run,
          href: undefined,
          hrefLabel: undefined,
          onClick: undefined,
        },
        {
          id: "convert",
          label: "02",
          eyebrow: "Converter",
          title: `Convert HDF5 to ${object}`,
          body: `Run the OpenMC export first, inspect or augment the HDF5 when needed, then use Converter to write ${object}.`,
          status: statuses.review,
          href: convertHref ?? inspectHref ?? undefined,
          hrefLabel: convertHref ? "Open Converter" : "Inspect HDF5",
          onClick: undefined,
        },
        {
          id: "bundle",
          label: "03",
          eyebrow: "Bundle",
          title: "Package the DONJON bundle",
          body:
            "Use a managed run directory to keep the MGXS input, ASCII output, summaries, and manifest together.",
          status: statuses.bundle,
          href: undefined,
          hrefLabel: "Open bundle builder",
          onClick: undefined,
        },
      ];

  return (
    <section className="mb-5 rounded-xl border border-cyan-300/20 bg-cyan-300/[0.035] p-4">
      <div className="flex flex-wrap items-start justify-between gap-4">
        <div>
          <div className="text-[10px] uppercase tracking-[0.14em] text-cyan-200/80">
            Main workflow
          </div>
          <h2 className="mt-1 text-base font-semibold tracking-tight">
            {isSphExport
              ? "Export the MGXS source before SPH"
              : isOpenmcSph
              ? `Prepare the input before Converter writes ${objectShort}`
              : "Prepare the input before Converter writes DONJON ASCII"}
          </h2>
          <p className="mt-1 max-w-3xl text-sm leading-relaxed text-[var(--fg-2)]">
            {isSphExport
              ? "This step ends at the MGXS HDF5. Inspect it, then build and apply SPH on the next page; Converter remains the single gateway that validates and writes DONJON ASCII."
              : isOpenmcSph
              ? "This route prepares one corrected Converter input from a heterogeneous fine reference and homogenized coarse model. Converter begins only after apply-sph; the project manifest decides how many other objects and which consumer are required."
              : "This route prepares the MGXS HDF5 from OpenMC inputs. Converter is enabled only after a real write and a fail-closed provenance check; a successful plan alone is not an artifact."}
          </p>
          {demo ? (
            <p className="mt-1 max-w-3xl text-[12px] leading-5 text-emerald-200/90">
              Demo: the {demo.mode === "mock" ? "bundled mock" : "live production"}{" "}
              two-region minicase can prefill this workflow; its per-step links
              load the evidence, converter, and bundle.
            </p>
          ) : null}
        </div>
        <div className="flex shrink-0 flex-wrap items-center gap-3">
          <span className="rounded border border-[var(--edge)] px-2 py-1 font-mono text-[11px] uppercase tracking-wider text-[var(--fg-2)]">
            {isOpenmcSph || isSphExport
              ? "CE fine + MG coarse → SPH → corrected HDF5 → Converter"
              : `OpenMC HDF5 → Converter · ${equivalenceLabel(equivalence)}`}
          </span>
          {demo ? (
            <>
              <span className="rounded border border-emerald-300/30 px-2 py-1 text-[10px] uppercase tracking-[0.14em] text-emerald-300">
                {demo.mode} demo
              </span>
              <button
                type="button"
                onClick={demo.onFill}
                className="btn btn-primary"
              >
                Fill SPH planner
              </button>
            </>
          ) : null}
        </div>
      </div>

      <div className="mt-4 grid gap-2 lg:grid-cols-3">
        {items.map((item) => (
          <article
            key={item.id}
            className={
              "rounded-lg border px-3 py-2 " +
              openmcWalkthroughStatusClass(item.status)
            }
          >
            <div className="flex items-center justify-between gap-2">
              <span className="rounded border border-current/25 px-1.5 py-0.5 font-mono text-[10px]">
                {item.label}
              </span>
              <span className="text-[10px] uppercase tracking-[0.14em] opacity-80">
                {item.status}
              </span>
            </div>
            <div className="mt-2 text-[10px] uppercase tracking-[0.14em] opacity-70">
              {item.eyebrow}
            </div>
            <h3 className="mt-2 text-sm font-semibold tracking-tight">
              {item.title}
            </h3>
            <p className="mt-1 text-[12px] leading-5 text-[var(--fg-2)]">
              {item.body}
            </p>
            {item.href && item.hrefLabel ? (
              <Link
                href={item.href}
                onClick={item.onClick}
                className="mt-3 inline-flex text-[12px] font-medium text-[var(--accent-2)] hover:underline"
              >
                {item.hrefLabel}
              </Link>
            ) : null}
          </article>
        ))}
      </div>

      {exported && !provenanceVerified ? (
        <div className="mt-4 rounded-lg border border-rose-300/25 bg-rose-300/[0.055] px-3 py-2 text-[12px] leading-5 text-rose-100">
          The HDF5 was written, but its embedded OpenMC provenance is not
          sufficient for this {production ? "production" : "engineering"}
          handoff. Inspect is available for diagnosis; Converter remains blocked.
        </div>
      ) : null}

      {planned && plan ? (
        <div className="mt-4 flex flex-wrap items-center gap-3">
          {planned && plan ? (
            <CopyCliButton
              value={plan.primary_command_text}
              label="Copy primary command"
              copiedLabel="Copied"
            />
          ) : null}
          {convertHref ? (
            <Link href={convertHref} className="btn btn-secondary">
              Open Converter
            </Link>
          ) : null}
          {inspectHref ? (
            <Link href={inspectHref} className="btn btn-secondary">
              Inspect HDF5
            </Link>
          ) : null}
        </div>
      ) : null}
    </section>
  );
}

function openmcWalkthroughRunFromState(
  state: OpenmcPlannerState,
): OpenmcWalkthroughRun {
  if (state.kind === "loading") return { kind: "loading" };
  if (state.kind === "ok") return { kind: "ok", ok: state.data.ok };
  return { kind: state.kind };
}

function openmcWalkthroughStatusClass(status: OpenmcWalkthroughStatus): string {
  if (
    status === "ready" ||
    status === "planning" ||
    status === "written" ||
    status === "verified"
  ) {
    return "border-cyan-300/25 bg-cyan-300/[0.06] text-cyan-100";
  }
  if (status === "planned") {
    return "border-amber-300/20 bg-amber-300/[0.045] text-amber-100";
  }
  if (status === "blocked") {
    return "border-rose-400/25 bg-rose-400/[0.06] text-rose-100";
  }
  if (status === "optional") {
    return "border-[var(--edge)] bg-white/[0.02] text-[var(--fg-3)]";
  }
  return "border-[var(--edge)] bg-white/[0.02] text-[var(--fg-2)]";
}

function equivalenceLabel(value: OpenmcEquivalenceMode): string {
  if (value === "adf") return "ADF/DF";
  if (value === "sph") return "SPH";
  if (value === "flux-ratio-adf") return "flux-ratio ADF";
  return "direct";
}
