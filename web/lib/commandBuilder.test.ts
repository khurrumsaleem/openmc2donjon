import { describe, expect, it } from "vitest";
import {
  builderCliIssues,
  builderValuesFromQuery,
  buildCommandCli,
  commandBuilderStage,
  commandBuilderSpec,
  defaultBuilderValues,
} from "./commandBuilder";

describe("commandBuilder", () => {
  it("builds repeated flags and numeric options for diff", () => {
    const spec = commandBuilderSpec("diff");
    expect(spec).not.toBeNull();
    const values = defaultBuilderValues(spec!);
    values.reference_h5 = "/runs/ref.h5";
    values.candidate_h5 = "/runs/candidate.h5";
    values.rtol = "1e-5";
    values.ignore_attr = "created_by, timestamp";
    values.no_fail = true;

    expect(buildCommandCli(spec!, values)).toBe(
      "openmc2donjon diff /runs/ref.h5 /runs/candidate.h5 --rtol 1e-5 " +
        "--ignore-attr created_by --ignore-attr timestamp --no-fail",
    );
  });

  it("builds PyGan writer comparison commands", () => {
    const spec = commandBuilderSpec("compare-writers");
    expect(spec).not.toBeNull();
    const values = defaultBuilderValues(spec!);
    expect(values.format).toBe("multicompo");
    values.input_h5 = "/runs/case/mgxs_library.h5";
    values.format = "multicompo";
    values.mixture = "ASM_01, ASM_02";
    values.summary_json = "/runs/case/writer_compare.json";
    values.keep_dir = "/runs/case/writer_compare_files";
    values.no_fail = true;

    expect(buildCommandCli(spec!, values)).toBe(
      "openmc2donjon compare-writers /runs/case/mgxs_library.h5 --format multicompo " +
        "--mixture ASM_01 --mixture ASM_02 --summary-json /runs/case/writer_compare.json " +
        "--keep-dir /runs/case/writer_compare_files --no-fail",
    );
  });

  it("keeps required placeholders visible for incomplete ADF driver commands", () => {
    const spec = commandBuilderSpec("make-homogeneous-face-flux");
    expect(spec).not.toBeNull();
    const cli = buildCommandCli(spec!, defaultBuilderValues(spec!));

    expect(cli).toContain("<mgxs_library.h5>");
    expect(cli).toContain("-o homogeneous_face_flux.h5");
    expect(cli).toContain("--volume-flux <volume_flux>");
    expect(cli).toContain("--net-current <net_current>");
  });

  it("builds OpenMC volume-flux export commands for SPH", () => {
    const spec = commandBuilderSpec("export-volume-flux");
    expect(spec).not.toBeNull();
    const values = defaultBuilderValues(spec!);
    values.statepoint = "/runs/case/mg_statepoint.h5";
    values.output = "/runs/case/openmc_mg_flux.h5";
    values.mgxs = "/runs/case/mgxs_library.h5";
    values.tally_name = "openmc_mg_volume_flux";
    values.dataset_name = "openmc_mg_flux";
    values.source_domain_ids = "91,17";
    values.summary_json = "/runs/case/openmc_mg_flux_summary.json";

    expect(buildCommandCli(spec!, values)).toBe(
      "openmc2donjon export-volume-flux /runs/case/mg_statepoint.h5 " +
        "-o /runs/case/openmc_mg_flux.h5 --mgxs /runs/case/mgxs_library.h5 " +
        "--tally-name openmc_mg_volume_flux --dataset-name openmc_mg_flux " +
        "--source-domain-ids 91,17 " +
        "--summary-json /runs/case/openmc_mg_flux_summary.json",
    );

    values.allow_zero_flux = true;

    expect(buildCommandCli(spec!, values)).toBe(
      "openmc2donjon export-volume-flux /runs/case/mg_statepoint.h5 " +
        "-o /runs/case/openmc_mg_flux.h5 --mgxs /runs/case/mgxs_library.h5 " +
        "--tally-name openmc_mg_volume_flux --dataset-name openmc_mg_flux " +
        "--source-domain-ids 91,17 " +
        "--allow-zero-flux " +
        "--summary-json /runs/case/openmc_mg_flux_summary.json",
    );

    const sourceIds = spec!.fields.find(
      (field) => field.name === "source_domain_ids",
    );
    expect(sourceIds?.required).not.toBe(true);
    expect(sourceIds?.help).toContain("canonical --mgxs mixture_names order");
    expect(sourceIds?.help).toContain("CE and MG geometries may use different");
    expect(sourceIds?.help).toContain(
      "derived from --mgxs source_domain_id metadata",
    );
  });

  it("builds SPH application commands for OpenMC MG reruns", () => {
    const spec = commandBuilderSpec("apply-sph");
    expect(spec).not.toBeNull();
    const values = defaultBuilderValues(spec!);
    values.input_h5 = "/runs/case/mg_case/mgxs_unapplied.h5";
    values.input_format = "openmc-mgxs";
    values.sph_source = "/runs/case/openmc_sph.h5";
    values.output = "/runs/case/mg_case/mgxs.h5";
    values.summary_json = "/runs/case/sph_apply_summary.json";
    values.force = true;

    expect(buildCommandCli(spec!, values)).toBe(
      "openmc2donjon apply-sph /runs/case/mg_case/mgxs_unapplied.h5 " +
        "--input-format openmc-mgxs --sph-source /runs/case/openmc_sph.h5 " +
        "-o /runs/case/mg_case/mgxs.h5 " +
        "--summary-json /runs/case/sph_apply_summary.json --force",
    );
  });

  it("defaults the SPH update table to production rate preservation", () => {
    const spec = commandBuilderSpec("make-sph-update-table");
    expect(spec).not.toBeNull();
    const values = defaultBuilderValues(spec!);

    expect(values.flux_normalization).toBe("auto");
    expect(values.sph_target).toBe("rate");
    expect(values.max_reference_flux_std_dev_rel).toBe("");
    expect(values.max_mg_flux_std_dev_rel).toBe("");
    expect(buildCommandCli(spec!, values)).toContain(
      "--flux-normalization auto --sph-target rate",
    );
    expect(buildCommandCli(spec!, values)).toContain(
      "--require-reference-flux-std-dev --max-reference-flux-std-dev-rel <CE_MAX_REL_STD_DEV>",
    );
    expect(buildCommandCli(spec!, values)).toContain(
      "--require-mg-flux-std-dev --max-mg-flux-std-dev-rel <MG_MAX_REL_STD_DEV>",
    );
    expect(builderCliIssues(spec!, values)).toEqual([
      "Physical route HOLD: CE max relative std dev is required when sph_target=rate.",
      "Physical route HOLD: MG max relative std dev is required when sph_target=rate.",
    ]);

    values.max_reference_flux_std_dev_rel = "0.025";
    values.max_mg_flux_std_dev_rel = "0.04";
    expect(buildCommandCli(spec!, values)).toContain(
      "--require-reference-flux-std-dev --max-reference-flux-std-dev-rel 0.025",
    );
    expect(buildCommandCli(spec!, values)).toContain(
      "--require-mg-flux-std-dev --max-mg-flux-std-dev-rel 0.04",
    );
    expect(builderCliIssues(spec!, values)).toEqual([]);

    values.flux_normalization = "none";
    values.sph_target = "flux";
    values.max_reference_flux_std_dev_rel = "";
    values.max_mg_flux_std_dev_rel = "";
    expect(buildCommandCli(spec!, values)).toContain(
      "--flux-normalization none --sph-target flux",
    );
    expect(buildCommandCli(spec!, values)).not.toContain(
      "--require-reference-flux-std-dev",
    );
    expect(buildCommandCli(spec!, values)).not.toContain(
      "--require-mg-flux-std-dev",
    );
    expect(builderCliIssues(spec!, values)).toEqual([]);

    values.max_reference_flux_std_dev_rel = "-0.01";
    expect(builderCliIssues(spec!, values)).toEqual([
      "CE max relative std dev must be a finite non-negative number.",
    ]);
  });

  it("builds the native DRAGON SPH physics validation command", () => {
    const spec = commandBuilderSpec("validate-native-sph");
    expect(spec).not.toBeNull();
    const values = defaultBuilderValues(spec!);
    values.reference_h5 = "/runs/case/reference.h5";
    values.reference_macrolib = "/runs/case/reference.macrolib.txt";
    values.sph_macrolib = "/runs/case/native_sph.macrolib.txt";
    values.verify_macrolib = "/runs/case/verify.macrolib.txt";
    values.result_listing = "/runs/case/donjon.result";
    values.execution_deck = "/runs/case/native_sph.x2m";
    values.energy_coverage = "/runs/case/energy_coverage.json";
    values.converter_receipt = "/runs/case/converter_receipt.json";
    values.summary_json = "/runs/case/physics_summary.json";

    expect(buildCommandCli(spec!, values)).toBe(
      "openmc2donjon validate-native-sph /runs/case/reference.h5 " +
        "--reference-macrolib /runs/case/reference.macrolib.txt " +
        "--sph-macrolib /runs/case/native_sph.macrolib.txt " +
        "--verify-macrolib /runs/case/verify.macrolib.txt " +
        "--result-listing /runs/case/donjon.result " +
        "--execution-deck /runs/case/native_sph.x2m " +
        "--energy-coverage /runs/case/energy_coverage.json " +
        "--converter-receipt /runs/case/converter_receipt.json " +
        "--summary-json /runs/case/physics_summary.json",
    );
    expect(commandBuilderStage("validate-native-sph").label).toBe(
      "Advanced · native DRAGON SPH",
    );
  });

  it("builds serve command with mock mode and repeated CORS origins", () => {
    const spec = commandBuilderSpec("serve");
    expect(spec).not.toBeNull();
    const values = defaultBuilderValues(spec!);
    values.host = "0.0.0.0";
    values.port = "8015";
    values.mock = true;
    values.cors_origin = "http://localhost:3000,http://127.0.0.1:3000";
    values.log_level = "DEBUG";

    expect(buildCommandCli(spec!, values)).toBe(
      "openmc2donjon serve --host 0.0.0.0 --port 8015 --mock " +
        "--cors-origin http://localhost:3000 --cors-origin http://127.0.0.1:3000 " +
        "--log-level DEBUG",
    );
  });

  it("prefills builder values from matching query parameters", () => {
    const spec = commandBuilderSpec("bundle");
    expect(spec).not.toBeNull();
    const values = builderValuesFromQuery(
      spec!,
      new URLSearchParams({
        mgxs: "/runs/case/mgxs_library.h5",
        mcompo: "/runs/case/out.mcompo.txt",
        output_dir: "/runs/case/bundle",
        force: "1",
        ignored: "nope",
      }),
    );

    expect(values.mgxs).toBe("/runs/case/mgxs_library.h5");
    expect(values.mcompo).toBe("/runs/case/out.mcompo.txt");
    expect(values.output_dir).toBe("/runs/case/bundle");
    expect(values.force).toBe(true);
    expect(buildCommandCli(spec!, values)).toContain(
      "--mgxs /runs/case/mgxs_library.h5 --mcompo /runs/case/out.mcompo.txt",
    );
  });

  it("labels SPH builders with the recommended CE/MG equivalence stage", () => {
    const stage = commandBuilderStage("export-volume-flux");

    expect(stage.label).toBe("Recommended OpenMC CE/MG SPH");
    expect(stage.summary).toContain("heterogeneous CE fine reference");
    expect(stage.summary).toContain("MG coarse model");
    expect(stage.reference).toContain("Different CE/MG geometries");
    expect(stage.reference).toContain("CE tallies use MG group boundaries");
    expect(stage.reference).toContain("state/BC/domain mapping");
  });

  it("emits the equals form for values that begin with a dash", () => {
    // Regression: `--mu-edges -1,-0.5,0.5,1` is rejected by argparse
    // (the leading dash classifies the value as an option string), and
    // physically valid mu bin edges always start at -1.
    const spec = commandBuilderSpec("export-surface-flux");
    expect(spec).not.toBeNull();
    const values = defaultBuilderValues(spec!);
    values.statepoint = "/runs/case/statepoint.h5";

    const cli = buildCommandCli(spec!, values);

    expect(cli).toContain("--mu-edges=-1,-0.5,0.5,1");
    expect(cli).not.toContain("--mu-edges -1");

    values.mu_edges = "-1,0,1";
    expect(buildCommandCli(spec!, values)).toContain("--mu-edges=-1,0,1");
  });

  it("flags argparse dependencies the doctor CLI enforces", () => {
    // Regression: the doctor form marked Statepoint freestanding-
    // optional, but the CLI rejects --statepoint without --recipe.
    const spec = commandBuilderSpec("doctor");
    expect(spec).not.toBeNull();

    const values = defaultBuilderValues(spec!);
    expect(builderCliIssues(spec!, values)).toEqual([]);

    values.statepoint = "/runs/case/statepoint.h5";
    expect(builderCliIssues(spec!, values)).toEqual([
      "Statepoint requires Recipe: the CLI rejects --statepoint without --recipe.",
    ]);

    values.recipe = "/runs/case/recipe.py";
    expect(builderCliIssues(spec!, values)).toEqual([]);

    const toggled = defaultBuilderValues(spec!);
    toggled.load_statepoint = true;
    expect(builderCliIssues(spec!, toggled)).toEqual([
      "Load statepoint requires Recipe: the CLI rejects --load-statepoint without --recipe.",
      "Load statepoint requires Statepoint: the CLI rejects --load-statepoint without --statepoint.",
    ]);
  });

  it("marks write-target paths as output-browse fields", () => {
    // The builder page keys its directory-select picker ("Browse for
    // output directory") off browse === "output"; input paths stay in
    // file-select mode.
    const spec = commandBuilderSpec("export-surface-flux");
    expect(spec).not.toBeNull();
    const field = (name: string) => spec!.fields.find((f) => f.name === name);

    expect(field("output")?.browse).toBe("output");
    expect(field("summary_json")?.browse).toBe("output");
    expect(field("statepoint")?.browse).toBe("file");
    expect(field("mgxs")?.browse).toBe("file");
  });
});
