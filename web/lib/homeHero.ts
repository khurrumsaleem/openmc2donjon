export interface HomeHero {
  kicker: string;
  heading: string;
  paragraph: string;
  supporting: string;
}

export const HOME_EQUIVALENCE_FLOW = [
  {
    label: "OpenMC CE fine + MG coarse",
    qualifier: "different geometries · matched physics contract",
  },
  {
    label: "Rate-preserving SPH",
    qualifier: "recommended when equivalence is required",
  },
  {
    label: "Corrected / Converter-ready HDF5",
    qualifier: "Converter-ready input",
  },
] as const;

export const HOME_HANDOFF_FLOW = [
  {
    label: "Converter-ready HDF5",
    qualifier: "direct users can start here",
  },
  {
    label: "Converter",
    qualifier: "required formal handoff boundary",
  },
  {
    label: "L_MULTICOMPO / L_MACROLIB",
    qualifier: "checked object + hash-linked receipt",
  },
  {
    label: "DONJON",
    qualifier: "downstream model and acceptance",
  },
] as const;

export const HOME_HERO: HomeHero = {
  kicker: "OpenMC → DRAGON / DONJON handoff",
  heading: "Convert OpenMC MGXS into a traceable DRAGON/DONJON object.",
  paragraph:
    "Converter is the formal handoff boundary after the HDF5 is ready. It checks the MGXS contract and declared mixture/domain ordering, then writes an L_MULTICOMPO or L_MACROLIB object with a hash-linked receipt. The conservative CE-to-MG map remains upstream SPH evidence.",
  supporting:
    "When fine-to-coarse equivalence is required, use the heterogeneous OpenMC CE reference and homogenized OpenMC MG coarse model to converge rate-preserving SPH, apply it to create the corrected HDF5, and only then enter Converter. Native DRAGON SPH remains available as an advanced, project-specific route. Built-in ASCII is the default writer; PyGan/LCM is optional.",
} as const;
