"""Converter vector names shared with pre-conversion transformations."""

INVERSE_VELOCITY_DATASETS = ("inverse_velocity", "inverse-velocity", "OVERV", "overv")
H_FACTOR_DATASETS = (
    "h_factor", "H-FACTOR", "H_FACTOR", "kappa_fission",
    "kappa_fission_xs", "kappa_fission_cross_section",
)

# Collapse writes one canonical dataset per physical quantity. These names
# retain the native exporter convention; no unit conversion is implied.
COLLAPSE_VECTOR_FIELDS = {
    name: (name,) for name in (
        "total", "transport_total", "absorption", "reduced_absorption",
        "fission", "nu_fission",
    )
} | {
    "kappa_fission": H_FACTOR_DATASETS,
    "inverse_velocity": INVERSE_VELOCITY_DATASETS,
}
