"""Shared stage order; no scientific libraries are imported by the launcher."""
STAGES = ("preprocess", "fc", "edges", "geometry", "inference", "pairwise", "descriptive", "stability", "longitudinal", "combined", "statistics", "report")
MODULES = {name: "srcs.workflow." + name for name in STAGES}
MODULES["seedfc"] = "plugins.seedfc"

PRODUCTS = {
    "preprocess": {"preprocessed_recordings": ("*",)},
    "fc": {"fc_tables": ("*.pkl",), "batch_report": ("*.json",)},
    "edges": {"edge_store": ("*",)},
    "geometry": {"geometry_metrics": ("basis.json","distances.csv","points.csv","dispersion.csv","displacement.csv"),
                 "discriminability": ("subjects.csv",), "crossval_difference": ("crossval.csv",)},
    "inference": {"inference_summary": ("*",)},
    "pairwise": {"pairwise_tensors": ("*.npz",), "pairwise_selection": ("*.json","*.csv")},
    "descriptive": {"descriptive": ("*",)},
    "stability": {"stability": ("*",)},
    "longitudinal": {"longitudinal": ("*",)},
    "combined": {"combined_identification": ("*",), "identifiability_aggregate": ("@aggregate",)},
    "statistics": {"statistics_tables": ("*",)},
    "report": {"unified_report": ("*",)},
}


def validate_plan(stages):
    from importlib.util import find_spec
    seen = set()
    for stage in stages:
        if find_spec(MODULES[stage]) is None:
            raise ValueError(f"Missing implementation: {stage}")
        declared = PRODUCTS.get(stage,{})
        conflict = seen.intersection(declared)
        if conflict:
            raise ValueError(f"Multiple producers: {sorted(conflict)}")
        seen.update(declared)
