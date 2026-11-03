"""Shared stage order; no scientific libraries are imported by the launcher."""
STAGES = ("preprocess", "fc", "edges", "geometry", "inference", "pairwise", "descriptive", "stability", "longitudinal", "combined", "statistics", "report")
MODULES = {name: "srcs.workflow." + name for name in STAGES}
MODULES["seedfc"] = "plugins.seedfc"

PRODUCTS = {
    "preprocess": {"preprocessed_recordings_v1": ("*",)},
    "fc": {"fc_tables_v1": ("*.pkl",), "batch_report_v1": ("*.json",)},
    "edges": {"edge_store_v1": ("*",)},
    "geometry": {"geometry_metrics_v1": ("basis.json","distances.csv","points.csv","dispersion.csv","displacement.csv"),
                 "discriminability_v1": ("subjects.csv",), "crossval_difference_v1": ("crossval.csv",)},
    "inference": {"inference_summary_v1": ("*",)},
    "pairwise": {"pairwise_tensors_v1": ("*.npz",), "pairwise_selection_v1": ("*.json","*.csv")},
    "descriptive": {"descriptive_v2": ("*",)},
    "stability": {"stability_v1": ("*",)},
    "longitudinal": {"longitudinal_v1": ("*",)},
    "combined": {"combined_identification_v1": ("*",), "identifiability_aggregate_v1": ("@aggregate",)},
    "statistics": {"statistics_tables_v1": ("*",)},
    "report": {"unified_report_v1": ("*",)},
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
