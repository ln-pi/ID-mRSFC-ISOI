"""Run-level amplitude/topology and same-day repeatability; no ICC/model fitting."""
from pathlib import Path

import pandas as pd

from srcs.analysis.connectome import save_arrays
from srcs.analysis.descriptive import describe
from srcs.utils.results import save_json, save_table
from srcs.utils.edge_store import load_store, store_path
from srcs.utils.provenance import sha256


def run(ctx):
    reliability, deviations = [], []
    baseline = ctx.section("descriptive").get("reference_week", "03w")
    for atlas in ["_".join(a) for a in ctx.atlases]:
        store = load_store(ctx, atlas)
        result = describe(store, baseline)
        reliability.extend(result["reliability"])
        deviations.extend(result["deviations"])
        variance, delta = result["edge_variance"], result["delta_from_reference"]
        save_arrays(Path(ctx.output_dir("results", "descriptive")) / f"topology_{atlas}.npz",
            dict(schema="descriptive_topology_v1", atlas=atlas, mouse_ids=store.mice,
                 week_ids=store.weeks, run_ids=store.runs, roi_labels=store.roi_labels,
                 edge_order="numpy.triu_indices(k=1)", scale="FC_edge_Fisher_z",
                 reference_week=baseline, reference_meaning="post_injection_reference_not_healthy_baseline",
                 edge_variance_ddof=1, code=store.metadata["code"],
                 edge_store=dict(path=str(store_path(ctx, atlas)), sha256=sha256(store_path(ctx, atlas)))),
            edge_variance=variance, delta_from_reference=delta, edge_index=store.edge_index,
            n_valid_mice=store.valid_record.sum(axis=0), valid_record=store.valid_record)
    folder = Path(ctx.output_dir("results", "descriptive"))
    save_table(folder / "reliability.csv", pd.DataFrame(reliability, columns=[
        "atlas", "week", "mouse", "run_1", "run_2", "pattern_r", "pattern_status", "rms_difference", "mean_signed_difference"]), index=False)
    save_table(folder / "deviations.csv", pd.DataFrame(deviations, columns=[
        "atlas", "mouse", "week", "run", "n_reference_mice", "reference_mice", "rms_from_group_loo",
        "mean_signed_group_deviation", "reference_week", "reference_status", "rms_from_reference", "mean_signed_reference_change"]), index=False)
    save_json(folder / "methods.json", dict(scale="FC_edge_Fisher_z", reference_week=baseline,
        group_reference="same_cohort_same_week_same_run_leave_one_mouse_out",
        inference="descriptive_only", icc="not_fitted", variance_model="not_fitted",
        missing="no_imputation_no_single_run_substitution"))
