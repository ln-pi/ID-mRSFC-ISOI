"""Persist combined numerical products; all rendering belongs to report."""
from pathlib import Path

import numpy as np
import pandas as pd

from srcs.analysis.edge_contribution import edge_contributions, edge_matrix
from srcs.analysis.identification import identify, permutation_test
from srcs.utils.dataset import select_complete_cohort, analysis_plan
from srcs.utils.edge_store import load_store, store_path
from srcs.utils.provenance import sha256
from srcs.utils.results import save_json, save_pickle, save_table


def run(ctx):
    mice, audit = select_complete_cohort(ctx)
    weeks = audit["retained_weeks"]
    metadata_dir = Path(ctx.output_dir("Metadata"))
    save_json(metadata_dir / "cohort_selection.json", audit)
    if len(mice) < 2:
        raise ValueError("Fewer than two complete mice remain; see cohort_selection.json")
    conditions = analysis_plan(ctx, weeks)
    groups = ctx.group_indices(weeks, ctx.runs)
    split_groups = ctx.split_groups(weeks, ctx.runs)
    stores = {"_".join(a): load_store(ctx, "_".join(a)) for a in ctx.atlases}
    metadata = dict(schema="week_mouse_run_v1", split_axis=ctx.split_axis,
        selection_mode=ctx.selection_mode, conditions=conditions, split_groups=split_groups,
        mice=mice, weeks=weeks, runs=ctx.runs, averaging="mean_F_z_then_Pearson_never_mean_S",
        similarity_scale="Fisher_z_capped_at_4_for_P1_compatibility", disjoint_halves=True,
        split_sizes=[len(g) for g in groups],
        edge_store_sha256={a: sha256(store_path(ctx, a)) for a in stores},
        code=next(iter(stores.values())).metadata["code"])
    audit["analysis"] = metadata
    save_json(metadata_dir / "cohort_selection.json", audit)
    ctx.save_analysis_metadata(metadata)
    aggregate = pd.DataFrame(index=list(stores), columns=list(conditions), dtype=object)
    aggregate.attrs.update(cohort_selection=audit, analysis=metadata)
    for condition, selection in conditions.items():
        selected_weeks, selected_runs = selection["weeks"], selection["runs"]
        summaries = []
        for atlas, store in stores.items():
            result = identify(*store.average_halves(mice, selected_weeks, selected_runs, ctx.split_axis, groups))
            result.update(mice=mice, roi_labels=store.roi_labels, weeks=selected_weeks,
                          runs=selected_runs, analysis=metadata, split_groups=split_groups)
            numeric = Path(ctx.output_dir("Identification", condition, "Numeric"))
            save_pickle(numeric / f"identification_{atlas}.pkl", result)
            save_table(numeric / f"similarity_{atlas}.csv", pd.DataFrame(result["similarity"], index=mice, columns=mice))
            save_table(numeric / f"subjects_{atlas}.csv", pd.DataFrame(
                {k: result[k] for k in ("group", "individual", "relative_effect")}, index=pd.Index(mice, name="mouse")))
            save_json(numeric / f"statistics_{atlas}.json", {**result["statistics"], "cohen_d": result["cohen_d"]})
            summaries.append(dict(atlas=atlas, condition=condition, split_axis=ctx.split_axis,
                selection_mode=ctx.selection_mode, week_count=len(selected_weeks), run_count=len(selected_runs),
                **result["statistics"]))
            aggregate.at[atlas, condition] = [result["group"], result["individual"], result["relative_effect"],
                                             result["statistics"], result["cohen_d"]]
            if ctx.section("ID").get("shuffle", True):
                shuffled = permutation_test(result["similarity"], ctx.section("ID").get("n_iter", 1000), ctx.session_seed)
                shuffled.update(analysis_selection=selection, split_axis=ctx.split_axis, split_groups=split_groups)
                save_pickle(numeric / f"shuffle_{atlas}.pkl", shuffled)
                save_json(numeric / f"shuffle_{atlas}.json", shuffled)
            opposite_axis = "week" if ctx.split_axis == "run" else "run"
            opposite_ids = selected_weeks if ctx.split_axis == "run" else selected_runs
            phi, dp = edge_contributions(store.subset(mice, selected_weeks, selected_runs, scale="r"),
                                         ctx.split_axis, groups, opposite_ids)
            edges = Path(ctx.output_dir("Identification", condition, "edges"))
            save_pickle(edges / f"edges_{atlas}_{condition}.pkl", dict(Phi=phi, DP=dp,
                roi_labels=store.roi_labels, edge_index=store.edge_index, mice=mice,
                weeks=selected_weeks, runs=selected_runs, split_axis=ctx.split_axis, split_groups=split_groups,
                pairing="shared_ID_groups", per_sample_axis=opposite_axis, effective_repeats=1,
                status="complete", reason=None))
            for name, values in (("phi", phi), ("dp", dp)):
                matrix = edge_matrix(np.mean(list(values.values()), axis=0), store.edge_index, len(store.roi_labels))
                save_table(edges / f"{name}_{atlas}_{condition}.csv",
                           pd.DataFrame(matrix, index=store.roi_labels, columns=store.roi_labels))
            ctx.logger.info("Combined %s/%s: %s mice; accuracy=%.3f", atlas, condition,
                            len(mice), result["statistics"]["identification_accuracy"])
        folder = Path(ctx.output_dir("Identification", condition))
        save_pickle(folder / ctx.statistics_file, aggregate[[condition]])
        save_table(folder / "statistics_summary.csv", pd.DataFrame(summaries), index=False)
    save_pickle(ctx.statistics_path, aggregate)
