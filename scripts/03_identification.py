import numpy as np
import os
import json
from pathlib import Path
import pandas as pd

from srcs.utils.dataset import (select_complete_cohort, analysis_plan, load_week_tables, axis_frames)
from srcs.utils.config import load_context
from srcs.utils.utils import build_half_dataframe
from srcs.utils.results import save_json, save_pickle, save_table
from srcs.comp.ID_Tools import Identify, Edgewise_Contribution
from srcs.vis.PlotSeries import PlotSeries
from srcs.vis.circularGraph import circularGraph

import warnings
warnings.filterwarnings(
    "ignore",
    message="Mean of empty slice"
)


def select_complete_mice(ctx):
    """Compatibility entry point for the canonical week-based cohort selector."""
    return select_complete_cohort(ctx)


ctx = load_context()
logger = ctx.logger
analysis_mice, cohort_audit = select_complete_mice(ctx)
analysis_weeks = cohort_audit["retained_weeks"]
save_json(os.path.join(ctx.output_dir("Identification"), "cohort_selection.json"), cohort_audit)
if len(analysis_mice) < 2:
    raise ValueError("Fewer than two complete mice remain; see Identification/cohort_selection.json")
conditions = analysis_plan(ctx, analysis_weeks)
group_indices = ctx.group_indices(analysis_weeks, ctx.runs)
split_groups = ctx.split_groups(analysis_weeks, ctx.runs)
analysis_meta = dict(schema="week_mouse_run_v1", split_axis=ctx.split_axis,
                     selection_mode=ctx.selection_mode, conditions=conditions,
                     split_groups=split_groups, mice=analysis_mice,
                     weeks=analysis_weeks, runs=ctx.runs)
cohort_audit["analysis"] = analysis_meta
save_json(os.path.join(ctx.output_dir("Identification"), "cohort_selection.json"), cohort_audit)
ctx.save_analysis_metadata(analysis_meta)
logger.info("Cohort=%s; weeks=%s; analysis=%s", analysis_mice, analysis_weeks, analysis_meta)
id_cfg = ctx.section("ID")
edge_cfg = ctx.section("Edge")
session_seed = ctx.session_seed
# The shared repeat ID resolves against the selected split axis, after cover.
repeat_id = id_cfg.get("repeat_id")
split_ids = analysis_weeks if ctx.split_axis == "week" else ctx.runs
if isinstance(repeat_id, str):
    repeat_id = split_ids.index(repeat_id) + 1

render_options = ctx.render_options()
rel_box_preset = ctx.viz_style("stat_box", "stat_rel_effect_box")
group_indiv_box_preset = ctx.viz_style("stat_box", "stat_group_indiv_box")
edge_phi_preset = ctx.viz_style("edge_phi")
edge_dp_preset = ctx.viz_style("edge_dp")

df = pd.DataFrame(
    index=["_".join(atlas) for atlas in ctx.atlases],
    columns=list(conditions),
)

df.attrs["cohort_selection"] = cohort_audit
df.attrs["analysis"] = analysis_meta

for tag, selection in conditions.items():
    data = [[], []]
    relE = []
    labels = []

    selected_weeks, selected_runs = selection["weeks"], selection["runs"]
    summary_rows = []
    trial_path = ctx.output_dir("Identification", tag)

    for atlas in ctx.atlases:
        label = "_".join(atlas)

        logger.info(f"\n>>> {tag} Atlas {label} excluded:")

        tables = load_week_tables(ctx, label, analysis_mice, selected_weeks, selected_runs)
        frames = axis_frames(tables, ctx.split_axis, analysis_mice, selected_weeks, selected_runs)
        df_out = build_half_dataframe(
            list(frames.values()), 0, ctx.group_mode,
            repeat_session=repeat_id, seed=session_seed,
            mouse_ids=analysis_mice,
        )
        logger.info("Split groups: %s", split_groups)

        # Identify
        ID = Identify(df_out, atlas_exclude=atlas)
        info, sim_mat, Indiv, Group, relEft, cohen_d = ID.ID_Test()
        stat = ID.statistics()

        shuffle_result = None
        if id_cfg["shuffle"]:
            shuffle_result = ID.shuffle(
                iter=id_cfg["n_iter"],
                seed=session_seed,
                tag=tag,
                save_path=ctx.output_dir("Identification", tag, "Similarity"),
                show=render_options["show"],
                log=render_options["log"],
                plot_config=ctx.viz_style("id_shuffle_hist"),
            )

        ID.output(
            save_path=ctx.output_dir("Identification", tag, "Similarity"),
            dpi=ctx.figure_dpi(),
            tag=tag,
            show=render_options["show"],
            log=render_options["log"],
            paired_config=ctx.viz_style("id_paired"),
            matrix_config=ctx.viz_style("id_similarity_matrix"),
        )

        numeric_path = ctx.output_dir("Identification", tag, "Numeric")
        save_pickle(os.path.join(numeric_path, f"identification_{label}.pkl"), {
            "similarity": sim_mat, "group": Group, "individual": Indiv,
            "relative_effect": relEft, "cohen_d": cohen_d, "statistics": stat,
            "mice": analysis_mice, "roi_labels": info[2],
            "weeks": selected_weeks, "runs": selected_runs, "analysis": analysis_meta,
            "split_groups": split_groups,
        })
        save_table(os.path.join(numeric_path, f"similarity_{label}.csv"),
                   pd.DataFrame(sim_mat, index=analysis_mice, columns=analysis_mice))
        save_table(os.path.join(numeric_path, f"subjects_{label}.csv"),
                   pd.DataFrame({"group": Group, "individual": Indiv, "relative_effect": relEft},
                                index=pd.Index(analysis_mice, name="mouse")))
        save_json(os.path.join(numeric_path, f"statistics_{label}.json"), {**stat, "cohen_d": cohen_d})
        summary_rows.append({"atlas": label, "condition": tag, "split_axis": ctx.split_axis, "selection_mode": ctx.selection_mode, "week_count": len(selected_weeks), "run_count": len(selected_runs), **stat, "cohen_d": cohen_d})
        if shuffle_result is not None:
            shuffle_result["analysis_selection"] = selection
            shuffle_result["split_axis"] = ctx.split_axis
            shuffle_result["split_groups"] = split_groups
            save_pickle(os.path.join(numeric_path, f"shuffle_{label}.pkl"), shuffle_result)
            save_json(os.path.join(numeric_path, f"shuffle_{label}.json"), shuffle_result)

        # Edgewise Contribution
        data[0].append(Group)
        data[1].append(Indiv)
        relE.append(relEft)
        df.at["_".join(atlas), tag] = [Group, Indiv, relEft, stat, cohen_d]

        # Edge holds each selected opposite-axis sample fixed and splits the same axis as ID.
        edge_axis = "week" if ctx.split_axis == "run" else "run"
        edge_frames = axis_frames(tables, edge_axis, analysis_mice, selected_weeks, selected_runs)
        EC = Edgewise_Contribution(
            path=ctx.correlation_path, mouse_ids=analysis_mice,
            sessions=list(edge_frames), seed=edge_cfg["seed"], atlas_exclude=atlas,
            frames=edge_frames, group_indices=group_indices,
        )
        Phi, DP = EC.run(repeat=edge_cfg["repeat"])

        save_path = ctx.output_dir("Identification", tag, "EdgeContri")

        save_pickle(os.path.join(save_path, f"edges_{label}_{tag}.pkl"), {
            "Phi": Phi, "DP": DP, "roi_labels": EC.roi_labels, "mice": analysis_mice,
            "weeks": selected_weeks, "runs": selected_runs,
            "split_axis": ctx.split_axis, "split_groups": split_groups,
            "pairing": "shared_ID_groups", "per_sample_axis": edge_axis,
            "effective_repeats": 1,
            "status": "complete" if Phi and DP else "unavailable",
            "reason": None if Phi and DP else "At least two samples on the split axis and two complete mice are required.",
        })
        if not Phi or not DP:
            logger.warning(f"{tag}: Edge unavailable (needs >=2 samples on split axis and >=2 mice).")
        else:
            # Group Consistency (Phi)
            avg_phi = np.mean(list(Phi.values()), axis=0)
            phi_mat = np.zeros((len(EC.roi_labels), len(EC.roi_labels)))
            phi_mat[EC.edge_index] = avg_phi
            phi_mat += phi_mat.T

            circularGraph(
                phi_mat, [lab.split("-")[1] for lab in EC.roi_labels],
                **edge_phi_preset,
                node_color_spec=EC.roi_labels,
                threshold=np.percentile(phi_mat, edge_cfg["threshold"]),
                gap=edge_cfg["gap"],
                cross_only=edge_cfg["cross_only"],
                save_path=os.path.join(save_path, f"EdgeContri_phi_{label}_excluded_thres{edge_cfg['threshold']}_{tag}"),
                save_mat=edge_cfg["save_mat"]
            )

            save_table(os.path.join(save_path, f"phi_{label}_{tag}.csv"),
                       pd.DataFrame(phi_mat, index=EC.roi_labels, columns=EC.roi_labels))
            # Differential Power (DP)
            avg_dp = np.mean(list(DP.values()), axis=0)
            dp_mat = np.zeros((len(EC.roi_labels), len(EC.roi_labels)))
            dp_mat[EC.edge_index] = avg_dp
            dp_mat += dp_mat.T

            circularGraph(
                dp_mat, [lab.split("-")[1] for lab in EC.roi_labels],
                **edge_dp_preset,
                node_color_spec=EC.roi_labels,
                threshold=np.percentile(dp_mat, edge_cfg["threshold"]),
                gap=edge_cfg["gap"],
                cross_only=edge_cfg["cross_only"],
                save_path=os.path.join(save_path, f"EdgeContri_dp_{label}_excluded_thres{edge_cfg['threshold']}_{tag}"),
                save_mat=edge_cfg["save_mat"]
            )

            save_table(os.path.join(save_path, f"dp_{label}_{tag}.csv"),
                       pd.DataFrame(dp_mat, index=EC.roi_labels, columns=EC.roi_labels))

        # labels
        if label == "SM_VIS_AUD":
            label = "Association"
        elif label in ("PF_LAT_MED", "PF_LAT_MED_AUDv"):
            label = "Sensory excluding AUDv" if "AUDv" in label else "Sensory"
        elif label == "None":
            label = "Whole Brain"
        elif label == "AUDv":
            label = "Whole Brain excluding AUDv"
        labels.append(label)
    
    rel_box_cfg = dict(rel_box_preset)
    rel_box = PlotSeries(
        figsize=tuple(rel_box_cfg.pop("figsize", (4, 4))),
        dpi=rel_box_cfg.pop("dpi", ctx.figure_dpi()),
        save_path=ctx.output_dir("Identification", tag, "Stat"),
    )
    rel_box.box(
        relE,
        name=f"boxmap_relEffect_{tag}",
        xticklabel=labels,
        **rel_box_cfg,
    ).render(**render_options)

    group_indiv_box_cfg = dict(group_indiv_box_preset)
    group_indiv_box = PlotSeries(
        figsize=tuple(group_indiv_box_cfg.pop("figsize", (4, 4))),
        dpi=group_indiv_box_cfg.pop("dpi", ctx.figure_dpi()),
        save_path=ctx.output_dir("Identification", tag, "Stat"),
    )
    group_indiv_box.box(
        data,
        name=f"boxmap_Group-Indiv_{tag}",
        xticklabel=labels,
        **group_indiv_box_cfg,
    ).render(**render_options)

    save_pickle(os.path.join(trial_path, ctx.statistics_file), df[[tag]])
    save_table(os.path.join(trial_path, "statistics_summary.csv"), pd.DataFrame(summary_rows), index=False)

df_save_path = ctx.statistics_path
save_pickle(df_save_path, df)
logger.info("\nStatistics done output.")
