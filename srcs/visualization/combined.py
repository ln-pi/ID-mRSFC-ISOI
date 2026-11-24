"""Draw combined figures exclusively from persisted numerical products."""
from pathlib import Path

import numpy as np
import pandas as pd

from srcs.utils.dataset import statistics_metadata
from srcs.visualization.plot_series import PlotSeries
from srcs.visualization.circular_graph import circularGraph
from srcs.visualization.scales import resolve

LABELS = {"SM_VIS_AUD": "Association", "PF_LAT_MED_AUDv": "Sensory excluding AUDv",
          "AUDv": "Whole Brain excluding AUDv"}


def plot(ctx, data, kind, name, folder, config):
    config = dict(config)
    series = PlotSeries(figsize=tuple(config.pop("figsize", (6, 5))),
                        dpi=config.pop("dpi", ctx.figure_dpi()), save_path=str(folder))
    series.add_plot(data, plot_type=kind, name=name, config=config).render(**ctx.render_options())


def render_identification(ctx, atlas, condition, result, folder):
    plot(ctx, [result["group"], result["individual"]], "paired_vertical",
         f"IDtest_paired_{atlas}_excluded_{condition}", folder, ctx.viz_style("id_paired"))
    style = ctx.viz_style("id_similarity_matrix")
    viz = ctx.viz()
    legacy_limits = style.pop("limits_by_atlas", {}).get(atlas, style.pop("default_limits", None))
    explicit_policy = ctx.config.get("visualization", {}).get("metrics", {}).get("id_similarity_matrix")
    scale = resolve(result["similarity"], viz, "id_similarity_matrix",
                    fixed=None if explicit_policy else legacy_limits)
    if scale is None:
        return
    tick_size = style.pop("tick_fontsize", 10)
    mice = result["mice"]
    style.update(**scale.image, cbar_extend=scale.extend, xlabel="2nd half", ylabel="1st half",
        xticks=dict(ticks=np.arange(len(mice)), labels=mice, rotation=90, fontsize=tick_size),
        yticks=dict(ticks=np.arange(len(mice)), labels=mice, rotation=0, fontsize=tick_size))
    plot(ctx, result["similarity"], "origin", f"IDtest_simMat_{atlas}_excluded_{condition}", folder, style)


def run(ctx):
    stat = pd.read_pickle(ctx.statistics_path)
    metadata = statistics_metadata(ctx, stat)
    files = [str(ctx.statistics_path)]
    enabled = ctx.render_options()["save"] or ctx.render_options()["show"]
    for condition in stat.columns:
        group_data, indiv_data, relative_data, labels = [], [], [], []
        for atlas in stat.index:
            numeric = Path(ctx.output_dir("Identification", str(condition), "Numeric", create=False))
            path = numeric / f"identification_{atlas}.pkl"
            result = pd.read_pickle(path)
            if result["analysis"] != metadata or result["mice"] != metadata["mice"]:
                raise ValueError(f"Combined result metadata mismatch: {path}")
            for key, position in (("group", 0), ("individual", 1), ("relative_effect", 2)):
                np.testing.assert_allclose(result[key], stat.at[atlas, condition][position], equal_nan=True)
            files.append(str(path))
            group_data.append(result["group"])
            indiv_data.append(result["individual"])
            relative_data.append(result["relative_effect"])
            labels.append(LABELS.get(atlas, atlas))
            folder = Path(ctx.output_dir("IdentificationFigures", str(condition), "Similarity"))
            if enabled:
                render_identification(ctx, atlas, str(condition), result, folder)
            shuffled = numeric / f"shuffle_{atlas}.pkl"
            if shuffled.exists():
                files.append(str(shuffled))
                payload = pd.read_pickle(shuffled)
                if enabled:
                    style = ctx.viz_style("id_shuffle_hist")
                    style["arrow_x"] = payload["real_accuracy"][0]
                    plot(ctx, payload["shuffle_acc"][0], "hist",
                         f"ID_shuffle_acc_{atlas}_excluded_{condition}", folder, style)
            edges = Path(ctx.output_dir("Identification", str(condition), "edges", create=False))
            for metric in ("phi", "dp"):
                path = edges / f"{metric}_{atlas}_{condition}.csv"
                files.append(str(path))
                frame = pd.read_csv(path, index_col=0)
                if list(frame.index) != result["roi_labels"] or list(frame.columns) != result["roi_labels"]:
                    raise ValueError(f"Edge ROI order mismatch: {path}")
                if enabled:
                    options = ctx.viz("edge_graph",ctx.section("Edge"))
                    target = Path(ctx.output_dir("IdentificationFigures", str(condition), "edges"))
                    matrix = frame.to_numpy()
                    circularGraph(matrix, [label.split("-", 1)[-1] for label in result["roi_labels"]],
                        visualization=ctx.viz(), **ctx.viz_style(f"edge_{metric}"), node_color_spec=result["roi_labels"],
                        threshold=np.percentile(matrix, options.get("threshold", 99)), gap=options.get("gap", .25),
                        cross_only=options.get("cross_only", False), show=ctx.render_options()["show"],
                        save_path=str(target / f"EdgeContri_{metric}_{atlas}_excluded_thres{options.get('threshold',99)}_{condition}")
                            if ctx.render_options()["save"] else None, save_mat=options.get("save_mat", True))
        if enabled:
            folder = Path(ctx.output_dir("IdentificationFigures", str(condition), "Stat"))
            relative = ctx.viz_style("stat_box", "stat_rel_effect_box")
            relative["xticklabel"] = labels
            plot(ctx, relative_data, "box", f"boxmap_relEffect_{condition}", folder, relative)
            combined = ctx.viz_style("stat_box", "stat_group_indiv_box")
            combined["xticklabel"] = labels
            plot(ctx, [group_data, indiv_data], "box", f"boxmap_Group-Indiv_{condition}", folder, combined)
    return files
