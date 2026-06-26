import numpy as np
import os
import pandas as pd

from srcs.utils.config import load_context
from srcs.utils.utils import build_half_dataframe
from srcs.comp.ID_Tools import Identify, Edgewise_Contribution
from srcs.vis.PlotSeries import PlotSeries
from srcs.vis.circularGraph import circularGraph

import warnings
warnings.filterwarnings(
    "ignore",
    message="Mean of empty slice"
)


ctx = load_context()
logger = ctx.logger
id_cfg = ctx.section("ID")
edge_cfg = ctx.section("Edge")
render_options = ctx.render_options()
rel_box_preset = ctx.viz_style("stat_box", "stat_rel_effect_box")
group_indiv_box_preset = ctx.viz_style("stat_box", "stat_group_indiv_box")
edge_phi_preset = ctx.viz_style("edge_phi")
edge_dp_preset = ctx.viz_style("edge_dp")

df = pd.DataFrame(
    index=["_".join(atlas) for atlas in ctx.atlases],
    columns=[f"{n}Trial" for n in range(1, len(ctx.trials) + 1)],
)

for trial_cnt in range(len(ctx.trials)):
    data = [[], []]
    relE = []
    labels = []

    trial_cnt += 1
    tag = f"{trial_cnt}Trial"

    for atlas in ctx.atlases:
        label = "_".join(atlas)

        logger.info(f"\n>>> {tag} Atlas {label} excluded:")

        session_paths = [
            os.path.join(ctx.processed_path, f"{session}_correlation_matrices_{label}_excluded.pkl")
            for session in ctx.sessions
        ]
        missing_paths = [path for path in session_paths if not os.path.exists(path)]
        if missing_paths:
            logger.raise_(FileNotFoundError(
                "Missing correlation matrix files from config sessions: "
                + ", ".join(missing_paths)
            ))

        df_out = build_half_dataframe(session_paths, trial_cnt, id_cfg["group_mode"])

        # Identify
        ID = Identify(df_out, atlas_exclude=atlas)
        info, sim_mat, Indiv, Group, relEft, cohen_d = ID.ID_Test()
        stat = ID.statistics()

        if id_cfg["shuffle"]:
            ID.shuffle(
                iter=id_cfg["n_iter"],
                save_path=ctx.output_dir("Similarity"),
                show=render_options["show"],
                log=render_options["log"],
                plot_config=ctx.viz_style("id_shuffle_hist"),
            )

        ID.output(
            save_path=ctx.output_dir("Similarity"),
            dpi=ctx.figure_dpi(),
            tag=tag,
            show=render_options["show"],
            log=render_options["log"],
            paired_config=ctx.viz_style("id_paired"),
            matrix_config=ctx.viz_style("id_similarity_matrix"),
        )

        # Edgewise Contribution
        data[0].append(Group)
        data[1].append(Indiv)
        relE.append(relEft)
        df.at["_".join(atlas), tag] = [Group, Indiv, relEft, stat, cohen_d]

        EC = Edgewise_Contribution(
            path=ctx.processed_path,
            mouse_ids=ctx.mice,
            sessions=ctx.sessions,
            seed=edge_cfg["seed"],
            atlas_exclude=atlas,
        )

        Phi, DP = EC.run(repeat=edge_cfg["repeat"])

        save_path = ctx.output_dir("EdgeContri")

        # Group Consistency (Phi)
        avg_phi = np.mean(list(Phi.values()), axis=0)
        phi_mat = np.zeros((len(info[2]), len(info[2])))
        phi_mat[EC.edge_index] = avg_phi
        phi_mat += phi_mat.T

        circularGraph(
            phi_mat, [lab.split("-")[1] for lab in info[2]],
            **edge_phi_preset,
            node_color_spec=info[2],
            threshold=np.percentile(phi_mat, edge_cfg["threshold"]),
            gap=edge_cfg["gap"],
            cross_only=edge_cfg["cross_only"],
            save_path=os.path.join(save_path, f"EdgeContri_phi_{label}_excluded_thres{edge_cfg['threshold']}_{tag}"),
            save_mat=edge_cfg["save_mat"]
        )

        # Differential Power (DP)
        avg_dp = np.mean(list(DP.values()), axis=0)
        dp_mat = np.zeros((len(info[2]), len(info[2])))
        dp_mat[EC.edge_index] = avg_dp
        dp_mat += dp_mat.T

        circularGraph(
            dp_mat, [lab.split("-")[1] for lab in info[2]],
            **edge_dp_preset,
            node_color_spec=info[2],
            threshold=np.percentile(dp_mat, edge_cfg["threshold"]),
            gap=edge_cfg["gap"],
            cross_only=edge_cfg["cross_only"],
            save_path=os.path.join(save_path, f"EdgeContri_dp_{label}_excluded_thres{edge_cfg['threshold']}_{tag}"),
            save_mat=edge_cfg["save_mat"]
        )

        # labels
        if label == "SM_VIS_AUD":
            label = "Association"
        elif label == "PF_LAT_MED":
            label = "Sensory"
        elif label == "None":
            label = "Whole Brain"
        elif label == "AUDv":
            label = "AUDv excluded"
        labels.append(label)
    
    rel_box_cfg = dict(rel_box_preset)
    rel_box = PlotSeries(
        figsize=tuple(rel_box_cfg.pop("figsize", (4, 4))),
        dpi=rel_box_cfg.pop("dpi", ctx.figure_dpi()),
        save_path=ctx.output_dir("Stat"),
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
        save_path=ctx.output_dir("Stat"),
    )
    group_indiv_box.box(
        data,
        name=f"boxmap_Group-Indiv_{tag}",
        xticklabel=labels,
        **group_indiv_box_cfg,
    ).render(**render_options)

df_save_path = os.path.join(ctx.output_path, "statistics.pkl")
df.to_pickle(df_save_path)
logger.info("\nStatistics done output.")
