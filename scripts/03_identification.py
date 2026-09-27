import numpy as np
import os
import pandas as pd

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


ctx = load_context()
logger = ctx.logger
id_cfg = ctx.section("ID")
edge_cfg = ctx.section("Edge")
session_seed = ctx.session_seed
repeat_session = ctx.repeat_session
ctx.save_analysis_metadata()
logger.info(f"Analysis output: {ctx.output_path}")
logger.info(f"Session groups: {ctx.session_groups}")
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
    summary_rows = []
    trial_path = ctx.output_dir(tag)

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

        df_out = build_half_dataframe(
            session_paths, trial_cnt, ctx.group_mode,
            repeat_session=repeat_session, seed=session_seed,
            mouse_ids=ctx.mice, trial_ids=ctx.trials,
        )
        logger.info(f"Session groups: {df_out.attrs['session_groups']}")

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
                save_path=ctx.output_dir(tag, "Similarity"),
                show=render_options["show"],
                log=render_options["log"],
                plot_config=ctx.viz_style("id_shuffle_hist"),
            )

        ID.output(
            save_path=ctx.output_dir(tag, "Similarity"),
            dpi=ctx.figure_dpi(),
            tag=tag,
            show=render_options["show"],
            log=render_options["log"],
            paired_config=ctx.viz_style("id_paired"),
            matrix_config=ctx.viz_style("id_similarity_matrix"),
        )

        numeric_path = ctx.output_dir(tag, "Numeric")
        save_pickle(os.path.join(numeric_path, f"identification_{label}.pkl"), {
            "similarity": sim_mat, "group": Group, "individual": Indiv,
            "relative_effect": relEft, "cohen_d": cohen_d, "statistics": stat,
            "mice": ctx.mice, "roi_labels": info[2],
            "trial_ids": ctx.trials[:trial_cnt], "session_groups": ctx.session_groups,
        })
        save_table(os.path.join(numeric_path, f"similarity_{label}.csv"),
                   pd.DataFrame(sim_mat, index=ctx.mice, columns=ctx.mice))
        save_table(os.path.join(numeric_path, f"subjects_{label}.csv"),
                   pd.DataFrame({"group": Group, "individual": Indiv, "relative_effect": relEft},
                                index=pd.Index(ctx.mice, name="mouse")))
        save_json(os.path.join(numeric_path, f"statistics_{label}.json"), {**stat, "cohen_d": cohen_d})
        summary_rows.append({"atlas": label, "trial_count": trial_cnt, **stat, "cohen_d": cohen_d})
        if shuffle_result is not None:
            save_pickle(os.path.join(numeric_path, f"shuffle_{label}.pkl"), shuffle_result)
            save_json(os.path.join(numeric_path, f"shuffle_{label}.json"), shuffle_result)

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
            trial_cnt=trial_cnt, trial_ids=ctx.trials,
        )

        Phi, DP = EC.run(repeat=edge_cfg["repeat"])

        save_path = ctx.output_dir(tag, "EdgeContri")

        save_pickle(os.path.join(save_path, f"edges_{label}_{tag}.pkl"), {
            "Phi": Phi, "DP": DP, "roi_labels": EC.roi_labels,
            "trial_ids": ctx.trials[:trial_cnt],
            "status": "complete" if Phi and DP else "unavailable",
            "reason": None if Phi and DP else "At least two valid trials and two mice per session are required.",
        })
        if not Phi or not DP:
            logger.warning(f"{tag}: Edge unavailable (needs >=2 trials and >=2 mice/session).")
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
        save_path=ctx.output_dir(tag, "Stat"),
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
        save_path=ctx.output_dir(tag, "Stat"),
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
