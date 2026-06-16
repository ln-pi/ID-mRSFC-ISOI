import glob
import numpy as np
import os
import pandas as pd

from srcs.io.readers import load_config
from srcs.utils.utils import parse_args, build_half_dataframe
from srcs.utils.logger import setup_logger
from srcs.comp.ID_Tools import Identify, Edgewise_Contribution
from srcs.vis.PlotSeries import PlotSeries
from srcs.vis.circularGraph import circularGraph

import warnings
warnings.filterwarnings(
    "ignore",
    message="Mean of empty slice"
)


args = parse_args()
config = load_config(args.config)
if config is None:
    raise ValueError("Config is empty (None). Exiting.")
logger = setup_logger(config["paths"]["log"])

processed_root = config["paths"]["processed"]
output_root = config["paths"]["output"]

dataset = config["variables"]["dataset"]
Session_ids = config["variables"]["Session_ids"]
Mouse_ids = config["variables"]["Mouse_ids"]
Trial_ids = config["variables"]["Trial_ids"]
Atlas = config["variables"]["Atlas"]

ID_group_mode = config["ID"]["group_mode"]
ID_shuffle = config["ID"]["shuffle"]
ID_n_iter = config["ID"]["n_iter"]

Edge_seed = config["Edge"]["seed"]
Edge_repeat = config["Edge"]["repeat"]
Edge_threshold = config["Edge"]["threshold"]
Edge_gap = config["Edge"]["gap"]
Edge_cross_only = config["Edge"]["cross_only"]
Edge_save_mat = config["Edge"]["save_mat"]

processed_path = os.path.join(processed_root, dataset)
output_path = os.path.join(output_root, dataset)
os.makedirs(output_path, exist_ok=True)

df = pd.DataFrame(
    index=["_".join(atlas) for atlas in Atlas],
    columns=[f"{n}Trial" for n in range(1, len(Trial_ids) + 1)],
)

for trial_cnt in range(len(Trial_ids)):
    data = [[], []]
    relE = []
    labels = []

    trial_cnt += 1
    tag = f"{trial_cnt}Trial"

    for atlas in Atlas:
        label = "_".join(atlas)

        logger.info(f"\n>>> {tag} Atlas {label} excluded:")

        session_paths = glob.glob(os.path.join(processed_path, f"Session*_correlation_matrices_{label}_excluded.pkl"))

        df_out = build_half_dataframe(session_paths, trial_cnt, ID_group_mode)

        # Identify
        ID = Identify(df_out, atlas_exclude=atlas)
        info, sim_mat, Indiv, Group, relEft, cohen_d = ID.ID_Test()
        stat = ID.statistics()

        if ID_shuffle:
            ID.shuffle(
                iter=1000,
                save_path=os.path.join(output_path, "Similarity"),
                dpi=180, show=False, log=False,
            )

        ID.output(
            figsize={"line": (4, 6), "mat": (6, 6)},
            save_path=os.path.join(output_path, "Similarity"),
            dpi=180, tag=tag, show=False, log=False,
        )

        # Edgewise Contribution
        data[0].append(Group)
        data[1].append(Indiv)
        relE.append(relEft)
        df.at["_".join(atlas), tag] = [Group, Indiv, relEft, stat, cohen_d]

        EC = Edgewise_Contribution(
            path=processed_path,
            mouse_ids=Mouse_ids,
            sessions=Session_ids,
            seed=Edge_seed,
            atlas_exclude=atlas,
        )

        Phi, DP = EC.run(repeat=Edge_repeat)

        save_path = os.path.join(output_path, "EdgeContri")
        os.makedirs(save_path, exist_ok=True)

        # Group Consistency (Phi)
        avg_phi = np.mean(list(Phi.values()), axis=0)
        phi_mat = np.zeros((len(info[2]), len(info[2])))
        phi_mat[EC.edge_index] = avg_phi
        phi_mat += phi_mat.T

        circularGraph(
            phi_mat, [lab.split("-")[1] for lab in info[2]],
            node_color_mode="fixed",
            node_color_spec=info[2],
            edge_color_mode="fixed",
            edge_fixed_color="blue",
            edge_width_mode="fixed",
            threshold=np.percentile(phi_mat, Edge_threshold),
            gap=Edge_gap,
            cross_only=Edge_cross_only,
            save_path=os.path.join(save_path, f"EdgeContri_phi_{label}_excluded_thres{Edge_threshold}_{tag}"),
            save_mat=Edge_save_mat
        )

        # Differential Power (DP)
        avg_dp = np.mean(list(DP.values()), axis=0)
        dp_mat = np.zeros((len(info[2]), len(info[2])))
        dp_mat[EC.edge_index] = avg_dp
        dp_mat += dp_mat.T

        circularGraph(
            dp_mat, [lab.split("-")[1] for lab in info[2]],
            node_color_mode="fixed",
            node_color_spec=info[2],
            edge_color_mode="fixed",
            edge_fixed_color="red",
            edge_width_mode="fixed",
            threshold=np.percentile(dp_mat, Edge_threshold),
            gap=Edge_gap,
            cross_only=Edge_cross_only,
            save_path=os.path.join(save_path, f"EdgeContri_dp_{label}_excluded_thres{Edge_threshold}_{tag}"),
            save_mat=Edge_save_mat
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
    
    PlotSeries(figsize=(4, 4), dpi=180, save_path=os.path.join(output_path, "Stat")).add_plot(
        relE,
        plot_type="box",
        name=f"boxmap_relEffect_{tag}",
        config = {
            "box_width": 0.1,
            "box_alpha": 0.7,
            "box_color": None,
            "showfliers": False,
            "fontsize": 7,
            "rotation": 0,
            "grid": False,
            "xticklabel": labels,
            "show_legend": False,
            "show_sub_label": True,
            "hide_ns": True,
        }
    ).render(show=False, save=True, log=False)

    PlotSeries(figsize=(4, 4), dpi=180, save_path=os.path.join(output_path, "Stat")).add_plot(
        data,
        plot_type="box",
        name=f"boxmap_Group-Indiv_{tag}",
        config={
            "box_width": 0.1,
            "box_alpha": 0.7,
            "box_color": None,
            "showfliers": False,
            "fontsize": 7,
            "rotation": 0,
            "group_gap": 2,
            "grid": False,
            "xticklabel": labels,
            "hide_ns": True,
        }
    ).render(show=False, save=True, log=False)

df_save_path = os.path.join(output_path, "statistics.pkl")
df.to_pickle(df_save_path)
logger.info("\nStatistics done output.")
