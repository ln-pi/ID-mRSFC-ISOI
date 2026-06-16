import numpy as np
import os
import pandas as pd

from tqdm import tqdm

from srcs.utils.logger import setup_logger
from srcs.io.readers import load_config
from srcs.utils.utils import parse_args, extract_roi
from srcs.vis.PlotSeries import PlotSeries

import warnings
warnings.filterwarnings(
    "ignore",
    message="invalid value encountered in divide"
)


args = parse_args()
config = load_config(args.config)
if config is None:
    raise ValueError("Config is empty (None). Exiting.")
logger = setup_logger(config["paths"]["log"])

processed_root = config["paths"]["processed"]
output_root = config["paths"]["output"]
mask_root = config["paths"]["mask"]

dataset = config["variables"]["dataset"]
Session_ids = config["variables"]["Session_ids"]
Mouse_ids = config["variables"]["Mouse_ids"]
Trial_ids = config["variables"]["Trial_ids"]
Atlas = config["variables"]["Atlas"]

processed_path = os.path.join(processed_root, dataset)
os.makedirs(processed_path, exist_ok=True)
output_path = os.path.join(output_root, dataset)
os.makedirs(output_path, exist_ok=True)

df = pd.DataFrame(index=Mouse_ids, columns=Trial_ids, dtype=object)

total = len(Atlas) * len(Session_ids) * len(Mouse_ids) * len(Trial_ids)
pbar = tqdm(total=total, desc="Correlating", leave=False)
for atlas in Atlas:
    atlas_tag = "_".join(atlas)
    for session in Session_ids:
        for mouse in Mouse_ids:
            for trial in Trial_ids:

                series_path = os.path.join(processed_root, session, f"{mouse}_{trial}_FiltDRR.npy")
                save_folder = os.path.join(processed_path, "FC_Matrix")
                os.makedirs(save_folder, exist_ok=True)

                matrix = np.load(series_path)
                roi_time_series, names = extract_roi(matrix, atlas_exclude=atlas, mask_root=mask_root, mode="symmetric")
                FC = np.corrcoef(roi_time_series.T)
                np.save(os.path.join(save_folder, f"{mouse}_{trial}_FC_Maticx_{atlas_tag}_excluded.npy"), FC)

                PS = PlotSeries(
                    figsize=(8, 8),
                    dpi=180,
                    save_path=os.path.join(output_path, "FC_Matrix")
                ).add_plot(
                    FC,
                    plot_type="heatmap",
                    name=f"{mouse}_{trial}_FC_Matrix_{atlas_tag}_excluded",
                    config={
                        "center": 0,
                        "vmin": -1,
                        "vmax": 1,
                        "square": True,
                        "cmap": "coolwarm",
                        "cbar": True,
                        "xticks": {
                            "ticks": np.array(range(len(names))) + 0.5,
                            "labels": [n.split("-", 1)[1] for n in names],
                            "rotation": 90,
                            "fontsize": 300 / len(names) if len(names) >=30 else 10,
                        },
                        "yticks": {
                            "ticks": np.array(range(len(names))) + 0.5,
                            "labels": [n.split("-", 1)[1] for n in names],
                            "rotation": 0,
                            "fontsize": 300 / len(names) if len(names) >=30 else 10,
                        },
                    }
                ).render(show=False, save=True, log=False)

                df.at[mouse, trial] = pd.DataFrame(FC, index=names, columns=names)

                pbar.update(1)
                pbar.set_postfix(({"atlas": atlas, "session": session, "mouse": mouse, "trial": trial}))
            
            tqdm.write(f"FC Matrix {mouse} in {session} saved.")
        
        df_output_path = os.path.join(processed_path, f"{session}_correlation_matrices_{atlas_tag}_excluded.pkl")
        df.to_pickle(df_output_path)
        tqdm.write(f">>> FC Matrices {session} packaged.")
