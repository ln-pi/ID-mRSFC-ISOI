import numpy as np
import os
import pandas as pd

from tqdm import tqdm

from srcs.utils.config import load_context
from srcs.utils.utils import extract_roi
from srcs.vis.PlotSeries import save_fc_heatmap

import warnings
warnings.filterwarnings(
    "ignore",
    message="invalid value encountered in divide"
)


ctx = load_context()
logger = ctx.logger
render_options = ctx.render_options()

total = len(ctx.atlases) * len(ctx.sessions) * len(ctx.mice) * len(ctx.trials)
pbar = tqdm(total=total, desc="Correlating", leave=False)
try:
    for atlas in ctx.atlases:
        atlas_tag = "_".join(atlas)
        for session in ctx.sessions:
            df = pd.DataFrame(index=ctx.mice, columns=ctx.trials, dtype=object)
            for mouse in ctx.mice:
                for trial in ctx.trials:
                    pbar.set_postfix({"Atlas": atlas, "Session": session, "Mouse": mouse, "Trial": trial})
                    logger.info(f"\n>>> Correlating Atlas: {atlas_tag}, Session: {session}, Mouse: {mouse}, Trial: {trial}")

                    series_path = os.path.join(ctx.processed_path, session, mouse, trial, f"{mouse}_{trial}_FiltDRR.npy")
                    save_folder = ctx.processed_dir("FC_Matrix")

                    matrix = np.load(series_path)
                    roi_time_series, names = extract_roi(matrix, atlas_exclude=atlas, mask_root=ctx.mask_root, mode="symmetric")
                    FC = np.corrcoef(roi_time_series.T)
                    np.save(os.path.join(save_folder, f"{session}_{mouse}_{trial}_FC_Matrix_{atlas_tag}_excluded.npy"), FC)

                    save_fc_heatmap(
                        FC,
                        names=names,
                        save_path=ctx.output_dir("FC_Matrix"),
                        name=f"{session}_{mouse}_{trial}_FC_Matrix_{atlas_tag}_excluded",
                        **render_options,
                        config=ctx.viz_style("fc_heatmap"),
                    )

                    df.at[mouse, trial] = pd.DataFrame(FC, index=names, columns=names)

                    pbar.update(1)

                logger.info(f"FC Matrix {mouse} in {session} saved.")

            df_output_path = os.path.join(ctx.processed_path, f"{session}_correlation_matrices_{atlas_tag}_excluded.pkl")
            df.to_pickle(df_output_path)
            logger.info(f">>> FC Matrices {session} packaged.")
finally:
    pbar.close()
