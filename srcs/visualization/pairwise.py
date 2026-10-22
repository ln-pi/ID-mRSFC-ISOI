"""Independently redraw saved pairwise numerical products on shared r scales."""
import json
from math import ceil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from srcs.analysis.connectome import load_arrays, SIM_SCHEMA
from srcs.utils.provenance import sha256
from srcs.utils.results import save_json


def run(ctx):
    source = Path(ctx.output_path) / "results/pairwise/selection.json"
    selection = json.loads(source.read_text(encoding="utf-8"))
    products = {(p["atlas"], p["run"]): p for p in selection["products"]}
    folder = Path(ctx.output_dir("report", "pairwise"))
    written = []
    def save(fig, name):
        fig.savefig(folder / f"{name}.png", dpi=ctx.figure_dpi(), bbox_inches="tight")
        fig.savefig(folder / f"{name}.pdf", bbox_inches="tight")
        plt.close(fig)
        written.extend([f"{name}.png", f"{name}.pdf"])
    for (atlas, run_id), product in products.items():
        path = Path(product["path"])
        if sha256(path) != product["sha256"]:
            raise ValueError(f"Saved pairwise tensor changed: {path}; rerun pairwise")
        meta, arrays = load_arrays(path, SIM_SCHEMA)
        mice, weeks = meta["mouse_ids"], meta["week_ids"]
        cols = min(4, len(mice))
        fig, axes = plt.subplots(ceil(len(mice)/cols), cols, figsize=(4*cols, 3.6*ceil(len(mice)/cols)), squeeze=False,
                                 constrained_layout=True)
        for i, ax in enumerate(axes.flat):
            if i >= len(mice):
                ax.set_visible(False)
                continue
            matrix = arrays["S_r"][i, i].copy()
            matrix[~arrays["statistical_mask"][i, i]] = np.nan
            image = ax.imshow(matrix, vmin=-1, vmax=1, cmap="coolwarm", interpolation="none")
            ax.set_xticks(range(len(weeks)), weeks, rotation=90, fontsize=7)
            ax.set_yticks(range(len(weeks)), weeks, fontsize=7)
            ax.set_title(mice[i])
            ax.set(xlabel="week v", ylabel="week u")
        fig.suptitle(f"{atlas} | {run_id} | within-mouse Pearson r; same-record diagonal masked")
        fig.colorbar(image, ax=axes.ravel().tolist(), label="Pattern Pearson r", shrink=.7)
        save(fig, f"stability_{atlas}_{run_id}")
        for view in selection["selections"]:
            if view["atlas"] != atlas or view["run"] != run_id or view["view"] == "all_weeks" or not view["pairs"]:
                continue
            pairs = view["pairs"]
            cols = min(3, len(pairs))
            fig, axes = plt.subplots(ceil(len(pairs)/cols), cols, figsize=(4*cols, 3.8*ceil(len(pairs)/cols)),
                                     squeeze=False, constrained_layout=True)
            for k, ax in enumerate(axes.flat):
                if k >= len(pairs):
                    ax.set_visible(False)
                    continue
                pair = pairs[k]
                u, v = weeks.index(pair["week_u"]), weeks.index(pair["week_v"])
                image = ax.imshow(arrays["S_r"][:, :, u, v], vmin=-1, vmax=1, cmap="coolwarm", interpolation="none")
                ax.set_xticks(range(len(mice)), mice, rotation=90, fontsize=7)
                ax.set_yticks(range(len(mice)), mice, fontsize=7)
                ax.set(xlabel=f"week v: {pair['week_v']}", ylabel=f"week u: {pair['week_u']}")
                ax.set_title(f"{pair['week_u']} to {pair['week_v']} | gap {pair['week_gap']} weeks")
            fig.suptitle(f"{atlas} | {run_id} | {view['view']}")
            fig.colorbar(image, ax=axes.ravel().tolist(), label="Pattern Pearson r", shrink=.7)
            save(fig, f"{view['view']}_{atlas}_{run_id}")
    save_json(folder / "index.json", dict(files=written, source=str(source), source_sha256=sha256(source),
        scale="pattern_Pearson_r_of_F_z", limits=[-1, 1], nan="masked", same_record="excluded",
        layout="standard_directed_panels; rows=mouse_i_at_u, columns=mouse_j_at_v"))
    ctx.logger.info("Report: %s saved figures", len(written))
