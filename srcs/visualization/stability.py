"""Stability reports consume saved numerical products and verify their hashes."""
from pathlib import Path
from math import ceil

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from srcs.analysis.connectome import load_arrays
from srcs.analysis.stability import SCHEMA
from srcs.analysis.topology import edge_matrix
from srcs.utils.analysis_products import verified_index
from srcs.utils.provenance import sha256
from srcs.utils.results import save_json
from srcs.visualization.descriptive import _heat, _upper


def run(ctx):
    source = Path(ctx.output_path)/"results/stability/index.json"
    index = verified_index(source)
    if index["schema"] != SCHEMA:
        raise ValueError("Unsupported stability schema")
    folder = Path(ctx.output_dir("report", "stability"))
    written = []
    def save(fig, name):
        for suffix in ("png", "pdf"):
            fig.savefig(folder/f"{name}.{suffix}", dpi=ctx.figure_dpi(), bbox_inches="tight")
            written.append(f"{name}.{suffix}")
        plt.close(fig)
    repeat = pd.read_csv(source.parent/"repeatability.csv")
    if "representation" in repeat:
        for atlas, data in repeat.groupby("atlas", sort=False):
            fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), constrained_layout=True)
            for label, part in data.groupby("representation", sort=False):
                axes[0].hist(part.pattern_r.dropna(), bins=20, alpha=.5, label=label)
                axes[1].scatter(part.pair_mean, part.mean_signed_difference, s=12, alpha=.5, label=label)
                axes[2].scatter(part.week, part.rms_difference, s=12, alpha=.5, label=label)
            axes[0].set(xlabel="Same-day pattern r", ylabel="Mouse-week pairs")
            axes[1].set(xlabel="Mean of Run edge means", ylabel="Run1 - Run2 edge mean")
            axes[1].axhline(0, color="gray", lw=.7)
            axes[2].set(xlabel="Week", ylabel="Same-day RMS difference")
            axes[2].tick_params(axis="x", rotation=90)
            for ax in axes:
                ax.legend(fontsize=7)
            fig.suptitle(f"{atlas} | paired Runs | descriptive repeatability, not independent observations")
            save(fig, f"repeatability_{atlas}")
    for product in index["products"]:
        meta, a = load_arrays(product["path"], SCHEMA)
        atlas = meta["atlas"]
        if product["kind"] == "residual":
            for r, view in enumerate(meta["view_ids"]):
                mice, weeks = meta["mouse_ids"], meta["week_ids"]
                cols = min(4, len(mice))
                fig, axes = plt.subplots(ceil(len(mice)/cols), cols, figsize=(4*cols, 3.5*ceil(len(mice)/cols)),
                                         squeeze=False, constrained_layout=True)
                for i, ax in enumerate(axes.flat):
                    if i >= len(mice):
                        ax.set_visible(False)
                        continue
                    im = _heat(ax, a["residual_similarity"][i, r], weeks, title=mice[i],
                               cmap="coolwarm", vmin=-1, vmax=1)
                fig.colorbar(im, ax=list(axes.flat), label="Delta FC pattern r", shrink=.7)
                fig.suptitle(f"{atlas} | {view} | {meta['cohort']} | same-record diagonal excluded")
                save(fig, f"delta_stability_{atlas}_{meta['cohort']}_{view}")
            continue
        fig, axes = plt.subplots(2, 3, figsize=(16, 10), constrained_layout=True)
        labels, edges = meta["roi_labels"], a["edge_index"]
        keys = ("stable_icc", "same_week_icc", "mouse_variance", "mouse_week_variance", "error_variance")
        for ax, key in zip(axes.flat, keys):
            matrix = edge_matrix(a[key], edges, len(labels))
            im = _heat(ax, matrix, labels, title=key, cmap="viridis", vmin=0,
                       vmax=1 if "icc" in key else _upper(matrix))
            fig.colorbar(im, ax=ax, shrink=.7)
        axes.flat[-1].scatter(a["mouse_variance"], a["stable_icc"], s=8, alpha=.5)
        axes.flat[-1].set(xlabel="Stable mouse variance (absolute)", ylabel="Stable mouse ICC", ylim=(-.02, 1.02))
        fig.suptitle(f"{atlas} | {meta['model_status']} | N={len(meta['mouse_ids'])} | fixed week + Run")
        save(fig, f"icc_variance_{atlas}")
        mice = meta["mouse_ids"]
        if mice:
            cols = min(4, len(mice))
            fig, axes = plt.subplots(ceil(len(mice)/cols), cols, figsize=(4*cols, 3.5*ceil(len(mice)/cols)),
                                     squeeze=False, constrained_layout=True)
            limit = _upper(np.abs(a["stable_mouse_blup"]))
            for i, ax in enumerate(axes.flat):
                if i >= len(mice):
                    ax.set_visible(False)
                    continue
                im = _heat(ax, edge_matrix(a["stable_mouse_blup"][i], edges, len(labels)), labels,
                           title=mice[i], cmap="coolwarm", vmin=-limit, vmax=limit)
            fig.colorbar(im, ax=list(axes.flat), label="Conditional stable mouse effect (Fisher-z FC)", shrink=.7)
            fig.suptitle(f"{atlas} | stable component BLUP, not pure neural identity")
            save(fig, f"stable_components_{atlas}")
    save_json(folder/"index.json", dict(schema=SCHEMA, files=written, source=str(source), source_sha256=sha256(source)))
    return [str(source)]+[item["path"] for item in index["files"]+index["products"]]
